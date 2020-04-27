#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
#ifdef TIPSYFORMAT
#include "tipsydefs.h"
#else
#ifdef TIPSYN2FORMAT
#include "tipsydefs_n2.h"
#else
#ifdef OWLSFORMAT 
#include "owls.h"
#else
#include "tipsydefs_n.h" 
#endif // OWLSFORMAT
#endif // TIPSYN2FORMAT 
#endif // TIPSYFORMAT
#include "specexbindefs.h"
#include "proto.h"
#include "extern.h"

#define RHOLOW -9.0
#define TEMPLOW 2.5

extern double unit_Mass,unit_Velocity,unit_Length,unit_Density;

#ifdef INTKERNELNHLIMIT
#define NINTERP 10000
#define NSRCHRAD 1.0
#define NHILIM 1.58489e17
//#define NHILIM 1.e18	// agrees better w/Faucher-Giguere,Keres 2010 Fig 3
#define P0BLITZ 3.5e+04 // Blitz & Rosolowski 2006
#define ALPHA0BLITZ 0.92 // Blitz & Rosolowski 2006
float KernIntTable[NINTERP+1][2];
#endif

#define clight		2.99792458e10
#define KPC 3.086e21

float findrandgauss();

int ContSmoothSpec()
{
  int i,j,k,m;
  struct spec_particle *cp;
  double irep[NDIM],part_pos[NDIM];
  float bound_min[NDIM],bound_max[NDIM];
  double vz;
  float CosmicTime();
  int cosmopar();
  float distnorm;
  float unit_vel;
  int comove=1;
  float kernel,vkernel;
  float zlower,zupper;
  float abs_zlower,abs_zupper; 
  float radius2,radius;
  float zo2,zo,zi2,zi;
  float d,d2;
  int bin,bin_min,bin_max;
  double zcoord;
  double rhocgs;
#ifndef PIPELINE
  FILE *binfile, *partfile;
  char binname[80], partname[80];
#ifndef OUTPUT_LOCAL_FOLDER
  char longbinname[400], longpartname[400];
  char spec_dir[100];
  strcpy(spec_dir, "/scratch/shuiyao/los/shortlos/");
#endif
#endif
  float IonFrac(),*ionfrac;
  int load_fraction_tables();
  double ion_weight;
#ifdef NONEQUIL
  int load_noneq_fraction_tables();
  float Noneq_Frac();
#endif
#ifdef VARGALBKGD
  float mgal_100kpc;
#endif
#ifdef PAINTAVERAGEMETALS		
  int iZ,jZ;
  float mu;
#endif
#ifdef INTKERNELNHLIMIT
  int InitKernIntTable();
  int ilo,ihi;
  float frh;
  float coldphasemassfrac,coldphasetemp;

  InitKernIntTable();
#endif
#ifdef PHEW
  int nvbinsnew;
  int ionid;
  int irepz;
  int vbin, vbin_min, vbin_max;
  
  double prob0, prob1, prob2;
  double nc_bin;
  float nc_this_bin;
  int Zcol;
  double *vbin_size, *vbin_coord;
  double vcoord, vstep;
  float hubble_expansion;
  float vlower,vupper;
  float abs_vlower,abs_vupper;
  float dvcol, unit_col;
  float b;
  double colcloud, rhocloud, tcloud, vcloud;
  double *voffset;
  ionStruct I;
#endif  

  i = floor(theta*180.0/PI+0.001);
  k = floor(redshift_begin*100+0.0001);
  //k = floor((redshift_low-delta_redshift+0.01)*100+0.0001);
  j = floor(redshift_end*100+0.0001);

  //#ifdef PARTIONFRAC
  redshift_hold = redshift;
  redshift_center = bin_redshift[(int)((nzbins)/2)];
  redshift = redshift_center;
  load_fraction_tables();
#ifdef NONEQUIL
  load_noneq_fraction_tables();
#endif

  redshift = redshift_hold;
  ionfrac = (float *) malloc(nions*sizeof(float));
  //#endif

#ifndef PIPELINE
#ifdef SHORTSPEC
  sprintf(binname,"binzfile.%s.%s",sim_id,namesuffix);
  sprintf(partname,"partzfile.%s.%s",sim_id,namesuffix);
#ifdef OUTPUT_LOCAL_FOLDER
  binfile = fopen(binname,"w");
  partfile = fopen(partname,"w");
#else
  strcat(strcpy(longbinname, spec_dir), binname);
  strcat(strcpy(longpartname, spec_dir), partname);
  binfile = fopen(longbinname,"w");
  partfile = fopen(longpartname,"w");
#endif
  fprintf(partfile,"#count = %d\n",count);
#else 
  if(theta>0){
    sprintf(binname,"binzfile.%s.%d.%d_%d",sim_id,i,j,k);
  }else{
    sprintf(binname,"binzfile.%s.%s.%d_%d",sim_id,id,j,k);
  }
  if(theta>0){
    sprintf(partname,"partzfile.%s.%d.%d_%d",sim_id,i,j,k);
  }else{
    sprintf(partname,"partzfile.%s.%s.%d_%d",sim_id,id,j,k);
  }
#ifdef OUTPUT_LOCAL_FOLDER
  binfile = fopen(binname,"a");
  partfile = fopen(partname,"a");
#else
  strcat(strcpy(longbinname, spec_dir), binname);
  strcat(strcpy(longpartname, spec_dir), partname);
  binfile = fopen(longbinname,"a");
  partfile = fopen(longpartname,"a");
#endif

  fprintf(binfile,"# xbegin = % 7.5f ybegin = % 7.5f zbegin = % 7.5f theta = %7.3f phi = %7.3f\n",xbegin,ybegin,zbegin,theta*180.0/PI,phi*180.0/PI);
  fprintf(stderr,"bin_size[0] = %5.3e nzbins = %d zlength %5.3e\n",bin_size[0],nzbins,zlength);

#endif
#endif

  //zout = redshift_hold;
#ifdef PHEW
  /* SH161020: Construct vbins for particle-by-particle tau calculation*/
  // --------------------------------
  // nzbins is are already assigned in GetSpecParticles.
  // On the other hand, nzloopbins and nvloopbins hasn't been updated at this point.
  // However, the IonTotal.redshift[nzloopbins ~ nzloopbins+nzbins] has not been constructed
  // The values are in bin_redshift, which is updated in getspecparticles.c
  // WARNING: In the first loop, IonTotal.redshift has not been updated.

  // voffset[nzbins+nzloopbins] corresponds to the hubble_expansion in tau.c
  voffset = malloc((nzbins+nzloopbins) * sizeof(double));
  if(IonTotal.redshift[0] > 0) // The first element
    voffset[0] = hubble*aex*unit_Velocity/1.e5*IonTotal.binsize[0];
  else
    voffset[0] = hubble*aex*unit_Velocity/1.e5*bin_size[0];    
  for(i = 1; i < nzloopbins; i ++){ // The old wraps
#ifndef SHORTSPEC
    cosmopar(CosmicTime(IonTotal.redshift[i]));
#endif
    voffset[i] = voffset[i-1] + hubble*aex*unit_Velocity/1.e5*IonTotal.binsize[i];
  }
  for(i = nzloopbins; i < nzloopbins + nzbins; i ++){ // The new wrap
#ifndef SHORTSPEC
    cosmopar(CosmicTime(bin_redshift[i]));
#endif
    voffset[i] = voffset[i-1] + hubble*aex*unit_Velocity/1.e5*bin_size[i];
  }
  if(IonTotal.redshift[0] > 0)
    redshift_track = IonTotal.redshift[0];
  else
    redshift_track = bin_redshift[0];
  vbin_size = malloc(sizeof(double));
  vbin_coord = malloc(sizeof(double));
  bin = 0; //
  vcoord = 0;
  zcoord = 0;
  /* while(redshift_track >= IonTotal.redshift[nzbins-1]){ */
  while(redshift_track >= bin_redshift[nzbins-1]){
    redshift_track -= vres; // vres = VRES is defined in specexbindefs.h
#ifndef SHORTSPEC  
    cosmopar(CosmicTime(redshift_track));
#endif
    zstep = vres/(BOXSIZE*hubble*unit_Velocity/clight);
    zcoord += zstep;
    vstep = zstep*aex*hubble*unit_Velocity;
    vcoord += vstep;
    bin++;
    vbin_size = realloc(vbin_size, bin*sizeof(double));
    vbin_coord = realloc(vbin_coord, bin*sizeof(double));
    vbin_size[bin-1] = vstep/1.e5;
    vbin_coord[bin-1] = vcoord/1.e5;
  }
  // bin does not necessarily equal to nzloopbins + nzbins, because
  // vres and zres are different thing.
  bin --;
  if(bin != nvloopbins + nvbins)
    fprintf(stderr, "Error: bin(%d) != nvloopbins(%d) + nvbins(%d)\n",
	    bin, nvloopbins, nvbins);
  nvbinsnew = nvbins + nvloopbins;
#endif // PHEW
  
#ifdef SHORTSPEC
  redshift = redshift_center;
  redshift_hold = redshift_center;
  cosmopar(CosmicTime(redshift_center));
#else
  cosmopar(CosmicTime(redshift));
#endif
  bound_min[0] = bound_min[1] = bound_min[2] = -HALFBOX;
  bound_max[0] = bound_max[1] = bound_max[2] = HALFBOX;
  for(i = nzloopbins; i < nzloopbins+nzbins; i++){
    IonTotal.mass[i] = IonTotal.vel[i] = IonTotal.temp[i] = IonTotal.rho[i] = 0.0;
    for(m=0;m<NMETALS;m++) IonTotal.metals[m][i] = 0;
#ifdef PHYSSPEC
    IonTotal.sfr[i] = IonTotal.wtmass[i] = IonTotal.mgal[i] = IonTotal.dgal[i] = IonTotal.age[i] = IonTotal.nrec[i] = IonTotal.vlaunch[i] = 0.0;
#endif
    for( k=0; k<nions; k++ ){ 
      Ion[k].mass[i] = Ion[k].vel[i] = Ion[k].temp[i] = Ion[k].rho[i] = 0.0;
      for(m=0;m<NMETALS;m++) Ion[k].metals[m][i] = 0;
#ifdef PHYSSPEC
      Ion[k].sfr[i] = Ion[k].wtmass[i] = Ion[k].mgal[i] = Ion[k].dgal[i] = Ion[k].age[i] = Ion[k].nrec[i] = Ion[k].vlaunch[i] = 0.0;
#endif	    
    }
  }
  for (i = 0 ;i < count ;i++) {
    cp = spec_particles+i;
    for(irep[0] = -1; irep[0] <= 1; irep[0]++) {
      for(irep[1] = -1; irep[1] <= 1; irep[1]++) {
	for(irep[2] = -1; irep[2] <= 1; irep[2]++) {
	  if(irep[0]*BOXSIZE + cp->pos[0] > bound_max[0]
	     && irep[0]*BOXSIZE+cp->pos[0]-2*cp->hsmooth
	     > bound_max[0])
	    continue;
	  if(irep[0]*BOXSIZE + cp->pos[0] < bound_min[0]
	     && irep[0]*BOXSIZE+cp->pos[0]+2*cp->hsmooth
	     < bound_min[0])
	    continue;
	  if(irep[1]*BOXSIZE + cp->pos[1] > bound_max[1]
	     && irep[1]*BOXSIZE+cp->pos[1]-2*cp->hsmooth
	     > bound_max[1])
	    continue;
	  if(irep[1]*BOXSIZE + cp->pos[1] < bound_min[1]
	     && irep[1]*BOXSIZE+cp->pos[1]+2*cp->hsmooth
	     < bound_min[1])
	    continue;
	  if(irep[2]*BOXSIZE + cp->pos[2] > bound_max[2]
	     && irep[2]*BOXSIZE+cp->pos[2]-2*cp->hsmooth
	     > bound_max[2])
	    continue;
	  if(irep[2]*BOXSIZE + cp->pos[2] < bound_min[2]
	     && irep[2]*BOXSIZE+cp->pos[2]+2*cp->hsmooth
	     < bound_min[2])
	    continue;
	  for( k=0; k<NDIM; k++ ) part_pos[k] = cp->pos[k]+irep[k]*BOXSIZE;
	  RotateCoords(part_pos[0], part_pos[1], part_pos[2], theta, phi);
	  distnorm = 1. / (cp->hsmooth * cp->hsmooth) ;
	  radius2 = (pow(xrot-xrotline,2)+pow(yrot-yrotline,2))*distnorm;
	  if(radius2 >= 4.0) continue;
	  radius = sqrt(radius2) ;
	  zo2 = 4. - radius2 ;
	  zo = sqrt(zo2) ;
	  if(radius2 < 1.){
	    zi2 = 1. - radius2 ;
	    zi = sqrt(zi2) ;
	  }
	  else{
	    zi2 = zi = 0.0;
	  }
	  zcoord = zrot; // hold zrot because we need to do velocity

	  RotateCoords(cp->vel[0],cp->vel[1],cp->vel[2],theta,phi);
	  vz = zrot;
#ifdef METALFLOOR
	  cp->metals[0] += 0.001*0.00213;
	  cp->metals[1] += 0.001*0.00541;
	  cp->metals[2] += 0.001*0.00067;
	  cp->metals[3] += 0.001*0.00117;
#endif

#ifdef PAINTAVEMETALLICITY
	  cp->metals[0] = PAINTAVEMETALLICITY*0.00213;
	  cp->metals[1] = PAINTAVEMETALLICITY*0.00541;
	  cp->metals[2] = PAINTAVEMETALLICITY*0.00067;
	  cp->metals[3] = PAINTAVEMETALLICITY*0.00117;
#endif

#ifdef INTKERNELNHLIMIT
	  coldphasemassfrac = 1;
	  coldphasetemp = cp->temp;
	  if(cp->sfr>0){
	    //coldphasemassfrac = (1e+08-cp->temp)/1e+08; // WRONG TO ASSUME TEMPERATURE 1E+08 FOR HOT PHASE, IS A COMPLICATED FUNCTION. 
	    coldphasemassfrac = 0.9; // MOTIVATED BY HONG ET AL.  TO GET AROUND COMPLICATED FUNCTION OF TEMPERATURE. 
	    coldphasetemp = 1e+03; /* Modified 7-28-11 to split SH03 two phase. 10^3 and 10^8 */
	  }
	  // SH161007: PHEW?
#endif

#ifdef PAINTAVERAGEMETALS // Mostly OFF
	  iZ = (log10(cp->rho*XH*unit_Density/(aex*aex*aex)/MHYDR)-RHOLOW)*10;
	  jZ = (log10(cp->temp)-TEMPLOW)*10;
	  for(m=0;m<NMETALS;m++){
	    //cp->metals[m] = rhoTZ[iZ][jZ][m]/53.32*pow(10,1.0*findrandgauss());
	    mu = log(rhoTZ[iZ][jZ][m])-pow(PAINTAVERAGEMETALS,2)/2;
	    cp->metals[m] = exp(mu+PAINTAVERAGEMETALS*findrandgauss());
	    //cp->metals[m] = rhoTZ[iZ][jZ][m];
	  }
	  printf("iZ= %3d jZ= %3d rho= % 5.3e temp= % 5.3e metals= % 5.3e %5.3e %5.3e %5.3e\n",iZ,jZ,log10(cp->rho*XH*unit_Density/(aex*aex*aex)/MHYDR),log10(cp->temp),cp->metals[0],cp->metals[1],cp->metals[2],cp->metals[3]);
#endif

#ifndef PIPELINE
	  fprintf(partfile,"% 6.4e % 6.4e % 6.4e %6.4e %6.4e %6.4e %6.4e ",cp->mass, cp->rho*XH*unit_Density/(aex*aex*aex)/MHYDR, cp->temp, cp->metals[0], cp->metals[1], cp->metals[2], cp->hsmooth);
#ifdef PHYSSPEC
	  fprintf(partfile,"%5.2f % 4.2f %5.2f ",cp->mgal,cp->dgal,cp->ageaway);
#else 
	  fprintf(partfile,"% 6.4e % 6.4e % 6.4e ",cp->pos[0],cp->pos[1],cp->pos[2]);
	  fprintf(partfile,"% 6.4e % 6.4e % 6.4e ",cp->vel[0],cp->vel[1],cp->vel[2]);

#endif
	  fprintf(partfile," %5.3f %6.4e\n", redshift, IonFrac(cp->temp,cp->rho*XH*unit_Density/(aex*aex*aex),0)); 
#endif


	  bin_min = binarysearch((zcoord - zo*cp->hsmooth - zbeginline),bin_coord,nzbins);
	  bin_max = binarysearch((zcoord + zo*cp->hsmooth - zbeginline),bin_coord,nzbins);

	  //fprintf(stdout,"BIN: irep= % 5.3f zcoord= %7.5f zhsmooth= %7.5f zbeginline= %7.5f bin_coord[0]= %5.3e bin_min= %5d bin_max= %5d nzbins= %5d %5d ratio= %5.3e\n",irep[2],zcoord, zo*cp->hsmooth, zbeginline,bin_coord[0],bin_min,bin_max,bin_max-bin_min,nzbins,(bin_max-bin_min)/(zo*cp->hsmooth));
	  fflush(stdout);
	  //#ifdef PARTIONFRAC
#ifndef SHORTSPEC
	  if (bin_min < 0) bin_min = 0;
	  if (bin_max >= nzbins) bin_max = nzbins - 1;

	  cosmopar(CosmicTime(bin_redshift[(int)((bin_max-bin_min)/2)]));
#endif

	  rhocgs = cp->rho*XH*unit_Density/(aex*aex*aex);

	  for( k=0; k<nions; k++ ){
#if defined(NONEQUIL) && defined(DO6IONS) 
	    if((k==1 && log10(rhocgs/MHYDR)>-3.20) || (k==2 && log10(rhocgs/MHYDR)>-3.30) || (k==4 && log10(rhocgs/MHYDR)>-4.10) || (k==5 && log10(rhocgs/MHYDR)>-2.35)){
	      ionfrac[k] = Noneq_Frac(cp->temp,(cp->metals[0]+cp->metals[1]+cp->metals[2]+cp->metals[3])*1.28,k);
	      printf("% 5.3f\n",log10(rhocgs/MHYDR));
	    }else{
	      ionfrac[k] = IonFrac(cp->temp,rhocgs,k);
	    }
#else		  
#ifdef VARGALBKGD
	    if(cp->mgal>6){
	      mgal_100kpc = log10(pow(10,cp->mgal)*pow(100/cp->dgal,2));
	    }else{
	      mgal_100kpc = 6.0;
	    }
	    if(mgal_100kpc>11.0){
	      mgal_100kpc = 11.0;
	    }

	    ionfrac[k] = IonFrac(cp->temp,rhocgs,mgal_100kpc,k);
#else
#ifdef TEMPOVER105
	    if(cp->temp>=1e+05){
	      ionfrac[k] = IonFrac(cp->temp,rhocgs,k);
	    }else{
	      ionfrac[k] = 0;
	    }
#else
#ifdef TEMPOVER1052
	    if(cp->temp>=1.585e+05){
	      ionfrac[k] = IonFrac(cp->temp,rhocgs,k);
	    }else{
	      ionfrac[k] = 0;
	    }
#else
	    ionfrac[k] = IonFrac(cp->temp,rhocgs,k);
#ifdef INTKERNELNHLIMIT
#ifdef DO9IONS
	    if(k==0 || k==7) // For HI and MgII and assuming i9 
#else
	      if(k==0 || k==3 || k==10 || k==20 || k==22 || k==24 || k==30) // For HI, CII, OI, MgII, SiII, AlII, FeII and using i31 
#endif
		{
		  ionfrac[k] = IonFrac(coldphasetemp,rhocgs,k);
		  ilo = 0; ihi = NINTERP-1;
		  frh = ionfrac[k]*rhocgs/MHYDR*cp->hsmooth*unit_Length*aex; // ???
		  while( ihi-ilo > 1 ) {
		    if( KernIntTable[(ilo+ihi)/2][1]*frh < NHILIM ) ihi = (ilo+ihi)/2;
		    else ilo = (ilo+ihi)/2;
		  }
		  if(ilo>0)ionfrac[k] = coldphasemassfrac*(ionfrac[k]*(KernIntTable[(ilo+ihi)/2][0]) + 1/(1+pow((rhocgs/MHYDR*coldphasetemp)/P0BLITZ,ALPHA0BLITZ)) *(1-KernIntTable[(ilo+ihi)/2][0]));		    
		  //if( KernIntTable[(ilo+ihi)/2][0] < 0.999999 ) 
		  //fprintf(stdout,"INTKERNELNHLIMIT: % 5.3f %5.3f % 5.3f %6.2f %5d %5d %5.3e %7.5f % 7.5f % 5.3e % 5.3e % 5.3e %5.3e % 5.3e % 5.3e % 5.3e %5.3e\n",log10(rhocgs/MHYDR),log10(cp->temp),log10(ionfrac[k]),log10(frh),ilo,ihi,KernIntTable[(ilo+ihi)/2][1]*frh,KernIntTable[(ilo+ihi)/2][0],KernIntTable[(ilo+ihi)/2][1],log10(1-KernIntTable[(ilo+ihi)/2][1]),cp->hsmooth*unit_Length,cp->hsmooth,cp->sfr,coldphasemassfrac,1/(1+pow((rhocgs/MHYDR*coldphasetemp)/P0BLITZ,ALPHA0BLITZ)),(rhocgs/MHYDR*cp->temp),IonFrac(coldphasetemp,rhocgs,k));
		}

#endif // INTKERNELNHLIMIT On

#ifdef NHLIMIT // Not used anymore.  
	    if(k==0 || k==7){ /* Now only H and MgII in i9 */
	      if(rhocgs/MHYDR*cp->temp>155 && rhocgs/MHYDR*cp->temp<810 && rhocgs/MHYDR>2.176e-08*pow(cp->temp,1.1)){
		ionfrac[k] = 1; 
	      }
	      //if(rhocgs/MHYDR*cp->temp>155)fprintf(stderr,"PARTSMOOTH: %d % 5.3e %5.3e %5.3f \n",i,rhocgs/MHYDR, cp->temp,rhocgs/MHYDR*cp->temp);
	    }
#endif // NHLIMIT On
#endif // TEMPOVER1052 Off
#endif // TEMPOVER105 Off
#endif // VARGALBKGD Off
#endif // NONEQUIL && DO6ION Off
	  }

	  if((bin_min <= 0 && bin_max == 0) || (bin_min >= nzbins-1 && bin_max >= nzbins-1)) continue;

	  for(bin = bin_min; bin <= bin_max; bin++){
	    if(bin_min >= bin_max) fprintf(stderr,"ALERT!!! bin_min (%d) >= bin_max (%d)\n",bin_min,bin_max);
#ifndef SHORTSPEC
	    if(bin >= nzbins) 
	      cosmopar(CosmicTime(bin_redshift[nzbins-1]));
	    else if(bin < bin_min) 
	      cosmopar(CosmicTime(bin_redshift[0]));
	    else
	      cosmopar(CosmicTime(bin_redshift[bin]));
#endif
	    zlower = (bin_coord[bin] + zbeginline - zcoord)/cp->hsmooth;
	    zupper = (bin_coord[bin+1] + zbeginline - zcoord)/cp->hsmooth;

	    /* SH161007: The kernel weight is an integration along the bin. The integration interval depends on whether the bin is fully submerged inside the SPH kernel. Note that h_specexbin = 0.5 h_gadget. A cubic spline kernel (3D) is used. Since it has two regions devided by u = 0.5, the integration must be done separately if the bin happens to cross over the transition. */

	    /* SH161007: The PhEW particles: Calculate the kernel as normal SPH particle. Then, assuming Poisson distribution of clouds, with average number density of cloud being nc_bin = ncloud * kernel * pi * (rcloud)^2. The determine the number of clouds inside the bin from the probabilistic distribution F_p(k, nc_bin). I guess nc_bin is very small so we only care about k=0 and k=1 case. Though I output k=2 prob just to make sure it is negligible.  */

	    kernel = 0.0 ;
	    vkernel = 0.0 ;
	    abs_zlower = fabs(zlower) ;
	    if(abs_zlower >= zo){
	      kernel += 1.3125*radius2*zo - 1.5*zo + 
		0.5*zo2*zo - (1.5*radius2 +
			      0.09375*radius2*radius2)*log(zo + 2.0);
	      /*kernel at zo in region 2 */
	      if(comove == 1){
		vkernel -= zo2 + 0.75*radius2*zo2 +
		  0.375*zo2*zo2 - 9.6 ;
	      }
	    }
	    else if(abs_zlower < zo && (abs_zlower > zi ||
					zlower == zi)){
	      d2 = radius2 + zlower*zlower ;
	      d = sqrt(d2) ;
	      kernel -= copysign(1.0,zlower)*((2.0 +
					       1.5*radius2)* abs_zlower + 
					      0.5*abs_zlower*abs_zlower*abs_zlower -
					      0.0625*abs_zlower*d*d2 - (1.5 + 0.09375*
									radius2)* abs_zlower*d - (1.5*radius2 +
												  0.09375*radius2* radius2)*log(abs_zlower + d));
	      /* sign(zlower)*kernel at abs_zlower in region 2 */
	      if(comove == 1){
		vkernel -= zlower*zlower*(1. +
					  0.75*radius2) + 0.375*zlower*
		  zlower*zlower*zlower - d*d2 -
		  0.05*d2*d2*d ;
	      }
	    }
	    else {
	      d2 = radius2 + zlower*zlower ;
	      d = sqrt(d2) ;
	      kernel -= copysign(1.0,zlower)*((1.0 -
					       1.5*radius2)*abs_zlower - 0.5*
					      abs_zlower*abs_zlower*abs_zlower +
					      0.1875*abs_zlower*d2*d + 0.28125*
					      radius2*abs_zlower*d + 0.28125*radius2*
					      radius2*log(abs_zlower + d)) ;
	      /* sign(zlower)* kernel at abs_zlower in region 1 */
	      if(comove == 1){
		vkernel -= (0.5 - 0.75*radius2)*zlower*
		  zlower - 0.375*zlower*zlower*
		  zlower*zlower + 0.15*d2*d2*d ;
	      }
	    }
	    abs_zupper = fabs(zupper) ;
	    if(zlower*zupper <= 0.0){  /* bin stradles zero */
	      if(radius2 < 1.){
		kernel -= 0.5625*radius2*radius2*
		  log(radius) ;
		/* 2.0 * kernel at zero in region 1 */
	      }
	      else {
		kernel -= -(3.0*radius2
			    + 0.1875*radius2*radius2)*
		  log(radius) ;
		/* 2.0 *kernel at zero in region 2 */
	      }
	    }

	    if(radius2 < 1.){
	      if(zlower*zupper < 0.0 && abs_zupper > zi &&
		 abs_zlower > zi) {
		kernel += 2.375*zi - 2.4375*radius2*zi -
		  zi2*zi + 0.5625*radius2*radius2*
		  log(zi+1) ;
		/* 2.0 * kernel at zi in region 1 */
		kernel -= 0.875*zi + 2.8125*radius2*zi +
		  zi2*zi - (3.0*radius2 + 0.1875*
			    radius2*radius2)*log(zi + 1) ;
		/* 2.0 * kernel at zi in region 2 */
	      }
	      else {
		if((zlower < -zi && zupper > -zi) ||
		   (zlower < zi && zupper > zi)){
		  kernel += 1.1875*zi - 1.21875*radius2*
		    zi - 0.5*zi2*zi + 0.28125*
		    radius2*radius2*log(zi+1) ;
		  /* kernel at zi in region 1 */
		  if(comove == 1){
		    vkernel += copysign(1.0, zlower)
		      *((0.5 - 0.75*radius2)*
			zi2 - 0.375*zi2*zi2 + 0.15) ;
		  }
		  kernel -= 0.4375*zi + 1.40625*radius2*
		    zi + 0.5*zi2*zi - (1.5*radius2 +
				       0.09375*radius2*radius2)*
		    log(zi + 1) ;
		  /* kernel at zi in region 2 */
		  if(comove == 1){
		    vkernel -= copysign(1.0, zlower)
		      *((1.0 + 0.75*radius2)*
			zi2 + 0.375*zi2*zi2 - 1.05) ;
		  }
		}
	      }
	    }
		  
	    if(abs_zupper >= zo){
	      kernel += 1.3125*radius2*zo - 1.5*zo + 
		0.5*zo2*zo - (1.5*radius2 +
			      0.09375*radius2*radius2)*log(zo + 2.0) ;
	      /* kernel at zo in region 2 */
	      if(comove == 1){
		vkernel += zo2 + 0.75*radius2*zo2 +
		  0.375*zo2*zo2 - 9.6 ;
	      }
	    }
	    else if(abs_zupper < zo && abs_zupper >= zi){
	      d2 = radius2 + zupper*zupper ;
	      d = sqrt(d2) ;
	      kernel += copysign(1.0,zupper)*((2.0 + 1.5*
					       radius2)*abs_zupper + 0.5*abs_zupper*
					      abs_zupper*abs_zupper - 0.0625*
					      abs_zupper*d*d2 - (1.5 + 0.09375*
								 radius2)* abs_zupper*d - (1.5*radius2 +
											   0.09375*radius2* radius2)*
					      log(abs_zupper + d)) ;
	      /* sign(zupper)*kernel at abs_zupper in region 2 */
	      if(comove == 1){
		vkernel += zupper*zupper*(1. +
					  0.75*radius2) + 0.375*zupper*
		  zupper*zupper*zupper - d*d2 -
		  0.05*d2*d2*d ;
	      }
	    }
	    else {
	      d2 = radius2 + zupper*zupper ;
	      d = sqrt(d2) ;
	      kernel += copysign(1.0,zupper)*((1.0 -
					       1.5*radius2)*abs_zupper - 0.5*
					      abs_zupper*abs_zupper*abs_zupper +
					      0.1875*abs_zupper*d2*d + 0.28125*
					      radius2*abs_zupper*d + 0.28125*radius2*
					      radius2*log(abs_zupper + d)) ;
	      /* sign(zupper)* kernel at abs_zupper in region 1 */
	      if(comove == 1){
		vkernel += (0.5 - 0.75*radius2)*zupper*
		  zupper - 0.375*zupper*zupper*
		  zupper*zupper + 0.15*d2*d2*d ;
	      }
	    }
	    kernel *= distnorm/PI ;
	    if(comove == 1){
	      vkernel /= PI*cp->hsmooth ;
	    }

	    //unit_vel = 1;
	    //#ifndef OWLSFORMAT /* Not sure about this... */
	    unit_vel = unit_Velocity*aex/1.e5;
	    //#endif
	    //if(kernel > -1.0){ /* It is an unexplained occurence why some kernel values are anomolously very negative... but they are rare so we will ignore them */
	    if(kernel > 0.0){ /* Really, negative kernel values are bad.  They should be ignored at all costs.  7-9-11 */
#ifdef PHEW
	      if(cp -> wind_flag <= 0){ /* Not a PhEW wind */
#endif		    
		IonTotal.mass[bin+nzloopbins] += kernel*cp->mass ;
		//if(bin==nzbins-2 || bin==nzbins-1){ 
		//fprintf(stdout,"LASTBIN: nzbins= %5d kernel= % 5.3e cp->mass= %5.3e mass= %5.3e tot= %5.3e zlower= % 5.3e zupper= % 5.3e radius2= % 5.3e radius= % 5.3e\n",bin,kernel,cp->mass,kernel*cp->mass,IonTotal.mass[bin+nzloopbins],zlower,zupper,radius2,radius);
		//}
		//if(bin%3000==0)fprintf(stderr,"(%d)= %g %g %g\n",bin+nzloopbins,IonTotal.vel[bin+nzloopbins],kernel,cp->mass);
#ifdef ZEROVEL
		IonTotal.vel[bin+nzloopbins] = 0;
#else 
		IonTotal.vel[bin+nzloopbins] += kernel*(cp->mass)*unit_vel*vz ;
#endif
		if(comove == 1){
#ifndef ZEROVEL
		  IonTotal.vel[bin+nzloopbins] += vkernel*(cp->mass)*hubble*unit_vel;
#endif
		}
		IonTotal.temp[bin+nzloopbins] += kernel*(cp->mass)*(cp->temp) ;
		IonTotal.rho[bin+nzloopbins] += kernel*(cp->mass)*(cp->rho*XH*unit_Density/(aex*aex*aex)) ;

		for(m=0;m<NMETALS;m++){
		  IonTotal.metals[m][bin+nzloopbins] += kernel*(cp->mass)*(cp->metals[m]) ;
		  if(cp->metals[m] > 10) {printf("BAD CP METAL!!\n"); exit(-1);}
		}
#ifdef PHYSSPEC
		IonTotal.sfr[bin+nzloopbins] += kernel*(cp->mass)*(cp->sfr);
		if(cp->mgal>0 && cp->dgal>0 && cp->ageaway>0){
		  IonTotal.wtmass[bin+nzloopbins] += kernel*cp->mass ;
		  IonTotal.age[bin+nzloopbins] += kernel*cp->mass*(cp->ageaway) ;
		  IonTotal.dgal[bin+nzloopbins] += kernel*cp->mass*(log10(cp->dgal)) ;
		  IonTotal.mgal[bin+nzloopbins] += kernel*cp->mass*(cp->mgal) ;
		  IonTotal.vlaunch[bin+nzloopbins] += kernel*cp->mass*(cp->vlaunch) ;
		  IonTotal.nrec[bin+nzloopbins] += kernel*cp->mass*(cp->nrec) ;
		}

#endif

		for( k=0; k<nions; k++ ) {		      
		  if(Ion[k].Zcolumn==-1){
		    ion_weight = ionfrac[k];
		  }else{
		    if(Ion[k].Zcolumn<-1){
		      //ion_weight = ionfrac[k]*cp->metals[3]/0.001267*Ion[k].fraction*pow(10,Ion[k].alpha);
		      ion_weight = ionfrac[k]*cp->metals[1]/0.009618*Ion[k].fraction*pow(10,Ion[k].alpha); /* 2-11-10 */
		    }else{
		      ion_weight = ionfrac[k]*cp->metals[Ion[k].Zcolumn];
		    }
		  }
		      
		      
		  Ion[k].mass[bin+nzloopbins] += kernel*ion_weight*cp->mass;


#ifdef ZEROVEL
		  Ion[k].vel[bin+nzloopbins] = 0;
#else
		  Ion[k].vel[bin+nzloopbins] += kernel*ion_weight*(cp->mass)*vz*unit_vel;
#endif
		  if(comove == 1){
		    Ion[k].vel[bin+nzloopbins] += vkernel*ion_weight*(cp->mass)*hubble*unit_vel;
		  }
		      
		  Ion[k].temp[bin+nzloopbins] += kernel*ion_weight*cp->mass*cp->temp ;
		  Ion[k].rho[bin+nzloopbins] += kernel*ion_weight*cp->mass*(cp->rho*XH*unit_Density/(aex*aex*aex)) ;
		  for(m=0;m<NMETALS;m++)Ion[k].metals[m][bin+nzloopbins] += kernel*ion_weight*(cp->mass)*(cp->metals[m]) ;
#ifdef PHYSSPEC
		  Ion[k].sfr[bin+nzloopbins] += kernel*ion_weight*(cp->mass)*(cp->sfr);
		  if(cp->mgal>0 && cp->dgal>0 && cp->ageaway>0){
		    Ion[k].wtmass[bin+nzloopbins] += kernel*ion_weight*cp->mass ;
		    Ion[k].age[bin+nzloopbins] += kernel*ion_weight*cp->mass*cp->ageaway ;
		    Ion[k].dgal[bin+nzloopbins] += kernel*ion_weight*cp->mass*(log10(cp->dgal)) ;
		    Ion[k].mgal[bin+nzloopbins] += kernel*ion_weight*cp->mass*cp->mgal ;
		    Ion[k].vlaunch[bin+nzloopbins] += kernel*ion_weight*cp->mass*(cp->vlaunch) ;
		    Ion[k].nrec[bin+nzloopbins] += kernel*ion_weight*cp->mass*(cp->nrec) ;
		  }
#endif
		} // FOR: 0 < k < nions
#ifdef PHEW
	      } // UP: wind_flag <= 0
	      else{ // IS a PhEW wind !
		nc_bin = kernel * cp->ncloud * (PI * cp->rcloud * cp->rcloud); // kernel has been normalized by PI * h^2
		/* prob0 = exp(-nc_bin); // k=0 term of the Poisson distribution		       */
		prob1 = nc_bin * exp(-nc_bin); // k=1 term of the Poisson distribution
		prob2 = prob1 * nc_bin / 2.0; // k=2 term of the Poisson distribution
#ifdef PHEW_VERBOSE
		fprintf(stderr, "PHEW: %d %g %5.3e %5.1f %4.2f %5.3f %4.2f %g %g %g\n",
			i, radius, kernel,
			cp->ncloud, log10(cp->temp),
			cp->rcloud * unit_Length * aex / KPC,
			cp->hsmooth * unit_Length * aex / KPC,
			nc_bin, prob1, prob2);
#endif		
		// Write prob0, prob1, prob2 to the partzfile ...
		if(get_random_number(cp->idx + 12) < prob1){
		  nc_this_bin = 1;
#ifdef PHEW_VERBOSE		  
		  fprintf(stderr, "GET ONE CLOUD!\n");
#endif		  
		  if(get_random_number(cp->idx + 12) < prob2){
		    nc_this_bin = 2;
#ifdef PHEW_VERBOSE		  
		  fprintf(stderr, "GET ONE MORE CLOUD!\n");
#endif		  
		  }
		  // I want to find out which zbins this particle really belongs to here, after vel shift
		  // I will still use *_interp variables, though they mean different!

		  // **************** Construction Begin ... ****************

		  // hubble_expansion and cosmopar() are already defined before, like normal SPH particles.

		  for(ionid = -1; ionid < nions; ionid ++){
		    if(ionid == -1){
		      I = IonTotal;
		      I.atomwt = Ion[0].atomwt; 
		      I.fraction = 0.0122; 
		      I.Zcolumn = -1; 
		      I.alpha = 0;
		      I.bsys = 1e-10;
		      I.Xsec = Ion[0].Xsec;
		      /* bin_cen = binarysearch(vcloud-vmin,vbin_coord,nvbinsnew); */
		      /* bin_cen = binarysearch(vcloud-vmin+vbin_size[bin_cen]/2,vbin_coord,nvbinsnew); */
		      /* if (bin_cen < 0 || bin_cen >= nvbinsnew) continue; // Moving on to the new ion. */
		    }
		    else{
		      I = Ion[ionid];
		    }
		    // colcloud: Average column density of a single cloudlet
		    Zcol = I.Zcolumn;
		    if(Zcol==-1){
		      if(ionid >= 0) ion_weight = ionfrac[ionid]; // HI, HeII. 
		      else ion_weight = 1.0; // IonTotal
		    }else{
		      if(Zcol<-1){
			ion_weight = ionfrac[ionid]*cp->metals[1]/0.009618*I.fraction*pow(10,I.alpha);
			// Not C, O, Si, Fe
		      }else{
			ion_weight = ionfrac[ionid]*cp->metals[Zcol];
		      }
		    } // IF: Zcol ==/</> -1
		    // note ionfrac[k] is calculated independently for each cp

		    colcloud = cp->rho * 4. / 3. * cp->rcloud * ion_weight;
		    rhocloud = colcloud * cp->rho * XH * unit_Density / (aex * aex * aex);
		    tcloud = colcloud * cp->temp;

		    vcloud = vz * unit_vel;
		    vcloud += vkernel / kernel * hubble * unit_vel;
		    // colcloud is equivalent to cp->mass * kernel * ion_weight
		    vcloud += voffset[bin+nzloopbins];
		    // voffset determines the cloud in velocity space
		    // in tau.c, it is hubble_expansion
		    
#ifdef ZEROVEL
		    b = 0.1; //Negligable temperature braodening in 0-velocity case.
		    vcloud = 0.0;  //no velocity case!
#else	
		    b = I.bsys*sqrt(cp->temp) ;
#endif
		    irepz = 0;
#ifdef SHORTSPEC
		    for(irepz = -1; irepz <= 1; irepz++) {
#endif
		      vbin_min = binarysearch((vcloud + irepz*vstep/1.e5*nvbinsnew - NBSMOOTH*b),vbin_coord,nvbinsnew);
		      vbin_max = binarysearch((vcloud + irepz*vstep/1.e5*nvbinsnew + NBSMOOTH*b),vbin_coord,nvbinsnew);

#ifdef PHEW_VERBOSE
		      if(PHEW_VERBOSE > 5)
			fprintf(stderr, "%g %g | %g [%g %g] %g %g %g %d %d %d %d\n",
				colcloud, kernel*cp->mass*ion_weight, vcloud, vbin_coord[0], vbin_coord[nvbinsnew-1], NBSMOOTH*b, vkernel/kernel*hubble*unit_vel, voffset[bin+nzloopbins], bin, vbin_min, vbin_max, nvbinsnew);
#endif		      

		      if (vbin_min < 0) vbin_min = 0;
		      if (vbin_max >= nvbinsnew) vbin_max = nvbinsnew-1;

		      if((vbin_min == 0 && vbin_max == 0) || (vbin_min >= nvbinsnew-1 && vbin_max >= nvbinsnew-1)){
		      continue;
		      }
		      for(vbin = vbin_min; vbin <= vbin_max; vbin++){ // Should we use another variable instead of "bin"
			// Integration over fac * Rc. The average should be fac = 4/3

			vlower = vbin_coord[vbin] - vcloud - irepz*vstep/1.e5*nvbinsnew ;
			vupper = vbin_coord[vbin+1] - vcloud - irepz*vstep/1.e5*nvbinsnew ;
			vlower /= b ;
			vupper /= b ;
			abs_vlower = fabs(vlower) ;
			abs_vupper = fabs(vupper) ;
	    
			if(vupper*vlower < 0){
			  dvcol = 0.5*(erf(abs_vlower) + erf(vupper)) ;
			}
			else{
			  if(abs_vlower < abs_vupper){
			    dvcol = 0.5*(erf(abs_vupper) - erf(abs_vlower)) ;
			  }
			  else{
			    dvcol = 0.5*(erf(abs_vlower) - erf(abs_vupper)) ;
			  }
			}
			if(colcloud > 0){
			  dvcol *= nc_this_bin;
			  I.vcbins[vbin] += dvcol*colcloud; // The I.fraction is moved to tau.c
			  // ionfrac is already in the colcloud
			  I.rhocbins[vbin] += dvcol*rhocloud;
			  I.tcbins[vbin] += dvcol*tcloud;
			  // norm[bin] = I.vcbins[bin]
			  // These bins are in v-space
			  /* fprintf(stderr, "I[%d]: %d %g %g %g %g %g %g\n", */
			  /* 	  ionid, vbin, I.vcbins[vbin], I.tcbins[vbin], I.rhocbins[vbin], dvcol, vlower, vupper); */
			} // colcloud > 0
		      } // bin_min < bin < bin_max
#ifdef SHORTSPEC
		    } // irepz
#endif
		    if(ionid==-1){
		      IonTotal = I;
		    }else{
		      Ion[ionid] = I;
		    }
		  } // Loop: -1 < ionid < nions
		  // **************** Construction End ... ****************		  
		  
		} // The bin contains ONE cloud
	      } // PhEW Wind
#endif	// PHEW	    
	    } // kernel > 0.0
		    
	      //if(bin==0) printf("MIN r = %9.7e k = % 9.7e m = % 9.7e r = % 9.7e b_max = %5d d = %5d cp->rho = %5.3e z = %7.5f aex = %7.5e\n",radius,kernel,IonTotal.mass[bin+nzloopbins],IonTotal.rho[bin+nzloopbins],bin_max,bin_max-bin_min,rhocgs,bin_redshift[bin],XH*unit_Density/(aex*aex*aex));
	      //if(bin==nzbins-1) printf("MAX r = %9.7e k = % 9.7e m = % 9.7e r = %9.7e b_min = %5d d = %5d cp->rho = %5.3e z = %7.5f aex = %7.5e\n",radius,kernel,IonTotal.mass[bin+nzloopbins],IonTotal.rho[bin+nzloopbins],bin_min,bin_max-bin_min,rhocgs,bin_redshift[bin],XH*unit_Density/(aex*aex*aex));

	    } // FOR: bin_min < bin < bin_max
	  }  // irep
	} // irep
      } // irep
    } // FOR: 0 < i < count

    //fprintf(binfile,"# MIN z = %9.7f coord = %9.7f vel = %7.5e temp = %7.5e rho = %7.5e metals = %7.5e mass = %7.5e\n",bin_redshift[nzloopbins],bin_redshift[nzloopbins],IonTotal.vel[nzloopbins],IonTotal.temp[nzloopbins],IonTotal.rho[nzloopbins],IonTotal.metals[0][nzloopbins],IonTotal.mass[nzloopbins]);
    //fprintf(binfile,"# MAX z = %9.7f coord = %9.7f vel = %7.5e temp = %7.5e rho = %7.5e metals = %7.5e mass = %7.5e\n",bin_redshift[nzloopbins+nzbins-1],bin_redshift[nzloopbins+nzbins-1],IonTotal.vel[nzloopbins+nzbins-1],IonTotal.temp[nzloopbins+nzbins-1],IonTotal.rho[nzloopbins+nzbins-1],IonTotal.metals[0][nzloopbins+nzbins-1],IonTotal.mass[nzloopbins+nzbins-1]);


    for(i = nzloopbins; i < nzloopbins+nzbins; i++){
      if(IonTotal.mass[i] != 0.0){
	IonTotal.vel[i] /= IonTotal.mass[i] ;
	IonTotal.temp[i] /= IonTotal.mass[i] ;
	IonTotal.rho[i] /= IonTotal.mass[i] ;
	for(m=0;m<NMETALS;m++){
	  IonTotal.metals[m][i] /= IonTotal.mass[i] ;
	  if(IonTotal.metals[m][i] < 0) IonTotal.metals[m][i] = 0.0e+00;
	}
#ifdef PHYSSPEC
	IonTotal.sfr[i] /= IonTotal.mass[i] ;
	if(IonTotal.wtmass[i] != 0.0){
	  IonTotal.age[i] /= IonTotal.wtmass[i];
	  IonTotal.dgal[i] /= IonTotal.wtmass[i];
	  IonTotal.mgal[i] /= IonTotal.wtmass[i];
	  IonTotal.vlaunch[i] /= IonTotal.wtmass[i];
	  IonTotal.nrec[i] /= IonTotal.wtmass[i];
	}
#endif
      } // IonTotal.mass[i] != 0
      /* nzloopbins <= i < nzloopbins+nzbins */
      IonTotal.redshift[i] = bin_redshift[i-nzloopbins]; // 0 <= idx < nzbins
      IonTotal.binsize[i] = bin_size[i-nzloopbins];
      IonTotal.bincoord[i] = bin_coord[i-nzloopbins]+hold_coord;
      IonExtra.x[i] = bin_x[i-nzloopbins];
      IonExtra.y[i] = bin_y[i-nzloopbins];
      IonExtra.z[i] = bin_z[i-nzloopbins];
      //#ifdef PARTIONFRAC
      for( k=0; k<nions; k++ ) {
	if(Ion[k].mass[i] > 0.0){
	  Ion[k].vel[i] /= Ion[k].mass[i] ;
	  Ion[k].temp[i] /= Ion[k].mass[i] ;
	  Ion[k].rho[i] /= Ion[k].mass[i] ;
	  for(m=0;m<NMETALS;m++){
	    Ion[k].metals[m][i] /= Ion[k].mass[i] ;
	    if(Ion[k].metals[m][i] < 0) Ion[k].metals[m][i] = 0.0e+00;
	  }
#ifdef PHYSSPEC
	  Ion[k].sfr[i] /= Ion[k].mass[i] ;
	  if(Ion[k].wtmass[i] != 0.0){
	    Ion[k].age[i] /= Ion[k].wtmass[i];
	    Ion[k].dgal[i] /= Ion[k].wtmass[i];
	    Ion[k].mgal[i] /= Ion[k].wtmass[i];
	    Ion[k].vlaunch[i] /= Ion[k].wtmass[i];
	    Ion[k].nrec[i] /= Ion[k].wtmass[i];
	  }
#endif
	} // Ion[k].mass[i] > 0
      } // 0 <= k < nions
      //#endif

#ifndef PIPELINE
      if(i==nzloopbins)fprintf(binfile,"# MIN z = %9.7f coord = %9.7f vel = %7.5e temp = %7.5e rho = %7.5e metals = %7.5e mass = %7.5e\n",bin_redshift[i-nzloopbins],bin_redshift[i-nzloopbins],IonTotal.vel[i],IonTotal.temp[i],IonTotal.rho[i],IonTotal.metals[0][i],IonTotal.mass[i]);
      if(i==nzloopbins+nzbins-1)fprintf(binfile,"# MAX z = %9.7f coord = %9.7f vel = %7.5e temp = %7.5e rho = %7.5e metals = %7.5e mass = %7.5e\n",bin_redshift[i-nzloopbins],bin_redshift[i-nzloopbins],IonTotal.vel[i],IonTotal.temp[i],IonTotal.rho[i],IonTotal.metals[0][i],IonTotal.mass[i]);
      /* 	  fprintf(binfile,"%5d  %9.7f %5.3e % 8.3f % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e % 9.7e %5.3e",i,IonTotal.redshift[i],IonTotal.mass[i], IonTotal.vel[i],IonTotal.temp[i],IonTotal.rho[i],IonTotal.metals[0][i],IonTotal.metals[1][i],IonTotal.metals[2][i],IonTotal.metals[3][i],IonTotal.bincoord[i],IonTotal.binsize[i]); */
      fprintf(binfile,"%5d  %9.7f %5.3e % 8.3f % 5.3e % 5.3e %g %g %g %g % 9.7e %5.3e",i,IonTotal.redshift[i],IonTotal.mass[i], IonTotal.vel[i],IonTotal.temp[i],IonTotal.rho[i],IonTotal.metals[0][i],IonTotal.metals[1][i],IonTotal.metals[2][i],IonTotal.metals[3][i],IonTotal.bincoord[i],IonTotal.binsize[i]);
      if(IonTotal.metals[3][i] > 10){
	printf("BAD METAL!! %d %g\n", i, IonTotal.metals[3][i]);
	exit(-1);
      }
      fprintf(binfile," % 7.5e % 7.5e % 7.5e",IonExtra.x[i],IonExtra.y[i],IonExtra.z[i]);
#ifdef DO6IONS 
      fprintf(binfile," ",Ion[4].mass[i],Ion[4].vel[i],Ion[4].rho[i],Ion[4].temp[i],Ion[4].metals[0][i]);
#endif
      fprintf(binfile," %5.3e %8.3f % 5.3e %5.3e %5.3e ",Ion[6].mass[i],Ion[6].vel[i],Ion[6].rho[i],Ion[6].temp[i],Ion[6].metals[0][i]);
#ifdef PHYSSPEC
      fprintf(binfile," %5.2f %5.2f %5.2f %6.1f %5.3f",IonTotal.mgal[i],IonTotal.age[i],IonTotal.dgal[i],IonTotal.vlaunch[i],IonTotal.nrec[i]);
#endif
      fprintf(binfile,"\n");
#endif

    } // nzloopbins <= i < nzloopbins + nzbins
#ifndef PIPELINE
    fprintf(binfile,"# edge redshift = %9.7f hold_coord = %10.7f\n",IonTotal.redshift[nzloopbins+nzbins-1],hold_coord);
#endif

    hold_coord += bin_coord[nzbins-1]+bin_size[nzbins-1];

#ifndef PIPELINE
    fclose(binfile);
    fclose(partfile);
#endif

    //free(bin_size);
    //free(bin_coord);
    //free(bin_redshift);
#ifndef SHORTSPEC	
    free(spec_particles); 
    //cosmopar(CosmicTime(zout)); /* Why is this necessary? BDO 10/28/09
#endif
#ifdef PHEW    
    free(voffset);
    free(vbin_size);
    free(vbin_coord);
#endif    
	
    return 0;
}


/* this function generates a random number in a Gaussian distribution 
   using the Box-Muller method */

#ifdef PAINTAVERAGEMETALS		

 float findrandgauss()
 {
  
   float x1, x2, w, y1, y2;
  
   do {
     x1 = 2.0 * drand48() - 1.0;
     x2 = 2.0 * drand48() - 1.0;
     w = x1 * x1 + x2 * x2;
   } while ( w >= 1.0 );
  
   w = sqrt( (-2.0 * log( w ) ) / w );
   y1 = x1 * w;
   y2 = x2 * w;

   return(y1);
 }
#endif

#ifdef INTKERNELNHLIMIT

 int InitKernIntTable()
 {
   int i;
   float xw,dx,sum,kint,kern;
   float Kernel();
  
   sum = kint = 0;
   dx = NSRCHRAD/NINTERP;
   for(i=NINTERP-1; i>=0; i--) {
     xw = i*dx;
     if(xw<=0.5){
       kern = 1-6*xw*xw+6*xw*xw*xw;
     }else{
       kern = 2*(1-xw)*(1-xw)*(1-xw);
     }
     kern *= 8./M_PI;
     sum += kern*4*M_PI*xw*xw*dx;
     kint += kern*dx;
     KernIntTable[i][0] = sum;
     KernIntTable[i][1] = kint;
     if( i%1000 == 0 ) fprintf(stdout,"KERNTABLE BUILD: %d %g %g %g\n",i,xw,sum,kint);
   }
   return 0;
 }

#endif
