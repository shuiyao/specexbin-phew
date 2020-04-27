#ifdef TIPSYFORMAT
#include "tipsydefs.h"
#else
#ifdef TIPSYN2FORMAT
#include "tipsydefs_n2.h"
#else
#include "tipsydefs_n.h" 
#endif
#endif
#include "specexbindefs.h"
#include "extern.h"

#define NRANSI

#define clight		2.99792458e10

int RotateCoords(float, float, float, float, float);
int InverseRotateCoords(float, float, float, float, float);
#ifdef CALCQSOBKGD
int poisson_QSO(float, long *);
float gaussrand(long *);
#endif 

extern double unit_Velocity,unit_Density,unit_Time;

int GetSpecParticles(char *finname){

#ifdef PHYSSPEC
  FILE *metaltrackfile;
  char metaltrackname[100];
#else
  FILE *binfile, *tabfile;
  char snapname[100], binname[100], tabname[30];
#ifdef TIPSYN2FORMAT
  FILE *auxfile;
  char auxname[100];
#endif
#endif
  
  int FindLOS(), CheckParticleLOS();
  int i; 
  Real hold,metalhold;
  int int_hold;
  float redshift_tab;
  float CosmicTime();
  int  cosmopar();

#ifdef CALCQSOBCKD
  int PlaceQSOs();
#endif

#ifdef PHYSSPEC
  sprintf(metaltrackname,"%s",finname);
  if( (metaltrackfile = fopen(metaltrackname,"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",metaltrackname);
    return 0;
  }  
  //sim_id = strtok(finnametmp,".");
  //fprintf(stderr,"SIM_ID= %s\n",sim_id);
#else
  sprintf(tabname,"%s.tab",finname);
  if( (tabfile = fopen(tabname,"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",tabname);
    return 0;
  }  
  while(!feof(tabfile)){
    fscanf(tabfile,"%s %f %f %f %s", snapname, &redshift_tab, &redshift_tab_begin, &redshift_tab_end, sim_id);
    if(redshift >= redshift_tab_begin && redshift <= redshift_tab_end){
      break;
    }
  }
  fclose(tabfile);
  fprintf(stderr,"ATTENTION: Opening %s corresponding to redshift %5.3f\n",snapname, redshift_tab);
  sprintf(binname,"%s.bin",snapname);
  if( (binfile = fopen(binname,"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",binname);
    return 0;
  }
#ifdef TIPSYN2FORMAT 
  sprintf(auxname,"%s.aux",snapname);
  if( (binfile = fopen(binname,"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",binname);
    return 0;
  }
  if( (auxfile = fopen(auxname,"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",auxname);
    return 0;
  }
#endif
#endif

  nzbins = FindLOS();
  
#ifdef CALCQSOBKGD
  PlaceQSOs();
#endif

#ifdef PHYSSPEC
  fprintf(stderr,"BEGIN READ METAL FILE\n");
  while(!feof(metaltrackfile)){
    if(feof(metaltrackfile))break;
    fread(hp,sizeof(struct spec_particle),1,metaltrackfile);
#else
    fread(&header,sizeof(struct dump),1,binfile);
    /*fread(&header.time, sizeof(header.time), 1, binfile);
    fread(&header.nbodies, sizeof(header.nbodies), 1, binfile);
    fread(&header.ndim, sizeof(header.ndim), 1, binfile);
    fread(&header.nsph, sizeof(header.nsph), 1, binfile);
    fread(&header.ndark, sizeof(header.ndark), 1, binfile);
    fread(&header.nstar, sizeof(header.nstar), 1, binfile);*/
    //#ifndef TIPSYFORMAT /* Goddamn tipsy header format, has 28 bytes and I'm trying to make it 32 bytes.*/
    //fread(&header.pad, sizeof(header.pad), 1, binfile);
    //#endif
  for(i=0; i < header.nsph ; i++) {
#ifdef TIPSYNFORMAT /* The format I used for my thesis that is now obslete (BDO) */
    fread(&hp->mass,sizeof(Real),1,binfile);
    fread(hp->pos,sizeof(Real),MAXDIM,binfile);
    fread(hp->vel,sizeof(Real),MAXDIM,binfile);
    fread(&hp->rho,sizeof(Real),1,binfile);
    fread(&hp->temp,sizeof(Real),1,binfile);
    fread(&hp->hsmooth,sizeof(Real),1,binfile);
    fread(hp->metals,sizeof(Real),NMETALS,binfile);
    fread(&hold,sizeof(Real),1,binfile);
    fread(&hold,sizeof(Real),1,binfile);
    fread(&hold,sizeof(Real),1,binfile);
#else
    fread(&hp->mass,sizeof(Real),1,binfile);
    fread(&hp->pos,sizeof(Real),MAXDIM,binfile);
    fread(&hp->vel,sizeof(Real),MAXDIM,binfile);
    fread(&hp->rho,sizeof(Real),1,binfile);
    fread(&hp->temp,sizeof(Real),1,binfile);
    fread(&hp->hsmooth,sizeof(Real),1,binfile);
    fread(&metalhold,sizeof(Real),1,binfile);    
    fread(&hold,sizeof(Real),1,binfile);
    /* Specexbin assumes 4 metals now: carbon, oxygen, silicon, iron.*/
    hp->metals[0] = metalhold*0.00213/0.0122;
    hp->metals[1] = metalhold*0.00541/0.0122;
    hp->metals[2] = metalhold*0.00067/0.0122;
    hp->metals[3] = metalhold*0.00117/0.0122;
#ifdef TIPSYN2FORMAT
    fread(&hp->metals,sizeof(Real),NMETALS,auxfile);
    fread(&hp->sfr,sizeof(Real),1,auxfile);
    fread(&hold,sizeof(Real),1,auxfile);
    fread(&hold,sizeof(Real),1,auxfile);
    fread(&hold,sizeof(Real),1,auxfile);
    fread(&hold,sizeof(Real),1,auxfile);
    fread(&int_hold,sizeof(int),1,auxfile);
#endif // TIPSYN2FORMAT
#endif // TIPSYNFORMAT
#endif // PHYSSPEC

    if(direction==1){
      hold = hp->pos[2];
      hp->pos[2] = hp->pos[1];
      hp->pos[1] = hp->pos[0];
      hp->pos[0] = hold;
      
      hold = hp->vel[2];
      hp->vel[2] = hp->vel[1];
      hp->vel[1] = hp->vel[0];
      hp->vel[0] = hold;
    }

    if(direction==0){
      hold = hp->pos[0];
      hp->pos[0] = hp->pos[1];
      hp->pos[1] = hp->pos[2];
      hp->pos[2] = hold;
      
      hold = hp->vel[0];
      hp->vel[0] = hp->vel[1];
      hp->vel[1] = hp->vel[2];
      hp->vel[2] = hold;
    }

    CheckParticleLOS();
    
  }
  
#ifdef PHYSSPEC
  fclose(metaltrackfile);
#else
  fclose(binfile);
#ifdef TIPSYN2FORMAT
  fclose(auxfile);
#endif
#endif
  return(count);
}


int FindLOS()
{
  int i,vi_hold;
  float CosmicTime();
  int  cosmopar();

  redshift_hold = redshift;

  cosmopar(CosmicTime(redshift));
  zout = redshift;
  zstep = ZRES/(BOXSIZE*hubble*unit_Velocity/clight);
  
  //fprintf(stderr,"Spectrum Particles Output: at z=%g, redshift param = %g\n",redshift,BOXSIZE*hubble*unit_Velocity/clight);

  xbegin = xspec;
  ybegin = yspec;
  zbegin = zspec;

  RotateCoords(xbegin, ybegin, zbegin, theta, phi);

  xrotline = xrot;
  yrotline = yrot;
  zbeginline = zrot+zstep;

  xorig = xbegin;
  yorig = ybegin;
  zorig = zbegin;

    /* find xend, yend, and zend */
  xbreak = 0;
  ybreak = 0;
  zbreak = 0;
  i = 0;

  bin_size = malloc(sizeof(double));
  bin_coord = malloc(sizeof(double));
  bin_redshift = malloc(sizeof(double));

  /* use this loop to run line-of-site through box and make bins in 
   box coordinates (bin_coord) and redshift coordinates (bin_redshift) */
  vi_hold = vi;
  while(xbreak == 0 && ybreak == 0 && zbreak == 0){
    zout -= ZRES;
    cosmopar(CosmicTime(zout));
    zstep = ZRES/(BOXSIZE*hubble*unit_Velocity/clight);
    zrot += zstep;
    InverseRotateCoords(xrot, yrot, zrot, theta, phi);
    i++;
    vi = floor((i+nzloopbins)*ZRES/VRES);
    if(vi > vi_hold){
      bin_x = realloc(bin_x,(vi+1)*sizeof(double));
      bin_y = realloc(bin_y,(vi+1)*sizeof(double));
      bin_z = realloc(bin_z,(vi+1)*sizeof(double));
      bin_x[vi] = xorig;
      bin_y[vi] = yorig;
      bin_z[vi] = zorig;
      vi_hold = vi;
      //if(i>1)printf("%d %d % f % f % f %f\n",nzloopbins+i,vi,bin_x[vi],bin_y[vi],bin_z[vi],bin_redshift[i-2]);
    }
    bin_size = realloc(bin_size,i*sizeof(double));
    bin_coord = realloc(bin_coord,i*sizeof(double));
    bin_redshift = realloc(bin_redshift,i*sizeof(double));
    bin_size[i-1] = zstep;
    bin_coord[i-1] = zrot-zbeginline;
    bin_redshift[i-1] = zout;
    if(xorig < -HALFBOX || xorig > HALFBOX) xbreak = xorig;
    if(yorig < -HALFBOX || yorig > HALFBOX) ybreak = yorig;
    if(zorig > HALFBOX) zbreak = zorig;
    //printf("zout = %5.3e zorig = %9.7f i = %5d,binsize = %5.3e, bincoord = %5.3e bin_redshift = %9.7f bin_hubble = %5.3e\n",zout,zorig,i-1,bin_size[i-1],bin_coord[i-1],bin_redshift[i-1], bin_hubble[i-1]);
  }        

  zrot -= zstep;
  i--;
  redshift += ZRES;
  zout = redshift;
  InverseRotateCoords(xrot, yrot, zrot, theta, phi);
  xspec = xorig;
  yspec = yorig;
  zspec = zorig;

  zlength = sqrt(pow(xbegin-xspec,2)+pow(ybegin-yspec,2)+pow(zbegin-zspec,2));

  fprintf(stderr,"zlength = %5.3e\n",zlength);
  fprintf(stderr,"xbegin = %9.7e, ybegin = %9.7e, zbegin = %9.7e\n",xbegin,ybegin,zbegin);
  fprintf(stderr,"xend = %9.7e, yend = %9.7e, zend = %9.7e\n",xspec,yspec,zspec);
  fprintf(stderr,"xbreak = %9.7e, ybreak = %9.7e, zbreak = %9.7e\n",xbreak,ybreak,zbreak);
  fprintf(stderr,"redshift = %5.3e redshift covered = %5.3e\n",redshift,redshift_hold-redshift);
  redshift_hold = redshift;
  count = 0;
  return(i);
} 


int CheckParticleLOS()
{
  double dx,dy,dr2;
  double irep[NDIM];
  int int_hold;

  //if( hp->mass == 0. ) continue;
  RotateCoords(hp->pos[0], hp->pos[1], hp->pos[2], theta, phi);
  
  dx = fabs(xrot-xrotline);
  dy = fabs(yrot-yrotline);

  for(irep[0] = -1; irep[0] <= 1; irep[0]++) {
    for(irep[1] = -1; irep[1] <= 1; irep[1]++) {
      for(irep[2] = -1; irep[2] <= 1; irep[2]++) {
	if(fabs(hp->pos[0]+irep[0]*BOXSIZE+2*hp->hsmooth)<HALFBOX || fabs(hp->pos[0]+irep[0]*BOXSIZE-2*hp->hsmooth)<HALFBOX){
	  if(fabs(hp->pos[1]+irep[1]*BOXSIZE+2*hp->hsmooth)<HALFBOX || fabs(hp->pos[1]+irep[1]*BOXSIZE-2*hp->hsmooth)<HALFBOX){
	    if(fabs(hp->pos[2]+irep[2]*BOXSIZE+2*hp->hsmooth)<HALFBOX || fabs(hp->pos[2]+irep[2]*BOXSIZE-2*hp->hsmooth)<HALFBOX){
	      RotateCoords(hp->pos[0]+irep[0]*BOXSIZE, hp->pos[1]+irep[1]*BOXSIZE, hp->pos[2]+irep[2]*BOXSIZE, theta, phi);
	      if(fabs(xrot-xrotline)<dx) dx = fabs(xrot-xrotline);
	      if(fabs(yrot-yrotline)<dy) dy = fabs(yrot-yrotline);
	    }
	  }
	}
      }
    }
  }
  
  dr2 = dx*dx+dy*dy;
  if( dr2 < NSRCHRAD*hp->hsmooth*NSRCHRAD*hp->hsmooth ){
#ifdef PHYSSPEC
    if(count==0){
      spec_particles = (struct spec_particle *) malloc(sizeof(struct spec_particle));
    }else{
      spec_particles = realloc(spec_particles,(count+2)*sizeof(struct spec_particle));
    }
#else    
    if(count==0){
      spec_particles = (struct spec_particle *) malloc(sizeof(struct spec_particle));
    }else{
      spec_particles = realloc(spec_particles,(count+2)*sizeof(struct spec_particle));
    }
#endif
    spec_particles[count] = hp[0];
    count++;
  }

  return(0);
}


#ifdef CALCQSOBKGD

int PlaceQSOs()
{

  float CosmicTime();
  int  cosmopar();
  int M,i,j,n_QSO;
  float z,xi,mu,dl,vol,phi,prob,R_0,R_M;
  float L_2500, L_228;
  time_t t1;
  long idum;
  float sum_vol=0, sum_L228=0, sum_L1216=0;
  float L, Lstar, alphas;
  //float ran1(long *idum);

  (void) time(&t1);
  idum = -time(&t1);

  for(i=0;i<nzbins;i++){
    z=bin_redshift[i];
    cosmopar(CosmicTime(z));
    xi = log10((1+z)/(1+2.45));   
    mu = -26+1.43*xi+36.63*xi*xi+34.39*xi*xi*xi;
    dl = ZRES/(BOXSIZE*hubble*unit_Velocity/clight)*boxsize;
    R_0 = 30*pow((1+z)/4,-3);
    R_M = R_0*5.0;
    vol = dl*PI*pow(R_M,2);
    sum_vol += vol;
    for(M=22;M<=30;M++){
      if(z>2.4) phi = pow(10,-5.70+(-M*1.-mu)*(0.83-0.11*(z-2.45))); //QSO Luminsoity frunction from Richards et al.2006 
      else phi = pow(10,-5.70+(-M*1.-mu)*0.84);
      //alphas = 0.3;
      //Lstar = pow(10,10.67)*pow((1+z),(alphas-1))*(exp(2.58*z)*(1+exp(3.16*1.9))/(exp(3.16*z)+exp(3.16*1.9))); // Madau et al. 1999, B-band (4860 Angstroms)
      //L = pow(10,-0.4*(-M-5.41));
      //phi = pow(10,-6.01)*L/Lstar/(pow(L/Lstar,3.41)+pow(L/Lstar,1.58));
      //printf("phi = % 5.3e L = % 5.3e L* = % 5.3e M = -%2d\n",phi, L, Lstar, M);

      //R_M = R_0*5.0*pow(2.5,0.5*(M-22));
      //R_M = 5.0*R_0;
      //vol = dl*PI*R_M*R_M;
      prob = phi*vol;
      n_QSO = poisson_QSO(prob, &idum);
      for(j=0;j<n_QSO;j++){
	QSOs.z[nQSO] = z;
	QSOs.r[nQSO] = sqrt(ran1(&idum))*R_M;
	QSOs.M[nQSO] = ran1(&idum)-0.5+M;
	if(QSOs.M[nQSO]>=22.0){
	  QSOs.R_M[nQSO] = R_M;
	  L_2500 = pow(10,-0.4*(-QSOs.M[nQSO]+48.60+2.5*log10(3.0)))*(4*PI*pow(3.08e+19,2)); //Richards et al. 2006 eqn 4
	  L_228 = L_2500*pow(1216.0/2500.0,0.5+0.3*gaussrand(&idum))*pow(228.0/1216.0,1.57+0.17*gaussrand(&idum));
	  QSOs.L228[nQSO] = L_228;
	  sum_L228 += L_228;
	  sum_L1216 += L_2500*pow(1216.0/2500.0,0.5+0.3*gaussrand(&idum));
	  nQSO++;
	  QSOs.z = realloc(QSOs.z,(nQSO+1)*sizeof(float));
	  QSOs.r = realloc(QSOs.r,(nQSO+1)*sizeof(float));
	  QSOs.M = realloc(QSOs.M,(nQSO+1)*sizeof(float));
	  QSOs.R_M = realloc(QSOs.R_M,(nQSO+1)*sizeof(float));
	  QSOs.L228 = realloc(QSOs.L228,(nQSO+1)*sizeof(float));
	}
      }
      if(M>26 && n_QSO > 0)fprintf(stderr,"z = %5.3f phi = %5.3e vol = %5.3e prob = %5.3e n_QSO = %2d M = %d i = %6d dl = %5.3e r_QSO = %5.1f, R_M = %5.1f M_QSO = %5.2f L_228 = % 5.3e\n",QSOs.z[nQSO-1],phi,vol,prob,n_QSO,M,i,dl,QSOs.r[nQSO-1],R_M,QSOs.M[nQSO-1],QSOs.L228[nQSO-1]);
    }	
  }    
  fprintf(stderr,"vol = %5.3e L228 = %5.3e L1216 = %5.3e eps_Q(228) = %5.3e eps_Q(1216) = %5.3e nQSO = %8d\n",sum_vol, sum_L228, sum_L1216, sum_L228/sum_vol, sum_L1216/sum_vol, nQSO);
  return(0);
}

int poisson_QSO(float f_QSO, long *idum)
{
  int nQSO;

  //fprintf(stderr,"f_QSO = %f idnum = %ld\n",f_QSO, *idum);
  nQSO = poidev(f_QSO, idum);
  //fprintf(stderr,"nQSO = %d\n",nQSO);
  return(nQSO);
}

float gaussrand(long *idum)
{
  
  float x1, x2, w, y1, y2;
  
  do {
    x1 = 2.0 * ran1(idum) - 1.0;
    x2 = 2.0 * ran1(idum) - 1.0;
    w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );
  
  w = sqrt( (-2.0 * log( w ) ) / w );
  y1 = x1 * w;
  y2 = x2 * w;

  return(y1);
}

#endif

