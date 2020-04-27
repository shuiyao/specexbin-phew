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
#ifdef HDF5FORMAT 
#include "hdf5.h"
#endif // HDF5FORMAT
#endif // TIPSYFORMAT
#include "proto.h"
#include "specexbindefs.h"
#include "extern.h"


//#define NSUBBINMIN	4 // changed from 4 to 20 on Feb 28,2011 by Amanda Ford, trying to make smaller pixels
#define NSUBBINMIN	4
#define NSUBBINMAX	1000000
#define SIGN(x) ((x)<0.0 ? -1.0 : 1.0)
#define clight		2.99792458e10
#define ZRES_ION_TABLE  0.05
extern double unit_Length,unit_Time,unit_Velocity,unit_Density,unit_Mass;

double vmin,vmax,vbin_size;
extern float H0;

int Tau()
/* Calculate optical depth for given ion along line of sight */
{
  int ionid;
  int i,j,k;
  int Zcol, imet;
  /* float z; */
  /* int l; */
  float b;
  double mass_interp, v_interp;
  double irepz;
  float t_interp,rho_interp,Z_interp;
  double voffset;
  int bin,bin_min,bin_max,bin_cen;
  float vlower,vupper;
  float abs_vlower,abs_vupper;
  float dvcol;
  double unit_col;
  float *norm;
  int *norm_field;
  double *vbin_size, *vbin_coord, *vbin_zsize, *vbin_zcoord;
  double vcoord, vstep, zcoord;
  float CosmicTime();
  float hubble_expansion;
  int cosmopar();
#ifdef SMOOTHSPH
  int load_fraction_tables();
  float IonFrac(), ionfrac_interp, ionize_redshift;
#endif
#ifdef NONEQUIL
  int load_noneq_fraction_tables();
  float Noneq_Frac();
#endif
#ifdef BTURB
  double bturb;
  double epsilon, lpart;
#endif
#ifdef VARGALBKGD
  float mgal_100kpc_interp;
#endif
  ionStruct I;
  double x_interp, y_interp, z_interp, weight_sub;
  int nsubbinvar;
  
  /* we are out of the loop so reassign nzbins and nvbins */
  nzbins = nzloopbins;
  nvbins = nvloopbins;

  fprintf(stderr,"Begin tau()\n");

#ifdef SHORTSPEC  
  cosmopar(CosmicTime(redshift_center));
#else
  cosmopar(CosmicTime(IonTotal.redshift[0]));
#endif

  redshift_track = IonTotal.redshift[0];
  vbin_size = malloc(sizeof(double));
  vbin_coord = malloc(sizeof(double));
  vbin_zsize = malloc(sizeof(double));
  vbin_zcoord = malloc(sizeof(double));
  norm = malloc(nvbins*sizeof(float));
  norm_field = malloc(nvbins*sizeof(int));

  i = 0;
  vcoord = 0;
  zcoord = 0;
  fprintf(stderr,"TAUBEGIN: redshift_track= %g  IonTotal.redshift[nzbins-1]= %g nzbins= %d\n",redshift_track,IonTotal.redshift[nzbins-1],nzbins);
  while(redshift_track >= IonTotal.redshift[nzbins-1]){
    redshift_track -= vres;
#ifndef SHORTSPEC  
    cosmopar(CosmicTime(redshift_track));
#endif
    zstep = vres/(BOXSIZE*hubble*unit_Velocity/clight);
    zcoord += zstep;
    vstep = zstep*aex*hubble*unit_Velocity;
    vcoord += vstep;
    i++;
    vbin_zsize = realloc(vbin_zsize,i*sizeof(double));
    vbin_zcoord = realloc(vbin_zcoord,i*sizeof(double));	
    vbin_size = realloc(vbin_size,i*sizeof(double));
    vbin_coord = realloc(vbin_coord,i*sizeof(double));
    
    vbin_size[i-1] = vstep/1.e5;
    vbin_coord[i-1] = vcoord/1.e5;
    vbin_zsize[i-1] = zstep;
    vbin_zcoord[i-1] = zcoord;
  }
  fprintf(stderr, "TAUBEGIN: i=%d nvbins=%d\n",
	  i-1, nvbins);

  voffset = 0;
  vmin = 0.0;

  for(ionid = -1; ionid < nions; ionid++){
    
    if( ionid==-1 ) fprintf(stderr,"Computing tau for ion ");
    fprintf(stderr,"%d ",ionid);
    if( ionid==nions-1 ) fprintf(stderr,"\n");
    fflush(stderr);    
    if(ionid==-1){ /* Outputing in physical space */
      I = IonTotal;
      I.atomwt = Ion[0].atomwt; 
      I.fraction = 0.0122; 
      I.Zcolumn = -1; 
      I.alpha = 0;
      I.bsys = 1e-10;
      I.Xsec = Ion[0].Xsec;
    }else{
      I = Ion[ionid];
    }
    for( i=0; i<nvbins; i++ ){
      I.vbins[i] = I.tbins[i] = I.rhobins[i] = I.Zbins[i] = norm[i] = 0.0;
    }
    Zcol = I.Zcolumn;
    hubble_expansion = 0;
    
    redshift = IonTotal.redshift[0];
#ifdef SMOOTHSPH
#ifdef NONEQUIL
    load_noneq_fraction_tables();
#endif
    load_fraction_tables();
    ionize_redshift = redshift - ZRES_ION_TABLE;
#endif    

    for(i = 0; i < nzbins; i++){
#ifndef SHORTSPEC
      cosmopar(CosmicTime(IonTotal.redshift[i]));
#endif
      //vi = floor(i*ZRES/VRES);
      //IonExtra.redshift[vi] = IonTotal.redshift[i];

      hubble_expansion += hubble*aex*unit_Velocity/1.e5*IonTotal.binsize[i];

      nsubbinvar = NSUBBINMIN;
      if(i>0 && i<nzbins-1){
	nsubbinvar = abs((int)((I.vel[i+1]-I.vel[i-1])/((IonTotal.redshift[i-1]-IonTotal.redshift[i+1])*clight/1e+05/(1+IonTotal.redshift[i]))));	
      }
      //if(i%100==0 && ionid==1)fprintf(stderr,"nsubbinvar= %3d %5.3e %5.3e",nsubbinvar,(I.vel[i+1]-I.vel[i-1]),(IonTotal.redshift[i-1]-IonTotal.redshift[i+1])*clight/1e+05/(1+IonTotal.redshift[i]));
      if(nsubbinvar<NSUBBINMIN) nsubbinvar=NSUBBINMIN;
      if(nsubbinvar>NSUBBINMAX) nsubbinvar=NSUBBINMAX;     
      //if(i%100==0 && ionid==1)fprintf(stderr," %3d\n",nsubbinvar);


      //if(i!=0)hubble_expansion += hubble*aex*unit_Velocity/1.e5*vbin_zsize[vi]*ZRES/VRES;
      if(IonTotal.mass[i] == 0.0)
	continue ;
      for(j = 0; j < nsubbinvar; j++){
	/* z = (IonTotal.bincoord[i]+((double)(j))/((double)(nsubbinvar))*IonTotal.binsize[i]); */
	if(2*j < nsubbinvar){
	  k = i - 1 ;
	}
	else {
	  k = i ;
	}
	if( k < 0 ) k=0;
	if( k > nzbins-2 ) k = nzbins-2;

#ifdef SMOOTHSPH /* We are going to smooth SPH particles and not treat them individually */
	I = IonTotal;
#endif

	//if(ionid==-1){
	//v_interp = 0;
	//}else{

	weight_sub = (i + ((double)(j))/((double)(nsubbinvar)) - (k + 0.5));
	if( I.mass[k+1] < I.mass[k] ) weight_sub *= 2*I.mass[k+1]/(I.mass[k+1] + I.mass[k]);
	else weight_sub = 1.-2*weight_sub*I.mass[k]/(I.mass[k+1] + I.mass[k]);
	if(ionid==-1){
	  x_interp = (IonExtra.x[k+1] - IonExtra.x[k])*weight_sub + IonExtra.x[k] ;
	  y_interp = (IonExtra.y[k+1] - IonExtra.y[k])*weight_sub +  IonExtra.y[k] ;
	  z_interp = (IonExtra.z[k+1] - IonExtra.z[k])*weight_sub + IonExtra.z[k] ;
	  v_interp = hubble_expansion + hubble*aex*unit_Velocity/1.e5*IonTotal.binsize[k]*((double)(j))/((double)(nsubbinvar)) + voffset;
	  bin_cen = binarysearch(v_interp-vmin,vbin_coord,nvbins);
	  bin_cen = binarysearch(v_interp-vmin+vbin_size[bin_cen]/2,vbin_coord,nvbins);
	  if (bin_cen < 0 || bin_cen >= nvbins) continue;
	  IonExtra.xbins[bin_cen]=x_interp;
	  IonExtra.ybins[bin_cen]=y_interp;
	  IonExtra.zbins[bin_cen]=z_interp;
	}
	mass_interp = (I.mass[k+1] - I.mass[k])*weight_sub + I.mass[k] ;
	// mass_interp must be the column density
	//if(ionid>=0 && mass_interp==0) continue;
	v_interp = (I.vel[k+1] - I.vel[k])*weight_sub + I.vel[k] ; /* No ionization weighting here: NOT TRUE NOW THERE IS */
	  //}
	t_interp = (I.temp[k+1] - I.temp[k])*weight_sub + I.temp[k] ;
	rho_interp = (I.rho[k+1] - I.rho[k])*weight_sub + I.rho[k] ;
	//if( ionid==2 && (k==5062||k==5061)) fprintf(stderr,"%d %g %g %g %g %g %g\n",k,log10(I.temp[k]),log10(I.temp[k+1]),I.mass[k],I.mass[k+1],weight_sub,log10(t_interp));

	if(Zcol==-1){ // H or He or IonTotal
	    Z_interp = 0;
	    for(imet=0;imet<NMETALS;imet++) Z_interp += (I.metals[imet][k+1] - I.metals[imet][k])*(i + ((double)(j))/((double)(nsubbinvar)) - (k + 0.5)) + I.metals[imet][k]; // Sum to total metals                                            
	    Z_interp *= 1.28;
	}else{
	  if(Zcol<-1){
	    Z_interp = (I.metals[1][k+1] - I.metals[1][k])*(i + ((double)(j))/((double)(nsubbinvar)) - (k + 0.5)) + I.metals[1][k]; 
	    //if(Z_interp>0)fprintf(stderr,"IONCALC %d: %g %g ",ionid,I.metals[3][k],Z_interp);
	    //Z_interp *= I.fraction/0.001267*pow(10,I.alpha); 
	    Z_interp *= I.fraction/0.009618*pow(10,I.alpha); /* 2-11-10 */
	    //if(Z_interp>0)fprintf(stderr," %g %g \n",I.fraction/0.001267*pow(10,I.alpha),Z_interp);
	  }else{
	    Z_interp = (I.metals[Zcol][k+1] - I.metals[Zcol][k])*(i + ((double)(j))/((double)(nsubbinvar)) - (k + 0.5)) + I.metals[Zcol][k]; // Fraction by mass of species
	  }
	}

#ifdef SMOOTHSPEC
	I = Ion[ionid];
#endif

#ifdef SMOOTHSPH /* Now run ionization fractions calculations */
	if(ionize_redshift > IonTotal.redshift[i]){
	  redshift = ionize_redshift;
	  load_fraction_tables();
	  ionize_redshift = redshift - ZRES_ION_TABLE;
	}

#ifdef VARGALBKGD
	if(ionid == -1){
	  ionfrac_interp = 1; /* No ion at 0 now.  */
	}else{
	  if(Z_interp>0 && t_interp>0 && rho_interp>0){
	    ionfrac_interp = IonFrac(t_interp,rho_interp,mgal_100kpc_interp,ionid);
	  }else{
	    ionfrac_interp = 0;
	  }
	}
#else
	if(ionid == -1){
	  ionfrac_interp = 1; /* No ion at 0 now.  */
	}else{
#if defined(NONEQUIL) && defined(DO6IONS) 
	  if((ionid==1 && log10(rho_interp/MHYDR)>-3.20) || (ionid==2 && log10(rho_interp/MHYDR)>-3.30) || (ionid==4 && log10(rho_interp/MHYDR)>-4.10) || (ionid==5 && log10(rho_interp/MHYDR)>-2.35)){
	    ionfrac_interp = Noneq_Frac(t_interp,Z_interp*(0.0189/I.fraction),ionid); /* This table is still in Anders & Grevesse Units: 0.189 */
	    printf("% 5.3f\n",log10(rho_interp/MHYDR));
	  }else{
#endif
	    ionfrac_interp = IonFrac(t_interp,rho_interp,ionid);
#if defined(NONEQUIL) && defined(DO6IONS) 
	  }
#endif
	}
	
#endif
	if(Zcol!=-1){
	  ionfrac_interp *= Z_interp; /* This is new */
	}
#endif // SMOOTHSPH ON

	

#ifdef ZEROVEL
	b = 0.1; //Negligable temperature braodening in 0-velocity case.
	v_interp = 0.0;  //no velocity case!
#else	
	b = I.bsys*sqrt(t_interp) ;
#endif

	v_interp += hubble_expansion + hubble*aex*unit_Velocity/1.e5*IonTotal.binsize[i]*((double)(j))/((double)(nsubbinvar)) + voffset;


#ifdef BTURB /* Add aritificial turbulence */
	
	bturb = 0;
	if(log10(rho_interp/MHYDR)>-5.5){
	  if(log10(rho_interp/MHYDR)<-4.5){
	    bturb = sqrt(1405.2*pow(log10(rho_interp/MHYDR),2)+15674.2*log10(rho_interp/MHYDR)+43609.8);
	  }else{
	    if(log10(rho_interp/MHYDR)<-3.0){
	      bturb = 13.93*log10(rho_interp/MHYDR)+101.8;
	    }else{
	      bturb = 60.01;
	    }
	  }
	}
	

	/*if((IonTotal.metals[1][k+1] - IonTotal.metals[1][k])*(i + ((double)(j))/((double)(nsubbinvar)) - (k + 0.5)) + IonTotal.metals[1][k]<0.0001){
	  bturb = 0;
	}
	if(bturb>0 && ionid>=0 && ionfrac_interp>0){
	  if((ionid==0 || ionid== 4)){
	    printf("redshift= %7.5f b(temp)= %5.1f b(turb)= %5.1f b(all)= %5.1f nh= % 5.3f T= %5.3f Z= % 5.3e ionfrac= %5.3e ionid= %2d\n",redshift,b,bturb,sqrt(b*b + bturb*bturb),log10(rho_interp/MHYDR),log10(t_interp),Z_interp,ionfrac_interp/Z_interp,ionid); 
	  }
	  b = sqrt(b*b + bturb*bturb);
	  }*/ // Taken out on 5/12/10- BDO.  
	b = sqrt(b*b + bturb*bturb);

#endif


	/* l = 0; */

	irepz = 0;
	/* SH161009: Distribute the broadened (thermal (+ turbulent)) elements from nsubbinvar to the nearby nzbins */
	
#ifdef SHORTSPEC
	for(irepz = -1; irepz <= 1; irepz++) {
#endif
	  bin_min = binarysearch((v_interp + irepz*vstep/1.e5*nvbins - NBSMOOTH*b - vmin),vbin_coord,nvbins);
	  bin_max = binarysearch((v_interp + irepz*vstep/1.e5*nvbins + NBSMOOTH*b - vmin),vbin_coord,nvbins);

	  /* if(ionid == 4) */
	  /*   fprintf(stderr, "SPH: %g %g %g | %g %g | %g %g %g %d %d %d\n", */
	  /* 	    mass_interp, t_interp, v_interp, */
	  /* 	    hubble_expansion, voffset, */
	  /* 	    vbin_coord[0], vbin_coord[nvbins-1], NBSMOOTH*b, */
	  /* 	    bin_min, bin_max, nvbins); */

	  if (bin_min < 0) bin_min = 0;
	  if (bin_max >= nvbins) bin_max = nvbins-1;
	  
	  //printf("ionid = %d i = %d vi = %5.3e bin_min = %d bin_max = %d nvbins = %d nzbins = %d\n",ionid,i,z,bin_min,bin_max,nvbins,nzbins);

	  if((bin_min == 0 && bin_max == 0) || (bin_min >= nvbins-1 && bin_max >= nvbins-1)) continue;

	  //if((bin_min == bin_max) && (irepz<-0.5 || irepz>0.5)) continue;

	  for(bin = bin_min; bin <= bin_max; bin++){

	    //if(bin == bin_min){
	    //vlower = -NBSMOOTH*b ; /* this appears to be screwing things up! 7-9-11 */
	    //}
	    //else{
	      vlower = vbin_coord[bin] + vmin - v_interp - irepz*vstep/1.e5*nvbins;
	      //}
	      //if(bin == bin_max){
	      //vupper = NBSMOOTH*b ;  /* along with this! 7-9-11 */
	      //}
	      //else{
	      vupper = vbin_coord[bin+1] + vmin - v_interp - irepz*vstep/1.e5*nvbins;
	      //}
	    vlower /= b ;
	    vupper /= b ;
	    abs_vlower = fabs(vlower) ;
	    abs_vupper = fabs(vupper) ;
	    
	    if(vupper*vlower < 0){
	      dvcol = 1./(double)nsubbinvar*0.5*(erf(abs_vlower) + erf(vupper)) ;
	    }
	    else{
	      if(abs_vlower < abs_vupper){
		dvcol = 1./(double)nsubbinvar*0.5*(erf(abs_vupper) - erf(abs_vlower)) ;
	      }
	      else{
		dvcol = 1./(double)nsubbinvar*0.5*(erf(abs_vlower) - erf(abs_vupper)) ;
	      }
	    }
	    if(mass_interp > 0){
#ifdef SMOOTHSPH
	      if(I.Zcolumn==-1){
		I.vbins[bin] += dvcol*IonTotal.mass[i]*ionfrac_interp*I.fraction;
	      }else{
		I.vbins[bin] += dvcol*IonTotal.mass[i]*ionfrac_interp;
	      }
#else 
	      if(I.Zcolumn==-1){ // H, He, or IonTotal
		// I.zolumn == -1 means IonTotal, fraction = 0.0122
		// Or means H, He ... 
		// Know that the mass_interp is from IonTotal.mass = SUM(cp->mass*kernel), no ion_weight
		// For H, He, we must have I.fraction, because M_H = Mass * I.fraction!!!!!
		I.vbins[bin] += dvcol*mass_interp*I.fraction;
		// In the case of wind, the I.fraction is multiplied at the very end.
		//    because I.vcbins is used to normalize rhocbins, tcbins, Zcbins later
	      }else{
		I.vbins[bin] += dvcol*mass_interp;
	      }
	      // Since tau depends on oscillator strength, which is col3.
	      //    This is done in the normalization later.
#endif
	      
	      I.rhobins[bin] += dvcol*rho_interp*mass_interp;
	      I.tbins[bin] += dvcol*t_interp*mass_interp;
	      I.Zbins[bin] += dvcol*Z_interp*mass_interp;
	      norm[bin] += dvcol*mass_interp;
	      //if(ionid==6 || ionid==1)fprintf(stdout,"WRAPAROUNDVEL: ionid= %2d irepz= % 5.3f bin_min= %5d bin_max= %5d bin= %5d v_interp= % 7.2f irepz*vstep/1.e5*nvbins= % 7.2f b= %7.2f vmin= %4.2f vbin_coord[bin]= %7.2f dvcol= % 5.3e vlower= % 7.2f vupper= % 7.2f b= %7.3f t_interp= %5.3e I.temp[k]= %5.3e I.temp[k+1]= %5.3e weight_sub= %5.3e (I.temp[k+1]-I.temp[k])*weight_sub= %5.3e\n",ionid,irepz,bin_min,bin_max,bin,v_interp,irepz*vstep/1.e5*nvbins,b,vmin,vbin_coord[bin],dvcol,vlower,vupper,I.bsys,t_interp,I.temp[k],I.temp[k+1],weight_sub,(I.temp[k+1] - I.temp[k])*weight_sub);
	      //if(weight_sub>1.0)fprintf(stdout,"WRAPAROUNDVELWRONG: i= %5d j= %5d k= %5d nsubbinvar= %5d mass= %5.3e %5.3e ionid= %2d irepz= % 5.3f bin_min= %5d bin_max= %5d bin= %5d v_interp= % 7.2f irepz*vstep/1.e5*nvbins= % 7.2f b= %7.2f vmin= %4.2f vbin_coord[bin]= %7.2f dvcol= % 5.3e vlower= % 7.2f vupper= % 7.2f b= %7.3f t_interp= %5.3e I.temp[k]= %5.3e I.temp[k+1]= %5.3e weight_sub= %5.3e (I.temp[k+1]-I.temp[k])*weight_sub= %5.3e\n",i,j,k,nsubbinvar,I.mass[k+1],I.mass[k],ionid,irepz,bin_min,bin_max,bin,v_interp,irepz*vstep/1.e5*nvbins,b,vmin,vbin_coord[bin],dvcol,vlower,vupper,I.bsys,t_interp,I.temp[k],I.temp[k+1],weight_sub,(I.temp[k+1] - I.temp[k])*weight_sub);



	      //if(ionid==6  && (irepz<-0.5 || irepz>0.5))fprintf(stdout,"WRAPAROUNDVEL: ionid= %2d irepz= % 5.3f bin_min= %5d bin_max= %5d bin= %5d v_interp= % 7.2f irepz*vstep/1.e5*nvbins= % 7.2f b= %7.2f vmin= %4.2f vbin_coord[bin]= %7.2f dvcol= % 5.3e vlower= % 7.2f vupper= % 7.2f b= %7.3f t_interp= %5.3e\n",ionid,irepz,bin_min,bin_max,bin,v_interp,irepz*vstep/1.e5*nvbins,b,vmin,vbin_coord[bin],dvcol,vlower,vupper,I.bsys,sqrt(t_interp));

	    } // mass_interp > 0	    
	  } // bin_min < bin < bin_max
#ifdef SHORTSPEC
	} // irepz
#endif
      } // 0 < j < nsubbinvar
    } // 0 < i < nzbins

    redshift_track = IonTotal.redshift[0];
    for(i = 0; i < nvbins; i++) {
      redshift_track -= vres;
#ifndef SHORTSPEC
      cosmopar(CosmicTime(redshift_track));
#endif      
      unit_col = I.Xsec/(vbin_size[i]) /(aex*aex*unit_Length*unit_Length) /MHYDR;
      // Xsec contains oscillator strength and lambda!! from col3 of the specions_*.dat
      // Ion[i].Xsec *= 2.648e-2*Ion[i].lambda*1.e-13

      I.vbins[i] *= unit_col;
      I.vbins[i] *= unit_Mass/I.atomwt;
      if(norm[i]>0){
	I.tbins[i] /= norm[i];
	I.rhobins[i] /= norm[i];
	I.Zbins[i] /= norm[i];
      }
    }
    if(ionid==-1){
      IonTotal = I;
    }else{
      Ion[ionid] = I;
    }
    
  }

#if defined(PHEW) || defined(PART_BY_PART)
    redshift_track = IonTotal.redshift[0];
    for(i = 0; i < nvbins; i++) {
      redshift_track -= vres;
#ifndef SHORTSPEC
      cosmopar(CosmicTime(redshift_track));
#endif      
      if(IonTotal.vcbins[i] > 0){
	IonTotal.tcbins[i] /= IonTotal.vcbins[i];
	IonTotal.rhocbins[i] /= IonTotal.vcbins[i];
	IonTotal.Zcbins[i] /= IonTotal.vcbins[i];	
      }
      unit_col = Ion[0].Xsec/(vbin_size[i]) /(aex*aex*unit_Length*unit_Length) /MHYDR;
      // Xsec contains oscillator strength and lambda!! from col3 of the specions_*.dat
      // Ion[i].Xsec *= 2.648e-2*Ion[i].lambda*1.e-13
      IonTotal.vcbins[i] *= 0.0122;
      IonTotal.vcbins[i] *= unit_col;
      IonTotal.vcbins[i] *= unit_Mass/Ion[0].atomwt;
      
      for(k = 0; k < nions; k++){
	if(Ion[k].vcbins[i] > 0){
	  Ion[k].tcbins[i] /= Ion[k].vcbins[i];
	  Ion[k].rhocbins[i] /= Ion[k].vcbins[i];
	  Ion[k].Zcbins[i] /= Ion[k].vcbins[i];	  
	}
	unit_col = Ion[k].Xsec/(vbin_size[i]) /(aex*aex*unit_Length*unit_Length) /MHYDR;
	// Xsec contains oscillator strength and lambda!! from col3 of the specions_*.dat
	// Ion[i].Xsec *= 2.648e-2*Ion[i].lambda*1.e-13

	if(Ion[k].Zcolumn == -1) Ion[k].vcbins[i] *= Ion[k].fraction;
	Ion[k].vcbins[i] *= unit_col;
	Ion[k].vcbins[i] *= unit_Mass/Ion[k].atomwt;
      }
    }  
#endif  // PHEW
  
  free(vbin_size);
  free(vbin_coord);
  free(vbin_zsize);
  free(vbin_zcoord);
  free(norm);
  free(norm_field);

  
  return 0;
}

int binarysearch(double key, double *array, int nbins)
{
  int low = 0;
  int high = nbins-1;
  int middle;

  if(key>=array[high]) return high;
  if(key<array[low]) return low;
  
  while( low < high ){
    middle = floor( low  + high ) / 2;
    if(key >= array[middle] && key < array[middle+1]) 
      return middle; 
    else if( key < array[middle] )
      high = middle;		//search low end of array
    else
      low = middle + 1;		//search high end of array    
  }
  return -1;		//search key not found
}

