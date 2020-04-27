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
#include "specexbindefs.h"
#include "proto.h"
#include "extern.h"


extern double unit_Density;

static float nhl[NHPTS],tl[TPTS];

#ifndef VARGALBKGD
static float fractab[MAXIONS][NHPTS][TPTS];
float IonFrac(T,rho_H,ionid)
float T,rho_H;
int ionid;
{
	float frac;
	int calc_fractions();

	rho_H /= flux_fac; /* This is new */
	calc_fractions(T,rho_H,ionid,&frac);
	if( frac < 1.e-20 ) frac = 0.0;

#ifdef NEUTRALGAS 
#ifdef DO6IONS
	if(ionid==0 || ionid==3){
	  frac = 1.0;
	}else{
	  frac = 0.0;
	}
#endif
#ifdef DOHIZIONS 
	if(ionid==0 || ionid==2 || ionid==3){
	  frac = 1.0;
	}else{
	  frac = 0.0;
	}
#endif
#endif

	return frac;
}
#else
static float fractab[MAXIONS][NHPTS][TPTS][GALPTS],gall[GALPTS];
float IonFrac(T,rho_H,galfac,ionid)
float T,rho_H,galfac;
int ionid;
{
	float frac;
	int calc_fractions();

	calc_fractions(T,rho_H,galfac,ionid,&frac);
	if( frac < 1.e-20 ) frac = 0.0;

#ifdef NEUTRALGAS 
#ifdef DO6IONS 
	if(ionid==0 || ionid==3){
	  frac = 1.0;
	}else{
	  frac = 0.0;
	}
#endif
#ifdef DOHIZIONS 
	if(ionid==0 || ionid==2 || ionid==3){
	  frac = 1.0;
	}else{
	  frac = 0.0;
	}
#endif
#endif


	return frac;
}

#endif



/* ------------------------------------------------------
   Calculate log(fractions) of HI,HeI,HeII,CII,CIV,O6,N5,SiIV
   by bilinear table interpolation in hydrogen density 
   and temperature. UH August 1996
   ------------------------------------------------------*/
#ifndef VARGALBKGD
int calc_fractions(temp, density, ionid, ionfraction)
  float temp;
  float density;
  int ionid;
  float *ionfraction;
{
  int inh,itemp;
  float y1,y2,y3,y4,t,u;
  float n_h;
  static float nhmax=0.000001,nhmin=999,tmax=11.,tmin=1.e+07;
  /* static int outbound=0;  */
  n_h = density/MHYDR;
  if(n_h<nhmin) nhmin=n_h;
  if(n_h>nhmax) nhmax=n_h;
  if(temp<tmin) tmin=temp;
  if(temp>tmax) tmax=temp;
/* --------------------------------------------------------
 * n_h is cgs number density of Hydrogen
 * temp is temperature in K
 * --------------------------------------------------------*/
 
 inh = (log10(n_h) - NHLOW)/DELTANH ;
 itemp = (log10(temp) - TLOW)/DELTAT ;
/* if((inh>NHPTS-2 || itemp>TPTS-2|| inh<0 || itemp <0) && ionid == 0 ){
 fprintf(stderr,"%d: n_h=%g %d, T=%g %d\n",outbound,n_h,inh,temp,itemp);
 outbound++;
 }*/
  if(inh>NHPTS-2) inh=NHPTS-2;
  if(itemp>TPTS-2)itemp=TPTS-2;
  if(inh<0)inh=0;
  if(itemp<0)itemp=0;
  t = (log10(n_h)-nhl[inh])/DELTANH;
  u = (log10(temp)-tl[itemp])/DELTAT;

  y1 = fractab[ionid][inh][itemp];
  y2 = fractab[ionid][inh+1][itemp];
  y3 = fractab[ionid][inh+1][itemp+1];
  y4 = fractab[ionid][inh][itemp+1];
  *ionfraction = pow(10.,(1.-t)*(1.-u)*y1+t*(1.-u)*y2+t*u*y3+(1.-t)*u*y4);
  if( *ionfraction < 1.e-10) fprintf(stderr,"%d: n_h=%f %d, T=%f %d  %g %g  %g\n",ionid,log10(n_h),inh,log10(temp),itemp,y1,y3,*ionfraction);

  return 0;
}
#else
int calc_fractions(temp, density, galfac, ionid, ionfraction)
  float temp;
  float density;
  float galfac;
  int ionid;
  float *ionfraction;
{
  int inh,itemp,igal;
  float y1,y2,y3,y4,y5,y6,y7,y8,t,u,v;
  float n_h;
  static float nhmax=0.000001,nhmin=999,tmax=11.0,tmin=1.e+07,galmax=-100,galmin=100;
  static int outbound=0; 
  n_h = density/MHYDR;
  if(n_h<nhmin) nhmin=n_h;
  if(n_h>nhmax) nhmax=n_h;
  if(temp<tmin) tmin=temp;
  if(temp>tmax) tmax=temp;
  if(galfac<galmin) galmin=galfac;
  if(galfac>galmax) galmax=galfac;

/* --------------------------------------------------------
 * n_h is cgs number density of Hydrogen
 * temp is temperature in K
 * --------------------------------------------------------*/
 
 inh = (log10(n_h) - NHLOW)/DELTANH ;
 itemp = (log10(temp) - TLOW)/DELTAT ;
 igal = (galfac - GALLOW)/DELTAGAL ;
/* if((inh>NHPTS-2 || itemp>TPTS-2|| inh<0 || itemp <0) && ionid == 0 ){
 fprintf(stderr,"%d: n_h=%g %d, T=%g %d\n",outbound,n_h,inh,temp,itemp);
 outbound++;
 }*/
  if(inh>NHPTS-2) inh=NHPTS-2;
  if(itemp>TPTS-2)itemp=TPTS-2;
  if(igal>GALPTS-1) igal=GALPTS-1;
  if(inh<0)inh=0;
  if(itemp<0)itemp=0;
  t = (log10(n_h)-nhl[inh])/DELTANH;
  u = (log10(temp)-tl[itemp])/DELTAT;
  v = (galfac-gall[igal])/DELTAGAL;

  y1 = fractab[ionid][inh][itemp][igal];
  y2 = fractab[ionid][inh+1][itemp][igal];
  y3 = fractab[ionid][inh+1][itemp+1][igal];
  y4 = fractab[ionid][inh][itemp+1][igal];
  y5 = fractab[ionid][inh][itemp][igal+1];
  y6 = fractab[ionid][inh+1][itemp][igal+1];
  y7 = fractab[ionid][inh+1][itemp+1][igal+1];
  y8 = fractab[ionid][inh][itemp+1][igal+1];
  
  if(igal>6)fprintf(stderr,"ionid=%d inh=%d (% 5.3e) itemp=%d (% 5.3e) igal=%d (% 5.3f) \n",ionid,inh,n_h,itemp,temp,igal,galfac);
  *ionfraction = pow(10.,(1.-t)*(1.-u)*(1.-v)*y1+t*(1.-u)*(1.-v)*y2+t*u*(1.-v)*y3+(1.-t)*u*(1.-v)*y4 + (1.-t)*(1.-u)*v*y5+t*(1.-u)*v*y6+t*u*v*y7+(1.-t)*u*v*y8);
  if( *ionfraction < 1.e-10 ) fprintf(stderr,"%d: n_h=%f %d, T=%f %d  %g %g  %g\n",outbound,log10(n_h),inh,log10(temp),itemp,y1,y3,*ionfraction);
  
  return 0;
}
#endif

#ifndef VARGALBKGD
void load_fraction_tables()
/*
  --------------------------------------------------------
  load_fraction_tables: reads in the tables of fractions
  of species (by number relative to the total number 
  of nuclei of the given element) as a function of
  density and temperature. At the present time, the
  following species are included: HI, HeI, HeII, CII, CIV,
  O6, NV, SiIV. UH August 1996 
  --------------------------------------------------------- 
*/
{
  int i;
  FILE *table1,*table2;
  int inh,itemp;
  int z10,zlo,zhi;
  float ftmp;
  float redzlo,redzhi;
  char fname1[80],fname2[80],prefix[50];
  char suffix[12];
  
	sprintf(prefix,"/home/shuiyao/temp/specexbin/ionfiles/");
#ifdef DOC4ONLY
	sprintf(suffix,"_c4");
#else
#ifdef DOLYAONLY
	sprintf(suffix,"_lya");
#else
#ifdef DOHANDHEONLY
	sprintf(suffix,"_hhe");
#else
#ifdef DO5IONS
	sprintf(suffix,"_i5");
#else
#ifdef DO6IONS
	sprintf(suffix,"_i6");
#else
#ifdef DO9IONS
	sprintf(suffix,"_i9");
#else
#ifdef DOXRAYIONS
	sprintf(suffix,"_xray");
#else
#ifdef PIPELINE
	sprintf(suffix,"_pipe");
#else
#ifdef DOHIZIONS
	sprintf(suffix,"_hiz");
#else
	sprintf(suffix,"_i31");
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#ifdef DOO6ONLY
        sprintf(suffix,"");
#endif

	z10 = (int) (redshift*10.);
	zlo = z10;
	sprintf(fname1,"%slt%.2d%s%s",prefix,z10,ionbkgd,suffix);
  	while ((table1=fopen(fname1,"r")) == NULL && zlo > 0) {
		zlo--;
		sprintf(fname1,"%slt%.2d%s%s",prefix,zlo,ionbkgd,suffix);
	}
	zhi = z10+1;
	sprintf(fname2,"%slt%.2d%s%s",prefix,zhi,ionbkgd,suffix);
  	while ((table2=fopen(fname2,"r")) == NULL && zhi < 80) {
	  zhi++;
		sprintf(fname2,"%slt%.2d%s%s",prefix,zhi,ionbkgd,suffix);
	}

  	if ( table1 == NULL && table2 == NULL ) {
     		printf("<unable to open files %s %s>\n",fname1,fname2);
     		exit(-1);
     	}
  	if ( table1 == NULL ) fprintf(stderr,"Using lookup table:%s\n",fname2);
	else if( table2 == NULL ) fprintf(stderr,"Using lookup tables:%s\n",fname1);
	else fprintf(stderr,"Using lookup tables:\n %s\n %s\n",fname1,fname2);
	if( table2 == NULL ) {
		table2 = table1;
		table1 = NULL;
	}
	redzlo = 0.1*zlo;
	redzhi = 0.1*zhi;
	fprintf(stderr,"Ionfrac: redshift = %5.3f redzlo = %5.3f redzhi = %5.3f\n",redshift,redzlo,redzhi);
	fprintf(stderr,"%g %g %g\n",redshift,redzlo,redzhi);


	if( table1 != NULL )
   	for(inh=0;inh<NHPTS;inh++) {
      		for(itemp=0;itemp<TPTS;itemp++)
         		for(i=0;i<nions;i++)
			  if(!fscanf(table1,"%g ", &fractab[i][inh][itemp]))
			    io_error_msg("ionfrac.c: table1 not read.");
	}
	if( table1 != NULL ) fclose(table1);

   	for(inh=0;inh<NHPTS;inh++) {
      		for(itemp=0;itemp<TPTS;itemp++) {
         		for(i=0;i<nions;i++) {
			  if(!fscanf(table2,"%g ", &ftmp))
			    io_error_msg("ionfrac.c: table1 not read.");
				if( table1 != NULL ) fractab[i][inh][itemp] = ((redzhi-redshift)*fractab[i][inh][itemp] + (redshift-redzlo)*ftmp)/(redzhi-redzlo);
				else fractab[i][inh][itemp] = ftmp;
			}
/*         		for(i=1;i<nions;i++) if(fractab[i][inh][itemp]<-6.) fractab[i][inh][itemp]=-6.;*/
		}
	} 
	fclose(table2);

/* Initialize gridpoint vectors for n_h and T */
	for(inh=0;inh<NHPTS;inh++) nhl[inh] = NHLOW + inh*DELTANH;
	for(itemp=0;itemp<TPTS;itemp++) tl[itemp] = TLOW + itemp*DELTAT;
      
   return;
}
#else
void load_fraction_tables()
/*
  --------------------------------------------------------
  load_fraction_tables: reads in the tables of fractions
  of species (by number relative to the total number 
  of nuclei of the given element) as a function of
  density and temperature. At the present time, the
  following species are included: HI, HeI, HeII, CII, CIV,
  O6, NV, SiIV. UH August 1996 
  --------------------------------------------------------- 
*/
{
  int i,j,k;
  FILE *table1,*table2;
  int inh,itemp,igal;
  int z10,zlo,zhi;
  int qfac;
  float gfac;
  float ftmp;
  float redzlo,redzhi;
  float hold;
  char fname1[80],fname2[80],prefix[50];
  char suffix[12];

  strcpy(prefix, "/data/collab/aford/ionfiles/");
#ifdef DOC4ONLY
  sprintf(suffix,"_c4");
#else
#ifdef DOLYAONLY
  sprintf(suffix,"_lya");
#else
#ifdef DOHANDHEONLY
  sprintf(suffix,"_hhe");
#else
#ifdef DO5IONS
  sprintf(suffix,"_i5");
#else
#ifdef DO6IONS
  sprintf(suffix,"_i6");
#else
#ifdef DO9IONS
  sprintf(suffix,"_i9");
#else
#ifdef DOXRAYIONS
  sprintf(suffix,"_xray");
#else
#ifdef PIPELINE
  sprintf(suffix,"_pipe");
#else
#ifdef DOHIZIONS
  sprintf(suffix,"_hiz");
#else
  sprintf(suffix,"_i31");
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#ifdef DOO6ONLY
	sprintf(suffix,"");
#endif

  for(igal=0;igal<GALPTS;igal++) {
    gfac = GALLOW + DELTAGAL*igal;
    z10 = (int) (redshift*10.);
    zlo = z10;
    sprintf(fname1,"%slt%.2dG%.3d%s",prefix,z10,(int)(gfac*10),suffix);
    while ((table1=fopen(fname1,"r")) == NULL && zlo > 0) {
      zlo--;
      sprintf(fname1,"%slt%.2dG%.3d%s",prefix,zlo,(int)(gfac*10),suffix);
    }
    zhi = z10+1;
    sprintf(fname2,"%slt%.2dG%.3d%s",prefix,zhi,(int)(gfac*10),suffix);
    while ((table2=fopen(fname2,"r")) == NULL && zhi < 80) {
      zhi++;
      sprintf(fname2,"%slt%.2dG%.3d%s",prefix,zhi,(int)(gfac*10),suffix);
    }
    if ( table1 == NULL && table2 == NULL ) {
      fprintf(stderr,"<unable to open files %s %s>\n",fname1,fname2);
      exit(-1);
    }
    if ( table1 == NULL ) fprintf(stderr,"Using lookup table:%s\n",fname2);
    else if( table2 == NULL ) fprintf(stderr,"Using lookup tables:%s\n",fname1);
    else fprintf(stderr,"Using lookup tables:\n %s\n %s\n",fname1,fname2);
    if( table2 == NULL ) {
      table2 = table1;
      table1 = NULL;
    }
    redzlo = 0.1*zlo;
    redzhi = 0.1*zhi;
    fprintf(stderr,"Ionfrac: redshift = %5.3f redzlo = %5.3f redzhi = %5.3f\n",redshift,redzlo,redzhi);
    fprintf(stderr,"%g %g %g\n",redshift,redzlo,redzhi);

    if( table1 != NULL )
      for(inh=0;inh<NHPTS;inh++) {
	for(itemp=0;itemp<TPTS;itemp++)
	  for(i=0;i<nions;i++)fscanf(table1,"%g ", &fractab[i][inh][itemp][igal]);
      }
    if( table1 != NULL ) fclose(table1);
    
    for(inh=0;inh<NHPTS;inh++) {
      for(itemp=0;itemp<TPTS;itemp++) {
	for(i=0;i<nions;i++) {
	  fscanf(table2,"%g ", &ftmp);
	  if( table1 != NULL ) fractab[i][inh][itemp][igal] = ((redzhi-redshift)*fractab[i][inh][itemp][igal] + (redshift-redzlo)*ftmp)/(redzhi-redzlo);
	  else fractab[i][inh][itemp][igal] = ftmp;
	}
	/*         		for(i=1;i<nions;i++) if(fractab[i][inh][itemp]<-6.) fractab[i][inh][itemp]=-6.;*/
      }
    } 
    fclose(table2);
  }

/* Initialize gridpoint vectors for n_h and T */
    for(inh=0;inh<NHPTS;inh++) nhl[inh] = NHLOW + inh*DELTANH;
    for(itemp=0;itemp<TPTS;itemp++) tl[itemp] = TLOW + itemp*DELTAT;
    for(igal=0;igal<GALPTS;igal++) gall[igal] = GALLOW + igal*DELTAGAL;
    return;
}
#endif

#ifdef NONEQUIL 

#define NZNEQ 25
#define NTNEQ 51
#define DELTAZNEQ 0.10
#define DELTATNEQ 0.05
#define ZNEQLOW -2.0
#define TNEQLOW 4.0

static float Znoneq[NZNEQ],tnoneq[NTNEQ];

static float noneqtab[MAXIONS][NZNEQ][NTNEQ];

float Noneq_Frac(T,Z,ionid)
     float T,Z;
     int ionid;
{
  float frac;
  int calc_noneq_fractions();

  if(Z>0){
    calc_noneq_fractions(T, Z/0.0189, ionid, &frac);
  }else{
    frac = 0.0;
  }
  if(frac<1e-10) frac = 0.0;
  if(Z>0){
    printf("NONEQ: %5.3f %5.3f % 5.3f %d %5.3e ",redshift,log10(T),log10((Z+1e-09)/0.0189),ionid,frac);
  }  
  
  return frac;

}

int calc_noneq_fractions(temp, Z, ionid, ionfraction)
     float temp;
     float Z;
     int ionid;
     float *ionfraction;
{
  int iZ,itemp;
  float y1,y2,y3,y4,t,u;
  static float Zmax=0.0001,Zmin=0.1,tmax=11,tmin=1.e+07;

  if(Z<Zmin) Zmin=Z;
  if(Z>Zmax) Zmax=Z;
  if(temp<tmin) tmin=temp;
  if(temp>tmax) tmax=temp;

  iZ = (log10(Z) - ZNEQLOW)/DELTAZNEQ ;
  itemp = (log10(temp) - TNEQLOW)/DELTATNEQ ;
  if(iZ>NZNEQ-2) iZ=NZNEQ-2;
  if(itemp>NTNEQ-2)itemp=NTNEQ-2;
  if(iZ<0)iZ=0;
  if(itemp<0)itemp=0;
  t = (log10(Z)-Znoneq[iZ])/DELTAZNEQ;
  u = (log10(temp)-tnoneq[itemp])/DELTATNEQ;

  y1 = noneqtab[ionid][iZ][itemp];
  y2 = noneqtab[ionid][iZ+1][itemp];
  y3 = noneqtab[ionid][iZ+1][itemp+1];
  y4 = noneqtab[ionid][iZ][itemp+1];
  *ionfraction = pow(10.,(1.-t)*(1.-u)*y1+t*(1.-u)*y2+t*u*y3+(1.-t)*u*y4);
  //if( *ionfraction > 0.5) fprintf(stderr,"%d: Z=%f %d, T=%f %d  %g %g  %g\n",ionid,log10(Z),iZ,log10(temp),itemp,y1,y3,*ionfraction);
  if( *ionfraction < 1.e-10) fprintf(stderr,"%d: Z=%f %d, T=%f %d  %g %g  %g\n",ionid,log10(Z),iZ,log10(temp),itemp,y1,y3,*ionfraction);

  return 0;
}

void load_noneq_fraction_tables()
{
  int i;
  FILE *table;
  int iZ,itemp;
  float tmp1,tmp2;
  char fname[80],prefix[50];

  sprintf(prefix,"/data/collab/aford/ionfiles/");

  sprintf(fname,"%slt_noneq",prefix);

  if((table=fopen(fname,"r")) == NULL){
    fprintf(stderr,"Cannot locate %s\n",fname);
    exit(-1);
  }

  if( table != NULL ){ 
    for(iZ=0;iZ<NZNEQ;iZ++) {
      for(itemp=0;itemp<NTNEQ;itemp++) {
	fscanf(table,"%g %g ",&tmp1,&tmp2);
	//printf("%5.3f % 5.3f ",tmp1, tmp2);
	for(i=0;i<nions;i++) {
	  fscanf(table,"%g ",&noneqtab[i][iZ][itemp]);
	  //printf("noneqtab[%d][%d][%d]= % 5.3e ",i,iZ,itemp,noneqtab[i][iZ][itemp]);
	}
	//printf("\n");
      }
    }
  }
  fclose(table);

  for(iZ=0;iZ<NZNEQ;iZ++) Znoneq[iZ] = ZNEQLOW + iZ*DELTAZNEQ;
  for(itemp=0;itemp<NTNEQ;itemp++) tnoneq[itemp] = TNEQLOW + itemp*DELTATNEQ;

}


#endif
