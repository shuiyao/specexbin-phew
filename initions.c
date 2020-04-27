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

/* Read from file specions.dat:
   ion_name  rest_freq  osc_strength  atom_wt  solar_abundance  rel_abundance
 * solar abundance mass fractions are taken from Arnett, "Supernovae 
	and Nucleosynthesis" (1996) Table A.1 column 4,
	which in turn is from Anders & Grevesse 1989 and others.
 * rel_abundance is in dex vs.solar (not used in specgen)
*/

int InitIons()
{
  int i,j;
  char line[80],prefix[50],specionfilename[80];
  int load_fraction_tables();
  FILE *specfile;

  sprintf(prefix,"/home/shuiyao/temp/specexbin/ionfiles/");
  /* nvbins = floor(nzbins/(VRES/ZRES)) + NVBINS_ADVANCED; */
  nvbins = floor(nzbins/(VRES/ZRES));

  fprintf(stderr,"nvbins = %d, nvloopbins = %d\n",nvbins, nvloopbins);
#ifdef DOC4ONLY
  sprintf(specionfilename,"%sspecions_c4.dat",prefix);
#else
#ifdef DOLYAONLY
  sprintf(specionfilename,"%sspecions_lya.dat",prefix);
#else
#ifdef DOHANDHEONLY
  sprintf(specionfilename,"%sspecions_hhe.dat",prefix);
#else
#ifdef DOO6ONLY
  sprintf(specionfilename,"%sspecions_o6.dat",prefix);
#else
#ifdef DO5IONS
  sprintf(specionfilename,"%sspecions_i5.dat",prefix);
#else
#ifdef DO6IONS
  sprintf(specionfilename,"%sspecions_i6.dat",prefix);
#else
#ifdef DO9IONS
  sprintf(specionfilename,"%sspecions_i9.dat",prefix);
#else
#ifdef DOXRAYIONS
  sprintf(specionfilename,"%sspecions_xray.dat",prefix);
#else
#ifdef PIPELINE
  sprintf(specionfilename,"%sspecions_pipe.dat",prefix);
#else
#ifdef DOHIZIONS
  sprintf(specionfilename,"%sspecions_hiz.dat",prefix);
#else
  sprintf(specionfilename,"%sspecions_i31.dat",prefix);
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
  if( (specfile = fopen(specionfilename,"r")) == NULL ) specfile = fopen(specionfilename,"r");

  //if( (specfile = fopen("specions_i5.dat","r")) == NULL ) specfile = fopen("specions_i5.dat","r");

  if( specfile == NULL ) {
    fprintf(stderr,"cannot find specion file anywhere\n");
    exit(-1);
  }
  if(loop==0){
    fprintf(stderr,"FIRST LOOP ALLOCATION\n");
    nvbins += nvloopbins; // In the case of NVBINS_ADVANCE
    Ion = (ionStruct *) malloc(MAXIONS*sizeof(ionStruct));
    for(i=0;i<MAXIONS;i++){
      Ion[i].mass = malloc(nzbins*sizeof(double));
      Ion[i].vel = malloc(nzbins*sizeof(double));
      Ion[i].temp = malloc(nzbins*sizeof(double));
      Ion[i].rho = malloc(nzbins*sizeof(double));
      for(j=0;j<NMETALS;j++)Ion[i].metals[j] = malloc(nzbins*sizeof(double));
#ifdef PHYSSPEC
      Ion[i].sfr = malloc(nzbins*sizeof(double));
      Ion[i].wtmass = malloc(nzbins*sizeof(double));
      Ion[i].mgal = malloc(nzbins*sizeof(double));
      Ion[i].dgal = malloc(nzbins*sizeof(double));
      Ion[i].age = malloc(nzbins*sizeof(double));
      Ion[i].nrec = malloc(nzbins*sizeof(double));
      Ion[i].vlaunch = malloc(nzbins*sizeof(double));
#endif
      Ion[i].vbins = malloc(nvbins*sizeof(double));
      Ion[i].tbins = malloc(nvbins*sizeof(double));
      Ion[i].rhobins = malloc(nvbins*sizeof(double));
      Ion[i].Zbins = malloc(nvbins*sizeof(double));      
#if defined(PHEW) || defined(PART_BY_PART)
      Ion[i].vcbins = malloc(nvbins*sizeof(double));
      Ion[i].tcbins = malloc(nvbins*sizeof(double));
      Ion[i].rhocbins = malloc(nvbins*sizeof(double));
      Ion[i].Zcbins = malloc(nvbins*sizeof(double));      
#endif      
#ifdef PHYSSPEC
      Ion[i].sfrbins = malloc(nvbins*sizeof(double));
      Ion[i].mgalbins = malloc(nvbins*sizeof(double));
      Ion[i].dgalbins = malloc(nvbins*sizeof(double));
      Ion[i].agebins = malloc(nvbins*sizeof(double));
      Ion[i].nrecbins = malloc(nvbins*sizeof(double));
      Ion[i].vlaunchbins = malloc(nvbins*sizeof(double));
#endif
    }
    IonTotal.mass = malloc(nzbins*sizeof(double));
    IonTotal.vel = malloc(nzbins*sizeof(double));
    IonTotal.temp = malloc(nzbins*sizeof(double));
    IonTotal.rho = malloc(nzbins*sizeof(double));
    for(j=0;j<NMETALS;j++)IonTotal.metals[j] = malloc(nzbins*sizeof(double));
#ifdef PHYSSPEC
    IonTotal.sfr = malloc(nzbins*sizeof(double));
    IonTotal.wtmass = malloc(nzbins*sizeof(double));
    IonTotal.mgal = malloc(nzbins*sizeof(double));
    IonTotal.dgal = malloc(nzbins*sizeof(double));
    IonTotal.age = malloc(nzbins*sizeof(double));
    IonTotal.nrec = malloc(nzbins*sizeof(double));
    IonTotal.vlaunch = malloc(nzbins*sizeof(double));
#endif
    IonTotal.redshift = malloc(nzbins*sizeof(double));
    IonTotal.binsize = malloc(nzbins*sizeof(double));
    IonTotal.bincoord = malloc(nzbins*sizeof(double));
    IonTotal.vbins = malloc(nvbins*sizeof(double));
    IonTotal.tbins = malloc(nvbins*sizeof(double));
    IonTotal.rhobins = malloc(nvbins*sizeof(double));
    IonTotal.Zbins = malloc(nvbins*sizeof(double));
#if defined(PHEW) || defined(PART_BY_PART)
    IonTotal.vcbins = malloc(nvbins*sizeof(double));
    IonTotal.tcbins = malloc(nvbins*sizeof(double));
    IonTotal.rhocbins = malloc(nvbins*sizeof(double));
    IonTotal.Zcbins = malloc(nvbins*sizeof(double));    
#endif      
#ifdef PHYSSPEC
    IonTotal.sfrbins = malloc(nvbins*sizeof(double));
    IonTotal.mgalbins = malloc(nvbins*sizeof(double));
    IonTotal.dgalbins = malloc(nvbins*sizeof(double));
    IonTotal.agebins = malloc(nvbins*sizeof(double));
    IonTotal.nrecbins = malloc(nvbins*sizeof(double));
    IonTotal.vlaunchbins = malloc(nvbins*sizeof(double));
#endif
    IonExtra.redshift = malloc(nvbins*sizeof(double));
    IonExtra.x = malloc(nzbins*sizeof(double));
    IonExtra.y = malloc(nzbins*sizeof(double));
    IonExtra.z = malloc(nzbins*sizeof(double));
    IonExtra.xbins = malloc(nvbins*sizeof(double));
    IonExtra.ybins = malloc(nvbins*sizeof(double));
    IonExtra.zbins = malloc(nvbins*sizeof(double));
    IonExtra.gal_field = malloc(nvbins*sizeof(float));
    nvbins -= nvloopbins;
  }else{
    fprintf(stderr,"SUBSEQUENT LOOP ALLOCATION: nzbins = %d, nzloopbins = %d\n",nzbins,nzloopbins);
    for(i=0;i<MAXIONS;i++){
      Ion[i].mass = realloc(Ion[i].mass,(nzloopbins+nzbins)*sizeof(double));
      Ion[i].vel = realloc(Ion[i].vel,(nzloopbins+nzbins)*sizeof(double));
      Ion[i].temp = realloc(Ion[i].temp,(nzloopbins+nzbins)*sizeof(double));
      Ion[i].rho = realloc(Ion[i].rho,(nzloopbins+nzbins)*sizeof(double));
      for(j=0;j<NMETALS;j++)Ion[i].metals[j] = realloc(Ion[i].metals[j],(nzloopbins+nzbins)*sizeof(double));
#ifdef PHYSSPEC
      Ion[i].sfr = realloc(Ion[i].sfr,(nzloopbins+nzbins)*sizeof(double));
      Ion[i].wtmass = realloc(Ion[i].wtmass,(nzloopbins+nzbins)*sizeof(double));
      Ion[i].mgal = realloc(Ion[i].mgal,(nzloopbins+nzbins)*sizeof(double));
      Ion[i].dgal = realloc(Ion[i].dgal,(nzloopbins+nzbins)*sizeof(double));
      Ion[i].age = realloc(Ion[i].age,(nzloopbins+nzbins)*sizeof(double));
      Ion[i].nrec = realloc(Ion[i].nrec,(nzloopbins+nzbins)*sizeof(double));
      Ion[i].vlaunch = realloc(Ion[i].vlaunch,(nzloopbins+nzbins)*sizeof(double));
#endif
      Ion[i].vbins = realloc(Ion[i].vbins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].tbins = realloc(Ion[i].tbins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].rhobins = realloc(Ion[i].rhobins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].Zbins = realloc(Ion[i].Zbins,(nvloopbins+nvbins)*sizeof(double));
#if defined(PHEW) || defined(PART_BY_PART)
      Ion[i].vcbins = realloc(Ion[i].vcbins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].tcbins = realloc(Ion[i].tcbins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].rhocbins = realloc(Ion[i].rhocbins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].Zcbins = realloc(Ion[i].Zcbins,(nvloopbins+nvbins)*sizeof(double));      
#endif      
#ifdef PHYSSPEC
      Ion[i].sfrbins = realloc(Ion[i].sfrbins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].mgalbins = realloc(Ion[i].mgalbins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].dgalbins = realloc(Ion[i].dgalbins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].agebins = realloc(Ion[i].agebins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].nrecbins = realloc(Ion[i].nrecbins,(nvloopbins+nvbins)*sizeof(double));
      Ion[i].vlaunchbins = realloc(Ion[i].vlaunchbins,(nvloopbins+nvbins)*sizeof(double));
#endif
    }
    IonTotal.mass = realloc(IonTotal.mass,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.vel = realloc(IonTotal.vel,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.temp = realloc(IonTotal.temp,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.rho = realloc(IonTotal.rho,(nzloopbins+nzbins)*sizeof(double));
    for(j=0;j<NMETALS;j++)IonTotal.metals[j] = realloc(IonTotal.metals[j],(nzloopbins+nzbins)*sizeof(double));
#ifdef PHYSSPEC
    IonTotal.sfr = realloc(IonTotal.sfr,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.wtmass = realloc(IonTotal.wtmass,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.mgal = realloc(IonTotal.mgal,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.dgal = realloc(IonTotal.dgal,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.age = realloc(IonTotal.age,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.nrec = realloc(IonTotal.nrec,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.vlaunch = realloc(IonTotal.vlaunch,(nzloopbins+nzbins)*sizeof(double));
#endif
    IonTotal.redshift = realloc(IonTotal.redshift,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.binsize = realloc(IonTotal.binsize,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.bincoord = realloc(IonTotal.bincoord,(nzloopbins+nzbins)*sizeof(double));
    IonTotal.vbins = realloc(IonTotal.vbins,(nvloopbins+nvbins)*sizeof(double));
    IonTotal.tbins = realloc(IonTotal.tbins,(nvloopbins+nvbins)*sizeof(double));
    IonTotal.rhobins = realloc(IonTotal.rhobins,(nvloopbins+nvbins)*sizeof(double));
#if defined(PHEW) || defined(PART_BY_PART)
    IonTotal.vcbins = realloc(IonTotal.vcbins,(nvloopbins+nvbins)*sizeof(double));
    IonTotal.tcbins = realloc(IonTotal.tcbins,(nvloopbins+nvbins)*sizeof(double));
    IonTotal.rhocbins = realloc(IonTotal.rhocbins,(nvloopbins+nvbins)*sizeof(double));
    IonTotal.Zcbins = realloc(IonTotal.Zcbins,(nvloopbins+nvbins)*sizeof(double));    
#endif      
    IonTotal.Zbins = realloc(IonTotal.Zbins,(nvloopbins+nvbins)*sizeof(double));
    IonExtra.redshift = realloc(IonExtra.redshift,(nvloopbins+nvbins)*sizeof(double));
    IonExtra.x = realloc(IonExtra.x,(nzloopbins+nzbins)*sizeof(double));
    IonExtra.y = realloc(IonExtra.y,(nzloopbins+nzbins)*sizeof(double));
    IonExtra.z = realloc(IonExtra.z,(nzloopbins+nzbins)*sizeof(double));
    IonExtra.xbins = realloc(IonExtra.xbins,(nvloopbins+nvbins)*sizeof(double));
    IonExtra.ybins = realloc(IonExtra.ybins,(nvloopbins+nvbins)*sizeof(double));
    IonExtra.zbins = realloc(IonExtra.zbins,(nvloopbins+nvbins)*sizeof(double));
    IonExtra.gal_field = realloc(IonExtra.gal_field,(nvloopbins+nvbins)*sizeof(double));
#ifdef PHYSSPEC
    IonTotal.sfrbins = realloc(IonTotal.sfrbins,(nvloopbins+nvbins)*sizeof(double));
    IonTotal.mgalbins = realloc(IonTotal.mgalbins,(nvloopbins+nvbins)*sizeof(double));
    IonTotal.dgalbins = realloc(IonTotal.dgalbins,(nvloopbins+nvbins)*sizeof(double));
    IonTotal.agebins = realloc(IonTotal.agebins,(nvloopbins+nvbins)*sizeof(double));
    IonTotal.nrecbins = realloc(IonTotal.nrecbins,(nvloopbins+nvbins)*sizeof(double));
    IonTotal.vlaunchbins = realloc(IonTotal.vlaunchbins,(nvloopbins+nvbins)*sizeof(double));
#endif
  } 

#if defined(PHEW) || defined(PART_BY_PART)
  if(nvloopbins == NVBINS_ADVANCED) // First call
    for( i=0; i<nvloopbins; i++ ){
      IonTotal.vcbins[i] = IonTotal.tcbins[i] = IonTotal.rhocbins[i] = IonTotal.Zcbins[i] = 0.0;
      for(j=0;j<MAXIONS;j++){
	Ion[j].vcbins[i] = Ion[j].tcbins[i] = Ion[j].rhocbins[i] = Ion[j].Zcbins[i] = 0.0;
      }
    }
  for( i=nvloopbins; i<nvloopbins+nvbins; i++ ){
    IonTotal.vcbins[i] = IonTotal.tcbins[i] = IonTotal.rhocbins[i] = IonTotal.Zcbins[i] = 0.0;
    for(j=0;j<MAXIONS;j++){
      Ion[j].vcbins[i] = Ion[j].tcbins[i] = Ion[j].rhocbins[i] = Ion[j].Zcbins[i] = 0.0;
    }
  }
  fprintf(stderr, "=== Update Ion.vcbins From %d to %d ===\n", nvloopbins, nvloopbins+nvbins);
#endif
  
  i = 0;
  while( fgets(line,80,specfile) != NULL ) {
    if( strstr(line,"#") != NULL ) continue;
    if( i >= MAXIONS ) break;
    sscanf(line,"%10s %g %g %g %g %d %g",Ion[i].name,&Ion[i].lambda,&Ion[i].Xsec,&Ion[i].atomwt,&Ion[i].fraction,&Ion[i].Zcolumn,&Ion[i].alpha);
    i++;
  }
  nions = i;
  fclose(specfile);
  
  fprintf(stderr,"Processing %d ions from specions.dat:\n",nions);
  for( i=0; i<nions; i++ ) {
    Ion[i].bsys = sqrt(2.*KBOLTZ/(MHYDR*Ion[i].atomwt))/1.e5;
    Ion[i].Xsec *= 2.648e-2*Ion[i].lambda*1.e-13;
    // 1.497e-2 * sqrt(pi) = 2.653e-2
    fprintf(stderr,"%5d %10s %12.6g %12.6g %10.5g %10.5g %10.5g % 3.1f % 2d\n",i,Ion[i].name,Ion[i].lambda,Ion[i].Xsec,Ion[i].atomwt,Ion[i].fraction,Ion[i].bsys*sqrt(1.e4),Ion[i].alpha,Ion[i].Zcolumn);
  }
  
  return 0;
}

void FreeIons(void)
{
  int i,j;

  for(i=0;i<MAXIONS;i++){
    free(Ion[i].mass);
    free(Ion[i].vel);
    free(Ion[i].temp);
    free(Ion[i].rho);
    for(j=0;j<NMETALS;j++)free(Ion[i].metals[j]);
    free(Ion[i].vbins);
    free(Ion[i].tbins);
    free(Ion[i].rhobins);
    free(Ion[i].Zbins);
#if defined(PHEW) || defined(PART_BY_PART)
    free(Ion[i].vcbins);
    free(Ion[i].rhocbins);
    free(Ion[i].tcbins);
    free(Ion[i].Zcbins);        
#endif    
  }
  free(IonTotal.mass);
  free(IonTotal.vel);
  free(IonTotal.temp);
  free(IonTotal.rho);
  for(j=0;j<NMETALS;j++)free(IonTotal.metals[j]);
  free(IonTotal.redshift);
  free(IonTotal.binsize);
  free(IonTotal.bincoord);
  free(IonTotal.vbins);
  free(IonTotal.tbins);
  free(IonTotal.rhobins);
#if defined(PHEW) || defined(PART_BY_PART)
  free(IonTotal.vcbins);
  free(IonTotal.rhocbins);
  free(IonTotal.tcbins);
  free(IonTotal.Zcbins);      
#endif    
  free(IonTotal.Zbins);
  free(IonExtra.redshift);
  free(IonExtra.x);
  free(IonExtra.y);
  free(IonExtra.z);
  free(IonExtra.gal_field);
}
