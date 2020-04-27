#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "specexbin.defs"

ionStruct *Ion;

int main(int argc,char **argv)
{



  
  if( (infile = fopen(argv[1],"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",argv[1]);
    return -1;
  }
    
  while(1){
    if(feof(infile)) break;
    fscanf(infile,"%f%f%f%f",&z[i],&rho[i],&T[i],&Z[i]);
    for(j=0;j<nions;j++){
      fscanf(infile,"%f",&tau[i][j]);
    }
    i++;
  }


  
}
 

int loadspecions()
{

  int i;
  char line[80],prefix[50],specionfilename[80];
  FILE *specfile;

  sprintf(prefix,"/data/collab/aford/ionfiles/");
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
#ifdef PIPELINE
  sprintf(specionfilename,"%sspecions_pipe.dat",prefix);
#else
#ifdef DOHIZIONS
  sprintf(specionfilename,"%sspecions_hiz.dat",prefix);
#else
  sprintf(specionfilename,"%sspecions_i29.dat",prefix);
#endif
#endif
#endif
#endif
#endif
#endif
#endif
  if( (specfile = fopen(specionfilename,"r")) == NULL ) specfile = fopen(specionfilename,"r");

  if( specfile == NULL ) {
    fprintf(stderr,"cannot find specion file anywhere\n");
    exit(-1);
  }


  Ion = (ionStruct *) malloc(MAXIONS*sizeof(ionStruct));

  i = 0;
  while( fgets(line,80,infile) != NULL ) {
    if( strstr(line,"#") != NULL ) continue;
    if( i >= MAXIONS ) break;
    sscanf(line,"%10s %g %g %g %g %d %g",Ion[i].name,&Ion[i].lambda,&Ion[i].Xsec,&Ion[i].atomwt,&Ion[i].fraction,&Ion[i].Zcolumn,&Ion[i].alpha);
    i++;
  }
  nions = i;
  fclose(infile);

  fprintf(stderr,"Processing %d ions from specions.dat:\n",nions);
  for( i=0; i<nions; i++ ) {
    Ion[i].bsys = sqrt(2.*KBOLTZ/(MHYDR*Ion[i].atomwt))/1.e5;
    Ion[i].Xsec *= 2.648e-2*Ion[i].lambda*1.e-13;
    fprintf(stderr,"%5d %10s %12.6g %12.6g %10.5g %10.5g %10.5g % 3.1f % 2d\n",i,Ion[i].name,Ion[i].lambda,Ion[i].Xsec,Ion[i].atomwt,Ion[i].fraction,Ion[i].bsys*sqrt(1.e4),Ion[i].alpha,Ion[i].Zcolumn);
  }

  return nions;
}
