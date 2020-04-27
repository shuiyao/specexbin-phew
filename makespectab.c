/* This program takes a list of tipsy snapshots and redshifts (see
   doIGM script) and outputs the tabular list of snapshots, redshifts,
   redshift ranges, and model base names, which is an input into
   specexbin.  BDO 11/16/08
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXFILES 200

int main(int argc,char **argv)
{
  int i, nbin;
  float zhi, zlow;
  float z[MAXFILES];
  float zhirange, zlowrange;
  char modelname[MAXFILES][100], gadgetname[MAXFILES][300];
  FILE *zfile;
  
 if( (zfile = fopen(argv[1],"r")) == NULL ) {
   fprintf(stderr,"Could not open file %s\n",argv[1]);
   return -1;
 }


  zhi = 100.0;
  zlow = 0.0;
  i = 0;
  while(1){
    if(feof(zfile)) break;
    fscanf(zfile,"%f%s%s",&z[i],&modelname[i],&gadgetname[i]);
    i++;
  }
  nbin = i-1;

  for(i=0;i<nbin;i++){
    if(i==0) 
      zhirange = zhi;
    else 
      zhirange = (z[i-1]+z[i])/2;
    
    if(i==nbin-1)
      zlowrange = zlow;
    else
      zlowrange = (z[i]+z[i+1])/2;

    printf("%s %f %f %f %s\n",gadgetname[i],z[i],zlowrange,zhirange,modelname[i]);
  } 
  
}
 
