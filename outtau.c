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

extern double unit_Velocity;

int OutTau()
{
  int i, k;
  /* int j, vi; */
  /* double z; */
  double rhomean;
  float CosmicTime();
  int cosmopar();
  FILE *outfile;
  char fname[80];
#ifndef OUTPUT_LOCAL_FOLDER
  char longfname[400];
  char spec_dir[100];
  strcpy(spec_dir, FOLDER_OUTPUT);
#endif

#ifndef SHORTSPEC
  int iname, jname, kname;
  iname = floor(theta*180.0/PI+0.001);
  kname = floor(redshift_begin*100+0.0001);
  //k = floor((redshift_begin-delta_redshift+0.01)*100+0.0001);
  jname = floor(redshift_end*100+0.0001);
#endif

#ifdef PIPELINE
  if(theta>0){
    sprintf(fname,"spectau.%s.%d.z%d_%d",sim_id,iname,jname,kname);
  }else{
    sprintf(fname,"spectau.%s.%s.z%d_%d",sim_id,id,jname,kname);
  }
#else
#ifdef SHORTSPEC
  sprintf(fname,"specaim.%s.%s",sim_id,namesuffix);
#else
  if(theta>0){
    sprintf(fname,"specztau.%s.%d.%d_%d",sim_id,iname,jname,kname);
  }else{
    sprintf(fname,"specztau.%s.%s.%d_%d",sim_id,id,jname,kname);
  }
#endif
#endif

#ifdef OUTPUT_LOCAL_FOLDER
  outfile = fopen(fname,"w");
#else
  strcat(strcpy(longfname, spec_dir), fname);
  outfile = fopen(longfname,"w");
#endif

  fprintf(stderr,"Outputting in Spectau!!!!\n");
  i = 0;
  redshift_track = IonTotal.redshift[0];
  fprintf(stderr,"redshift_track = %f\n",redshift_track);
  //delta_redshift -= 0.01;
  //while(redshift_track > redshift_begin && redshift_track < redshift_end){
  while(redshift_track>redshift_begin-0.01){
    redshift_track -= vres;
    //fprintf(stderr,"redshift_track = %f redshift_begin = %f redshift_end = %f\n",redshift_track,redshift_begin,redshift_end);
    if( redshift_track > redshift_begin && redshift_track < redshift_end){
    /* if(1>0){       */
      fprintf(outfile, "%10.7f ",redshift_track);
      /* fprintf(outfile, "%d ",i);       */
#ifdef SHORTSPEC
      cosmopar(CosmicTime(redshift_center));
#else
      cosmopar(CosmicTime(redshift_track));
#endif
      /* Convert rho to overdensities for output */
      rhomean = XH*1.88e-29*Omega_b*h*h/(aex*aex*aex);
      for(k=-1; k<nions; k++) {
	if(k==-1){
          IonTotal.rhobins[i] = log10(IonTotal.rhobins[i]/rhomean);
          IonTotal.tbins[i] = log10(IonTotal.tbins[i]);
	  if(isnan(IonTotal.rhobins[i]) || isinf(IonTotal.rhobins[i])) IonTotal.rhobins[i] = 0.0;
	  if(isnan(IonTotal.tbins[i]) || isinf(IonTotal.tbins[i])) IonTotal.tbins[i] = 0.0;
	}else{
	  Ion[k].rhobins[i] = log10(Ion[k].rhobins[i]/rhomean);
	  Ion[k].tbins[i] = log10(Ion[k].tbins[i]);
	  if(isnan(Ion[k].rhobins[i]) || isinf(Ion[k].rhobins[i])) Ion[k].rhobins[i] = 0.0;
	  if(isnan(Ion[k].tbins[i]) || isinf(Ion[k].tbins[i])) Ion[k].tbins[i] = 0.0;
	  if(isnan(Ion[k].Zbins[i]) || isinf(Ion[k].Zbins[i])) Ion[k].Zbins[i] = 0.0;
	  if(isnan(Ion[k].vbins[i]) || isinf(Ion[k].vbins[i])) Ion[k].vbins[i] = 0.0;
	}
	if(k==-1){
	  fprintf(outfile, "% 6.3f %6.3f %5.3e ", IonTotal.rhobins[i],IonTotal.tbins[i],IonTotal.Zbins[i]);
	}else{
#ifdef PIPELINE
	  fprintf(outfile, "%5.3e ",Ion[k].vbins[i]);
#else
	  fprintf(outfile, "% 5.2f %5.2f %5.3e %5.3e ", Ion[k].rhobins[i],Ion[k].tbins[i],Ion[k].Zbins[i],Ion[k].vbins[i]);
#endif
	}
      }
#ifndef PIPELINE
      fprintf(outfile,"% 7.5f % 7.5f % 7.5f",IonExtra.xbins[i],IonExtra.ybins[i],IonExtra.zbins[i]);
#endif
      fprintf(outfile,"\n");
    }
    i++;
  } // while(redshift_track>redshift_begin-0.01)
  fclose(outfile); // specztau or specaim

#ifdef PHEW
  // The format is exactly the same as specztau or specaim. Only contains wind information
#ifdef SHORTSPEC
  sprintf(fname,"specaimc.%s.%s",sim_id,namesuffix);
#else
  if(theta>0){
    sprintf(fname,"specztauc.%s.%d.%d_%d",sim_id,iname,jname,kname);
  }else{
    sprintf(fname,"specztauc.%s.%s.%d_%d",sim_id,id,jname,kname);
  }
#endif
#ifdef OUTPUT_LOCAL_FOLDER
  outfile = fopen(fname,"w");
#else
  strcat(strcpy(longfname, spec_dir), fname);
  outfile = fopen(longfname,"w");
#endif

  fprintf(stderr,"Outputting in Spectau (Clouds)!!!!\n");
  i = 0;
  redshift_track = IonTotal.redshift[0];
  while(redshift_track>redshift_begin-0.01){
    redshift_track -= vres;
    //fprintf(stderr,"redshift_track = %f redshift_begin = %f redshift_end = %f\n",redshift_track,redshift_begin,redshift_end);
    if( redshift_track > redshift_begin && redshift_track < redshift_end){
      fprintf(outfile, "%10.7f ",redshift_track);
#ifdef SHORTSPEC
      cosmopar(CosmicTime(redshift_center));
#else
      cosmopar(CosmicTime(redshift_track));
#endif
      /* Convert rho to overdensities for output */
      rhomean = XH*1.88e-29*Omega_b*h*h/(aex*aex*aex);
      for(k=-1; k<nions; k++) {
	if(k==-1){
          IonTotal.rhocbins[i] = log10(IonTotal.rhocbins[i]/rhomean);
          IonTotal.tcbins[i] = log10(IonTotal.tcbins[i]);
	  if(isnan(IonTotal.rhocbins[i]) || isinf(IonTotal.rhocbins[i])) IonTotal.rhocbins[i] = 0.0;
	  if(isnan(IonTotal.tcbins[i]) || isinf(IonTotal.tcbins[i])) IonTotal.tcbins[i] = 0.0;
	}else{
	  Ion[k].rhocbins[i] = log10(Ion[k].rhocbins[i]/rhomean);
	  Ion[k].tcbins[i] = log10(Ion[k].tcbins[i]);
	  if(isnan(Ion[k].rhocbins[i]) || isinf(Ion[k].rhocbins[i])) Ion[k].rhocbins[i] = 0.0;
	  if(isnan(Ion[k].tcbins[i]) || isinf(Ion[k].tcbins[i])) Ion[k].tcbins[i] = 0.0;
	  /* if(isnan(Ion[k].Zcbins[i]) || isinf(Ion[k].Zcbins[i])) Ion[k].Zcbins[i] = 0.0; */
	  if(isnan(Ion[k].vcbins[i]) || isinf(Ion[k].vcbins[i])) Ion[k].vcbins[i] = 0.0;
	}
	if(k==-1){
	  fprintf(outfile, "% 6.3f %6.3f 0.0 ", IonTotal.rhocbins[i],IonTotal.tcbins[i]);
	}else{
	  fprintf(outfile, "% 5.2f %5.2f 0.0 %5.3e ", Ion[k].rhocbins[i],Ion[k].tcbins[i],Ion[k].vcbins[i]);
	}
      }
      fprintf(outfile,"% 7.5f % 7.5f % 7.5f",IonExtra.xbins[i],IonExtra.ybins[i],IonExtra.zbins[i]);
      fprintf(outfile,"\n");
    }
    i++;
  } // while(redshift_track>redshift_begin-0.01)
  fclose(outfile); // specztau or specaim
#endif // PHEW

#ifdef PART_BY_PART 
  // The format is exactly the same as specztau or specaim. Only contains wind information
#ifdef SHORTSPEC
  sprintf(fname,"specaimw.%s.%s",sim_id,namesuffix);
#else
  if(theta>0){
    sprintf(fname,"specztauw.%s.%d.%d_%d",sim_id,iname,jname,kname);
  }else{
    sprintf(fname,"specztauw.%s.%s.%d_%d",sim_id,id,jname,kname);
  }
#endif
#ifdef OUTPUT_LOCAL_FOLDER
  outfile = fopen(fname,"w");
#else
  strcat(strcpy(longfname, spec_dir), fname);
  outfile = fopen(longfname,"w");
#endif

  fprintf(stderr,"Outputting in Spectau (Winds)!!!!\n");
  i = 0;
  redshift_track = IonTotal.redshift[0];
  while(redshift_track>redshift_begin-0.01){
    redshift_track -= vres;
    //fprintf(stderr,"redshift_track = %f redshift_begin = %f redshift_end = %f\n",redshift_track,redshift_begin,redshift_end);
    if( redshift_track > redshift_begin && redshift_track < redshift_end){
      /* Only write within the z boundary. (Exclude z < 0 for example) */
      fprintf(outfile, "%10.7f ",redshift_track);
      /* fprintf(outfile, "%d ",i); */
#ifdef SHORTSPEC
      cosmopar(CosmicTime(redshift_center));
#else
      cosmopar(CosmicTime(redshift_track));
#endif
      /* Convert rho to overdensities for output */
      rhomean = XH*1.88e-29*Omega_b*h*h/(aex*aex*aex);
      for(k=-1; k<nions; k++) {
	if(k==-1){
          IonTotal.rhocbins[i] = log10(IonTotal.rhocbins[i]/rhomean);
          IonTotal.tcbins[i] = log10(IonTotal.tcbins[i]);
	  if(isnan(IonTotal.rhocbins[i]) || isinf(IonTotal.rhocbins[i])) IonTotal.rhocbins[i] = 0.0;
	  if(isnan(IonTotal.tcbins[i]) || isinf(IonTotal.tcbins[i])) IonTotal.tcbins[i] = 0.0;
	}else{
	  Ion[k].rhocbins[i] = log10(Ion[k].rhocbins[i]/rhomean);
	  Ion[k].tcbins[i] = log10(Ion[k].tcbins[i]);
	  if(isnan(Ion[k].rhocbins[i]) || isinf(Ion[k].rhocbins[i])) Ion[k].rhocbins[i] = 0.0;
	  if(isnan(Ion[k].tcbins[i]) || isinf(Ion[k].tcbins[i])) Ion[k].tcbins[i] = 0.0;
	  if(isnan(Ion[k].Zcbins[i]) || isinf(Ion[k].Zcbins[i])) Ion[k].Zcbins[i] = 0.0;
	  if(isnan(Ion[k].vcbins[i]) || isinf(Ion[k].vcbins[i])) Ion[k].vcbins[i] = 0.0;
	}
	if(k==-1){
	  fprintf(outfile, "% 6.3f %6.3f %5.3e ", IonTotal.rhocbins[i],IonTotal.tcbins[i],IonTotal.Zcbins[i]);
	}else{
	  fprintf(outfile, "% 5.2f %5.2f %5.3e %5.3e ", Ion[k].rhocbins[i],Ion[k].tcbins[i],Ion[k].Zcbins[i],Ion[k].vcbins[i]);
	}
      /* if(k == 0) */
      /* 	if(Ion[k].vcbins[i] < 0.0 || (Ion[k].vcbins[i] > 0.0 && Ion[k].vcbins[i] < 1.e-30)){ */
      /* 	  fprintf(stderr, "STRANGE BIN: Ion[0].vcbins[%d] = %g\n", */
      /* 		  i, Ion[0].vcbins[i]); */
      /* 	} */
      }
      fprintf(outfile,"% 7.5f % 7.5f % 7.5f",IonExtra.xbins[i],IonExtra.ybins[i],IonExtra.zbins[i]);
      fprintf(outfile,"\n");
    }
    i++;
  } // while(redshift_track>redshift_begin-0.01)
  fclose(outfile); // specztau or specaim
  fprintf(stderr, "SUCCESSFUL! Done Writting....\n");
#endif // PART_BY_PART

  return 0;
}
