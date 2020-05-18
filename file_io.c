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
#include "loadhdf5.h"
#endif // HDF5FORMAT
#endif // TIPSYFORMAT
#include "specexbindefs.h"
#include "proto.h"
#include "extern.h"

#define RHOLOW -9.0
#define TEMPLOW 2.5

#define NRANSI

#define clight		2.99792458e10

extern double unit_Velocity,unit_Density,unit_Time;


int Check_Z_File(char *finname){

  FILE *tabfile;
  char tabname[300],snapname[300];
  double redshift_tab;
#ifdef PAINTAVERAGEMETALS
  int iZ,jZ,k;
  double barZ;
  char line[80],rhotZname[300];
  FILE *rhotZfile;
#endif

  sprintf(tabname,"%s.tab",finname);
  if( (tabfile = fopen(tabname,"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",tabname);
    return 0;
  }  
  while(!feof(tabfile)){
    if(!fscanf(tabfile,"%s %lf %lf %lf %s", snapname, &redshift_tab, &redshift_tab_begin, &redshift_tab_end, sim_id))
      io_error_msg("file_io.c: tabfile not read.");
    fprintf(stderr,"CHECK REDSHIFTS: %g %g %g %g \n",redshift_tab,redshift_tab_begin, redshift_tab_end, redshift);
    if(redshift >= redshift_tab_begin && redshift <= redshift_tab_end){
      break;
    }
  }
  fclose(tabfile);

  if(redshift_tab != redshift_tab_now){
    if(redshift_tab_now>=0) free(gp);
    fprintf(stderr,"ATTENTION: Opening %s corresponding to redshift %5.3f\n",snapname, redshift_tab);
    Open_Snapshot(snapname);
    fprintf(stderr,"ATTENTION: Finished opening %s corresponding to redshift %5.3f\n",snapname, redshift_tab);

#ifdef PAINTAVERAGEMETALS
    sprintf(rhotZname,"%s.ionrhot",snapname);
    if( (rhotZfile = fopen(rhotZname,"r")) == NULL ) {
      fprintf(stderr,"Could not open file %s\n",rhotZname);
      return 0;
    }
    fprintf(stderr,"ATTENTION: Opened %s corresponding to redshift %5.3f\n",rhotZname, redshift_tab);
    while(fgets(line,80,rhotZfile)){
      sscanf(line,"%lg %lg %lg %lg %lg %lg %lg",&rhoZ,&tempZ,&barZ,&Z[0],&Z[1],&Z[2],&Z[3]);
      iZ = (rhoZ+1e-03-RHOLOW)*10;
      jZ = (tempZ+1e-03-TEMPLOW)*10;
      for(k=0;k<NMETALS;k++){
	rhoTZ[iZ][jZ][k] = Z[k]/barZ;
      }
    }
    
#endif

    redshift_tab_now = redshift_tab;
  }

  return 1;
}


int Open_Snapshot(char *snapname){
  int i, j;
#ifdef TIPSYFORMAT  
  Real hold;
  FILE *binfile;
  char binname[300];
  FILE *auxfile;
  char auxname[300];
  struct gas_particle *gps;
  struct aux_gas_data *auxgp;
#ifdef PHEW
  FILE *awfile;
  char awname[300];
  struct aw_gas_data *awgp;
#endif  // PHEW
#endif // TIPSYFORMAT

  int Rotate90Box();

#ifdef HDF5FORMAT
  // Do Nothing
#else
  sprintf(binname,"%s.bin",snapname);
  if( (binfile = fopen(binname,"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",binname);
    return 0;
  }
#ifdef TIPSYFORMAT 
  sprintf(auxname,"%s.aux",snapname);
  if( (auxfile = fopen(auxname,"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",auxname);
    return 0;
  }
#ifdef PHEW
  sprintf(awname,"%s.aw",snapname);
  if( (awfile = fopen(awname,"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",awname);
    return 0;
  }
#endif  // PHEW
#endif // TIPSYFORMAT
#endif // HDF5FORMAT

#ifdef HDF5FORMAT 
  double MeanWeight;

  nsph = load_hdf5(snapname, 0);
  // NOTE: read_hdf5 should automatically load global variable gheader and *P
  tipsyunits(); // Need gheader from HDF5 file
  gp = malloc(sizeof(struct spec_particle)*nsph);

  for(i=0; i < nsph ; i++){
    gp[i].mass = P[i].Mass * UNIT_M / unit_Tipsy_Mass;
    for(j=0;j<MAXDIM;j++){
      gp[i].pos[j] = P[i].Pos[j] * UNIT_L / unit_Tipsy_Length - 0.5;
      if(gp[i].pos[j] >= 0.5) gp[i].pos[j] = 0.499999;
      if(gp[i].pos[j] <= -0.5) gp[i].pos[j] = -0.499999;
      gp[i].vel[j] = P[i].Vel[j] * UNIT_V / unit_Tipsy_Velocity / sqrt(gheader.time);
#ifdef VELOCITY_UNIT_CORRECTION
      gp[i].vel[j] = gp[i].vel[j] * gheader.HubbleParam;
#endif
      // note: vel_corr = sqrt(a^3) / h; vel /= vel_corr; but when converting from hdf5 to tipsy, I ignore the *= sqrt(a^3) factor, which is canceled out.
    } // j
    gp[i].rho = P[i].Rho * UNIT_M / (pow(UNIT_L, 3) * unit_Tipsy_Density);
#ifdef DENSITY_H2_FACTOR
    gp[i].rho = gp[i].rho * gheader.HubbleParam * gheader.HubbleParam;
#endif

    MeanWeight = (1 + 4 * XHE) / (1 + P[i].Ne + XHE);
    gp[i].temp = P[i].Temp * unit_Tipsy_Temp;
    gp[i].temp *= GAMMAM1 * MeanWeight * MHYDR / KBOLTZ;
    if(i<5)
      fprintf(stdout, "gp[%d].temp = %g\n", i, gp[i].temp);

    gp[i].hsmooth = P[i].Hsml * UNIT_L / unit_Tipsy_Length * 0.5;
#ifdef QUINTIC_KERNEL
    gp[i].hsmooth = gp[i].hsmooth / 1.2275;
#endif

    for(j=0;j<NMETALS;j++) gp[i].metals[j] = P[i].metal[j];
    gp[i].sfr = P[i].Sfr;
    if(gp[i].metals[3] > 10){
      printf("BAD METAL: %d %g %lu %lu %lu\n", i, gp[i].metals[3],
	     sizeof(Real), sizeof(int), sizeof(short int));
      exit(-1);
    }
    if(i<10)
      printf("Check Format: P[%d]->ne = %g\n", i, P[i].Ne);

#ifdef PHEW
    gp[i].wind_flag = (P[i].Mcloud > 0) ? 1 : 0;
    if(gp[i].wind_flag){ // Surely a PhEW particle
      gp[i].idx = i; // Debug
      gp[i].rcloud = P[i].Rcloud;
#ifdef PHEW_RCLOUD_CORRECTION 
      gp[i].rcloud /= (gheader.time * gheader.time);
#endif	  
#ifdef PHEW_NCLOUD	  
      gp[i].ncloud = PHEW_NCLOUD;
#else
      gp[i].ncloud = 0;
#endif
#ifdef PHEW_HSMOOTH
      gp[i].hsmooth = gp[i].hsmooth / 5.04;
      // 5.04 = 128 ** (1./3.)
      // m.gad = 4./3. * rho.gad * l.gad ** 3 ?
#endif	  
    }
    else{
      gp[i].idx = -1;
      gp[i].rcloud = 0.0;
      gp[i].ncloud = 0.0;
    }
#endif    // PHEW
  } // i loop: 0 ~ nsph
  free(P);
  
#else // ----------------> Tipsy Format

#ifdef TIPSYFORMAT
  if(!fread(&header,sizeof(struct dump),1,binfile))
    io_error_msg("file_io.c: binfile not read.");
  gp = malloc(sizeof(struct spec_particle)*header.nsph);
  auxgp = (struct aux_gas_data *)malloc(sizeof(struct aux_gas_data));
  gps = (struct gas_particle *)malloc(sizeof(struct gas_particle));

#ifdef VELOCITY_UNIT_CORRECTION
  double vel_corr;
  vel_corr = sqrt(header.time * header.time * header.time) / h;
#endif

  nsph = header.nsph;
  for(i=0; i < header.nsph ; i++){
    if(!fread((char *)gps, sizeof(struct gas_particle),1,binfile))
      io_error_msg("file_io.c: gps not read.");
    gp[i].mass = gps->mass;
    for(j=0;j<MAXDIM;j++){
      gp[i].pos[j] = gps->pos[j];
#ifndef VELOCITY_UNIT_CORRECTION
      gp[i].vel[j] = gps->vel[j];
#else
      gp[i].vel[j] = gps->vel[j] / vel_corr;
#endif
    }
#ifndef DENSITY_H2_FACTOR
    gp[i].rho = gps->rho;
#else
    gp[i].rho = gps->rho * h * h;
#endif
    gp[i].temp = gps->temp;
#ifndef QUINTIC_KERNEL
    gp[i].hsmooth = gps->hsmooth * 0.5;
#else
    gp[i].hsmooth = gps->hsmooth * 0.5 / 1.2275;
#endif
    /* On 7-28-11, it became apparent I was always assuming solar abundance ratios since late 2009, because of these parentheses. -- only added this wrong way after adapted for HDF5, only affects Opp. et al 2011... phew, just in time. */ 
    if(!fread((char *)auxgp, sizeof(struct aux_gas_data),1,auxfile))
      io_error_msg("file_io.c: auxgp not read.");
    for(j=0;j<NMETALS;j++) gp[i].metals[j] = auxgp->metal[j];
    gp[i].sfr = auxgp->sfr;
    if(gp[i].metals[3] > 10){
      printf("BAD METAL: %d %g %lu %lu %lu\n", i, gp[i].metals[3],
	     sizeof(Real), sizeof(int), sizeof(short int));
      exit(-1);
    }
    //	fread(&gp[i].sfr,sizeof(Real),1,auxfile);
    if(i<10)
      printf("Check Format: auxgp[%d]->ne = %g\n", i, auxgp->ne);
#ifdef PART_BY_PART
    gp[i].delaytime = auxgp->delaytime;
#endif	

#ifdef PHEW
    awgp = (struct aw_gas_data *)malloc(sizeof(struct aw_gas_data));    
    fread((char *)awgp, sizeof(struct aw_gas_data),1,awfile);
    gp[i].wind_flag = awgp -> wind_flag;
    if(awgp -> wind_flag > 0 && auxgp -> delaytime > 0){ // Surely a PhEW particle
      gp[i].idx = i;
      gp[i].rcloud = awgp -> rcloud;
#ifdef PHEW_RCLOUD_CORRECTION // TEMPORARY, SHOULD BE REMOVED LATER
      gp[i].rcloud /= (header.time * header.time);
#endif	  
#ifndef DENSITY_H2_FACTOR
      gp[i].rho = awgp -> rho;
#else	  
      gp[i].rho = awgp -> rho * h * h; // Be careful of the units!
#endif	  
      gp[i].temp = awgp -> temp; // Be careful of the units!
#ifdef PHEW_NCLOUD	  
      gp[i].ncloud = PHEW_NCLOUD;
#else
      gp[i].ncloud = gps -> mass / awgp -> mass_cloud;
#endif
#ifdef PHEW_HSMOOTH
      gp[i].hsmooth = gps -> hsmooth * 0.5 / 5.04;
      // 5.04 = 128 ** (1./3.)
      // m.gad = 4./3. * rho.gad * l.gad ** 3 ?
#endif	  
    }
    else{
      gp[i].idx = -1;
      gp[i].rcloud = 0.0;
      gp[i].ncloud = 0.0;
    }
#endif    // PHEW
  } // i loop (tipsyformat)
  free(auxgp);
  free(gps);
#ifdef PHEW
  free(awgp);
#endif    

#endif // TIPSYFORMAT
#endif // HDF5FORMAT

  Rotate90Box();
  return 1;
}

int Rotate90Box()
{
  int i,ri;
  double hold;

  ri = 0;

  for(i=0;i<nsph;i++){
    if(direction==1){
      hold = gp[i].pos[2];
      gp[i].pos[2] = gp[i].pos[1];
      gp[i].pos[1] = gp[i].pos[0];
      gp[i].pos[0] = hold;
      
      hold = gp[i].vel[2];
      gp[i].vel[2] = gp[i].vel[1];
      gp[i].vel[1] = gp[i].vel[0];
      gp[i].vel[0] = hold;
      ri++;
    }
    
    if(direction==0){
      hold = gp[i].pos[0];
      gp[i].pos[0] = gp[i].pos[1];
      gp[i].pos[1] = gp[i].pos[2];
      gp[i].pos[2] = hold;
      
      hold = gp[i].vel[0];
      gp[i].vel[0] = gp[i].vel[1];
      gp[i].vel[1] = gp[i].vel[2];
      gp[i].vel[2] = hold;
      ri++;
    }
  }
  fprintf(stderr,"Physically rotated %d particles, because direction= %d\n", ri, direction);
  return(ri);
}

void io_error_msg(char *errmsg){
  fprintf(stderr, "Run-time Error. %s\n", errmsg);
  exit(5431);
}
