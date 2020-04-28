/* 180630 [CRITICAL] Bug Fix. Fix the mass and other properties for stars. */

#include <stdio.h>
#include "defs.h"
#include "loadhdf5.h"
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Unit Conversion for .tipsy */
/*   - Read in .hdf file, write it into .tipsy file. */
/*   - Usage: hdf52tipsy snap_* (no .hdf5) */

/* WARNING: */
/* Tipsy metallicity array contains only C, O, Si, Fe. Make sure the corresponding field is loaded from the HDF5 output. Default: [C, O, Si, Fe] = [2, 4, 7, 10] (NMETALS = 11) */

#define  CM_PER_MPC  3.085678e24

static long NumPart;

double unit_Tipsy_Time;
double unit_Tipsy_Length;
double unit_Tipsy_Density;
double unit_Tipsy_Mass;
double unit_Tipsy_Velocity;
double unit_Tipsy_Temp;

struct particle_data *P;
struct gadget_header gheader;

void tipsyunits(void)
{
  double HubbleParam, BoxSize;
  HubbleParam = gheader.HubbleParam;
  BoxSize = gheader.BoxSize;
  unit_Tipsy_Time=sqrt(8*M_PI/3)*CM_PER_MPC/(100*HubbleParam*1.e5);
  unit_Tipsy_Density=1.8791E-29*HubbleParam*HubbleParam;
  unit_Tipsy_Length=BoxSize*CM_PER_MPC*1.e-3;
  unit_Tipsy_Mass=unit_Tipsy_Density*unit_Tipsy_Length*unit_Tipsy_Length*unit_Tipsy_Length/(HubbleParam*HubbleParam);
  unit_Tipsy_Velocity=unit_Tipsy_Length/unit_Tipsy_Time;
  unit_Tipsy_Temp = pow(UNIT_V, 2);
  fprintf(stdout, "unit_Tipsy_Length = %g\n", unit_Tipsy_Length);
  fprintf(stdout, "unit_Tipsy_Mass = %g\n", unit_Tipsy_Mass);
  fprintf(stdout, "unit_Tipsy_Velocity = %g\n", unit_Tipsy_Velocity);
  // = pow(UNIT_TIPSY_L, 2) / pow((UNIT_TIPSY_L/Unit_Tipsy_V), 2); 
  return;
}

int load_hdf5(char *basename, int itype)
{
  char hdf5base[MAX_LEN_FILENAME];
  char hdf5File[MAX_LEN_FILENAME];
  FILE *fp;
  int multipart = 0;
 char itypebuf[20];
 sprintf(itypebuf, "/PartType%d", itype);

  int allocate_memory();

 sprintf(hdf5base, "%s", basename);
 sprintf(hdf5File,"%s.hdf5",hdf5base);

 if(!(fp=fopen(hdf5File,"r"))){ // More than one file
    fprintf(stdout, "%s not found.\n", hdf5File);
    sprintf(hdf5File,"%s.0.hdf5",hdf5base);
    multipart = 1;
    if(!(fp=fopen(hdf5File,"r"))){
      fprintf(stderr,"Error opening file '%s' \n",hdf5File);
      exit(0);
    }
  }
 fclose(fp);

  long i, cnt;
  int j, k;
  long noffset;
  hid_t hdf5_file, hdf5_headergrp, hdf5_attribute, hdf5_grp, hdf5_dataset;
  /* int ngas  = 1; // Why? */
  long ngas  = 0; 
  long ndark  = 0; 
  long nstar  = 0; 
  float *posvel, *single, *metals;
  int *intsingle;

  hdf5_file = H5Fopen(hdf5File, H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5_headergrp = H5Gopen1(hdf5_file, "/Header");
  
  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumFilesPerSnapshot");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &gheader.num_files);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Redshift");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.redshift);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"BoxSize");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.BoxSize);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Omega0");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.Omega0);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"OmegaLambda");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.OmegaLambda);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Time");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.time);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"HubbleParam");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gheader.HubbleParam);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumPart_Total");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, gheader.npartTotal);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumPart_ThisFile");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, gheader.npart);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"MassTable");
  H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, gheader.mass);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Metals");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_metals);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_StellarAge");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_stellarage);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Sfr");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_sfr);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Cooling");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_cooling);
  H5Aclose(hdf5_attribute);

  hdf5_attribute = H5Aopen_name(hdf5_headergrp,"Flag_Feedback");
  H5Aread(hdf5_attribute, H5T_NATIVE_INT, &gheader.flag_feedback);
  H5Aclose(hdf5_attribute);

  /* Missing several fields including flag_sfr, flag_feedback, nparttotal, etc. */

  H5Gclose(hdf5_headergrp);
  H5Fclose(hdf5_file);

  NumPart = gheader.npartTotal[0] + gheader.npartTotal[1] + gheader.npartTotal[4];

  gheader.flag_metals = NMETALSHDF5;

  fprintf(stderr,"NumFiles = %d\n", gheader.num_files);  
  fprintf(stderr,"NumPart = %ld\n", NumPart);
  fprintf(stderr,"Time = %g; Redshift = %g\n", gheader.time, gheader.redshift);  
  
  allocate_memory();

  for(k=0; k<gheader.num_files; k++){
    if(multipart)
      sprintf(hdf5File,"%s.%d.hdf5",basename,k);

    hdf5_file = H5Fopen(hdf5File, H5F_ACC_RDONLY, H5P_DEFAULT);
    if(multipart){
      hdf5_headergrp = H5Gopen1(hdf5_file, "/Header");
      hdf5_attribute = H5Aopen_name(hdf5_headergrp,"NumPart_ThisFile");
      H5Aread(hdf5_attribute, H5T_NATIVE_INT, gheader.npart);
      H5Aclose(hdf5_attribute);
      hdf5_attribute = H5Aopen_name(hdf5_headergrp,"MassTable");
      H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, gheader.mass);
      H5Aclose(hdf5_attribute);
      H5Gclose(hdf5_headergrp);
    }

    //GAS
    if(itype == 0){
    fprintf(stdout, "Reading GAS for file #%d\n", k);    
    if(!(posvel = (float *)malloc(sizeof(float)*gheader.npart[0]*3))){
      fprintf(stderr, "Failed to allocate memory for posvel\n");
      exit(-1);
    }
    if(!(single = (float *)malloc(sizeof(float)*gheader.npart[0]))){
      fprintf(stderr, "Failed to allocate memory for single\n");
      exit(-1);
    }
    if(!(metals = (float *)malloc(sizeof(float)*gheader.npart[0]*gheader.flag_metals))){
      fprintf(stderr, "Failed to allocate memory for single\n");
      exit(-1);
    }
    if(!(intsingle = (int *)malloc(sizeof(int)*gheader.npart[0]))){
      fprintf(stderr, "Failed to allocate memory for intsingle\n");
      exit(-1);
    }
    
    hdf5_grp = H5Gopen1(hdf5_file, "/PartType0");
    
    hdf5_dataset = H5Dopen1(hdf5_grp, "Coordinates");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      for(j=0; j<3; j++)
	P[i].Pos[j] = posvel[cnt*3 + j];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "Velocities");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      for(j=0; j<3; j++)
	P[i].Vel[j] = posvel[cnt*3 + j];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
    H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].ID = intsingle[cnt];
      cnt += 1;
    }

    if(gheader.mass[0] == 0 && gheader.npart[0] > 0){
      hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
	P[i].Mass = single[cnt];
	cnt += 1;
      }
    }
    else{
      for(i=ngas; i<gheader.npart[0]+ngas; i++)
	P[i].Mass = gheader.mass[0];
    }
    
    hdf5_dataset = H5Dopen1(hdf5_grp, "InternalEnergy");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].Temp = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "Density");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].Rho = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "ElectronAbundance");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].Ne = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "NeutralHydrogenAbundance");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].Nh = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "SmoothingLength");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].Hsml = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "StarFormationRate");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++) {
      P[i].Sfr = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "TemperatureMax");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].Tmax = single[cnt];
      cnt += 1;
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "DelayTime");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].DelayTime = single[cnt];
      cnt += 1;
    }

    if(H5Lexists(hdf5_grp,"FractionH2",H5P_DEFAULT)==0){
      for(i=ngas; i<gheader.npart[0]+ngas; i++){
	P[i].fH2 = 0.0;
      }
    }
    else{
      hdf5_dataset = H5Dopen1(hdf5_grp, "FractionH2");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
	P[i].fH2 = single[cnt];
	cnt += 1;
      }
    }

    hdf5_dataset = H5Dopen1(hdf5_grp, "Metallicity");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, metals);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      for(j=0; j<gheader.flag_metals; j++)
	P[i].metal[j] = metals[cnt*gheader.flag_metals + j];
      cnt += 1;
    }

#ifdef PHEW
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWKey");
    H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].Key = intsingle[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWMcloud");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].Mcloud = single[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWRcloud");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].Mcloud = single[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWWindMass");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].WindMass = single[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWLastSFTime");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].LastSFTime = single[cnt];
      cnt += 1;
    }
    hdf5_dataset = H5Dopen1(hdf5_grp, "PhEWVinit");
    H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
    H5Dclose(hdf5_dataset);
    for(i=ngas, cnt=0; i<gheader.npart[0]+ngas; i++){
      P[i].Vinit = single[cnt];
      cnt += 1;
    }
#endif    

    H5Gclose(hdf5_grp);
    free(metals);
    free(single);
    free(posvel);
    ngas += gheader.npart[0];
    }

    // DARK
    if(itype == 1){
      fprintf(stdout, "Reading DARK from file #%d\n", k);    
      noffset = gheader.npartTotal[0];
      if(!(posvel = (float *)malloc(sizeof(float)*gheader.npart[1]*3))){
	fprintf(stderr, "Failed to allocate memory for posvel\n");
	exit(-1);
      }
      if(!(single = (float *)malloc(sizeof(float)*gheader.npart[1]))){
	fprintf(stderr, "Failed to allocate memory for single\n");
	exit(-1);
      }
      if(!(intsingle = (int *)malloc(sizeof(int)*gheader.npart[1]))){
	fprintf(stderr, "Failed to allocate memory for intsingle\n");
	exit(-1);
      }
    
      hdf5_grp = H5Gopen1(hdf5_file, "/PartType1");

      hdf5_dataset = H5Dopen1(hdf5_grp, "Coordinates");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+ndark, cnt=0; i<noffset+gheader.npart[1]+ndark; i++){
	for(j=0; j<3; j++)
	  P[i].Pos[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "Velocities");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+ndark, cnt=0; i<noffset+gheader.npart[1]+ndark; i++){
	for(j=0; j<3; j++)
	  P[i].Vel[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
      H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
      H5Dclose(hdf5_dataset);
      for(i=noffset+ndark, cnt=0; i<noffset+gheader.npart[1]+ndark; i++){
	P[i].ID = intsingle[cnt];
	cnt += 1;
      }

      if(gheader.mass[1] == 0 && gheader.npart[1] > 0){
	hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
	H5Dclose(hdf5_dataset);
	for(i=noffset+ndark, cnt=0; i<noffset+gheader.npart[1]+ndark; i++){
	  P[i].Mass = single[cnt];
	  cnt += 1;
	}
      }
      else
	for(i=noffset+ndark; i<noffset+gheader.npart[1]+ndark; i++)
	  P[i].Mass = gheader.mass[1];

      H5Gclose(hdf5_grp);
      free(single);
      free(posvel);
      ndark += gheader.npart[1];
    } // flag_read_dark

    if(itype == 4){
      fprintf(stdout, "Reading STAR from file #%d\n", k);
      // STAR
      noffset = gheader.npartTotal[0] + gheader.npartTotal[1];
      if(!(posvel = (float *)malloc(sizeof(float)*gheader.npart[4]*3))){
	fprintf(stderr, "Failed to allocate memory for posvel\n");
	exit(-1);
      }
      if(!(single = (float *)malloc(sizeof(float)*gheader.npart[4]))){
	fprintf(stderr, "Failed to allocate memory for single\n");
	exit(-1);
      }
      if(!(metals = (float *)malloc(sizeof(float)*gheader.npart[4]*gheader.flag_metals))){
	fprintf(stderr, "Failed to allocate memory for single\n");
	exit(-1);
      }
      if(!(intsingle = (int *)malloc(sizeof(int)*gheader.npart[4]))){
	fprintf(stderr, "Failed to allocate memory for intsingle\n");
	exit(-1);
      }
    
      hdf5_grp = H5Gopen1(hdf5_file, "/PartType4");

      hdf5_dataset = H5Dopen1(hdf5_grp, "Coordinates");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+gheader.npart[4]+nstar; i++){
	for(j=0; j<3; j++)
	  P[i].Pos[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "Velocities");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, posvel);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+gheader.npart[4]+nstar; i++){
	for(j=0; j<3; j++)
	  P[i].Vel[j] = posvel[cnt*3 + j];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "ParticleIDs");
      H5Dread(hdf5_dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, intsingle);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+gheader.npart[4]+nstar; i++){
	P[i].ID = intsingle[cnt];
	cnt += 1;
      }

      if(gheader.mass[4] == 0 && gheader.npart[4] > 0){
	hdf5_dataset = H5Dopen1(hdf5_grp, "Masses");
	H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
	H5Dclose(hdf5_dataset);
	for(i=noffset+nstar, cnt=0; i<noffset+gheader.npart[4]+nstar; i++){
	  P[i].Mass = single[cnt];
	  cnt += 1;
	}
      }
      else
	for(i=noffset+nstar; i<noffset+gheader.npart[4]+nstar; i++)
	  P[i].Mass = gheader.mass[4];

      hdf5_dataset = H5Dopen1(hdf5_grp, "StellarFormationTime");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+gheader.npart[4]+nstar; i++){
	P[i].Sfr = single[cnt];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "TemperatureMax");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, single);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+gheader.npart[4]+nstar; i++){
	P[i].Tmax = single[cnt];
	cnt += 1;
      }

      hdf5_dataset = H5Dopen1(hdf5_grp, "Metallicity");
      H5Dread(hdf5_dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, metals);
      H5Dclose(hdf5_dataset);
      for(i=noffset+nstar, cnt=0; i<noffset+gheader.npart[4]+nstar; i++){
	for(j=0; j<gheader.flag_metals; j++)
	  P[i].metal[j] = metals[cnt*gheader.flag_metals + j];
	cnt += 1;
      }
    
      H5Gclose(hdf5_grp);
      free(posvel);
      free(single);
      free(metals);
      free(intsingle);
      nstar += gheader.npart[4];
    } // flag_read_star

    H5Fclose(hdf5_file);
    fprintf(stderr, "File: %d ngas = %d(%5.3f)\n", k, gheader.npart[0],
	    (float)(gheader.npart[0])/(float)(NumPart));
    fprintf(stderr, "File: %d ndark = %d(%5.3f)\n", k, gheader.npart[1],
	    (float)(gheader.npart[1])/(float)(NumPart));
    fprintf(stderr, "File: %d nstar = %d(%5.3f)\n", k, gheader.npart[4],
	    (float)(gheader.npart[4])/(float)(NumPart));
  }
  return gheader.npart[itype];
}

int allocate_memory(void)
{
  fprintf(stdout, "Allocating %6.3f GB Memory for %ld particles.\n",
	  NumPart * sizeof(struct particle_data) / (1024. * 1024. * 1024.),
	  NumPart);

  if(!(P=malloc(NumPart * sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }

  fprintf(stdout, "Memory allocated.\n");
  /* P--;    */
  /* start with offset 1 */
  return 0;
}
