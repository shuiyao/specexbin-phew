#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
#include <hdf5.h>
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

#define OWLSHSMOOTHCONVERT  2 //1//1.0575 //2.165

extern double unit_Mass,unit_Velocity,unit_Length,unit_Density; /* Temporary until OWLS is integrated*/


int readhdf5(char* basename, int itype){


 hid_t hdf5_file, hdf5_grp, hdf5_attribute, hdf5_dataset;
 hid_t hdf5_dataspace, hdf5_memspace;
 hid_t hdf5_holder;
 H5E_auto_t old_func;
 void *old_client_data;

 int gadget_num_files; // number of files in multifile mode
 unsigned int  gadget_npart[6];  // number of particles of each type in a single file
 int gadget_npart_tot[6]; //total number of particels of each type (all files)
 int dims[2];
 double gadget_mass[6]; // mass code; see description below
 double gadget_time; //for cosm. sim.: scale factor, else: time
 double gadget_boxsize;
 double gadget_omega0;
 double gadget_omegab0;
 double gadget_lambda0;
 double gadget_redshift;
 double gadget_h100;
 int i;
 long n;
 long ipos;
 int ifile;
 char itypebuf[20];
 sprintf(itypebuf, "/PartType%d", itype);

 char hdf5File[300];
 char hdf5base[300];
 FILE *fp;

 sprintf(hdf5base, "%s", basename);
 sprintf(hdf5File,"%s.hdf5",hdf5base);

 if( (fp=fopen(hdf5File,"r")) != NULL ){
   gadget_num_files = 1;
 }else{
   sprintf(hdf5File,"%s.0.hdf5",hdf5base);
   if( (fp=fopen(hdf5File,"r")) == NULL ){
     fprintf(stderr, "File not found: %s\n", hdf5File);
     exit(-1);
   }
 } 



 // physical parameters
 //----------------------------------------------------------------------------
 const double Mpc = 3.086e24;                         // [cm]
 const double H0 = 1e7/Mpc;                           // [s^{-1}] (100 km/s/Mpc)


 /* Open hdf5 file. */
 //---------------------------------------------------------------------

 /* Save old error handler */
 // H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
 H5Eget_auto(&old_func, &old_client_data);

 /* Turn off error handling */
 // H5Eset_auto(H5E_DEFAULT, NULL, NULL);
 H5Eset_auto(NULL, NULL);

 // hdf5_file = H5Fopen(hdf5File, H5F_ACC_RDONLY, H5P_DEFAULT);
 hdf5_file = H5Fopen(hdf5File, H5F_ACC_RDONLY, H5P_DEFAULT);

 /* Restore previous error handler */
 // H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);
 H5Eset_auto(old_func, old_client_data);


 if(hdf5_file < 0){
   fprintf(stderr,"Could not open file %s \n", hdf5File);
   exit(0);
 }

 /* Read header information */
 //---------------------------------------------------------------------
 // hdf5_grp = H5Gopen2(hdf5_file, "/Header", hdf5_holder);
 hdf5_grp = H5Gopen(hdf5_file, "/Header");

 hdf5_attribute = H5Aopen_name(hdf5_grp, "NumPart_Total");
 H5Aread(hdf5_attribute, H5T_NATIVE_UINT, gadget_npart_tot);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_grp, "NumPart_ThisFile");
 H5Aread(hdf5_attribute, H5T_NATIVE_UINT, gadget_npart);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_grp, "NumFilesPerSnapshot");
 H5Aread(hdf5_attribute, H5T_NATIVE_UINT, &gadget_num_files);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_grp, "MassTable");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gadget_mass);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_grp, "ExpansionFactor");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gadget_time);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_grp, "Redshift");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gadget_redshift);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_grp, "BoxSize");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gadget_boxsize);
 H5Aclose(hdf5_attribute);


 hdf5_attribute = H5Aopen_name(hdf5_grp, "Omega0");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gadget_omega0);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_grp, "OmegaBaryon");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gadget_omegab0);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_grp, "OmegaLambda");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gadget_lambda0);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_grp, "HubbleParam");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &gadget_h100);
 H5Aclose(hdf5_attribute);


 H5Gclose(hdf5_grp);

 totMass = gadget_lambda0;
 Lambda = gadget_omega0;
 h = gadget_h100;
 Omega_b = gadget_omegab0;


 fprintf(stderr,"NumPart_Total= %d NumPart_ThisFile= %d ExpansionFactor= %g Redshift= %g BoxSize= %g Omega0= %g OmegaBaryon= %g OmegaLambda= %g HubbleParam= %g\n",gadget_npart_tot[itype],gadget_npart[itype],gadget_time,gadget_redshift,gadget_boxsize,gadget_omega0,gadget_omegab0,gadget_lambda0,gadget_h100);

 if(gadget_npart_tot[itype] == 0){
   fprintf(stderr,"File %s contains no particles! Exiting.\n", hdf5File);
   exit(1);
 }

 for(i = 0; i < 6; i++){
     npart[i] = gadget_npart_tot[i];
     mass_array[i] = gadget_mass[i];
 }
 atime = gadget_time;
 //redshift = gadget_redshift;
 h100 = gadget_h100;
 hubble0 = H0 * h100; //[1/s]
 boxsize = gadget_boxsize;
 omega0 = gadget_omega0;
 omegab0 = gadget_omegab0;
 lambda0 = gadget_lambda0;

  totMass = gadget_omega0;
  Lambda = gadget_lambda0;
  Omega_b = gadget_omegab0;

 /* allocate memory */
 //---------------------------------------------------------------------



 gp = (struct spec_particle *) malloc(gadget_npart_tot[itype]*sizeof(struct spec_particle));
 

 pos_tmp = (double*) malloc(3 * gadget_npart[itype] * sizeof(double));
 vel_tmp = (double*) malloc(3 * gadget_npart[itype] * sizeof(double));
 rho = (double*) malloc(gadget_npart[itype] * sizeof(double));
 temp = (double*) malloc(gadget_npart[itype] * sizeof(double));
 hsmooth = (double*) malloc(gadget_npart[itype] * sizeof(double));
 mass = (double*) malloc(gadget_npart[itype] * sizeof(double));
 carbon = (double*) malloc(gadget_npart[itype] * sizeof(double));
 oxygen = (double*) malloc(gadget_npart[itype] * sizeof(double));
 silicon = (double*) malloc(gadget_npart[itype] * sizeof(double));
 iron = (double*) malloc(gadget_npart[itype] * sizeof(double));

 /* read in the particle properties */
 //---------------------------------------------------------------------
 ipos = 0;
 for(ifile = 0; ifile < gadget_num_files; ifile++){ // loop over all subfiles

   int dotpos = strlen(hdf5File) - 7;
   char prefix[200];
   strncpy(prefix, hdf5File, dotpos);
   prefix[dotpos] = '\0';
   char buf[200];
   if(gadget_num_files > 1)
     sprintf(buf, "%s.%d.hdf5", prefix, ifile);
   else
     sprintf(buf, "%s", hdf5File);

   hdf5_file = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);

   fprintf(stderr,"Reading ptype %d from file %s\n", itype, buf);

   // read number of particles
   //---------------------------------------------------------------------
   //  hdf5_grp = H5Gopen2(hdf5_file, "/Header", hdf5_holder);
   hdf5_grp = H5Gopen(hdf5_file, "/Header");
   hdf5_attribute = H5Aopen_name(hdf5_grp, "NumPart_ThisFile");
   H5Aread(hdf5_attribute, H5T_NATIVE_UINT, gadget_npart);
   H5Aclose(hdf5_attribute);
   H5Gclose(hdf5_grp);

   pos_tmp = (double*) realloc(pos_tmp, 3 * gadget_npart[itype] * sizeof(double));
   vel_tmp = (double*) realloc(vel_tmp, 3 * gadget_npart[itype] * sizeof(double));
   mass = (double*) realloc(mass, gadget_npart_tot[itype] * sizeof(double));
   rho = (double*) realloc(rho, gadget_npart_tot[itype] * sizeof(double));
   temp = (double*) realloc(temp, gadget_npart_tot[itype] * sizeof(double));
   hsmooth = (double*) realloc(hsmooth, gadget_npart_tot[itype] * sizeof(double));
   carbon = (double*) realloc(carbon, gadget_npart_tot[itype] * sizeof(double));
   oxygen = (double*) realloc(oxygen, gadget_npart_tot[itype] * sizeof(double));
   silicon = (double*) realloc(silicon, gadget_npart_tot[itype] * sizeof(double));
   iron = (double*) realloc(iron, gadget_npart_tot[itype] * sizeof(double));

   fprintf(stderr,"-> reading positions ...\n");
   //   hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
   hdf5_grp = H5Gopen(hdf5_file, itypebuf);
   //   hdf5_dataset = H5Dopen2(hdf5_grp, "Coordinates", hdf5_holder);
   hdf5_dataset = H5Dopen(hdf5_grp, "Coordinates");
   H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pos_tmp[0]);
   H5Dclose(hdf5_dataset);
   H5Gclose(hdf5_grp);
   for(i=0;i<gadget_npart[itype];i++){
     gp[ipos+i].pos[0] = pos_tmp[3*i]/gadget_boxsize-0.5;
     gp[ipos+i].pos[1] = pos_tmp[3*i+1]/gadget_boxsize-0.5;
     gp[ipos+i].pos[2] = pos_tmp[3*i+2]/gadget_boxsize-0.5;
   }


   fprintf(stderr,"-> reading velocities ...\n");
   //   hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
   hdf5_grp = H5Gopen(hdf5_file, itypebuf);
   //   hdf5_dataset = H5Dopen2(hdf5_grp, "Velocity", hdf5_holder);
   hdf5_dataset = H5Dopen(hdf5_grp, "Velocity");
   H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vel_tmp[0]);
   H5Dclose(hdf5_dataset);
   H5Gclose(hdf5_grp);
   for(i=0;i<gadget_npart[itype];i++){
     gp[ipos+i].vel[0] = vel_tmp[3*i]/sqrt(gadget_time)*1e+05/unit_Velocity;
     gp[ipos+i].vel[1] = vel_tmp[3*i+1]/sqrt(gadget_time)*1e+05/unit_Velocity;
     gp[ipos+i].vel[2] = vel_tmp[3*i+2]/sqrt(gadget_time)*1e+05/unit_Velocity;
   }


   fprintf(stderr,"-> reading masses ...\n");
   if(gadget_mass[itype] == 0){
     //     hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
     hdf5_grp = H5Gopen(hdf5_file, itypebuf);
     //     hdf5_dataset = H5Dopen2(hdf5_grp, "Mass", hdf5_holder);
     hdf5_dataset = H5Dopen(hdf5_grp, "Mass");
     H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mass[ipos]);
     H5Dclose(hdf5_dataset);
     H5Gclose(hdf5_grp);
   }
   else{
       for(i=0; i<gadget_npart[itype]; i++)
         gp[ipos+i].mass = gadget_mass[itype];
   }

   fprintf(stderr,"-> reading hsml ...\n");
   //   hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
   hdf5_grp = H5Gopen(hdf5_file, itypebuf);
   //   hdf5_dataset = H5Dopen2(hdf5_grp, "SmoothingLength", hdf5_holder);
   hdf5_dataset = H5Dopen(hdf5_grp, "SmoothingLength");
   H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &hsmooth[ipos]);
   H5Dclose(hdf5_dataset);
   H5Gclose(hdf5_grp);

   fprintf(stderr,"-> reading density ...\n");
   //   hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
   hdf5_grp = H5Gopen(hdf5_file, itypebuf);
   //   hdf5_dataset = H5Dopen2(hdf5_grp, "Density", hdf5_holder);
   hdf5_dataset = H5Dopen(hdf5_grp, "Density");
   H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rho[ipos]);
   H5Dclose(hdf5_dataset);
   H5Gclose(hdf5_grp);
   
   fprintf(stderr,"-> reading temperature ...\n");
   //   hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
   hdf5_grp = H5Gopen(hdf5_file, itypebuf);
   //   hdf5_dataset = H5Dopen2(hdf5_grp, "Temperature", hdf5_holder);
   hdf5_dataset = H5Dopen(hdf5_grp, "Temperature");
   H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &temp[ipos]);
   H5Dclose(hdf5_dataset);
   H5Gclose(hdf5_grp);

   fprintf(stderr,"-> reading carbon ...\n");
   //   hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
   hdf5_grp = H5Gopen(hdf5_file, itypebuf);
   //   hdf5_dataset = H5Dopen2(hdf5_grp, "ElementAbundance/Carbon", hdf5_holder);
   hdf5_dataset = H5Dopen(hdf5_grp, "ElementAbundance/Carbon");
   H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &carbon[ipos]);
   H5Dclose(hdf5_dataset);
   H5Gclose(hdf5_grp);

   fprintf(stderr,"-> reading oxygen ...\n");
   //   hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
   hdf5_grp = H5Gopen(hdf5_file, itypebuf);
   //   hdf5_dataset = H5Dopen2(hdf5_grp, "ElementAbundance/Oxygen", hdf5_holder);
   hdf5_dataset = H5Dopen(hdf5_grp, "ElementAbundance/Oxygen");
   H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &oxygen[ipos]);
   H5Dclose(hdf5_dataset);
   H5Gclose(hdf5_grp);

   fprintf(stderr,"-> reading silicon ...\n");
   //   hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
   hdf5_grp = H5Gopen(hdf5_file, itypebuf);
   //   hdf5_dataset = H5Dopen2(hdf5_grp, "ElementAbundance/Silicon", hdf5_holder);
   hdf5_dataset = H5Dopen(hdf5_grp, "ElementAbundance/Silicon");
   H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &silicon[ipos]);
   H5Dclose(hdf5_dataset);
   H5Gclose(hdf5_grp);
   
   fprintf(stderr,"-> reading iron ...\n");
   //   hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
   hdf5_grp = H5Gopen(hdf5_file, itypebuf);
   //   hdf5_dataset = H5Dopen2(hdf5_grp, "ElementAbundance/Iron", hdf5_holder);
   hdf5_dataset = H5Dopen(hdf5_grp, "ElementAbundance/Iron");
   H5Dread(hdf5_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &iron[ipos]);
   H5Dclose(hdf5_dataset);
   H5Gclose(hdf5_grp);
   
   for(i=0;i<gadget_npart[itype];i++){
     gp[ipos+i].mass = mass[i]*6.77e-31/(1.88e-29*pow(gadget_boxsize,3));//1.88e-29*Omega_b*pow(Mpc*gadget_boxsize/gadget_boxsize,3)/unit_Mass;  //3.94e+05;
     gp[ipos+i].rho = rho[i]*(6.77e-31/1.88e-29);    //Omega_b;
     gp[ipos+i].temp = temp[i];
     gp[ipos+i].hsmooth = hsmooth[i]/gadget_boxsize/(OWLSHSMOOTHCONVERT);
     gp[ipos+i].metals[0] = carbon[i];
     gp[ipos+i].metals[1] = oxygen[i];
     gp[ipos+i].metals[2] = silicon[i];
     gp[ipos+i].metals[3] = iron[i];
     //fprintf(stderr,"%g %g %g %g\n",gp[ipos+i].hsmooth,gp[ipos+i].rho,gp[ipos+i].temp,gp[ipos+i].metals[0]);
   }
   
   /* Close hdf5 file. */
   H5Fclose(hdf5_file);

   ipos += gadget_npart[itype];

 }
 
 free(pos_tmp);
 free(vel_tmp);
 free(mass);
 free(rho);
 free(temp);
 free(hsmooth);
 free(carbon);
 free(oxygen);
 free(silicon);
 free(iron);

 // Read unit information
 //----------------------------------------------------------------------

 if(itype == 0){
     lunit = readunits(hdf5File, itype, "Coordinates");
     vunit = readunits(hdf5File, itype, "Velocity");
     if(gadget_mass[itype] != 0){
         hdf5_file = H5Fopen(hdf5File, H5F_ACC_RDONLY, H5P_DEFAULT);
	 //         hdf5_grp = H5Gopen2(hdf5_file, "/Units", hdf5_holder);
         hdf5_grp = H5Gopen(hdf5_file, "/Units");
         hdf5_attribute = H5Aopen_name(hdf5_grp, "UnitMass_in_g");
         H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &munit);
         H5Aclose(hdf5_attribute);
         H5Gclose(hdf5_grp);
         H5Fclose(hdf5_file);
         munit /= h100;

     }
     else
         munit = readunits(hdf5File, itype, "Mass");

     dunit = readunits(hdf5File, itype, "Density");
     tunit = readunits(hdf5File, itype, "Temperature");

 }

 if(itype == 1){
     lunit = readunits(hdf5File, itype, "Coordinates");
     vunit = readunits(hdf5File, itype, "Velocity");
     if(gadget_mass[itype] != 0){
         hdf5_file = H5Fopen(hdf5File, H5F_ACC_RDONLY, H5P_DEFAULT);
	 //         hdf5_grp = H5Gopen2(hdf5_file, "/Units",hdf5_holder);
         hdf5_grp = H5Gopen(hdf5_file, "/Units");
         hdf5_attribute = H5Aopen_name(hdf5_grp, "UnitMass_in_g");
         H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &munit);
         H5Aclose(hdf5_attribute);
         H5Gclose(hdf5_grp);
         H5Fclose(hdf5_file);
         munit /= h100;
     }
     else
         munit = readunits(hdf5File, itype, "Mass");
 }

 if(itype == 4){
     lunit = readunits(hdf5File, itype, "Coordinates");
     if(gadget_mass[itype] != 0){
         hdf5_file = H5Fopen(hdf5File, H5F_ACC_RDONLY, H5P_DEFAULT);
	 //         hdf5_grp = H5Gopen2(hdf5_file, "/Units",hdf5_holder);
	 hdf5_grp = H5Gopen(hdf5_file, "/Units");
	 hdf5_attribute = H5Aopen_name(hdf5_grp, "UnitMass_in_g");
         H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &munit);
         H5Aclose(hdf5_attribute);
         H5Gclose(hdf5_grp);
         H5Fclose(hdf5_file);
         munit /= h100;
     }
     else
       munit = readunits(hdf5File, itype, "Mass");
 }
 fprintf(stderr,"munit= %g lunit= %g dunit= %g tunit= %g vunit= %g\n",munit,lunit,dunit,tunit,vunit);

 return gadget_npart_tot[itype]; //number of particles
}


double readunits(char* hdf5File, int itype, char* datasetname){


  hid_t hdf5_file, hdf5_grp, hdf5_attribute, hdf5_dataset, hdf5_holder;
 char itypebuf[20];
 sprintf(itypebuf, "/PartType%d", itype);


 hdf5_file = H5Fopen(hdf5File, H5F_ACC_RDONLY, H5P_DEFAULT);
 // hdf5_grp = H5Gopen2(hdf5_file, itypebuf, hdf5_holder);
 hdf5_grp = H5Gopen(hdf5_file, itypebuf);
 // hdf5_dataset = H5Dopen2(hdf5_grp, datasetname, hdf5_holder);
 hdf5_dataset = H5Dopen(hdf5_grp, datasetname);

 hdf5_attribute = H5Aopen_name(hdf5_dataset, "CGSConversionFactor");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &cgs);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_dataset, "h-scale-exponent");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &hscale);
 H5Aclose(hdf5_attribute);

 hdf5_attribute = H5Aopen_name(hdf5_dataset, "aexp-scale-exponent");
 H5Aread(hdf5_attribute, H5T_NATIVE_DOUBLE, &ascale);
 H5Aclose(hdf5_attribute);

 H5Dclose(hdf5_dataset);
 H5Gclose(hdf5_grp);
 H5Fclose(hdf5_file);


 return cgs * pow(atime, ascale) * pow(h100, hscale);


}
