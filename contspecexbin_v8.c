/* contspecexbin"_v8 - 

   contspecexbin_v8 infile redshift_begin redshift_end Omega Lambda boxsize H_0 Omega_b flux_fac theta (x_coord) (y_coord)
*/
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
#ifdef PHEW
#include "extern.h"
#endif

#define NINTERP 1000

#ifdef PHYSSPEC
//struct metal_track_particle *spec_particles;
char *sim_id;
#else
//struct gas_particle *spec_particles;
char sim_id[20];
#endif
struct spec_particle *spec_particles;
struct spec_particle *gp;
double xspec,yspec,zspec;
int direction;
double xrot,yrot,zrot;
double aex,hubble;
int nions;
char ionbkgd[15];
ionStruct *Ion,IonTotal;
ionExtra IonExtra;
double taufact0,flux_fac;
double totMass,Lambda,boxsize,h,redshift,Omega_b,theta,zstep;
double redshift_end,redshift_hold,redshift_begin,redshift_center,redshift_track;
double hold_coord;
int NSpecTot,NSpecQuad,count;
double deltaz;
double xend, yend, zend;
double xorig, yorig, zorig, phi;
double xbegin, ybegin, zbegin;
double xbreak, ybreak, zbreak;
double xrotline, yrotline, zbeginline, zlength;
int nzbins,nvbins;
int tempcounter;
int loop,nzloopbins,nvloopbins;
double *bin_size, *bin_coord, *bin_redshift, *bin_hubble;
double *bin_x, *bin_y, *bin_z;
double *gal_field;
double redshift_tab_begin, redshift_tab_end, redshift_tab_now;
char id[10];
int vi;
int nsph;
char namesuffix[200];
double zres;
double vres;

//#ifdef OWLSFORMAT
double *carbon, *oxygen, *silicon, *iron, *temp, *rho, *hsmooth, *mass;
double *pos_tmp, *vel_tmp;

double gadget_mass, atime, h100, hubble0, box100, omega0, omegab0, lambda0;
double lunit, vunit, munit, dunit, tunit, cgs, hscale, ascale;
int mass_array[6];
unsigned int npart[6];
//#endif

double unit_Velocity;

#ifdef PAINTAVERAGEMETALS
double rhoZ, tempZ, rhoTZ[120][70][4],Z[4];
#endif

#ifdef PHEW
double RndTable[RNDTABLE];
gsl_rng *random_generator;

double get_random_number(unsigned int id)
{
  return RndTable[(id % RNDTABLE)];
}

void set_random_numbers(void)
{
  int i;

  for(i = 0; i < RNDTABLE; i++)
    RndTable[i] = gsl_rng_uniform(random_generator);
}
#endif

int main(int argc,char **argv)
{
  int i;
  double H_0;
  float CosmicTime();
  int cosmopar();
  int ContSmoothSpec(),Tau(),OutTau(),InitIons(),FreeIons(),Check_Z_File(),GetSpecParticles();
#ifdef SHORTSPEC
  double hold;
  FILE *LOSfile;
#endif
#ifdef PHYSSPEC
  char fname[100];
#endif

#if defined(WIND_BY_WIND) && defined(PHEW)
  fprintf(stderr, "ERROR: WIND_BY_WIND and PHEW are mutually exclusive!! QUIT!\n");
  exit(-1);
#endif

  zres = ZRES;
  vres = VRES;

  tempcounter = 0;

  if(argc<6){
    //fprintf(stderr,"Usage: contspecexbin_v8 infile redshift_begin redshift_end Omega Lambda boxsize H_0 Omega_b flux_factor theta ('x' coord) ('y' coord) (direction)\n");
#ifdef SHORTSPEC
    fprintf(stderr,"Usage: specexbin_short sim_infile LOSfile redshift_center boxsize flux_factor[where f is Gamma_corr=f*Gamma_input] (direction)\n");
#else
    fprintf(stderr,"Usage: contspecexbin_v8 infile redshift_begin redshift_end boxsize flux_factor[where f is Gamma_corr=f*Gamma_input] theta ('x' coord) ('y' coord) (direction)\n");
#endif
    return 0;
  }

  strcpy(ionbkgd,"HMQG");
#ifdef DOH157
  strcpy(ionbkgd,"H157");
#endif
#ifdef DOHM12
  strcpy(ionbkgd,"HM12");
#endif
#ifdef DOH12H
  strcpy(ionbkgd,"H12H");
#endif
  direction = 2;  /* Going in the z-direction unless otherwise noted */
#ifdef SHORTSPEC
  if( (LOSfile = fopen(argv[2],"r")) == NULL ) {
    fprintf(stderr,"Could not open file %s\n",argv[2]);
    return 0;
  }  
  sscanf(argv[3],"%lg",&redshift_center) ;
  sscanf(argv[4],"%lg",&boxsize) ;
  sscanf(argv[5],"%lg",&flux_fac);
  theta = 0;
  phi = 0;
  redshift_begin = redshift_center;
  redshift_end = redshift_center;
  redshift = redshift_center;
  if(argc==7) sscanf(argv[6],"%d",&direction) ;
#else
  sscanf(argv[2],"%lg",&redshift_begin) ;
  sscanf(argv[3],"%lg",&redshift_end) ;
  sscanf(argv[4],"%lg",&boxsize) ;
  sscanf(argv[5],"%lg",&flux_fac);
  sscanf(argv[6],"%lg",&theta) ;
#endif
  totMass = 0.30; //0.28; //0.238; //0.28;
  Lambda = 0.70; //0.72; //0.762; //0.72;
  Omega_b = 0.045; //0.046; //0.0418; //0.046;
  H_0 = 70; //73;
  h = 0.01*H_0;

  //cosmopar(CosmicTime(argv[3]));

  //  sscanf(argv[4],"%lg",&totMass) ;
  //  sscanf(argv[5],"%lg",&Lambda) ;
  //sscanf(argv[4],"%lg",&boxsize) ;
  //  sscanf(argv[7],"%lg",&H_0) ;
  //  sscanf(argv[8],"%lg",&Omega_b) ;

#ifdef PHEW
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42);	/* start-up seed */
  set_random_numbers();
#endif

#ifdef PHYSSPEC
  strcpy(fname,argv[1]);
  sim_id = strtok(fname,".");
  //strcpy(sim_id,"test"); 
  fprintf(stderr,"SIM_ID= %s\n",sim_id);
#endif

  if(argc>7){
    fprintf(stderr,"Using input xspec yspec zspec\n");
    sscanf(argv[7],"%lg",&xspec);
    sscanf(argv[8],"%lg",&yspec);
    sscanf(argv[9],"%lg",&zspec);
    fprintf(stderr,"Using inputted xspec= %g yspec= %g direciton= %d\n",xspec,yspec,direction);
  }else{
#ifndef SHORTSPEC
    fprintf(stderr,"Using random xspec yspec (argc=%d)\n",argc);
    i = (int) 1000000*redshift_end;
    i = (int) i*theta;
    srand48(i+100);
    i = 0;
    xspec = BOXSIZE*(drand48()-HALFBOX);
    yspec = BOXSIZE*(drand48()-HALFBOX);
    zspec = -HALFBOX;
    theta *= PI/180.0;
    phi = 2*PI*drand48();
#endif
  }
  
  
  fprintf(stderr,"xspec = %5.3e, yspec = %5.3e, zspec = %5.3e, theta = %5.3e, phi = %5.3e\n",xspec,yspec,zspec,theta,phi);

  //OpenFile(argv[1]);
  //getspecparticles(argv[1]);
  
  loop = 0;
  nzloopbins = 0;
  nvloopbins = 0;
  hold_coord = 0;
  vi = 0;
  redshift_tab_now = -1;

#ifdef SHORTSPEC					
  redshift = redshift_center;
  Check_Z_File(argv[1]); /* See to open up new snapshot */
  while(!feof(LOSfile)){
    nzloopbins = 0;
    nvloopbins = 0;
    hold_coord = 0;
    count = 0;

    fscanf(LOSfile,"%lf%lf%lf%d%s",&xspec,&yspec,&zspec,&direction,namesuffix);
    if(xspec>HALFBOX) xspec -= BOXSIZE; // New as of 10/12/12
    if(xspec<-HALFBOX) xspec += BOXSIZE;
    if(yspec>HALFBOX) yspec -= BOXSIZE;
    if(yspec<-HALFBOX) yspec += BOXSIZE;
    if(zspec>HALFBOX) zspec -= BOXSIZE;
    if(zspec<-HALFBOX) zspec += BOXSIZE; 
    if(direction==1){
      hold = zspec;
      zspec = yspec;
      yspec = xspec;
      xspec = hold;
    }
    if(direction==0){
      hold = xspec;
      xspec = yspec;
      yspec = zspec;
      zspec = hold;
    }
    GetSpecParticles();
    fprintf(stderr,"nzbins= %d, count= %d direction= %d\n",nzbins,count,direction); 
    InitIons();
    taufact0 = 1.0;
    ContSmoothSpec();
    nzloopbins += nzbins;
    nvloopbins += nvbins;
    Tau();
    OutTau();
    free(bin_x);
    free(bin_y);
    free(bin_z);
    free(bin_size);
    free(bin_coord);
    free(bin_redshift);
    FreeIons();
  }
  fclose(LOSfile);
#else // LONGSPEC
#if defined(PHEW) || defined(WIND_BY_WIND)
  nvloopbins = NVBINS_ADVANCED;
#endif  
  bin_x = malloc(sizeof(double));
  bin_y = malloc(sizeof(double));
  bin_z = malloc(sizeof(double));

  redshift = redshift_end+0.01; /* so we don't have edge effects!!! */
  redshift_track = redshift;
  while(redshift >= redshift_begin-0.01){
    //if(redshift < redshift_tab_end || redshift > redshift_tab_begin) getspecparticles(argv[1]);
    Check_Z_File(argv[1]); /* See to open up new snapshot */
    fprintf(stderr, "Begin: GetSpecParticles()\n");
    GetSpecParticles(); // nzbins is defined
    //nzbins = OutContSpectra();
    fprintf(stderr,"------------------------------------------------\n");     
    fprintf(stderr,"loop: %d nzbins = %d, count = %d, nzloopbins = %d\n",loop,nzbins,count,nzloopbins); 
    InitIons();
    // nvbins is obtained here as floor(nzbins(VRES/ZRES));
    // Also, *vbins, tcbins, etc. are realloc'ed.
    taufact0 = 1.0;
    ContSmoothSpec();
    fprintf(stderr, "ContSmoothSpec() Done.\n");

    // PhEW particles within this loop have been added to the spectra
    if(xbreak > 0) xspec -= BOXSIZE;
    if(xbreak < 0) xspec += BOXSIZE;
    if(ybreak > 0) yspec -= BOXSIZE;
    if(ybreak < 0) yspec += BOXSIZE;
    if(zbreak > 0) zspec -= BOXSIZE;
    loop++;
    nzloopbins += nzbins;
    nvloopbins += nvbins;    
    //free(spec_particles);
  }
  
  /* free(spec_particles); */

  Tau(); // nzbins, nvbins are updated as nzloopbins and nvloopbins
  OutTau();
#endif
  return 0;
}


int RotateCoords(double x, double y, double z, double theta, double phi)
{
  xrot = x*cos(phi) + y*sin(phi);
  yrot = -x*sin(phi)*cos(theta) + y*cos(phi)*cos(theta) + z*sin(theta);
  zrot = x*sin(phi)*sin(theta) - y*sin(theta)*cos(phi) + z*cos(theta);
  return 0;	
}


int InverseRotateCoords(double x, double y, double z, double theta, double phi)
{

  xorig =  x*cos(phi) - y*sin(phi)*cos(theta) + z*sin(phi)*sin(theta);
  yorig =  x*sin(phi) + y*cos(phi)*cos(theta) - z*cos(phi)*sin(theta);
  zorig =  y*sin(theta) + z*cos(theta);
  return 0;	
}

