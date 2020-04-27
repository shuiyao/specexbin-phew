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

#define NRANSI

#define clight		2.99792458e10

#ifdef CALCQSOBKGD
int poisson_QSO(float, long *);
float gaussrand(long *);
#endif 

extern double unit_Velocity,unit_Density,unit_Time,unit_Mass;

int GetSpecParticles(){
  int FindLOS(), ShortLOS(), CheckParticleLOS();
  int i; 
  /* Real hold; */
  /* Real metalhold; */
  /* int int_hold; */
  /* float redshift_tab; */

  fprintf(stderr, "Begin: FindLOS()\n");
#ifdef SHORTSPEC
  nzbins = ShortLOS();
#else
  nzbins = FindLOS();
#endif
  // nvbins is first initialized in initions.c as floor(nzbins/(VRES/ZRES))

  //Check_Z_File(finname); /* See to open up new snapshot */
  
  for(i=0;i<nsph;i++){
    //if(i%10000==0)fprintf(stderr,"%9d %9d %5.3e % 5.3e % 5.3e % 5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n",i,count,gp[i].mass,gp[i].vel[0], gp[i].vel[1],gp[i].vel[2],gp[i].rho,gp[i].temp,gp[i].metals[0],gp[i].hsmooth,unit_Mass);
    CheckParticleLOS(i);
  }

  fprintf(stderr,"nsph= %d count= %d\n",nsph,count);

  return(count);
}

int ShortLOS()
{
  int i, nbins;
  double dz, zstep;
  float CosmicTime();
  int  cosmopar();

  cosmopar(CosmicTime(redshift_center));
  redshift = redshift_center;
  dz = BOXSIZE*hubble*unit_Velocity/clight;
  zstep = zres/(BOXSIZE*hubble*unit_Velocity/clight);
  nbins = (BOXSIZE*hubble*unit_Velocity/clight)/zres;
  zbeginline = 0;


  //bin_x = realloc(bin_x,i*sizeof(double));
  //bin_y = realloc(bin_y,i*sizeof(double));
  //bin_z = realloc(bin_z,i*sizeof(double));
  bin_x = malloc((nbins+1)*sizeof(double));
  bin_y = malloc((nbins+1)*sizeof(double));
  bin_z = malloc((nbins+1)*sizeof(double));

  bin_size = malloc((nbins+1)*sizeof(double));
  bin_coord = malloc((nbins+1)*sizeof(double));
  bin_redshift = malloc((nbins+1)*sizeof(double));
  for(i=0;i<=nbins;i++){  
    bin_size[i] = zstep;
    bin_coord[i] =  ((float)i- (float)nbins/2)/((float)nbins); // (float)i/((float)nbins);
    bin_redshift[i] = (nbins/2-i)*zres+redshift_center;

    bin_x[i] = xspec;
    bin_y[i] = yspec;
    bin_z[i] = zspec + ((float)i- (float)nbins/2)/((float)nbins);
    if(bin_z[i]>HALFBOX) bin_z[i] -= BOXSIZE;
    if(bin_z[i]<-HALFBOX) bin_z[i] += BOXSIZE;
  }
  redshift_end = bin_redshift[0];
  redshift_begin = bin_redshift[nbins-1];

  fprintf(stderr,"Straight LOS: dz= %7.5e zstep= %7.5e nbins= %d from %7.5e to %7.5e\n",dz,zstep,nbins,bin_coord[0],bin_coord[nbins-1]);
  fprintf(stderr,"Redshift_begin= %7.5e Redshift_end= %7.5e Redshift_center = %7.5e\n",redshift_begin, redshift_end, redshift_center);
  return(nbins);
}

int FindLOS()
{
  int i;
  /* int vi_hold; */
  float CosmicTime();
  int  cosmopar();

  redshift_hold = redshift;

  cosmopar(CosmicTime(redshift));
  //redshift_track = redshift;
  zstep = zres/(BOXSIZE*hubble*unit_Velocity/clight);
 
  /* 161005: {x,y,z}spec, and also phi, is defined in contspecexbin_v8.c */
  /* Option 1: Given when the program is called */
  /* Option 2: Use a random coord as the origin */
  xbegin = xspec;
  ybegin = yspec;
  zbegin = zspec;

  RotateCoords(xbegin, ybegin, zbegin, theta, phi);

  /* 161005: Now the {x,y}rotline is the projected coords of the sightline */
  xrotline = xrot;
  yrotline = yrot;
  zbeginline = zrot+zstep;

  xorig = xbegin;
  yorig = ybegin;
  zorig = zbegin;

  /* find xend, yend, and zend */
  xbreak = 0;
  ybreak = 0;
  zbreak = 0;
  i = 0;

  bin_size = malloc(sizeof(double));
  bin_coord = malloc(sizeof(double));
  bin_redshift = malloc(sizeof(double));

  /* use this loop to run line-of-site through box and make bins in 
   box coordinates (bin_coord) and redshift coordinates (bin_redshift) */
  /* vi_hold = vi; */
  while(xbreak == 0 && ybreak == 0 && zbreak == 0){
    redshift_track -= zres;
    cosmopar(CosmicTime(redshift_track));
    zstep = zres/(BOXSIZE*hubble*unit_Velocity/clight);
    zrot += zstep;
    InverseRotateCoords(xrot, yrot, zrot, theta, phi);
    i++;
    bin_size = realloc(bin_size,i*sizeof(double));
    bin_coord = realloc(bin_coord,i*sizeof(double));
    bin_redshift = realloc(bin_redshift,i*sizeof(double));
    bin_x = realloc(bin_x,i*sizeof(double));
    bin_y = realloc(bin_y,i*sizeof(double));
    bin_z = realloc(bin_z,i*sizeof(double));
    bin_size[i-1] = zstep;
    bin_coord[i-1] = zrot-zbeginline;
    bin_redshift[i-1] = redshift_track;
    bin_x[i-1] = xorig;
    bin_y[i-1] = yorig;
    bin_z[i-1] = zorig;
    if(xorig < -HALFBOX || xorig > HALFBOX) xbreak = xorig;
    if(yorig < -HALFBOX || yorig > HALFBOX) ybreak = yorig;
    if(zorig > HALFBOX) zbreak = zorig;
  }        

  zrot -= zstep;
  i--;
  redshift_track += zres;
  InverseRotateCoords(xrot, yrot, zrot, theta, phi);
  xspec = xorig;
  yspec = yorig;
  zspec = zorig;

  zlength = sqrt(pow(xbegin-xspec,2)+pow(ybegin-yspec,2)+pow(zbegin-zspec,2));

  fprintf(stderr,"i = %d\n",i);  
  fprintf(stderr,"zlength = %5.3e\n",zlength);
  fprintf(stderr,"xbegin = %9.7e, ybegin = %9.7e, zbegin = %9.7e\n",xbegin,ybegin,zbegin);
  fprintf(stderr,"xend = %9.7e, yend = %9.7e, zend = %9.7e\n",xspec,yspec,zspec);
  fprintf(stderr,"xbreak = %9.7e, ybreak = %9.7e, zbreak = %9.7e\n",xbreak,ybreak,zbreak);
  fprintf(stderr,"redshift = %9.7f redshift covered = %9.7f\n",redshift,redshift_hold-redshift);
  count = 0;
  return(i);
} 

int CheckParticleLOS(int i)
{
  double dx,dy,dr2;
  double irep[NDIM];
  /* int int_hold; */

#ifdef SHORTSPEC  
  dx = fabs(gp[i].pos[0]-xspec);
  dy = fabs(gp[i].pos[1]-yspec);
#else
  RotateCoords(gp[i].pos[0], gp[i].pos[1], gp[i].pos[2], theta, phi);
  dx = fabs(xrot-xrotline);
  dy = fabs(yrot-yrotline);
#endif

  for(irep[0] = -1; irep[0] <= 1; irep[0]++) {
    for(irep[1] = -1; irep[1] <= 1; irep[1]++) {
      for(irep[2] = -1; irep[2] <= 1; irep[2]++) {
	if(fabs(gp[i].pos[0]+irep[0]*BOXSIZE+2*gp[i].hsmooth)<HALFBOX || fabs(gp[i].pos[0]+irep[0]*BOXSIZE-2*gp[i].hsmooth)<HALFBOX){
	  if(fabs(gp[i].pos[1]+irep[1]*BOXSIZE+2*gp[i].hsmooth)<HALFBOX || fabs(gp[i].pos[1]+irep[1]*BOXSIZE-2*gp[i].hsmooth)<HALFBOX){
#ifdef SHORTSPEC
	    if(fabs(gp[i].pos[0]+irep[0]*BOXSIZE-xspec)<dx){ 
	      dx = fabs(gp[i].pos[0]+irep[0]*BOXSIZE-xspec);
	      // gp[i].pos[0] += irep[0]*BOXSIZE;
	    }
	    if(fabs(gp[i].pos[1]+irep[1]*BOXSIZE-yspec)<dy){
	      dy = fabs(gp[i].pos[1]+irep[1]*BOXSIZE-yspec);
	      // gp[i].pos[1] += irep[1]*BOXSIZE;
	    }
#else
	    if(fabs(gp[i].pos[2]+irep[2]*BOXSIZE+2*gp[i].hsmooth)<HALFBOX || fabs(gp[i].pos[2]+irep[2]*BOXSIZE-2*gp[i].hsmooth)<HALFBOX){
	      RotateCoords(gp[i].pos[0]+irep[0]*BOXSIZE, gp[i].pos[1]+irep[1]*BOXSIZE, gp[i].pos[2]+irep[2]*BOXSIZE, theta, phi);
	      if(fabs(xrot-xrotline)<dx) dx = fabs(xrot-xrotline);
	      if(fabs(yrot-yrotline)<dy) dy = fabs(yrot-yrotline);
	    }
#endif 
	  }
	}
      }
    }
  }
  
  dr2 = dx*dx+dy*dy;
  if( dr2 < NSRCHRAD*gp[i].hsmooth*NSRCHRAD*gp[i].hsmooth ){
    if(count==0){
      spec_particles = (struct spec_particle *) malloc(sizeof(struct spec_particle));
    }else{
      spec_particles = realloc(spec_particles,(count+2)*sizeof(struct spec_particle));
    }
    spec_particles[count] = gp[i];
#ifdef SHORTSPEC
    spec_particles[count].pos[0] -= xspec;
    spec_particles[count].pos[1] -= yspec;
    spec_particles[count].pos[2] -= zspec;
#endif
    count++;
  }
  return(0);
}


