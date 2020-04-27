#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>

typedef float Real;

#define NDIM		3
#define BOXSIZE         1.0
#define HALFBOX         0.5
#define HBEXTRA         0.55
#define NSRCHRAD	2.0

#define NSUBBIN		10

#define NMETALS 4
#define NIONS 29
#define SOLAR_METALS 0.0189


/* redshift and velocity resolution */
#define ZRES            3.0e-06
#define VRES            1.5e-05

/* Physical constants */
#define XH 	0.76
#define XHE 	(1.0 - XH) / (4.0 * XH)
#define KBOLTZ	1.3806e-16
#define MHYDR	1.6726e-24
#define GAMMAM1 (2.0/3)
#define CLIGHT	2.99792458e10
#define PI	3.14159265

/* Ion lookup table definitions */
#include "iontab.h"

struct opt_tau {
  float z;
  float rho;
  float temp;
  float metals[NIONS];
  float ions[NIONS];
} ;

struct dark_particle {
    Real mass;
    Real pos[NDIM];
    Real vel[NDIM];
    Real eps;
    Real phi ;
} ;

struct dark_particle dp[1];

struct star_particle {
    Real mass;
    Real pos[NDIM];
    Real vel[NDIM];
  //    Real metals[NMETALS] ;
    Real metals;
    Real tform ;
    Real eps;
    Real phi ;
} ;

struct star_particle sp[1];


struct dump {
  double time ;
  int nbodies ;
  int ndim ;
  int nsph ;
  int ndark ;
  int nstar ;
} ;
struct dump header ;

typedef struct ionStruct {
  char name[10];
  double *mass;
  double *vel;
  double *temp;
  double *rho;
  double *metals[NMETALS];
#ifdef MAXPARTMEASURE
  double *onepart;
#endif
  double *redshift;
  double *binsize;
  double *bincoord;
  double *vbins;
  double *tbins;
  double *rhobins;
  double *Zbins;
  float lambda,fraction,Xsec,atomwt,bsys,alpha;
  int Zcolumn;
} ionStruct;

typedef struct ionExtra {
  double *x;
  double *y;
  double *z;
  double *gal_field;
} ionExtra;

typedef struct QSO {
  float *z;
  float *r;
  float *M;
  float *R_M;
  float *L228;
} QSO;

#define ZGAP 0.237734 //2.9321736e-01 //0.195145                                                                                      
#define COEFF_1 2.51904 //2.03298 //3.03446                                                                                                           
#define COEFF_2 -0.0418026 //-0.0339018 //-0.0508473                                                                                                                   
#define COEFF_3 0.00035573 //0.000256754 //-0.000453879                                                                                                                    
#define COEFF_A 61.1189 //88.4521 //93.0455                                                                                                       
#define COEFF_B -24.2613 //-57.5255 //-41.6596                                                                       
#define COEFF_C -0.00137742 //6.89475 //3.62392                                                                                              

#define MASS_CONST 1.52408e+16
#define VEL_CONST 1166.04


