#define NDIM		3
#define BOXSIZE         1.0
#define HALFBOX         0.5
#define HBEXTRA         0.55
#define NSRCHRAD	2.0
#define NBSMOOTH	5.0

#if defined(PHEW) || defined(WIND_BY_WIND)
#define NVBINS_ADVANCED 250 // About 200 km/s
#endif

#define NIONS 29
#define SOLAR_METALS 0.0189

/* redshift and velocity resolution */
//#define ZRES            6.0e-07
//#define VRES            3.0e-06
#define ZRES            3.0e-06
#define VRES            1.5e-05
//#define ZRES            3.0e-07


/* Physical constants */
#define XH 	0.76
#define KBOLTZ	1.381e-16
#define MHYDR	1.673e-24
#define CLIGHT	2.99792458e10
#define PI	3.14159265

#define NMETALS 4

/* Ion lookup table definitions */
#include "iontab.h"

struct spec_particle {
  int ID;
  Real mass;
  Real rho;
  Real temp;
  Real metals[NMETALS];
  Real hsmooth;
  Real pos[NDIM];
  Real vel[NDIM];
  Real sfr;
  Real vlaunch;
  Real ageaway;
  int nrec;
  Real mgal_launch;
  Real mgal;
  Real dtravel;
  Real dgal;
  Real dpeculiar;
  Real vrel;
#ifdef WIND_BY_WIND
  Real delaytime;
#endif  
#ifdef PHEW
  int idx;
  int wind_flag;
  Real rcloud;
  Real ncloud;
#endif  
};

struct spec_particle hp[1];

typedef struct ionStruct {
  char name[10];
  double *mass;
  double *vel;
  double *temp;
  double *rho;
  double *metals[NMETALS];
#ifdef PHYSSPEC
  double *sfr;
  double *wtmass;
  double *mgal;
  double *dgal;
  double *age;
  double *nrec;
  double *vlaunch;
#endif
  double *redshift;
  double *binsize;
  double *bincoord;
  double *vbins;
  double *tbins;
  double *rhobins;
  double *Zbins;
#if defined(PHEW) || defined(WIND_BY_WIND)
  // Make sure whatever defined here are properly:
  // malloc'ed, realloc'ed, free'd
  double *vcbins;
  double *tcbins;
  double *rhocbins;
  double *Zcbins;
#endif
#ifdef PHYSSPEC
  double *sfrbins;
  double *mgalbins;
  double *dgalbins;
  double *agebins;
  double *nrecbins;
  double *vlaunchbins;
#endif
  float lambda,fraction,Xsec,atomwt,bsys,alpha;
  int Zcolumn;
} ionStruct;

typedef struct ionExtra {
  double *redshift;
  double *x;
  double *y;
  double *z;
  double *xbins;
  double *ybins;
  double *zbins;
  double *gal_field;
} ionExtra;
