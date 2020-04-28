typedef float Real;

#define FOLDER_IONS "/scratch/shuiyao/specexbin/ionfiles/"

#define NDIM		3
#define BOXSIZE         1.0
#define HALFBOX         0.5
#define HBEXTRA         0.55
#ifdef INTKERNELNHLIMIT
#define NSRCHRAD	1.0
#else
#define NSRCHRAD	2.0
#endif

#define NSUBBIN		10
#define NBSMOOTH	5.0

#if defined(PHEW) || defined(PART_BY_PART)
#define NVBINS_ADVANCED 250 // About 200 km/s
#endif

#define NIONS 29
#define SOLAR_METALS 0.0189
#ifdef HDF5FORMAT
#define NMETALS 11
#else
#ifdef TIPSYFORMAT
#define NMETALS 4
#endif
#endif

/* redshift and velocity resolution */
//#define ZRES            6.0e-07
//#define VRES            3.0e-06
#define ZRES            3.0e-06
#define VRES            1.5e-05
//#define ZRES            3.0e-07

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
#ifdef PART_BY_PART
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
  double *redshift;
  double *binsize;
  double *bincoord;
  double *vbins;
  double *tbins;
  double *rhobins;
  double *Zbins;
#if defined(PHEW) || defined(PART_BY_PART)
  // Make sure whatever defined here are properly:
  // malloc'ed, realloc'ed, free'd
  double *vcbins;
  double *tcbins;
  double *rhocbins;
  double *Zcbins;
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
