#ifdef PHEW
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
/* The Random Number Generator */
/* SH161008: I grabbed it from the Gadget-3 */
#define  RNDTABLE 32768
extern double RndTable[RNDTABLE];
extern gsl_rng *random_generator;	/*!< the random number generator used */

extern double phew_mcinit; 	/* Initialized in cmd */
#endif // PHEW

//extern struct gas_particle *spec_particles;
char sim_id[20];
extern struct spec_particle *gp;
extern struct spec_particle *spec_particles;

extern double xspec,yspec,zspec;
extern int direction;
extern double xrot,yrot,zrot;
extern double aex,hubble;
extern char ionbkgd[15];
extern int nions;
extern ionStruct *Ion,IonTotal;
extern ionExtra IonExtra;
extern double taufact0,flux_fac;
extern double totMass,Lambda,boxsize,h,redshift,Omega_b,theta,zstep;
extern double redshift_end,redshift_hold,redshift_begin,redshift_center,redshift_track;
extern int NSpecTot,NSpecQuad,count;
extern double hold_coord;
extern double deltaz;
extern double xend, yend, zend;
extern double xorig, yorig, zorig, phi;
extern double xbegin, ybegin, zbegin;
extern double xbreak, ybreak, zbreak;
extern double xrotline, yrotline, zbeginline, zlength;
extern int nzbins,nvbins;
extern int tempcounter;
extern int loop,nzloopbins,nvloopbins;
extern double *bin_size, *bin_coord, *bin_redshift, *bin_hubble;
extern double *bin_x, *bin_y, *bin_z;
extern double *gal_field;
extern double redshift_tab_begin, redshift_tab_end, redshift_tab_now;
extern char id[10];
extern int vi;
extern int nsph;
extern char namesuffix[200];
extern double zres, vres;

//#ifdef HDF5FORMAT
extern double *carbon, *oxygen, *silicon, *iron, *temp, *rho, *hsmooth, *mass;
extern double *pos_tmp, *vel_tmp;

extern double gadget_mass, atime, h100, hubble0, box100, omega0, omegab0, lambda0;
extern double lunit, vunit, munit, dunit, tunit, cgs, hscale, ascale;
extern int mass_array[6];
extern unsigned int npart[6];
//#endif

#ifdef PAINTAVERAGEMETALS
extern double rhoZ, tempZ, rhoTZ[120][70][4],Z[4];
#endif
