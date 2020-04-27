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

static int startflag=1;

double unit_Time,unit_Density,unit_Length,unit_Mass,unit_Velocity,unit_Temperature,unit_DUDT;
double H0,t0,aex3,aexhub;

static double etaold;

int cosmounits()
{
	float L;
    double Pi=3.14159265358979323846;
    double km=1.E5;
    double Mpc=3.086E24;
    double m_p=1.6726231E-24;     /* proton mass */
    double k_B=1.380622E-16;      /* Boltzman constant */

	boxsize /= h;
	L = boxsize;
	if( Lambda > 0.01 && totMass < 1.0 ) {
		if( Lambda != 1.-totMass ) {
/*			fprintf(stderr,"Setting Lambda = %g\n",1.-totMass);*/
			Lambda = 1.-totMass;
		}
	}
	else Lambda = 0.0;

    H0=sqrt(8*Pi/3);
	t0 = 2./(3*H0);

    unit_Time=H0*Mpc/(100*h*km);
    unit_Density=1.8791E-29*h*h;
    unit_Length=L*Mpc;
    unit_Mass=unit_Density*unit_Length*unit_Length*unit_Length;
    unit_Velocity=unit_Length/unit_Time;
    unit_Temperature=m_p*unit_Velocity*unit_Velocity/k_B;
    unit_DUDT=unit_Density*unit_Velocity*unit_Velocity/unit_Time;
	fprintf(stderr,"COSMO PARAMS:  L=%g Mpc, h=%g, Omega=%g, Lambda=%g Omega_b=%g\n",L,h,totMass,Lambda,Omega_b);
	fprintf(stderr,"UNITS: T=%g rho=%g L=%g M=%g v=%g\n",unit_Time,unit_Density,unit_Length,unit_Mass,unit_Velocity);

	return 0;
}

int cosmopar(t)
float t;
{
	float t1;
	double tol=1.e-6;
	double a0,astar,eta,etalast;
	double f,fprime;
	int it;

	if( startflag ) {
		cosmounits();
		startflag = 0;
		etaold = 1.;
	}
	if( fabs(totMass - 1.0) < tol ) {
		t1 = t/t0;
		aex3 = t1*t1;
		aex = pow((aex3),1./3.);
		hubble = 2.0/3.0/t;
		aexhub = aex*hubble;
		redshift = 1./aex - 1.;
	}
	else if( totMass < 1.0 ) {
		if( Lambda > 0.0 ) {	/* Flat, low-density universe */
                eta = sqrt(1.-totMass)*1.5*H0*t;
                aex = pow(sqrt(totMass/(1.-totMass))*sinh(eta),2./3);
		aex3 = aex*aex*aex;
		hubble = H0*sqrt(totMass/aex3+Lambda);
		aexhub = aex*hubble;
		redshift = 1./aex - 1.;
		}
		else {			/* Open universe */
        	a0=1./H0/sqrt(1.-totMass);
        	astar=.5*H0*H0*totMass*a0*a0*a0;
        	it=0;
        	eta=etaold;
		do {
        		f=astar*(sinh(eta)-eta)-t;
        		fprime=astar*(cosh(eta)-1.);
			etalast=eta;
        		eta=eta-f/fprime;
        		if( (it++) > 20 ) {
				fprintf(stderr,"Overiterated in cosmopar %d %g %g\n",it,eta,etalast);
				break;
			}
        	} while( fabs(eta-etalast)/etalast > tol );

		aex = astar*(cosh(eta)-1.)/a0;
		aex3 = aex*aex*aex;
		etaold = eta;
		redshift = 1./aex - 1.;
		hubble=H0*(1.+redshift)*sqrt(1.+totMass*redshift);
		aexhub = aex*hubble;
		}
	}
	else if( totMass > 1.0 ) {
        	a0=1./H0/sqrt(totMass-1.);
        	astar=.5*H0*H0*totMass*a0*a0*a0;
        	it=0;
        	eta=etaold;
		do {
        		f=astar*(eta-sin(eta))-t;
        		fprime=astar*(1.-cos(eta));
			etalast=eta;
        		eta=eta-f/fprime;
        		if( (it++) > 20 ) {
				fprintf(stderr,"Overiterated in cosmopar %d %g %g\n",it,eta,etalast);
				break;
			}
        	} while( fabs(eta-etalast)/etalast > tol );

		aex = astar*(1.-cos(eta))/a0;
		aex3 = aex*aex*aex;
		etaold = eta;
		redshift = 1./aex - 1.;
		hubble=H0*(1.+redshift)*sqrt(1.+totMass*redshift);
		aexhub = aex*hubble;
	}

	return 0;
}

float CosmicTime(z)	/* returns system time at redshift z */
float z;
{
	double tol=1.e-6;
	double a0,astar,eta,etalast,t=1.;
	double f,fprime,aextemp;
	int it;

	if( startflag ) {
		cosmounits();
		startflag = 0;
		etaold = 1.;
	}
	if( fabs(totMass - 1.0) < tol ) t = t0*pow(1.+z,-1.5);
	else if( totMass < 1.0 ) {
		if( Lambda > 0.0 ) {
		it = 0;
		aextemp = 1./(1.+z);
		t = etaold;
		do {
			f = sqrt(Lambda/totMass)*pow(aextemp,1.5) +
			    sqrt(aextemp*aextemp*aextemp*Lambda/totMass+1) -
			    exp(1.5*sqrt(Lambda)*H0*t);
			fprime =-1.5*sqrt(Lambda)*H0*exp(1.5*sqrt(Lambda)*H0*t);
			etalast = t;
			t = t-f/fprime;
        		if( (it++) > 20 ) break;
        	} while( fabs(t-etalast)/etalast > tol );
		etaold = t;
		}
		else {
        	a0=1./H0/sqrt(1.-totMass);
        	astar=.5*H0*H0*totMass*a0*a0*a0;
		aextemp = 1./(1.+z);
        	it=0;
        	eta=etaold;
		do {
        		f=astar*(cosh(eta)-1.)/a0-aextemp;
        		fprime=sinh(eta)*astar/a0;
			etalast=eta;
        		eta=eta-f/fprime;
        		if( (it++) > 20 ) {
				fprintf(stderr,"Overiterated in CosmicTime %d %g %g\n",it,eta,etalast);
				break;
			}
        	} while( fabs(eta-etalast)/etalast > tol );

		t = astar*(sinh(eta)-eta);
		etaold = eta;
		}
	}
	else if( totMass > 1.0 ) {
        	a0=1./H0/sqrt(totMass-1.);
        	astar=.5*H0*H0*totMass*a0*a0*a0;
		aextemp = 1./(1.+z);
        	it=0;
        	eta=etaold;
		do {
        		f=astar*(1.-cos(eta))/a0-aextemp;
        		fprime=sin(eta)*astar/a0;
			etalast=eta;
        		eta=eta-f/fprime;
        		if( (it++) > 20 ) {
				fprintf(stderr,"Overiterated in CosmicTime %d %g %g\n",it,eta,etalast);
				break;
			}
        	} while( fabs(eta-etalast)/etalast > tol );

		t = astar*(eta-sin(eta));
		etaold = eta;
	}

	return t;
}

