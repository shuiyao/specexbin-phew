#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "malloc.h"
#include <math.h>

#define MAXIONS		50
#define MAXPIX		135000
#define NOISEFACT	1.0	// Factor to multiply noise by

#define ZMETALINIT	0
#define SQRT2PI		2.506628275
#define SQRTPI		1.7724539  
#define ckms		2.99792458e5

#define SOLARMET	0.0122


int main(int argc,char **argv)
{
	int i,j,nl;
	double *z,*rho,*T,*Zmet;
	double *tau,*lambda;
	double *flux,*vel,*vin;
	double *intensity, *noise;
	double *rhobin,*Tbin,*Zbin;
	double dvres,norm;
	double prob,select,deltaint;
	double taufact0; /* overall tau multiplication factor */
	double COS_G130M_LSF[101][10];	// pix=-50,50; lam=1150,1450,50
	double COS_G160M_LSF[101][10];	// pix=-50,50; lam=1450,1750,50
        double lambda0[MAXIONS],ftrans[MAXIONS],Zover[MAXIONS];
	int ncos_G130M,ncos_G160M;
	int npix,ns,nions,nwidth,ionnum;
	double zbegin,zend,dz,zquasar;
	double lion,maxflux,junk;
	double a[100];
	char ion[80],fname[80],line[2048];
	int ilo,ihi;
	double df,dv;
	double SN;
	FILE *outspec,*inions=NULL,*inspec,*coslsf;
	int scanline();


	z = malloc(MAXPIX*sizeof(double));
        rho = malloc(MAXPIX*sizeof(double));
        T = malloc(MAXPIX*sizeof(double));
        Zmet = malloc(MAXPIX*sizeof(double));
        tau = malloc(MAXPIX*sizeof(double));
        flux = malloc(MAXPIX*sizeof(double));
        lambda = malloc(MAXPIX*sizeof(double));
        vel = malloc(MAXPIX*sizeof(double));
        vin = malloc(MAXPIX*sizeof(double));
        intensity = malloc(MAXPIX*sizeof(double));
        noise = malloc(MAXPIX*sizeof(double));
        rhobin = malloc(MAXPIX*sizeof(double));
        Tbin = malloc(MAXPIX*sizeof(double));
        Zbin = malloc(MAXPIX*sizeof(double));

	if( argc != 8 ) {
		fprintf(stderr,"usage: mkspec infile S/N redshift taufact dvpix dvres ion\n");
		exit(-1);
	}
	if( (inspec=fopen(argv[1],"r")) == NULL ) {
		fprintf(stderr,"Could not open %s\n",argv[1]);
		exit(-1);
	}
	SN = atof(argv[2]);
	zquasar = atof(argv[3]);
	taufact0 = atof(argv[4]);
	dv = atof(argv[5]);
	dvres = atof(argv[6]);
	strcpy(ion,argv[7]);
	if( taufact0 <= 0 ) {
		fprintf(stdout,"taufact input nonsense: setting taufact=1\n");
		taufact0 = 1.0;
	}

/* read in COS LSF */
	if( (coslsf=fopen("./COS_G130M_LSF.dat","r")) == NULL ) {
		fprintf(stderr,"Could not open COS_G130M_LSF.dat\n");
		exit(-1);
	}
	fgets(line,2048,coslsf);
	nl=0; while( fgets(line,2048,coslsf) != NULL ) {
		sscanf(line,"%lg %lg %lg %lg %lg %lg %lg %lg",&junk,&COS_G130M_LSF[nl][0],&COS_G130M_LSF[nl][1],&COS_G130M_LSF[nl][2],&COS_G130M_LSF[nl][3],&COS_G130M_LSF[nl][4],&COS_G130M_LSF[nl][5],&COS_G130M_LSF[nl][6]);
		nl++;
	}
	ncos_G130M=nl;
	fclose(coslsf);

/* input ion info */
#ifdef DO9IONS
	fprintf(stdout, "-------- DO 9 IONS! --------\n");
        if( (inions=fopen("./specions_i9.dat","r")) == NULL ) {
                fprintf(stderr,"specions_i9.dat not found :-(\n");
                exit(-1);
        }
#else
        if( (inions=fopen("./specions_i31.dat","r")) == NULL ) {
                fprintf(stderr,"specions_i31.dat not found :-(\n");
                exit(-1);
        }
#endif
        i = 0;
	ionnum = 0;
        while ( fgets(line,2048,inions) != NULL && i < MAXIONS ) {
                sscanf(line,"%s %lg %lg %lg %lg %lg",fname,&lambda0[i],&ftrans[i],&df,&df,&Zover[i]);
                //fprintf(stderr,"%s %s %lg %lg %lg %lg %lg\n",fname,ion,lambda0[i],ftrans[i],df,df,Zover[i]);
                if( strcmp(ion,fname) == 0 ) ionnum = i;
                i++;
        }
	//	if( strstr(ion,"HI") == NULL ) taufact0 = 1.0; // no tau adjustment if not HI
        fprintf(stdout,"Doing ion %s %d lam= %g f= %g Zover= %g ftau= %g\n",ion,ionnum,lambda0[ionnum],ftrans[ionnum],Zover[ionnum],taufact0);
        fclose(inions);
	
/* read in simulated optical depths versus redshift */
	nl = 0;
        while ( fgets(line,2048,inspec) != NULL && nl < MAXPIX ) {
                nions = scanline(line,a)/4;
                if( nions > MAXIONS ) nions = MAXIONS;
		z[nl] = a[0];
                rho[nl] = a[4*(ionnum+2)-4];
                T[nl] = a[4*(ionnum+2)-3];
                Zmet[nl] = a[4*(ionnum+2)-2];
		if(Zmet[nl]<0) Zmet[nl] = 0;
                tau[nl] = a[4*(ionnum+2)-1];
		if(nl < 10) fprintf(stdout, "tau[%d] = %g\n", nl, tau[nl]);
	//if( Zmet[nl]>0 ) fprintf(stderr,"%d %g %g %g %g %g\n",nl,z[nl],rho[nl],T[nl],Zmet[nl],tau[nl]);
		nl++;
        }

/*	nl=0; while( fgets(line,2048,inspec) != NULL ) {
		sscanf(line,"%lg %lg %lg %lg %lg %lg %lg %lg",&z[nl],&junk,&junk,&junk,&rho[nl],&T[nl],&Zmet[nl],&tau[nl]);
		nl++;
	}*/
	fclose(inspec);
	npix = nl;
	fprintf(stdout,"zq=%g npix=%d zbegin=%g zend=%g\n",zquasar,npix,z[0],z[npix-1]);

/* optical depth to flux */
	vin[0] = 0.;
	for( i=0; i<npix; i++ ) {
		//if( i>0 ) vin[i] = (lambda[i]-lambda[i-1])*ckms/lambda0+vin[i-1];
		//vin[i] = ckms*((1+z[i])*(1+z[i])-1)/((1+z[i])*(1+z[i])+1);
		vin[i] = ckms*log(1+z[i]);
		intensity[i] = exp(-tau[i]*taufact0);
	}
	lion=0; for( i=0; i<npix; i++ ) lion += intensity[i];
	fprintf(stdout,"<Int>= %g (npixin=%d) %g %g\n",lion/npix,npix,vin[0],vin[npix-1]);

/* Smooth spectrum to desired velocity resolution */
	vel[0] = vin[npix-1];
	if( dvres > 0 ) {
	    dvres /= 2.35482;	// FWHM to sigma
	    for( ns=0; 1; ns++ ) {
		vel[ns] = vin[npix-1]+dv*ns;
		for( i=0; i<npix; i++ ) if( fabs(vin[i]-vel[ns]) < 3*dvres ) break;
		ilo = i-1; if( ilo<0 ) ilo=0;
		for( i=npix-1; i>=0; i-- ) if( fabs(vel[ns]-vin[i]) < 3*dvres ) break;
		ihi = i+1; if( ihi>npix-1 ) ihi=npix-1;
		norm = flux[ns] = rhobin[ns] = Tbin[ns] = Zbin[ns] = 0.;
		for( j=ilo; j<=ihi; j++ ) {
			df = exp(-0.5*(vin[j]-vel[ns])*(vin[j]-vel[ns])/(dvres*dvres));
			flux[ns] += intensity[j]*df;
			rhobin[ns] += pow(10.,rho[j])*df;
			Tbin[ns] += pow(10.,T[j])*df;
			Zbin[ns] += Zmet[j]*df;
			norm += df;
		}
		flux[ns] /= (norm+1.e-30);
		rhobin[ns] /= (norm+1.e-30);
		Tbin[ns] /= (norm+1.e-30);
		Zbin[ns] /= (norm+1.e-30);
		if( vel[ns]+dv > vin[0] ) break;
	    }
	//for( i=0; i<10; i++ ) fprintf(stdout,"=1= %d %g %g\n",i,vel[i],flux[i]);
	}
	else if (dvres==0 ) {	// Use COS LSF
	    fprintf(stderr,"Using COS LSF ... **IF YOU ARE SEAN YOU SHOULDN'T BE SEEING THIS LINE!!!!*****\n");
	    for( ns=0; 1; ns++ ) {

		vel[ns] = vin[npix-1]+dv*ns;
		for( i=0; i<npix; i++ ) if( vin[i] > vel[ns] && vin[i+1] < vel[ns] ) nl = i;
		ilo = nl-50; if( ilo < 0 ) ilo=0;
		ihi = nl+50; if( ihi > npix-1 ) ihi=npix-1;
		norm = flux[ns] = rhobin[ns] = Tbin[ns] = Zbin[ns] = 0.;
		for( j=ilo; j<=ihi; j++ ) {
			df = COS_G130M_LSF[j-nl+50][6];  // use G130M, 1450A
			flux[ns] += intensity[j]*df;
			rhobin[ns] += pow(10.,rho[j])*df;
			Tbin[ns] += pow(10.,T[j])*df;
			Zbin[ns] += Zmet[j]*df;
			norm += df;
		}
		flux[ns] /= (norm+1.e-30);
		rhobin[ns] /= (norm+1.e-30);
		Tbin[ns] /= (norm+1.e-30);
		Zbin[ns] /= (norm+1.e-30);
//		if( ns < 10 ) fprintf(stderr,"=+= %d %d %d %d %g %g %g %g\n",ns,nl,ilo,ihi,vel[ns],vin[nl],flux[ns],rhobin[ns]);
		if( vel[ns]+dv > vin[0] ) break;
	    }
	}
	else {
		fprintf(stderr,"dvres %g must be >=0; if 0, uses COS LSF\n",dvres);
		exit(-1);
	}

	lion=0; for( i=0; i<ns; i++ ) lion += flux[i];
	fprintf(stdout,"<Int>= %g (npixout=%d)\n",lion/ns,ns);

/* Add noise to spectrum and output */
	sprintf(fname,"%s.raw",ion);
	outspec = fopen(fname,"w");
	for(i=-20;i<0;i++){
	  fprintf(outspec,"%12.4f %12.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",lambda0[ionnum]*(exp((vel[0]+(i*dv))/ckms)),vel[0]+(i*dv),1.0,1./SN,log10(rhobin[0]),log10(Tbin[0]),log10(Zbin[0]+1.e-20));
	}
	for( i=0,df=0; i<ns; i++ ) {
	//	if(vel[i]>=0 ) continue;
		lion = lambda0[ionnum]*exp(vel[i]/ckms);
		if( lion/lambda0[ionnum]-1 > zquasar ) continue;
		noise[i] = 1./SN;
		do {
			deltaint = (drand48()-0.5)*10.*noise[i];
			prob = 1./(noise[i]*SQRT2PI) * exp( -0.5*(deltaint/noise[i])*(deltaint/noise[i]) );
			select = drand48()/(noise[i]*SQRT2PI);
		} while ( prob < select );
#ifndef NONOISE		
		flux[i] += deltaint;
#endif
		df += deltaint;
		fprintf(outspec,"%12.4f %12.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",lion,vel[i],flux[i],noise[i],log10(rhobin[i]),log10(Tbin[i]),log10(Zbin[i]+1.e-20));
	}
	for(i=1;i<=20;i++){
	  fprintf(outspec,"%12.4f %12.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",lambda0[ionnum]*(exp((vel[ns-1]+(i*dv))/ckms)),vel[ns-1]+(i*dv),1.0,noise[ns-1],log10(rhobin[ns-1]),log10(Tbin[ns-1]),log10(Zbin[ns-1]+1.e-20));
	}
	lion=0; for( i=0; i<ns; i++ ) lion += flux[i];
	fprintf(stdout,"<Int>= %g <Noise>= %g\n",lion/ns,df/ns);
	fclose(outspec);

	exit(0);
}

int scanline(line,a)
char *line;
double *a;
{
	int i,ip;
	int count;
	char word[80];

	count = 0;
	i=0; while( line[i++] == ' ' ); if( i>0 ) i--;
	ip = i;
	for(; i<strlen(line); i++ ) {
		if( line[i] == '\n' ) break;
		word[i-ip] = line[i];
		if( line[i] == ' ' ) {
			word[i-ip] = '\0';
			a[count] = atof(word);
                        if( strlen(word) != 0 ) count++;
			ip = i+1;
		}
	}
	return count;
}

