#include "specexbin_defs.h"
#include "extern.h"

int InTauFact(fname,ntf,taufact)
char *fname;
int ntf;
float *taufact;
{
	int i;
	FILE *infile;

	if( (infile = fopen(fname,"r")) == NULL ) {
		fprintf(stderr,"cannot find %s\n",fname);
		exit(-1);
	}

	fprintf(stderr,"Reading taufact file %s %d\n",fname,ntf);
	for( i=0; i<ntf; i++ ) {
		if( fscanf(infile,"%g",&taufact[i]) == EOF ) {
			fprintf(stderr,"%s ended prematurely: %d\n",fname,i);
			break;
		}
	}
	fclose(infile);

	return 0;
}
