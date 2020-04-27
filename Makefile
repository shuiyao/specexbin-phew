OPT     +=  -DTIPSYN2FORMAT
OPT     +=  -DDENSITY_H2_FACTOR
OPT	+=  -DOUTPUT_LOCAL_FOLDER
#OPT	+=  -DQUINTIC_KERNEL
OPT	+=  -DVELOCITY_UNIT_CORRECTION

OPT	+=  -DWIND_BY_WIND # Mutually exclusive with PHEW
OPT	+=  -DSINGLE_VOFFSET_PER_PARTICLE # Not recommended. The Pygad way.

#OPT	+=  -DPHEW # Mutually exclusive with WIND_BY_WIND
#OPT	+=  -DPHEW_MOREINFO
#OPT	+=  -DPHEW_VERBOSE=1
#OPT	+=  -DPHEW_HSMOOTH
#OPT	+=  -DPHEW_RCLOUD_CORRECTION # Should be taken out then
#OPT	+=  -DPHEW_NCLOUD=1000.0 # Should be taken out in the end

#OPT	+=  -DNO_WIND_NGB_STAT
#OPT    +=  -DVARYGALBKG
#OPT    +=  -DHOTWINDS
#OPT    +=  -DCALCQSOBKGD
#OPT	+=  -DZEROVEL
#OPT	+=  -DDOLYAONLY
#OPT    +=  -DDOHANDHEONLY
#OPT    +=  -DDOC4ONLY
#OPT	+=  -DDOO6ONLY
#OPT     +=  -DDO5IONS
#OPT     +=  -DDO6IONS
OPT     +=  -DDO9IONS
#OPT     +=  -DDOHIZIONS
#OPT     +=  -DDOXRAYIONS
#OPT     +=  -DDOH157
OPT     +=  -DDOHM12
#OPT     +=  -DDOH12H
#OPT    +=  -DPIPELINE
#OPT     +=  -DNEUTRALGAS
#OPT     +=  -DCALCDENSFIELD
#OPT	 +=  -DNONEQUIL
#OPT	 +=  -DBTURB
#OPT     +=  -DSMOOTHSPH
#OPT     +=  -DPHYSSPEC
#OPT     +=  -DVARGALBKGD
#OPT	+=  -DMETALFLOOR
OPT     +=  -DSHORTSPEC
#OPT     +=  -DTEMPOVER105
#OPT	+=  -DTEMPOVER1052
#OPT     +=  -DNHLIMIT

OPT     +=  -DINTKERNELNHLIMIT

#OPT     +=  -DHALFSMOOTH
#OPT     +=  -DPAINTAVERAGEMETALS=2.0
#OPT	+=  -DPAINTAVEMETALLICITY=0.009574

# OBSOLETE - SH161007
# --------------------------------
#OPT    +=  -DTIPSYFORMAT
#OPT	+=  -DOWLSFORMAT

CC= gcc
FC= f77
CLINK=gcc
FLINK=f77
GSL_INCL=-I/home/shuiyao/include
GSL_LIBS=-L/home/shuiyao/lib
CFLAGS= ${OPT} -g -Wall $(GSL_INCL) -O2
FFLAGS= -O
CLIB= -lm -lgsl -lgslcblas $(GSL_LIBS) #-lhdf5
FLIB= 

ifeq (PHEW,$(findstring PHEW,$(OPT)))
	SMOOTHER = contsmoothspec_phew.o
else
	SMOOTHER = contsmoothspec.o
endif

all: contspecexbin_v8

OBJSCONTV6= contspecexbin_v6.o cosmo.o getspecparticles.o $(SMOOTHER) ionfrac.o initions.o tau.o outtau.o
OBJSBIN2SPECV6 = bin2spec_v6.o cosmo.o ionfrac.o tau.o outtau.o 
OBJSCONTV8= contspecexbin_v8.o cosmo.o getspecparticles.o file_io.o $(SMOOTHER) ionfrac.o initions.o tau.o outtau.o # read_hdf5.o
OBJSBIN2SPECV8 = bin2spec_v8.o cosmo.o ionfrac.o tau.o outtau.o 


NROBJ= poidev.o gammln.o ran1.o

contspecexbin_v8:  $(OBJSCONTV8) $(NROBJ) Makefile
	$(CLINK) $(CFLAGS) -o contspecexbin_v8 $(NROBJ) $(OBJSCONTV8) $(CLIB)

contspecexbin_v7:  $(OBJSCONTV7) $(NROBJ) 
	$(CLINK) $(CFLAGS) -o contspecexbin_v7 $(NROBJ) $(OBJSCONTV7) $(CLIB)

specexbin_pipeline:  $(OBJSCONTV6) $(NROBJ) Makefile
	 $(CLINK) -DN2FORMAT -DPIPELINE -O2 -g -Wall -lm -o specexbin_pipeline $(NROBJ) $(OBJSCONTV6) $(CLIB)

bin2spec_v6:	$(OBJSBIN2SPECV6) Makefile
	$(CLINK) $(CFLAGS) -o bin2spec_v6  $(OBJSBIN2SPECV6) $(CLIB)

clean:
	rm -f *.o
	rm -f contspecexbin_v8

