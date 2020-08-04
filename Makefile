#OPT     +=  -DTIPSYFORMAT
OPT     +=  -DHDF5FORMAT # See appendix. on the metallicity fields

OPT     +=  -DDENSITY_H2_FACTOR
OPT	+=  -DOUTPUT_LOCAL_FOLDER
OPT    +=   -DNO_EXTRA_OUTPUT # no binzfile partfile etc.
OPT	+=  -DQUINTIC_KERNEL
OPT	+=  -DVELOCITY_UNIT_CORRECTION

OPT	+=  -DPART_BY_PART # Mutually exclusive with PHEW
# OPT	+=  -DSINGLE_VOFFSET_PER_PARTICLE # Not recommended. The Pygad way.

OPT	+=  -DPHEW # Mutually exclusive with PART_BY_PART
#OPT     +=  -DPHEW_IGNORE_PHEWS
# OPT	+=  -DPHEW_MOREINFO
# OPT	+=  -DPHEW_VERBOSE=1
#OPT	+=  -DPHEW_HSMOOTH
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
#OPT     +=  -DVARGALBKGD
#OPT	+=  -DMETALFLOOR
#OPT     +=  -DSHORTSPEC
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

ifeq (HDF5FORMAT,$(findstring HDF5FORMAT,$(OPT)))
	HDF5_INCL = -I/usr/include/hdf5/serial -DH5_USE_16_API
	HDF5_LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_serial -lz
	# HDF5_LIBS= -lhdf5
endif

CC= gcc
FC= f77
CLINK=gcc
FLINK=f77
GSL_INCL=-I/home/shuiyao/local/include
GSL_LIBS=-L/home/shuiyao/local/lib
CFLAGS= ${OPT} -g $(GSL_INCL) $(HDF5_INCL) #-fPIC -shared
FFLAGS= -O
CLIB= -lm -lgsl -lgslcblas $(GSL_LIBS) $(HDF5_LIBS)
FLIB= 


ifeq (PHEW,$(findstring PHEW,$(OPT)))
	SMOOTHER = contsmoothspec.o
else
	SMOOTHER = contsmoothspec.o
endif

all: contspecexbin

OBJSCONT= contspecexbin.o cosmo.o getspecparticles.o file_io.o $(SMOOTHER) ionfrac.o initions.o tau.o outtau.o loadhdf5.o

NROBJ= poidev.o gammln.o ran1.o

contspecexbin:  $(OBJSCONT) $(NROBJ) Makefile
	$(CLINK) $(CFLAGS) -o contspecexbin $(NROBJ) $(OBJSCONT) $(CLIB)

specexbin_pipeline:  $(OBJSCONTV6) $(NROBJ) Makefile
	 $(CLINK) -DN2FORMAT -DPIPELINE -O2 -g -Wall -lm -o specexbin_pipeline $(NROBJ) $(OBJSCONTV6) $(CLIB)

clean:
	rm -f *.o
	rm -f contspecexbin


# Appendix 1. The GIZMO Metallicities. (galaxy_sf/metals.c)
#     All.SolarAbundances[1]=0.28;    // He  (10.93 in units where log[H]=12, so photospheric mass fraction -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314)                                                                
#     All.SolarAbundances[2]=3.26e-3; // C   (8.43 -> 2.38e-3, AG=3.18e-3)        
#     All.SolarAbundances[3]=1.32e-3; // N   (7.83 -> 0.70e-3, AG=1.15e-3)        
#     All.SolarAbundances[4]=8.65e-3; // O   (8.69 -> 5.79e-3, AG=9.97e-3)        
#     All.SolarAbundances[5]=2.22e-3; // Ne  (7.93 -> 1.26e-3, AG=1.72e-3)        
#     All.SolarAbundances[6]=9.31e-4; // Mg  (7.60 -> 7.14e-4, AG=6.75e-4)        
#     All.SolarAbundances[7]=1.08e-3; // Si  (7.51 -> 6.71e-4, AG=7.30e-4)        
#     All.SolarAbundances[8]=6.44e-4; // S   (7.12 -> 3.12e-4, AG=3.80e-4)        
#     All.SolarAbundances[9]=1.01e-4; // Ca  (6.34 -> 0.65e-4, AG=0.67e-4)        
#     All.SolarAbundances[10]=1.73e-3; // Fe (7.50 -> 1.31e-3, AG=1.92e-3)   
