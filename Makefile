#OPT     +=  -DTIPSYFORMAT
OPT     +=  -DHDF5FORMAT
OPT     +=  -DDENSITY_H2_FACTOR
OPT	+=  -DOUTPUT_LOCAL_FOLDER
OPT	+=  -DQUINTIC_KERNEL
OPT	+=  -DVELOCITY_UNIT_CORRECTION

# OPT	+=  -DPART_BY_PART # Mutually exclusive with PHEW
# OPT	+=  -DSINGLE_VOFFSET_PER_PARTICLE # Not recommended. The Pygad way.

OPT	+=  -DPHEW # Mutually exclusive with PART_BY_PART
# OPT	+=  -DPHEW_MOREINFO
# OPT	+=  -DPHEW_VERBOSE=1
OPT	+=  -DPHEW_HSMOOTH
OPT	+=  -DPHEW_RCLOUD_CORRECTION # Should be taken out then
OPT	+=  -DPHEW_NCLOUD=1000.0 # Should be taken out in the end

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

ifeq (HDF5FORMAT,$(findstring HDF5FORMAT,$(OPT)))
	HDF5_INCL = -I/usr/include/hdf5/serial -DH5_USE_16_API
	HDF5_LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_serial -lz
	# HDF5_LIBS= -lhdf5
endif

CC= gcc
FC= f77
CLINK=gcc
FLINK=f77
GSL_INCL=-I/home/shuiyao/include
GSL_LIBS=-L/home/shuiyao/lib
CFLAGS= ${OPT} -g $(GSL_INCL) $(HDF5_INCL) -O2
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

