CC = gcc
FC = gfortran

default:	autovp

FFLAGS = -O

COMPILE = gfortran
OPTIONS = -O
SOURCE = autofit.f ew_finder.f smooth.f fitregion.f model.f combine.f fitsat.f inparam.f readdata.f polyfit.f doublet.f
OBJECTS = $(SOURCE:.f=.o)

autofit:	$(OBJECTS)
	$(COMPILE) $(OPTIONS) $(OBJECTS) -o autofit

$(OBJECTS): Makefile vpdefs.h const.dek

SOURC2 = minfit.f minimize.f model.f ew_finder.f oneparmin.f tossline.f splitreg.f readdata.f initparm.f fixparm.f loadline.f funcs.f nparmin.f inparam.f polyfit.f piksrt.f piksr2.f addline.f getlines.f doublet.f
OBJECT2 = $(SOURC2:.f=.o)

minfit:	$(OBJECT2)
	$(COMPILE) $(OPTIONS) $(OBJECT2) -o minfit

$(OBJECT2): Makefile vpdefs.h const.dek

SOURC3 = fitcont.f
OBJECT3 = $(SOURC3:.f=.o)

autovp:	autofit minfit 

clean:
	rm -f *.o
	rm -f autofit minfit