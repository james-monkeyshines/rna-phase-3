include ../make.libs

#######################################################################
#  This makefile creates a cutdown version of the BLAS library
#  http://www.netlib.org
#
#---------------------------------------------------------------------
#  Edward Anderson, University of Tennessee
#  March 26, 1990
#  Susan Ostrouchov, Last updated September 30, 1994
#######################################################################

BLAS = $(shell ls ../blas/*.f | sed "s/\.f/\.o/" | sed "s/\/blas/\/lib/")

all: $(BLAS)
	$(ARCH) $(ARCHFLAGS) $(BLASLIB) $(BLAS)
	$(RANLIB) $(BLASLIB)

clean:
	rm -f ../lib/*.o $(BLASLIB)

../lib/%.o: %.f
	$(FORTRAN) $(OPTS) $(FLAGS) -o $@ -c $<
