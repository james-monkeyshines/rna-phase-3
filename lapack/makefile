include ../make.libs

#######################################################################
#  This the makefile to create a cutdown version of the LAPACK library
#  http://www.netlib.org
#
#  ***Note***
#  The functions lsame, second, dsecnd, slamch, and dlamch may have
#  to be installed before compiling the library.  Refer to the
#  installation guide, LAPACK Working Note 41, for instructions.
#
#######################################################################


LAPACK = $(shell ls ../lapack/*.f | sed "s/\.f/\.o/" | sed "s/\/lapack/\/lib/")

all: $(LAPACK) #$(DLASRC) $(ALLAUX)
	$(ARCH) $(ARCHFLAGS) $(LAPACKLIB) $(LAPACK) #$(DLASRC) $(ALLAUX)
	$(RANLIB) $(LAPACKLIB)

clean:
	rm -f ../lib/*.o $(LAPACKLIB)

../lib/slamch.o: slamch.f
	$(FORTRAN) $(NOOPT) -o $@ -c $<
../lib/dlamch.o: dlamch.f
	$(FORTRAN) $(NOOPT) -o $@ -c $<


../lib/%.o: %.f
	$(FORTRAN) $(OPTS) -o $@ -c $<
