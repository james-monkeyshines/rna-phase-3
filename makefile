########################################################################
#  PHASE (Phylogenies and Sequence Evolution) makefile	         ######
#  PHASE version 3.0                                            ######
#                                                              ######
####################################################################

include makephase

CC = g++
LD = g++

# Change these two settings if you already have optimized versions of blas and
# lapack in your system (probably the case for Mac/Linux - see README).
#OPTLIBS = false
#LIBS = -llapack -lblas -lg2c -lm

# Probable best Linux settings
OPTLIBS = true
LIBS = -llapack -lblas -lm

# Probable best Mac settings
#OPTLIBS = true
#LIBS = -framework vecLib

##### Set the compiler options

CFLAGS = $(MACHINE) $(OPT) $(MACHINEOPT) -DNDEBUG
LFLAGS = $(MACHINE) $(OPT) $(MACHINEOPT) -L./lib -DNDEBUG

RM=/bin/rm -f

TOOLS = $(shell ls src/Tools/*.cc | sed "s/\.cc/\.o/" | sed "s/src/obj/")
MODELS = $(shell ls src/Models/*.cc | sed "s/\.cc/\.o/" | sed "s/src/obj/")
TREE = $(shell ls src/Tree/*.cc | sed "s/\.cc/\.o/" | sed "s/src/obj/")
UTILS = $(shell ls src/Util/*.cc | sed "s/\.cc/\.o/" | sed "s/src/obj/")
SEQUENCE = $(shell ls src/Sequence/*.cc | sed "s/\.cc/\.o/" | sed "s/src/obj/")
PROG = $(shell ls src/*.cc | sed "s/\.cc/\.o/" | sed "s/src/obj/")


SRCS = $(shell find src -name "*.cc")


#####  Compilation rules
all : dir libs mcmcphase analyzer optimizer mcmcsummarize mlphase likelihood simulate tools distphase

tools : dir splitdataset

#####  Create object/bin directories
dir:
	mkdir -p obj/
	mkdir -p obj/Models/
	mkdir -p obj/Tree/
	mkdir -p obj/Sequence/
	mkdir -p obj/Util/
	mkdir -p obj/Tools/
	mkdir -p lib/
	mkdir -p bin/

#####  BLAS/LAPACK libraries (Mac users see readme)
lib/liblapack.a : lapack/*.f
ifneq ($(OPTLIBS),true)
	(cd lapack ; make)
endif
lib/libblas.a : blas/*.f
ifneq ($(OPTLIBS),true)
	(cd blas ; make)
endif


#####  Set the location of the include files
INCLUDE = -I./include

libs : lib/liblapack.a lib/libblas.a

likelihood : obj/Likelihood.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE)
	$(CC) $(INCLUDE) $(LFLAGS) -o likelihood obj/Likelihood.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE) $(LIBS) -o bin/likelihood
optimizer : obj/Optimizer.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE)
	$(CC) $(INCLUDE) $(LFLAGS) -o optimizer obj/Optimizer.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE) $(LIBS) -o bin/optimizer
analyzer : obj/Analyzer.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE)
	$(CC) $(INCLUDE) $(LFLAGS) -o analyzer obj/Analyzer.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE) $(LIBS) -o bin/analyzer
mcmcsummarize : obj/MCMCsummarize.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE)
	$(CC) $(INCLUDE) $(LFLAGS) -o mcmcsummarize obj/MCMCsummarize.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE) $(LIBS) -o bin/mcmcsummarize
simulate : obj/Simulate.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE)
	$(CC) $(INCLUDE) $(LFLAGS) -o simulate obj/Simulate.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE) $(LIBS) -o bin/simulate
mlphase : obj/MLphase.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE)
	$(CC) $(INCLUDE) $(LFLAGS) -o mlphase obj/MLphase.o obj/Phase.o $(MODELS) $(TREE) $(UTILS) $(SEQUENCE) $(LIBS) -o bin/mlphase
mcmcphase : obj/MCMCphase.o obj/Phase.o $(MODELS) $(UTILS) $(TREE) $(SEQUENCE)
	$(CC) $(INCLUDE) $(LFLAGS) -o mcmcphase obj/MCMCphase.o obj/Phase.o $(MODELS) $(UTILS) $(TREE) $(SEQUENCE) $(LIBS) -o bin/mcmcphase
distphase : obj/DistPhase.o obj/Phase.o $(MODELS) $(UTILS) $(TREE) $(SEQUENCE)
	$(CC) $(INCLUDE) $(LFLAGS) -o distphase obj/DistPhase.o obj/Phase.o $(MODELS) $(UTILS) $(TREE) $(SEQUENCE) $(LIBS) -o bin/distphase

splitdataset : obj/Tools/splitDataSet.o obj/Util/FileParser.o obj/Util/ParametersSet.o obj/Sequence/SequenceTable.o
	$(CC) $(INCLUDE) $(LFLAGS) obj/Tools/splitDataSet.o obj/Util/FileParser.o obj/Util/ParametersSet.o obj/Sequence/SequenceTable.o -o bin/splitdataset


obj/%.o : src/%.cc
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<


### Clean up
clean:
	$(RM) $(PROG) $(MODELS) $(TREE) $(UTILS) $(SEQUENCE) $(TOOLS)
	$(RM) bin/*
	$(RM) *.a

libclean:
	(cd lapack; make clean);
	(cd blas; make clean);
