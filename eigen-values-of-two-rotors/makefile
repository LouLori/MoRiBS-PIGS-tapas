#.SUFFIXES: .o .cxx

#options= -Ofast -march=native -fopenmp
options= -Ofast -fopenmp

CFLAGS =-I./sprng/include -I/usr/local/include -DEXCLUDEROTATION #-DTESTNORM #-DIOREAD

#Below is the LDFLAGS Toby Zeng use on nlogn
#LDFLAGS=-L/home/pnroy/Dev/lib64/ -lm -L./sprng/lib -llcg -L/home/pnroy/Dev/lib64/ -lgfortran  -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#below is the LDFLAGS with minimum flags
LDFLAGS= -lm -L./sprng/lib -lgfortran -lblas -llapack

#-------------------------------------------------------------------------
#  Compilers
#-------------------------------------------------------------------------

#CC=mpic++
#CC=g++
CC=g++
FC=gfortran
GCC=gcc
#FC=/home/pnroy/Dev/bin/gfortran

#-------------------------------------------------------------------------

runOBJS= plgndr.o gauleg.o vh2h2.o input.o

run:	eigen.o $(runOBJS)
		$(CC) -o $@  eigen.o $(runOBJS) $(LDFLAGS) -fopenmp

eigen.o:	eigen.cc eigen.h setup.h input.h
			$(CC) $(options) -c  $(CFLAGS) -o $@ $*.cc

input.o:	input.cc input.h
			$(CC) $(options) -c  $(CFLAGS) -o $@ $*.cc

.PHONY: clean
            
clean:
	rm -f *.o core run

wipe:
	rm -f mc00* mc00*.*old *.dump

.f.o:
	${FC} -O3 -c $<

.c.o:
	${GCC} -O3 -c $<
#---------------------------------------------------------------------------
# %.o: %.cxx
# $(CC) $(options) -c  $(CFLAGS) $(CPPFLAGS) $< -o $@
#---------------------------------------------------------------------------
