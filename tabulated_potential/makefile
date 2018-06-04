FC=gfortran
all: hfc60.x
hfc60.x:	hfc60.o
		$(FC) -o $@ hfc60.o

hfc60.o:	hfc60.f
		$(FC) -c hfc60.f

clean:
		rm -f *.o hfc60.x
