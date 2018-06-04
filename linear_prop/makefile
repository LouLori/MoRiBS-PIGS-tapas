FC=gfortran
all: qpr.x linden.x
linden.x:	linden.o
		$(FC) -o $@ linden.o

linden.o:	linden.f
		$(FC) -c linden.f
qpr.x:	qpr.o
		$(FC) -o $@ qpr.o

qpr.o:	qpr.f
		$(FC) -c qpr.f

clean:
		rm -f *.o linden.x
