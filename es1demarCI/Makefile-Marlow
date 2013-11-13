build:
	gfortran src/Methods_noprint.f90 -c  -o Methods.o
	gfortran src/NelderMead.f90 -c  -o NelderMead.o
	gfortran src/demarcationsCI_noprint.f90 Methods.o NelderMead.o -cpp -fbounds-check -fopenmp -g -lgfortran -lm -O3 -Wall -Wextra -o demarcationsCI.i386
	rm -f methods.mod
	rm -f neldermead.mod
	rm -f Methods.o
	rm -f Neldermead.o
