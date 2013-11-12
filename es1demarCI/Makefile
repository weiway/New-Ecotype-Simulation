build:
	gfortran src/Methods_print.f90 -c  -o Methods.o
	gfortran src/NelderMead.f90 -c  -o NelderMead.o
	gfortran src/demarcationsCI_print.f90 Methods.o NelderMead.o  -o demarcationsCI.app
	rm -f methods.mod
	rm -f neldermead.mod
	rm -f Methods.o
	rm -f Neldermead.o