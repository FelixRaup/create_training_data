FLAGS = -DCHEMISTRYNETWORK=5 -g -fbacktrace -fcheck=all -std=legacy -O2 -fdefault-real-8 -ffree-line-length-0 -fallow-argument-mismatch -fdefault-double-8 -fbounds-check

FC = gfortran $(FLAGS)

OBJS  = coolinmo.o cheminmo.o spline.o cool_func.o photoinit_ism.o \
		         dvode.o evolve_abundances.o rate_eq.o jac.o cool_util.o \
				          const_rates.o validate_rates.o test_harness.o

test: $(OBJS)
	    $(FC) -o test $(OBJS)

