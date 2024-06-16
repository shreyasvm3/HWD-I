bindir = ../bin/
srcdir = ./src
FC = mpif90 -O2
LIB = -llapack -lblas
.SUFFIXES: .f90 .o

#--------------------------------------------------

MOD = traj_mD.f90 \
      supply_mD.f90 \
      potential_mD.f90 \
      MonteCarlo_mD.f90 \
      params_mD.o

SRC = LSC_timecorr_mD.f90

OBJ = $(MOD) $(SRC)

# rule to make a .o file from a .f90 file
.f90.o:
	$(FC) -c $<
.f.o:
	$(FC) -c $<

dyn.x :  $(OBJ)
	 $(FC) $(OPTION) -o $@ $(OBJ) $(LIB)

clear   :
	rm *.x
	rm *.o
	rm *.mod
