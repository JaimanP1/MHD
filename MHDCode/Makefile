#FC = sxmpif90
#FC = /SX/opt/mpisx/LC310106/bin/sxmpif90_313
FC = mpif90
#FC = mpifrt

EXEC = mhd_run
OBJ = constants.o common.o mpiut.o pset.o mhd.o job.o iset.o rkg.o out.o r2f.o dcft2.o 
INC = constants.f90 

#OPT = -ftrace
#OPT = -ftrace -Wf"-pvctl fullmsg"
#OPT = -P auto -ftrace
#OPT = -mcmodel=large -i-dynamic -CB -traceback
OPT = -mcmodel=medium
#FLIB = -lasl

all : main.f90 $(OBJ)
	$(FC) $(OPT) -o $(EXEC) main.f90 $(OBJ) $(FLIB)

constants.o : constants.f90
	$(FC) $(OPT) -c constants.f90

common.o : constants.o common.f90
	$(FC) $(OPT) -c common.f90

mpiut.o : common.o mpiut.f90
	$(FC) $(OPT) -c mpiut.f90

pset.o : common.o mpiut.o pset.f90
	$(FC) $(OPT) -c pset.f90

job.o : constants.o common.f90 job.f90
	$(FC) $(OPT) -c job.f90

iset.o : common.o mhd.o iset.f90 
	$(FC) $(OPT) -c iset.f90

rkg.o : common.o rkg.f90
	     $(FC) $(OPT) -c rkg.f90

out.o : common.o mpiut.o out.f90
	$(FC) $(OPT) -c out.f90

mhd.o : common.o mpiut.o mhd.f90
	$(FC) $(OPT) -c mhd.f90

r2f.o : common.o r2f.f90
	$(FC) $(OPT) -c r2f.f90

dcft2.o : dcft2.f
	$(FC) $(OPT) -c dcft2.f

clean : 
	rm *.o mhd_run

cleanall : 
	rm *.o *.mod $(EXEC)






