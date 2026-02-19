!FC = ifort
FC = mpif90
!FC = f90
!FC = esf90
!FC = esmpif90

OPT_IFC = -mcmodel=large 
#OPT_IFC = -mcmodel=large -CB -traceback -convert -big_endian
EXEC = MergeCode
FILES = main subroutines ana common constants
OBJS = ${FILES:=.o}
INC = 

FFLAGS = $(OPT_IFC)

all : $(OBJS) $(INC)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS)

### rules ###
.f90.o:
	$(FC) -c $(FFLAGS) $< 

# main.o : subroutines.o main.f90
#	$(FC) -c $(FFLAGS) main.f90
main.o : subroutines.o ana.o main.f90
	$(FC) -c $(FFLAGS) main.f90
subroutines.o : common.o subroutines.f90
	$(FC) -c $(FFLAGS) subroutines.f90
ana.o  : common.o ana.f90
	$(FC) -c $(FFLAGS) ana.f90
common.o : constants.o common.f90
	$(FC) -c $(FFLAGS) common.f90
constants.o : constants.f90
	$(FC) -c $(FFLAGS) constants.f90

clean : 
	rm -f *.o *.mod $(EXEC)

