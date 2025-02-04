# Compiler and executable
FC = mpif90
EXEC = mhd_run.o

# Files and dependencies
SRCS = constants.f90 common.f90 mpiut.f90 pset.f90 mhd.f90 job.f90 iset.f90 rkg.f90 out.f90 r2f.f90 dcft2.f main.f90
OBJS = $(SRCS:.f90=.o)
MODS = $(SRCS:.f90=.mod)

# Compilation options
OPT = -mcmodel=medium
DEBUG = -g -O0 -fcheck=all -Wall

# Default target
all: $(EXEC)

# Debug build
debug: OPT += $(DEBUG)
debug: $(EXEC)

# Linking
$(EXEC): $(OBJS)
	$(FC) $(OPT) -o $@ $(OBJS)

# Pattern rule for object files
%.o: %.f90
	$(FC) $(OPT) -c $< -o $@

%.o: %.f
	$(FC) $(OPT) -c $< -o $@

# Clean targets
clean:
	rm -f *.o $(EXEC)

cleanall:
	rm -f *.o *.mod $(EXEC)

