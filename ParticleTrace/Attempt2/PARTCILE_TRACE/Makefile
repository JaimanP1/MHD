FORT = mpif90
#OPT = -O3 -parallel -mcmodel=large -i-dynamic -CB -traceback
OPT = -mcmodel=large 

#LIBS =  -L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11 -lm
SRC = main.f90 bline.f90 sint.f90 out.f90  
INC = paramete.f90 common.f90
OBJ = $(SRC:%.f90=%.o)
EXEC = bline

main : $(OBJ) 
	$(FORT) $(OPT) $(LIBS) $(OBJ) -o $(EXEC)

%.o : %.f90 $(INC)
	$(FORT) -c $< $(OPT) 

clean:
	rm -f $(OBJ) $(EXEC)

