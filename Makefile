F90       = ifort
EXEC      = Monte_Carlo
all: $(EXEC)
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
	@echo "Successfully created $(EXEC) binary file"
	@echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

OBJS   =  date_mod.o eigen_mod.o eof_mod.o inout_mod.o params_mod.o utils_mod.o critical_statistic_mod.o dynamical_analogue_mod.o Monte_Carlo.o

.SUFFIXES:

.SUFFIXES: .f90 .o

$(EXEC): $(OBJS)
	$(F90) -o $@ $(OBJS)

.f90.o:
	$(F90) -o $@ -c $<

clean:
	rm -f *.mod *.o $(EXEC)

neat:
	-rm -f *.mod *.o
