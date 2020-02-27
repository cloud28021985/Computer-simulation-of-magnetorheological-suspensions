# 2020 Computer simulations of anisotropic structures in magnetorheological elastomers


CFLAGS = -c -Ofast -Wall
RUN = mpirun -np # run command of programm 'prog'
N_PROC = 2 # number of the processors
SRCPATH = src/
OBJPATH = obj/
VMDPATH = vmd/
DATAPATH = data/
FIGPATH = figs/
SCRIPT = script.tcl
PLOTSRC = plot.gp
SOURCES = $(wildcard $(SRCPATH)*.f90)
OBJECTS = $(patsubst $(SRCPATH)%.f90, $(OBJPATH)%.o, $(SOURCES))
EXECUTABLE = $(OBJPATH)prog


all:
	$(RUN) $(N_PROC) $(EXECUTABLE)
	gnuplot plot.gp


all: $(EXECUTABLE)


$(EXECUTABLE): $(OBJECTS)
	mpif77 $(OBJECTS) -lm -o $@


$(OBJPATH)%.o: $(SRCPATH)%.f90
	mpif77 $(CFLAGS) $< -o $@


$(OBJECTS): makefile


# command $make clean
clean:
	rm -rf $(OBJPATH)
	rm -rf $(VMDPATH)
	rm -rf $(DATAPATH)
	rm -rf $(FIGPATH)
	mkdir $(OBJPATH)
	mkdir $(VMDPATH)
	mkdir $(DATAPATH)
	mkdir $(FIGPATH)
	cp $(SCRIPT) $(VMDPATH)$(SCRIPT)
