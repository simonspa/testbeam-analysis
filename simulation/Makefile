#
# D. Pitzl, 2013
# original Fortran code converted with f2c -C++
#	-std=c++0x \
#
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs) # with Gui

pixelav: pixelav.c Makefile
	gcc  pixelav.c \
	-Wall -O2 -o pixelav -L/usr/lib64/cernlib/2006/lib -lmathlib -lf2c -lm  -lgsl -L/usr/lib64/atlas -lcblas
	@echo 'done: pixelav'

evrd:	evrd.C Makefile
	@echo 'root C flags = ' $(ROOTCFLAGS)
	g++ $(ROOTCFLAGS) evrd.C \
	-Wall -O2 -o evrd $(ROOTLIBS)
	@echo 'done: evrd'

sim:	sim.C Makefile
	@echo 'root C flags = ' $(ROOTCFLAGS)
	g++ $(ROOTCFLAGS) sim.C \
	-Wall -O2 -o sim $(ROOTLIBS)
	@echo 'done: sim'

mobil: mobil.cpp Makefile
	g++ $(ROOTCFLAGS) mobil.cpp  \
	-Wall -Wextra -O2 -o mobil $(ROOTGLIBS)
	@echo 'done: mobil'
