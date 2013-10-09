# This is a very simple makefile just for demonstration purposes
# you need to build libGeant.so before making this executable
# you need to make sure that the libtbb is in path ${TBBROOT}/lib/ and ${TBBROOT} needs to be set as an environment variable

SYSTEM := $(shell uname)

CXXFLAGS      = $(shell root-config --cflags)

MAINSRC  = runmain.cxx
MAINEXE  = runmain.x

# determine root libs
ROOTLIBS    = $(shell root-config --libs)
ROOTLIBS    += "-lGeom" # for geometry
ROOTLIBS    += "-lVMC" # for TParticleDatabase
GEANTINCDIR = ../inc


#------------------------------------------------------------------------------

all:            $(MAINEXE)

$(MAINEXE):  $(MAINSRC)
		$(CXX) $(CXXFLAGS) -I$(GEANTINCDIR) $< -o $@ $(ROOTLIBS)  -L../ -lGeant -L$(TBBROOT)/lib -ltbb
		@echo "$@ done"
