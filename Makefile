SrcSuf        = cc
ObjSuf        = o
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

OPT           = -O2

CXX           = g++
CXXFLAGS      = $(OPT) -Wall -fPIC
LD            = g++
LDFLAGS       = $(OPT)
SOFLAGS       = -shared -Wl -m64 -g
INCDIR        = include
SRCDIR        = src
CXXFLAGS     += -I./$(INCDIR) -I$(ROOTSYS)/include


GEANTSOURCES  = PhysicsProcess.cxx GeantTrack.cxx GeantOutput.cxx GeantVolumeBasket.cxx GeantPropagator.cxx GeantBasket.cxx WorkloadManager.cxx sync_objects.cxx
GEANTDICTS    = G__Geant.cxx
GEANTO        = PhysicsProcess.o GeantTrack.o GeantOutput.o GeantBasket.o GeantVolumeBasket.o GeantPropagator.o WorkloadManager.o sync_objects.o G__Geant.o
GEANTSL       = $(patsubst %.o,$(SRCDIR)/%.o,$(GEANTO))
GEANTDEPS     = $(patsubst %.cxx,$(INCDIR)/%.h,$(GEANTSOURCES))

GEANTSO       = libGeant.$(DllSuf)
GEANTLIB      = $(shell pwd)/$(GEANTSO)
OBJS          = $(GEANTSO)


#------------------------------------------------------------------------------

all:            $(GEANTSO)

$(GEANTSO):     $(GEANTDICTS) $(GEANTO)
		$(LD) $(LDFLAGS) $(SOFLAGS) $(GEANTSL) $(OutPutOpt)$@
		@echo "$@ done"
PhysicsProcess.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/PhysicsProcess.cxx
GeantTrack.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantTrack.cxx
GeantBasket.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantBasket.cxx
GeantOutput.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantOutput.cxx
GeantVolumeBasket.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantVolumeBasket.cxx
GeantPropagator.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantPropagator.cxx
WorkloadManager.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/WorkloadManager.cxx
sync_objects.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/sync_objects.cxx
$(GEANTDICTS):
		@echo "Generating dictionary $@"
#		cp $(INCBRIDGE)/*.h $(INCDIR)
		$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/$@ \
                -c $(INCDIR)/PhysicsProcess.h $(INCDIR)/GeantVolumeBasket.h $(INCDIR)/GeantBasket.h $(INCDIR)/GeantPropagator.h \
		$(INCDIR)/GeantTrack.h $(INCDIR)/GeantOutput.h $(INCDIR)/WorkloadManager.h $(INCDIR)/LinkDef.h
G__Geant.o:
		$(ROOTSYS)/bin/rmkdepend -R -f$(SRCDIR)/G__Geant.d -Y -w 1000 -- \
                pipe -m64 -Wshadow -Wall -W -Woverloaded-virtual -fPIC \
                -I$(INCDIR) -pthread  -D__cplusplus -I$(ROOTSYS)cint/cint/lib/prec_stl \
                -I$(ROOTSYS)cint/cint/stl -I$(ROOTSYS)/cint/cint/inc -- $(SRCDIR)/$(GEANTDICTS)
		$(CXX)  $(CXXFLAGS) -I. -I$(ROOTSYS)/include -o $(SRCDIR)/$@ \
                -c $(SRCDIR)/G__Geant.cxx
clean:
		@rm -rf $(GEANTSO) src/*.o src/*.d src/G__*
                

