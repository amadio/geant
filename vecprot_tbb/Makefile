
SYSTEM := $(shell uname)

SrcSuf        = cc
ObjSuf        = o
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

OPT           = -O2

CXX			  = $(shell root-config --cxx)
LD			  = $(shell root-config --ld)
CXXFLAGS      = $(shell root-config --cflags)

ifneq (,$(findstring icc,$(CXX)))
CXXFLAGS	  += -fPIC
endif
ifeq ($(CXX),g++)
CXXFLAGS      += $(OPT) -fPIC # -Wall
LD            = g++
LDFLAGS       = $(OPT)
endif
ifeq ($(CXX),clang++)
CXXFLAGS      += $(OPT) -Wall
LD            = $(CXX)
LDFLAGS       = $(OPT)
endif
ifeq ($(SYSTEM),Darwin)
SOFLAGS       = -shared -m64 -dynamiclib -undefined dynamic_lookup -single_module
else
ifneq (,$(findstring icc,$(CXX)))
SOFLAGS	      = -shared -m64
else
SOFLAGS	      = -shared -m64 -g
endif
endif
INCDIR        = inc
SRCDIR        = src
CXXFLAGS     += -I./$(INCDIR) -I$(ROOTSYS)/include

GEANTSOURCES  = GeantThreadData.cxx PhysicsProcess.cxx GeantTrack.cxx \
                GeantOutput.cxx GeantVolumeBasket.cxx GeantBasket.cxx \
                GeantMainScheduler.cxx \
                GeantPropagator.cxx \
                GeantMainPropagator.cxx \
                Event_Hit.cxx \
                PropTask.cxx CollDispTask.cxx \
                WorkloadManager.cxx

GEANTDICTS    = G__Geant.cxx
GEANTO        = GeantThreadData.o PhysicsProcess.o GeantTrack.o \
                GeantOutput.o GeantVolumeBasket.o GeantBasket.o \
                GeantMainScheduler.o \
                GeantPropagator.o \
                GeantMainPropagator.o \
                Event_Hit.o \
                PropTask.o CollDispTask.o \
                WorkloadManager.o G__Geant.o

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
GeantThreadData.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantThreadData.cxx
PhysicsProcess.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/PhysicsProcess.cxx
GeantTrack.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantTrack.cxx
GeantOutput.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantOutput.cxx
GeantVolumeBasket.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantVolumeBasket.cxx
GeantBasket.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantBasket.cxx
GeantMainScheduler.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantMainScheduler.cxx
GeantPropagator.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantPropagator.cxx
GeantMainPropagator.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantMainPropagator.cxx
Event_Hit.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/Event_Hit.cxx
PropTask.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/PropTask.cxx
CollDispTask.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/CollDispTask.cxx
WorkloadManager.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/WorkloadManager.cxx

$(GEANTDICTS):
		@echo "Generating dictionary $@"
#		cp $(INCBRIDGE)/*.h $(INCDIR)
		$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/$@ \
             -c $(INCDIR)/PhysicsProcess.h \
                $(INCDIR)/GeantOutput.h \
                $(INCDIR)/GeantVolumeBasket.h \
                $(INCDIR)/GeantBasket.h \
                $(INCDIR)/GeantMainScheduler.h \
                $(INCDIR)/GeantMainPropagator.h \
                $(INCDIR)/LinkDef.h

G__Geant.o:
		$(ROOTSYS)/bin/rmkdepend -R -f$(SRCDIR)/G__Geant.d -Y -w 1000 -- \
                pipe -m64 -Wshadow -Wall -W -Woverloaded-virtual -fPIC \
                -I$(INCDIR) -pthread  -D__cplusplus -I$(ROOTSYS)cint/cint/lib/prec_stl \
                -I$(ROOTSYS)cint/cint/stl -I$(ROOTSYS)/cint/cint/inc -- $(SRCDIR)/$(GEANTDICTS)
		$(CXX)  $(CXXFLAGS) -I. -I$(ROOTSYS)/include -o $(SRCDIR)/$@ \
                -c $(SRCDIR)/G__Geant.cxx
clean:
		@rm -rf $(GEANTSO) src/*.o src/*.d src/G__*
		@rm -f src/*~ inc/*~ *~
