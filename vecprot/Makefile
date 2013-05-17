-include Config.mk

SYSTEM := $(shell uname)

SrcSuf        = cc
ObjSuf        = o
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

OPT           = -g 
# -O2

CXX			  = $(shell root-config --cxx)
LD			  = $(shell root-config --ld)
CXXFLAGS      = $(shell root-config --cflags)

ifneq (,$(findstring icc,$(CXX)))
CXXFLAGS	  += -fPIC
endif
ifeq ($(CXX),g++) 
CXXFLAGS      += $(OPT) -Wall -fPIC
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
SOFLAGS	      = -shared -Wl -m64 -g
endif
endif
INCDIR        = inc
SRCDIR        = src
CXXFLAGS     += -I./$(INCDIR) -I$(ROOTSYS)/include


GEANTSOURCES  = GeantThreadData.cxx PhysicsProcess.cxx GeantTrack.cxx \
                GeantOutput.cxx GeantVolumeBasket.cxx GeantPropagator.cxx \
                GeantBasket.cxx WorkloadManager.cxx sync_objects.cxx
GEANTDICTS    = $(SRCDIR)/G__Geant.cxx
GEANTO        = GeantThreadData.o PhysicsProcess.o GeantTrack.o \
                GeantOutput.o GeantBasket.o GeantVolumeBasket.o GeantPropagator.o \
                WorkloadManager.o sync_objects.o G__Geant.o
GEANTSOURCES := $(patsubst %.cxx,$(SRCDIR)/%.cxx,$(GEANTSOURCES))
GEANTO       := $(patsubst %.o,$(SRCDIR)/%.o,$(GEANTO))
GEANTDEPS     = $(patsubst %.cxx,%.d,$(GEANTSOURCES))

GEANTSO       = libGeant.$(DllSuf)
GEANTLIB      = $(shell pwd)/$(GEANTSO)
OBJS          = $(GEANTSO)


#------------------------------------------------------------------------------

all:            $(GEANTSO)

-include $(GEANTDEPS)

ifeq ($(EnableCuda),yes)
include Makefile.cuda
endif

$(GEANTSO):     $(GEANTDICTS) $(GEANTO)
		$(LD) $(LDFLAGS) $(SOFLAGS) $(GEANTO) $(OutPutOpt)$@
		@echo "$@ done"

$(SRCDIR)/%.o: $(SRCDIR)/%.cxx
		$(ROOTSYS)/bin/rmkdepend -R -f$(SRCDIR)/$*.d -Y -w 1000 -- \
	             	$(CXXFLAGS) \
		 	 -- $<
		$(CXX)  $(CXXFLAGS) -o $@ -c $<

$(GEANTDICTS):
		@echo "Generating dictionary $@"
#		cp $(INCBRIDGE)/*.h $(INCDIR)
		$(ROOTSYS)/bin/rootcint -f $@ \
                -c $(INCDIR)/GeantThreadData.h $(INCDIR)/PhysicsProcess.h \
                $(INCDIR)/GeantVolumeBasket.h $(INCDIR)/GeantBasket.h \
                $(INCDIR)/GeantPropagator.h $(INCDIR)/GeantTrack.h $(INCDIR)/GeantOutput.h \
                $(INCDIR)/WorkloadManager.h $(INCDIR)/CoprocessorBroker.h $(INCDIR)/LinkDef.h
src/G__Geant.o:
		$(ROOTSYS)/bin/rmkdepend -R -f$(SRCDIR)/G__Geant.d -Y -w 1000 -- \
                pipe -m64 -Wshadow -Wall -W -Woverloaded-virtual -fPIC \
                -I$(INCDIR) -pthread  -D__cplusplus -I$(ROOTSYS)cint/cint/lib/prec_stl \
                -I$(ROOTSYS)cint/cint/stl -I$(ROOTSYS)/cint/cint/inc -- $(SRCDIR)/$(GEANTDICTS)
		$(CXX)  $(CXXFLAGS) -I. -I$(ROOTSYS)/include -o $@ \
                -c $(SRCDIR)/G__Geant.cxx
clean:
		@rm -rf $(GEANTSO) src/*.o src/*.d src/G__*


# DO NOT DELETE
