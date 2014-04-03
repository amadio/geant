-include Config.mk

SYSTEM := $(shell uname)

SrcSuf        = cc
ObjSuf        = o
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

OPT           = -g

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
SOFLAGS	      = -shared -m64 -g
endif
endif
INCDIR        = inc
SRCDIR        = src
ifeq ($(GEANT_DEBUG),1)
CXXFLAGS     += -D__STAT_DEBUG
endif
CXXFLAGS     += -I./$(INCDIR) -I$(ROOTSYS)/include


GEANTSOURCES  = GeantThreadData.cxx GeantTrack.cxx GeantBasket.cxx  sync_objects.cxx \
                GeantOutput.cxx GeantPropagator.cxx GeantEvent.cxx GeantScheduler.cxx \
                PhysicsProcess.cxx WorkloadManager.cxx GeantVApplication.cxx GeantFactoryStore.cxx \
                GeantTrackStat.cxx MyHit.cxx MyApplication.cxx
GEANTH        = GeantThreadData.h PhysicsProcess.h GeantBasket.h GeantScheduler.h GeantPropagator.h \
                GeantEvent.h GeantTrack.h GeantOutput.h WorkloadManager.h GeantVApplication.h \
                GeantTrackStat.h MyHit.h MyApplication.h LinkDef.h
GEANTHI       = $(patsubst %.h,inc/%.h,$(GEANTH))
GEANTDICTS    = $(SRCDIR)/G__Geant.cxx
GEANTO        = $(patsubst %.cxx,src/%.o,$(GEANTSOURCES)) $(patsubst %.cxx,%.o,$(GEANTDICTS))
GEANTDEPS     = $(patsubst %.cxx,%.d,$(GEANTSOURCES)) $(SRCDIR)/G__Geant.d 

GEANTSO       = libGeant.$(DllSuf)
GEANTLIB      = $(shell pwd)/$(GEANTSO)

OBJS          = $(GEANTSO)


#------------------------------------------------------------------------------

all:            $(GEANTSO)
geant:          $(GEANTSO)

ifeq ($(EnableCuda),yes)
include cuda/Makefile.cuda
endif

$(GEANTSO):     $(GEANTDICTS) $(GEANTO)
		$(LD) $(LDFLAGS) $(SOFLAGS) $(GEANTO) $(OutPutOpt)$@
		@echo "$@ done"

src/%.o: src/%.cxx inc/%.h compiler_flags
		$(ROOTSYS)/bin/rmkdepend -R -f$(SRCDIR)/$*.d -Y -w 1000 -- \
	             	$(CXXFLAGS) \
		 	 -- $<
		$(CXX)  $(CXXFLAGS) -MMD -o $@ -c $<

-include $(GEANTDEPS)

$(GEANTDICTS): $(GEANTHI) 
		@echo "Generating dictionary $@"
		$(ROOTSYS)/bin/rootcint -f $@ -c -p $(GEANTHI)

$(SRCDIR)/G__Geant.o:     $(GEANTDICTS)
		$(ROOTSYS)/bin/rmkdepend -R -f $(SRCDIR)/G__Geant.d -Y -w 1000 -- \
                   pipe -m64 -Wshadow -Wall -W -Woverloaded-virtual -fPIC \
                   -I$(INCDIR) -pthread  -D__cplusplus -I$(ROOTSYS)cint/cint/lib/prec_stl \
                   -I$(ROOTSYS)cint/cint/stl -I$(ROOTSYS)/cint/cint/inc -- $(GEANTDICTS)
		$(CXX)  $(CXXFLAGS) -I. -I$(ROOTSYS)/include -o $@ -c $(SRCDIR)/G__Geant.cxx

clean:
		@rm -rf $(GEANTSO) $(USERSO) $(GEANTDEPS) src/*.o src/G__*

.PHONY: force
compiler_flags: force
		@echo '$(CXXFLAGS)' | cmp -s - $@ || echo '$(CXXFLAGS)' > $@

