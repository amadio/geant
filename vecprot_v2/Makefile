
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
CXXFLAGS     += -I./$(INCDIR) -I$(ROOTSYS)/include


GEANTSOURCES  = GeantThreadData.cxx GeantTrack.cxx GeantBasket.cxx  sync_objects.cxx \
                GeantOutput.cxx GeantPropagator.cxx GeantEvent.cxx GeantScheduler.cxx \
                PhysicsProcess.cxx WorkloadManager.cxx GeantVApplication.cxx GeantFactoryStore.cxx
GEANTH        = GeantThreadData.h PhysicsProcess.h GeantBasket.h GeantScheduler.h GeantPropagator.h \
                GeantEvent.h GeantTrack.h GeantOutput.h WorkloadManager.h GeantVApplication.h LinkDef.h
GEANTHI       = $(patsubst %.h,inc/%.h,$(GEANTH))
GEANTDICTS    = G__Geant.cxx
GEANTO        = $(patsubst %.cxx,src/%.o,$(GEANTSOURCES)) G__Geant.o
GEANTDEPS     = $(patsubst %.cxx,$(INCDIR)/%.h,$(GEANTSOURCES)) $(patsubst %.cxx,%.h,$(GEANTDICTS))
GEANTSO       = libGeant.$(DllSuf)
GEANTLIB      = $(shell pwd)/$(GEANTSO)

USERSOURCES   = MyHit.cxx MyApplication.cxx
USERH         = MyHit.h MyApplication.h UserLinkDef.h
USERHI        = $(patsubst %.h,inc/%.h,$(USERH))
USERDICTS     = G__User.cxx
USERO        = $(patsubst %.cxx,src/%.o,$(USERSOURCES)) G__User.o
USERDEPS      = $(patsubst %.cxx,$(INCDIR)/%.h,$(USERSOURCES)) $(patsubst %.cxx,%.h,$(USERDICTS))
USERSO        = libUser.$(DllSuf)
USERLIB       = $(shell pwd)/$(USERSO)

OBJS          = $(GEANTSO) $(USERSO)


#------------------------------------------------------------------------------

all:            $(GEANTSO) $(USERSO)
geant:          $(GEANTSO)
user:           $(USERSO)

$(GEANTSO):     $(GEANTDICTS) $(GEANTO)
		$(LD) $(LDFLAGS) $(SOFLAGS) $(GEANTO) $(OutPutOpt)$@
		@echo "$@ done"
$(USERSO):      $(USERDICTS) $(USERO)
		$(LD) $(LDFLAGS) $(SOFLAGS) $(USERO) $(OutPutOpt)$@
		@echo "$@ done"

src/%.o: src/%.cxx inc/%.h
		$(CXX)  $(CXXFLAGS) -o $@ -c $<

$(GEANTDICTS): $(GEANTHI) 
		@echo "Generating dictionary $@"
		$(ROOTSYS)/bin/rootcint -f $@ -c -p $(GEANTHI)
$(USERDICTS):
		@echo "Generating dictionary $@"
		$(ROOTSYS)/bin/rootcint -f $@ -c $(USERHI)
G__Geant.o:     $(GEANTDICTS)
#		$(ROOTSYS)/bin/rmkdepend -R -f G__Geant.d -Y -w 1000 -- \
                pipe -m64 -Wshadow -Wall -W -Woverloaded-virtual -fPIC \
                -I$(INCDIR) -pthread  -D__cplusplus -I$(ROOTSYS)cint/cint/lib/prec_stl \
                -I$(ROOTSYS)cint/cint/stl -I$(ROOTSYS)/cint/cint/inc -- $(GEANTDICTS)
		$(CXX)  $(CXXFLAGS) -I. -I$(ROOTSYS)/include -o $@ -c G__Geant.cxx
G__User.o:      $(USERDICTS)
#		$(ROOTSYS)/bin/rmkdepend -R -f G__User.d -Y -w 1000 -- \
                pipe -m64 -Wshadow -Wall -W -Woverloaded-virtual -fPIC \
                -I$(INCDIR) -pthread  -D__cplusplus -I$(ROOTSYS)cint/cint/lib/prec_stl \
                -I$(ROOTSYS)cint/cint/stl -I$(ROOTSYS)/cint/cint/inc -- $(USERDICTS)
		$(CXX)  $(CXXFLAGS) -I. -I$(ROOTSYS)/include -o $@ -c G__User.cxx
clean:
		@rm -rf $(GEANTSO) $(USERSO) src/*.o G__*
                

