
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


GEANTSOURCES  = GeantThreadData.cxx GeantTrack.cxx GeantBasket.cxx  sync_objects.cxx \
                GeantOutput.cxx GeantPropagator.cxx \
                PhysicsProcess.cxx WorkloadManager.cxx GeantVApplication.cxx GeantFactoryStore.cxx
USERSOURCES   = MyHit.cxx MyApplication.cxx
GEANTDICTS    = G__Geant.cxx
USERDICTS     = G__User.cxx
GEANTO        = GeantThreadData.o GeantTrack.o GeantBasket.o sync_objects.o \
                PhysicsProcess.o GeantOutput.o GeantPropagator.o \
                WorkloadManager.o GeantVApplication.o G__Geant.o
USERO         = MyHit.o  MyApplication.o  G__User.o
GEANTSL       = $(patsubst %.o,$(SRCDIR)/%.o,$(GEANTO))
USERSL        = $(patsubst %.o,$(SRCDIR)/%.o,$(USERO))
GEANTDEPS     = $(patsubst %.cxx,$(INCDIR)/%.h,$(GEANTSOURCES))
USERDEPS      = $(patsubst %.cxx,$(INCDIR)/%.h,$(USERSOURCES))

GEANTSO       = libGeant.$(DllSuf)
USERSO        = libUser.$(DllSuf)
GEANTLIB      = $(shell pwd)/$(GEANTSO)
USERLIB       = $(shell pwd)/$(USERSO)
OBJS          = $(GEANTSO) $(USERSO)


#------------------------------------------------------------------------------

all:            $(GEANTSO) $(USERSO)
geant:          $(GEANTSO)
user:           $(USERSO)

$(GEANTSO):     $(GEANTDICTS) $(GEANTO)
		$(LD) $(LDFLAGS) $(SOFLAGS) $(GEANTSL) $(OutPutOpt)$@
		@echo "$@ done"
$(USERSO):      $(USERDICTS) $(USERO)
		$(LD) $(LDFLAGS) $(SOFLAGS) $(USERSL) $(OutPutOpt)$@
		@echo "$@ done"
GeantThreadData.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantThreadData.cxx
PhysicsProcess.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/PhysicsProcess.cxx
GeantTrack.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantTrack.cxx
GeantBasket.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantBasket.cxx
GeantOutput.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantOutput.cxx
GeantPropagator.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantPropagator.cxx
WorkloadManager.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/WorkloadManager.cxx
GeantVApplication.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantVApplication.cxx
GeantFactoryStore.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/GeantFactoryStore.cxx
sync_objects.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/sync_objects.cxx
MyHit.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/MyHit.cxx
MyApplication.o:
		$(CXX)  $(CXXFLAGS) -o src/$@ -c src/MyApplication.cxx
$(GEANTDICTS):
		@echo "Generating dictionary $@"
		$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/$@ \
                -c -p $(INCDIR)/GeantThreadData.h $(INCDIR)/PhysicsProcess.h \
                $(INCDIR)/GeantBasket.h $(INCDIR)/GeantScheduler.h \
                $(INCDIR)/GeantPropagator.h $(INCDIR)/GeantTrack.h $(INCDIR)/GeantOutput.h \
                $(INCDIR)/WorkloadManager.h $(INCDIR)/GeantVApplication.h $(INCDIR)/LinkDef.h
$(USERDICTS):
		@echo "Generating dictionary $@"
		$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/$@ \
                -c  $(INCDIR)/MyHit.h $(INCDIR)/MyApplication.h \
                $(INCDIR)/UserLinkDef.h
G__Geant.o:
		$(ROOTSYS)/bin/rmkdepend -R -f$(SRCDIR)/G__Geant.d -Y -w 1000 -- \
                pipe -m64 -Wshadow -Wall -W -Woverloaded-virtual -fPIC \
                -I$(INCDIR) -pthread  -D__cplusplus -I$(ROOTSYS)cint/cint/lib/prec_stl \
                -I$(ROOTSYS)cint/cint/stl -I$(ROOTSYS)/cint/cint/inc -- $(SRCDIR)/$(GEANTDICTS)
		$(CXX)  $(CXXFLAGS) -I. -I$(ROOTSYS)/include -o $(SRCDIR)/$@ \
                -c $(SRCDIR)/G__Geant.cxx
G__User.o:
		$(ROOTSYS)/bin/rmkdepend -R -f$(SRCDIR)/G__User.d -Y -w 1000 -- \
                pipe -m64 -Wshadow -Wall -W -Woverloaded-virtual -fPIC \
                -I$(INCDIR) -pthread  -D__cplusplus -I$(ROOTSYS)cint/cint/lib/prec_stl \
                -I$(ROOTSYS)cint/cint/stl -I$(ROOTSYS)/cint/cint/inc -- $(SRCDIR)/$(USERDICTS)
		$(CXX)  $(CXXFLAGS) -I. -I$(ROOTSYS)/include -o $(SRCDIR)/$@ \
                -c $(SRCDIR)/G__User.cxx
clean:
		@rm -rf $(GEANTSO) $(USERSO) src/*.o src/*.d src/G__*
                

