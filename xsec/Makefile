# Simple Makefile to build "stand alone" this lib

ifeq ($(shell uname),Darwin)
  SO:=so
  SOFLAGS:=-dynamiclib -single_module
endif

CXXFLAGS := $(shell root-config --cflags) -I$(PWD)/inc

OBJS := $(subst .cxx,.o,$(wildcard $(PWD)/src/*.cxx))

%.o: %.cxx
	$(CXX) -c -o $@ $(CXXFLAGS)  $^

CWD := $(shell pwd)

all: $(CWD)/lib/libXsec.$(SO)

$(CWD)/lib/libXsec.$(SO): $(OBJS) $(CWD)/src/xsecCint.o
	if [ ! -r $(CWD)/lib ] ; then mkdir $(CWD)/lib; fi
	echo $^
	$(CXX) -o $@ $(SOFLAGS) $(shell root-config --libs) $^

$(CWD)/src/xsecCint.cxx: $(filter-out %xsecLinkDef.h,$(wildcard $(CWD)/inc/*.h))
	rootcint -f $@ -c $(CXXFLAGS) $^ $(CWD)/inc/xsecLinkDef.h

clean: 
	rm -rf $(OBJS) $(CWD)/src/xsecCint.* lib
