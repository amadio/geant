ifeq ($(CXX),g++) 
OPT=-O3 -mfpmath=sse -mavx -funroll-loops -flax-vector-conversions -fivopts  \
-ffast-math -funsafe-math-optimizations -msse4.1 -Wattributes -ftree-vectorizer-verbose=3 -fPIC
LDOPTS=-shared -Wl,--hash-style=sysv -Wl,-Bsymbolic
endif

ifeq ($(CXX),clang++) 
OPT=-O3 -mllvm -unroll-allow-partial -mllvm -unroll-runtime \
-mllvm -vectorize -mllvm -bb-vectorize-aligned-only
endif

ifeq ($(CXX),icc) 
OPT=-O3 -DHAVE_AVX -mkl -mavx -xAVX -Ofast -vec-report3 -restrict -fPIC
LDOPTS=-shared -Wl,--hash-style=sysv -Wl,-Bsymbolic
endif


ROOTLIBS:=$(shell root-config --libs)
ROOTCFLAGS:=$(shell root-config --cflags)

tests := $(shell ls TGeoBBox_Test*.cxx | xargs -I {} basename {} .cxx)

looplines:=$(shell grep -n "EXPECT" TGeoBBox_v.cxx | sed 's/:.*//')
gitcommit:=$(shell git reflog | head -n 1 | awk '{print $1}')
compilerversion:=$(shell $(CXX) --version)


TGeoBBox_loopvecreport.txt: TGeoBBox_v.o
	echo $< > $@
	echo "GIT COMMIT $(gitcommit)" >> $@ 
	echo "COMPILER VERSION $(compilerversion)" >> $@ 
	echo "COMPILER FLAGS $(OPT)" >> $@ 
	grep "LOOP" compilelog >> $@
	echo "EXPECTED $(looplines)" >> $@

all: libtest.so

test : $(tests)

impvect: libtest.so simpvect.o
	$(CXX) -o simpvect simpvect.o $(ROOTLIBS) -lGeom -L. -ltest -o $@

simpvect.o: simpvect.cxx
	$(CXX) -m64 $(OPT) -c simpvect.cxx $(ROOTCFLAGS) 

libtest.so: TGeoBBox_v.o testCint.o
	$(CXX) -m64 $(LDOPTS) $(ROOTLIBS) -lGeom -o $@ $^

testCint.o: testCint.cxx
	$(CXX) -m64 $(OPT) -c -I $(ROOTSYS)/include $*.cxx

testCint.cxx: LinkDef.h TGeoBBox_v.h
	rootcint -f testCint.cxx -c TGeoBBox_v.h LinkDef.h

TGeoBBox_v.o: TGeoBBox_v.cxx TGeoBBox_v.h
	$(CXX) -m64 $(OPT) -c -I $(ROOTSYS)/include $< 2> compilelog

TGeoBBox_Test_%: TGeoBBox_Test_%.o
	$(CXX) -o $@ $<  $(ROOTLIBS) -lGeom -L$(TBBROOT)/lib -ltbb -L. -ltest

TGeoBBox_Test_%.o: TGeoBBox_Test_%.cxx
	$(CXX) -m64 $(OPT) $(ROOTCFLAGS) -I$(TBBROOT)/include -c $<

clean:
	rm TGeoBBox_v.o libtest.so testCint.* 

cleantest:
	rm $(tests)