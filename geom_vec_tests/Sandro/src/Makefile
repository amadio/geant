#

#gcc: -std=c99 -O3 -funroll-loops -vopts -
#ax-vector-conversions
#-ast-math -funsafe-math-optimizations -msse4.1 (with
#-fno-tree-vectorize to turn o autovectorization)
#clang: -O3 -mllvm -unroll-allow-partial -mllvm -unroll-runtime
#-funsafe-math-optimizations -ast-math (with -mllvm -vectorize
#-mllvm -bb-vectorize-aligned-only for autovectorization)

#OPT=-O3 -mllvm -unroll-allow-partial -mllvm -unroll-runtime \
-mllvm -vectorize -mllvm -bb-vectorize-aligned-only
#CC=clang++

ifeq ($(CXX),g++) 
OPT=-O3 -mfpmath=sse -mmmx -msse -msse2 -funroll-loops -flax-vector-conversions -fivopts  \
-ffast-math -funsafe-math-optimizations -msse4.1 -Wattributes -ftree-vectorizer-verbose=3 -fPIC
LDOPTS=-shared -Wl,--hash-style=sysv -Wl,-Bsymbolic
endif

ifeq ($(CXX),clang++) 
OPT=-O3 -mllvm -unroll-allow-partial -mllvm -unroll-runtime \
-mllvm -vectorize -mllvm -bb-vectorize-aligned-only
endif

ifeq ($(CXX),icc) 
OPT=-O3 -DHAVE_AVX -mkl -mavx -vec-report3 -restrict -fPIC
LDOPTS=-shared -Wl,--hash-style=sysv -Wl,-Bsymbolic
endif


ROOTLIBS:=$(shell root-config --libs)
ROOTCFLAGS:=$(shell root-config --cflags)

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
	$(CXX) -m64 $(OPT) -c -I $(ROOTSYS)/include $*.cxx

TGeoBBox_Test_Safety: TGeoBBox_Test_Safety.o
	$(CXX) -m64 $(OPT) $(ROOTCFLAGS) -I$(TBBROOT)/include -c TGeoBBox_Test_Safety.cxx
	$(CXX) -o TGeoBBox_Test_Safety TGeoBBox_Test_Safety.o  $(ROOTLIBS) -lGeom -L$(TBBROOT)/lib -ltbb -L. -ltest

TGeoBBox_Test_Safety.o: TGeoBBox_Test_Safety.cxx
	$(CXX) -m64 $(OPT) $(ROOTCFLAGS) -I$(TBBROOT)/include -c TGeoBBox_Test_Safety.cxx

clean:
	rm TGeoBBox_v.o libtest.so testCint.* 