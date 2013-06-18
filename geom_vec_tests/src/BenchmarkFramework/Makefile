ifeq ($(CXX),g++)
OPT=-O3 -mavx -mfpmath=sse -funroll-loops -ftree-vectorize -flax-vector-conversions -fivopts  \
-ffast-math -funsafe-math-optimizations -msse4.1 -Wattributes -ftree-vectorizer-verbose=7 -fPIC
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

tests := $(shell ls BenchTGeo*.cxx | xargs -I {} basename {} .cxx | xargs -I {} echo {}_$(CXX))

test : $(tests) PerformanceTester.o

PerformanceTester.o: PerformanceTester.cxx
	$(CXX) -m64 $(OPT) $(ROOTCFLAGS) -I$(TBBROOT)/include -c $<


BenchTGeo%_$(CXX):  BenchTGeo%.o PerformanceTester.o
	$(CXX) -o $@ $< PerformanceTester.o $(ROOTLIBS) -lGeom -L$(TBBROOT)/lib -ltbb 

BenchTGeo%.o: BenchTGeo%.cxx
	$(CXX) -m64 $(OPT) $(ROOTCFLAGS) -c $<

clean:
	rm *.o 

cleantest:
	rm $(tests)