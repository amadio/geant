# relatice include dir
INC=../inc/
VCINC=${VCROOT}include
VCLIB=${VCROOT}lib/libVc.a

ifeq ($(CXX),g++)
OPT=-std=c++11 -fabi-version=6 -O3 -mavx -mfpmath=sse -funroll-loops -ftree-vectorize -flax-vector-conversions -fivopts  \
-ffast-math -funsafe-math-optimizations -msse4.1 -Wattributes -ftree-vectorizer-verbose=7 -fPIC -I ${INC} -I ${VCINC} # -DVEC_EXTENSIONS
LDOPTS=-shared -Wl,--hash-style=sysv -Wl,-Bsymbolic
endif

# ifeq ($(CXX),g++)
# OPT=-ggdb -std=c++11 -fabi-version=6 -O2 -mavx -mfpmath=sse -funroll-loops -ftree-vectorize -flax-vector-conversions -fivopts  \
# -ffast-math -funsafe-math-optimizations -msse4.1 -Wattributes -ftree-vectorizer-verbose=7 -fPIC -I ${INC} -I ${VCINC} # -DVEC_EXTENSIONS
# LDOPTS=-shared -Wl,--hash-style=sysv -Wl,-Bsymbolic
# endif

ifeq ($(CXX),clang++) 
OPT=-O3 -mllvm -unroll-allow-partial -mllvm -unroll-runtime \
-mllvm -vectorize -mllvm -bb-vectorize-aligned-only -I ${INC}
endif

ifeq ($(CXX),icc) 
OPT=-O3 -DHAVE_AVX -mkl -mavx -xAVX -Ofast -vec-report3 -restrict -fPIC -I ${INC}
LDOPTS=-shared -Wl,--hash-style=sysv -Wl,-Bsymbolic
endif


ROOTLIBS:=$(shell root-config --libs)
ROOTCFLAGS:=$(shell root-config --cflags)

tests := $(shell ls TGeo*_Test*.cxx | xargs -I {} basename {} .cxx | xargs -I {} echo {}_$(CXX))

looplines:=$(shell grep -n "EXPECT" TGeoBBox_v.cxx | sed 's/:.*//')
gitcommit:=$(shell git reflog | head -n 1 | awk '{print $1}')
compilerversion:=$(shell $(CXX) --version)


TGeoBBox_$(CXX)_loopvecreport.txt: libtest_$(CXX).so
	echo $< > $@
	echo "GIT COMMIT $(gitcommit)" >> $@ 
	echo "COMPILER VERSION $(compilerversion)" >> $@ 
	echo "COMPILER FLAGS $(OPT)" >> $@ 
	grep "LOOP" compilelog_$(CXX) >> $@
	echo "EXPECTED $(looplines)" >> $@


all: libtest_$(CXX).so

test : $(tests)

impvect: libtest_$(CXX).so simpvect.o
	$(CXX) -o simpvect simpvect.o $(ROOTLIBS) -lGeom -L. -ltest_$(CXX) -o $@

simpvect.o: simpvect.cxx
	$(CXX) -m64 $(OPT) -c simpvect.cxx $(ROOTCFLAGS) 

libtest_$(CXX).so: TGeoBBox_v_$(CXX).o # testCint.o
	$(CXX) -m64 $(LDOPTS)  $(ROOTLIBS) -lGeom -o $@ $^ ${VCLIB}

#testCint.o: testCint.cxx
#	$(CXX) -m64 $(OPT) -c -I $(ROOTSYS)/include $*.cxx

#testCint.cxx: LinkDef.h TGeoBBox_v.h
#	rootcint -f testCint.cxx -c TGeoBBox_v.h LinkDef.h

TGeoBBox_v_$(CXX).o: TGeoBBox_v.cxx TGeoBBox_v.h
	$(CXX) -m64 $(OPT) -c -I $(ROOTSYS)/include $< -o $@ 2> compilelog_$(CXX)

TGeoBBox_Test_%_$(CXX): TGeoBBox_Test_%.o libtest_$(CXX).so
	$(CXX) -o $@ $<  $(ROOTLIBS) ${VCLIB} -lGeom -L$(TBBROOT)/lib -ltbb -L. -ltest_$(CXX) 

TGeoBBox_Test_%.o: TGeoBBox_Test_%.cxx
	$(CXX) -m64 $(OPT) $(ROOTCFLAGS) -I$(TBBROOT)/include -c $<

TGeoCone_Test_%_$(CXX): TGeoCone_Test_%.o
	$(CXX) -o $@ $<  $(ROOTLIBS) ${VCLIB} -lGeom -L$(TBBROOT)/lib -ltbb 

TGeoCone_Test_%.o: TGeoCone_Test_%.cxx
	$(CXX) -m64 $(OPT) $(ROOTCFLAGS) -I$(TBBROOT)/include -c $<

clean:
	rm TGeoBBox_v_$(CXX).o libtest_$(CXX).so # testCint.* 

cleantest:
	rm $(tests)