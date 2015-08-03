CXX           = icc
ifeq ($(CXX),g++) 
CXXFLAGS      = -O3 -fPIC -m64 -pthread -march=amdfam10 -msse -msse2 -msse3 -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=2
endif

ifeq ($(CXX),icc) 
CXXFLAGS      = -O3 -fPIC -m64 -wd1476 -wd1572 -wd279 -pthread -vec-report2
endif

LIBS          := $(shell root-config --libs)

vectors:   vectors.o Makefile
		$(CXX) -o $@ vectors.o $(LIBS)
vectors.o: vectors.C Makefile
		$(CXX) -I$(ROOTSYS)/include -c vectors.C $(CXXFLAGS)
clean:
		rm vectors vectors.o