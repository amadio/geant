#ifndef SHAPE_PERF_TEST
#define SHAPE_PERF_TEST

#include "TGeoBBox.h"
//#include "TGeoMatrix.h"
#include "TGeoMatrixPlain.h"
#include <vector>
#include "tbb/tick_count.h" // timing from Intel TBB 
#include <cassert>

struct StopWatch 
{
  tbb::tick_count t1;
  tbb::tick_count t2;
  void Start(){ t1=tbb::tick_count::now(); }
  void Stop(){ t2=tbb::tick_count::now(); }
  void Reset(){ /* */ ;}
  //  void Print(){  std::cerr << (t2 - t1).seconds() << std::endl; }
  double getDeltaSecs() { return (t2-t1).seconds(); }
};

#define N 14
//template<typename ShapeT>
class TrafoBenchmarker
{
 private:
  StopWatch timer;

  TGeoBBox * testshape;
  TGeoMatrix const* trafo;

  unsigned int vecsizes[14]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8182};
  static const unsigned int MAXSIZE=8182;
  static const unsigned int NREPS = 10000; // number of timing repetitions;

  // containers for the timings
  std::vector<double> TdO; // dist out
  std::vector<double> Tt; // matrix trafo

  // containers for testdata ( only one vector prepared for each function from which we will resample )
  double *points; // distanceIn
  double *dirs; // directions distanceIn

 public:
  // ctor
 TrafoBenchmarker( TGeoBBox *s, TGeoMatrix const *m ) : testshape(s), trafo(m), TdO(N,0.), Tt(N,0.)
    {
      points = new double[3*MAXSIZE];
      dirs = new double[3*MAXSIZE];
      // calculate Bounding Box for shape
      testshape->ComputeBBox();
    }

  ~TrafoBenchmarker()
    {
      delete[] points;
      delete[] dirs;
    }

  // these functions prepare the test data such that we test all possible cases
  // for example: distanceFromOut should have 1/3. of particles actually hitting the object etc.
  void initData();

  // actual timing functions which call the shapes routine
  double timeTrafo( double &, unsigned int );
  void timeDistanceFromOutside( double &, unsigned int );
  void timeIt();
};

#endif // SHAPE_PERF_TEST
