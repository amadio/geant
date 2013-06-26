#ifndef SHAPE_PERF_TEST
#define SHAPE_PERF_TEST

#include "TGeoBBox.h"
#include <vector>
#include "tbb/tick_count.h" // timing from Intel TBB 
#include <cassert>
#include "Util.h"

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
class ShapeBenchmarker
{
 private:
  StopWatch timer;

  TGeoBBox * testshape;
  unsigned int vecsizes[14]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8182};
  static const unsigned int MAXSIZE=8182;
  static const unsigned int NREPS = 2000; // number of timing repetitions;

  // containers for the timings
  std::vector<double> TdO; // dist out
  std::vector<double> TdI; // dist in
  std::vector<double> Tc;  //  contains
  std::vector<double> Ts;  // safety

  // containers for testdata ( only one vector prepared for each function from which we will resample )
  double *points_C; // contains
  double *points_dO; // distanceIn
  double *points_dI; // distanceOut
  double *points_s; // safety

  double *dirs_dO; // directions distanceIn
  double *dirs_dI; // directions distanceOut
 
  double *steps; // steps for distanceIn and distanceOut


 
 public:
  // ctor
  ShapeBenchmarker( TGeoBBox *s );

  ~ShapeBenchmarker()
    {
      delete[] points_C;
      delete[] points_dO;
      delete[] points_dI;
      delete[] points_s;
      delete[] dirs_dO;
      delete[] dirs_dI;
      delete[] steps;
    }

  // these functions prepare the test data such that we test all possible cases
  // for example: distanceFromOut should have 1/3. of particles actually hitting the object etc.
  void initDataDistanceFromInside();
  void initDataDistanceFromOutside();
  void initDataSafety();
  void initDataContains();

  // actual timing functions which call the shapes routine
  void timeDistanceFromInside( double &, unsigned int );
  void timeDistanceFromOutside( double &, unsigned int );
  void timeContains( double &, unsigned int );
  void timeSafety( double &, unsigned int );
  void timeIt();
};


#endif // SHAPE_PERF_TEST
