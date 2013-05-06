//
// tests and benchmarks for simple vectorized math functions
//
#include "VecUtils.h"
#include "TStopwatch.h" // timer from root
#include "TRandom.h" // random num from root
#include <iostream>
#include "tbb/tick_count.h"" // timing from Intel TBB

main(int argc, char *argv[])
{
  int npoints=100000000;

  double *x=new double[npoints];
  double *y=new double[npoints];
  double *z=new double[npoints];

  // fill the input vectors
  for(int i=0; i<npoints; ++i) {
    x[i]=(1-2*gRandom->Rndm());
    y[i]=(1-2*gRandom->Rndm());
  }

  tbb::tick_count ts;
  tbb::tick_count te;

  // testing vectorized min version
  /*
  tbb::tick_count ts = tbb::tick_count::now();
  VecUtils::min_l( x, y, z, npoints );
  tbb::tick_count te = tbb::tick_count::now();
  std::cerr << (te - ts).seconds() << std::endl;
  */

  // testing loop/ blocked min
  // tt.Start();

  /*
  ts=tbb::tick_count::now();
  VecUtils::min_v( x, y, z, npoints );
  te=tbb::tick_count::now();
  std::cerr << (te - ts).seconds() << std::endl;
  */

  //  tt.Stop();
  //  tt.Print();
  //  tt.Reset();

  // testing serial min
  // testing loop/ blocked min
  //  tt.Start();
  
  ts=tbb::tick_count::now();
  for(int i=0;i<npoints;i++)
    {
      z[i]=VecUtils::min( x[i], y[i]);
    }
  te=tbb::tick_count::now();
  std::cerr << (te - ts).seconds() << std::endl;
  
  //  tt.Stop();
  //  tt.Print();
  //  tt.Reset();

  ts=tbb::tick_count::now();
  double d=-1.23;
  for(unsigned int i= 1;i<10000000;i++)
    {
      volatile double v = VecUtils::abs2( d );
      d*=-1.02;
    }
  te=tbb::tick_count::now();
  std::cerr << (te - ts).seconds() << std::endl;

  ts=tbb::tick_count::now();
  d=-1.23;
  for(unsigned int i= 1;i<10000000;i++)
    {
      volatile double v = VecUtils::abs1( d );
      d*=-1.02;
    }
  te=tbb::tick_count::now();
  std::cerr << (te - ts).seconds() << std::endl;

  delete[] x;
  delete[] y;
  delete[] z;
 
 return 0;
}
