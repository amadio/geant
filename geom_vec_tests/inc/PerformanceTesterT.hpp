// a set of generic classes that to obtain benchmark information
// on the most important TGeo methods for both old and (possible new) vector interfaces
// author: Sandro Wenzel, CERN, June 2013

#ifndef SHAPE_PERF_TEST
#define SHAPE_PERF_TEST

#include "TGeoBBox.h"
#include <vector>
#include "tbb/tick_count.h" // timing from Intel TBB 
#include <iostream>
#include "TRandom.h"
#include "PointStruct.h" // for the SOA data container
#include <cmath>
#include "Util.h" // for sampling points and dirs
#include "PointStruct.h"
#include "mm_malloc.h" // for aligned malloc

#define N 14

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


// T is supposed to be a TGeoShape inherited from TGeoBBox
template <typename T>
class ShapeBenchmarker{
 protected:
  StopWatch timer;

  T * testshape;
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

  // containers for results
  double *results_dI; // distance In
  double *results_dO; // distance from Out
  double *results_s;  // safety
  bool *results_C;  // contains

 public:
  // ctor
  ShapeBenchmarker( T *s ) : testshape(s), TdO(N,0.), TdI(N,0.), Tc(N,0.), Ts(N,0.) {
  points_C = (double*) _mm_malloc(sizeof(double)*3*MAXSIZE,32);
  points_dO = (double*) _mm_malloc(sizeof(double)*3*MAXSIZE,32);
  points_dI = (double*) _mm_malloc(sizeof(double)*3*MAXSIZE,32);
  points_s = (double*) _mm_malloc(sizeof(double)*3*MAXSIZE,32);
  dirs_dO = (double*) _mm_malloc(sizeof(double)*3*MAXSIZE,32);
  dirs_dI = (double*) _mm_malloc(sizeof(double)*3*MAXSIZE,32);
  steps = (double*) _mm_malloc(sizeof(double)*MAXSIZE,32);
  results_dI = (double*) _mm_malloc(sizeof(double)*MAXSIZE,32);
  results_dO = (double*) _mm_malloc(sizeof(double)*MAXSIZE,32);
  results_s = (double*) _mm_malloc(sizeof(double)*MAXSIZE,32);
  results_C = (bool*) _mm_malloc(sizeof(bool)*MAXSIZE,32);

  // need to fill steps
  for(unsigned int i=0;i<MAXSIZE;i++)
    {
      steps[i]=TGeoShape::Big();
    }

  // calculate Bounding Box for shape
  testshape->T::ComputeBBox();
  } 

  ~ShapeBenchmarker()
    {
      _mm_free(points_C);
      _mm_free(points_dO);
      _mm_free(points_dI);
      _mm_free(points_s);
      _mm_free(dirs_dO);
      _mm_free(dirs_dI);
      _mm_free(steps);
      _mm_free(results_dI);
      _mm_free(results_dO);
      _mm_free(results_s);
      _mm_free(results_C);
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

template <typename T>
void ShapeBenchmarker<T>::initDataContains()
{
  // store box sizes instead of doing many virtual function calls
  double dx = testshape->T::GetDX();
  double dy = testshape->T::GetDX();
  double dz = testshape->T::GetDX();

  // initialize the data first off all within twice bounding box
  for( unsigned int i=0; i<MAXSIZE; ++i) 
    {
      Util::samplePoint(&points_C[3*i],dx,dy,dz,2);
    }

  int insidecounter=0;

  // now try to do some resampling
  for(int i=0; i< MAXSIZE; ++i) 
    {
      insidecounter+=testshape->T::Contains( &points_C[3*i] );
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
  while(insidecounter < MAXSIZE/3.)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( ! testshape->T::Contains( &points_C[ 3*index ] ))
	{
	  do{
	    Util::samplePoint(&points_C[3*index],dx,dy,dz, 2);
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are inside
	  }while( ! testshape->T::Contains( &points_C[ 3*index ] ) );
	  insidecounter++;
	}	
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
}

template <typename T>
void ShapeBenchmarker<T>::initDataSafety()
{
  // store box sizes instead of doing many virtual function calls
  double dx = testshape->T::GetDX();
  double dy = testshape->T::GetDX();
  double dz = testshape->T::GetDX();

  // initialize the data first off all within twice bounding box
  for( unsigned int i=0; i<MAXSIZE; ++i) 
    {
      Util::samplePoint(&points_s[3*i],dx,dy,dz,1);
    }

  int insidecounter=0;

  // now try to do some resampling
  for(int i=0; i< MAXSIZE; ++i) 
    {
      insidecounter+=testshape->T::Contains( &points_s[3*i] );
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
  while(insidecounter < MAXSIZE)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( ! testshape->T::Contains( &points_s[ 3*index ] ))
	{
	  do{
	    Util::samplePoint(&points_s[3*index],dx,dy,dz,1);
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are inside
	    }while( ! testshape->T::Contains( &points_s[ 3*index ] ) );
	}	
      insidecounter++;
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
}


template <typename T>
void ShapeBenchmarker<T>::initDataDistanceFromInside()
{
  // store box sizes instead of doing many virtual function calls
  double dx = testshape->T::GetDX();
  double dy = testshape->T::GetDX();
  double dz = testshape->T::GetDX();

  // initialize the data first off all within twice bounding box
  for( unsigned int i=0; i<MAXSIZE; ++i) 
    {
      Util::samplePoint(&points_dI[3*i],dx,dy,dz,1);
      Util::sampleDir(&dirs_dI[3*i]);
    }

  int insidecounter=0;

  // now try to do some resampling
  for(int i=0; i< MAXSIZE; ++i) 
    {
      insidecounter+=testshape->T::Contains( &points_dI[3*i] );
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
  while(insidecounter < MAXSIZE)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( ! testshape->T::Contains( &points_dI[ 3*index ] ))
	{
	  do{
	    Util::samplePoint(&points_dI[3*index],dx,dy,dz,1);
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are inside
	    }while( ! testshape->T::Contains( &points_dI[ 3*index ] ) );
	  insidecounter++;
      	}	
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
}

template <typename T>
void ShapeBenchmarker<T>::initDataDistanceFromOutside()
{
  // store box sizes instead of doing many virtual function calls
  double dx = testshape->T::GetDX();
  double dy = testshape->T::GetDX();
  double dz = testshape->T::GetDX();

  // initialize the data first off all within twice bounding box and arbitrary direction
  for( unsigned int i=0; i<MAXSIZE; ++i) 
    {
      Util::samplePoint(&points_dO[3*i],dx,dy,dz,5);
      Util::sampleDir(&dirs_dO[3*i]);
    }

  int insidecounter=0;

  // now try to do some adjustments ( we want all points outside and a fraction of 1./3. of hitting particles )
  for(int i=0; i< MAXSIZE; ++i) 
    {
      insidecounter+=testshape->T::Contains( &points_dO[3*i] );
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;

  // move all particles outside
  while(insidecounter > 0)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( testshape->T::Contains( &points_dO[ 3*index ] ))
	{
	  do{
	    Util::samplePoint(&points_dO[3*index],dx,dy,dz,5);
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are outside
	    }while( testshape->T::Contains( &points_dO[ 3*index ] ) );
	  insidecounter--;
	}	
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
  
  // now check how many particles hit this shape
  int hitcounter=0;
  for(unsigned int i=0; i< MAXSIZE; ++i) 
    {
      hitcounter+=testshape->T::CouldBeCrossed( &points_dO[3*i], &dirs_dO[3*i] );
    }

  std::cerr << " have " << hitcounter << " points hitting " << std::endl;

  while(hitcounter < MAXSIZE/3.)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( ! testshape->T::CouldBeCrossed( &points_dO[ 3*index ], &dirs_dO[3*index] ))
	{
	  do{
	    Util::sampleDir( &dirs_dO[3*index] );
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are outside
	  }while( ! testshape->T::CouldBeCrossed( &points_dO[ 3*index ], &dirs_dO[3*index]) );
	  hitcounter++;
	}	
    }
    std::cerr << " have " << hitcounter << " points hitting " << std::endl;
}



// elemental function do one time measurement for a given number of points
template <typename T>
void ShapeBenchmarker<T>::timeContains(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*MAXSIZE-1.*vecsize);
  if(vecsize==1)
    {
      timer.Start();
      results_C[startindex]=testshape->T::Contains( &points_C[3*startindex] );
      timer.Stop();
    }
  else
    {
      timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  results_C[index]=testshape->T::Contains( &points_C[3*index] );
	}
      timer.Stop();
    }
  Tacc+=timer.getDeltaSecs();
}



// elemental function do one time measurement for a given number of points
template <typename T>
void ShapeBenchmarker<T>::timeSafety(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*MAXSIZE-1.*vecsize);
  if(vecsize==1)
    {
      timer.Start();
      results_s[startindex]=testshape->T::Safety( &points_s[3*startindex], kTRUE );
      timer.Stop();
    }
  else
    {
      timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  results_s[index]=testshape->T::Safety( &points_s[3*index], kTRUE );
	}
      timer.Stop();
    }
  Tacc+=timer.getDeltaSecs();
}



// elemental function do one time measurement for a given number of points
template <typename T>
void ShapeBenchmarker<T>::timeDistanceFromInside(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*MAXSIZE-1.*vecsize);
  volatile double distance;
  if(vecsize==1)
    {
      timer.Start();
      results_dI[startindex]=testshape->T::DistFromInside( &points_dI[3*startindex], &dirs_dI[3*startindex], 3, TGeoShape::Big(), 0);
      timer.Stop();
    }
  else
    {
      timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  results_dI[index]=testshape->T::DistFromInside( &points_dI[3*index], &dirs_dI[3*index], 3, TGeoShape::Big(), 0);
	}
      timer.Stop();
    }
  Tacc+=timer.getDeltaSecs();
}



// elemental function do one time measurement for a given number of points
template <typename T>
void ShapeBenchmarker<T>::timeDistanceFromOutside(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*MAXSIZE-1.*vecsize);
  if(vecsize==1)
    {
      timer.Start();
      results_dO[startindex]=testshape->T::DistFromOutside( &points_dO[3*startindex], &dirs_dO[3*startindex], 3, TGeoShape::Big(), 0);
      timer.Stop();
    }
  else
    {
      timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  results_dO[index]=testshape->T::DistFromOutside( &points_dO[3*index], &dirs_dO[3*index], 3, TGeoShape::Big(), 0);
	}
      timer.Stop();
    }
  Tacc+=timer.getDeltaSecs();
}


template <typename T>
void ShapeBenchmarker<T>::timeIt( )
{
  // init all data
  initDataContains();
  initDataSafety();
  initDataDistanceFromInside();
  initDataDistanceFromOutside();

  // to avoid measuring the same function over and over again we interleave calls to different functions and different data
  for(unsigned int rep=0; rep< NREPS; rep++)
    {
      for(unsigned int vectype =0 ; vectype < N; ++vectype )
	{
	  // Safety
	  timeSafety ( Ts[vectype], vecsizes[vectype] );

	  // Contains
	  timeContains( Tc[vectype], vecsizes[vectype] );

	  // distFromInside
	  timeDistanceFromInside ( TdI[vectype], vecsizes[vectype] );

	  // distFromOutside
	  timeDistanceFromOutside ( TdO[vectype], vecsizes[vectype] );
	}
    }
  
  // print result
  for(unsigned int vectype =0 ; vectype < N; ++vectype )
    {
      std::cerr << vecsizes[vectype] 
		<< " " << Tc[vectype]/NREPS  /* timing for Contains method */
		<< " " << Tc[0]/(Tc[vectype]/vecsizes[vectype]) /* speedup with respect to 1 particle */
		<< " " <<  Ts[vectype]/NREPS   /* timing for safety method */
		<< " " << Ts[0]/(Ts[vectype]/vecsizes[vectype]) 
		<< " " <<  TdI[vectype]/NREPS 
		<< " " << TdI[0]/(TdI[vectype]/vecsizes[vectype]) 
		<< " " <<  TdO[vectype]/NREPS 
		<< " " << TdO[0]/(TdO[vectype]/vecsizes[vectype]) 
		<< std::endl;
    }  
} 

// extension to benchmark vectorized interfaces
template <typename T>
class ShapeBenchmarker_v : ShapeBenchmarker<T> 
{
private:
 // additional containers for the timings
  std::vector<double> TdO_v; // dist out
  std::vector<double> TdI_v; // dist in
  std::vector<double> Tc_v;  //  contains
  std::vector<double> Ts_v;  // safety

  double *results_v_dI; // distance In
  double *results_v_dO; // distance from Out
  double *results_v_s;  // safety
  bool *results_v_C;  // contains

  // data in SOA form
  StructOfCoord points_C_SOA; // contains
  StructOfCoord points_dO_SOA; // distanceIn
  StructOfCoord points_dI_SOA; // distanceOut
  StructOfCoord points_s_SOA; // safety
  StructOfCoord dirs_dO_SOA; // directions distanceIn
  StructOfCoord dirs_dI_SOA; // directions distanceOut


public:
  ShapeBenchmarker_v( T  *s ) : ShapeBenchmarker<T>(s), TdO_v(N,0.), TdI_v(N,0.), Tc_v(N,0.), Ts_v(N,0.) {
    points_C_SOA.alloc(ShapeBenchmarker<T>::MAXSIZE); // contains
    points_dO_SOA.alloc(ShapeBenchmarker<T>::MAXSIZE); // distanceIn
    points_dI_SOA.alloc(ShapeBenchmarker<T>::MAXSIZE); // distanceOut
    points_s_SOA.alloc(ShapeBenchmarker<T>::MAXSIZE); // safety
    
    dirs_dO_SOA.alloc(ShapeBenchmarker<T>::MAXSIZE); // directions distanceIn
    dirs_dI_SOA.alloc(ShapeBenchmarker<T>::MAXSIZE); // directions distanceOut
    
    results_v_dI=(double *)_mm_malloc(sizeof(double)*this->MAXSIZE,32);
    results_v_dO=(double *)_mm_malloc(sizeof(double)*this->MAXSIZE,32);
    results_v_s=(double *)_mm_malloc(sizeof(double)*this->MAXSIZE,32);
    results_v_C=(bool *)_mm_malloc(sizeof(bool)*this->MAXSIZE,32);
  }

  // actual timing functions which call the shapes routine
  void timeDistanceFromInside_v( double &, unsigned int );
  void timeDistanceFromOutside_v( double &, unsigned int );
  void timeContains_v( double &, unsigned int );
  void timeSafety_v( double &, unsigned int );
  void timeIt();
};

template<typename T> 
void ShapeBenchmarker_v<T>::timeIt()
{
  // call Base class timeIt()
  this->ShapeBenchmarker<T>::timeIt();

  // need to bring original data in SOA Form();
  points_C_SOA.fill(ShapeBenchmarker<T>::points_C); // contains
  points_dO_SOA.fill(ShapeBenchmarker<T>::points_dO); // distanceIn
  points_dI_SOA.fill(ShapeBenchmarker<T>::points_dI); // distanceOut
  points_s_SOA.fill(ShapeBenchmarker<T>::points_s); // safety
  
  dirs_dO_SOA.fill(ShapeBenchmarker<T>::dirs_dO); // directions distanceIn
  dirs_dI_SOA.fill(ShapeBenchmarker<T>::dirs_dI); // directions distanceOut

  // do additional vector interface measurements
  // to avoid measuring the same function over and over again we interleave calls to different functions and different data
  for(unsigned int rep=0; rep< ShapeBenchmarker<T>::NREPS; rep++)
    {
      for(unsigned int vectype =0 ; vectype < N; ++vectype )
	{
	  // Safety
	  timeSafety_v ( Ts_v[vectype], ShapeBenchmarker<T>::vecsizes[vectype] );

	  // Contains
	  timeContains_v ( Tc_v[vectype], ShapeBenchmarker<T>::vecsizes[vectype] );

	  // distFromInside
	  timeDistanceFromInside_v ( TdI_v[vectype], ShapeBenchmarker<T>::vecsizes[vectype] );

	  // distFromOutside
	  timeDistanceFromOutside_v ( TdO_v[vectype], ShapeBenchmarker<T>::vecsizes[vectype] );
	}
    }

  // compare results ( good place for tests and asserts )
  /*
  for(unsigned int i =0 ; i < this->MAXSIZE; ++i )
    {
      // std::cerr << results_v_C[i] << " " << this->results_C[i] << std::endl;
      //      std::cerr << results_v_s[i] - this->results_s[i] << std::endl;
      //      std::cerr << results_v_dO[i] - this->results_v_dO[i] << std::endl;
      //      std::cerr << results_v_dI[i] << " " << this->results_dI[i] << " " << results_v_dI[i] - this->results_dI[i] << std::endl;
    }
  */

  // print result
  for(unsigned int vectype =0 ; vectype < N; ++vectype )
    {
      std::cerr << ShapeBenchmarker<T>::vecsizes[vectype] 
		<< " " << Tc_v[vectype]/ShapeBenchmarker<T>::NREPS  /* timing for Contains method */
		<< " " << Tc_v[0]/(Tc_v[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) /* speedup with respect to 1 particle */
		<< " " <<  Ts_v[vectype]/ShapeBenchmarker<T>::NREPS   /* timing for safety method */
		<< " " << Ts_v[0]/(Ts_v[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) 
		<< " " <<  TdI_v[vectype]/ShapeBenchmarker<T>::NREPS 
		<< " " << TdI_v[0]/(TdI_v[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) 
		<< " " <<  TdO_v[vectype]/ShapeBenchmarker<T>::NREPS 
		<< " " << TdO_v[0]/(TdO_v[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) 
		<< std::endl;
    }  
}

template <typename T>
void ShapeBenchmarker_v<T>::timeSafety_v(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*ShapeBenchmarker<T>::MAXSIZE-1.*vecsize);
  while(startindex % 4 !=0) startindex--; // for alignment reasons
  points_s_SOA.setstartindex(startindex);
  this->timer.Start();
     
  // we should call the SOA version here
  // so we pass a reference to Struct of Point
  this->testshape->T::Safety_v( points_s_SOA, kTRUE, results_v_s, vecsize );
  this->timer.Stop();
  Tacc+=this->timer.getDeltaSecs();
}


template <typename T>
void ShapeBenchmarker_v<T>::timeContains_v(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*ShapeBenchmarker<T>::MAXSIZE-1.*vecsize);
  while(startindex % 4 !=0) startindex--; // for alignment reasons
  points_C_SOA.setstartindex(startindex);
  this->timer.Start();
  this->testshape->T::Contains_v( points_C_SOA, results_v_C, vecsize );
  this->timer.Stop();
  Tacc+=this->timer.getDeltaSecs();
}


template <typename T>
void ShapeBenchmarker_v<T>::timeDistanceFromOutside_v(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*ShapeBenchmarker<T>::MAXSIZE-1.*vecsize);
  while(startindex % 4 !=0) startindex--; // for alignment reasons
  points_dO_SOA.setstartindex(startindex);
  dirs_dO_SOA.setstartindex(startindex);
  this->timer.Start();
  this->testshape->T::DistFromOutside_v( points_dO_SOA, dirs_dO_SOA, 3, this->steps, 0, results_v_dO, vecsize );
  this->timer.Stop();
  Tacc+=this->timer.getDeltaSecs();
}


template <typename T>
void ShapeBenchmarker_v<T>::timeDistanceFromInside_v(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*ShapeBenchmarker<T>::MAXSIZE-1.*vecsize);
  while(startindex % 4 !=0) startindex--; // for alignment reasons  
  points_dI_SOA.setstartindex(startindex);
  dirs_dI_SOA.setstartindex(startindex);
  
  this->timer.Start();
  this->testshape->T::DistFromInside_v( points_dI_SOA, dirs_dI_SOA, 3, this->steps, 0, results_v_dI, vecsize );
  this->timer.Stop();
  Tacc+=this->timer.getDeltaSecs();
}


#endif // SHAPE_PERF_TEST
