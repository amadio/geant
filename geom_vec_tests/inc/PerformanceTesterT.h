// a set of generic classes that to obtain benchmark information
// on the most important TGeo methods for both old and (possible new) vector interfaces
// author: Sandro Wenzel, CERN, June 2013

#ifndef SHAPE_PERF_TEST
#define SHAPE_PERF_TEST

#include "TGeoBBox.h"
#include <vector>
#include "TBBStopWatch.h" // timing from Intel TBB 
#include "RDTSCStopWatch.h" // RDTSC timer
#include <iostream>
#include "TRandom.h"
#include "PointStruct.h" // for the SOA data container
#include <cmath>
#include "Util.h" // for sampling points and dirs
#include "PointStruct.h"
#include "mm_malloc.h" // for aligned malloc
#include <fstream>

#define NS 14


// T is supposed to be a TGeoShape inherited from TGeoBBox
template <typename T>
class ShapeBenchmarker{
 protected:
#ifdef USE_RDTSC // to choose timer at compile time
  RDTSCStopWatch timer;
#else
  StopWatch timer;
#endif

  // for the timer overhead
  double Toverhead;

  T * testshape;
  unsigned int vecsizes[14]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8182};
  static const unsigned int MAXSIZE=8182;
  static const unsigned int NREPS = 1000; // number of timing repetitions;

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

  void heatup();

 public:
  // ctor
  ShapeBenchmarker( T *s ) : testshape(s), TdO(NS,0.), TdI(NS,0.), Tc(NS,0.), Ts(NS,0.) {
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

  // get timer overhead 
  Toverhead=timer.getOverhead( 10000 );
  
  std::cerr << "# TIMER OVERHEAD IS :" << Toverhead << std::endl;

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

  void printTimings( char const * filename ) const;


  // correct times by overhead and normalize
  // we expect here the timings to be accumulated over NREP repetitions
  void correctTimingAndNormalize( std::vector<double> & timing )
  {
    for(unsigned int vectype =0 ; vectype < NS; ++vectype )
      {
	// take out overhead
	timing[ vectype ] -= Toverhead*NREPS;
	// normalize
	timing[ vectype ] /= NREPS;
      }
  }



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
      Util::samplePoint(&points_C[3*i],dx,dy,dz, testshape->T::GetOrigin(), 1);
    }

  int insidecounter=0;

  // now try to do some resampling
  for(int i=0; i< MAXSIZE; ++i) 
    {
      insidecounter+=testshape->T::Contains( &points_C[3*i] );
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
  while(insidecounter < MAXSIZE/3. )
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( ! testshape->T::Contains( &points_C[ 3*index ] ))
	{
	  do{
	    Util::samplePoint(&points_C[3*index],dx,dy,dz, testshape->T::GetOrigin(), 1);
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are inside
	  } while( ! testshape->T::Contains( &points_C[ 3*index ] ) );
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
      Util::samplePoint(&points_s[3*i],dx,dy,dz,testshape->T::GetOrigin(),1);
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
	    Util::samplePoint(&points_s[3*index],dx,dy,dz, testshape->T::GetOrigin(), 1);
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
      Util::samplePoint(&points_dI[3*i],dx,dy,dz,testshape->T::GetOrigin(),1);
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
	    Util::samplePoint(&points_dI[3*index],dx,dy,dz,testshape->T::GetOrigin(),1);
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
      Util::samplePoint(&points_dO[3*i],dx,dy,dz,testshape->T::GetOrigin(),5);
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
	    Util::samplePoint(&points_dO[3*index],dx,dy,dz, testshape->T::GetOrigin(), 5);
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
  timer.HeatUp();

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
  // heat up the timer/caches
  timer.HeatUp();

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
  timer.HeatUp();

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
  timer.HeatUp();

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

  heatup();

  // to avoid measuring the same function over and over again we interleave calls to different functions and different data
  for(unsigned int rep=0; rep< NREPS; rep++)
    {
      for(unsigned int vectype =0 ; vectype < NS; ++vectype )
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

  correctTimingAndNormalize(Tc);
  correctTimingAndNormalize(Ts);
  correctTimingAndNormalize(TdI);
  correctTimingAndNormalize(TdO);
  
  // print result
  for(unsigned int vectype =0 ; vectype < NS; ++vectype )
    {
      std::cout << vecsizes[vectype] 
		<< " " << Tc[vectype]  /* timing for Contains method */

		<< " " << Tc[0]/(Tc[vectype]/vecsizes[vectype]) /* speedup with respect to 1 particle */

		<< " " <<  Ts[vectype]  /* timing for safety method */

		<< " " << Ts[0]/(Ts[vectype]/vecsizes[vectype]) 

		<< " " <<  TdI[vectype]

		<< " " << TdI[0]/(TdI[vectype]/vecsizes[vectype]) 

		<< " " <<  TdO[vectype]

		<< " " << TdO[0]/(TdO[vectype]/vecsizes[vectype]) 
		<< std::endl;
    }  
} 

template <typename T>
void ShapeBenchmarker<T>::printTimings( char const * filename ) const
{
  ofstream outstr(filename);
  // print result
  for(unsigned int vectype =0 ; vectype < NS; ++vectype )
    {
      outstr << this->vecsizes[vectype] 
		<< " " << this->Tc[vectype]  /* timing for Contains method */
		<< " " << this->Tc[0]/(Tc[vectype]/vecsizes[vectype]) /* speedup with respect to 1 particle */
		<< " " << this->Ts[vectype]   /* timing for safety method */
		<< " " << this->Ts[0]/(Ts[vectype]/vecsizes[vectype]) 
		<< " " <<  this->TdI[vectype] 
		<< " " << this->TdI[0]/(TdI[vectype]/vecsizes[vectype]) 
		<< " " <<  this->TdO[vectype] 
		<< " " << this->TdO[0]/(TdO[vectype]/vecsizes[vectype]) 
		<< std::endl;
    }  
  outstr.close();
}


template <typename T>
void ShapeBenchmarker<T>::heatup()
{
  results_dI[0] = this->testshape->T::DistFromInside( &points_dI[3*0], &dirs_dI[3*0], 3, TGeoShape::Big(), 0 ); 
  results_dO[0] = this->testshape->T::DistFromOutside( &points_dO[3*0], &dirs_dO[3*0], 3, TGeoShape::Big(), 0 ); 
  results_C[0] = this->testshape->T::Contains( points_C );
  results_s[0] = this->testshape->T::Safety( points_s, kTRUE );
}


#endif // SHAPE_PERF_TEST
