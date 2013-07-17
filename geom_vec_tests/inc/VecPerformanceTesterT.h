#ifndef SHAPE_PERF_TEST_V
#define SHAPE_PERF_TEST_V

#include "PerformanceTesterT.h"

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
  ShapeBenchmarker_v( T  *s ) : ShapeBenchmarker<T>(s), TdO_v(NS,0.), TdI_v(NS,0.), Tc_v(NS,0.), Ts_v(NS,0.) {
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

  // heatup ( to avoid overhead from first call to dynamic library )
  void heatup();

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

  // make contact to dynamic library
  heatup();

  // do additional vector interface measurements
  // to avoid measuring the same function over and over again we interleave calls to different functions and different data
  for(unsigned int rep=0; rep< ShapeBenchmarker<T>::NREPS; rep++)
    {
      for(unsigned int vectype =0 ; vectype < NS; ++vectype )
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
  for(unsigned int i =0 ; i < this->MAXSIZE; ++i )
    {
      if( results_v_C[i] != this->results_C[i] )
	{
	  std::cerr << i << " #difference in contains function" << std::endl;
	}
      if( fabs(results_v_s[i] - this->results_s[i]) > 1e-10 )
	{
	  //	  VUSolid::EnumInside inside = testusolid->U::Inside( &this->points_C[3*index] );
	  std::cerr << i << " #difference in safety function " << results_v_s[i] << "  " <<  this->results_s[i] << " CONTAINS: " << results_v_C[i] << " " << this->results_C[i] << std::endl;
	}
      if( ( fabs(results_v_dO[i] - this->results_dO[i]) > 1e-10 ) && ( fabs(results_v_dO[i] - this->results_dO[i]) < TGeoShape::Big() ) )
	{
	  std::cerr << i << " #difference in DistFromOut function " << results_v_dO[i] << " " <<  this->results_dO[i] << std::endl;
	}
      if( ( fabs(results_v_dI[i] - this->results_dI[i]) > 1e-10 ) && ( fabs(results_v_dI[i] - this->results_dI[i]) < TGeoShape::Big() ) ) 
	{
	  std::cerr << i << " #difference in DistFromIn function " <<  results_v_dI[i] << " " <<  this->results_dI[i] << std::endl;
	}
    }

  /*
  for(unsigned int i =0 ; i < this->MAXSIZE; ++i )
    {
      // std::cerr << results_v_C[i] << " " << this->results_C[i] << std::endl;
      //      std::cerr << results_v_s[i] - this->results_s[i] << std::endl;
      //      std::cerr << results_v_dO[i] - this->results_v_dO[i] << std::endl;
      //      std::cerr << results_v_dI[i] << " " << this->results_dI[i] << " " << results_v_dI[i] - this->results_dI[i] << std::endl;
    }
  */

  this->correctTimingAndNormalize(Tc_v);
  this->correctTimingAndNormalize(Ts_v);
  this->correctTimingAndNormalize(TdI_v);
  this->correctTimingAndNormalize(TdO_v);

  // print result
  for(unsigned int vectype =0 ; vectype < NS; ++vectype )
    {
      std::cout << ShapeBenchmarker<T>::vecsizes[vectype] 
		<< " " << Tc_v[vectype]  /* timing for Contains method */
		<< " " << Tc_v[0]/(Tc_v[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) /* speedup with respect to 1 particle */
		<< " " <<  Ts_v[vectype]   /* timing for safety method */
		<< " " << Ts_v[0]/(Ts_v[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) 
		<< " " <<  TdI_v[vectype] 
		<< " " << TdI_v[0]/(TdI_v[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) 
		<< " " <<  TdO_v[vectype] 
		<< " " << TdO_v[0]/(TdO_v[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) 
		<< std::endl;
    }  
}

template <typename T>
void ShapeBenchmarker_v<T>::timeSafety_v(double & Tacc, unsigned int vecsize)
{
  this->timer.HeatUp();

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
  this->timer.HeatUp();

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
  this->timer.HeatUp();

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
  this->timer.HeatUp();

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

template <typename T>
void ShapeBenchmarker_v<T>::heatup()
{
  this->testshape->T::DistFromInside_v( points_dI_SOA, dirs_dI_SOA, 3, this->steps, 0, results_v_dI, this->MAXSIZE );
  this->testshape->T::DistFromOutside_v( points_dO_SOA, dirs_dO_SOA, 3, this->steps, 0, results_v_dO, this->MAXSIZE );
  this->testshape->T::Contains_v( points_C_SOA, results_v_C, this->MAXSIZE );
  this->testshape->T::Safety_v( points_s_SOA, kTRUE, results_v_s, this->MAXSIZE );
}


#endif
