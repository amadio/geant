#ifndef SHAPE_PERF_TEST_V
#define SHAPE_PERF_TEST_V

#include "PerformanceTesterT.h"

// for the SolidConversion
#include "TGeoToUSolidConvertor.h"

// 
#include "UVector3.hh"
#include "VUSolid.hh"

// extension to benchmark and check USolid classes
template <typename U, typename T>
class USolidBenchmarker : ShapeBenchmarker<T> 
{
private:
  U * testusolid;

 // additional containers for the timings
  std::vector<double> TdO_USolid; // dist out
  std::vector<double> TdI_USolid; // dist in
  std::vector<double> Tc_USolid;  // contains
  std::vector<double> Ts_USolid;  // safety

  double *results_USolid_dI; // distance In
  double *results_USolid_dO; // distance from Out
  double *results_USolid_s;  // safety
  bool *results_USolid_C;  // contains

public:
 USolidBenchmarker( T *s ) : ShapeBenchmarker<T>(s), TdO_USolid(N,0.), TdI_USolid(N,0.), Tc_USolid(N,0.), Ts_USolid(N,0.) 
    {
      results_USolid_dI=(double *)_mm_malloc(sizeof(double)*this->MAXSIZE,32);
      results_USolid_dO=(double *)_mm_malloc(sizeof(double)*this->MAXSIZE,32);
      results_USolid_s=(double *)_mm_malloc(sizeof(double)*this->MAXSIZE,32);
      results_USolid_C=(bool *)_mm_malloc(sizeof(bool)*this->MAXSIZE,32);

      // get the usolid from the TGeoShape ( T )
      testusolid = (U*) TGeoToUSolidConvertor::Convert( this->testshape ); 

      if( testusolid == 0 )
	{
	  std::cerr << " SHAPE CONVERSION FAILED " << std::endl;
	}
    }

  // actual timing functions which call the shapes routine
  void timeDistanceFromInside_USolid( double &, unsigned int );
  void timeDistanceFromOutside_USolid( double &, unsigned int );
  void timeContains_USolid( double &, unsigned int );
  void timeSafety_USolid( double &, unsigned int );
  void timeIt();
};

template<typename U, typename T> 
  void USolidBenchmarker<U,T>::timeIt()
{
  // call Base class timeIt()
  this->ShapeBenchmarker<T>::timeIt();

  // do additional vector interface measurements
  // to avoid measuring the same function over and over again we interleave calls to different functions and different data
  for(unsigned int rep=0; rep< ShapeBenchmarker<T>::NREPS; rep++)
    {
      for(unsigned int vectype =0 ; vectype < N; ++vectype )
	{
	  // Safety
	  timeSafety_USolid ( Ts_USolid[vectype], ShapeBenchmarker<T>::vecsizes[vectype] );

	  // Contains
	  timeContains_USolid ( Tc_USolid[vectype], ShapeBenchmarker<T>::vecsizes[vectype] );

	  // distFromInside
	  timeDistanceFromInside_USolid ( TdI_USolid[vectype], ShapeBenchmarker<T>::vecsizes[vectype] );

	  // distFromOutside
	  timeDistanceFromOutside_USolid ( TdO_USolid[vectype], ShapeBenchmarker<T>::vecsizes[vectype] );
	}
    }

  // compare results ( good place for tests and asserts )
  for(unsigned int i =0 ; i < this->MAXSIZE; ++i )
    {
      if( results_USolid_C[i] != this->results_C[i] )
	{
	  std::cerr << i << " #difference in contains function" << std::endl;
	}
      if( fabs(results_USolid_s[i] - this->results_s[i]) > 1e-10 )
	{
	  //	  VUSolid::EnumInside inside = testusolid->U::Inside( &this->points_C[3*index] );
	  std::cerr << i << " #difference in safety function " << results_USolid_s[i] << "  " <<  this->results_s[i] << " CONTAINS: " << results_USolid_C[i] << " " << this->results_C[i] << std::endl;
	}
      if( ( fabs(results_USolid_dO[i] - this->results_dO[i]) > 1e-10 ) && ( fabs(results_USolid_dO[i] - this->results_dO[i]) < TGeoShape::Big() ) )
	{
	  std::cerr << i << " #difference in DistFromOut function " << results_USolid_dO[i] << " " <<  this->results_dO[i] << std::endl;
	}
      if( ( fabs(results_USolid_dI[i] - this->results_dI[i]) > 1e-10 ) && ( fabs(results_USolid_dI[i] - this->results_dI[i]) < TGeoShape::Big() ) ) 
	{
	  std::cerr << i << " #difference in DistFromIn function " <<  results_USolid_dI[i] << " " <<  this->results_dI[i] << std::endl;
	}

      //       std::cerr << results_USolid_C[i] << " " << this->results_C[i] << std::endl;
      //       std::cerr << results_USolid_s[i] - this->results_s[i] << std::endl;
      //       std::cerr << results_USolid_dO[i] - this->results_USolid_dO[i] << std::endl;
      //       std::cerr << results_USolid_dI[i] << " " << this->results_USolid_dI[i] << " " << results_USolid_dI[i] - this->results_USolid_dI[i] << std::endl;
    }
  
  this->correctTimingAndNormalize(Tc_USolid);
  this->correctTimingAndNormalize(Ts_USolid);
  this->correctTimingAndNormalize(TdI_USolid);
  this->correctTimingAndNormalize(TdO_USolid);

  // print result
  for(unsigned int vectype =0 ; vectype < N; ++vectype )
    {
      std::cout << ShapeBenchmarker<T>::vecsizes[vectype] 
		<< " " << Tc_USolid[vectype]  /* timing for Contains method */
		<< " " << Tc_USolid[0]/(Tc_USolid[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) /* speedup with respect to 1 particle */
		<< " " <<  Ts_USolid[vectype]   /* timing for safety method */
		<< " " << Ts_USolid[0]/(Ts_USolid[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) 
		<< " " <<  TdI_USolid[vectype] 
		<< " " << TdI_USolid[0]/(TdI_USolid[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) 
		<< " " <<  TdO_USolid[vectype] 
		<< " " << TdO_USolid[0]/(TdO_USolid[vectype]/ShapeBenchmarker<T>::vecsizes[vectype]) 
		<< std::endl;
    }  
}

template<typename U, typename T> 
  void USolidBenchmarker<U,T>::timeSafety_USolid(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*this->MAXSIZE-1.*vecsize);
  if(vecsize==1)
    {
      this->timer.Start();
      results_USolid_s[startindex]=testusolid->U::SafetyFromInside( &this->points_s[3*startindex], true );
      this->timer.Stop();
    }
  else
    {
      this->timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  results_USolid_s[index]=testusolid->U::SafetyFromInside( &this->points_s[3*index], true ); // true for precise mode
	}
      this->timer.Stop();
    }
  Tacc+=this->timer.getDeltaSecs();
}


template<typename U, typename T> 
  void USolidBenchmarker<U,T>::timeContains_USolid(double & Tacc, unsigned int vecsize)
{
  this->timer.HeatUp();

  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*this->MAXSIZE-1.*vecsize);
  if(vecsize==1)
    {
      this->timer.Start();
      VUSolid::EnumInside inside = testusolid->U::Inside( &this->points_C[3*startindex] );
      this->timer.Stop();
      results_USolid_C[startindex] = (inside == VUSolid::eOutside)? false : true;
    }
  else
    {
      this->timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  VUSolid::EnumInside inside = testusolid->U::Inside( &this->points_C[3*index] );
	  results_USolid_C[index] = (inside == VUSolid::eOutside)? false : true;
	}
      this->timer.Stop();
    }
  Tacc+=this->timer.getDeltaSecs();
}


template<typename U, typename T> 
  void USolidBenchmarker<U,T>::timeDistanceFromOutside_USolid(double & Tacc, unsigned int vecsize)
{
  this->timer.HeatUp();

  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*this->MAXSIZE-1.*vecsize);
  volatile double distance;
  if(vecsize==1)
    {
      this->timer.Start();

      results_USolid_dO[startindex]=testusolid->U::DistanceToIn( &this->points_dO[3*startindex], &this->dirs_dO[3*startindex], TGeoShape::Big());
      this->timer.Stop();
    }
  else
    {
      this->timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  results_USolid_dO[index]=testusolid->U::DistanceToIn( &this->points_dO[3*index], &this->dirs_dO[3*index], TGeoShape::Big());
	}
      this->timer.Stop();
    }
  Tacc+=this->timer.getDeltaSecs();
}


template <typename U, typename T>
  void USolidBenchmarker<U,T>::timeDistanceFromInside_USolid(double & Tacc, unsigned int vecsize)
{
  this->timer.HeatUp();

  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*this->MAXSIZE-1.*vecsize);
  if(vecsize==1)
    {
      this->timer.Start();
      UVector3 normal;
      bool convex;   
      results_USolid_dI[startindex]=testusolid->U::DistanceToOut( &this->points_dI[3*startindex], &this->dirs_dI[3*startindex], normal, convex, TGeoShape::Big() );
      this->timer.Stop();
    }
  else
    {
      this->timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  UVector3 normal;
	  bool convex;   
	  results_USolid_dI[index]=testusolid->U::DistanceToOut( &this->points_dI[3*index], &this->dirs_dI[3*index], normal, convex, TGeoShape::Big() );
	}
      this->timer.Stop();
    }
  Tacc+=this->timer.getDeltaSecs();
}

#endif
