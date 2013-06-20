#include "PerformanceTesterTGeoMatrix.hpp"
#include "TRandom.h"
//#include "TMath.h"
#include <iostream>
#include <cmath>
#include "Util.h"


void TrafoBenchmarker::initData()
{
  // store box sizes instead of doing many virtual function calls
  double dx = testshape->GetDX();
  double dy = testshape->GetDX();
  double dz = testshape->GetDX();

  // initialize the data first off all within twice bounding box and arbitrary direction
  for( unsigned int i=0; i<MAXSIZE; ++i) 
    {
      Util::samplePoint(&points[3*i],dx,dy,dz,5);
      Util::sampleDir(&dirs[3*i]);
    }

  int insidecounter=0;

  // now try to do some adjustments ( we want all points outside and a fraction of 1./3. of hitting particles )
  for(int i=0; i< MAXSIZE; ++i) 
    {
      insidecounter+=testshape->Contains( &points[3*i] );
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;

  // move all particles outside
  while(insidecounter > 0)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( testshape->Contains( &points[ 3*index ] ))
	{
	  do{
	    Util::samplePoint(&points[3*index],dx,dy,dz,5);
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are outside
	    }while( testshape->Contains( &points[ 3*index ] ) );
	  insidecounter--;
	}	
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
  
  // now check how many particles hit this shape
  int hitcounter=0;
  for(unsigned int i=0; i< MAXSIZE; ++i) 
    {
      hitcounter+=testshape->CouldBeCrossed( &points[3*i], &dirs[3*i] );
    }

  std::cerr << " have " << hitcounter << " points hitting " << std::endl;

  while(hitcounter < MAXSIZE/3.)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( ! testshape->CouldBeCrossed( &points[ 3*index ], &dirs[3*index] ))
	{
	  do{
	    Util::sampleDir( &dirs[3*index] );
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are outside
	  }while( ! testshape->CouldBeCrossed( &points[ 3*index ], &dirs[3*index]) );
	  hitcounter++;
	}	
    }
    std::cerr << " have " << hitcounter << " points hitting " << std::endl;
}



// elemental function do one time measurement for a given number of points
void TrafoBenchmarker::timeDistanceFromOutside(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*MAXSIZE-1.*vecsize);
  volatile double distance;
  if(vecsize==1)
    {
      timer.Start();
      distance=testshape->DistFromOutside( &points[3*startindex], &dirs[3*startindex], 3, TGeoShape::Big(), 0);
      timer.Stop();
    }
  else
    {
      timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  distance=testshape->DistFromOutside( &points[3*index], &dirs[3*index], 3, TGeoShape::Big(), 0);
	}
      timer.Stop();
    }
  Tacc+=timer.getDeltaSecs();
}

double TrafoBenchmarker::timeTrafo( double & Tacc, unsigned int vecsize )
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*MAXSIZE-1.*vecsize);
  double local[3];
  double accum=0.;
  if(vecsize==1)
    {
      timer.Start();
      trafo->MasterToLocal(&points[3*startindex], local );
      timer.Stop();
    }
  else
    {
      timer.Start();
#pragma ivdep
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  trafo->MasterToLocal(&points[3*index], local );
	}
      timer.Stop();
      accum+=local[0]+local[1]+local[2];
    }
  Tacc+=timer.getDeltaSecs();
  // if(vecsize==8182) std::cerr << accum << std::endl;
}

void TrafoBenchmarker::timeIt( )
{
  // init all data
  initData();

  // to avoid measuring the same function over and over again we interleave calls to different functions and different data
  for(unsigned int rep=0; rep< NREPS; rep++)
    {
      for(unsigned int vectype =0 ; vectype < N; ++vectype )
	{
	  // distFromOutside
	  timeDistanceFromOutside ( TdO[vectype], vecsizes[vectype] );

	  // time trafo
	  timeTrafo( Tt[vectype], vecsizes[vectype] );
	}
    }

  // print result
  for(unsigned int vectype =0 ; vectype < N; ++vectype )
    {
      std::cerr << vecsizes[vectype] 
		<< " " <<  TdO[vectype]/NREPS 
		<< " " << TdO[0]/(TdO[vectype]/vecsizes[vectype]) 
		<< " " <<  Tt[vectype]/NREPS 
		<< " " << Tt[0]/(Tt[vectype]/vecsizes[vectype]) 
		<< std::endl;
    }
} 
