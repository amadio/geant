#include "PerformanceTester.h"
#include "TRandom.h"
//#include "TMath.h"
#include <iostream>
#include <cmath>
#include "Util.h"



ShapeBenchmarker::ShapeBenchmarker( TGeoBBox *s ) : testshape(s), TdO(N,0.), TdI(N,0.), Tc(N,0.), Ts(N,0.) {
  points_C = new double[3*MAXSIZE];
  points_dO = new double[3*MAXSIZE];
  points_dI = new double[3*MAXSIZE];
  points_s = new double[3*MAXSIZE];
  dirs_dO = new double[3*MAXSIZE];
  dirs_dI = new double[3*MAXSIZE];
  steps = new double[MAXSIZE];

  // calculate Bounding Box for shape
  testshape->ComputeBBox();
} 


void ShapeBenchmarker::initDataContains()
{
  // store box sizes instead of doing many virtual function calls
  double dx = testshape->GetDX();
  double dy = testshape->GetDX();
  double dz = testshape->GetDX();

  // initialize the data first off all within twice bounding box
  for( unsigned int i=0; i<MAXSIZE; ++i) 
    {
      Util::samplePoint(&points_C[3*i],dx,dy,dz,2);
    }

  int insidecounter=0;

  // now try to do some resampling
  for(int i=0; i< MAXSIZE; ++i) 
    {
      insidecounter+=testshape->Contains( &points_C[3*i] );
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
  while(insidecounter < MAXSIZE/3.)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( ! testshape->Contains( &points_C[ 3*index ] ))
	{
	  do{
	    Util::samplePoint(&points_C[3*index],dx,dy,dz, 2);
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are inside
	    }while( ! testshape->Contains( &points_C[ 3*index ] ) );
	  insidecounter++;
	}	
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
}


void ShapeBenchmarker::initDataSafety()
{
  // store box sizes instead of doing many virtual function calls
  double dx = testshape->GetDX();
  double dy = testshape->GetDX();
  double dz = testshape->GetDX();

  // initialize the data first off all within twice bounding box
  for( unsigned int i=0; i<MAXSIZE; ++i) 
    {
      Util::samplePoint(&points_s[3*i],dx,dy,dz,1);
    }

  int insidecounter=0;

  // now try to do some resampling
  for(int i=0; i< MAXSIZE; ++i) 
    {
      insidecounter+=testshape->Contains( &points_s[3*i] );
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
  while(insidecounter < MAXSIZE)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( ! testshape->Contains( &points_s[ 3*index ] ))
	{
	  do{
	    Util::samplePoint(&points_s[3*index],dx,dy,dz,1);
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are inside
	    }while( ! testshape->Contains( &points_s[ 3*index ] ) );
	}	
      insidecounter++;
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
}



void ShapeBenchmarker::initDataDistanceFromInside()
{
  // store box sizes instead of doing many virtual function calls
  double dx = testshape->GetDX();
  double dy = testshape->GetDX();
  double dz = testshape->GetDX();

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
      insidecounter+=testshape->Contains( &points_dI[3*i] );
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
  while(insidecounter < MAXSIZE)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( ! testshape->Contains( &points_dI[ 3*index ] ))
	{
	  do{
	    Util::samplePoint(&points_dI[3*index],dx,dy,dz,1);
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are inside
	    }while( ! testshape->Contains( &points_dI[ 3*index ] ) );
	  insidecounter++;
      	}	
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
}


void ShapeBenchmarker::initDataDistanceFromOutside()
{
  // store box sizes instead of doing many virtual function calls
  double dx = testshape->GetDX();
  double dy = testshape->GetDX();
  double dz = testshape->GetDX();

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
      insidecounter+=testshape->Contains( &points_dO[3*i] );
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;

  // move all particles outside
  while(insidecounter > 0)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( testshape->Contains( &points_dO[ 3*index ] ))
	{
	  do{
	    Util::samplePoint(&points_dO[3*index],dx,dy,dz,5);
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are outside
	    }while( testshape->Contains( &points_dO[ 3*index ] ) );
	  insidecounter--;
	}	
    }
  std::cerr << " prepared data with " << insidecounter << " points inside " << std::endl;
  
  // now check how many particles hit this shape
  int hitcounter=0;
  for(unsigned int i=0; i< MAXSIZE; ++i) 
    {
      hitcounter+=testshape->CouldBeCrossed( &points_dO[3*i], &dirs_dO[3*i] );
    }

  std::cerr << " have " << hitcounter << " points hitting " << std::endl;

  while(hitcounter < MAXSIZE/3.)
    {
      // pick a point randomly
      int index = gRandom->Rndm()*MAXSIZE;
      if( ! testshape->CouldBeCrossed( &points_dO[ 3*index ], &dirs_dO[3*index] ))
	{
	  do{
	    Util::sampleDir( &dirs_dO[3*index] );
	      // this might be totally inefficient and we should replace this by some importance sampling or cloning of points which are outside
	  }while( ! testshape->CouldBeCrossed( &points_dO[ 3*index ], &dirs_dO[3*index]) );
	  hitcounter++;
	}	
    }
    std::cerr << " have " << hitcounter << " points hitting " << std::endl;
}



// elemental function do one time measurement for a given number of points
void ShapeBenchmarker::timeContains(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*MAXSIZE-1.*vecsize);
  volatile bool isin;
  if(vecsize==1)
    {
      timer.Start();
      isin=testshape->Contains( &points_C[3*startindex] );
      timer.Stop();
    }
  else
    {
      timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  isin=testshape->Contains( &points_C[3*index] );
	}
      timer.Stop();
    }
  Tacc+=timer.getDeltaSecs();
}



// elemental function do one time measurement for a given number of points
void ShapeBenchmarker::timeSafety(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*MAXSIZE-1.*vecsize);
  volatile double safety;
  if(vecsize==1)
    {
      timer.Start();
      safety=testshape->Safety( &points_s[3*startindex], kTRUE );
      timer.Stop();
    }
  else
    {
      timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  safety=testshape->Safety( &points_s[3*index], kTRUE );
	}
      timer.Stop();
    }
  Tacc+=timer.getDeltaSecs();
}



// elemental function do one time measurement for a given number of points
void ShapeBenchmarker::timeDistanceFromInside(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*MAXSIZE-1.*vecsize);
  volatile double distance;
  if(vecsize==1)
    {
      timer.Start();
      distance=testshape->DistFromInside( &points_dI[3*startindex], &dirs_dI[3*startindex], 3, TGeoShape::Big(), 0);
      timer.Stop();
    }
  else
    {
      timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  distance=testshape->DistFromInside( &points_dI[index], &dirs_dI[3*index], 3, TGeoShape::Big(), 0);
	}
      timer.Stop();
    }
  Tacc+=timer.getDeltaSecs();
}



// elemental function do one time measurement for a given number of points
void ShapeBenchmarker::timeDistanceFromOutside(double & Tacc, unsigned int vecsize)
{
  // choose a random start point in vector
  int startindex=gRandom->Rndm()*(1.*MAXSIZE-1.*vecsize);
  volatile double distance;
  if(vecsize==1)
    {
      timer.Start();
      distance=testshape->DistFromOutside( &points_dO[3*startindex], &dirs_dO[3*startindex], 3, TGeoShape::Big(), 0);
      timer.Stop();
    }
  else
    {
      timer.Start();
      for(unsigned int index=startindex; index < startindex+vecsize; ++index)
	{
	  distance=testshape->DistFromOutside( &points_dO[3*index], &dirs_dO[3*index], 3, TGeoShape::Big(), 0);
	}
      timer.Stop();
    }
  Tacc+=timer.getDeltaSecs();
}



void ShapeBenchmarker::timeIt( )
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
