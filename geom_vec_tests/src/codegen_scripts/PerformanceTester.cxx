#include "PerformanceTester.h"
#include "TRandom.h"
//#include "TMath.h"
#include <iostream>

class Util{
  public:
  static void samplePoint( double *point, double dx, double dy, double dz, double scale )
  {
    point[0]=scale*(1-2.*gRandom->Rndm())*dx;
    point[1]=scale*(1-2.*gRandom->Rndm())*dy;
    point[2]=scale*(1-2.*gRandom->Rndm())*dz;
  }
  
  void sampleDir( double * dir)
  {
    
  }
  
};


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
	}	
      insidecounter++;
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
	  isin=testshape->Contains( &points_C[3*startindex] );
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
	  safety=testshape->Safety( &points_s[3*startindex], kTRUE );
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


  // to avoid measuring the same function over and over again we interleave calls to different functions and different data
  for(unsigned int rep=0; rep< NREPS; rep++)
    {
      for(unsigned int vectype =0 ; vectype < N; ++vectype )
	{
	  // Contains
	  timeContains( Tc[vectype], vecsizes[vectype] );

	  // Safety
	  timeSafety ( Ts[vectype], vecsizes[vectype] );
	}
    }

  
  // print result
  for(unsigned int vectype =0 ; vectype < N; ++vectype )
    {
      std::cerr << vecsizes[vectype] << " " << Tc[vectype]/NREPS << " " << Tc[vectype]/NREPS/vecsizes[vectype] << " " <<  Ts[vectype]/NREPS << " " << Ts[vectype]/NREPS/vecsizes[vectype] << std::endl;
    }  
} 
