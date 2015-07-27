#include "TGeoManager.h"
#include "TGeoBBox_v.h"
#include "TRandom.h"

#include <iostream>
#include "tbb/tick_count.h"" // timing from Intel TBB 
#include <cassert>

struct TStopWatch 
{
  tbb::tick_count t1;
  tbb::tick_count t2;
  void Start(){ t1=tbb::tick_count::now(); }
  void Stop(){ t2=tbb::tick_count::now(); }
  void Reset(){ /* */ ;}
  void Print(){  std::cerr << (t2 - t1).seconds() << std::endl; }
  double getDeltaSecs() { return (t2-t1).seconds(); }
};

/*
#ifndef __INTEL_COMPILER
void * _mm_malloc(int s, int p)
{
  void *v;
  posix_memalign(&v,s,p);
  return v;
  // return malloc(s);
}

void _mm_free(void *p)
{
  free(p);
}
#endif
*/

#define NREP 1000

main(int argc, char *argv[])
{
  int npoints=100;
  if(argc>1) sscanf(argv[1],"%d",&npoints);
  printf("npoints = %d\n",npoints);

  const Double_t dx=10; // these are half-distances
  const Double_t dy=20;
  const Double_t dz=30;

  TGeoManager *testvec = new TGeoManager("Test","This is a naive test");
  TGeoMaterial *vacmat = new TGeoMaterial("vacuum",0,0,0);
  TGeoMedium *vacmed = new TGeoMedium("vacuum",0,vacmat);

  TGeoVolume *world = testvec->MakeBox("world",vacmed,100,100,100);
  testvec->SetTopVolume(world);

  TGeoVolume *tbox = testvec->MakeBox("tbox",vacmed,dx,dy,dz);
  tbox->SetLineColor(kRed);
  tbox->SetFillColor(kRed);
  tbox->SetVisibility(1);
  world->AddNode(tbox,1,0);
  
  testvec->CloseGeometry();

  Double_t origin[3]={0,0,0};
  TGeoBBox_v *box = new TGeoBBox_v(dx, dy, dz,origin);
  const Double_t r3two = pow(2,1./3.);

  npoints=1;
  for(int i = 0 ;i < 14; i++) 
    {
      Double_t *points = new Double_t[3*npoints];
      Double_t *dir = new Double_t[3*npoints];

      StructOfCoord p,d;
      p.alloc(npoints);
      d.alloc(npoints);

      TStopWatch tt;

      points[3*i] = -10*dx;
      points[3*i+1] = -10*dy;
      points[3*i+2] = -10*dz;

      for(int i=1; i<npoints; ++i) {
	/*
	points[3*i  ]=0;
	points[3*i+1]=0;
	points[3*i+2]=0;
	*/

	// they are all completely in
	points[3*i  ]=(1-2.*gRandom->Rndm())*dx;
	points[3*i+1]=(1-2.*gRandom->Rndm())*dy;
	points[3*i+2]=(1-2.*gRandom->Rndm())*dz;

	dir[3*i  ]=(1-2.*gRandom->Rndm());
	dir[3*i+1]=(1-2.*gRandom->Rndm());
	dir[3*i+2]=(1-2.*gRandom->Rndm());
      }
      p.fill(points);
      d.fill(dir);

      // TODO: introduce some special boundary cases for tests

      double *distance_v = new double[npoints];

      double DeltaT=0., DeltaT_v=0., DeltaT_l=0., DeltaT_b=0.;
      for ( unsigned int repetitions = 0; repetitions < NREP; repetitions ++ ) 
	{
	  // assert correctness of result (simple checksum check)
	  {
	    double checksum=0., checksum_l=0., checksum_box=0., checksum_v=0.;
	    #pragma novector 
	    for(int i=0; i<npoints; ++i) {
	      distance_v[i]=TGeoBBox_v::DistFromInsideS(&points[3*i],&dir[3*i], dx,dy,dz, origin, TGeoShape::Big());
 	      checksum+=distance_v[i];
	    }
	    
	    TGeoBBox_v::DistFromInsideS_l(points, dir, dx, dy, dz, origin, TGeoShape::Big(), distance_v, npoints);
	    for(int i=0; i<npoints; ++i) {
	      checksum_l+=distance_v[i];
	    }

	    TGeoBBox_v::DistFromInsideS_v(p, d, dx, dy, dz, origin, TGeoShape::Big(), distance_v, npoints);
	    for(int i=0; i<npoints; ++i) {
	      checksum_v+=distance_v[i];
	    }

	    // crosscheck with result from box instance ( non - statically )
	    for(int i=0; i<npoints; ++i) {
	      distance_v[i]=box->DistFromInside(&points[3*i],&dir[3*i], 3, TGeoShape::Big(), 0);
 	      checksum_box+=distance_v[i];
	    }

	    if(checksum_v != checksum || checksum_box != checksum );
	    {
	      //  std::cerr << "#" << checksum_box << " " << checksum_v <<  " " << checksum << " " <<  checksum_l << std::endl;
	    }	  
	  }

	  
	  tt.Start();
	  TGeoBBox_v::DistFromInsideS_v(p, d, dx, dy, dz, origin, TGeoShape::Big(), distance_v, npoints);
	  tt.Stop();
	  DeltaT_v+= tt.getDeltaSecs(); //      tt.Print();
	  tt.Reset();

	  tt.Start();
	  TGeoBBox_v::DistFromInsideS_l(points, dir, dx, dy, dz, origin, TGeoShape::Big(), distance_v, npoints);
	  tt.Stop();
	  DeltaT_l+= tt.getDeltaSecs(); //      tt.Print();
	  tt.Reset();


	  tt.Start();
	  #pragma novector
	  for(int i=0; i<npoints; ++i) {
	    distance_v[i] = TGeoBBox_v::DistFromInsideS(&points[3*i],&dir[3*i], dx,dy,dz, origin, TGeoShape::Big());
	  }
	  tt.Stop();
	  DeltaT+= tt.getDeltaSecs();
	  tt.Reset();

	  // crosscheck with result from box instance ( non - statically )
	  tt.Start();
	  for(int i=0; i<npoints; ++i) {
	    distance_v[i]=box->DistFromInside(&points[3*i],&dir[3*i], 3, TGeoShape::Big(), 0);
	  }
	  tt.Stop();
	  DeltaT_b+= tt.getDeltaSecs();
	  tt.Reset();
	}

      std::cerr << "#P " << npoints << " " << DeltaT/NREP << " " << DeltaT_l/NREP << " " << DeltaT_v/NREP << " " << DeltaT/DeltaT_l << " " << DeltaT/DeltaT_v << " " << DeltaT_b/DeltaT_v << std::endl;
      
      delete[] distance_v;
      delete[] dir;
      delete[] points;

      p.dealloc();
      d.dealloc();

      npoints*=2;
    }
  return 0;
}
