#include "TGeoManager.h"
#include "TGeoBBox_v.h"
#include "TRandom.h"

#include <iostream>
#include "tbb/tick_count.h" // timing from Intel TBB 
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
  for(int s = 0 ;s < 14; s++) 
    {
        Double_t *points = new Double_t[3*npoints];
        Double_t *dir = new Double_t[3*npoints];

	StructOfCoord p,d;
	p.alloc(npoints);
	d.alloc(npoints);
        
        Double_t *step=(double*)_mm_malloc(npoints*sizeof(double),32);//.;new Double_t[npoints];
        TStopWatch tt;
	points[3*0  ]=10*dx;
	points[3*0+1]=0;
	points[3*0+2]=0;
	dir[3*0  ]=-1.;
	dir[3*0+1]=0.;
	dir[3*0+2]=0.;
	step[0]=TGeoShape::Big();

        for(int i=1; i<npoints; ++i) {
	  //            points[3*i  ]=1.4*(1-2.*gRandom->Rndm())*dx;
	  //            points[3*i+1]=1.4*(1-2.*gRandom->Rndm())*dy;
	  //            points[3*i+2]=1.4*(1-2.*gRandom->Rndm())*dz;
	  
	  // most of test particles should be outside to be fair
	  int signx = (gRandom->Rndm() < 0.5)? -1 : 1;
	  int signy = (gRandom->Rndm() < 0.5)? -1 : 1;
	  int signz = (gRandom->Rndm() < 0.5)? -1 : 1;

	  points[3*i  ]=signx*(dx + 0.1*(gRandom->Rndm()-0.1));
	  points[3*i+1]=signy*(dy + 0.1*(gRandom->Rndm()-0.1));
	  points[3*i+2]=signz*(dz + 0.1*(gRandom->Rndm()-0.1));
	  
	  dir[3*i  ]=-1.;//(1-2.*gRandom->Rndm());
	  dir[3*i+1]=0.;///(1-2.*gRandom->Rndm());
	  dir[3*i+2]=0.;//(1-2.*gRandom->Rndm());
          
	  step[i]=0.001; // TGeoShape::Big();
	}
	p.fill(points);
	d.fill(dir);
      
	Double_t *distance = new Double_t[npoints];
	Double_t *distance_v = (double*)_mm_malloc(npoints*sizeof(double),32);//.;new Double_t[npoints];new Double_t[npoints];

	// assert correctness of result (simple checksum check)
	double checksum=0., checksum_v=0.;
	// perform a check
	for(int i=0; i<npoints; ++i)
	  distance[i]=TGeoBBox_v::DistFromOutsideS(&points[3*i], &dir[3*i],dx,dy,dz,origin, step[i]);
	
	TGeoBBox_v::DistFromOutsideS_v(p, d,dx,dy,dz,origin,step, distance_v, npoints);
	
	for(int i=0; i<npoints; i++)
	  {	
	    // std::cerr << i << " " << distance[i] << " " << distance_v[i] << std::endl;
	    assert(distance[i] == distance_v[i]);
	  }

	double DeltaT=0., DeltaT_l=0., DeltaT_v=0.;
	for ( unsigned int repetitions = 0; repetitions < NREP; repetitions ++ ) 
	  {
	    tt.Start();
	    for(int i=0; i<npoints; i++)
	      distance[i]=TGeoBBox_v::DistFromOutsideS( &points[3*i], &dir[3*i], dx, dy, dz, origin, step[i] );
	    tt.Stop();
	    DeltaT+= tt.getDeltaSecs(); //      tt.Print();
	    tt.Reset();
	    
	    tt.Start();
	    TGeoBBox_v::DistFromOutsideS_v(p, d, dx, dy, dz, origin, step, distance_v, npoints);
	    tt.Stop();
	    DeltaT_v+= tt.getDeltaSecs(); //      tt.Print();
	    tt.Reset();

	    tt.Start();
	    //	    TGeoBBox_v::DistFromOutside_l(p, d, dx, dy, dz, origin, step, distance_v, npoints);
	    TGeoBBox_v::DistFromOutsideS_l(points, dir, dx, dy, dz, origin, step, distance_v, npoints);
	    tt.Stop();
	    DeltaT_l+= tt.getDeltaSecs(); //      tt.Print();
	    tt.Reset();
	  }
	std::cerr <<  "#P " << npoints << " " << DeltaT/NREP << " " << DeltaT_l/NREP << " " << DeltaT_v/NREP << " "  << DeltaT/DeltaT_l << " "  << DeltaT/DeltaT_v  << std::endl;

	delete[] dir;
	_mm_free(step);//delete[] step;
	delete[] points;
	delete[] distance;
	_mm_free(distance_v);
	p.dealloc();
	d.dealloc();
	
	npoints*=2;
    }
  return 0;
}
