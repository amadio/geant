#include "TGeoManager.h"
#include "TGeoBBox_v.h"
#include "TRandom.h"

#include <iostream>
#include "tbb/tick_count.h" " // timing from Intel TBB 
#include <cassert>

struct TStopWatch {
  tbb::tick_count t1;
  tbb::tick_count t2;
  void Start() { t1 = tbb::tick_count::now(); }
  void Stop() { t2 = tbb::tick_count::now(); }
  void Reset() { /* */
    ;
  }
  void Print() { std::cerr << (t2 - t1).seconds() << std::endl; }
  double getDeltaSecs() { return (t2 - t1).seconds(); }
};

#define NREP 1000

main(int argc, char *argv[]) {
  int npoints = 100;
  if (argc > 1)
    sscanf(argv[1], "%d", &npoints);
  printf("npoints = %d\n", npoints);

  const Double_t dx = 10; // these are half-distances
  const Double_t dy = 20;
  const Double_t dz = 30;

  TGeoManager *testvec = new TGeoManager("Test", "This is a naive test");
  TGeoMaterial *vacmat = new TGeoMaterial("vacuum", 0, 0, 0);
  TGeoMedium *vacmed = new TGeoMedium("vacuum", 0, vacmat);

  TGeoVolume *world = testvec->MakeBox("world", vacmed, 100, 100, 100);
  testvec->SetTopVolume(world);

  TGeoVolume *tbox = testvec->MakeBox("tbox", vacmed, dx, dy, dz);
  tbox->SetLineColor(kRed);
  tbox->SetFillColor(kRed);
  tbox->SetVisibility(1);
  world->AddNode(tbox, 1, 0);

  testvec->CloseGeometry();

  Double_t origin[3] = {0, 0, 0};
  TGeoBBox_v *box = new TGeoBBox_v(dx, dy, dz, origin);
  const Double_t r3two = pow(2, 1. / 3.);

  npoints = 1;
  for (int i = 0; i < 14; i++) {
    Double_t *points = new Double_t[3 * npoints];

    StructOfCoord p;
    p.alloc(npoints);

    TStopWatch tt;

    for (int i = 0; i < npoints; ++i) {
      /*
      points[3*i  ]=0;
      points[3*i+1]=0;
      points[3*i+2]=0;
      */

      points[3 * i] = (1 - 2. * gRandom->Rndm()) * dx;
      points[3 * i + 1] = (1 - 2. * gRandom->Rndm()) * dy;
      points[3 * i + 2] = (1 - 2. * gRandom->Rndm()) * dz;
    }
    p.fill(points);

    double *safety_v = new double[npoints];
    double DeltaT = 0., DeltaT_v = 0., DeltaT_l = 0.;
    for (unsigned int repetitions = 0; repetitions < NREP; repetitions++) {
      // assert correctness of result (simple checksum check)
      {
        double checksum = 0., checksum_v = 0.;
        for (int i = 0; i < npoints; ++i) {
          safety_v[i] = box->Safety(&points[3 * i], true);
          //	      std::cerr << "Ax " << points[3*i+0] << " y " << points[3*i+1] << " z " << points[3*i+2] << " safet
          //" << safety_v[i] << std::endl;
          checksum += safety_v[i];
        }

        box->Safety_v(p, safety_v, npoints, true);
        for (int i = 0; i < npoints; ++i) {
          //	      std::cerr << "Bx " << points[3*i+0] << " y " << points[3*i+1] << " z " << points[3*i+2] << " safet
          //" << safety_v[i] << std::endl;
          checksum_v += safety_v[i];
        }
        assert(checksum_v == checksum);
      }

      // measure timings here separately
      tt.Start();
      for (int i = 0; i < npoints; ++i) {
        safety_v[i] = box->Safety(&points[3 * i], true);
      }
      tt.Stop();
      DeltaT += tt.getDeltaSecs();
      tt.Reset();

      tt.Start();
      box->Safety_v(p, safety_v, npoints, true);
      tt.Stop();
      DeltaT_v += tt.getDeltaSecs(); //      tt.Print();
      tt.Reset();

      tt.Start();
      box->Safety_l(points, safety_v, npoints, true);
      tt.Stop();
      DeltaT_l += tt.getDeltaSecs(); //      tt.Print();
      tt.Reset();
    }

    std::cerr << "#P " << npoints << " " << DeltaT / NREP << " " << DeltaT_l / NREP << " " << DeltaT_v / NREP << " "
              << DeltaT / DeltaT_l << " " << DeltaT / DeltaT_v << std::endl;

    delete[] safety_v;
    delete[] points;
    p.dealloc();
    npoints *= 2;
  }
  return 0;
}
