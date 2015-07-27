#include "TGeoManager.h"
#include "TGeoBBox_v.h"
#include "TRandom.h"

#include <iostream>
#include "tbb/tick_count.h" // timing from Intel TBB
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

  const double dx = 10; // these are half-distances
  const double dy = 20;
  const double dz = 30;

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

  double origin[3] = {0, 0, 0};
  TGeoBBox_v *box = new TGeoBBox_v(dx, dy, dz, origin);
  const double r3two = pow(2, 1. / 3.);

  npoints = 10;
  for (int i = 0; i < 14; i++) {
    double *points = new double[3 * npoints];
    double *dir = new double[3 * npoints];
    double *norm = new double[3 * npoints];
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

      dir[3 * i] = (1 - 2. * gRandom->Rndm()) * dx; // correct?
      dir[3 * i + 1] = (1 - 2. * gRandom->Rndm()) * dy;
      dir[3 * i + 2] = (1 - 2. * gRandom->Rndm()) * dz;
    }

    double DeltaT = 0., DeltaT_v = 0., DeltaT_l = 0.;
    for (unsigned int repetitions = 0; repetitions < NREP; repetitions++) {
      // assert correctness of result (simple checksum check)
      {
        double checksum = 0., checksum_v = 0.;
        for (int i = 0; i < npoints; ++i) {
          box->ComputeNormal(&points[3 * i], &dir[3 * i], &norm[3 * i]);
          //	      std::cerr << "Ax " << points[3*i+0] << " y " << points[3*i+1] << " z " << points[3*i+2] << " safet
          //" << safety_v[i] << std::endl;
        }

        for (int i = 0; i < npoints * 3; i++)
          checksum += norm[i];

        box->ComputeNormal_l(points, dir, norm, npoints);
        for (int i = 0; i < npoints * 3; ++i) {
          //	      std::cerr << "Bx " << points[3*i+0] << " y " << points[3*i+1] << " z " << points[3*i+2] << " safet
          //" << safety_v[i] << std::endl;
          checksum_v += norm[i];
        }
        assert(checksum_v == checksum);
      }

      // measure timings here separately
      tt.Start();
      for (int i = 0; i < npoints; ++i) {
        box->ComputeNormal(&points[3 * i], &dir[3 * i], &norm[3 * i]);
      }
      tt.Stop();
      DeltaT += tt.getDeltaSecs();
      tt.Reset();

      tt.Start();
      box->ComputeNormal_v(points, dir, norm, npoints);
      tt.Stop();
      DeltaT_v += tt.getDeltaSecs(); //      tt.Print();
      tt.Reset();

      tt.Start();
      box->ComputeNormal_l(points, dir, norm, npoints);
      tt.Stop();
      DeltaT_l += tt.getDeltaSecs(); //      tt.Print();
      tt.Reset();
    }

    std::cerr << "P# " << npoints << " " << DeltaT / NREP << " " << DeltaT_l / NREP << " " << DeltaT_v / NREP << " "
              << DeltaT / DeltaT_l << " " << DeltaT / DeltaT_v << std::endl;

    delete[] dir;
    delete[] norm;
    delete[] points;
    npoints *= 2;
  }
  return 0;
}
