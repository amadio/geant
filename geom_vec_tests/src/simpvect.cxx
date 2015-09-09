#include "TGeoManager.h"
#include "TGeoBBox_v.h"
#include "TRandom.h"
#include "TStopwatch.h"

int
#ifdef __CINT__
simpvect(int npoints=10000)
{
#else
main(int argc, char *argv[])
{
  int npoints = 10000;
  if (argc > 1)
    sscanf(argv[1], "%d", &npoints);
  printf("npoints = %d\n", npoints);
#endif

  const double dx = 10;
  const double dy = 20;
  const double dz = 30;

  TGeoManager *testvec = new TGeoManager("Test", "This is a naive test");
  TGeoMaterial *vacmat = new TGeoMaterial("vacuum", 0, 0, 0);
  TGeoMedium *vacmed = new TGeoMedium("vacuum", 0, vacmat);

  TGeoVolume *world = testvec->MakeBox("world", vacmed, 100, 100, 100);
  testvec->SetTopVolume(world);

  TGeoVolume *tbox = testvec->MakeBox("tbox", vacmed, 10, 20, 30);
  tbox->SetLineColor(kRed);
  tbox->SetFillColor(kRed);
  tbox->SetVisibility(1);
  world->AddNode(tbox, 1, 0);

  testvec->CloseGeometry();

  double origin[3] = {0, 0, 0};

  TGeoBBox_v *box = new TGeoBBox_v(dx, dy, dz, origin);

  double *points = new double[3 * npoints];

  const double r3two = pow(2, 1. / 3.);

  TStopwatch tt;
  tt.Start();
  for (int i = 0; i < npoints; ++i) {
    points[3 * i] = r3two * (1 - 2 * gRandom->Rndm()) * dx;
    points[3 * i + 1] = r3two * (1 - 2 * gRandom->Rndm()) * dy;
    points[3 * i + 2] = r3two * (1 - 2 * gRandom->Rndm()) * dz;
  }
  tt.Stop();
  tt.Print();
  tt.Reset();

  int ipin = 0;
  bool *isin = new bool[npoints];

  tt.Start();
  for (int i = 0; i < npoints; ++i) {
    isin[i] = box->Contains(&points[3 * i]);
  }
  tt.Stop();
  tt.Print();
  tt.Reset();

  tt.Start();
  box->Contains_v(points, isin, npoints);
  tt.Print();
  tt.Reset();

  return 0;
}
