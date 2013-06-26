#include "TGeoCone.h"
#include "PerformanceTester.h"
#include <iostream>

int main(int argc, char *argv[])
{
  // make a Box
  const Double_t dx=10; // these are half-distances
  const Double_t dy=20;
  const Double_t dz=30;
  TGeoBBox *box = new TGeoBBox(dx, dy, dz);

  // test vector performance of cone
  ShapeBenchmarker test( box );
  test.timeIt();// initDataContains();

  delete box;
}
