#include "TGeoBBox_v.h"
#include "VecPerformanceTesterT.h" // for the vectorized interface
#include <iostream>

int main(int argc, char *argv[])
{
  // make a Box
  const Double_t dx=10; // these are half-distances
  const Double_t dy=20;
  const Double_t dz=30;
  TGeoBBox_v *box = new TGeoBBox_v(dx, dy, dz);

  // test vector performance of cone
  //  ShapeBenchmarker_v<TGeoBBox_v> test( box );
  ShapeBenchmarker_v<TGeoBBox_v> test( box );
  test.timeIt();// initDataContains();

  delete box;
}
