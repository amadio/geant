#include "TGeoBBox_v.h"
#include "MaskedVecPerformanceTesterT.h" // for the vectorized interface
#include <iostream>

int main(int argc, char *argv[])
{
  // make a Box
  const double dx=10; // these are half-distances
  const double dy=20;
  const double dz=30;
  TGeoBBox_v *box = new TGeoBBox_v(dx, dy, dz);

  // test vector performance of box with masked particles
  ShapeBenchmarker_masked_v<TGeoBBox_v> test( box );
  test.timeIt();

  delete box;
}
