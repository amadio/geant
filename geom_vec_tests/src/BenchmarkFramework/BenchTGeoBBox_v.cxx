#include "TGeoBBox_v.h"
#include "VecPerformanceTesterT.h" // for the vectorized interface
#include <iostream>

int main(int argc, char *argv[])
{
  // make a Box
  const double dx=10; // these are half-distances
  const double dy=20;
  const double dz=30;
  TGeoBBox_v *box = new TGeoBBox_v(dx, dy, dz);

  // test vector performance of cone
  //  ShapeBenchmarker_v<TGeoBBox_v> test( box );
  ShapeBenchmarker_v<TGeoBBox_v> test( box );
  test.timeIt();// initDataContains();

  delete box;
}
