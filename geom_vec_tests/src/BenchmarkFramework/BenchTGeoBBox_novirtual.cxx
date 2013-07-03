#include "TGeoCone.h"
#include "PerformanceTesterT.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
  // make a Box
  const Double_t dx=10; // these are half-distances
  const Double_t dy=20;
  const Double_t dz=30;
  TGeoBBox *box = new TGeoBBox(dx, dy, dz);
  box->Draw();

  // test vector performance of cone
  ShapeBenchmarker<TGeoBBox> test( box );
  test.timeIt();// initDataContains();

  delete box;
}
