#include "TGeoCone.h"
#include "PerformanceTester.h"
#include <iostream>

int main(int argc, char *argv[])
{
  // make a cone
  const double dz=10; // these are half-distances
  const double Rmin1=0;
  const double Rmax1=30;
  const double Rmin2=0;
  const double Rmax2=50;
  TGeoCone *cone = new TGeoCone(dz, Rmin1, Rmax2, Rmin2, Rmax2);


  // test vector performance of cone
  ShapeBenchmarker test( cone );
  test.timeIt();// initDataContains();

  delete cone;
}
