#include "TGeoBBox.h"
#include "TGeoMatrixPlain.h"
//#include "TGeoMatrix.h"
#include "PerformanceTesterTGeoMatrix.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
  // make a Box
  const double dx=10; // these are half-distances
  const double dy=20;
  const double dz=30;
  TGeoBBox *box = new TGeoBBox(dx, dy, dz);

  // setup a trafo
  TGeoMatrix r1(0,0,0,90,0,30);    
  // TGeoTranslation t1(-10,10,0);
  // TGeoRotation r1;r1.SetAngles(90,0,30);

  // test vector performance of cone
  TrafoBenchmarker test( box, &r1 );
  test.timeIt();// initDataContains();

  delete box;
}
