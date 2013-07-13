#include "TGeoCone_v.h"
#include "VecPerformanceTesterT.h" // for the vectorized interface
#include "TGeoManager.h"
#include "TObjArray.h"
#include <iostream>
#include <cassert>

//  this tests SIMD vector implementations of TGeoCone
//  volume 374 from CMS.root = TGeoCone 

int main(int argc, char *argv[])
{
  TGeoManager::Import("cms.root");
  TObjArray *vlist= gGeoManager->GetListOfVolumes();
  TGeoCone_v *acone = (TGeoCone_v *)( ((TGeoVolume*) vlist->UncheckedAt(374))->GetShape());


  acone->InspectShape();
  ShapeBenchmarker_v<TGeoCone_v> testcone( acone );
  testcone.timeIt();
}
