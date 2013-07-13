#include "TGeoCone_v.h"
#include "VecPerformanceTesterT.h" // for the vectorized interface
#include "TGeoManager.h"
#include "TObjArray.h"
#include <iostream>
#include <cassert>

//  this tests SIMD vector implementations of TGeoCone
//  volume 258 from CMS.root = TGeoConeSeg 

int main(int argc, char *argv[])
{
  TGeoManager::Import("cms.root");
  TObjArray *vlist= gGeoManager->GetListOfVolumes();
  TGeoConeSeg_v *acone = (TGeoConeSeg_v *)( ((TGeoVolume*) vlist->UncheckedAt(258))->GetShape());


  acone->InspectShape();
  ShapeBenchmarker_v<TGeoConeSeg_v> testcone( acone );
  testcone.timeIt();
}
