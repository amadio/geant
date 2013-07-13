#include "TGeoTube_v.h"
#include "VecPerformanceTesterT.h" // for the vectorized interface
#include "TGeoManager.h"
#include "TObjArray.h"
#include <iostream>
#include <cassert>

//  this tests SIMD vector implementations of TGeoTube
//  volume 18 from CMS.root = TGeoTube with inner and outer radius

int main(int argc, char *argv[])
{
  TGeoManager::Import("cms.root");
  TObjArray *vlist= gGeoManager->GetListOfVolumes();
  TGeoTube_v *atube = (TGeoTube_v *)( ((TGeoVolume*) vlist->UncheckedAt(18))->GetShape());

  atube->InspectShape();
  ShapeBenchmarker_v<TGeoTube_v> testtube( atube );
  testtube.timeIt();
}
