#include "TGeoPcon_v.h"
#include "VecPerformanceTesterT.h" // for the vectorized interface
#include "TGeoManager.h"
#include "TObjArray.h"
#include "TGeoPcon.h"
#include <iostream>

//  this takes a couple of shapes from the CMS detector description 
//  and performs the standard benchmark
//  we know here ( study of volume list ), that
//  volume 312 = TGeoPcone (Segment)


int main(int argc, char *argv[])
{
  TGeoManager::Import("cms.root");
  TObjArray *vlist= gGeoManager->GetListOfVolumes();
  TGeoPcon_v *apcon = (TGeoPcon_v *)( ((TGeoVolume*) vlist->UncheckedAt(312))->GetShape());

  // test vector performance for Pcon
  apcon->InspectShape();
  ShapeBenchmarker_v<TGeoPcon_v> testpcon( apcon );
  testpcon.timeIt();
}
