#include "TGeoPcon_v.h"
#include "VecPerformanceTesterT.h" // for the vectorized interface
#include <iostream>

//  this takes a couple of shapes from the CMS detector description 
//  and performs the standard benchmark
//  we know here ( study of volume list ), that
//  volume 2 = TGeoPcone


int main(int argc, char *argv[])
{
  TGeoManager::Import("cms.root");
  TObjArray *vlist= gGeoManager->GetListOfVolumes();
  TGeoPcon *apcon = (TGeoPcon *)( ((TGeoVolume*) vlist->UncheckedAt(2))->GetShape());

  // test vector performance for Pcon
  apcon->InspectShape();
  ShapeBenchmarker_v<TGeoPcon_v> testpcon( apcon );
  testpcon.timeIt();
  //  testpcon.printTimings( ss.str().c_str() );
}
