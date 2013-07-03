#include "TGeoCone.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "TGeoPcon.h"
#include "TGeoArb8.h"
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "PerformanceTesterT.h"
#include <iostream>
#include <sstream>

//  this takes a couple of shapes from the CMS detector description 
//  and performs the standard benchmark
//  we know here ( study of volume list ), that
//  volume 2 = TGeoPcone
//  volume 261 = TGeoConeSeg
//  volume 374 = TGeoCone
//  volume 281 = TGeoTubeSeg
//  volume 285 = TGeoTrap
//  volume 40 = TGeoBBox


int main(int argc, char *argv[])
{
  // get number of this run
  int run=0;
  if(argc>1) run=atoi(argv[1]);
  // output filenames
  std::stringstream ss;

  TGeoManager::Import("cms.root");
  TObjArray *vlist= gGeoManager->GetListOfVolumes();
  TGeoPcon *apcon = (TGeoPcon *)( ((TGeoVolume*) vlist->UncheckedAt(2))->GetShape());
  TGeoCone *acone = (TGeoCone *) ( ((TGeoVolume*) vlist->UncheckedAt(374))->GetShape());
  TGeoTubeSeg *atubeseg= (TGeoTubeSeg *) ( ((TGeoVolume*) vlist->UncheckedAt(281))->GetShape());
  TGeoConeSeg *aconeseg= (TGeoConeSeg *) ( ((TGeoVolume*) vlist->UncheckedAt(261))->GetShape());
  TGeoTrap *atrap= (TGeoTrap *)( ((TGeoVolume*) vlist->UncheckedAt(285))->GetShape());
  TGeoBBox *abox = (TGeoBBox *) ( ((TGeoVolume*) vlist->UncheckedAt(40))->GetShape());
  // test vector performance of cone

  apcon->InspectShape();
  ShapeBenchmarker<TGeoPcon> testpcon( apcon );
  testpcon.timeIt();
  ss << "CMStimingPcon_run"<<run<<".dat";
  testpcon.printTimings( ss.str().c_str() );

  acone->InspectShape();
  ShapeBenchmarker<TGeoCone> testcone( acone );
  testcone.timeIt();
  ss.str(""); ss << "CMStimingCone_run"<<run<<".dat";
  testcone.printTimings( ss.str().c_str() );

  atubeseg->InspectShape();
  ShapeBenchmarker<TGeoTubeSeg> testtubeseg( atubeseg );
  testtubeseg.timeIt();
  ss.str("");ss << "CMStimingTubeSeg_run"<<run<<".dat";
  testtubeseg.printTimings( ss.str().c_str() );

  aconeseg->InspectShape();
  ShapeBenchmarker<TGeoConeSeg> testconeseg( aconeseg );
  testconeseg.timeIt();
  ss.str(""); ss<< "CMStimingConeSeg_run"<<run<<".dat";
  testconeseg.printTimings( ss.str().c_str() );

  atrap->InspectShape();
  ShapeBenchmarker<TGeoTrap> testtrap( atrap );
  testtrap.timeIt();
  ss.str(""); ss << "CMStimingTrap_run"<<run<<".dat";
  testtrap.printTimings( ss.str().c_str() );

  abox->InspectShape();
  ShapeBenchmarker<TGeoBBox> testbox( abox );
  testbox.timeIt();
  ss.str("");ss << "CMStimingBox_run"<<run<<".dat";
  testbox.printTimings( ss.str().c_str() );
}
