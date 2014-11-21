#include <iostream>
#include "GUXSectionKleinNishina.h"
#include "GUAliasSampler.h"
#include "GURandom.h"
#include "GUTrack.h"
#include "GUComptonKleinNishina.h"
#include "GUTrackHandler.h"

using namespace vecgeom;

int main(int argc, char* argv[]) {

  int ntrack = 10;
  std::cout << "GUTrackingTest Processing Ntracks  = " << ntrack << std::endl;

  GUTrackHandler *handler_in = new GUTrackHandler();
  handler_in->GenerateRandomTracks(ntrack);
  GUTrack_v track_in = handler_in->GetSoATracks();

  GUTrackHandler *handler_out = new GUTrackHandler(ntrack);
  GUTrack_v track_out = handler_out->GetSoATracks();

  //put GUComptonProcess here which will call
  GUComptonKleinNishina *model = new GUComptonKleinNishina();

  int *targetElements = new int(ntrack);
  for(int i = 0 ; i < ntrack ; ++i) {
    targetElements[i] = i ;
  }

  model->Interact(track_in,targetElements,&track_out);

  for(int j = 0; j < ntrack ; ++j) {
    std::cout << "secondary E " <<  (track_out.E)[j] << std::endl;
  }

  delete handler_in;
  delete handler_out;
  delete model;

  return 0;
}
