#include <iostream>
#include "GUXSectionKleinNishina.h"
#include "GUAliasSampler.h"
#include "GURandom.h"
#include "GUTrack.h"
#include "GUComptonKleinNishina.h"
#include "GUTrackHandler.h"

#include "base/Stopwatch.h"

using namespace vecgeom;

int main(int argc, char* argv[]) {

  //argument
  int ntrack   = 1000;
  int testType = 1; // 0 test run

  if(argc >= 2) ntrack = atoi(argv[1]);
  if(argc >= 3) testType = atoi(argv[2]);

  if(testType < 0 || testType > 2) {
    std::cout << "Usage: GUTrackingTest [ntrack=1000] [testType=1] " 
	      << std::endl;
    return 0;
  }
  if( testType == 0 ) ntrack = 10;

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

  Stopwatch timer; 
  timer.Start();

  model->Interact(track_in,targetElements,&track_out);

  timer.Stop();
  double vtime =  timer.Elapsed();

  if(testType==0) {
    for(int j = 0; j < ntrack ; ++j) {
      std::cout << "vector primary secondary E " <<  (track_in.E)[j] << " " 
		<< (track_out.E)[j] << std::endl;
    }
    std::cout << " " << std::endl;
  }
  //scalar
  GUTrack* track_aos = handler_in->GetAoSTracks();
  GUTrack* track_aos_out = (GUTrack*) malloc(ntrack*sizeof(GUTrack));

  timer.Start();
  for(int i = 0 ; i < ntrack ; ++i) {
    model->Interact(track_aos[i],targetElements[i],&track_aos_out[i]);
  }
  timer.Stop();
  double stime =  timer.Elapsed();

  if(testType==0) {
    for(int i = 0 ; i < ntrack ; ++i) {
      std::cout << "scalar primary secondary E " <<  track_aos[i].E << " " 
		<< track_aos_out[i].E << std::endl;
    }
  }

  //original method
  
  timer.Start();
  for(int i = 0 ; i < ntrack ; ++i) {
    model->InteractG4(track_aos[i],targetElements[i],&track_aos_out[i]);
  }
  timer.Stop();
  double g4time =  timer.Elapsed();

  std::cout << " vtime stime g4time " << vtime << " " << stime << " " 
	    << g4time << std::endl;

  //  delete handler_in;
  //  delete handler_out;
  //  delete model;

  free(track_aos_out);

  return 0;
}
