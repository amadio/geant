#include "MCTruthMgr.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
void MCTruthMgr::OpenEvent(int evID) {
  
   MCEvent* current_event = new MCEvent();
   current_event->event_id = evID;
   events_map[evID] = current_event;
   
}

//______________________________________________________________________________
void MCTruthMgr::AddTrack(geant::GeantTrack &gtrack) {

  // get the event from the map  
  MCEvent* current_event = events_map.find(gtrack.Event());
  
  if(CheckTrack(gtrack, current_event))
    {
      /*
	Printf("MCTruthMgr: Adding to event %i track ID %i mother %i energy %f pdg %i X %f Y %f Z %f time %f",
	gtrack.fEvent, gtrack.fParticle,
	gtrack.fMother, gtrack.fE, gtrack.fPDG, gtrack.fXpos, gtrack.fYpos, gtrack.fZpos, gtrack.fTime);
      */ 
      
      MCParticle* part = new MCParticle();
      part->pid = gtrack.GVcode();
      part->motherid = gtrack.Mother();
      part->fPx = gtrack.Px(); 
      part->fPy = gtrack.Py();
      part->fPz = gtrack.Pz();
      part->fXpos = gtrack.X(); 
      part->fYpos = gtrack.Y();
      part->fZpos = gtrack.Z();
      part->fTime = gtrack.Time();
      part->fE = gtrack.E();
      part->has_end = false;

      // add to the map GeantTrackID -> MCParticle*
      (current_event->particles)[gtrack.Particle()] = part;
    }
  /*
    else Printf("MCTruthMgr: NOT Adding to event %i track ID %i mother %i energy %f pdg %i X %f Y %f Z %f time %f",
    gtrack.fEvent, gtrack.fParticle,
    gtrack.fMother, gtrack.fE, gtrack.fPDG, gtrack.fXpos, gtrack.fYpos, gtrack.fZpos, gtrack.fTime);
  */
}

//______________________________________________________________________________
void MCTruthMgr::EndTrack(GeantTrack *track) {

  MCEvent* current_event = 0;
  
  // get the event to which the track belongs  
  if(events_map.find(track->Event(), current_event))
    {
      // check of the track was stored
      // if not, do nothing

      MCParticle* current_particle = 0;
      if(current_event->particles.find(track->Particle(), current_particle))
	{
	  current_particle->fXend = track->X();
	  current_particle->fYend = track->Y();
	  current_particle->fZend = track->Z();
	  current_particle->fTend = track->Time();
	  current_particle->has_end = true;
	}
      else return;      
    }
  else return; 
}

} // GEANT_IMPL_NAMESPACE
} // Geant
