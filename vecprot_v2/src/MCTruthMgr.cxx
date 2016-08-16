#include "MCTruthMgr.h"

//______________________________________________________________________________
void MCTruthMgr::OpenEvent(int evID) {
  
   MCEvent* current_event = new MCEvent();
   current_event->event_id = evID;
   events_map[evID+1] = current_event;
   
}

//______________________________________________________________________________
void MCTruthMgr::AddTrack(Geant::GeantTrack &gtrack) {
  
  // get the event from the map
  MCEvent* current_event = events_map.find(gtrack.fEvent+1);;
  
  if(CheckTrack(gtrack, current_event))
    {
      /*
	Printf("MCTruthMgr: Adding to event %i track ID %i mother %i energy %f pdg %i X %f Y %f Z %f time %f",
	gtrack.fEvent, gtrack.fParticle,
	gtrack.fMother, gtrack.fE, gtrack.fPDG, gtrack.fXpos, gtrack.fYpos, gtrack.fZpos, gtrack.fTime);
      */ 
      
      MCParticle* part = new MCParticle();
      part->pid = gtrack.fPDG;
      part->motherid = gtrack.fMother;
      part->fPx = gtrack.fXdir * gtrack.fP; 
      part->fPy = gtrack.fYdir * gtrack.fP;
      part->fPz = gtrack.fZdir * gtrack.fP;
      part->fXpos = gtrack.fXpos; 
      part->fYpos = gtrack.fYpos;
      part->fZpos = gtrack.fZpos;
      part->fTime = gtrack.fTime;
      part->fE = gtrack.fE;
      part->has_end = false;

      // add to the map GeantTrackID -> MCParticle*
      (current_event->particles)[gtrack.fParticle] = part;
    }
  /*
    else Printf("MCTruthMgr: NOT Adding to event %i track ID %i mother %i energy %f pdg %i X %f Y %f Z %f time %f",
    gtrack.fEvent, gtrack.fParticle,
    gtrack.fMother, gtrack.fE, gtrack.fPDG, gtrack.fXpos, gtrack.fYpos, gtrack.fZpos, gtrack.fTime);
  */
}

//______________________________________________________________________________
void MCTruthMgr::EndTrack(const Geant::GeantTrack_v &tracks, int itr) {
    
  // get the event to which the track belongs  
  if(events_map.contains(tracks.fEventV[itr]+1))
    {
      MCEvent* current_event = events_map.find(tracks.fEventV[itr]+1);

      // check of the track was stored
      // if not, do nothing
      if(current_event->particles.contains(tracks.fParticleV[itr]) )
	{
	  //  Printf("MCTruthMgr: Ending track ID %i in event %i", tracks.fParticleV[itr], tracks.fEventV[itr]);

	  MCParticle* current_particle = current_event->particles.find(tracks.fParticleV[itr]);

	  current_particle->fXend = tracks.fXposV[itr];
	  current_particle->fYend = tracks.fYposV[itr];
	  current_particle->fZend = tracks.fZposV[itr];
	  current_particle->fTend = tracks.fTimeV[itr];
	  current_particle->has_end = true;
	}
      else return;      
    }
  else return; 
}
