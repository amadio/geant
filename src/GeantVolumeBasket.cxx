#include "GeantVolumeBasket.h"
#include "globals.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "PhysicsProcess.h"

#include "TThread.h"
#include "TArrayI.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TGeoBranchArray.h"

ClassImp(GeantVolumeBasket)

const Double_t gTolerance = TGeoShape::Tolerance();

//______________________________________________________________________________
GeantVolumeBasket::GeantVolumeBasket(TGeoVolume *vol)
                  :TObject(),
                   fVolume(vol),
                   fNtracks(0),
                   fNpending(0),
                   fFirstFree(0),
                   fMaxTracks(10),
                   fMaxPending(10),
                   fIndex(0),
                   fPending(0)
{
// Constructor
}                   

//______________________________________________________________________________
GeantVolumeBasket::~GeantVolumeBasket()
{
// Clean up
   delete [] fIndex;
}   

//______________________________________________________________________________
TGeoBranchArray *GeantVolumeBasket::GetBranchArray(Int_t itrack) const
{
// Return path for the track.
   return gPropagator->fTracks[fIndex[itrack]]->path;
}

//______________________________________________________________________________
Int_t GeantVolumeBasket::GetNchunks(Int_t nthreads) const
{
// Get optimum number of work chunks for the basket for nthreads. The basket
// should be flushed.
   Int_t nchunk = fNtracks/nthreads;
   Int_t nworkers = nthreads;
   if (!nchunk) nworkers = fNtracks;
   return nworkers;   
}

//______________________________________________________________________________
GeantTrack *GeantVolumeBasket::GetTrack(Int_t itrack) const
{
// Get track from the basket.
   return gPropagator->fTracks[fIndex[itrack]];
}
   
//______________________________________________________________________________
void GeantVolumeBasket::GetWorkload(Int_t &indmin, Int_t &indmax)
{
// Get a range of indices to work with for a given thread.
   Int_t nthreads = gPropagator->fNthreads;
   TThread::Lock();
   indmin = indmax = fFirstFree;
   if (!fNtracks || fFirstFree==fNtracks) {
      TThread::UnLock();
      return;
   }
   Int_t fair_share = fNtracks/nthreads;
   Int_t remaining = fNtracks%nthreads;
   indmax = indmin+fair_share;
   if (remaining) indmax++;
   if (indmax > fNtracks) indmax = fNtracks;
   fFirstFree = indmax;
   TThread::UnLock();
}   

//______________________________________________________________________________
void GeantVolumeBasket::AddPendingTrack(Int_t itrack)
{
// Add a new track to the pending tracks list. This can be done in parallel
// by different threads.
   TThread::Lock();   
   if (!fPending) fPending = new Int_t[fMaxPending];
   if (fNpending==fMaxPending) {
      Int_t *pending = new Int_t[2*fMaxPending];
      memcpy(pending, fPending, fNpending*sizeof(Int_t));
      delete [] fPending;
      fPending = pending;
      fMaxPending *= 2;
   }
   fPending[fNpending++] = itrack;
   TThread::UnLock();     
}      

//______________________________________________________________________________
void GeantVolumeBasket::AddTrack(Int_t itrack)
{
// Add a track and its path to the basket.
   TThread::Lock();
   if (!fNtracks) {
      fIndex = new Int_t[fMaxTracks];
      fIndex[fNtracks] = itrack;
   }
   fIndex[fNtracks] = itrack;
   fNtracks++;
   // Increase arrays of tracks and path indices if needed
   if (fNtracks == fMaxTracks) {
      Int_t *newindex = new Int_t[2*fMaxTracks];
      memcpy(newindex, fIndex, fNtracks*sizeof(Int_t));
      delete [] fIndex;
      fIndex = newindex;
      fMaxTracks *= 2;
   }   
   TThread::UnLock();   
}     

//______________________________________________________________________________
void GeantVolumeBasket::Clear(Option_t *)
{
// Clear all particles and paths 
//   TThread::Lock();
   fNtracks = 0;
   fFirstFree = 0;
//   TThread::UnLock();   
}   

//______________________________________________________________________________
void GeantVolumeBasket::ComputeTransportLength(Int_t ntracks, Int_t *trackin)
{
// Computes snext and safety for an array of tracks. This is the transportation 
// process. Tracks are assumed to be inside gVolume.
   static Int_t icalls = 0;
   Double_t pdir[3];
   Int_t itr;
   Bool_t isOnBoundary = kFALSE;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   nav->SetOutside(kFALSE);
   
   for (itr=0; itr<ntracks; itr++) {
      GeantTrack *track = gPropagator->fTracks[trackin[itr]];
      track->Direction(pdir);
      track->path->UpdateNavigator(nav);
      nav->SetCurrentPoint(&track->xpos);
      nav->SetCurrentDirection(pdir);
      isOnBoundary = track->frombdr;
      Double_t pstep = TMath::Min(1.E20, track->pstep);
      nav->FindNextBoundary(pstep,"",isOnBoundary);
//      gGeoManager->SetVerboseLevel(0);
      track->safety = nav->GetSafeDistance();
//      if (!isOnBoundary && track->safety<gTolerance) {
//         nav->FindNextBoundary(track->pstep,"",isOnBoundary);
//      }
//      track->snext = nav->GetStep();
      track->snext = TMath::Max(gTolerance,nav->GetStep());
      if (gPropagator->fUseDebug && (gPropagator->fDebugTrk==trackin[itr] || gPropagator->fDebugTrk<0)) {
//         Printf("    %d   track %d: %s  snext=%19.15f safe=%19.15f pstep=%f", icalls,trackin[itr], nav->GetPath(), track->snext, track->safety, track->pstep);
         Printf("       track %d: %s  snext=%19.15f safe=%19.15f pstep=%f", trackin[itr], nav->GetPath(), track->snext, track->safety, track->pstep);
         track->Print(trackin[itr]);
      }   
   }
   icalls++;  
}            

//______________________________________________________________________________
void GeantVolumeBasket::ComputeTransportLengthSingle(Int_t *trackin)
{
// Computes snext and safety for a given track. This is the transportation 
// process. Track is assumed to be inside gVolume. Lighter version of normal 
// navigation, using only volume data and ignoring MANY option.
   static Int_t icalls = 0;
   Double_t pdir[3];
   Bool_t isOnBoundary = kFALSE;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   nav->SetOutside(kFALSE);
   GeantTrack *track = gPropagator->fTracks[trackin[0]];
   track->Direction(pdir);
   nav->SetCurrentPoint(&track->xpos);
   nav->SetCurrentDirection(pdir);
   isOnBoundary = track->frombdr;
   if (gPropagator->fUseDebug && (gPropagator->fDebugTrk==trackin[0] || gPropagator->fDebugTrk<0)) {
//         gGeoManager->SetVerboseLevel(10);
   }   
   Double_t pstep = TMath::Min(1.E20, track->pstep);
   nav->FindNextBoundary(pstep,"",isOnBoundary);
   gGeoManager->SetVerboseLevel(0);
   track->safety = nav->GetSafeDistance();
//   if (!isOnBoundary && track->safety<gTolerance) {
//      nav->FindNextBoundary(track->pstep,"",isOnBoundary);
//   }
//      track->snext = nav->GetStep();
   track->snext = TMath::Max(gTolerance,nav->GetStep());
   if (gPropagator->fUseDebug && (gPropagator->fDebugTrk==trackin[0] || gPropagator->fDebugTrk<0)) {
      Printf("       track %d: %s  snext=%19.15f safe=%19.15f pstep=%f", trackin[0], nav->GetPath(), track->snext, track->safety, track->pstep);
//      Printf("    %d   track %d: %s  snext=%19.15f safe=%19.15f pstep=%f", icalls,trackin[0], nav->GetPath(), track->snext, track->safety, track->pstep);
      track->Print(trackin[0]);
   }   
   icalls++;  
}            

//______________________________________________________________________________
void GeantVolumeBasket::PropagateTracks(Int_t ntracks, Int_t *trackin, Int_t &nout, Int_t *trackout, Int_t &ntodo, Int_t *tracktodo, Int_t &ncross, Int_t *trackcross)
{
// Propagate the ntracks with their selected steps. If a boundary is
// found in the way, the track is stopped. Nout must be initialized from outside.
//     trackin = array of <ntracks> input tracks
//     trackout = array of <nout> tracks propagated to physics processes
//     tracktodo = array of <ntodo> tracks propagated with a safe step or crossing 
//                 inside the same volume. These have to be propagated  again.
//     trackcross = array of <ncross> tracks that crossed the boundary. For these tracks
//                 the continuous processes have to be applied after propagation
   GeantTrack *track;
   Double_t step, snext, safety, c;
   ntodo = 0;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   GeantVolumeBasket *basket = 0;
//   Printf("===== PropagateTracks: ntracks=%d nout=%d", ntracks, nout);
   for (Int_t itr=0; itr<ntracks; itr++) {
      track = gPropagator->fTracks[trackin[itr]];
      // Skip neutral tracks for the time being (!!!)
      if (!track->charge) continue;
      if (!track->IsAlive()) continue;
      track->nsteps++;
      if (track->nsteps > gPropagator->fMaxSteps) {
         track->Kill();
         continue;
      }   
      step = track->pstep;
      snext = track->snext;
      safety = track->safety;
      c = track->Curvature();
      // If proposed step less than safety, just propagate to physics process
      if (step<safety) {
         track->izero = 0;
         track->PropagateInField(step, kFALSE, trackin[itr]);
         gPropagator->fNsafeSteps++; // increment-only, thread safe
         // track transported to physics process
         if (gPropagator->fUseDebug && (gPropagator->fDebugTrk==trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d process: %s", trackin[itr], gPropagator->Process(track->process)->GetName());
         trackout[nout++] = trackin[itr]; // <- survives geometry, stopped due to physics
         continue; // -> to next track
      }
      // Check if we can propagate to boundary
      if (0.25*c*snext<1E-6 && snext<1E-3 && snext<step-1E-6) {
         // Propagate with snext and check if we crossed
         //   backup track position and momentum
         if (track->izero>10) snext = 1.E-3;
         basket = track->PropagateInField(snext+10*gTolerance, kTRUE, trackin[itr]);
         if (snext<1.E-6) track->izero++;
         else track->izero = 0;
         gPropagator->fNsnextSteps++;
         if (!basket) {
            // track exiting
            if (nav->IsOutside()) {
               if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d exiting geometry", trackin[itr]);
               trackcross[ncross++] = trackin[itr];
               continue;
            }   
            // these tracks do not cross
//            if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d propagated with snext=%19.15f", trackin[itr], snext);
            tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
            continue; // -> next track
         }
         // These tracks are reaching boundaries
         trackcross[ncross++] = trackin[itr];
//         if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d pushed to boundary of %s", trackin[itr], basket->GetName());
//         basket = track->PropagateStraight(snext, trackin[itr]);
         continue; // -> to next track
      }
      // Track has safety<pstep but next boundary not close enough.
      // We propagate in field with the safety value.
      if (safety<gTolerance) {
         // Track getting away from boundary. Work to be done here
         // In principle we need a safety value for the content of the current volume only
         // This does not behave well on corners...
         // ... so we peek a small value and chech if this crosses, than recompute safety
         safety = 1.E-3;
         track->izero++;
         if (track->izero > 10) safety = 0.5*snext;
         basket = track->PropagateInField(safety, kTRUE, trackin[itr]);
      } else {
         if (track->izero > 10) {
            // Propagate with snext
            basket = track->PropagateInField(snext+10*gTolerance, kTRUE, trackin[itr]);
            track->izero = 0; 
            gPropagator->fNsnextSteps++;
            if (!basket) {
               // track exiting geometry
               if (nav->IsOutside()) {
                  if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d exiting geometry", trackin[itr]);
                  trackcross[ncross++] = trackin[itr];
                  continue;
               }   
               // these tracks do not cross
//               if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d propagated with snext=%19.15f", trackin[itr], snext);
               tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
               continue; // -> next track
            }
            // These tracks are reaching boundaries
            trackcross[ncross++] = trackin[itr];
            continue;
         }
         if (safety<1.E-3) track->izero++;
         // Propagate with safety without checking crossing
         basket = track->PropagateInField(safety, kFALSE, trackin[itr]);
      }   
      gPropagator->fNsafeSteps++;
      if (!basket) {
         // check if track exiting
         if (nav->IsOutside()) {
            if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d exiting geometry", trackin[itr]);
            trackcross[ncross++] = trackin[itr];
            continue;
         }   
         // these tracks do not cross
//         if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d propagated with safety=%19.15f", trackin[itr], safety);
         tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
         continue; // -> next track
      }
      // These tracks are reaching boundaries
      trackcross[ncross++] = trackin[itr];
   }   
   // Recompute snext and safety for todo tracks
   if (ntodo) ComputeTransportLength(ntodo, tracktodo);
}

//______________________________________________________________________________
Bool_t GeantVolumeBasket::PropagateTrack(Int_t *trackin)
{
// Propagate a single track with its selected physics step. If a boundary is
// found in the way return. 
   GeantTrack *track;
   Double_t step = 0;
   Double_t snext, safety, c;
//   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   track = gPropagator->fTracks[trackin[0]];
   // Skip neutral tracks for the moment (!!!)
   Bool_t crossed = kFALSE;
   // If proposed step less than safety, just propagate to physics process
   while (!crossed) {
      step = track->pstep;
      snext = track->snext;
      safety = track->safety;
      c = track->Curvature();
      if (step<safety) {
         track->izero = 0;
         crossed = track->PropagateInFieldSingle(step, kFALSE, trackin[0]);
         gPropagator->fNsafeSteps++;
         return crossed;
      }
      // Check if we can propagate to boundary
      if (0.25*c*snext<1E-6 && snext<1E-3 && snext<step-1E-6) {
         // Propagate with snext and check if we crossed
         if (track->izero>10) snext = 1.E-3;
         crossed = track->PropagateInFieldSingle(snext+10*gTolerance, kTRUE, trackin[0]);
         if (snext<1.E-6) track->izero++;
         else track->izero = 0;
         gPropagator->fNsnextSteps++;
         if (crossed) return crossed;
         ComputeTransportLengthSingle(trackin);
         continue;
      }   

      // Track has safety<pstep but next boundary not close enough.
      // We propagate in field with the safety value.
      if (safety<gTolerance) {
         // Track getting away from boundary. Work to be done here
         // In principle we need a safety value for the content of the current volume only
         // This does not behave well on corners...
         // ... so we peek a small value and chech if this crosses, than recompute safety
         safety = 1.E-3;
         track->izero++;
         if (track->izero > 10) safety = 0.5*snext;
         crossed = track->PropagateInFieldSingle(safety, kTRUE, trackin[0]);
      } else {        
         if (track->izero > 10) {
            // Propagate with snext
            crossed = track->PropagateInFieldSingle(snext+10*gTolerance, kTRUE, trackin[0]);
            track->izero = 0; 
            gPropagator->fNsnextSteps++;
            if (crossed) return crossed;
            ComputeTransportLengthSingle(trackin);
            continue;
         }
         if (safety<1.E-3) track->izero++;
         crossed = track->PropagateInField(safety, kFALSE, trackin[0]);
      } 
      gPropagator->fNsafeSteps++;
      if (crossed) return crossed;  
      // Recompute snext and safety for todo tracks
      ComputeTransportLengthSingle(trackin);
   }
   return crossed;
}

//______________________________________________________________________________
void GeantVolumeBasket::Prepare()
{
// Prepare basket for transport.
   if (!fNpending) return;
   if (!fIndex || fNpending > fMaxTracks) {
      fMaxTracks = fNpending;
      delete [] fIndex;
      fIndex = new Int_t[fMaxTracks];
   }
   memcpy(fIndex, fPending, fNpending*sizeof(Int_t));
   fNtracks = fNpending;
   fNpending = 0;   
}

//______________________________________________________________________________
void GeantVolumeBasket::Print(Option_t *) const
{
// Print info about the basket content.
   if (gPropagator->fUseDebug) {
      for (Int_t i=0; i<fNtracks; i++) {
//         if (gDebug && (gPropagator->fDebugTrk==fIndex[i] || gPropagator->fDebugTrk<0)) {
//            Printf("   %d - track %d: ", i, fIndex[i]);
            gPropagator->fTracks[fIndex[i]]->Print();
            if (gPropagator->fTracks[fIndex[i]]->path->GetNode(gPropagator->fTracks[fIndex[i]]->path->GetLevel())->GetVolume() != fVolume) {
               Printf("ERROR: Wrong path for basket %s", GetName());
               *((int*)0)=0;
            }  
//            GetBranchArray(i)->Print();
//         }   
      }
   }      
}

//______________________________________________________________________________
void GeantVolumeBasket::ResetStep(Int_t ntracks, Int_t *array)
{
// Reset current step for a list of tracks.
   for (Int_t i=0; i<ntracks; i++) {
      GeantTrack *track = gPropagator->fTracks[array[i]];
      track->step = 0.;
   }
}

//______________________________________________________________________________
void GeantVolumeBasket::Suspend(Int_t ntracks, Int_t *trackin)
{
// Suspend transport for this basket and put remaining tracks as pending.
   if (!ntracks) return;
   TThread::Lock();
   Clear();
   if (!fPending || fMaxPending<ntracks+fNpending) {
      Int_t *pending = new Int_t[ntracks+fNpending];
      if (fNpending) memcpy(pending, fPending, fNpending*sizeof(Int_t));
      delete [] fPending;
      fPending = pending;
   }
   memcpy(&fPending[fNpending], trackin, ntracks*sizeof(Int_t));
   fNpending += ntracks;
   TThread::UnLock();
}

//______________________________________________________________________________
void GeantVolumeBasket::TransportSingle()
{
// Transport all particles in this basket one by one (empty basket).
   if (!fNtracks) return;
   // Main loop
   Int_t itrack, nout;   
   Int_t nprocesses = gPropagator->fNprocesses;
   Int_t *particles = gPropagator->fPartInd[0]->GetArray();
   Int_t *partnext  = gPropagator->fPartNext[0]->GetArray();
   memcpy(particles, fIndex, fNtracks*sizeof(Int_t));
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
//   Int_t tid = nav->GetThreadId();
//   gPushed[tid]->ResetAllBits();
   Bool_t transported = kFALSE;
   Bool_t crossed = kFALSE;
   Bool_t usePhysics = gPropagator->fUsePhysics;
   Int_t generation = 0;
   TGeoBranchArray a;
   TGeoBranchArray b;
   Int_t n10 = fNtracks/10;
   for (itrack=0; itrack<fNtracks; itrack++) {
      if (n10) {
         if ((itrack%n10) == 0) Printf("%i percent", Int_t(100*itrack/fNtracks));
      }
      // Skip neutral particles for the moment (!!!)
      if (!gPropagator->fTracks[particles[itrack]]->charge) continue;
      transported = kFALSE;
      // Init navigator
      GetBranchArray(gPropagator->fTracks[particles[itrack]])->UpdateNavigator(nav);
      gPropagator->fVolume[0] = gGeoManager->GetCurrentVolume();
      // Physics step      
      if (usePhysics) gPropagator->PhysicsSelect(1,&particles[itrack],0);
      while (!transported) {
         crossed = kFALSE;
         a.InitFromNavigator(nav);
         // Loop all tracks to generate physics/geometry steps
         generation++;
         // Geometry snext and safety
         ComputeTransportLengthSingle(&particles[itrack]);
         crossed = PropagateTrack(&particles[itrack]);
         if (!crossed) {
            // Do post-step actions
            if (usePhysics) {
               // Apply continuous processes first
               for (Int_t iproc=0; iproc<nprocesses; iproc++) {
                  if (gPropagator->Process(iproc)->IsType(PhysicsProcess::kDiscrete)) continue;
                  nout = 0;
                  GeantPropagator::Instance()->Process(iproc)->PostStep(gPropagator->fVolume[0],1, &particles[itrack], nout, partnext,0);
                  if (!nout) {
                     transported = kTRUE;
                     break;
                  }   
               }
               // Apply discrete process if selected
               if (!transported && gPropagator->Process(gPropagator->fTracks[particles[itrack]]->process)->IsType(PhysicsProcess::kDiscrete)) {
                  nout = 0;
                  GeantPropagator::Instance()->Process(gPropagator->fTracks[particles[itrack]]->process)->PostStep(gPropagator->fVolume[0],1, &particles[itrack], nout, partnext,0);
//                  PostStepAction[gPropagator->fTracks[particles[itrack]]->process](1, &particles[itrack], nout, partnext,0);
                  if (!nout) {
                     transported = kTRUE;
                     break;
                  } 
               }     
               if (!transported) gPropagator->PhysicsSelect(1,&particles[itrack],0);
            }   
            continue;
         } else {
            // Energy deposition up to boundary
            if (usePhysics) {
               for (Int_t iproc=0; iproc<nprocesses; iproc++) {
                  if (gPropagator->Process(iproc)->IsType(PhysicsProcess::kDiscrete)) continue;
                  // In future take into account that particles may be generated 
                  // here (like dray)
                  Int_t nafter = 0;
                  GeantPropagator::Instance()->Process(iproc)->PostStep(gPropagator->fVolume[0],1, &particles[itrack], nafter, NULL,0);
//                  PostStepAction[iproc](1, &particles[itrack], nafter, NULL,0);
                  ResetStep(1, &particles[itrack]);
               }   
            }
         }   
         // Exit setup ?
         if (nav->IsOutside()) break;
         // Particle crossed boundary
         b.InitFromNavigator(nav);
         if (b==a) continue; // Particle entered the same volume, just continue the step
         // Physics step      
         gPropagator->fVolume[0] = gGeoManager->GetCurrentVolume();
         if (usePhysics) gPropagator->PhysicsSelect(1,&particles[itrack],0);
      }
   }
}   
