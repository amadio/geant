#include "globals.h"
#include "TGeoBranchArray.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoHelix.h"
#include "GeantTrack.h"
#include "GeantVolumeBasket.h"
#include "WorkloadManager.h"

const Double_t gTolerance = TGeoShape::Tolerance();

//______________________________________________________________________________
GeantTrack::GeantTrack(Int_t ipdg) 
           :event(-1),
            particle(-1),
            pdg(ipdg),
            species(kHadron),
            status(kAlive),
            charge(0),
            mass(0),
            process(-1),
            xpos(0),
            ypos(0),
            zpos(0),
            px(0),
            py(0),
            pz(0),
            e(0), 
            pstep(1.E20), 
            step(0), 
            snext(0), 
            safety(0), 
            frombdr(false), 
            izero(0), 
            nsteps(0),
            path(0),
            nextpath(0),
            pending(false)
{
// Constructor
   path = new TGeoBranchArray(30);
   nextpath = new TGeoBranchArray(30);
}

//______________________________________________________________________________
GeantTrack::~GeantTrack()
{
// Destructor.
   delete path;
   delete nextpath;
}   

//______________________________________________________________________________
void GeantTrack::Direction(Double_t dir[3]) {
   dir[0] = px; dir[1] = py; dir[2] = pz;
   TMath::Normalize(dir);
}

//______________________________________________________________________________
void GeantTrack::Print(Int_t) const {
   TString spath;
//   if (path) path->GetPath(spath);
   Printf("=== Track %d (%s): Process=%d, pstep=%g Charge=%d  Position:(%f,%f,%f) Mom:(%f,%f,%f) P:%g E:%g snext=%g safety=%g nsteps=%d",
           particle,spath.Data(), process,pstep,charge,xpos,ypos,zpos,px,py,pz,TMath::Sqrt(px*px+py*py+pz*pz),e,snext,safety,nsteps);
}

//______________________________________________________________________________
GeantVolumeBasket *GeantTrack::PropagateStraight(Double_t crtstep, Int_t itr)
{
// Propagate along a straight line for neutral particles, for B=0 or for last tiny step.
// The method adds the particle to the next volume basket. 
// Returns the basket pointer, null if exiting geometry.
// Find next volume
   Double_t dir[3];
//   frombdr = kTRUE;
   Direction(dir);
   pstep -= crtstep;
   safety = 0;
   // Change path to reflect the physical volume for the current track;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   Int_t tid = nav->GetThreadId();
   // Particle crossed?
   if (frombdr) nextpath->UpdateNavigator(nav);
   else path->UpdateNavigator(nav);
   nav->SetOutside(kFALSE);
   nav->SetStep(crtstep);
   xpos += crtstep*dir[0];
   ypos += crtstep*dir[1];
   zpos += crtstep*dir[2];
   TGeoVolume *vol = nav->GetCurrentVolume();
   if (vol->IsAssembly()) Printf("### ERROR ### Entered assembly %s", vol->GetName());
   GeantVolumeBasket *basket = (GeantVolumeBasket*)vol->GetField();
   if (frombdr) {
      gPropagator->fCollections[tid]->AddTrack(itr,basket);
      path->InitFromNavigator(nav);
   }   
   // Signal that the transport is still ongoing if the particle entered a new basket
   gPropagator->fTransportOngoing = kTRUE;
   if (gPropagator->fUseDebug && (gPropagator->fDebugTrk==itr || gPropagator->fDebugTrk<0)) {
      Printf("   track %d: entering %s at:(%f, %f, %f)", itr, nav->GetCurrentNode()->GetName(), xpos,ypos,zpos);
      if (gPropagator->fTracks[itr]->path) gPropagator->fTracks[itr]->path->Print();
   }
//   basket->Print();
   return basket;  
}   

//______________________________________________________________________________
GeantVolumeBasket *GeantTrack::PropagateInField(Double_t crtstep, Bool_t checkcross, Int_t itr)
{
// Propagate with step using the helix propagator. Returns a basket if a 
// boundary was crossed. In such case, the track position and step will reflect
// the boundary crossing point.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   Int_t tid = nav->GetThreadId();
   Bool_t useDebug = gPropagator->fUseDebug;
   Int_t debugTrk = gPropagator->fDebugTrk;
   nav->ResetState();
   if (checkcross) {
      path->UpdateNavigator(nav);
      nav->SetLastSafetyForPoint(safety, &xpos);
   }   
   // Reset relevant variables
   frombdr = kFALSE;
   pstep -= crtstep;
   safety -= crtstep;
   if (safety<0.) safety = 0.;
   step += crtstep;
   // Set curvature, charge
   Double_t c = 0.;
   Double_t dir[3];
   Double_t ptot = 0;
   const Double_t *point = 0;
   const Double_t *newdir = 0;
   Direction(dir);
   if (charge) {
      c = Curvature();
      gPropagator->fFieldPropagator[tid]->SetXYcurvature(c);
      gPropagator->fFieldPropagator[tid]->SetCharge(charge);
      gPropagator->fFieldPropagator[tid]->SetHelixStep(TMath::Abs(TMath::TwoPi()*pz/(c*Pt())));
      gPropagator->fFieldPropagator[tid]->InitPoint(xpos,ypos,zpos);
      gPropagator->fFieldPropagator[tid]->InitDirection(dir);
      gPropagator->fFieldPropagator[tid]->UpdateHelix();
      gPropagator->fFieldPropagator[tid]->Step(crtstep);
      point = gPropagator->fFieldPropagator[tid]->GetCurrentPoint();
      newdir = gPropagator->fFieldPropagator[tid]->GetCurrentDirection();
      xpos = point[0]; ypos = point[1]; zpos = point[2];
      ptot = P();
      px = ptot*newdir[0];
      py = ptot*newdir[1];
      pz = ptot*newdir[2];
   } else {
      xpos += crtstep*dir[0];
      ypos += crtstep*dir[1];
      zpos += crtstep*dir[2];
   }
   if (!checkcross) return 0;
   if (charge) {
      if (nav->IsSameLocation(xpos,ypos,zpos,kTRUE)) return 0;
   } else {
      frombdr = kTRUE;
      TGeoNode *skip = nav->GetCurrentNode();
      TGeoNode *next = 0;
      next = nav->CrossBoundaryAndLocate(kTRUE, skip);
      if (!next && !nav->IsOutside()) {
         path->UpdateNavigator(nav);
         next = nav->FindNextBoundaryAndStep(TGeoShape::Big(),kFALSE);
      }   
   }   
      
   // Boundary crossed
   TGeoNode *checked = nav->GetCurrentNode();
   TGeoVolume *vol = checked->GetVolume();
   Double_t ldir[3], ld[3];
   Double_t local[3], lp[3];
   Double_t delta;
   Bool_t outside = nav->IsOutside();
   // Swap track direction and compute distance back to boundary
   dir[0] = -newdir[0]; dir[1] = -newdir[1]; dir[2] = -newdir[2];
   Int_t level = nav->GetLevel();
   Bool_t entering = kTRUE;
   TGeoNode *node1 = 0;
   TGeoNode *node2 = 0;
   if (level < path->GetLevel() && !outside) {
      for (Int_t lev=0; lev<=level; lev++) {
         node1 = nav->GetMother(level-lev);
         node2 = path->GetNode(lev);
         if (node1 == node2) {
            if (lev==level) entering = kFALSE;
         } else {
         // different nodes at some level -> entering current node
            break;
         }   
      }
   }   
   if (!entering) {
      checked = path->GetNode(level+1);
      if (!checked) return 0;
   }   
   nav->MasterToLocal(&xpos, local);
   nav->MasterToLocalVect(dir, ldir);
   if (entering) {
      if (outside) delta = vol->GetShape()->DistFromOutside(local,ldir,3);
      else         delta = vol->GetShape()->DistFromInside(local,ldir,3);
   } else {
      checked->MasterToLocal(local,lp);
      checked->MasterToLocalVect(ldir,ld);
      delta = checked->GetVolume()->GetShape()->DistFromOutside(lp,ld,3);
   }   
   if (useDebug && (debugTrk<0 || itr==debugTrk)) {
      if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
      else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
   }   
   if (delta>crtstep) {
      if (useDebug && (debugTrk<0 || itr==debugTrk)) {
         if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
         else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
         Printf("Error propagating back %19.15f (forward %19.15f) track for track %d", delta, crtstep, itr);
         if (entering) Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", vol->GetName(), local[0],local[1],local[2],ldir[0],ldir[1],ldir[2]);
         else          Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", checked->GetName(), lp[0],lp[1],lp[2],ld[0],ld[1],ld[2]);
      }   
      delta = 10*gTolerance;
   }
   delta -= 10*gTolerance;
   // Propagate back to boundary and store position/direction
   xpos += delta*dir[0];
   ypos += delta*dir[1];
   zpos += delta*dir[2]; 
   px = -ptot*dir[0];
   py = -ptot*dir[1];
   pz = -ptot*dir[2];
   xpos = point[0]; ypos = point[1]; zpos = point[2];
   // Mark track as "on boundary" and update step/pstep
   frombdr = kTRUE;
   safety = 0.;
   pstep += delta;
   step -= delta;
   if (outside) {
      nav->SetOutside(kTRUE);
      return 0;
   }   
   // Create a new branch array
   path->InitFromNavigator(nav);
   if (vol->IsAssembly()) Printf("### ERROR ### Entered assembly %s", vol->GetName());
   GeantVolumeBasket *basket = (GeantVolumeBasket*)vol->GetField();
   gPropagator->fCollections[tid]->AddTrack(itr, basket);
   // Signal that the transport is still ongoing if the particle entered a new basket
   gPropagator->fTransportOngoing = kTRUE;
//   if (gUseDebug && (gDebugTrk==itr || gDebugTrk<0)) {
//      Printf("   track %d: entering %s at:(%f, %f, %f)   ", itr, nav->GetCurrentNode()->GetName(), xpos,ypos,zpos);
//      basket->GetBranchArray(gTracks[itr])->Print();
//   }
   return basket;

}   
//______________________________________________________________________________
Bool_t GeantTrack::PropagateInFieldSingle(Double_t crtstep, Bool_t checkcross, Int_t itr)
{
// Propagate with step using the helix propagator. Navigation must point to
// current particle location. Returns true if crossed next boundary.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   Bool_t useDebug = gPropagator->fUseDebug;
   Int_t debugTrk = gPropagator->fDebugTrk;
   Int_t tid = nav->GetThreadId();
   nav->ResetState();
   TGeoBranchArray a;
   a.InitFromNavigator(nav);
   if (checkcross) {
      nav->SetLastSafetyForPoint(safety, &xpos);
   }   
   // Reset relevant variables
   frombdr = kFALSE;
   pstep -= crtstep;
   safety -= crtstep;
   if (safety<0.) safety = 0.;
   step += crtstep;
   // Set curvature, charge
   Double_t c = Curvature();
   gPropagator->fFieldPropagator[tid]->SetXYcurvature(c);
   gPropagator->fFieldPropagator[tid]->SetCharge(charge);
   gPropagator->fFieldPropagator[tid]->SetHelixStep(TMath::Abs(TMath::TwoPi()*pz/(c*Pt())));
   gPropagator->fFieldPropagator[tid]->InitPoint(xpos,ypos,zpos);
   Double_t dir[3];
   Direction(dir);
   gPropagator->fFieldPropagator[tid]->InitDirection(dir);
   gPropagator->fFieldPropagator[tid]->UpdateHelix();
   gPropagator->fFieldPropagator[tid]->Step(crtstep);
   const Double_t *point = gPropagator->fFieldPropagator[tid]->GetCurrentPoint();
   const Double_t *newdir = gPropagator->fFieldPropagator[tid]->GetCurrentDirection();
   xpos = point[0]; ypos = point[1]; zpos = point[2];
   Double_t ptot = P();
   px = ptot*newdir[0];
   py = ptot*newdir[1];
   pz = ptot*newdir[2];
   if (!checkcross) return kFALSE;
   if (nav->IsSameLocation(xpos,ypos,zpos,kTRUE)) return kFALSE;
   // Boundary crossed
   TGeoNode *checked = nav->GetCurrentNode();
   TGeoVolume *vol = checked->GetVolume();
   Double_t ldir[3], ld[3];
   Double_t local[3], lp[3];
   Double_t delta;
   Bool_t outside = nav->IsOutside();
   // Swap track direction and compute distance back to boundary
   dir[0] = -newdir[0]; dir[1] = -newdir[1]; dir[2] = -newdir[2];
   Int_t level = nav->GetLevel();
   Bool_t entering = kTRUE;
   TGeoNode *node1 = 0;
   TGeoNode *node2 = 0;
   if (level < a.GetLevel() && !outside) {
      for (Int_t lev=0; lev<=level; lev++) {
         node1 = nav->GetMother(level-lev);
         node2 = a.GetNode(lev);
         if (node1 == node2) {
            if (lev==level) entering = kFALSE;
         } else {
         // different nodes at some level -> entering current node
            break;
         }   
      }
   }
   if (!entering) checked = a.GetNode(level+1);
   nav->MasterToLocal(&xpos, local);
   nav->MasterToLocalVect(dir, ldir);
   if (entering) {
      if (outside) delta = vol->GetShape()->DistFromOutside(local,ldir,3);
      else         delta = vol->GetShape()->DistFromInside(local,ldir,3);
   } else {
      checked->MasterToLocal(local,lp);
      checked->MasterToLocalVect(ldir,ld);
      delta = checked->GetVolume()->GetShape()->DistFromOutside(lp,ld,3);
   }   
   if (useDebug && (debugTrk<0 || itr==debugTrk)) {
      if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
      else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
   }   
   if (delta>crtstep) {
      if (useDebug && (debugTrk<0 || itr==debugTrk)) {
         if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
         else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
         Printf("Error propagating back %19.15f (forward %19.15f) track for track %d", delta, crtstep, itr);
         if (entering) Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", vol->GetName(), local[0],local[1],local[2],ldir[0],ldir[1],ldir[2]);
         else          Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", checked->GetName(), lp[0],lp[1],lp[2],ld[0],ld[1],ld[2]);
      }   
      delta = 10*gTolerance;
   }
   delta -= 10*gTolerance;
   // Propagate back to boundary and store position/direction
   xpos += delta*dir[0];
   ypos += delta*dir[1];
   zpos += delta*dir[2]; 
   px = -ptot*dir[0];
   py = -ptot*dir[1];
   pz = -ptot*dir[2];
   xpos = point[0]; ypos = point[1]; zpos = point[2];
   // Mark track as "on boundary" and update step/pstep
   frombdr = kTRUE;
   safety = 0.;
   pstep += delta;
   step -= delta;
   if (outside) {
      nav->SetOutside(kTRUE);
   }
   // Boundary crossed, navigation points to new location
   return kTRUE;
}   

//______________________________________________________________________________
void GeantEvent::AddTrack()
{
// Thread safe track addition
   the_mutex.Lock();
   ntracks++;
   the_mutex.UnLock();
}

//______________________________________________________________________________
void GeantEvent::StopTrack()
{
// Thread safe track addition
   the_mutex.Lock();
   ndone++;
   the_mutex.UnLock();
}

   
