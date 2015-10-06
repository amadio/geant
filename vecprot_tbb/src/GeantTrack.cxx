#include "GeantTrack.h"

#include "TGeoBranchArray.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoHelix.h"

#include "globals.h"
#include "GeantPropagator.h"
#include "GeantVolumeBasket.h"
#include "GeantThreadData.h"
#include "WorkloadManager.h"

#include "base/Global.h"
using vecgeom::kTwoPi;

#include <iostream>

const double gTolerance = TGeoShape::Tolerance();

//______________________________________________________________________________
GeantTrack::GeantTrack(int ipdg)
    : event(-1), evslot(-1), particle(-1), pdg(ipdg), species(kHadron), status(kAlive), charge(0), mass(0), process(-1),
      xpos(0), ypos(0), zpos(0), px(0), py(0), pz(0), e(0), pstep(1.E20), step(0), snext(0), safety(0), frombdr(false),
      izero(0), nsteps(0), path(0), nextpath(0), pending(false) {
  // Constructor
  path = new TGeoBranchArray(30);
  nextpath = new TGeoBranchArray(30);
}

//______________________________________________________________________________
GeantTrack::GeantTrack(const GeantTrack &other)
    : event(other.event), evslot(other.evslot), particle(other.particle), pdg(other.pdg), species(other.species),
      status(other.status), charge(other.charge), mass(other.mass), process(other.process), xpos(other.xpos),
      ypos(other.ypos), zpos(other.zpos), px(other.px), py(other.py), pz(other.pz), e(other.e), pstep(other.pstep),
      step(other.step), snext(other.snext), safety(other.safety), frombdr(other.frombdr), izero(other.izero),
      nsteps(other.nsteps), path(new TGeoBranchArray(*other.path)), nextpath(new TGeoBranchArray(*other.nextpath)),
      pending(other.pending) {
  // Copy constructor
}

//______________________________________________________________________________
GeantTrack &GeantTrack::operator=(const GeantTrack &other) {
  // Assignment
  if (&other != this) {
    event = other.event;
    evslot = other.evslot;
    particle = other.particle;
    pdg = other.pdg;
    species = other.species;
    status = other.status;
    charge = other.charge;
    mass = other.mass;
    process = other.process;
    xpos = other.xpos;
    ypos = other.ypos;
    zpos = other.zpos;
    px = other.px;
    py = other.py;
    pz = other.pz;
    e = other.e;
    pstep = other.pstep;
    step = other.step;
    snext = other.snext;
    safety = other.safety;
    frombdr = other.frombdr;
    izero = other.izero;
    nsteps = other.nsteps;
    path = new TGeoBranchArray(*other.path);
    ;
    nextpath = new TGeoBranchArray(*other.nextpath);
    pending = other.pending;
  }
  return *this;
}

//______________________________________________________________________________
GeantTrack::~GeantTrack() {
  // Destructor.
  delete path;
  delete nextpath;
}

double GeantTrack::Curvature() const {
  GeantPropagator *propagator = GeantPropagator::Instance();
  return (charge) ? fabs(kB2C * propagator->fBmag / Pt()) : 0.;
}

//______________________________________________________________________________
void GeantTrack::Reset() {
  // Resets track content.
  event = -1;
  evslot = -1;
  particle = -1;
  pdg = 0;
  species = kHadron;
  status = kAlive;
  charge = 0;
  mass = 0;
  process = -1;
  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  px = 0.;
  py = 0.;
  pz = 0.;
  e = 0;
  pstep = 1.E20;
  step = 0.;
  snext = 0.;
  safety = 0.;
  frombdr = false;
  izero = 0;
  nsteps = 0;
  pending = false;
}

//______________________________________________________________________________
void GeantTrack::Direction(double dir[3]) {
  dir[0] = px;
  dir[1] = py;
  dir[2] = pz;
  double mdir = px * px + py * py + pz * pz;
  if (mdir > 0) {
    mdir = sqrt(1. / mdir);
    dir[0] * / mdir;
    dir[1] * / mdir;
    dir[2] * / mdir;
  }
}

//______________________________________________________________________________
void GeantTrack::Print(int) const {
  TString spath;
  //   if (path) path->GetPath(spath);
  Printf("=== Track %d (ev=%d): Process=%d, pstep=%g Charge=%d  Position:(%f,%f,%f) Mom:(%f,%f,%f) P:%g E:%g snext=%g "
         "safety=%g nsteps=%d",
         particle, event, process, pstep, charge, xpos, ypos, zpos, px, py, pz, sqrt(px * px + py * py + pz * pz), e,
         snext, safety, nsteps);
}

//______________________________________________________________________________
GeantVolumeBasket *GeantTrack::PropagateStraight(double crtstep, int itr) {
  // Propagate along a straight line for neutral particles, for B=0 or for last tiny step.
  // The method adds the particle to the next volume basket.
  // Returns the basket pointer, null if exiting geometry.

  GeantPropagator *gPropagator = GeantPropagator::Instance();
  PerThread::reference TBBperThread = gPropagator->fTBBthreadData.local();

  // Find next volume
  double dir[3];
  //   frombdr = kTRUE;
  Direction(dir);
  pstep -= crtstep;
  safety = 0;
  // Change path to reflect the physical volume for the current track;
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();

  if (!nav) {
    nav = gGeoManager->AddNavigator();
    std::cerr << "[PropagateStraight] Added navigator" << std::endl;
  }

  /*int tid = nav->GetThreadId();*/
  // Particle crossed?
  if (frombdr)
    nextpath->UpdateNavigator(nav);
  else
    path->UpdateNavigator(nav);
  nav->SetOutside(kFALSE);
  nav->SetStep(crtstep);
  xpos += crtstep * dir[0];
  ypos += crtstep * dir[1];
  zpos += crtstep * dir[2];
  TGeoVolume *vol = nav->GetCurrentVolume();
  if (vol->IsAssembly())
    Printf("### ERROR ### Entered assembly %s", vol->GetName());
  GeantVolumeBasket *basket = (GeantVolumeBasket *)vol->GetField();
  if (frombdr) {
    TBBperThread.fCollection->AddTrack(itr, basket);
    path->InitFromNavigator(nav);
  }
  if (gPropagator->fUseDebug && (gPropagator->fDebugTrk == itr || gPropagator->fDebugTrk < 0)) {
    Printf("   track %d: entering %s at:(%f, %f, %f)", itr, nav->GetCurrentNode()->GetName(), xpos, ypos, zpos);
    if (gPropagator->fTracks[itr]->path)
      gPropagator->fTracks[itr]->path->Print();
  }
  //   basket->Print();
  return basket;
}

//______________________________________________________________________________
GeantVolumeBasket *GeantTrack::PropagateInField(double crtstep, bool checkcross, int itr) {
  // Propagate with step using the helix propagator. Returns a basket if a
  // boundary was crossed. In such case, the track position and step will reflect
  // the boundary crossing point.
  GeantPropagator *gPropagator = GeantPropagator::Instance();
  PerThread::reference TBBperThread = gPropagator->fTBBthreadData.local();

  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  /*int tid = nav->GetThreadId();*/

  if (!nav) {
    nav = gGeoManager->AddNavigator();
    std::cerr << "[PropagateInField] Added navigator" << std::endl;
  }

  bool useDebug = gPropagator->fUseDebug;
  int debugTrk = gPropagator->fDebugTrk;
  nav->ResetState();
  if (checkcross) {
    path->UpdateNavigator(nav);
    nav->SetLastSafetyForPoint(safety, &xpos);
  }
  // Reset relevant variables
  frombdr = kFALSE;
  pstep -= crtstep;
  safety -= crtstep;
  if (safety < 0.)
    safety = 0.;
  step += crtstep;
  // Set curvature, charge
  double c = 0.;
  double dir[3];
  double ptot = 0;
  const double *point = 0;
  const double *newdir = 0;
  Direction(dir);
  if (charge) {
    c = Curvature();
    TBBperThread.fFieldPropagator->SetXYcurvature(c);
    TBBperThread.fFieldPropagator->SetCharge(charge);
    TBBperThread.fFieldPropagator->SetHelixStep(fabs(kTwoPi * pz / (c * Pt())));
    TBBperThread.fFieldPropagator->InitPoint(xpos, ypos, zpos);
    TBBperThread.fFieldPropagator->InitDirection(dir);
    TBBperThread.fFieldPropagator->UpdateHelix();
    TBBperThread.fFieldPropagator->Step(crtstep);
    point = TBBperThread.fFieldPropagator->GetCurrentPoint();
    newdir = TBBperThread.fFieldPropagator->GetCurrentDirection();
    xpos = point[0];
    ypos = point[1];
    zpos = point[2];
    ptot = P();
    px = ptot * newdir[0];
    py = ptot * newdir[1];
    pz = ptot * newdir[2];
  } else {
    xpos += crtstep * dir[0];
    ypos += crtstep * dir[1];
    zpos += crtstep * dir[2];
  }
  if (!checkcross)
    return 0;
  if (charge) {
    if (nav->IsSameLocation(xpos, ypos, zpos, kTRUE))
      return 0;
  } else {
    frombdr = kTRUE;
    TGeoNode *skip = nav->GetCurrentNode();
    TGeoNode *next = 0;
    next = nav->CrossBoundaryAndLocate(kTRUE, skip);
    if (!next && !nav->IsOutside()) {
      path->UpdateNavigator(nav);
      next = nav->FindNextBoundaryAndStep(TGeoShape::Big(), kFALSE);
    }
  }

  // Boundary crossed
  TGeoNode *checked = nav->GetCurrentNode();
  TGeoVolume *vol = checked->GetVolume();
  double ldir[3], ld[3];
  double local[3], lp[3];
  double delta;
  bool outside = nav->IsOutside();
  // Swap track direction and compute distance back to boundary
  dir[0] = -newdir[0];
  dir[1] = -newdir[1];
  dir[2] = -newdir[2];
  int level = nav->GetLevel();
  bool entering = kTRUE;
  TGeoNode *node1 = 0;
  TGeoNode *node2 = 0;
  if (level < path->GetLevel() && !outside) {
    for (int lev = 0; lev <= level; lev++) {
      node1 = nav->GetMother(level - lev);
      node2 = path->GetNode(lev);
      if (node1 == node2) {
        if (lev == level)
          entering = kFALSE;
      } else {
        // different nodes at some level -> entering current node
        break;
      }
    }
  }
  if (!entering) {
    checked = path->GetNode(level + 1);
    if (!checked)
      return 0;
  }
  nav->MasterToLocal(&xpos, local);
  nav->MasterToLocalVect(dir, ldir);
  if (entering) {
    if (outside)
      delta = vol->GetShape()->DistFromOutside(local, ldir, 3);
    else
      delta = vol->GetShape()->DistFromInside(local, ldir, 3);
  } else {
    checked->MasterToLocal(local, lp);
    checked->MasterToLocalVect(ldir, ld);
    delta = checked->GetVolume()->GetShape()->DistFromOutside(lp, ld, 3);
  }
  if (useDebug && (debugTrk < 0 || itr == debugTrk)) {
    if (entering)
      Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr,
             vol->GetName(), xpos, ypos, zpos, crtstep, delta);
    else
      Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(),
             vol->GetName(), crtstep, delta);
  }
  if (delta > crtstep) {
    if (useDebug && (debugTrk < 0 || itr == debugTrk)) {
      if (entering)
        Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr,
               vol->GetName(), xpos, ypos, zpos, crtstep, delta);
      else
        Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(),
               vol->GetName(), crtstep, delta);
      Printf("Error propagating back %19.15f (forward %19.15f) track for track %d", delta, crtstep, itr);
      if (entering)
        Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", vol->GetName(), local[0],
               local[1], local[2], ldir[0], ldir[1], ldir[2]);
      else
        Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", checked->GetName(), lp[0],
               lp[1], lp[2], ld[0], ld[1], ld[2]);
    }
    delta = 10 * gTolerance;
  }
  delta -= 10 * gTolerance;
  // Propagate back to boundary and store position/direction
  xpos += delta * dir[0];
  ypos += delta * dir[1];
  zpos += delta * dir[2];
  px = -ptot * dir[0];
  py = -ptot * dir[1];
  pz = -ptot * dir[2];
  xpos = point[0];
  ypos = point[1];
  zpos = point[2];
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
  if (vol->IsAssembly())
    Printf("### ERROR ### Entered assembly %s", vol->GetName());
  GeantVolumeBasket *basket = (GeantVolumeBasket *)vol->GetField();
  TBBperThread.fCollection->AddTrack(itr, basket);
  //   if (gUseDebug && (gDebugTrk==itr || gDebugTrk<0)) {
  //      Printf("   track %d: entering %s at:(%f, %f, %f)   ", itr, nav->GetCurrentNode()->GetName(), xpos,ypos,zpos);
  //      basket->GetBranchArray(gTracks[itr])->Print();
  //   }
  return basket;
}
