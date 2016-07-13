#include "GeantTrack.h"

#include "globals.h"
#include "Geant/Error.h"
#include <execinfo.h>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GeantTrack::GeantTrack()
    : fEvent(-1), fEvslot(-1), fParticle(-1), fMother(0), fPDG(0), fGVcode(0), fEindex(0), fCharge(0), fProcess(-1),
      fNsteps(0), fSpecies(kHadron), fStatus(kAlive), fMass(0), fXpos(0), fYpos(0), fZpos(0), fXdir(0), fYdir(0),
      fZdir(0), fP(0), fE(0), fTime(0), fEdep(0), fPstep(1.E20), fStep(0), fSnext(0), fSafety(0), fNintLen(0), fIntLen(0), 
      fBoundary(false), fPending(false), fOwnPath(true), fPath(0), fNextpath(0) {
  // Dummy constructor
}

/* Obtain a backtrace and print it to stdout. */
void printrace(void) {
  void *array[10];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace(array, 10);
  strings = backtrace_symbols(array, size);

  printf("Obtained %zd stack frames.\n", size);

  for (i = 0; i < size; i++)
    printf("%s\n", strings[i]);

  free(strings);
}

//______________________________________________________________________________
GeantTrack::GeantTrack(int ipdg)
    : fEvent(-1), fEvslot(-1), fParticle(-1), fMother(0), fPDG(ipdg), fGVcode(0), fEindex(0), fCharge(0), fProcess(-1),
      fNsteps(0), fSpecies(kHadron), fStatus(kAlive), fMass(0), fXpos(0), fYpos(0), fZpos(0), fXdir(0), fYdir(0),
      fZdir(0), fP(0), fE(0), fTime(0), fEdep(0), fPstep(1.E20), fStep(0), fSnext(0), fSafety(0), fNintLen(0), fIntLen(0),
      fBoundary(false), fPending(false), fOwnPath(true), fPath(0), fNextpath(0) {
  // Constructor
  int maxdepth = GeantPropagator::Instance()->fMaxDepth;
  fPath = VolumePath_t::MakeInstance(maxdepth);
  fNextpath = VolumePath_t::MakeInstance(maxdepth);
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
GeantTrack::GeantTrack(int ipdg, int maxdepth)
    : fEvent(-1), fEvslot(-1), fParticle(-1), fMother(0), fPDG(ipdg), fGVcode(0), fEindex(0), fCharge(0), fProcess(-1),
      fNsteps(0), fSpecies(kHadron), fStatus(kAlive), fMass(0), fXpos(0), fYpos(0), fZpos(0), fXdir(0), fYdir(0),
      fZdir(0), fP(0), fE(0), fTime(0), fEdep(0), fPstep(1.E20), fStep(0), fSnext(0), fSafety(0), fNintLen(0), fIntLen(0),
      fBoundary(false), fPending(false), fOwnPath(true), fPath(0), fNextpath(0) {
  // Constructor
  fPath = VolumePath_t::MakeInstance(maxdepth);
  fNextpath = VolumePath_t::MakeInstance(maxdepth);
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
GeantTrack::GeantTrack(void *addr, int maxdepth)
    : fEvent(-1), fEvslot(-1), fParticle(-1), fPDG(0), fGVcode(0), fEindex(0), fCharge(0), fProcess(-1),
      fNsteps(0), fSpecies(kHadron), fStatus(kAlive), fMass(0), fXpos(0), fYpos(0), fZpos(0), fXdir(0), fYdir(0),
      fZdir(0), fP(0), fE(0), fTime(0), fEdep(0), fPstep(1.E20), fStep(0), fSnext(0), fSafety(0), fNintLen(0), fIntLen(0),
      fBoundary(false), fPending(false), fOwnPath(true), fPath(nullptr), fNextpath(nullptr) {
  // In place private constructor
  char *path_addr = round_up_align((char*)addr + sizeof(GeantTrack));
  fPath = VolumePath_t::MakeInstanceAt(maxdepth, path_addr);
  path_addr += round_up_align(VolumePath_t::SizeOfInstance(maxdepth));
  fNextpath = VolumePath_t::MakeInstanceAt(maxdepth, path_addr);
}

//______________________________________________________________________________
GeantTrack::GeantTrack(const GeantTrack &other)
    : fEvent(other.fEvent), fEvslot(other.fEvslot), fParticle(other.fParticle), fMother(other.fMother), fPDG(other.fPDG),
      fGVcode(other.fGVcode), fEindex(other.fEindex), fCharge(other.fCharge), fProcess(other.fProcess),
      fNsteps(other.fNsteps), fSpecies(other.fSpecies), fStatus(other.fStatus),
      fMass(other.fMass), fXpos(other.fXpos), fYpos(other.fYpos), fZpos(other.fZpos), fXdir(other.fXdir),
      fYdir(other.fYdir), fZdir(other.fZdir), fP(other.fP), fE(other.fE), fTime(other.fTime), fEdep(other.fEdep),
      fPstep(other.fPstep), fStep(other.fStep), fSnext(other.fSnext), fSafety(other.fSafety), fNintLen(other.fNintLen), fIntLen(other.fIntLen),
      fBoundary(other.fBoundary), fPending(other.fPending), fOwnPath(true), fPath(0), fNextpath(0) {
  // Copy constructor
  int maxdepth = GeantPropagator::Instance()->fMaxDepth;
  fPath = VolumePath_t::MakeInstance(maxdepth);
  fNextpath = VolumePath_t::MakeInstance(maxdepth);
  *fPath = *other.fPath;
  *fNextpath = *other.fNextpath;
}

//______________________________________________________________________________
GeantTrack &GeantTrack::operator=(const GeantTrack &other) {
  // Assignment
  if (&other != this) {
    fEvent = other.fEvent;
    fEvslot = other.fEvslot;
    fParticle = other.fParticle;
    fMother = other.fMother;
    fPDG = other.fPDG;
    fGVcode = other.fGVcode;
    fEindex = other.fEindex;
    fCharge = other.fCharge;
    fProcess = other.fProcess;
    fNsteps = other.fNsteps;
    fSpecies = other.fSpecies;
    fStatus = other.fStatus;
    fMass = other.fMass;
    fXpos = other.fXpos;
    fYpos = other.fYpos;
    fZpos = other.fZpos;
    fXdir = other.fXdir;
    fYdir = other.fYdir;
    fZdir = other.fZdir;
    fP = other.fP;
    fE = other.fE;
    fTime = other.fTime;
    fEdep = other.fEdep;
    fPstep = other.fPstep;
    fStep = other.fStep;
    fSnext = other.fSnext;
    fSafety = other.fSafety;
    fNintLen = other.fNintLen;
    fIntLen = other.fIntLen;
    fBoundary = other.fBoundary;
    fPending = other.fPending;
    int maxdepth = GeantPropagator::Instance()->fMaxDepth;
    fPath = VolumePath_t::MakeInstance(maxdepth);
    fNextpath = VolumePath_t::MakeInstance(maxdepth);
    *fPath = *other.fPath;
    *fNextpath = *other.fNextpath;
  }
  return *this;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
GeantTrack::~GeantTrack() {
  // Destructor.
  if (fOwnPath) {
    VolumePath_t::ReleaseInstance(fPath);
    VolumePath_t::ReleaseInstance(fNextpath);
  }
}

//______________________________________________________________________________
void GeantTrack::Clear(const char *) {
  // Resets track content.
  fEvent = -1;
  fEvslot = -1;
  fParticle = -1;
  fPDG = 0;
  fGVcode = 0;
  fEindex = 0;
  fCharge = 0;
  fProcess = -1;
  fNsteps = 0;
  fSpecies = kHadron;
  fStatus = kAlive;
  fMass = 0.;
  fXpos = 0.;
  fYpos = 0.;
  fZpos = 0.;
  fXdir = 0.;
  fYdir = 0.;
  fZdir = 0.;
  fP = 0.;
  fE = 0.;
  fTime = 0.;
  fEdep = 0;
  fPstep = 1.E20;
  fStep = 0.;
  fSnext = 0.;
  fSafety = 0.;
  fNintLen = 0.;
  fIntLen = 0.;
  fBoundary = false;
  fPending = false;
#ifdef USE_VECGEOM_NAVIGATOR
  fPath->Clear();
  fNextpath->Clear();
#endif
}

//______________________________________________________________________________
Volume_t const*GeantTrack::GetVolume() const {
#ifdef USE_VECGEOM_NAVIGATOR
  return fPath->Top()->GetLogicalVolume();
#else
  return (fPath->GetCurrentNode()->GetVolume());
#endif
}

//______________________________________________________________________________
Volume_t const*GeantTrack::GetNextVolume() const {
#ifdef USE_VECGEOM_NAVIGATOR
  // Next volume the track is entering
  return fNextpath->Top()->GetLogicalVolume();
#else
  // Next volume the track is entering
  return fNextpath->GetCurrentNode()->GetVolume();
#endif
}

//______________________________________________________________________________
Material_t *GeantTrack::GetMaterial() const {
   // Current material the track is into
#ifdef USE_VECGEOM_NAVIGATOR
   auto med = (Medium_t *) GetVolume()->GetTrackingMediumPtr();
#else
   auto med = GetVolume()->GetMedium();
#endif
   if (!med)
      return 0;
   return med->GetMaterial();
}


//______________________________________________________________________________
bool GeantTrack::IsNormalized(double tolerance) const {
  // Check if track direction is normalized within tolerance
  double norm = fXdir * fXdir + fYdir * fYdir + fZdir * fZdir;
  if (fabs(1. - norm) > tolerance)
    return false;
  return true;
}

//______________________________________________________________________________
void GeantTrack::SetPath(VolumePath_t const *const path) {
  // Set path.
  *fPath = *path;
}

//______________________________________________________________________________
void GeantTrack::SetNextPath(VolumePath_t const *const path) {
  // Set next path.
  *fNextpath = *path;
}

//______________________________________________________________________________
void GeantTrack::Print(const char *msg) const {
  const char *status[8] = {"alive", "killed", "inflight", "boundary", "exitSetup", "physics", "postponed", "new"};
#ifdef USE_VECGEOM_NAVIGATOR
  Geant::Print(msg,
      ": evt=%d slt=%d part=%d pdg=%d gvc=%d eind=%d bind=%d chg=%d proc=%d nstp=%d spc=%d status=%s mass=%g "
      "xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g time=%g pstp=%g stp=%g snxt=%g saf=%g nil=%g ile=%g bdr=%d\n",
      fEvent, fEvslot, fParticle, fPDG, fGVcode, fEindex, fBindex,
      fCharge, fProcess, fNsteps, (int)fSpecies, status[int(fStatus)],
      fMass, fXpos, fYpos, fZpos, fXdir, fYdir, fZdir, fP, fE,
      fTime, fPstep, fStep, fSnext, fSafety, fNintLen, fIntLen, fBoundary);
  
#ifndef GEANT_NVCC
  fPath->Print();
  fNextpath->Print();
#endif
#else
  TString path;
  fPath->GetPath(path);
  TString nextpath;
  fNextpath->GetPath(nextpath);

  Geant::Print(msg, "== Track %d: evt=%d slt=%d part=%d pdg=%d gvc=%d eind=%d bind=%d chg=%d proc=%d nstp=%d "
         "spc=%d status=%s mass=%g xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g "
         "time=%g edep=%g pstp=%g stp=%g snxt=%g saf=%g nil=%g ile=%g bdr=%d\n pth=%s npth=%s\n",
         itr, fEvent, fEvslot, fParticle, fPDG, fGVcode, fEindex, fBindex,
         fCharge, fProcess, fNsteps, (int)fSpecies, status[int(fStatus)],
         fMass, fXpos, fYpos, fZpos, fXdir, fYdir, fZdir, fP, fE,
         fTime, fEdep, fPstep, fStep, fSnext, fSafety, fNintLen, fIntLen, fBoundary, path.Data(),
         nextpath.Data());
#endif
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
size_t GeantTrack::SizeOfInstance(size_t maxdepth) {
  // return the contiguous memory size needed to hold a GeantTrack
  // The VecGeom navigation state requires alignment, so we need to account for the
  // worst possible scenario:
  //
  // |--padding--|--padding--|--padding--|--padding--|--padding--| (adresses)
  // |--*start----*empty_spc.*fPath-------*empty_spc.*fNextpath--  (obj. layout)
  return ( sizeof(GeantTrack) + 
           2 * VolumePath_t::SizeOfInstance(maxdepth) + 
           2 * GEANT_ALIGN_PADDING );
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
GeantTrack *GeantTrack::MakeInstanceAt(void *addr, int maxdepth) {
  return new (addr) GeantTrack(addr, maxdepth);
}

} // GEANT_IMPL_NAMESPACE
} // Geant
