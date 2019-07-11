#include "Geant/Track.h"

#include "Geant/Error.h"
#include <execinfo.h>
#include "Geant/Propagator.h"
#include "Geant/TaskData.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
TrackDataMgr::TrackDataMgr(size_t maxdepth) : fMaxDepth(maxdepth)
{
// Private constructor. Make sure that the instance was initialized with a depth.
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  if (fMaxDepth == 0) std::runtime_error("Track data manager was not provided a geometry depth");
#endif
  // Compute the total track size in the assumption that there is no user data
  fTrackSize = sizeof(Track) + 2 * VolumePath_t::SizeOfInstance(maxdepth) + 2 * GEANT_ALIGN_PADDING;
  // Each registered user data of type T will top up round_up_align(sizeof(T)) bytes.
}

//______________________________________________________________________________
TrackDataMgr *TrackDataMgr::GetInstance(size_t maxdepth)
{
  // Use a static local variable to store the singleton. This avoids a global
  // data member while being thread safe. The instance is initialized after
  // invoking the user detector construction.
  static TrackDataMgr *fgInstance = new TrackDataMgr(maxdepth);
  return fgInstance;
}

/* Obtain a backtrace and print it to stdout. */
void printrace(void)
{
  void *array[10];
  size_t size;
  char **strings;
  size_t i;

  size    = backtrace(array, 10);
  strings = backtrace_symbols(array, size);

  printf("Obtained %zd stack frames.\n", size);

  for (i = 0; i < size; i++)
    printf("%s\n", strings[i]);

  free(strings);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Track::Track(void *addr)
{
  // In place private constructor. The provided address does NOT need to be aligned.
  size_t maxdepth = TrackDataMgr::GetInstance()->GetMaxDepth();
  // The start address of the extra data block has to be aligned
  fExtraData = round_up_align((char *)addr + sizeof(Track));
  // Initialize all user data calling the in-place constructors
  TrackDataMgr::GetInstance()->InitializeTrack(*this);
  // Geometry paths follow
  char *path_addr = fExtraData + TrackDataMgr::GetInstance()->GetDataSize();
  fPath           = VolumePath_t::MakeInstanceAt(maxdepth, path_addr);
  path_addr += round_up_align(VolumePath_t::SizeOfInstance(maxdepth));
  fNextpath = VolumePath_t::MakeInstanceAt(maxdepth, path_addr);

  // this will be changed/removed after we already have the previous step length stored in the track
  for (size_t i = 0; i < kNumPhysicsProcess; ++i) {
    fPhysicsNumOfInteractLengthLeft[i] = -1.0;
    fPhysicsInteractLength[i]          = 1.0;
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Track &Track::operator=(const Track &other)
{
  // Assignment
  if (&other != this) {
    fEvent              = other.fEvent;
    fEvslot             = other.fEvslot;
    fParticle           = other.fParticle;
    fPrimaryIndx        = other.fPrimaryIndx;
    fMother             = other.fMother;
    fGVcode             = other.fGVcode;
    fEindex             = other.fEindex;
    fCharge             = other.fCharge;
    fProcess            = other.fProcess;
    fNsteps             = other.fNsteps;
    fMaxDepth           = other.fMaxDepth;
    fStage              = other.fStage;
    fGeneration         = other.fGeneration;
    fSpecies            = other.fSpecies;
    fStatus             = other.fStatus;
    fMass               = other.fMass;
    fXpos               = other.fXpos;
    fYpos               = other.fYpos;
    fZpos               = other.fZpos;
    fXdir               = other.fXdir;
    fYdir               = other.fYdir;
    fZdir               = other.fZdir;
    fP                  = other.fP;
    fE                  = other.fE;
    fLogEkin            = other.fLogEkin;
    fLocalTime          = other.fLocalTime;
    fGlobalTime         = other.fGlobalTime;
    fVelocity           = other.fVelocity;
    fEdep               = other.fEdep;
    fPstep              = other.fPstep;
    fStep               = other.fStep;
    fSnext              = other.fSnext;
    fSafety             = other.fSafety;
    fNintLen            = other.fNintLen;
    fIntLen             = other.fIntLen;
    fBoundary           = other.fBoundary;
    fPending            = other.fPending;
    fStage              = other.fStage;
    fIsOnBoundaryPreStp = other.fIsOnBoundaryPreStp;
    fVolume             = other.fVolume;

    // Copy user data
    memcpy(fExtraData, other.fExtraData, TrackDataMgr::GetInstance()->GetDataSize());

    *fPath     = *other.fPath;
    *fNextpath = *other.fNextpath;

    fPhysicsProcessIndex = other.fPhysicsProcessIndex;
    for (size_t i = 0; i < kNumPhysicsProcess; ++i) {
      fPhysicsNumOfInteractLengthLeft[i] = other.fPhysicsNumOfInteractLengthLeft[i];
      fPhysicsInteractLength[i]          = other.fPhysicsInteractLength[i];
    }
  }
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Track::Clear(const char *)
{
  // Resets track content.
  fEvent       = -1;
  fEvslot      = -1;
  fParticle    = -1;
  fPrimaryIndx = -1;
  fMother      = 0;
  fGVcode      = 0;
  fEindex      = 0;
  fCharge      = 0;
  fProcess     = -1;
  fNsteps      = 0;
  fSpecies     = kHadron;
  fStatus      = kAlive;
  fMass        = 0.;
  fXpos        = 0.;
  fYpos        = 0.;
  fZpos        = 0.;
  fXdir        = 0.;
  fYdir        = 0.;
  fZdir        = 0.;
  fP           = 0.;
  fE           = 0.;
  fLogEkin     = -1;
  fLocalTime   = 0.;
  fGlobalTime  = 0.;
  fVelocity    = 0.;
  fEdep        = 0;
  fPstep       = 1.E20;
  fStep        = 0.;
  fSnext       = 0.;
  fSafety      = 0.;
  fNintLen     = 0.;
  fIntLen      = 0.;
  fBoundary    = false;
  fPending     = false;
  fMaxDepth    = 0;
  fStage       = 0;
  fGeneration  = 0;
  fVolume      = nullptr;
  fPath->Clear();
  fNextpath->Clear();

  // this will be changed/removed after we already have the previous step length stored in the track
  fPhysicsProcessIndex = -1;
  for (size_t i = 0; i < kNumPhysicsProcess; ++i) {
    fPhysicsNumOfInteractLengthLeft[i] = -1.0;
    fPhysicsInteractLength[i]          = 1.0;
  }

  // Clear user data
  TrackDataMgr::GetInstance()->InitializeTrack(*this);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Track::Reset(Track const &blueprint)
{
  // Fast reset of the content of an existing track, avoiding in-place construction

  // Copy main content from blueprint, except the pointers to geometry states and user data
  memcpy(&fEvent, &blueprint.fEvent, sizeof(Track) - 3 * sizeof(void *));

  // Copy user data
  memcpy(fExtraData, blueprint.fExtraData, TrackDataMgr::GetInstance()->GetDataSize());

  // Clear Geometry path
  fPath->Clear();
  fNextpath->Clear();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Track *Track::Clone(TaskData *td)
{
  Track &track = td->GetNewTrack();
  track        = *this;
  return &track;
}

//______________________________________________________________________________
Material_t *Track::GetMaterial() const
{
  // Current material the track is into
  return ((Material_t *)fVolume->GetMaterialPtr());
}

//______________________________________________________________________________p
VECCORE_ATT_HOST_DEVICE
void Track::SetPath(VolumePath_t const *const path)
{
  // Set path.
  *fPath = *path;
  UpdateVolume();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Track::SetNextPath(VolumePath_t const *const path)
{
  // Set next path.
  *fNextpath = *path;
}

//______________________________________________________________________________
void Track::Print(const char *msg) const
{
  const char *status[8] = {"alive", "killed", "inflight", "boundary", "exitSetup", "physics", "postponed", "new"};

  printf("%s: evt=%d slt=%d part=%d prim=%d mth=%d gvc=%d eind=%d bind=%d chg=%d proc=%d nstp=%d spc=%d "
         "status=%s mass=%g "
         "xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g time=%g pstp=%g stp=%g snxt=%g saf=%g nil=%g "
         "ile=%g bdr=%d\n",
         msg, fEvent, fEvslot, fParticle, fPrimaryIndx, fMother, fGVcode, fEindex, fBindex, fCharge, fProcess, fNsteps,
         (int)fSpecies, status[int(fStatus)], fMass, fXpos, fYpos, fZpos, fXdir, fYdir, fZdir, fP, fE, fLocalTime,
         fPstep, fStep, fSnext, fSafety, fNintLen, fIntLen, fBoundary);

  TrackDataMgr::GetInstance()->PrintUserData(*this);
#ifndef VECCORE_CUDA
  fPath->Print();
  fNextpath->Print();
#endif
}

//______________________________________________________________________________
void Track::PrintTracks(TrackVec_t &tracks)
{
  for (auto track : tracks)
    track->Print("xxx");
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
size_t Track::SizeOfInstance()
{
  // return the contiguous memory size needed to hold a Track
  return (TrackDataMgr::GetInstance()->GetTrackSize());
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Track *Track::MakeInstanceAt(void *addr)
{
  return new (addr) Track(addr);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Track *Track::MakeInstance()
{
  char *buffer    = new char[Track::SizeOfInstance()];
  Track *track    = new (buffer) Track(buffer);
  track->fOwnPath = true;
  return track;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Track::ReleaseInstance(Track *track)
{
  if (track->fOwnPath) delete[](char *) track;
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
