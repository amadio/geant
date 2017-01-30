#include "Basket.h"
#include "Geant/Error.h"
#include "GeantTaskData.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Basket::Basket(int size, ELocality locality, int threshold, int node)
      : fLocality(locality), fThreshold(threshold), fNode(node) {
  assert(size > 0 && threshold <= size);
  fTracks.reserve(size);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool Basket::Contains(int evstart, int nevents) const {
  // Checks if any of the tracks in the input array belongs to the given event
  // range.
  for (auto track : fTracks) {
    if (track->fEvent >= evstart && track->fEvent < evstart + nevents)
      return true;
  }
  return false;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
const char *Basket::GetLocalityString() const
{
  switch (fLocality) {
    case kNone:
      return ( "mixed" );
    case kVolume:
      return ( "volume" );
    case kParticle:
      return ( "particle" );
    case kProcess:
      return ( "process" );
  }
  return ( "unknown" );
}

//______________________________________________________________________________
void Basket::Print(const char *msg) const {
  // Print basket content.
  Printf("*** %s basket %p: ntr = %d  locality = %s  NUMA node = %d",
                msg, this, GetLocalityString(), GetNtracks(), fNode);
#ifndef VECCORE_CUDA
  std::string sno;
  for (unsigned int itr=0; itr<fTracks.size(); ++itr) {
    sno = std::to_string(itr);
    fTracks[itr]->Print(sno.c_str());
  }
#endif
}

//______________________________________________________________________________
void Basket::PrintTrack(int itr) const {
  // Print a given track.
#ifndef VECCORE_CUDA
  assert(itr < fTracks.size());
  std::string sno = std::to_string(itr);
  fTracks[itr]->Print(sno.c_str());
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void Basket::Recycle(GeantTaskData */*td*/) {
  // Recycle the basket to the volume scheduler.
//  td->RecycleBasket(this);
}

//______________________________________________________________________________
void Basket::SetThreshold(int threshold) {
  // Set transport threshold for the basket
  fThreshold = threshold;
  fTracks.reserve(threshold);
}

} // GEANT_IMPL_NAMESPACE
} // Geant
