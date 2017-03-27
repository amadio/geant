#include "Basket.h"
#include "Geant/Error.h"
#include "GeantTaskData.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Basket::Basket(int size, int threshold, int node)
      : fThreshold(threshold), fNode(node) {
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
void Basket::Print(const char *msg) const {
  // Print basket content.
  Printf("*** %s basket %p: ntr = %d   NUMA node = %d",
                msg, this, GetNtracks(), fNode);
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
  assert(itr < int(fTracks.size()));
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
VECCORE_ATT_HOST_DEVICE
void Basket::SetThreshold(int threshold) {
  // Set transport threshold for the basket
  fThreshold = threshold;
  fTracks.reserve(threshold);
}

} // GEANT_IMPL_NAMESPACE
} // Geant
