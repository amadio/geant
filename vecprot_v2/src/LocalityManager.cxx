#include "LocalityManager.h"

#include <cassert>
#include <iostream>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

LocalityManager *LocalityManager::fgInstance = nullptr;

//______________________________________________________________________________
LocalityManager *LocalityManager::Instance() {
  if (fgInstance) return fgInstance;
  return (new LocalityManager());
}

//______________________________________________________________________________
LocalityManager::LocalityManager()
    :fInitialized(false),
     fPolicy(NumaPolicy::kSysDefault),
     fNnodes(0),
     fTrackMgr(nullptr),
     fNblocks(0),
     fBlockSize(0),
     fMaxdepth(0) {
  // Private constructor.
  fgInstance = this;
}

//______________________________________________________________________________
LocalityManager::~LocalityManager() {
  if (fTrackMgr) {
    for (auto i=0; i<fNnodes; ++i) delete fTrackMgr[i];
    delete [] fTrackMgr;
  }
  fgInstance = nullptr;
}

//______________________________________________________________________________
void LocalityManager::Init() {
  if (fInitialized) return;
  // If number of nodes is not specified, select according the policy
  assert(fNblocks > 0 && "Number of initial blocks not set.");
  assert(fBlockSize > 0 && "Block size for tracks not set.");
  assert(fMaxdepth > 0 && "Maximum geometry depth not set.");
  std::cout << *fPolicy.GetTopology() << std::endl;
  if (!fNnodes) fNnodes = fPolicy.GetTopology()->fNodes;
  if (!fNnodes) fNnodes = 1;
  fTrackMgr = new TrackManager*[fNnodes];
  for (int inode=0; inode<fNnodes; ++inode) {
    fTrackMgr[inode] = new TrackManager(fNblocks, fBlockSize, fMaxdepth, inode);
  }
  std::cout << "=== Locality manager allocated " << fNblocks << " blocks of " 
            << fBlockSize << " tracks each on " << fNnodes << " locality nodes ===" << std::endl;
  fInitialized = true;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
