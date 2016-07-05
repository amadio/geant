#include "LocalityManager.h"

#include <cassert>

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
     fBlockSize(0) {
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
  if (!fNnodes) fNnodes = fPolicy.GetTopology()->fNodes;
  fTrackMgr = new TrackManager*[fNnodes];
  for (int inode=0; inode<fNnodes; ++inode) {
    fTrackMgr[inode] = new TrackManager(fNblocks, fBlockSize, inode);
  }
  
}

} // GEANT_IMPL_NAMESPACE
} // Geant
