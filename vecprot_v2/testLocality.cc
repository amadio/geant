#include "LocalityManager.h"
#include <iostream>
#include <random>
#include <sys/time.h>
#include <thread>
#include <VecCore/VecMath.h>
#include "GeantNuma.h"
#include "NumaAllocator.h"

#include "Basketizer.h"
#include "GeantTrack.h"

#include "TGeoManager.h"
#include "TGeoBBox.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "management/RootGeoManager.h"
#include "management/GeoManager.h"
#include "navigation/VNavigator.h"
#include "navigation/SimpleNavigator.h"
#include "navigation/NewSimpleNavigator.h"
#include "navigation/SimpleABBoxNavigator.h"
#include "navigation/SimpleABBoxLevelLocator.h"
#include "navigation/HybridNavigator2.h"
#endif

//______________________________________________________________________________
void help() { std::cout << "Usage: testLocality [nthreads] [geom]>\n"; }

//______________________________________________________________________________
double get_wall_time() {
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

//______________________________________________________________________________
double get_cpu_time() { return (double)clock() / CLOCKS_PER_SEC; }

//______________________________________________________________________________
struct LocalData {
  using Basketizer = Geant::Basketizer<Geant::cxx::GeantTrack*>;
  size_t bsize;
  size_t buffer_size;
  std::atomic_flag fLock;
  Basketizer *basketizer;

  LocalData() : bsize(0), buffer_size(0), fLock(), basketizer(0) { fLock.clear(); }
  void Lock() {
    while (fLock.test_and_set(std::memory_order_acquire))
      ;
  }

  void Unlock() { fLock.clear(std::memory_order_release); }
};

//______________________________________________________________________________
inline void InitTrack(Geant::cxx::GeantTrack &track, double dx, double dy, double dz) {
// Initialize track position within the range (-dx,dx), (-dy,dy), (-dz,dz)
// Compute initial path and initialize direction randomly.
  using namespace Geant;
  using namespace vecgeom;
  using namespace vecCore::math;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> disx(-dx, dx), disy(-dy, dy), disz(-dz, dz);
  std::uniform_real_distribution<> disphi(0., 2.*M_PI), disrnd(0., 1.);
  track.fXpos = disx(gen);
  track.fYpos = disy(gen);
  track.fZpos = disz(gen);
  double phi = disphi(gen);
  double theta = ACos(1. - 2*disrnd(gen));
  track.fXdir = Sin(theta) * Cos(phi);
  track.fYdir = Sin(theta) * Sin(phi);
  track.fZdir = Cos(theta);
#ifdef USE_VECGEOM_NAVIGATOR
  SimpleNavigator nav;
  nav.LocatePoint(GeoManager::Instance().GetWorld(),
                    Vector3D<Precision>(track.fXpos, track.fYpos, track.fZpos), *track.Path(), true);
#else
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  nav->FindNode(track.fXpos, track.fYpos, track.fZpos);
#endif
}

//______________________________________________________________________________
void ConsumeTracks(size_t tid, LocalData *ldata) {
  using namespace Geant;
  LocalityManager *loc_mgr = LocalityManager::Instance();
  int node = loc_mgr->GetPolicy().AllocateNextThread();
  ldata->Lock();
  std::cout << "thread " << tid << " on node " << node << std::endl;
  ldata->Unlock();
  
  using track_allocator = NumaAllocator<GeantTrack*>;
  // The basket collecting the tracks
  std::vector<GeantTrack*, track_allocator> basket(ldata[node].bsize, nullptr, track_allocator(node));

  TrackManager &trk_mgr = loc_mgr->GetTrackManager(node);

  TGeoBBox* box = (TGeoBBox*)gGeoManager->GetTopVolume()->GetShape();
  double dx = box->GetDX();
  double dy = box->GetDY();
  double dz = box->GetDZ();
#ifndef USE_VECGEOM_NAVIGATOR
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif
  for (int itr=0; itr<1000000; ++itr) {
    GeantTrack &track = trk_mgr.GetTrack();
    InitTrack(track, dx, dy, dz);
    // Basketize the track
    loc_mgr->ReleaseTrack(track);
  }
}

//______________________________________________________________________________
int main(int argc, char *argv[]) {
  using namespace Geant;
  using namespace vecgeom;

  if (argc < 3) {
    help();
    exit(0);
  }
  size_t nthreads = atoi(argv[1]);
  std::string file = argv[2];  
  
  // Load geometry and convert to VecGeom
  TGeoManager::Import(file.c_str());
  if (!gGeoManager) {
    std::cout << "### Cannot load geometry from file " << file << std::endl;
    return 1;
  }
  int maxdepth = TGeoManager::GetMaxLevels();
  TrackDataMgr::GetInstance(maxdepth);
  if (nthreads > 1) gGeoManager->SetMaxThreads(nthreads);

#ifdef USE_VECGEOM_NAVIGATOR
  RootGeoManager::Instance().LoadRootGeometry();
  maxdepth = GeoManager::Instance().getMaxDepth();
  for (auto &lvol : GeoManager::Instance().GetLogicalVolumesMap()) {
    if (lvol.second->GetDaughtersp()->size() < 4) {
      lvol.second->SetNavigator(NewSimpleNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 5) {
      lvol.second->SetNavigator(SimpleABBoxNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 10) {
      lvol.second->SetNavigator(HybridNavigator<>::Instance());
      HybridManager2::Instance().InitStructure((lvol.second));
    }
    lvol.second->SetLevelLocator(SimpleABBoxLevelLocator::GetInstance());
  }
#endif
  
  // Configure the locality manager
  LocalityManager *mgr = LocalityManager::Instance();
  mgr->SetPolicy(NumaPolicy::kCompact);
  mgr->SetNblocks(100);
  mgr->SetBlockSize(10000);
  mgr->Init();
  int nnodes = mgr->GetNnodes();
  
  using Basketizer = Geant::Basketizer<GeantTrack*>;
  LocalData *ldata = new LocalData[nnodes];

  // Use NumaAllocator for the thread basket
  for (auto node=0; node<nnodes; node++) {
    ldata[node].bsize = 16;
    ldata[node].buffer_size = 1<<16;
    size_t needed = Basketizer::SizeofInstance(ldata[node].buffer_size);
    ldata[node].basketizer = Basketizer::MakeInstanceAt(NumaUtils::NumaAlignedMalloc(needed, node, 64), ldata[node].buffer_size, ldata[node].bsize);
  }

  std::vector<std::thread> v;
  double cpu0 = get_cpu_time();
  double rt0 = get_wall_time();
  for (size_t n = 0; n < nthreads; ++n) {
    v.emplace_back(ConsumeTracks, n, ldata);
  }
  for (auto &t : v) {
    t.join();
  }
  double rt1 = get_wall_time();
  double cpu1 = get_cpu_time();
  delete mgr;
//  delete gGeoManager;
  for (auto node=0; node<nnodes; node++) {
    NumaUtils::NumaAlignedFree(ldata[node].basketizer);
  }
  delete [] ldata;
  delete gGeoManager;
  std::cout << "run time/thread: " << (rt1 - rt0)/nthreads << "   cpu time/thread: " << (cpu1 - cpu0)/nthreads << std::endl;
  return 0;
}
