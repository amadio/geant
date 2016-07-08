#include "LocalityManager.h"
#include <iostream>
#include <sys/time.h>
#include <thread>

//______________________________________________________________________________
void help() { std::cout << "Usage: testLocality <Nthreads>\n"; }

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
void ConsumeTracks(size_t) {
  using namespace Geant;
  LocalityManager *loc_mgr = LocalityManager::Instance();
  int node = loc_mgr->GetPolicy().AllocateNextThread();
  TrackManager &trk_mgr = loc_mgr->GetTrackManager(node);
  for (int itr=0; itr<1000000; ++itr) {
    GeantTrack &track = trk_mgr.GetTrack();
//    printf("%p\n", &track);
//    track.Print("");
    loc_mgr->ReleaseTrack(track);
  }
}

//______________________________________________________________________________
int main(int argc, char *argv[]) {
  using namespace Geant;

  if (argc == 1) {
    help();
    exit(0);
  }
  size_t nthreads = atoi(argv[1]);
  
  LocalityManager *mgr = LocalityManager::Instance();
  mgr->SetPolicy(NumaPolicy::kCompact);
  mgr->SetNblocks(100);
  mgr->SetBlockSize(10000);
  mgr->SetMaxDepth(10);
  mgr->Init();
  std::vector<std::thread> v;
  double cpu0 = get_cpu_time();
  double rt0 = get_wall_time();
  for (size_t n = 0; n < nthreads; ++n) {
    v.emplace_back(ConsumeTracks, n);
  }
  for (auto &t : v) {
    t.join();
  }
  double rt1 = get_wall_time();
  double cpu1 = get_cpu_time();
  delete mgr;

  std::cout << "run time/thread: " << (rt1 - rt0)/nthreads << "   cpu time/thread: " << (cpu1 - cpu0)/nthreads << std::endl;
}
