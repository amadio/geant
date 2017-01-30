#include "Basketizer.h"
#include <iostream>
#include <sys/time.h>
#include <thread>
#include <cmath>
#include <future>
#include <random>

#include "GeantNuma.h"

struct test_track {
  int id_;
  int numa_;
  int type_;
  double x_;
  double y_;
  double z_;
  double px_;
  double py_;
  double pz_;
  test_track() : id_(0), numa_(-1), type_(-1), x_(0.), y_(0.), z_(0.), px_(0.), py_(0.), pz_(0.) {}
};

//______________________________________________________________________________
void help() { std::cout << "Usage: testBasketizer <Nthreads>\n"; }

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
void Process(std::vector<test_track *> &basket) {
  // Emulate CPU time on a basket
  const int load = 35;
  for (size_t itr = 0; itr < basket.size(); ++itr) {
    double x = 1.;
    for (int j = 0; j < 100 * load; ++j) {
      x *= ((*basket[itr]).id_) % 100;
      double y = std::sqrt(x) * std::tanh(x);
      x = std::sin(y) / std::atan(x);
    }
    (*basket[itr]).x_ = x;
    (*basket[itr]).y_ = x;
    (*basket[itr]).z_ = x;
  }
}

using namespace Geant;
//______________________________________________________________________________
struct Workload {
  NumaPolicy policy_; // kSysDefault, kCompact, kScatter
  size_t nnodes_;
  size_t bsize_;
  size_t buf_size_;
  size_t nthreads_;
  std::promise<size_t> *promise_;
  std::future<size_t> *future_;
  size_t checksum_ref_;
  test_track *tracks_;
  size_t ntracks_;
  Basketizer<test_track> **basketizers_;
  void *allocated_;
  std::atomic_flag fLock;

  void Lock() {
    while (fLock.test_and_set(std::memory_order_acquire))
      ;
  }

  void Unlock() { fLock.clear(std::memory_order_release); }

  Workload(size_t nthreads)
      :
        policy_(Geant::NumaPolicy::kCompact),
        nnodes_(1), nthreads_(nthreads), fLock() {
    using Basketizer = Geant::Basketizer<test_track>;
    std::cout << *policy_.GetTopology() << std::endl;
    nnodes_ = policy_.fTopo.fNodes;
    if (nnodes_ < 1) nnodes_ = 1;
    //numa_set_localalloc();
    fLock.clear();
    promise_ = new std::promise<size_t>[nthreads];
    future_ = new std::future<size_t>[nthreads];
    basketizers_ = new Basketizer *[nnodes_];
    for (size_t n = 0; n < nnodes_; ++n)
      basketizers_[n] = nullptr;
    for (size_t n = 0; n < nthreads; ++n)
      future_[n] = promise_[n].get_future();
    ntracks_ = 1 << 20;
    const size_t nfilters = 1000;
    bsize_ = 16;
    buf_size_ = 1 << 14;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> rnd(0, nfilters - 1);
// Create a pool of numbers
    allocated_ = NumaAlignedMalloc(ntracks_ * sizeof(test_track), 0 /*numa_node*/, 64);
    tracks_ = new (allocated_) test_track[ntracks_];
    //Lock();
    //    std::cout << "Allocated data " << allocated_ << " on node: " << numa_node_addr(allocated_) << std::endl;
    //Unlock();
    for (size_t i = 0; i < ntracks_; ++i) {
      tracks_[i].id_ = i;
      // Generate randomly the filtering type
      tracks_[i].type_ = rnd(gen);
    }
    checksum_ref_ = ntracks_ * (tracks_[0].id_ + tracks_[ntracks_ - 1].id_) / 2;
  };

  ~Workload() {
    NumaAlignedFree(allocated_);
    for (size_t i = 0; i < nnodes_; ++i)
      NumaAlignedFree(basketizers_[i]);
    delete[] basketizers_;
  }

  void InitBasketizers(size_t node) {
    // Create basketizers for each numa node. To be call by workers.
    using Basketizer = Geant::Basketizer<test_track>;
    Lock();
    size_t basket_size = Basketizer::SizeofInstance(buf_size_);
    if (!basketizers_[node]) {
      basketizers_[node] = Basketizer::MakeInstanceAt(NumaAlignedMalloc(basket_size, node, 64), buf_size_, bsize_);
//      basketizers_[node] = new Basketizer(buf_size_, bsize_);
      std::cout << "basketizer[" << node << "] allocated on NUMA node " << NumaNodeAddr(basketizers_[node])
                << std::endl;
    }
    Unlock();
  }
  
  void CheckBasketizers() {
    for (size_t i=0; i<nnodes_; ++i) {
      std::cout << "Basketizer #" << i << ":\n";
      if (basketizers_[i]) basketizers_[i]->CheckBaskets();
    }
  }
};

//______________________________________________________________________________
void AddTracks(int tid, Workload *work, size_t nchunk, size_t ntracks) {
  // Code run by each thread
  // Pin thread according to the NUMA policy
  int numa_node = 0;
  numa_node = work->policy_.AllocateNextThread();
  //  numa_node = (numa_node+1)%2;
  work->InitBasketizers(numa_node);
//  work->Lock();
//  std::cout << "thread " << tid << " allocated on NUMA node: " << numa_node << std::endl;
//  work->Unlock();
  test_track *tracks = &work->tracks_[tid * nchunk];
  std::vector<test_track *> basket;
  basket.reserve(work->bsize_);
  size_t checksum = 0;
  //  size_t checksum_ref = ntracks * (tracks[0].id_ + tracks[ntracks-1].id_) / 2;
  //  int npending, npmax=0;
  for (size_t i = 0; i < ntracks; ++i) {
    //    int type = tracks[i].type_;
    basket.clear();
    if (work->basketizers_[numa_node]->AddElement(&tracks[i], basket)) {
      // We got a basket = add to partial checksum
      for (size_t itr = 0; itr < basket.size(); ++itr) {
        checksum += (*basket[itr]).id_;
        //        npending = basketizers[0].GetNpending();
        //        if (npending > npmax) npmax = npending;
      }
      Process(basket);
    }
  }
  basket.clear();

  if (work->basketizers_[numa_node]->Flush(basket)) {
    for (size_t itr = 0; itr < basket.size(); ++itr) {
      checksum += (*basket[itr]).id_;
      Process(basket);
    }
  }

  work->promise_[tid].set_value(checksum);
  //  std::cout << "thread " << tid << " npmax: " << npmax << std::endl;
  //  std::cout << "thread " << tid << " checksum: " << checksum << " ref: " <<
  //            checksum_ref << std::endl;
}

//______________________________________________________________________________
int main(int argc, char *argv[]) {
  using namespace Geant;

  if (argc == 1) {
    help();
    exit(0);
  }
  size_t nthreads = atoi(argv[1]);
  Workload work(nthreads);
  size_t nchunk = work.ntracks_ / nthreads;
  std::vector<std::thread> v;
  double cpu0 = get_cpu_time();
  double rt0 = get_wall_time();
  for (size_t n = 0; n < nthreads; ++n) {
    size_t ntoprocess = nchunk;
    if (n == nthreads - 1)
      ntoprocess = work.ntracks_ - n * nchunk;
    v.emplace_back(AddTracks, n, &work, nchunk, ntoprocess);
  }
  for (auto &t : v) {
    t.join();
  }
  double rt1 = get_wall_time();
  double cpu1 = get_cpu_time();
  size_t checksum = 0;
  for (size_t n = 0; n < nthreads; ++n)
    checksum += work.future_[n].get();

  std::cout << "run time: " << rt1 - rt0 << "   cpu time: " << cpu1 - cpu0 << "  checksum: " << checksum
            << " ref: " << work.checksum_ref_ << std::endl;
  if (  checksum > work.checksum_ref_) std::cout << "### Data overwritten! ###\n";
  else if (checksum < work.checksum_ref_) std::cout << "### Not all tracks collected! ###\n";
//  work.CheckBasketizers();
}
