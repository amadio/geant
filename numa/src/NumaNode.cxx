#include "NumaNode.h"
#include "NumaUtils.h"

#ifdef USE_NUMA
#include <numa.h>

constexpr size_t MByte = 1024*1024;
constexpr size_t kByte = 1024;
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
NumaNode::NumaNode(int id, int maxcpus)
         :fId(id), fNcpus(0), fTotcpus(maxcpus), fNthreads(0),
          fMemTotal(0), fMemFree(0), fCpus(0), fMutex()
{
  // Constructor.
#ifdef USE_NUMA
  NumaUtils *utils = NumaUtils::Instance();
  hwloc_topology_t &topology = utils->fTopology;
  hwloc_obj_t obj = hwloc_get_numanode_obj_by_os_index(topology, id);
  fMemTotal = obj->memory.local_memory;

//  fMemTotal = numa_node_size(id, &fMemFree);
  fCpus = new int[maxcpus];  
  bitmask* bm = numa_bitmask_alloc(maxcpus);
  numa_node_to_cpus(id, bm);
  for(size_t i=0;i<bm->size;++i) {
    if(numa_bitmask_isbitset(bm, i))
      fCpus[fNcpus++] = i;
  }      
  numa_bitmask_free(bm);
#endif
}

//______________________________________________________________________________
int NumaNode::PinThread()
{
// Pin caller thread to this NUMA node. Pins actually to next free core.
// If all cores are given, restart from first core for this NUMA node.
  int core = 0;
#ifdef USE_NUMA
  std::lock_guard<std::mutex> lock(fMutex);
  core = fCpus[fNthreads%fNcpus];
  NumaUtils::Instance()->PinToCore(core);
  fNthreads++;
#endif
  return core;
}

std::ostream& operator<<(std::ostream& os, const NumaNode& node)
{
#ifdef USE_NUMA
  os << "Node id: " << node.fId << "   phys. cores: " << node.fNphysical << "   logical cores: " << node.fNcpus;
  os << " (";
  for (auto i=0; i<node.fNcpus; ++i)
    os << node.fCpus[i] << " ";
  os << ")\n";
  os << "  Total memory: " << node.fMemTotal/MByte << " MB";
#endif
  (void)node;
  return os;                            
}

} // GEANT_IMPL_NAMESPACE
} // Geant
