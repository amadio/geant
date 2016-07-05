#include "NumaTopology.h"

#ifdef USE_NUMA
#include <numa.h>
#include <hwloc.h>
#endif

#include "NumaNode.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
NumaTopology::NumaTopology()
             :fAvailable(false), fNodes(0), fNcpus(0), fNphysical(0),
              fPageSize(0), fListNodes(0), fNthreads(0)
{
// Default constructor
#ifdef USE_NUMA
  fAvailable = (numa_available() < 0) ? false : true;
  if (fAvailable) {
    fNodes = numa_num_configured_nodes();
    fNcpus = numa_num_task_cpus();
    fPageSize = numa_pagesize();
    // Get CPU configuration per numa node
    fNphysical = FindPhysCores(fHT);
    fListNodes = new NumaNode*[fNodes];
    fNthreads = new int[fNodes];
    for (auto i=0; i<fNodes; ++i) {
      fListNodes[i] = new NumaNode(i, fNcpus);
      fListNodes[i]->fNphysical = fNphysical/fHT;
      fNthreads[i] = 0;
    }
  }
#endif  
}

NumaTopology::~NumaTopology()
{
// Destructor
  for (auto i=0; i<fNodes; ++i) delete fListNodes[i];
  delete [] fListNodes;
}

int NumaTopology::FindPhysCores(int &ht) const
{
  // Allocate, initialize, and perform topology detection
  int cores = 0;
  ht = 0;
#ifdef USE_NUMA
  hwloc_topology_t topology;
  hwloc_topology_init(&topology);
  hwloc_topology_load(topology);

  // Try to get the number of CPU cores from topology
  int depth = hwloc_get_type_depth(topology, HWLOC_OBJ_CORE);
  if(depth == HWLOC_TYPE_DEPTH_UNKNOWN)
    std::cout << "*** The number of cores is unknown\n";
  else {
    cores = hwloc_get_nbobjs_by_depth(topology, depth);
    ht = hwloc_get_nbobjs_by_depth(topology, depth+1);
    ht /= cores;
    if (ht == 0) ht = 1;
  }  

  // Destroy topology object and return
  hwloc_topology_destroy(topology);
#endif
  return cores;
}

int NumaTopology::PinToNode(int 
#ifdef USE_NUMA
                              node
#endif
                           )
{
// Pin thread to NUMA node
#ifdef USE_NUMA
  if (node > fNodes-1) {
    std::cout << "Error: trying to pin to a NUMA node beyond max range\n";
    return -1;
  }
  fNthreads[node]++;
  return fListNodes[node]->PinThread();
#else
  return -1;
#endif  
}  

std::ostream& operator<<(std::ostream& os, const NumaTopology& 
#ifdef USE_NUMA
                           topo
#endif
                         )
{
#ifdef USE_NUMA
  os << "NUMA nodes: " << topo.fNodes << "   cpu's: " << topo.fNcpus
     << "   physical cpu's: " << topo.fNphysical 
     << "   hyperthreading: " << topo.fHT << std::endl;
  for (auto i=0; i<topo.fNodes; ++i) os << *(topo.fListNodes[i]) << std::endl;
#else
  os << "NUMA not available.\n";
#endif
  return os;
}      

} // GEANT_IMPL_NAMESPACE
} // Geant

