#ifdef USE_NUMA
#include <unistd.h>
#endif

#include "NumaUtils.h"
#include "NumaTopology.h"
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
  NumaUtils *utils = NumaUtils::Instance();
  hwloc_topology_t &topology = utils->fTopology;
  fAvailable = utils->fAvailable;
  if (fAvailable) {
    hwloc_const_bitmap_t cset;
    // retrieve the entire set of NUMA nodes and count them
    cset = hwloc_topology_get_topology_nodeset(topology);
    fNodes = hwloc_bitmap_weight(cset);
    int err = FindCores();
    if (err < 0) {
      std::cout << "*** NUMA disabled.\n";
      fAvailable = false;
      return;
    }
    fPageSize = sysconf(_SC_PAGESIZE);
    fListNodes = new NumaNode*[fNodes];
    fNthreads = new int[fNodes];
    for (auto i=0; i<fNodes; ++i) {
      fListNodes[i] = new NumaNode(i, fNcpus);
      fListNodes[i]->fNphysical = fNphysical/fHT;
      fNthreads[i] = 0;
    }
    std::cout << "*** NUMA enabled.\n";
  }
#endif  
}

NumaTopology::~NumaTopology()
{
// Destructor
  for (auto i=0; i<fNodes; ++i) delete fListNodes[i];
  delete [] fListNodes;
}

int NumaTopology::FindCores()
{
  // Allocate, initialize, and perform topology detection
  fNphysical = 1;
  fNcpus = 2;
  fHT = 2;
#ifdef USE_NUMA
  // Try to get the number of CPU cores from topology
  NumaUtils *utils = NumaUtils::Instance();
  hwloc_topology_t &topology = utils->fTopology;
  int depth_core = hwloc_get_type_depth(topology, HWLOC_OBJ_CORE);
  int depth_pu = hwloc_get_type_depth(topology, HWLOC_OBJ_PU);
  if((depth_core == HWLOC_TYPE_DEPTH_UNKNOWN) || (depth_pu == HWLOC_TYPE_DEPTH_UNKNOWN)) {
    std::cout << "*** The number of core/cpu is unknown - problems with hwloc detection\n";
    return -1;
  } else {
    fNphysical = hwloc_get_nbobjs_by_depth(topology, depth_core);
    fNcpus = hwloc_get_nbobjs_by_depth(topology, depth_pu);
    if (fNphysical*fNcpus == 0) {
      std::cout << "*** Problem counting CPUs by depth\n";
      fNphysical = 1;
      fNcpus = 2;
      return -1;
    }
    fHT = fNcpus / fNphysical;
    return 0;
  }
#endif
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
     << "   hyperthreading: " << topo.fHT
     << "   page size: " << topo.fPageSize << std::endl;
  for (auto i=0; i<topo.fNodes; ++i) os << *(topo.fListNodes[i]) << std::endl;
#else
  os << "NUMA not available.\n";
#endif
  return os;
}      

} // GEANT_IMPL_NAMESPACE
} // Geant

