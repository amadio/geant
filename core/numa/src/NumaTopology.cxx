#ifdef GEANT_USE_NUMA
#include <unistd.h>
#endif

#include "Geant/NumaUtils.h"
#include "Geant/NumaTopology.h"
#include "Geant/NumaNode.h"


namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
NumaTopology::NumaTopology()
{
// Default constructor
#ifdef GEANT_USE_NUMA
  hwloc_topology_t const &topology = NumaUtils::Topology();
  fAvailable = NumaUtils::NumaAvailable();
  if (fAvailable) {
    hwloc_const_bitmap_t cset;
    // retrieve the entire set of NUMA nodes and count them
    cset = hwloc_topology_get_topology_nodeset(topology);
    fNodes = hwloc_bitmap_weight(cset);
    fPageSize = sysconf(_SC_PAGESIZE);
    fListNodes = new NumaNode*[fNodes];
    for (auto i=0; i<fNodes; ++i) {
      fListNodes[i] = new NumaNode(i);
      fNcores += fListNodes[i]->fNcores;
      fNcpus += fListNodes[i]->fNcpus;
    }
    if (fNcores*fNcpus == 0) {
      std::cout << "*** Problem counting CPUs by depth\n";
      fNcores = 1;
      fNcpus = 2;
      fAvailable = false;
      std::cout << "*** NUMA disabled.\n";
      return;
    }
    fHT = fNcpus/fNcores;
    
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

//______________________________________________________________________________
int NumaTopology::NumaNodeOfCpu(int cpu) const
{
// Find the numa node to which the cpu belongs.
  for (int i=0; i<fNodes; ++i)
    if (fListNodes[i]->HasCpu(cpu)) return i;
  return -1;
}

//______________________________________________________________________________
int NumaTopology::BindToNode(int node)
{
// Pin thread to NUMA node
#ifdef GEANT_USE_NUMA
  if (node  < 0 || node > fNodes-1) {
    std::cout << "Error: trying to bind to NUMA node " << node << " outside range\n";
    return -1;
  }
  int cpu = fListNodes[node]->BindThread();
  if (cpu < 0) {
    std::cout << "Error: could not bind to NUMA node " << node << "\n";
    return -1;
  }
  return cpu;
#else
  (void)node;
  return -1;
#endif  
}  

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& os, const NumaTopology& 
#ifdef GEANT_USE_NUMA
                           topo
#endif
                         )
{
#ifdef GEANT_USE_NUMA
  os << "NUMA nodes: " << topo.fNodes << "   cores: " << topo.fNcores
     << "   cpus: " << topo.fNcpus 
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

