#include "NumaNode.h"
#include "NumaCore.h"
#include "NumaUtils.h"

#ifdef USE_NUMA
constexpr size_t MByte = 1024*1024;
constexpr size_t kByte = 1024;
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
NumaNode::NumaNode(int id)
         :fId(id)
{
  // Constructor.
#ifdef USE_NUMA
  hwloc_topology_t const &topology = NumaUtils::Topology();
  hwloc_obj_t onode = hwloc_get_numanode_obj_by_os_index(topology, id);
  fMemTotal = onode->memory.local_memory;

  // Loop all physical cores
  fNcores = hwloc_get_nbobjs_inside_cpuset_by_type(topology, onode->cpuset, HWLOC_OBJ_CORE);
  fCores = new NumaCore*[fNcores];
  for (int icore = 0; icore < fNcores; ++icore) {
    fCores[icore] = new NumaCore(icore, id);
    fNcpus += fCores[icore]->fNcpus;
  }
#endif
}

//______________________________________________________________________________
NumaNode::~NumaNode()
{
  for (int i=0; i<fNcores; ++i)
    delete fCores[i];
  delete [] fCores;
}

//______________________________________________________________________________
bool NumaNode::HasCpu(int cpu) const
{
  for (int i=0; i<fNcores; ++i)
    if (fCores[i]->HasCpu(cpu)) return true;
  return false;
}

//______________________________________________________________________________
int NumaNode::BindThread()
{
// Pin caller thread to this NUMA node. Pins actually to next free core.
// If all cores are given, restart from first core for this NUMA node.
  int cpu = 0;
#ifdef USE_NUMA
  std::lock_guard<std::mutex> lock(fMutex);
  NumaCore *core = fCores[fNthreads%fNcores];
  cpu = core->BindThread();
  fNthreads++;
#endif
  return cpu;
}

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& os, const NumaNode& node)
{
#ifdef USE_NUMA
  os << "NUMA node: " << node.fId << " cores: " << node.fNcores
     << " cpus: " << node.fNcpus << " memory: " << node.fMemTotal/MByte << "MB\n";
  for (auto i=0; i<node.fNcores; ++i)
    os << "    " << *node.fCores[i];
#endif
  (void)node;
  return os;                            
}

} // GEANT_IMPL_NAMESPACE
} // Geant
