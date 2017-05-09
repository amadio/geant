#include "NumaCore.h"
#include "NumaUtils.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
NumaCore::NumaCore(int id, int node)
         :fId(id), fNode(node)
{
  // Constructor.
#if defined(USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
  hwloc_topology_t const &topology = NumaUtils::Topology();
  hwloc_obj_t onode = hwloc_get_numanode_obj_by_os_index(topology, node);
  // Object for the core with index 'id'
  fObjCore = hwloc_get_obj_inside_cpuset_by_type(topology, onode->cpuset, HWLOC_OBJ_CORE, id);
  // Number of logical PUs in the core
  fNcpus = hwloc_get_nbobjs_inside_cpuset_by_type(topology, fObjCore->cpuset, HWLOC_OBJ_PU);
  fCpus = new int[fNcpus];
  // now loop over all PU objects
  for (int ipu = 0; ipu < fNcpus;  ++ipu) {
    hwloc_obj_t objpu = hwloc_get_obj_inside_cpuset_by_type(topology, fObjCore->cpuset, HWLOC_OBJ_PU, ipu);
    fCpus[ipu] = objpu->os_index;
  }
#endif
}

//______________________________________________________________________________
bool NumaCore::HasCpu(int cpu) const
{
  for (int i=0; i<fNcpus; ++i)
    if (fCpus[i] == cpu) return true;
  return false;
}

//______________________________________________________________________________
int NumaCore::BindThread()
{
// Pin caller thread to the next cpu on this core.
// If all cores are given, restart from first cpu.
  int cpu = -1;
#if defined(USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
  std::lock_guard<std::mutex> lock(fMutex);
  int cpuindex = fNthreads % fNcpus;
  cpu = fCpus[cpuindex];
  // Now bind the thread to the selected cpu
  hwloc_topology_t const &topology = NumaUtils::Topology();
  // Check current binding
  int binding = NumaUtils::GetCpuBinding();
//  std::cout << ">>> try to bind thread to cpu# " << cpu << "   currently bound to cpu# " << binding << std::endl;
  hwloc_obj_t objpu = hwloc_get_obj_inside_cpuset_by_type(topology, fObjCore->cpuset, HWLOC_OBJ_PU, cpuindex);
  assert((int)objpu->os_index == cpu);
  assert(hwloc_bitmap_weight(objpu->cpuset) == 1);
  hwloc_set_cpubind(topology, objpu->cpuset, HWLOC_CPUBIND_THREAD);
  binding = NumaUtils::GetCpuBinding();
  if (binding != cpu)
    std::cout << "### Cannot bind thread to cpu# " << cpu << std::endl;
//  std::cout << "    thread now bound to cpu# " << binding << std::endl;
  fNthreads++;
#endif
  return cpu;
}

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& os, const NumaCore& core)
{
#if defined(USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
  os << "core#" << core.fId << ": ";
  for (auto i=0; i<core.fNcpus; ++i)
    os << core.fCpus[i] << " ";
  os << "\n";
#endif
  (void)core;
  return os;                            
}

} // GEANT_IMPL_NAMESPACE
} // Geant
