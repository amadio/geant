#include "NumaPolicy.h"
#include "NumaUtils.h"
#include "NumaNode.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

/* NUMA policy class  */
//______________________________________________________________________________
int NumaPolicy::NextNumaNode()
{
  // Find node with smallest number of pinned threads
  int nnodes = fTopo.fNodes;
  int nodemin = -1;
  int minthr = 1e8;
  for (int inode=0; inode<nnodes; ++inode) {
    if (fTopo.GetNode(inode)->fNthreads < minthr) {
      minthr = fTopo.GetNode(inode)->fNthreads;
      nodemin = inode;
    }
  }
  return nodemin;
}

//______________________________________________________________________________
int NumaPolicy::AllocateNextThread(int node)
{
// Pin the caller thread according to the NUMA policy. Return the logical core id.
// Returns NUMA node id
  fNthreads++;
#ifdef GEANT_USE_NUMA
  auto crt_cpu = NumaUtils::GetCpuBinding();
  auto crt_node = fTopo.NumaNodeOfCpu(crt_cpu);
  if (fPolicy == kSysDefault)
    return crt_node;
  if (node >=0) {
    fTopo.BindToNode(node);
    return node;
  }
  
  // Find next NUMA node
  int nnodes = fTopo.fNodes;
  int nodemin = NextNumaNode();
  
  if (fPolicy & kCompact) {
    // Fill current NUMA node
    for (int inode=0; inode<nnodes; ++inode) {
      int npernode = fTopo.GetNode(inode)->fNcores;
      if (fPolicy & kHTcompact) npernode = fTopo.GetNode(inode)->fNcpus;
      if (fTopo.GetNode(inode)->fNthreads < npernode) {
        // Pin to this NUMA node
        fTopo.BindToNode(inode);
        return inode;
      }
    }
    // All NUMA nodes are full: allocate on the node having minimum nthreads
    fTopo.BindToNode(nodemin);
    return nodemin;
  }

  if (fPolicy & kScatter) {
    // Fill evenly NUMA nodes
    fTopo.BindToNode(nodemin);
    return nodemin;
  }     
  return crt_node;
#else
  (void)node;
  return 0;
#endif
}

//______________________________________________________________________________
int NumaPolicy::MembindNode(int node)
{
  // Use hwloc to set the default binding policy of current thread to prefer
  // the NUMA node specified.
  if (node < 0) return -2;
#ifdef GEANT_USE_NUMA
  hwloc_topology_t const &topology = NumaUtils::Topology();
  hwloc_nodeset_t nodeset = hwloc_bitmap_alloc();
  hwloc_bitmap_only(nodeset, unsigned(node));
  //assert( hwloc_bitmap_isset(nodeset, node) == 0 );

  int result = hwloc_set_membind_nodeset(topology,
                                         nodeset,
                                         HWLOC_MEMBIND_BIND,
                                         HWLOC_MEMBIND_THREAD);
  hwloc_bitmap_free(nodeset);
  return result;
#endif
  // NUMA not enabled
  return -2;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
