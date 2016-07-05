#include "NumaPolicy.h"

#ifdef USE_NUMA
#include <numa.h>
#include <numaif.h>
#endif

#include "NumaNode.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/* NUMA policy class  */
//______________________________________________________________________________
int NumaPolicy::AllocateNextThread()
{
// Pin the caller thread according to the NUMA policy. Return the logical core id.
// Returns NUMA node id
  fNthreads++;
#ifdef USE_NUMA
  auto crt_cpu = sched_getcpu();
  auto crt_node = numa_node_of_cpu(crt_cpu);
  if (fPolicy == kSysDefault)
    return (crt_node);
  int nodemin = 0;
  int minthr = fTopo.fNthreads[0];
  // Find node with smallest number of pinned threads
  int nnodes = fTopo.fNodes;
  for (int inode=0; inode<nnodes; ++inode) {
    if (fTopo.fNthreads[inode] < minthr) {
      minthr = fTopo.fNthreads[inode];
      nodemin = inode;
    }
  }  

  if (fPolicy & kCompact) {
    // Fill current NUMA node
    int nnodes = fTopo.fNodes;
    for (int inode=0; inode<nnodes; ++inode) {
      int npernode = fTopo.GetNode(inode)->fNphysical;
      if (fPolicy & kHTcompact) npernode = fTopo.GetNode(inode)->fNcpus;
      if (fTopo.fNthreads[inode] < npernode) {
        // Pin to this NUMA node
        fTopo.PinToNode(inode);
        return ( numa_node_of_cpu(sched_getcpu()) );
      }
    }
    // All NUMA nodes are full: allocate on the node having minimum nthreads
    fTopo.PinToNode(nodemin);
    return ( numa_node_of_cpu(sched_getcpu()) );
  }

  if (fPolicy & kScatter) {
    // Fill evenly NUMA nodes
    fTopo.PinToNode(nodemin);
    return ( numa_node_of_cpu(sched_getcpu()) );      
  }     
  return (crt_node);
#else
  return 0;
#endif
}

} // GEANT_IMPL_NAMESPACE
} // Geant
