//===--- NumaTopology.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaTopology.h
 * @brief Class describing the NUMA topology of the system, allowing to set NUMA policies
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_TOPOLOGY
#define GEANT_NUMA_TOPOLOGY

#include <numa.h>
#include <numaif.h>
#include <hwloc.h>
#include <thread>
#include <mutex>
#include <vector>
#include <iostream>
#include <cmath>
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/syscall.h>
#include <pthread.h>
#include <sched.h>
#include <random>

namespace Geant {

constexpr size_t MByte = 1024*1024;
constexpr size_t kByte = 1024;

/** Pin a tread to a core */
//______________________________________________________________________________
void pin_to_core(size_t core)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core, &cpuset);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
}

/** Class describung NUMA node properties */
class NumaNode {
public:
  int        fId;            /* NUMA node id */
  int        fNcpus;         /* Number of CPU's allowed for the node */
  int        fTotcpus;       /* Total number of CPU's */
  int        fNphysical;     /* Number of physical cores */
  int        fNthreads;      /* Number of assigned threads for the node */ 
  long       fMemTotal;      /* Total memory */
  long       fMemFree;       /* Free memoru */
  int       *fCpus;          /* List of CPU's for the node */
  std::mutex fMutex;         /* Mutex for the node */
  
public:
  NumaNode(int id, int maxcpus) : fId(id), fNcpus(0), fTotcpus(maxcpus), fNthreads(0), 
               fMemTotal(0), fMemFree(0), fCpus(0), fMutex()
  {
    fMemTotal = numa_node_size(id, &fMemFree);
    fCpus = new int[maxcpus];
    bitmask* bm = numa_bitmask_alloc(maxcpus);
    numa_node_to_cpus(id, bm);
    for(size_t i=0;i<bm->size;++i) {
      if(numa_bitmask_isbitset(bm, i))
        fCpus[fNcpus++] = i;
    }      
    numa_bitmask_free(bm);
  }
  
  ~NumaNode() { delete [] fCpus; }
  
  int PinThread() {
  /* Pin caller thread to this NUMA node. Pins actually to next free core.
     If all cores are given, restart from first core for this NUMA node. */
    std::lock_guard<std::mutex> lock(fMutex);
    int core = fCpus[fNthreads%fNcpus];
    pin_to_core(core);
    fNthreads++;
    return core;
  }
};

/* NUMA topology class storing NUMA-specific data detected for the system */
//______________________________________________________________________________
class NumaTopology {
public:
  bool       fAvailable;     /* Is NUMA available on the system */
  int        fNodes;         /* Number of allowed NUMA nodes */
  int        fNcpus;         /* Number of allowed CPU's */
  int        fNphysical;     /* Number of physical cores */
  int        fHT;            /* Number of threads per core */
  size_t     fPageSize;      /* NUMA page size on system */
  NumaNode **fListNodes;     /* List of nodes */
  int       *fNthreads;      /* Number of threads currently allocated per NUMA node */
  
public:
  NumaTopology() : fAvailable(false), fNodes(0), fNcpus(0), fNphysical(0), 
                   fPageSize(0), fListNodes(0), fNthreads(0)
  {
    fAvailable = (numa_available() < 0) ? false : true;
    if (!fAvailable) return;
    fNodes = numa_num_configured_nodes();
    if (fNodes != numa_max_node()+1) {
      std::cout << "ERROR: numa_num_max_node (" << numa_max_node()
        << ") not matching numa_num_configured_nodes ( " << fNodes << ")\n";
      return;
    }
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
  
  ~NumaTopology() 
  {
    for (auto i=0; i<fNodes; ++i) delete fListNodes[i];
    delete [] fListNodes;
  }
  
  int FindPhysCores(int &ht) const
  {
    // Allocate, initialize, and perform topology detection
    hwloc_topology_t topology;
    hwloc_topology_init(&topology);
    hwloc_topology_load(topology);

    // Try to get the number of CPU cores from topology
    int cores = 0;
    ht = 0;
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
    return cores;
  }
  
  NumaNode *GetNode(int node) { return fListNodes[node]; }
  
  int PinToNode(int node)
  {
    if (node > fNodes-1) {
      std::cout << "Error: trying to pin to a NUMA node beyond max range\n";
      return -1;
    }
    fNthreads[node]++;
    return fListNodes[node]->PinThread();
  }  
};

/* NUMA policy class  */
//______________________________________________________________________________
class NumaPolicy {
public:
  enum EType {
    kSysDefault = 0,       /** System default scheduling policy */
    kCompact    = 1 << 0,  /** Compact allocation until filling every NUMA node */
    kScatter    = 1 << 1,  /** Scatter threads evenly across the system */
    kHTcompact  = 1 << 2   /** Use HT in compact mode before pinning to next NUMA node */
  };

  int          fNthreads;        /* Number of threads to be pinned */
  NumaTopology fTopo;            /* Numa topology of the machine */
  EType        fPolicy;          /* NUMA policy */
  
  NumaPolicy(EType policy) : fNthreads(0), fTopo(), fPolicy(policy) {
  // Default constructor
  }  
  
  NumaTopology *GetTopology() { return &fTopo; }
  void SetPolicy(unsigned int policy);

//______________________________________________________________________________
  int AllocateNextThread() {
  // Pin the caller thread according to the NUMA policy. Return the logical core id.
  // Returns NUMA node id
    fNthreads++;
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
  }
};
  
//______________________________________________________________________________
std::ostream& operator<<(std::ostream& os, const bitmask& bm)
{
    for(size_t i=0;i<bm.size;++i)
    {
        os << numa_bitmask_isbitset(&bm, i);
    }
    return os;
}

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& os, const NumaNode& node)
{
  os << "Node id: " << node.fId << "   phys. cores: " << node.fNphysical << "   logical cores: " << node.fNcpus;
  os << " (";
  for (auto i=0; i<node.fNcpus; ++i)
    os << node.fCpus[i] << " ";
  os << ")\n";
  os << "  Total memory: " << node.fMemTotal/MByte << " MB  free: " << 
                              node.fMemFree/MByte << " MB";
  return os;                            
}

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& os, const NumaTopology& topo)
{
    os << "NUMA nodes: " << topo.fNodes << "   cpu's: " << topo.fNcpus
       << "   physical cpu's: " << topo.fNphysical 
       << "   hyperthreading: " << topo.fHT << std::endl;
    for (auto i=0; i<topo.fNodes; ++i) os << *(topo.fListNodes[i]) << std::endl;
    return os;
}      

/* Get the numa node for a given page address */
//______________________________________________________________________________
int numa_node_addr(void *ptr)
{
  int numa_node = -1;
  get_mempolicy(&numa_node, NULL, 0, ptr, MPOL_F_NODE | MPOL_F_ADDR);
  return numa_node;
}

/** Basic allocator for aligned memory on a given NUMA node.
    Input: number of bytes required and alignment boundaru
    Output: nullptr on error, valid aligned address on success */
//______________________________________________________________________________
void *numa_aligned_malloc(size_t bytes, int node, size_t alignment)
{
  constexpr size_t twosize = 2*sizeof(size_t);
  void *p1 ,*p2; // basic pointer needed for computation.

  /* We need to use numa_malloc provided by C. First we need to allocate memory
  of size bytes + alignment + 2*sizeof(size_t) . We need 'bytes' because 
  user requested it. We need to add 'alignment' because malloc can give us 
  any address and we need to find multiple of 'alignment', so at maximum multiple
  of alignment will be 'alignment' bytes away from any location. We need 
  '2*sizeof(size_t)' for implementing 'aligned_free', since we are returning modified 
  memory pointer, not given by malloc, to the user, we must free the memory 
  allocated by malloc not anything else. So I am storing address given by malloc just above 
  pointer returning to user. Thats why I need extra space to store that address. 
  Then I am checking for error returned by malloc, if it returns NULL then 
  aligned_malloc will fail and return NULL.
  */
  if((p1 =(void *) numa_alloc_onnode(bytes + alignment + twosize, node))==NULL)
  return NULL;

  /*	Next step is to find aligned memory address multiples of alignment.
  By using basic formule I am finding next address after p1 which is 
  multiple of alignment.I am storing new address in p2.
  */
  size_t addr=(size_t)p1+alignment+twosize;

  p2=(void *)(addr - (addr%alignment));

  /*	Final step, I am storing the address returned by malloc 'p1' just "size_t"
  bytes above p2, which will be useful while calling aligned_free.
  */
  *((size_t *)p2-2)=(size_t)p1;
  *((size_t *)p2-1)=bytes + alignment + twosize;

  return p2;
}

/** Frees memory allocated using numa_aligned_malloc*/
//______________________________________________________________________________
void numa_aligned_free(void *p )
{
/**	Find the address stored by aligned_malloc ,"size_t" bytes above the 
current pointer then free it using normal free routine provided by C.
*/
  numa_free((void *)(*((size_t *)p - 2)), *((size_t *)p - 1));
}

} // Geant

#endif
