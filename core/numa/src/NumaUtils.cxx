#include "NumaUtils.h"

#include <iostream>
#ifndef GEANT_USE_NUMA
  #ifdef __INTEL_COMPILER
  #include <immintrin.h>
  #else
  #include "mm_malloc.h"
  #endif
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)

NumaUtilsStruct *NumaUtilsStruct::fgInstance = nullptr;

//______________________________________________________________________________
NumaUtilsStruct *NumaUtilsStruct::Instance() {
  if (fgInstance) return fgInstance;
  return (new NumaUtilsStruct());
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
NumaUtilsStruct::NumaUtilsStruct()
{
  fgInstance = this;
  fAvailable = false;

  std::lock_guard<std::mutex> lock(fLock);
  int err;
  hwloc_bitmap_t set = hwloc_bitmap_alloc();
  if (!set) {
    std::cout << "*** Failed to allocate a bitmap\n*** NUMA disabled.\n";
    return;
  }
  err = hwloc_topology_init(&fTopology);
  if (err < 0) {
    std::cout << "*** Failed to initialize the topology\n*** NUMA disabled.\n";
    hwloc_bitmap_free(set);
    return;
  }
  err = hwloc_topology_load(fTopology);
  if (err < 0) {
    std::cout << "*** Failed to load the topology.\n*** NUMA disabled.\n";
    hwloc_bitmap_free(set);
    hwloc_topology_destroy(fTopology);
    return;
  }
  hwloc_const_bitmap_t cset;
  // retrieve the entire set of NUMA nodes and count them
  cset = hwloc_topology_get_topology_nodeset(fTopology);
  int nnodes = hwloc_bitmap_weight(cset);
  if (nnodes <= 0) {
    std::cout << "*** This machine is not NUMA, nothing to do\n";
    hwloc_bitmap_free(set);
    hwloc_topology_destroy(fTopology);
    return;
  }
  // get the process memory binding as a nodeset
  hwloc_membind_policy_t policy;
  err = hwloc_get_membind_nodeset(fTopology, set, &policy, 0);
  if (err < 0) {
    std::cout << "*** Failed to retrieve my memory binding and policy.\n*** NUMA disabled.\n";
    hwloc_bitmap_free(set);
    hwloc_topology_destroy(fTopology);
    return;
  }
  hwloc_obj_t obj;
  unsigned i;
  char *s = nullptr;
  hwloc_bitmap_asprintf(&s, set);
  std::cout << "*** Bound to nodeset " << s << " with content:\n";
  delete s;
  hwloc_bitmap_foreach_begin(i, set) {
    obj = hwloc_get_numanode_obj_by_os_index(fTopology, i);
    std::cout << "  node #" << obj->logical_index << " (OS index " << i
              << ") with " << obj->memory.local_memory << " bytes of memory\n";
  } hwloc_bitmap_foreach_end();

  /* Test memory binding */
  const struct hwloc_topology_support *support = hwloc_topology_get_support(fTopology);
  char *buffer = nullptr;
  if (support->membind->replicate_membind) {
    std::cout << "***    replicated memory binding is supported\n";
    buffer = (char*)hwloc_alloc_membind_nodeset(fTopology, 4096, cset, HWLOC_MEMBIND_REPLICATE, HWLOC_MEMBIND_STRICT);
  } else {
    std::cout << "***   replicated memory binding is NOT supported\n";
  }
    
  obj = nullptr;
  while ((obj = hwloc_get_next_obj_by_type(fTopology, HWLOC_OBJ_NUMANODE, obj)) != nullptr) {
    //obj = hwloc_get_numanode_obj_by_os_index(fTopology, obj->logical_index);
    buffer = (char*)hwloc_alloc_membind_nodeset(fTopology, 4096, obj->nodeset, HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_STRICT);
    if (!buffer) {
      std::cout << "*** Failed to manually allocate memory on node " << obj->os_index
                << "\n*** NUMA disabled.\n";
      hwloc_bitmap_free(set);
      hwloc_topology_destroy(fTopology);
      return;
    }
    // check where the memory was allocated
    hwloc_bitmap_zero(set);
    err = hwloc_get_area_membind_nodeset(fTopology, buffer, 4096, set, &policy, 0);
    if (err < 0) {
      std::cout << "*** Failed to retrieve the buffer binding and policy.\n*** NUMA disabled.\n";
      hwloc_bitmap_free(set);
      hwloc_topology_destroy(fTopology);
      return;
    }
    /* check the binding policy, it should be what we requested above,
     * but may be different if the implementation of different policies
     * is identical for the current operating system.
     */
    std::cout << "***   buffer membind policy is " << policy << " while we requested "
              << HWLOC_MEMBIND_REPLICATE << " or " << HWLOC_MEMBIND_BIND << std::endl;
    // print the corresponding NUMA nodes
    hwloc_bitmap_foreach_begin(i, set) {
    obj = hwloc_get_numanode_obj_by_os_index(fTopology, i);
    std::cout << "***   buffer" << i << " bound to node #" << obj->logical_index << std::endl;
    } hwloc_bitmap_foreach_end();
  }
  hwloc_bitmap_free(set);
  fAvailable = true;  
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
NumaUtilsStruct::~NumaUtilsStruct()
{
  // Destroy topology object
  if (fAvailable)
    hwloc_topology_destroy(fTopology);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void *NumaUtilsStruct::NumaAlignedMalloc(size_t bytes, int node, size_t alignment)
{
// Fallback to 
// Basic allocator for aligned memory on a given NUMA node.
//      Input: number of bytes required and alignment boundaru
//      Output: nullptr on error, valid aligned address on success
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

  hwloc_obj_t obj = hwloc_get_numanode_obj_by_os_index(fTopology, node);
  p1 = hwloc_alloc_membind_nodeset(fTopology, bytes + alignment + twosize,
                                   obj->nodeset, HWLOC_MEMBIND_BIND, HWLOC_MEMBIND_STRICT);
  if (p1 == nullptr)
    return nullptr;
//  if((p1 =(void *) numa_alloc_onnode(bytes + alignment + twosize, node))==NULL)
//  return NULL;

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
  
//  assert(NumaNodeAddr(p2) == node);
  return p2;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void NumaUtilsStruct::NumaAlignedFree(void *p )
{
// Find the address stored by aligned_malloc ,"size_t" bytes above the
// current pointer then free it using normal free routine provided by C.
  if (p) hwloc_free(fTopology, (void *)(*((size_t *)p - 2)), *((size_t *)p - 1));
//  if (p) numa_free((void *)(*((size_t *)p - 2)), *((size_t *)p - 1));
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int NumaUtilsStruct::NumaNodeAddr(void *ptr)
{
// Returns the NUMA node id associated with the memory pointer.
  // check where the memory was allocated
  hwloc_bitmap_t set = hwloc_bitmap_alloc();
  hwloc_membind_policy_t policy;
  int err = hwloc_get_membind_nodeset(fTopology, set, &policy, 0);
  if (err < 0) {
    std::cout << "*** Failed to retrieve my memory binding and policy.\n*** NUMA disabled.\n";
    return -1;
  }
  err = hwloc_get_area_membind_nodeset(fTopology, ptr, 1 /* *((size_t *)ptr - 1)*/, set, &policy, 0);
  if (err < 0) {
    std::cout << "*** Failed to retrieve the buffer binding and policy.\n*** NUMA disabled.\n";
    return -1;
  }
  int index = hwloc_bitmap_first(set);
  if (index < 0) {
    std::cout << "*** Failed to determine memory area binding\n";
    return -1;  
  }
  hwloc_obj_t obj = hwloc_get_numanode_obj_by_os_index(fTopology, index);
  return ( obj->logical_index );
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int NumaUtilsStruct::GetCpuBinding() const
{
// Check the current binding of the thread.
  hwloc_bitmap_t set = hwloc_bitmap_alloc();
  int err = hwloc_get_last_cpu_location(fTopology, set, HWLOC_CPUBIND_THREAD);
  if (err < 0) {
    std::cout << "*** Failed to get last cpu location\n";
    hwloc_bitmap_free(set);
    return err;
  }
  assert(hwloc_bitmap_weight(set) == 1);
  // extract the PU OS index from the bitmap
  unsigned i = hwloc_bitmap_first(set);
  hwloc_obj_t obj = hwloc_get_pu_obj_by_os_index(fTopology, i);
  hwloc_bitmap_free(set);
  
  return obj->os_index;
}

#endif // GEANT_USE_NUMA

namespace NumaUtils {
  VECCORE_ATT_HOST_DEVICE
  void *NumaAlignedMalloc(std::size_t bytes, int node, std::size_t alignment)
  {
  #ifndef VECCORE_CUDA_DEVICE_COMPILATION
  #ifdef GEANT_USE_NUMA
    return ( NumaUtilsStruct::Instance()->NumaAlignedMalloc(bytes, node, alignment) );
  #else
    (void)node;
    return _mm_malloc(bytes, alignment);
  #endif
  #else
    (void)node;
    (void)alignment;
    return malloc(bytes);
  #endif
  }
  
  VECCORE_ATT_HOST_DEVICE
  void  NumaAlignedFree(void *p)
  {
  #ifndef VECCORE_CUDA_DEVICE_COMPILATION
  #ifdef GEANT_USE_NUMA
    NumaUtilsStruct::Instance()->NumaAlignedFree(p);
  #else
    _mm_free(p);
  #endif
  #else
    free(p);
  #endif
  }
  
  VECCORE_ATT_HOST_DEVICE
  int NumaNodeAddr(void *ptr)
  {
  #if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
    return NumaUtilsStruct::Instance()->NumaNodeAddr(ptr);
  #else
    (void)ptr;
    return 0;
  #endif
  }
  
  VECCORE_ATT_HOST_DEVICE
  int GetCpuBinding()
  {
  #if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
    return NumaUtilsStruct::Instance()->GetCpuBinding();
  #else
    return -1;
  #endif
  }
  
  VECCORE_ATT_HOST_DEVICE
  bool NumaAvailable()
  {
#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
    return NumaUtilsStruct::Instance()->fAvailable;
#else
    return false;
#endif
  }

#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
  hwloc_topology_t const &Topology()
  {
    return NumaUtilsStruct::Instance()->fTopology;
  }
#endif
}

} // GEANT_IMPL_NAMESPACE
} // Geant
