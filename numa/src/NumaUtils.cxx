#include "NumaUtils.h"

#ifdef USE_NUMA
  #include <numa.h>
  #include <numaif.h>
#else
  #ifdef __INTEL_COMPILER
  #include <immintrin.h>
  #else
  #include "mm_malloc.h"
  #endif
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void *NumaAlignedMalloc(size_t bytes, int node, size_t alignment)
{
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
#ifdef USE_NUMA
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
#else
  (void)node;
  return _mm_malloc(bytes, alignment);
#endif
#else
  return malloc(bytes);
#endif
}

VECCORE_ATT_HOST_DEVICE
void NumaAlignedFree(void *p )
{
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
// Find the address stored by aligned_malloc ,"size_t" bytes above the
// current pointer then free it using normal free routine provided by C.
#ifdef USE_NUMA
  if (p) numa_free((void *)(*((size_t *)p - 2)), *((size_t *)p - 1));
#else
  _mm_free(p);
#endif
#else
  free(p);
#endif
}

int NumaNodeAddr(void *ptr)
{
// Returns the NUMA node id associated with the memory pointer.
#ifdef USE_NUMA
  (void)ptr;
  int numa_node = -1;
  get_mempolicy(&numa_node, NULL, 0, ptr, MPOL_F_NODE | MPOL_F_ADDR);
  return numa_node;
#else
  (void)ptr;
  return 0;
#endif  
}

void PinToCore(size_t core)
{
// Pin caller thread to a given physical core
#ifdef USE_NUMA
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(core, &cpuset);
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#else
  (void)core;
#endif
}

} // GEANT_IMPL_NAMESPACE
} // Geant
