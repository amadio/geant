/// \file AlignedBase.h
/// \author Johannes de Fine Licht (johannes.definelicht@cern.ch)

#ifndef ALIGNEDBASE_H_
#define ALIGNEDBASE_H_

#include "base/Global.h"
#include <Vc/Vc>


#ifdef VECGEOM_VC
  // unfortunately the version macros have changed in Vc over time
  // so I am checking which one exist
#ifdef Vc_VERSION_NUMBER
    class AlignedBase : public Vc::VectorAlignedBase {
#endif
#ifdef VC_VERSION_NUMBER
#if VC_VERSION_NUMBER >= VC_VERSION_CHECK(0,99,71)
    class AlignedBase : public Vc::VectorAlignedBase<Vc::Vector<double> > {
#else
    class AlignedBase : public Vc::VectorAlignedBase {
#endif
#endif
    public:
      virtual ~AlignedBase() {}
};
#elif !defined(VECGEOM_NVCC)
class AlignedBase {

public:

  inline
  void *operator new(size_t size) {
    return _mm_malloc(size, kAlignmentBoundary);
  }

  inline
  void *operator new(size_t, void *p) {
    return p;
  }

  inline
  void *operator new[](size_t size) {
    return _mm_malloc(size, kAlignmentBoundary);
  }

  inline
  void *operator new[](size_t , void *p) {
    return p;
  }
  
  inline
  void operator delete(void *ptr, size_t) {
    _mm_free(ptr);
  }

  inline
  void operator delete(void *, void *) {}

  inline
  void operator delete[](void *ptr, size_t) {
    _mm_free(ptr);
  }

  inline
  void operator delete[](void *, void *) {}

};
#else
class AlignedBase {};
#endif


#endif // ALIGNEDBASE_H_
