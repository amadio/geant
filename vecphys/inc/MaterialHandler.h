#ifndef MaterialHandler_H
#define MaterialHandler_H 1

#include "base/VPGlobal.h"
#include "GUConstants.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class MaterialHandler
{
public:

  VECCORE_CUDA_HOST
  static MaterialHandler* Instance();

  VECCORE_CUDA_HOST
  MaterialHandler();

  VECCORE_CUDA_HOST
  ~MaterialHandler();

public:
  VECCORE_CUDA_HOST_DEVICE
  int GetNumberOfElements() { return fNumberOfElements; }

  VECCORE_CUDA_HOST_DEVICE
    int* GetElementArray() { return &fElementArray[0]; }

  //a temporary method for the purpose of validation/benchmarking
  VECCORE_CUDA_HOST
  void PrepareTargetElements(int *targetElements, int ntracks, int elementMode = 0);

private:
  VECCORE_CUDA_HOST
  void BuildElementTable();

  VECCORE_CUDA_HOST
  void AddElement(int element);

private:
  static MaterialHandler* fInstance;
  int fElementMode;
  int fNumberOfElements;
  int fElementArray[maximumZ];
};

} // end namespace impl
} // end namespace vecphys

#endif
