#ifndef MaterialHandler_H
#define MaterialHandler_H 1

#include "base/VecPhys.h"
#include "GUConstants.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class MaterialHandler
{
public:
  VECCORE_CUDA_HOST
  static MaterialHandler *Instance();

  VECCORE_CUDA_HOST
  MaterialHandler();

  VECCORE_CUDA_HOST
  ~MaterialHandler();

public:
  VECCORE_CUDA_HOST
  int GetNumberOfElements() { return fNumberOfElements; }

  VECCORE_CUDA_HOST
    int* GetElementArray() { return &fElementArray[0]; }

  //temporary methods for the purpose of validation/benchmarking
  VECCORE_CUDA_HOST
  void PrepareTargetElements(int *targetElements, int ntracks, int elementMode = 0);

  VECCORE_CUDA_HOST
  void PrepareMaterialIndex(int *materialIndex, int ntracks, int materialMode = 0);

  VECCORE_CUDA_HOST
  void BuildElementTable();

  VECCORE_CUDA_HOST
  void BuildMaterialTable();

private:
  VECCORE_CUDA_HOST
  void AddElement(int element);

private:
  static MaterialHandler *fInstance;
  int fElementMode;
  int fNumberOfElements;
  int fElementArray[maximumZ];

};

} // end namespace impl
} // end namespace vecphys

#endif
