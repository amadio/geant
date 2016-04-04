#ifndef MaterialHandler_H
#define MaterialHandler_H 1

#include "backend/Backend.h"
#include "GUConstants.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class MaterialHandler 
{
public:

  VECPHYS_CUDA_HEADER_HOST
  static MaterialHandler* Instance();

  VECPHYS_CUDA_HEADER_HOST
  MaterialHandler();

  VECPHYS_CUDA_HEADER_HOST
  ~MaterialHandler();

public:
  VECPHYS_CUDA_HEADER_BOTH
  int GetNumberOfElements() { return fNumberOfElements; }

  VECPHYS_CUDA_HEADER_BOTH
    int* GetElementArray() { return &fElementArray[0]; }

  //a temporary method for the purpose of validation/benchmarking
  VECPHYS_CUDA_HEADER_HOST
  void PrepareTargetElements(int *targetElements, int ntracks, int elementMode = 0);

private:
  VECPHYS_CUDA_HEADER_HOST
  void BuildElementTable();

  VECPHYS_CUDA_HEADER_HOST
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
