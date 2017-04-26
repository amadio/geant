#ifndef MaterialHandler_H
#define MaterialHandler_H 1

#include "GUConstants.h"
#include "base/VecPhys.h"
#include <vector>

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class MaterialHandler {
public:
  VECCORE_ATT_HOST
  static MaterialHandler *Instance();

  VECCORE_ATT_HOST
  MaterialHandler();

  VECCORE_ATT_HOST
  ~MaterialHandler();

public:
  VECCORE_ATT_HOST
  int GetNumberOfElements() { return fNumberOfElements; }

  VECCORE_ATT_HOST
  int *GetElementArray() { return &fElementArray[0]; }

  // temporary methods for the purpose of validation/benchmarking
  VECCORE_ATT_HOST
  void PrepareTargetElements(int *targetElements, int ntracks, int elementMode = 0);

  VECCORE_ATT_HOST
  void PrepareMaterialIndex(int *materialIndex, int ntracks, int materialMode = 0);

  VECCORE_ATT_HOST
  void BuildElementTable();

  VECCORE_ATT_HOST
  void BuildElementTable(std::vector<int> list);

  VECCORE_ATT_HOST
  void BuildMaterialTable();

private:
  VECCORE_ATT_HOST
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
