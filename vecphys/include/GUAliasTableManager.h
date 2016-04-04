#ifndef GUALIASTABLEMANAGER_H
#define GUALIASTABLEMANAGER_H 1

#include "GUAliasTable.h"
#include "GUConstants.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

//should be template for precision choices
class GUAliasTableManager
{
public:
  VECCORE_CUDA_HOST
  GUAliasTableManager(int nelements, int ngrid);

  VECCORE_CUDA_HOST
  ~GUAliasTableManager();

  VECCORE_CUDA_HOST
  void AddAliasTable(int Z, GUAliasTable* table);

  VECPHYS_CUDA_HEADER_BOTH
  GUAliasTable* GetAliasTable(int Z);

  VECPHYS_CUDA_HEADER_BOTH
  int GetTableIndex(int nelement);

  VECPHYS_CUDA_HEADER_BOTH
  int GetNumberOfElements();

  VECCORE_CUDA_HOST
  int SizeOfManager();

#ifdef VECCORE_NVCC
  void Relocate(void *devPtr);
#endif

private:
  VECCORE_CUDA_HOST
  void SetTableIndex(int Z);

private:
  //  MaterialHandler *fMaterialHandler;
  GUAliasTable** fAliasTables;
  int            fIndex[maximumZ];
  int            fNElement;        // the number of elements
} ;

} // end namespace impl
} // end namespace vecphys

#endif
