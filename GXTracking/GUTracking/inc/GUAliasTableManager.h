#ifndef GUALIASTABLEMANAGER_H
#define GUALIASTABLEMANAGER_H 1

#include "backend/Backend.h"
#include "GUAliasTable.h"
#include "GUConstants.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

//should be template for precision choices
class GUAliasTableManager 
{ 
public:
  VECPHYS_CUDA_HEADER_HOST
  GUAliasTableManager();

  VECPHYS_CUDA_HEADER_HOST
  ~GUAliasTableManager();

  VECPHYS_CUDA_HEADER_HOST
  void AddAliasTable(int Z, GUAliasTable* table);

  VECPHYS_CUDA_HEADER_BOTH
  GUAliasTable* GetAliasTable(int Z);

  VECPHYS_CUDA_HEADER_BOTH
  int GetTableIndex(int nelement);

  VECPHYS_CUDA_HEADER_BOTH
  int GetNumberOfElements();

  VECPHYS_CUDA_HEADER_HOST
  int SizeOfManager();

#ifdef VECPHYS_NVCC
  void Relocate(void *devPtr);
#endif

private:
  VECPHYS_CUDA_HEADER_HOST
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
