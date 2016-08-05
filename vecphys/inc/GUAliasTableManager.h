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
  VECCORE_ATT_HOST
  GUAliasTableManager(int nelements, int ngrid);

  VECCORE_ATT_HOST
  ~GUAliasTableManager();

  VECCORE_ATT_HOST
  void AddAliasTable(int Z, GUAliasTable* table);

  VECCORE_ATT_HOST_DEVICE
  GUAliasTable* GetAliasTable(int Z);

  VECCORE_ATT_HOST_DEVICE
  int GetTableIndex(int nelement);

  VECCORE_ATT_HOST_DEVICE
  int GetNumberOfElements();

  VECCORE_ATT_HOST
  int SizeOfManager();

#ifdef VECCORE_CUDA
  void Relocate(void *devPtr);
#endif

private:
  VECCORE_ATT_HOST
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
