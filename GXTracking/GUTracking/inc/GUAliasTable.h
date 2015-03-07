#ifndef GUALIASTABLE_H
#define GUALIASTABLE_H 1

#include "backend/Backend.h"

namespace vecphys {

inline namespace VECPHYS_IMPL_NAMESPACE {

//should be template for precision choices
class GUAliasTable
{ 
 public:
  VECPHYS_FUNC_QUALIFIER
  GUAliasTable(int ngrid);

  VECPHYS_FUNC_QUALIFIER
  GUAliasTable(const GUAliasTable& table);

  VECPHYS_FUNC_QUALIFIER
  GUAliasTable& operator=(const GUAliasTable& table);

  VECPHYS_FUNC_QUALIFIER
  ~GUAliasTable();

  VECPHYS_FUNC_QUALIFIER
  void Allocate(int ngrid);

  VECPHYS_FUNC_QUALIFIER
  void Deallocate();

  VECPHYS_FUNC_QUALIFIER
  void CopyData(const GUAliasTable& table);

#ifdef VECPHYS_NVCC
  void Relocate(void *devPtr);
#endif

  VECPHYS_FUNC_QUALIFIER
  int SizeOfGrid() { return fNGrid; }

  VECPHYS_FUNC_QUALIFIER
  int SizeOfTable();

  VECPHYS_FUNC_QUALIFIER
  void PrintInfo();

  Precision* fpdf;         // original p.d.f distribution
  Precision* fProbQ;       // non-alias probability
  int*       fAlias;       // alias index
  int        fNGrid;       // the number of bins that can be stored
} ;

} // end namespace impl
} // end namespace vecphys

#endif
