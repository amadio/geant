#ifndef GUALIASTABLE_H
#define GUALIASTABLE_H 1

#include "base/VecPhys.h"

namespace vecphys {

inline namespace VECPHYS_IMPL_NAMESPACE {

// should be template for precision choices
class GUAliasTable {
public:
  VECCORE_ATT_HOST_DEVICE
  GUAliasTable(int ngrid);

  VECCORE_ATT_HOST_DEVICE
  GUAliasTable(const GUAliasTable &table);

  VECCORE_ATT_HOST_DEVICE
  GUAliasTable &operator=(const GUAliasTable &table);

  VECCORE_ATT_HOST_DEVICE
  ~GUAliasTable();

  VECCORE_ATT_HOST_DEVICE
  void Allocate(int ngrid);

  VECCORE_ATT_HOST_DEVICE
  void Deallocate();

  VECCORE_ATT_HOST_DEVICE
  void CopyData(const GUAliasTable &table);

#ifdef VECCORE_CUDA
  void Relocate(void *devPtr);
#endif

  VECCORE_ATT_HOST_DEVICE
  int SizeOfGrid() { return fNGrid; }

  VECCORE_ATT_HOST_DEVICE
  int SizeOfTable();

  VECCORE_ATT_HOST_DEVICE
  void PrintInfo();

  Real_t *fpdf;   // original p.d.f distribution
  Real_t *fProbQ; // non-alias probability
  int *fAlias;    // alias index
  int fNGrid;     // the number of bins that can be stored
};

} // end namespace impl
} // end namespace vecphys

#endif
