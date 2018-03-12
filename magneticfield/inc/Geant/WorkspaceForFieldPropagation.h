#ifndef Workspace_For_Field_Propagation_h
#define Workspace_For_Field_Propagation_h

#include "base/SOA3D.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

struct WorkspaceForFieldPropagation {
public:
  WorkspaceForFieldPropagation(size_t bufferSize);
  ~WorkspaceForFieldPropagation();

  GEANT_FORCE_INLINE
  size_t capacity() { return fCapacity; }

  /** @brief  Resize - likely to ignore old content */
  GEANT_FORCE_INLINE
  void Resize(size_t size);

  /** @brief  Throw away old content */
  GEANT_FORCE_INLINE
  void Clear();

  /** @brief  Enlarge the buffers, throwing away old content */
  GEANT_FORCE_INLINE
  void ClearAndResize(size_t numEntries);

  GEANT_FORCE_INLINE
  bool CheckSize(size_t numNeeded);
  // Invariant:
  //    All SOA3D containers must have the same size.
  //
  // Consequence:
  //  -> the size of fPositionInpForFieldProp with represent others
public:
  size_t                  fCapacity      = 0;        // container capacity
  double                 *fChargeInp     = nullptr;  // charge array
  double                 *fMomentumInp   = nullptr;  // momentum array
  double                 *fStepsInp      = nullptr;  // steps array
  vecgeom::SOA3D<double> *fPositionInp   = nullptr;  // position array
  vecgeom::SOA3D<double> *fDirectionInp  = nullptr;  // direction array
  vecgeom::SOA3D<double> *fPositionOutp  = nullptr;  // new position array
  vecgeom::SOA3D<double> *fDirectionOutp = nullptr;  // new direction array
};

GEANT_FORCE_INLINE
WorkspaceForFieldPropagation::WorkspaceForFieldPropagation(size_t bufferSize)
{
  // using SoADd = vecgeom::SOA3D<type>;
  assert(bufferSize > 0);

  fChargeInp     = vecgeom::AlignedAllocate<double>(bufferSize);
  fMomentumInp   = vecgeom::AlignedAllocate<double>(bufferSize);
  fStepsInp      = vecgeom::AlignedAllocate<double>(bufferSize);
  fPositionInp   = new vecgeom::SOA3D<double>(bufferSize);
  fDirectionInp  = new vecgeom::SOA3D<double>(bufferSize);
  fPositionOutp  = new vecgeom::SOA3D<double>(bufferSize);
  fDirectionOutp = new vecgeom::SOA3D<double>(bufferSize);
}

GEANT_FORCE_INLINE
WorkspaceForFieldPropagation::~WorkspaceForFieldPropagation()
{
  vecgeom::AlignedFree(fChargeInp);
  vecgeom::AlignedFree(fMomentumInp);
  vecgeom::AlignedFree(fStepsInp);
  delete fPositionInp;
  delete fDirectionInp;
  delete fPositionOutp;
  delete fDirectionOutp;
}

GEANT_FORCE_INLINE
void WorkspaceForFieldPropagation::Clear()
{
  fPositionInp->clear();
  fPositionOutp->clear();
  fDirectionInp->clear();
  fDirectionOutp->clear();
}

GEANT_FORCE_INLINE
void WorkspaceForFieldPropagation::ClearAndResize(size_t numEntries)
{
  Clear();

  if (numEntries > fCapacity) {
    fCapacity = std::max(2 * fCapacity, ((numEntries + 1) / 16 + 1) * 16);

    vecgeom::AlignedFree(fChargeInp);
    vecgeom::AlignedFree(fMomentumInp);
    vecgeom::AlignedFree(fStepsInp);
    fChargeInp   = vecgeom::AlignedAllocate<double>(fCapacity);
    fMomentumInp = vecgeom::AlignedAllocate<double>(fCapacity);
    fStepsInp    = vecgeom::AlignedAllocate<double>(fCapacity);
    fPositionInp->reserve(fCapacity);
    fPositionOutp->reserve(fCapacity);
    fDirectionInp->reserve(fCapacity);
    fDirectionOutp->reserve(fCapacity);
  }
}

GEANT_FORCE_INLINE
void WorkspaceForFieldPropagation::Resize(size_t size)
{
  assert(size < fCapacity);
  fPositionInp->resize(size);
  fPositionOutp->resize(size);
  fDirectionInp->resize(size);
  fDirectionOutp->resize(size);
}

GEANT_FORCE_INLINE
bool WorkspaceForFieldPropagation::CheckSize(size_t numNeeded)
{
  bool goodInp = fPositionInp && (fPositionInp->capacity() >= numNeeded);
  assert(goodInp && "Bad capacity of PositionInp in Workspace for Field Propagation.");

  bool goodOutp = fPositionOutp && (fPositionOutp->capacity() >= numNeeded);
  assert(goodOutp && "Bad capacity of PositionOutp in Workspace for Field Propagation.");

  bool goodInpDir = fDirectionInp && (fDirectionInp->capacity() >= numNeeded);
  assert(goodInpDir && "Bad capacity of DirectionInp in Workspace for Field Propagation.");

  bool goodOutpDir = fDirectionOutp && (fDirectionOutp->capacity() >= numNeeded);
  assert(goodOutpDir && "Bad capacity of DirectionOutp in Workspace for Field Propagation.");

  return goodInp && goodOutp && goodInpDir && goodOutpDir;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
#endif
