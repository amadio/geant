#ifndef GEANT_VECTOR_TYPES_H
#define GEANT_VECTOR_TYPES_H

#include <base/Global.h>

namespace Geant {

  typedef vecgeom::VectorBackend::Float_v Float_v;
  typedef vecgeom::VectorBackend::Real_v Double_v;
  typedef vecgeom::VectorBackend::Int_v Int_v;

  typedef vecCore::Mask_v<Float_v> MaskF_v;
  typedef vecCore::Mask_v<Double_v> MaskD_v;
  typedef vecCore::Mask_v<Int_v> MaskI_v;

  const auto kVecLenF = vecCore::VectorSize<Float_v>();
  const auto kVecLenD = vecCore::VectorSize<Double_v>();

}

#endif
