#ifndef GEANT_VECTOR_TYPES_H
#define GEANT_VECTOR_TYPES_H

#include <base/Global.h>
#include <base/Vector3D.h>
#include <Geant/Config.h>

namespace Geant {

  typedef vecgeom::VectorBackend::Float_v Float_v;
  typedef vecgeom::VectorBackend::Real_v Double_v;
  typedef vecgeom::VectorBackend::Int_v Int_v;

  typedef vecCore::Mask_v<Float_v> MaskF_v;
  typedef vecCore::Mask_v<Double_v> MaskD_v;
  typedef vecCore::Mask_v<Int_v> MaskI_v;

  const auto kVecLenF = vecCore::VectorSize<Float_v>();
  const auto kVecLenD = vecCore::VectorSize<Double_v>();

  GEANT_FORCE_INLINE
  void CopyFltToDbl(Float_v const &flt_v, Double_v &dbl1_v, Double_v &dbl2_v)
  {
  	// Copy the float SIMD lanes into 2 Double_v variables
  	for (size_t lane = 0; lane < kVecLenD; ++lane) {
  	  vecCore::Set(dbl1_v, lane, (double)vecCore::Get(flt_v, lane));
  	  vecCore::Set(dbl2_v, lane, (double)vecCore::Get(flt_v, lane + kVecLenD));
  	}
  }
 
  GEANT_FORCE_INLINE
  void CopyFltToDbl(vecgeom::Vector3D<Float_v> const &flt_v,
  	                vecgeom::Vector3D<Double_v> &dbl1_v,
  	                vecgeom::Vector3D<Double_v> &dbl2_v)
  {
  	// Copy the float SIMD lanes into 2 Double_v variables
  	for (size_t lane = 0; lane < kVecLenD; ++lane) {
  	  for (size_t i = 0; i < 3; ++i) {  	  	
  	    vecCore::Set(dbl1_v[i], lane, (double)vecCore::Get(flt_v[i], lane));
  	    vecCore::Set(dbl2_v[i], lane, (double)vecCore::Get(flt_v[i], lane + kVecLenD));
  	  }
  	}
  }

  GEANT_FORCE_INLINE
  void CopyDblToFlt(Double_v const &dbl1_v, Double_v const &dbl2_v, Float_v &flt_v)
  {
 	// Copy the 2 Double_v SIMD lanes into one Float_v variable
   	for (size_t lane = 0; lane < kVecLenD; ++lane) {
  	  vecCore::Set(flt_v, lane, (float)vecCore::Get(dbl1_v, lane));
  	  vecCore::Set(flt_v, lane + kVecLenD, (double)vecCore::Get(dbl2_v, lane));
  	}
  }

  GEANT_FORCE_INLINE
  void CopyDblToFlt(vecgeom::Vector3D<Double_v> const &dbl1_v,
  	                vecgeom::Vector3D<Double_v> const &dbl2_v,
  	                vecgeom::Vector3D<Float_v> &flt_v)
  {
 	// Copy the 2 Double_v SIMD lanes into one Float_v variable
   	for (size_t lane = 0; lane < kVecLenD; ++lane) {
   	  for (size_t i = 0; i < 3; ++i) {  	  	
 	    vecCore::Set(flt_v[i], lane, (float)vecCore::Get(dbl1_v[i], lane));
  	    vecCore::Set(flt_v[i], lane + kVecLenD, (double)vecCore::Get(dbl2_v[i], lane));
      }
  	}
  }

}

#endif
