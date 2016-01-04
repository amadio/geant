//===--- GeantTrack.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GUFieldAux.h
 * @brief Implementation of GetFieldValue static method
 * @author John Apostolakis
 */
//===----------------------------------------------------------------------===//

#ifndef FIELD_LOOKUP_H
#define FIELD_LOOKUP_H

#include "Geant/Typedefs.h"
#include "GUVField.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTaskData;
   
/**
 * @brief  Enable case of constant field without field access method
 *
 */
class FieldLookup
{
   public:
     FieldLookup() {}
     ~FieldLookup() {}
   
   /**
    * @brief Function that return magnetic field at a coordinates 'Position'
    * @param  Position      Location ( in global coordinates )
    * @param  BfieldOut[3]  Output magnetic field vector value (global coordinates)
    * @param  bmagOut       Output (optional) field magnitude
    */   

   static
      VECCORE_ATT_HOST_DEVICE
      void GetFieldValue(GeantTaskData *td, vecgeom::Vector3D<double> Position, double B[3], double *bmag);
};

}
} // namespace Geant 

#endif
