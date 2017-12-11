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
#include "VScalarField.h"

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
    FieldLookup() {} //  = delete; ??
    ~FieldLookup() {}

   /**
    * @brief Function that return magnetic field at a coordinates 'Position'
    * @param  Position       Location ( in global coordinates )
    * @param  MagFieldValue  Output magnetic field vector value (global coordinates)
    * @param  bmag           Output field magnitude
    */
   static
   VECCORE_ATT_HOST_DEVICE   
   void GetFieldValue( const vecgeom::Vector3D<double>& Position,
                             vecgeom::Vector3D<double>& MagFieldValue, // Out
                             double                   & bmag,
                       const GeantTaskData            * td );

#if 0   
   /**
    * @brief Function that return magnetic field at a coordinates 'Position'
    * @param  Position      Location ( in global coordinates )
    * @param  BfieldOut[3]  Output magnetic field vector value (global coordinates)
    * @param  bmag          Output field magnitude
    */   
   static
   VECCORE_ATT_HOST_DEVICE
   void GetFieldValue( const vecgeom::Vector3D<double> & Position,
                             double                      BfieldOut[3],
                             double                    & bmag,
                       const Geant::GeantTaskData      * td                          
         );
#endif

   /* @brief Ensure that either a uniform field is set or a field class is registered. */
   static   
      bool CheckConfig( const Geant::GeantTaskData * td );
};

}
} // namespace Geant 

#endif
