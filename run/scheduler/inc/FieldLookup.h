//===--- Track.h - GeantV ---------------------------------*- C++ -*-===//
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
#include "VVectorField.h"

#include "FieldConfig.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class TaskData;
   
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
                             double                   & bmag // ,
                       // const TaskData            * td
      );

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
                             double                    & bmag // ,
                       // const geant::TaskData      * td                          
         );
#endif

   static void SetFieldConfig( FieldConfig* fldCfg ) { fFieldConfig = fldCfg; }

   static FieldConfig* /*const*/ GetFieldConfig() { return fFieldConfig; }
   
private:
   // static VVectorField               *fFieldObj;         // To get value of the field!
   // static vecgeom::Vector3D<double>  fConstFieldValue;   // Value - if field is constant.
   // static bool                       fBfieldIsConst;      /** Flag - is the B field constant ?  */
   static FieldConfig* fFieldConfig;
};

}
} // namespace geant 

#endif
