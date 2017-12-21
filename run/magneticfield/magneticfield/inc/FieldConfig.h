//===--- GeantTrack.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file   FieldConfig.h
 * @brief  Class to hold configuration of field, including pointer to global field
 * @author John Apostolakis
 * @date   January 2018
 */
//===----------------------------------------------------------------------===//

#ifndef FIELD_CONFIG_H
#define FIELD_CONFIG_H

#include "Geant/Typedefs.h"
#include "VVectorField.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

// This class owns the field Object which it is passed.
          
class FieldConfig
{
   public:

   FieldConfig() : fFieldObj(nullptr), fConstFieldValue( vecgeom::Vector3D<double>(0., 0., 0.) ),
                   fBfieldMag(0.0), fBfieldIsConst(false)
      { }

   /* @brief Register an existing field (object).  FieldConfig assumes ownership */
   inline
   FieldConfig( VVectorField* vf, bool isUniform );

   /* @brief Create a uniform field. If its magnitude is 0.0, the field will be ignored. */
   FieldConfig( vecgeom::Vector3D<double> fieldValue ) :FieldConfig()  { SetUniformField( fieldValue ); }

   ~FieldConfig() { delete fFieldObj; }

   /* @brief Overwrite with an existing field class.  FieldConfig assumes ownership */
   inline void SetField( VVectorField *fieldObj );

   /* @brief Overwrite with a uniform (constant) field. 
              Note: If its magnitude is 0.0, the field will be ignored. */
   inline void SetUniformField( vecgeom::Vector3D<double> fieldValue );

   /* @brief Overwrite with a uniform (constant) field.  Use field object' value. */ 
   inline void SetUniformField( VVectorField *fieldObject );   

   /* @brief Get the value of the Uniform field  (if configured - else it not define. ) */
   vecgeom::Vector3D<double> GetUniformFieldValue() { return fConstFieldValue; }
   
   /* @brief Check whether a field exists, i.e. is one registered. */   
   bool FieldExists() const;

   /* @brief Ensure that either a uniform field is set or a field class is registered. */
   bool CheckConfig() const;

   /* @brief Check whether field is uniform. */   
   bool IsFieldUniform() const { return fBfieldIsConst; }

   /* @brief Obtain pointer to global field object (if any) */
   VVectorField       *GetFieldObject()       { return fFieldObj; }
   const VVectorField *GetFieldObject() const { return fFieldObj; }
   
   double    GetUniformFieldMag() const { return fBfieldMag; }

private:
   VVectorField            * fFieldObj= nullptr; /** point to field class object */
   vecgeom::Vector3D<double> fConstFieldValue = {0. , 0., 0.};
   double                    fBfieldMag = 0.0;  /** Magnitude of field in case of const field [kiloGauss] */
   bool                      fBfieldIsConst= false;
};

//______________________________________________________________________________          
inline
FieldConfig::FieldConfig( VVectorField* vf, bool isUniform )
   : FieldConfig()
{
   if( isUniform ) {
      SetUniformField( vf );
   } else {
      SetField( vf );            
   }
}
          
//______________________________________________________________________________          
inline
bool FieldConfig::FieldExists() const
{
   bool hasField= false;
   if( fBfieldIsConst ) { hasField = ( fBfieldMag != 0.0 ); }
   else {                 hasField = ( fFieldObj != nullptr); }
   return hasField;
}

//______________________________________________________________________________                    
inline
void FieldConfig::SetField( VVectorField *fieldObj )
{
   delete fFieldObj;
   fFieldObj = fieldObj;
   fBfieldIsConst = false;
   fConstFieldValue = vecgeom::Vector3D<double>(0., 0., 0.);  // Not relevant
   fBfieldMag = fConstFieldValue.Mag();
}

//______________________________________________________________________________          
inline
void FieldConfig::SetUniformField( vecgeom::Vector3D<double> fieldValue )
{
   delete fFieldObj;
   fConstFieldValue= fieldValue;
   fBfieldMag= fieldValue.Mag();  // If its value is 0.0, the field will be ignored
   fFieldObj = nullptr;
   fBfieldIsConst= true;
}

//______________________________________________________________________________          
inline
void FieldConfig::SetUniformField( VVectorField *fldUnif )
{
   vecgeom::Vector3D<double> position( 0.0, 0.0, 0.0);   
   vecgeom::Vector3D<double> fieldValue( 0.0, 0.0, 0.0);
   delete fFieldObj;
   fFieldObj = fldUnif;
   
   if( fFieldObj ) {
      fFieldObj->ObtainFieldValue( position, fieldValue );
   }
   fConstFieldValue= fieldValue;
   fBfieldIsConst= true;
   fBfieldMag= fieldValue.Mag();  // If its value is 0.0, the field will be ignored
}

}
} // namespace Geant 

#endif
