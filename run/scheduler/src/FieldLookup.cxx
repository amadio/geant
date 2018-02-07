#include <cassert>

#include "VVectorField.h"
#include "FieldLookup.h"

#include "FieldConfig.h"

#include "Geant/Error.h"
#include "GeantTaskData.h"
#include "GeantTrackVec.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

FieldConfig* FieldLookup::fFieldConfig= nullptr;

VECCORE_ATT_HOST_DEVICE
void FieldLookup::GetFieldValue( const vecgeom::Vector3D<double>& Position,
                                       vecgeom::Vector3D<double>& MagFieldOut,
                                       double                   & bmagOut
                                 // , const Geant::GeantTaskData     * td  => Not needed !!
   )
{
//   auto tkp = td->fPropagator;
//   auto config = tkp ? tkp->fConfig : nullptr;
   double bmag= 0.0;
   assert( fFieldConfig != nullptr );

   auto pField= fFieldConfig->GetFieldObject();
   
   if( fFieldConfig->IsFieldUniform() ) {
      MagFieldOut= fFieldConfig->GetUniformFieldValue();
      bmag=        fFieldConfig->GetUniformFieldMag();
   }
   else // if( pField )
   {
      assert( pField != nullptr );
      
      pField->ObtainFieldValue(Position, MagFieldOut);
      bmag= MagFieldOut.Mag();
      // printf(" GeantTrack_v::GetFieldValue> Field at ( %f %f %f ) is (%f %f %f) kGauss - mag = %f \n",
         //       Position.x(), Position.y(), Position.z(), MagFldD.x(), MagFldD.y(), MagFldD.z(), *bmag );
   }
   bmagOut= bmag;
}

#if 0
//*****          
//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldLookup::GetFieldValue( const vecgeom::Vector3D<double>& Position,
                                       double                     Bfield[3],   // Out
                                       double                   & bmag,
                                 const Geant::GeantTaskData     * td
   )
{
   using ThreeVector_d = vecgeom::Vector3D<double>;
   ThreeVector_d  MagFldD(0., 0., 0.);
   GetFieldValue( td, Position, MagFldD, bmag );
   B[0]= MagFldD.x();
   B[1]= MagFldD.y();
   B[2]= MagFldD.z(); 
}
#endif
//*****

}
} // End of namespace Geant
