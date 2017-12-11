#include "VScalarField.h"
#include "FieldLookup.h"

#include "Geant/Error.h"
#include "GeantTaskData.h"
#include "GeantTrackVec.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void FieldLookup::GetFieldValue( const vecgeom::Vector3D<double>& Position,
                                       vecgeom::Vector3D<double>& MagFieldOut,
                                       double                   & bmag,
                                 const Geant::GeantTaskData     * td
   )
{
   // using ThreeVector_f = vecgeom::Vector3D<float>;
   // using ThreeVector_d = vecgeom::Vector3D<double>;
   
   bmag= 0.0;
   if( td->fBfieldIsConst ) {
      MagFieldOut= td->fConstFieldValue;
      bmag=        td->fBfieldMag;
   }
   else
   {
      td->fFieldObj->GetFieldValue(Position, MagFieldOut);
      bmag= MagFieldOut.Mag();
      // printf(" GeantTrack_v::GetFieldValue> Field at ( %f %f %f ) is (%f %f %f) kGauss - mag = %f \n",
         //       Position.x(), Position.y(), Position.z(), MagFldD.x(), MagFldD.y(), MagFldD.z(), *bmag );
   }
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

//______________________________________________________________________________          
bool
FieldLookup::CheckConfig( const Geant::GeantTaskData * td )
{
   bool ok= ( td->fBfieldIsConst || ( td->fFieldObj != nullptr ) );
   if( !ok )
      Error("FieldLookup::CheckConfig",
            "Field configuration incorrect - neither uniform field nor is a field object set.");
   return ok;
}
          
}
} // End of namespace Geant
