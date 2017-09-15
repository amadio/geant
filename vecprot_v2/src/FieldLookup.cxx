#include "GUVField.h"
#include "FieldLookup.h"

#include "GeantTaskData.h"
#include "GeantTrackVec.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE
void FieldLookup::GetFieldValue(Geant::GeantTaskData *td, vecgeom::Vector3D<double> Position, double B[3], double *bmag)
{
      // using ThreeVector_f = vecgeom::Vector3D<float>;
      using ThreeVector_d = vecgeom::Vector3D<double>;
      
      if( bmag ) *bmag= 0.0;
      ThreeVector_d MagFldD;  //  Transverse wrt direction of B
      
      if( td->fBfieldIsConst ) {
         MagFldD= td->fConstFieldValue;
        if( bmag ) *bmag=   td->fBfieldMag;
      }
      else
      {
         td->fFieldObj->GetFieldValue(Position, MagFldD);
         if( bmag ) *bmag= MagFldD.Mag();
         
         // printf(" GeantTrack_v::GetFieldValue>  Field at x,y,z= ( %f %f %f ) is (%f %f %f) kGauss - mag = %f \n",
         //       fXposV[i], fZposV[i], fZposV[i], MagFldF.x(), MagFldF.y(), MagFldD.z(), *bmag );
         /***
             int oldPrec= cout.precision(3);
             std::cout << " GeantTrack_v::GetFieldValue>  Field at x,y,z= ( "
             << std::setw(8) << fXposV[i]<< " , " << std::setw(8) << fYposV[i]
             << " , " << std::setw(8) << fZposV[i] << " ) is ( "
             << std::setw(8) << MagFldD.x() << " , " << std::setw(8) << MagFldD.y()
             << " , " << std::setw(8) << MagFldD.z() <<  " ) kGauss - "
             << "mag = " << *bmag << std::endl;
             cout.precision(oldPrec);
         ****/
      }
      B[0]= MagFldD.x();
      B[1]= MagFldD.y();
      B[2]= MagFldD.z();   
}

}
} // End of namespace Geant
