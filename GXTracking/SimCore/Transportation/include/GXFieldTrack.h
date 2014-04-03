#ifndef GXFieldTrack_HH
#define GXFieldTrack_HH

#include "GPTypeDef.h"
#include "GPConstants.h"
#include "GPThreeVector.h"

struct GXFieldTrack 
{
  G4double  SixVector[GXFieldTrack_ncompSVEC]; //position and momentum

  G4double  fDistanceAlongCurve;  // distance along curve of point
  G4double  fKineticEnergy;
  G4double  fRestMass_c2;
  G4double  fCharge;

};


extern "C" {

FQUALIFIER
void GXFieldTrack_Constructor( GXFieldTrack *This, 
                               GPThreeVector pPosition, 
                               GPThreeVector pMomentum,
                               G4double      kineticEnergy,
                               G4double      restMass_c2,
                               G4double      charge, 
                               G4double      curve_length=0.0 );

FQUALIFIER
void GXFieldTrack_UpdateFourMomentum( GXFieldTrack *This,
				      G4double kineticEnergy, 
				      GPThreeVector momentum );

FQUALIFIER
void GXFieldTrack_UpdateState(GXFieldTrack *This, 
			      GPThreeVector position, 
			      GPThreeVector momentum,
			      G4double      kineticEnergy);
                              
FQUALIFIER
GXFieldTrack& GXFieldTrack_SetCurvePnt(GXFieldTrack *This, 
				       GPThreeVector pPosition, 
				       GPThreeVector pMomentum,  
				       G4double      s_curve );

FQUALIFIER
GPThreeVector  GXFieldTrack_GetMomentum(const GXFieldTrack *This);

FQUALIFIER
GPThreeVector  GXFieldTrack_GetPosition(const GXFieldTrack *This);

FQUALIFIER
G4double       GXFieldTrack_GetCurveLength(const GXFieldTrack *This);

FQUALIFIER
void GXFieldTrack_SetPosition(GXFieldTrack *This, GPThreeVector pPosition) ;

FQUALIFIER
void GXFieldTrack_SetMomentum(GXFieldTrack *This, GPThreeVector pMomentum );

FQUALIFIER
void GXFieldTrack_SetCurveLength(GXFieldTrack *This, G4double nCurve_s);

FQUALIFIER
G4double GXFieldTrack_GetKineticEnergy(GXFieldTrack *This);

FQUALIFIER
void GXFieldTrack_SetKineticEnergy(GXFieldTrack *This, G4double newKinEnergy);

FQUALIFIER
void GXFieldTrack_DumpToArray(GXFieldTrack *This,
			      G4double valArr[GXFieldTrack_ncompSVEC] );

FQUALIFIER
void GXFieldTrack_LoadFromArray(GXFieldTrack *This,
			        G4double valArrIn[GXFieldTrack_ncompSVEC]);

}

#endif 
