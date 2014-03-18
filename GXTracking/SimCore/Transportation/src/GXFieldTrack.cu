#include "GXFieldTrack.h"

FQUALIFIER
void GXFieldTrack_Constructor( GXFieldTrack *This, 
			       GPThreeVector pPosition, 
			       GPThreeVector pMomentum,
			       G4double      kineticEnergy,
			       G4double      restMass_c2,
			       G4double      charge, 
			       G4double      curve_length )
{
  GXFieldTrack_SetCurvePnt(This, pPosition, pMomentum, curve_length );

  This->fKineticEnergy = kineticEnergy;
  This->fRestMass_c2 = restMass_c2;
  This->fCharge = charge;

}

FQUALIFIER
void GXFieldTrack_UpdateFourMomentum( GXFieldTrack *This,
				      G4double kineticEnergy, 
				      GPThreeVector momentum )
{
  GXFieldTrack_SetMomentum(This, momentum ); 
  This->fKineticEnergy= kineticEnergy;
}

FQUALIFIER
void GXFieldTrack_UpdateState(GXFieldTrack *This, 
			      GPThreeVector position, 
			      GPThreeVector momentum,
			      G4double      kineticEnergy )
{ 
  // SetCurvePnt( position, momentumVector, s_curve=0.0);     
  GXFieldTrack_SetPosition(This, position); 
  This->fDistanceAlongCurve= 0.0;

  GXFieldTrack_UpdateFourMomentum(This, kineticEnergy, momentum); 
}

FQUALIFIER
GXFieldTrack& GXFieldTrack_SetCurvePnt(GXFieldTrack *This, 
				       GPThreeVector pPosition, 
				       GPThreeVector pMomentum,  
				       G4double       s_curve )
{
  This->SixVector[0] = GPThreeVector_x(pPosition); 
  This->SixVector[1] = GPThreeVector_y(pPosition); 
  This->SixVector[2] = GPThreeVector_z(pPosition); 

  This->SixVector[3] = GPThreeVector_x(pMomentum); 
  This->SixVector[4] = GPThreeVector_y(pMomentum); 
  This->SixVector[5] = GPThreeVector_z(pMomentum); 

  This->fDistanceAlongCurve= s_curve;

  return *This;
} 

// Assignment operator

FQUALIFIER
GPThreeVector  GXFieldTrack_GetMomentum(const GXFieldTrack *This)
{
  return GPThreeVector_create( This->SixVector[3], 
			       This->SixVector[4], 
			       This->SixVector[5]);

}   

FQUALIFIER
GPThreeVector  GXFieldTrack_GetPosition(const GXFieldTrack *This)
{
  return GPThreeVector_create( This->SixVector[0], 
			       This->SixVector[1], 
			       This->SixVector[2]);
}

FQUALIFIER
G4double       GXFieldTrack_GetCurveLength(const GXFieldTrack *This)
{
  return This->fDistanceAlongCurve;  

}

FQUALIFIER
G4double GXFieldTrack_GetKineticEnergy(GXFieldTrack *This)
{
  return This->fKineticEnergy;
}

// Accessors.

FQUALIFIER
void GXFieldTrack_SetPosition(GXFieldTrack *This, GPThreeVector pPosition) 
{
   This->SixVector[0] = GPThreeVector_x(pPosition); 
   This->SixVector[1] = GPThreeVector_y(pPosition); 
   This->SixVector[2] = GPThreeVector_z(pPosition); 
} 

FQUALIFIER
void GXFieldTrack_SetMomentum(GXFieldTrack *This, GPThreeVector pMomentum )
{
  This->SixVector[3] = GPThreeVector_x(pMomentum); 
  This->SixVector[4] = GPThreeVector_y(pMomentum); 
  This->SixVector[5] = GPThreeVector_z(pMomentum); 
}

FQUALIFIER
void GXFieldTrack_SetCurveLength(GXFieldTrack *This, G4double nCurve_s)
{
  This->fDistanceAlongCurve = nCurve_s;  
}

FQUALIFIER
void GXFieldTrack_SetKineticEnergy(GXFieldTrack *This, G4double newKinEnergy)
{
  This->fKineticEnergy=newKinEnergy;
}

FQUALIFIER
void GXFieldTrack_DumpToArray(GXFieldTrack *This,
			      G4double valArr[GXFieldTrack_ncompSVEC] )
{
  for(G4int i = 0 ; i < GXFieldTrack_ncompSVEC ; i++) {
    valArr[i]=This->SixVector[i];
  }

}

FQUALIFIER
void GXFieldTrack_LoadFromArray(GXFieldTrack *This,
			        G4double valArrIn[GXFieldTrack_ncompSVEC])
{

  for(G4int i = 0 ; i< GXFieldTrack_ncompSVEC ; i++){
     This->SixVector[i] = valArrIn[i];
  }

  GPThreeVector Momentum = 
    GPThreeVector_create(valArrIn[3],valArrIn[4],valArrIn[5]);

  G4double momentum_square= GPThreeVector_mag2(Momentum);

  This->fKineticEnergy = momentum_square / 
    (sqrt(momentum_square+This->fRestMass_c2*This->fRestMass_c2)
     + This->fRestMass_c2 ); 
}
