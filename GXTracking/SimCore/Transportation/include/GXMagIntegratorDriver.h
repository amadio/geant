#ifndef GXMagInt_Driver_HH
#define GXMagInt_Driver_HH

#include "GPTypeDef.h"
#include "GXFieldTrack.h"
#include "GXClassicalRK4.h"

CONSTTYPE const G4double max_stepping_increase = 5.0;
CONSTTYPE const G4double max_stepping_decrease = 0.1;

struct GXMagInt_Driver
{
  G4double  fMinimumStep;
  GXClassicalRK4 *pIntStepper;

  G4double safety;
  G4double pshrnk;   //  exponent for shrinking
  G4double pgrow;    //  exponent for growth
  G4double errcon;

  G4double  fSmallestFraction;  //   Expected range 1e-12 to 5e-15;  
  G4int   fMaxNoSteps;
};

extern "C" {

FQUALIFIER
void GXMagInt_Driver_Constructor(GXMagInt_Driver *This, 
                                 G4double         hminimum, 
                                 GXClassicalRK4  *pStepper);

FQUALIFIER
G4bool GXMagInt_Driver_AccurateAdvance( GXMagInt_Driver *This,
                                        GXFieldTrack& y_current,
                                        G4double     hstep,
                                        G4double     eps,
                                        G4double hinitial );

FQUALIFIER
G4bool GXMagInt_Driver_AccurateAdvance2( GXMagInt_Driver *This,
					 GXFieldTrack& y_current,
					 G4double     hstep,
					 G4double     eps,
					 G4double hinitial );

FQUALIFIER
void GXMagInt_Driver_OneGoodStep( GXMagInt_Driver *This,
                                  G4double y[],      
                                  const G4double dydx[],
                                  G4double& x,       
                                  G4double htry,
                                  G4double eps_rel_max,
                                  G4double& hdid,    
                                  G4double& hnext ); 

FQUALIFIER
void GXMagInt_Driver_OneGoodStep2( GXMagInt_Driver *This,
                                  G4double y[],      
                                  const G4double dydx[],
                                  G4double& x,       
                                  G4double htry,
                                  G4double eps_rel_max,
                                  G4double& hdid,    
                                  G4double& hnext ); 

FQUALIFIER
G4bool  GXMagInt_Driver_QuickAdvance(GXMagInt_Driver *This,        
                                     GXFieldTrack& y_posvel,
                                     const G4double     dydx[],  
                                     G4double     hstep, 
                                     G4double&    dchord_step,
                                     G4double&    dyerr );
FQUALIFIER
G4double GXMagInt_Driver_ComputeNewStepSize(GXMagInt_Driver *This, 
					    G4double  errMaxNorm, 
					    G4double  hstepCurrent);

FQUALIFIER
G4double GXMagInt_Driver_GetHmin( GXMagInt_Driver *This );

FQUALIFIER
G4double GXMagInt_Driver_GetSafety( GXMagInt_Driver *This );

FQUALIFIER
G4double GXMagInt_Driver_GetPshrnk( GXMagInt_Driver *This );

FQUALIFIER
G4double GXMagInt_Driver_GetPgrow( GXMagInt_Driver *This );

FQUALIFIER
G4double GXMagInt_Driver_GetErrcon( GXMagInt_Driver *This );

FQUALIFIER
const GXClassicalRK4* GXMagInt_Driver_GetStepper(GXMagInt_Driver *This);

FQUALIFIER
void GXMagInt_Driver_GetDerivatives( GXMagInt_Driver *This,
                                     GXFieldTrack &y_curr,
                                     G4double     dydx[]);
}

#endif
