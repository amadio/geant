#ifndef GPIonisParamMat_HH
#define GPIonisParamMat_HH

#include "GPTypeDef.h"
#include "GPMaterial.h"

struct GPIonisParamMat
{
  GPMaterial* fMaterial;                    // this material

  // parameters for mean energy loss calculation
  G4double  fMeanExcitationEnergy;          // 
  G4double  fLogMeanExcEnergy;              // 

  // parameters of the density correction
  G4double fCdensity;                      // mat.constant
  G4double fMdensity;                      // exponent
  G4double fAdensity;                      //
  G4double fX0density;                     //
  G4double fX1density;                     //
  G4double fD0density;

  // parameters of the energy loss fluctuation model
  G4double fF1fluct;                       
  G4double fF2fluct;                       
  G4double fEnergy1fluct;
  G4double fLogEnergy1fluct;
  G4double fEnergy2fluct;
  G4double fLogEnergy2fluct;
  G4double fEnergy0fluct;

  G4double twoln10;
};

extern "C" {

  FQUALIFIER 
  void GPIonisParamMat_Constructor(GPIonisParamMat *This,
				   GPMaterial* material); 
  FQUALIFIER    
  void GPIonisParamMat_ComputeParameters(GPIonisParamMat *This);

  FQUALIFIER
  G4double GPIonisParamMat_GetIonisationPotential(GPIonisParamMat *This,
						  G4double Z);   

  //get  parameters
  FQUALIFIER
  G4double  GPIonisParamMat_GetMeanExcitationEnergy(GPIonisParamMat *This);   

  FQUALIFIER
  G4double  GPIonisParamMat_GetLogMeanExcEnergy(GPIonisParamMat *This); 

  FQUALIFIER
  G4double  GPIonisParamMat_GetF1fluct(GPIonisParamMat *This);      

  FQUALIFIER
  G4double  GPIonisParamMat_GetF2fluct(GPIonisParamMat *This);               

  FQUALIFIER
  G4double  GPIonisParamMat_GetEnergy1fluct(GPIonisParamMat *This);           

  FQUALIFIER
  G4double  GPIonisParamMat_GetLogEnergy1fluct(GPIonisParamMat *This);        

  FQUALIFIER
  G4double  GPIonisParamMat_GetEnergy2fluct(GPIonisParamMat *This);           

  FQUALIFIER
  G4double  GPIonisParamMat_GetLogEnergy2fluct(GPIonisParamMat *This);        

  FQUALIFIER
  G4double  GPIonisParamMat_GetEnergy0fluct(GPIonisParamMat *This);           

  //density correction
  FQUALIFIER
  void GPIonisParamMat_ComputeDensityEffect(GPIonisParamMat *This);

  FQUALIFIER
  G4double GPIonisParamMat_DensityCorrection(GPIonisParamMat *This,
					     G4double x);

}

#endif
