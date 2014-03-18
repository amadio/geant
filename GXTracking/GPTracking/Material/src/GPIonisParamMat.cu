#include "GPIonisParamMat.h"
#include "GPElement.h"
#include "GPPhysicalConstants.h"
#include "GPSystemOfUnits.h"

#include "GPConstantArrays.h"

void GPIonisParamMat_Constructor(GPIonisParamMat *This,
				 GPMaterial* material)
{
  This->fMaterial = material;

  This->fMeanExcitationEnergy = 0.0;
  This->fLogMeanExcEnergy = 0.0;

  This->fCdensity = 0.0; 
  This->fMdensity = 0.0; 
  This->fAdensity = 0.0; 
  This->fX0density = 0.0;
  This->fX1density = 0.0;
  This->fD0density = 0.0;

  This->fF1fluct = 0.0;          
  This->fF2fluct = 0.0;                       
  This->fEnergy1fluct = 0.0;
  This->fLogEnergy1fluct = 0.0;
  This->fEnergy2fluct = 0.0;
  This->fLogEnergy2fluct = 0.0;
  This->fEnergy0fluct = 0.0;

  This->twoln10 = 2.*log(10.);

  // compute parameters
  GPIonisParamMat_ComputeParameters(This);
  GPIonisParamMat_ComputeDensityEffect(This);

  // void G4NistMaterialBuilder::NistSimpleMaterials()
  //  GPIonisParamMat_NistSimpleMaterials(This);
}

//G4IonisParamMat::ComputeMeanParameters + G4IonisParamMat::ComputeFluctModel
void GPIonisParamMat_ComputeParameters(GPIonisParamMat *This)
{
  // compute mean excitation energy

  size_t nElements = GPMaterial_GetNumberOfElements(This->fMaterial);
  GPElement* elmVector = GPMaterial_GetElementVector(This->fMaterial);
  G4double* nAtomsPerVolume = 
    GPMaterial_GetVecNbOfAtomsPerVolume(This->fMaterial);
  G4double* fracVector = GPMaterial_GetMassFractionVector(This->fMaterial);
 
  G4double Zeff = 0.;

  for (size_t i=0; i < nElements; i++) {
    GPElement* elm = &(elmVector[i]);
    G4double Z = GPElement_GetZ(elm);
    This->fLogMeanExcEnergy += nAtomsPerVolume[i]*Z
      *log(GPIonisParamMat_GetIonisationPotential(This,Z));
    //      *log(elm->GetIonisation()->GetMeanExcitationEnergy());
    Zeff += fracVector[i]*Z;
  }

  This->fLogMeanExcEnergy /= 
    GPMaterial_GetTotNbOfElectPerVolume(This->fMaterial);
  This->fMeanExcitationEnergy = exp(This->fLogMeanExcEnergy);

  // compute parameters for the energy loss fluctuation model
  // needs an 'effective Z' 

  if (Zeff > 2.) { This->fF2fluct = 2./Zeff; }
  else           { This->fF2fluct = 0.; }

  This->fF1fluct         = 1. - This->fF2fluct;
  This->fEnergy2fluct    = 10.*Zeff*Zeff*eV;
  This->fLogEnergy2fluct = log(This->fEnergy2fluct);
  This->fLogEnergy1fluct = (This->fLogMeanExcEnergy - 
			 This->fF2fluct*This->fLogEnergy2fluct)/This->fF1fluct;
  This->fEnergy1fluct    = exp(This->fLogEnergy1fluct);
  This->fEnergy0fluct    = 10.*eV;
}

//G4IonisParamElm->fMeanExcitationEnergy = 
//G4NistManager::Instance()->GetMeanIonisationEnergy(Z) =
//G4NistMaterialBuilder::GetMeanIonisationEnergy(Z-1) =
//G4NistMaterialBuilder::ionPotentials[Z-1]
G4double GPIonisParamMat_GetIonisationPotential(GPIonisParamMat *This,
						G4double Z) 
{
  G4int index=(G4int)(Z-1);
  G4double res = 10*index;
  if(index >= 0 && index < GPConstantArrays::numElements) { 
    //    res = (This->ionPotentials)[index]; 
    res = GPConstantArrays::IonisationPotential[index]*eV; 
  }
  return res;  
}

FQUALIFIER
G4double  GPIonisParamMat_GetMeanExcitationEnergy(GPIonisParamMat *This)
{
  return This->fMeanExcitationEnergy;
}   

FQUALIFIER
G4double  GPIonisParamMat_GetLogMeanExcEnergy(GPIonisParamMat *This)
{
  return This->fLogMeanExcEnergy;
}

FQUALIFIER
G4double  GPIonisParamMat_GetF1fluct(GPIonisParamMat *This)
{
  return This->fF1fluct;
} 

FQUALIFIER
G4double  GPIonisParamMat_GetF2fluct(GPIonisParamMat *This)
{
  return This->fF2fluct;
}               

FQUALIFIER
G4double  GPIonisParamMat_GetEnergy1fluct(GPIonisParamMat *This)
{
  return This->fEnergy1fluct;
}

FQUALIFIER
G4double  GPIonisParamMat_GetLogEnergy1fluct(GPIonisParamMat *This)
{
  return This->fLogEnergy1fluct;
}

FQUALIFIER
G4double  GPIonisParamMat_GetEnergy2fluct(GPIonisParamMat *This)
{
  return This->fEnergy2fluct;
}

FQUALIFIER
G4double  GPIonisParamMat_GetLogEnergy2fluct(GPIonisParamMat *This)
{
  return This->fLogEnergy2fluct;
}

FQUALIFIER
G4double  GPIonisParamMat_GetEnergy0fluct(GPIonisParamMat *This)
{ 
  return This->fEnergy0fluct;
}

FQUALIFIER
G4double GPIonisParamMat_DensityCorrection(GPIonisParamMat *This,
					   G4double x)
{
  G4double y = 0.0;
  if(x < This->fX0density) {
    if(This->fD0density > 0.0) { 
      y = This->fD0density*pow(10.,2*(x - This->fX0density)); 
    }
  } 
  else if(x >= This->fX1density) { 
    y = This->twoln10*x - This->fCdensity; 
  }
  else {
    y = This->twoln10*x - This->fCdensity 
      + This->fAdensity*pow(This->fX1density - x, This->fMdensity);
  }
  return y;
}

FQUALIFIER
void GPIonisParamMat_ComputeDensityEffect(GPIonisParamMat *This)
{
  //@@@ use the original parameterization
  //@@@ assume that fMaterial is not gas

  //@  G4State State = fMaterial->GetState();

  // Check if density effect data exist in the table
  // R.M. Sternheimer, Atomic Data and Nuclear Data Tables, 30: 261 (1984)
  //@  G4int idx = fDensityData->GetIndex(fMaterial->GetName());

  G4int nelm= GPMaterial_GetNumberOfElements(This->fMaterial);
  //  G4int Z0  = G4int((*(fMaterial->GetElementVector()))[0]->GetZ()+0.5);
  G4int Z0  = G4int(GPElement_GetZ( 
		    &(GPMaterial_GetElementVector(This->fMaterial))[0])+0.5);

  /*
  if(idx < 0 && 1 == nelm) {
    idx = fDensityData->GetElementIndex(Z0, fMaterial->GetState());
  }

  if(idx >= 0) {

    // Take parameters for the density effect correction from
    // R.M. Sternheimer et al. Density Effect For The Ionization Loss 
    // of Charged Particles in Various Substances. 
    // Atom. Data Nucl. Data Tabl. 30 (1984) 261-271. 

    fCdensity = fDensityData->GetCdensity(idx); 
    fMdensity = fDensityData->GetMdensity(idx);
    fAdensity = fDensityData->GetAdensity(idx);
    fX0density = fDensityData->GetX0density(idx);
    fX1density = fDensityData->GetX1density(idx);
    fD0density = fDensityData->GetDelta0density(idx);
    fPlasmaEnergy = fDensityData->GetPlasmaEnergy(idx);
    fAdjustmentFactor = fDensityData->GetAdjustmentFactor(idx);

    // Correction for base material
    const G4Material* bmat = fMaterial->GetBaseMaterial();
    if(bmat) {
      G4double corr = std::log(bmat->GetDensity()/fMaterial->GetDensity());
      fCdensity  += corr;
      fX0density += corr/twoln10;
      fX1density += corr/twoln10;
    }

  } else {
  */

  //  const G4double Cd2 = 4*pi*hbarc_squared*classic_electr_radius;
  const G4double Cd2 = 4*pi*hbarc*hbarc*classic_electr_radius;
  G4double fPlasmaEnergy = 
    sqrt(Cd2*GPMaterial_GetTotNbOfElectPerVolume(This->fMaterial));
  
  // Compute parameters for the density effect correction in DE/Dx formula.
  // The parametrization is from R.M. Sternheimer, Phys. Rev.B,3:3681 (1971)
  G4int icase;
  
  This->fCdensity = 1. + 2*log(This->fMeanExcitationEnergy/fPlasmaEnergy);
  //
  // condensed materials
  //  
  //  if ((State == kStateSolid)||(State == kStateLiquid)) {
    
    const G4double E100eV  = 100.*eV; 
    const G4double ClimiS[2] = {3.681 , 5.215 };
    const G4double X0valS[2] = {1.0   , 1.5   };
    const G4double X1valS[2] = {2.0   , 3.0   };
    
    if(This->fMeanExcitationEnergy < E100eV) { icase = 0; }
    else                               { icase = 1; } 
    
    if(This->fCdensity < ClimiS[icase])    { This->fX0density = 0.2; }
    else { This->fX0density = 0.326*This->fCdensity - X0valS[icase]; }
    
    This->fX1density = X1valS[icase]; 
    This->fMdensity = 3.0;
    
    //special: Hydrogen
    if (1 == nelm && 1 == Z0) {
      This->fX0density = 0.425; 
      This->fX1density = 2.0; 
      This->fMdensity = 5.949;
    }
    //  }
  
  //
  // gases
  //
  /*
    if (State == kStateGas) { 
    
    const G4double ClimiG[] = { 10. , 10.5 , 11. , 11.5 , 12.25 , 13.804};
    const G4double X0valG[] = { 1.6 , 1.7 ,  1.8 ,  1.9 , 2.0   ,  2.0 };
    const G4double X1valG[] = { 4.0 , 4.0 ,  4.0 ,  4.0 , 4.0   ,  5.0 };
    
    icase = 5;
    fX0density = 0.326*fCdensity-2.5 ; fX1density = 5.0 ; fMdensity = 3. ; 
    while((icase > 0)&&(fCdensity < ClimiG[icase])) { icase-- ; }
    fX0density = X0valG[icase]; fX1density = X1valG[icase];
    
    //special: Hydrogen
    if (1 == nelm && 1 == Z0) {
    fX0density = 1.837; fX1density = 3.0; fMdensity = 4.754;
    }
    
    //special: Helium
    if (1 == nelm && 2 == Z0) {
    fX0density = 2.191; fX1density = 3.0; fMdensity = 3.297;
    }
    }
  */
  //  }

  // change parameters if the gas is not in STP.
  // For the correction the density(STP) is needed. 
  // Density(STP) is calculated here : 
  
  /*    
    if (State == kStateGas) { 
    G4double Density  = fMaterial->GetDensity();
    G4double Pressure = fMaterial->GetPressure();
    G4double Temp     = fMaterial->GetTemperature();
      
    G4double DensitySTP = Density*STP_Pressure*Temp/(Pressure*STP_Temperature);

    G4double ParCorr = std::log(Density/DensitySTP);
  
    fCdensity  -= ParCorr;
    fX0density -= ParCorr/twoln10;
    fX1density -= ParCorr/twoln10;
  }
  */
  // fAdensity parameter can be fixed for not conductive materials 
  if(0.0 == This->fD0density) {
    G4double Xa = This->fCdensity/This->twoln10;
    This->fAdensity = This->twoln10*(Xa-This->fX0density)
      /pow((This->fX1density-This->fX0density),This->fMdensity);
  }
}
