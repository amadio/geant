#include "GPVEmModel.h"
#include "GPVParticleChange.h"
//#include "G4ParticleChangeForGamma.hh"
#include "GPRandom.h"

FQUALIFIER
GPVEmModel::GPVEmModel(curandState* devStates, int threadId)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  lowLimit = 0.1*keV; 
  highLimit = 100.0*TeV;
  eMinActive = 0.0;
  eMaxActive = DBL_MAX;
  polarAngleLimit = pi;
  secondaryThreshold = DBL_MAX;
  theLPMflag = false;
  pParticleChange = 0;
  xSectionTable = 0;
  fCurrentElement = 0;

}

FQUALIFIER
GPVEmModel::~GPVEmModel()
{
  ;
}

FQUALIFIER 
void GPVEmModel::SampleSecondaries(GXTrack* track,
				   GPMaterial* material,
				   G4double tmin,
				   G4double tmax)
{
  //pure virtual
  ; 
}

GXTrack& GPVEmModel::GetSecondary()
{
  //secondary electron
  return theSecondaryElectron;
}


FQUALIFIER
GPVParticleChange* GPVEmModel::GetParticleChange()
{
   return  pParticleChange;
}

FQUALIFIER
G4double GPVEmModel::G4VEmModel_CrossSectionPerVolume(GPMaterial* material,
						      G4double ekin,
						      G4double emin,
						      G4double emax)
{
  SetupForMaterial(material, ekin);
  G4double cross = 0.0;
  GPElement* theElementVector = GPMaterial_GetElementVector(material);
  G4double* theAtomNumDensityVector = 
    GPMaterial_GetVecNbOfAtomsPerVolume(material);
  G4int nelm = GPMaterial_GetNumberOfElements(material); 

  for (G4int i=0; i<nelm; i++) {
    G4double Z = GPElement_GetZ(&(theElementVector[i]));
    cross += theAtomNumDensityVector[i]*
      ComputeCrossSectionPerAtom(Z,ekin,emin,emax);
    xsec[i] = cross;
  }
  return cross;
}

FQUALIFIER
G4double GPVEmModel::CrossSectionPerVolume(GPMaterial* material,
					   G4double ekin,
					   G4double emin,
					   G4double emax)
{
  //@@@virtual implementation in a derived class
  G4double cross = 0.0;
  return cross;
}

FQUALIFIER
GPElement* GPVEmModel::SelectRandomAtom(GPMaterial* material,
					G4double kinEnergy,
					G4double tcut,
					G4double tmax)
{
  GPElement* theElementVector = GPMaterial_GetElementVector(material);
  G4int n = GPMaterial_GetNumberOfElements(material) - 1;
  fCurrentElement = &(theElementVector[n]);

  if (n > 0) {
    G4double x = GPUniformRand(fDevStates,fThreadId)*
                 G4VEmModel_CrossSectionPerVolume(material,kinEnergy,
						  tcut,tmax);
    for(G4int i=0; i<n; ++i) {
      if (x <= xsec[i]) {
	fCurrentElement = &(theElementVector[i]);
	break;
      }
    }
  }
  return fCurrentElement;
}

FQUALIFIER
G4double GPVEmModel::ComputeCrossSectionPerAtom(G4double, G4double, G4double,
						G4double, G4double)
{
  return 0.0;
}

FQUALIFIER
G4double GPVEmModel::ChargeSquareRatio(GXTrack* track)
{
  return GetChargeSquareRatio(track->q);
}

FQUALIFIER
G4double GPVEmModel::GetChargeSquareRatio(G4double q)
{
  return q*q;
}

FQUALIFIER
G4double GPVEmModel::GetParticleCharge()
{
  return 0;
}

FQUALIFIER
G4double GPVEmModel::Value(GPMaterial* material,
			   G4double e)
{
  return e*e*CrossSectionPerVolume(material,e,0.0,DBL_MAX);
}

FQUALIFIER
G4double GPVEmModel::MinPrimaryEnergy()
{
  return 0.0;
}

FQUALIFIER
G4double GPVEmModel::MaxSecondaryEnergy(G4double kineticEnergy)
{
  return kineticEnergy;
}

FQUALIFIER
void GPVEmModel::SetupForMaterial(GPMaterial*, 
				  G4double)
{
  ; 
}

FQUALIFIER
void GPVEmModel::SetParticleChange(GPVParticleChange* p)
{
  if(p && pParticleChange != p) { pParticleChange = p; }
}

FQUALIFIER
void GPVEmModel::SetCrossSectionTable(GPPhysicsTable* p)
{
  xSectionTable = p;
}

//inline methods

FQUALIFIER 
void GPVEmModel::SetCurrentElement(GPElement* elm)
{
  fCurrentElement = elm;
}

FQUALIFIER 
GPElement* GPVEmModel::GetCurrentElement()
{
  return fCurrentElement;
}


FQUALIFIER 
G4double GPVEmModel::MaxSecondaryKinEnergy(G4double kineticEnergy)
{
  return MaxSecondaryEnergy(kineticEnergy);
}

FQUALIFIER
G4double GPVEmModel::ComputeDEDXPerVolume(GPMaterial* material,
					  G4double kinEnergy,
					  G4double cutEnergy)
{
  ///@@@virtual implementation in a drived class
  return 0.0;
}

FQUALIFIER
G4double GPVEmModel::ComputeDEDX(GPMaterial* material,
				 G4double kinEnergy,
				 G4double cutEnergy)
{
  return ComputeDEDXPerVolume(material,kinEnergy,cutEnergy);
}


FQUALIFIER
G4double GPVEmModel::CrossSection(GPMaterial* material,
				  G4double kinEnergy,
				  G4double cutEnergy,
				  G4double maxEnergy)
{
  return CrossSectionPerVolume(material,kinEnergy,cutEnergy,maxEnergy);
}

FQUALIFIER
G4double GPVEmModel::ComputeMeanFreePath(G4double ekin,
					 GPMaterial* material,
					 G4double emin,
					 G4double emax)
{
  G4double mfp = DBL_MAX;
  G4double cross = CrossSectionPerVolume(material,ekin,emin,emax);
  if (cross > DBL_MIN) { mfp = 1./cross; }
  return mfp;
}

FQUALIFIER
G4double GPVEmModel::ComputeCrossSectionPerAtom(GPElement* elm,
						G4double kinEnergy, 
						G4double cutEnergy,
						G4double maxEnergy)
{
  fCurrentElement = elm;
  return ComputeCrossSectionPerAtom(kinEnergy,
				    GPElement_GetZ(elm),
				    GPElement_GetN(elm),
				    cutEnergy,
				    maxEnergy);
}

FQUALIFIER
G4double GPVEmModel::HighEnergyLimit()
{
  return highLimit;
}

FQUALIFIER
G4double GPVEmModel::LowEnergyLimit()
{
  return lowLimit;
}

FQUALIFIER
G4double GPVEmModel::HighEnergyActivationLimit()
{
  return eMaxActive;
}

FQUALIFIER
G4double GPVEmModel::LowEnergyActivationLimit()
{
  return eMinActive;
}

FQUALIFIER
G4double GPVEmModel::PolarAngleLimit()
{
  return polarAngleLimit;
}

FQUALIFIER
G4double GPVEmModel::SecondaryThreshold()
{
  return secondaryThreshold;
}

FQUALIFIER 
G4bool GPVEmModel::LPMFlag() 
{
  return theLPMflag;
}

FQUALIFIER
void GPVEmModel::SetHighEnergyLimit(G4double val)
{
  highLimit = val;
}

FQUALIFIER
void GPVEmModel::SetLowEnergyLimit(G4double val)
{
  lowLimit = val;
}

FQUALIFIER
void GPVEmModel::SetActivationHighEnergyLimit(G4double val)
{
  eMaxActive = val;
}

FQUALIFIER
void GPVEmModel::SetActivationLowEnergyLimit(G4double val)
{
  eMinActive = val;
}

FQUALIFIER
G4bool GPVEmModel::IsActive(G4double kinEnergy)
{
  return (kinEnergy >= eMinActive && kinEnergy <= eMaxActive);
}

FQUALIFIER
void GPVEmModel::SetPolarAngleLimit(G4double val)
{
  polarAngleLimit = val;
}

FQUALIFIER
void GPVEmModel::SetSecondaryThreshold(G4double val) 
{
  secondaryThreshold = val;
}

FQUALIFIER
void GPVEmModel::SetLPMFlag(G4bool val) 
{
  theLPMflag = val;
}

FQUALIFIER
GPPhysicsTable* GPVEmModel::GetCrossSectionTable()
{
  return xSectionTable;
}
