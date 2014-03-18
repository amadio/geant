#include "GPPhotonModel.h"
#include "GPAtomicShells.h"
#include "GPVParticleChange.h"
#include "GPRandom.h"

FQUALIFIER
GPPhotonModel::GPPhotonModel(curandState* devStates, G4int threadId)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  //Process
  theProcessType = kNullType;

  //GPVEmModel
  lowLimit = 0.1*keV; 
  highLimit = 100.0*TeV;
  eMinActive = 0.0;
  eMaxActive = DBL_MAX;
  polarAngleLimit = pi;
  secondaryThreshold = DBL_MAX;
  theLPMflag = false;
  fParticleChange = 0;
  xSectionTable = 0;
  fCurrentElement = 0;

  //G4PEEffectFluoModel
  fSandiaTable= 0;
}

FQUALIFIER
GPPhotonModel::~GPPhotonModel()
{
}

FQUALIFIER
void GPPhotonModel::SetProcessType(GPEmProcessType aType)
{
  theProcessType = aType;
}

FQUALIFIER 
void GPPhotonModel::SampleSecondaries(GXTrack* track,
				      GPMaterial* material,
				      G4double tmin,
				      G4double tmax)
{
  //pure virtual
  if(theProcessType == kCompton) {
    G4KleinNishinaCompton_SampleSecondaries(track);
  }
  else if (theProcessType == kPhotoElectric) {
    G4PEEffectFluoModel_SampleSecondaries(material,track);    
  }
  else if (theProcessType == kConversion) {
    G4BetheHeitlerModel_SampleSecondaries(material,track);
  }
}

FQUALIFIER
GPVParticleChange* GPPhotonModel::GetParticleChange()
{
   return  fParticleChange;
}

FQUALIFIER
G4double GPPhotonModel::G4VEmModel_CrossSectionPerVolume(GPMaterial* material,
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
      ComputeCrossSectionPerAtom(ekin,Z,emin,emax);
    xsec[i] = cross;
  }
  return cross;
}

FQUALIFIER
G4double GPPhotonModel::CrossSectionPerVolume(GPMaterial* material,
					   G4double ekin,
					   G4double emin,
					   G4double emax)
{
  //virtual implementation
  G4double cross = 0.0;

  if (theProcessType == kPhotoElectric) {
    cross = G4PEEffectFluoModel_CrossSectionPerVolume(material,ekin);
  }
  else {
    cross = G4VEmModel_CrossSectionPerVolume(material,ekin,emin,emax);
  }

  return cross;
}

FQUALIFIER
GPElement* GPPhotonModel::SelectRandomAtom(GPMaterial* material,
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
G4double GPPhotonModel::ComputeCrossSectionPerAtom(G4double kinEnergy, 
						   G4double Z, 
						   G4double, 
						   G4double, 
						   G4double)
{
  //virtual
  G4double cross = 0.0;
  if(theProcessType == kCompton) {
    cross = G4KleinNishinaCompton_ComputeCrossSectionPerAtom(kinEnergy,Z);
  }
  else if (theProcessType == kPhotoElectric) {
    cross = G4PEEffectFluoModel_ComputeCrossSectionPerAtom(kinEnergy,Z);
  }
  else if (theProcessType == kConversion) {
    cross = G4BetheHeitlerModel_ComputeCrossSectionPerAtom(kinEnergy,Z);
  }
  return cross;
}

FQUALIFIER
G4double GPPhotonModel::ChargeSquareRatio(GXTrack* track)
{
  return GetChargeSquareRatio(track->q);
}

FQUALIFIER
G4double GPPhotonModel::GetChargeSquareRatio(G4double q)
{
  return q*q;
}

FQUALIFIER
G4double GPPhotonModel::GetParticleCharge()
{
  return 0;
}

FQUALIFIER
G4double GPPhotonModel::Value(GPMaterial* material,
			   G4double e)
{
  return e*e*CrossSectionPerVolume(material,e,0.0,DBL_MAX);
}

FQUALIFIER
G4double GPPhotonModel::MinPrimaryEnergy()
{
  return 0.0;
}

FQUALIFIER
G4double GPPhotonModel::MaxSecondaryEnergy(G4double kineticEnergy)
{
  return kineticEnergy;
}

FQUALIFIER
void GPPhotonModel::SetupForMaterial(GPMaterial*, 
				  G4double)
{
  ; 
}

FQUALIFIER
void GPPhotonModel::SetParticleChange(GPVParticleChange* p)
{
  if(p && fParticleChange != p) { fParticleChange = p; }
}

FQUALIFIER
void GPPhotonModel::SetCrossSectionTable(GPPhysicsTable* p)
{
  xSectionTable = p;
}

//inline methods

FQUALIFIER 
void GPPhotonModel::SetCurrentElement(GPElement* elm)
{
  fCurrentElement = elm;
}

FQUALIFIER 
GPElement* GPPhotonModel::GetCurrentElement()
{
  return fCurrentElement;
}


FQUALIFIER 
G4double GPPhotonModel::MaxSecondaryKinEnergy(G4double kineticEnergy)
{
  return MaxSecondaryEnergy(kineticEnergy);
}

FQUALIFIER
G4double GPPhotonModel::ComputeDEDXPerVolume(GPMaterial* material,
					  G4double kinEnergy,
					  G4double cutEnergy)
{
  ///@@@virtual implementation in a drived class
  return 0.0;
}

FQUALIFIER
G4double GPPhotonModel::ComputeDEDX(GPMaterial* material,
				 G4double kinEnergy,
				 G4double cutEnergy)
{
  return ComputeDEDXPerVolume(material,kinEnergy,cutEnergy);
}


FQUALIFIER
G4double GPPhotonModel::CrossSection(GPMaterial* material,
				  G4double kinEnergy,
				  G4double cutEnergy,
				  G4double maxEnergy)
{
  return CrossSectionPerVolume(material,kinEnergy,cutEnergy,maxEnergy);
}

FQUALIFIER
G4double GPPhotonModel::ComputeMeanFreePath(G4double ekin,
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
G4double GPPhotonModel::ComputeCrossSectionPerAtom(GPElement* elm,
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
G4double GPPhotonModel::HighEnergyLimit()
{
  return highLimit;
}

FQUALIFIER
G4double GPPhotonModel::LowEnergyLimit()
{
  return lowLimit;
}

FQUALIFIER
G4double GPPhotonModel::HighEnergyActivationLimit()
{
  return eMaxActive;
}

FQUALIFIER
G4double GPPhotonModel::LowEnergyActivationLimit()
{
  return eMinActive;
}

FQUALIFIER
G4double GPPhotonModel::PolarAngleLimit()
{
  return polarAngleLimit;
}

FQUALIFIER
G4double GPPhotonModel::SecondaryThreshold()
{
  return secondaryThreshold;
}

FQUALIFIER 
G4bool GPPhotonModel::LPMFlag() 
{
  return theLPMflag;
}

FQUALIFIER
void GPPhotonModel::SetHighEnergyLimit(G4double val)
{
  highLimit = val;
}

FQUALIFIER
void GPPhotonModel::SetLowEnergyLimit(G4double val)
{
  lowLimit = val;
}

FQUALIFIER
void GPPhotonModel::SetActivationHighEnergyLimit(G4double val)
{
  eMaxActive = val;
}

FQUALIFIER
void GPPhotonModel::SetActivationLowEnergyLimit(G4double val)
{
  eMinActive = val;
}

FQUALIFIER
G4bool GPPhotonModel::IsActive(G4double kinEnergy)
{
  return (kinEnergy >= eMinActive && kinEnergy <= eMaxActive);
}

FQUALIFIER
void GPPhotonModel::SetPolarAngleLimit(G4double val)
{
  polarAngleLimit = val;
}

FQUALIFIER
void GPPhotonModel::SetSecondaryThreshold(G4double val) 
{
  secondaryThreshold = val;
}

FQUALIFIER
void GPPhotonModel::SetLPMFlag(G4bool val) 
{
  theLPMflag = val;
}

FQUALIFIER
GPPhysicsTable* GPPhotonModel::GetCrossSectionTable()
{
  return xSectionTable;
}

//-----------------------------------------------------------------------------
// Utility Functions for Seoncdaries
//-----------------------------------------------------------------------------

GXTrack& GPPhotonModel::GetSecondaryElectron()
{
  //secondary electron from all models
  return theSecondaryElectron;
}

GXTrack& GPPhotonModel::GetSecondaryPositron()
{
  //secondary positron from G4BetheHeitlerModel
  return theSecondaryPositron;
}

FQUALIFIER 
void GPPhotonModel::FillSecondary(GXTrack* secondary,
				  GXTrack* track, 
				  GPThreeVector direction,
				  G4double energy,
				  G4double charge)
{
  G4double p = sqrt(energy*(energy+2.0*electron_mass_c2*charge*charge));
  secondary->x  = track->x;
  secondary->y  = track->y;
  secondary->z  = track->z;
  secondary->s  = 0;
  secondary->px = p*direction.x;
  secondary->py = p*direction.y;
  secondary->pz = p*direction.z;
  secondary->E  = energy;
  secondary->q  = charge;
}

//-----------------------------------------------------------------------------
// G4KleinNishinaCompton for G4ComptonScattering
//-----------------------------------------------------------------------------

G4double GPPhotonModel::G4KleinNishinaCompton_ComputeCrossSectionPerAtom(
                                              G4double GammaEnergy,
                                              G4double Z)
{
  G4double xSection = 0.0 ;
  if ( Z < 0.9999 )                 return xSection;
  if ( GammaEnergy < 0.1*keV      ) return xSection;
  //  if ( GammaEnergy > (100.*GeV/Z) ) return xSection;

  const G4double a = 20.0 , b = 230.0 , c = 440.0;
  
  const G4double
    d1= 2.7965e-1*barn, d2=-1.8300e-1*barn, d3= 6.7527   *barn, d4=-1.9798e+1*barn,
    e1= 1.9756e-5*barn, e2=-1.0205e-2*barn, e3=-7.3913e-2*barn, e4= 2.7079e-2*barn,
    f1=-3.9178e-7*barn, f2= 6.8241e-5*barn, f3= 6.0480e-5*barn, f4= 3.0274e-4*barn;
     
  G4double p1Z = Z*(d1 + e1*Z + f1*Z*Z), p2Z = Z*(d2 + e2*Z + f2*Z*Z),
           p3Z = Z*(d3 + e3*Z + f3*Z*Z), p4Z = Z*(d4 + e4*Z + f4*Z*Z);

  G4double T0  = 15.0*keV; 
  if (Z < 1.5) T0 = 40.0*keV; 

  G4double X   = fmax(GammaEnergy, T0) / electron_mass_c2;
  xSection = p1Z*log(1.+2.*X)/X
               + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
                
  //  modification for low energy. (special case for Hydrogen)
  if (GammaEnergy < T0) {
    G4double dT0 = 1.*keV;
    X = (T0+dT0) / electron_mass_c2 ;
    G4double sigma = p1Z*log(1.+2*X)/X
                    + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
    G4double   c1 = -T0*(sigma-xSection)/(xSection*dT0);             
    G4double   c2 = 0.150; 
    if (Z > 1.5) c2 = 0.375-0.0556*log(Z);
    G4double    y = log(GammaEnergy/T0);
    xSection *= exp(-y*(c1+c2*y));          
  }
  return xSection;
}

void GPPhotonModel::G4KleinNishinaCompton_SampleSecondaries(GXTrack* 
							    aDynamicGamma)
{
  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used 
  // (Nuc Phys 20(1960),15).
  // Note : Effects due to binding of atomic electrons are negliged.
 
  G4double gamEnergy0 = aDynamicGamma->E; //aDynamicGamma->GetKineticEnergy();

  const G4double lowestGammaEnergy = 1.0*eV;

  // extra protection
  if(gamEnergy0 < lowestGammaEnergy) {
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->ProposeLocalEnergyDeposit(gamEnergy0);
    fParticleChange->SetProposedKineticEnergy(0.0);
    return;
  }

  G4double E0_m = gamEnergy0 / electron_mass_c2 ;

  GPThreeVector gamDirection0 = // aDynamicGamma->GetMomentumDirection();
    GPThreeVector_unit(GPThreeVector_create(aDynamicGamma->px,
					    aDynamicGamma->py,
					    aDynamicGamma->pz));
  //
  // sample the energy rate of the scattered gamma 
  //

  G4double epsilon, epsilonsq, onecost, sint2, greject ;

  G4double eps0       = 1./(1. + 2.*E0_m);
  G4double epsilon0sq = eps0*eps0;
  G4double alpha1     = - log(eps0);
  G4double alpha2     = 0.5*(1.- epsilon0sq);

  do {
    if ( alpha1/(alpha1+alpha2) > GPUniformRand(fDevStates,fThreadId) ) {
      epsilon   = exp(-alpha1*GPUniformRand(fDevStates,fThreadId)); // eps0**r
      epsilonsq = epsilon*epsilon; 

    } else {
      epsilonsq = epsilon0sq 
	        + (1.- epsilon0sq)*GPUniformRand(fDevStates,fThreadId);
      epsilon   = sqrt(epsilonsq);
    };

    onecost = (1.- epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);
    greject = 1. - epsilon*sint2/(1.+ epsilonsq);

  } while (greject < GPUniformRand(fDevStates,fThreadId));
  //
  // scattered gamma angles. ( Z - axis along the parent gamma)
  //

  if(sint2 < 0.0) { sint2 = 0.0; }
  G4double cosTeta = 1. - onecost; 
  G4double sinTeta = sqrt (sint2);
  G4double Phi     = twopi * GPUniformRand(fDevStates,fThreadId);

  //
  // update G4VParticleChange for the scattered gamma
  //
   
  GPThreeVector gamDirection1 = 
    GPThreeVector_create(sinTeta*cos(Phi), sinTeta*sin(Phi), cosTeta);

  //  gamDirection1.rotateUz(gamDirection0);
  GPThreeVector_rotateUz(&gamDirection1,gamDirection0);
  G4double gamEnergy1 = epsilon*gamEnergy0;

  if(gamEnergy1 > lowestGammaEnergy) {
    fParticleChange->SetProposedMomentumDirection(gamDirection1);
    fParticleChange->SetProposedKineticEnergy(gamEnergy1);
  } else { 
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->ProposeLocalEnergyDeposit(gamEnergy1);
    fParticleChange->SetProposedKineticEnergy(0.0);
  }

  //
  // kinematic of the scattered electron
  //

  G4double eKinEnergy = gamEnergy0 - gamEnergy1;

  if(eKinEnergy > DBL_MIN) {
    GPThreeVector eDirection = GPThreeVector_unit(
             GPThreeVector_sub(GPThreeVector_mult(gamDirection0,gamEnergy0), 
			       GPThreeVector_mult(gamDirection1,gamEnergy1)));

    // create G4DynamicParticle object for the electron.
    // G4DynamicParticle* dp = 
    // new G4DynamicParticle(theElectron,eDirection,eKinEnergy);
    // fvect->push_back(dp);
    FillSecondary(&theSecondaryElectron,aDynamicGamma,eDirection,eKinEnergy,-1);
  }
}

//-----------------------------------------------------------------------------
// G4PEEffectFluoModel for G4PhotoElectricEffect
//-----------------------------------------------------------------------------

FQUALIFIER
G4double GPPhotonModel::G4PEEffectFluoModel_ComputeCrossSectionPerAtom(
						   G4double energy,
						   G4double Z)
{
  G4double* SandiaCof = GPSandiaTable_GetSandiaCofPerAtom(fSandiaTable,
							  (G4int)Z, energy);

  G4double energy2 = energy*energy;
  G4double energy3 = energy*energy2;
  G4double energy4 = energy2*energy2;

  return SandiaCof[0]/energy  + SandiaCof[1]/energy2 +
    SandiaCof[2]/energy3 + SandiaCof[3]/energy4;
}

FQUALIFIER
G4double GPPhotonModel::G4PEEffectFluoModel_CrossSectionPerVolume(
                                                 GPMaterial* material,
					         G4double energy)
{
  //  G4double* SandiaCof = 
  //    material->GetSandiaTable()->GetSandiaCofForMaterial(energy);
  //@@@ this will return an array of zeros - this is not used for now
  G4double* SandiaCof = GPSandiaTable_GetSandiaCofForMaterial(fSandiaTable,
							      energy);
  G4double energy2 = energy*energy;
  G4double energy3 = energy*energy2;
  G4double energy4 = energy2*energy2;
          
  return SandiaCof[0]/energy  + SandiaCof[1]/energy2 +
    SandiaCof[2]/energy3 + SandiaCof[3]/energy4; 
}

FQUALIFIER
void 
GPPhotonModel::G4PEEffectFluoModel_SampleSecondaries(GPMaterial* aMaterial,
						     GXTrack* aDynamicPhoton)
{
  //  const G4Material* aMaterial = couple->GetMaterial();

  G4double energy = aDynamicPhoton->E; //aDynamicPhoton->GetKineticEnergy();

  // select randomly one element constituing the material.
  GPElement* anElement = SelectRandomAtom(aMaterial,energy);
  
  //
  // Photo electron
  //

  // Select atomic shell
  G4int Z = G4int(GPElement_GetZ(anElement));

  //  G4int nShells = anElement->GetNbOfAtomicShells();
  G4int nShells = GPAtomicShells_GetNumberOfShells(Z);

  G4int i = 0;  
  for(; i<nShells; ++i) {
    //    if(energy >= anElement->GetAtomicShell(i)) { break; }
    if(energy >= GPAtomicShells_GetBindingEnergy(Z,i)) { break; }
  }

  G4double edep = energy;

  // Normally one shell is available 
  if (i < nShells) { 
  
    //    G4double bindingEnergy  = anElement->GetAtomicShell(i);
    G4double bindingEnergy  = GPAtomicShells_GetBindingEnergy(Z,i);
    G4double elecKineEnergy = energy - bindingEnergy;

    // create photo electron
    //
    G4double fminimalEnergy = 1.0*eV;
    if (elecKineEnergy > fminimalEnergy) {
      edep = bindingEnergy;
      GPThreeVector elecDirection = 
        SauterGavrilaAngularDistribution_SampleDirection(aDynamicPhoton);
      //    GetAngularDistribution()->SampleDirection(aDynamicPhoton, 
      //                                              elecKineEnergy,
      //                                              i, 
      //                                              couple->GetMaterial());
      //    G4DynamicParticle* aParticle = 
      //      new G4DynamicParticle(theElectron, elecDirection, elecKineEnergy);
      //      fvect->push_back(aParticle);
      FillSecondary(&theSecondaryElectron,aDynamicPhoton,elecDirection,
		    elecKineEnergy,-1);
      
    } 

    // sample deexcitation
    //
    /*
    if(fAtomDeexcitation) {
      G4int index = couple->GetIndex();
      if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
        G4int Z = G4lrint(anElement->GetZ());
        G4AtomicShellEnumerator as = G4AtomicShellEnumerator(i);
        const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
        size_t nbefore = fvect->size();
        fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
        size_t nafter = fvect->size();
        if(nafter > nbefore) {
          for (size_t j=nbefore; j<nafter; ++j) {
            edep -= ((*fvect)[j])->GetKineticEnergy();
          } 
        }
      }
    }
    */
  }
  // kill primary photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  if(edep > 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

FQUALIFIER GPThreeVector
GPPhotonModel::SauterGavrilaAngularDistribution_SampleDirection(GXTrack*
								aDynamicGamma)
{
  //  G4double tau = dp->GetKineticEnergy()/electron_mass_c2;
  G4double tau = aDynamicGamma->E/electron_mass_c2;
  const G4double taulimit = 30.0;

  GPThreeVector fLocalDirection;

  if (tau > taulimit) {
    fLocalDirection = // dp->GetMomentumDirection(); 
      GPThreeVector_unit(GPThreeVector_create(aDynamicGamma->px,
					      aDynamicGamma->py,
					      aDynamicGamma->pz));

    // Bugzilla 1120
    // SI on 05/09/2010 as suggested by JG 04/09/10 
  } else {
 
    G4double invgamma  = 1.0/(tau + 1.0);
    G4double beta      = sqrt(tau*(tau + 2.0))*invgamma;
    G4double b         = 0.5*tau*(tau*tau - 1.0);
    G4double invgamma2 = invgamma*invgamma;
   
    G4double rndm,term,greject,grejsup,costeta,sint2;
    if (tau < 1.) { grejsup = (1.+b-beta*b)/invgamma2; }
    else          { grejsup = (1.+b+beta*b)/invgamma2; }

    do { 
      rndm    = 1 - 2*GPUniformRand(fDevStates,fThreadId);
      costeta = (rndm + beta)/(rndm*beta + 1);
      term    = invgamma2/(1 + beta*rndm);
      sint2   = (1 - costeta)*(1 + costeta);
      greject = sint2*(1 + b*term)/(term*term);

    } while(greject < GPUniformRand(fDevStates,fThreadId)*grejsup);
       
    G4double sinteta = sqrt(sint2);
    G4double phi  = twopi*GPUniformRand(fDevStates,fThreadId); 

    //    fLocalDirection.set(sinteta*cos(phi), sinteta*sin(phi), costeta);
    fLocalDirection = GPThreeVector_create(sinteta*cos(phi), 
					   sinteta*sin(phi),
					   costeta);

    //    fLocalDirection.rotateUz(dp->GetMomentumDirection());
    GPThreeVector refDirection = // dp->GetMomentumDirection(); 
      GPThreeVector_unit(GPThreeVector_create(aDynamicGamma->px,
					      aDynamicGamma->py,
					      aDynamicGamma->pz));
    GPThreeVector_rotateUz(&fLocalDirection,refDirection);
  }
  return fLocalDirection;
}

FQUALIFIER
void GPPhotonModel::SetaSandiaTable(GPSandiaTable* table)
{
   fSandiaTable = table;
}

//-----------------------------------------------------------------------------
// G4BetheHeitlerModel for G4GammaConversion
//-----------------------------------------------------------------------------

FQUALIFIER
G4double 
GPPhotonModel::G4BetheHeitlerModel_ComputeCrossSectionPerAtom(
                                                G4double GammaEnergy, 
						G4double Z)
// Calculates the microscopic cross section in GEANT4 internal units.
// A parametrized formula from L. Urban is used to estimate
// the total cross section.
// It gives a good description of the data from 1.5 MeV to 100 GeV.
// below 1.5 MeV: sigma=sigma(1.5MeV)*(GammaEnergy-2electronmass)
//                                   *(GammaEnergy-2electronmass) 
{
  //  static 
  const G4double GammaEnergyLimit = 1.5*MeV;
  G4double xSection = 0.0 ;
  if ( Z < 0.9 || GammaEnergy <= 2.0*electron_mass_c2 ) { return xSection; }

  //  static 
  const G4double
    a0= 8.7842e+2*microbarn, a1=-1.9625e+3*microbarn, a2= 1.2949e+3*microbarn,
    a3=-2.0028e+2*microbarn, a4= 1.2575e+1*microbarn, a5=-2.8333e-1*microbarn;

  //  static 
  const G4double
    b0=-1.0342e+1*microbarn, b1= 1.7692e+1*microbarn, b2=-8.2381   *microbarn,
    b3= 1.3063   *microbarn, b4=-9.0815e-2*microbarn, b5= 2.3586e-3*microbarn;

  //  static 
  const G4double
    c0=-4.5263e+2*microbarn, c1= 1.1161e+3*microbarn, c2=-8.6749e+2*microbarn,
    c3= 2.1773e+2*microbarn, c4=-2.0467e+1*microbarn, c5= 6.5372e-1*microbarn;

  G4double GammaEnergySave = GammaEnergy;
  if (GammaEnergy < GammaEnergyLimit) { GammaEnergy = GammaEnergyLimit; }

  G4double X=log(GammaEnergy/electron_mass_c2), X2=X*X, X3=X2*X, X4=X3*X, X5=X4*X;

  G4double F1 = a0 + a1*X + a2*X2 + a3*X3 + a4*X4 + a5*X5,
           F2 = b0 + b1*X + b2*X2 + b3*X3 + b4*X4 + b5*X5,
           F3 = c0 + c1*X + c2*X2 + c3*X3 + c4*X4 + c5*X5;     

  xSection = (Z + 1.)*(F1*Z + F2*Z*Z + F3);

  if (GammaEnergySave < GammaEnergyLimit) {

    X = (GammaEnergySave  - 2.*electron_mass_c2)
      / (GammaEnergyLimit - 2.*electron_mass_c2);
    xSection *= X*X;
  }

  if (xSection < 0.) { xSection = 0.; }
  return xSection;
}

FQUALIFIER
void GPPhotonModel::G4BetheHeitlerModel_SampleSecondaries(
                                            GPMaterial* aMaterial,
                                            GXTrack* aDynamicGamma)
// The secondaries e+e- energies are sampled using the Bethe - Heitler
// cross sections with Coulomb correction.
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : Effects due to the breakdown of the Born approximation at
//          low energy are ignored.
// Note 2 : The differential cross section implicitly takes account of 
//          pair creation in both nuclear and atomic electron fields.
//          However triplet prodution is not generated.
{
  //  const G4Material* aMaterial = couple->GetMaterial();

  G4double GammaEnergy = aDynamicGamma->E;//aDynamicGamma->GetKineticEnergy();

  //  G4ParticleMomentum 
  GPThreeVector GammaDirection = // aDynamicGamma->GetMomentumDirection();
    GPThreeVector_unit(GPThreeVector_create(aDynamicGamma->px,
                                            aDynamicGamma->py,
                                            aDynamicGamma->pz));

  G4double epsil ;
  G4double epsil0 = electron_mass_c2/GammaEnergy ;
  if(epsil0 > 1.0) { return; }

  // do it fast if GammaEnergy < 2. MeV
  //  static 
  const G4double Egsmall=2.*MeV;

  // select randomly one element constituing the material
  //  G4Element* anElement = SelectRandomAtom(aMaterial, theGamma, GammaEnergy);
  GPElement* anElement = SelectRandomAtom(aMaterial, GammaEnergy);

  if (GammaEnergy < Egsmall) {

    epsil = epsil0 + (0.5-epsil0)*GPUniformRand(fDevStates,fThreadId);

  } else {
    // now comes the case with GammaEnergy >= 2. MeV

    // Extract Coulomb factor for this Element
    G4int logZ3 = log(1.0*G4int(GPElement_GetZ(anElement) + 0.5))/3.0; 
    G4double FZ = 8.*logZ3 ; //8.*(anElement->GetIonisation()->GetlogZ3());

    if (GammaEnergy > 50.*MeV) { FZ += 8.*(GPElement_GetfCoulomb(anElement)); }

    // limits of the screening variable
    G4int Z3 = pow(1.0*G4int(GPElement_GetZ(anElement) + 0.5),1/3.0);
    G4double screenfac = 136.*epsil0/Z3; //(anElement->GetIonisation()->GetZ3());
    G4double screenmax = exp ((42.24 - FZ)/8.368) - 0.952 ;
    G4double screenmin = fmin(4.*screenfac,screenmax);

    // limits of the energy sampling
    G4double epsil1 = 0.5 - 0.5*sqrt(1. - screenmin/screenmax) ;
    G4double epsilmin = fmax(epsil0,epsil1) , epsilrange = 0.5 - epsilmin;

    //
    // sample the energy rate of the created electron (or positron)

    //G4double epsil, screenvar, greject ;
    G4double  screenvar, greject ;

    G4double F10 = G4BetheHeitlerModel_ScreenFunction1(screenmin) - FZ;
    G4double F20 = G4BetheHeitlerModel_ScreenFunction2(screenmin) - FZ;

    G4double NormF1 = fmax(F10*epsilrange*epsilrange,0.); 
    G4double NormF2 = fmax(1.5*F20,0.);

    do {
      if ( NormF1/(NormF1+NormF2) > GPUniformRand(fDevStates,fThreadId) ) {
        epsil = 0.5 - 
	  epsilrange*pow(GPUniformRand(fDevStates,fThreadId), 0.333333);
        screenvar = screenfac/(epsil*(1-epsil));
	greject = (G4BetheHeitlerModel_ScreenFunction1(screenvar) - FZ)/F10;
              
      } else { 
        epsil = epsilmin + epsilrange*GPUniformRand(fDevStates,fThreadId);
        screenvar = screenfac/(epsil*(1-epsil));
	greject = (G4BetheHeitlerModel_ScreenFunction2(screenvar) - FZ)/F20;
      }

    } while( greject < GPUniformRand(fDevStates,fThreadId) );
  }   //  end of epsil sampling
   
  //
  // fixe charges randomly
  //

  G4double ElectTotEnergy, PositTotEnergy;
  if (GPUniformRand(fDevStates,fThreadId) > 0.5) {

    ElectTotEnergy = (1.-epsil)*GammaEnergy;
    PositTotEnergy = epsil*GammaEnergy;
     
  } else {
    
    PositTotEnergy = (1.-epsil)*GammaEnergy;
    ElectTotEnergy = epsil*GammaEnergy;
  }

  //
  // scattered electron (positron) angles. ( Z - axis along the parent photon)
  //
  //  universal distribution suggested by L. Urban 
  // (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))
  G4double u;
  const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  if (9./(9.+d) >GPUniformRand(fDevStates,fThreadId)) {
    u= - log(GPUniformRand(fDevStates,fThreadId)*
	     GPUniformRand(fDevStates,fThreadId))/a1;
  }
  else {                            
    u= - log(GPUniformRand(fDevStates,fThreadId)*
	     GPUniformRand(fDevStates,fThreadId))/a2;
  }

  G4double TetEl = u*electron_mass_c2/ElectTotEnergy;
  G4double TetPo = u*electron_mass_c2/PositTotEnergy;
  G4double Phi  = twopi * GPUniformRand(fDevStates,fThreadId);
  G4double dxEl= sin(TetEl)*cos(Phi),dyEl= sin(TetEl)*sin(Phi),dzEl=cos(TetEl);
  G4double dxPo=-sin(TetPo)*cos(Phi),dyPo=-sin(TetPo)*sin(Phi),dzPo=cos(TetPo);

  //
  // kinematic of the created pair
  //
  // the electron and positron are assumed to have a symetric
  // angular distribution with respect to the Z axis along the parent photon.

  G4double ElectKineEnergy = fmax(0.,ElectTotEnergy - electron_mass_c2);

  GPThreeVector ElectDirection = GPThreeVector_create(dxEl, dyEl, dzEl);
  GPThreeVector_rotateUz(&ElectDirection,GammaDirection);   

  // create G4DynamicParticle object for the particle1  
  //  G4DynamicParticle* aParticle1= new G4DynamicParticle(
  //                     theElectron,ElectDirection,ElectKineEnergy);
  
  // the e+ is always created (even with Ekine=0) for further annihilation.

  G4double PositKineEnergy = fmax(0.,PositTotEnergy - electron_mass_c2);

  GPThreeVector PositDirection = GPThreeVector_create(dxPo, dyPo, dzPo);
  GPThreeVector_rotateUz(&PositDirection,GammaDirection);   

  // create G4DynamicParticle object for the particle2 
  //  G4DynamicParticle* aParticle2= new G4DynamicParticle(
  //                      thePositron,PositDirection,PositKineEnergy);

  // Fill output vector
  //  fvect->push_back(aParticle1);
  //  fvect->push_back(aParticle2);
  FillSecondary(&theSecondaryElectron,aDynamicGamma,ElectDirection,
		ElectKineEnergy,-1);
  FillSecondary(&theSecondaryPositron,aDynamicGamma,PositDirection,
		PositKineEnergy,1);

  // kill incident photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);   
}

FQUALIFIER G4double 
GPPhotonModel::G4BetheHeitlerModel_ScreenFunction1(G4double ScreenVariable)
{
  // compute the value of the screening function 3*PHI1 - PHI2
  G4double screenVal;
  
  if (ScreenVariable > 1.)
    screenVal = 42.24 - 8.368*log(ScreenVariable+0.952);
  else
    screenVal = 42.392 - ScreenVariable*(7.796 - 1.961*ScreenVariable);
  
  return screenVal;
}

FQUALIFIER G4double 
GPPhotonModel::G4BetheHeitlerModel_ScreenFunction2(G4double ScreenVariable)
{
  // compute the value of the screening function 1.5*PHI1 - 0.5*PHI2
  G4double screenVal;
  
  if (ScreenVariable > 1.)
    screenVal = 42.24 - 8.368*log(ScreenVariable+0.952);
  else
    screenVal = 41.405 - ScreenVariable*(5.828 - 0.8945*ScreenVariable);

  return screenVal;
}
