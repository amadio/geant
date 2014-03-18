#include "GPMollerBhabhaModel.h"
#include "GPPhysicalConstants.h"
#include "GPSystemOfUnits.h"
//#include "G4Electron.hh"
//#include "G4Positron.hh"
//#include "Randomize.hh"
//#include "G4ParticleChangeForLoss.hh"

#include "GPIonisParamMat.h"

#include "GPRandom.h"

//using namespace std;

#include "stdio.h"
// #ifdef GPUPLOTS
// #ifndef __CUDA_ARCH__
//   #include "GPHistoManager.hh"
// #endif
// #endif

FQUALIFIER
GPMollerBhabhaModel::GPMollerBhabhaModel(curandState* devStates, int threadId)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;

  isElectron = true;
  twoln10 = 2.0*log(10.0);

  //secondary from GPMollerBhabhaModel
  //  theDeltaRay = (GXTrack*) malloc(sizeof(GXTrack));

  //G4VEmModel
  pParticleChange = 0;

  lowLimit = 0.02*keV;
  highLimit = 100.0*TeV;
  eMinActive = 0.0;
  eMaxActive = DBL_MAX;
}

FQUALIFIER
GPMollerBhabhaModel::~GPMollerBhabhaModel()
{}

FQUALIFIER
G4double GPMollerBhabhaModel::MaxSecondaryEnergy(G4double kinEnergy) 
{
  G4double tmax = kinEnergy;
  if(isElectron) { tmax *= 0.5; }
  return tmax;
}

//void GPMollerBhabhaModel::SetParticle(const G4ParticleDefinition* p)
FQUALIFIER
void GPMollerBhabhaModel::SetElectron(G4bool eflag)
{
  //  particle = p;
  //  if(p != theElectron) { isElectron = false; }
  isElectron = eflag;
}

FQUALIFIER G4double 
GPMollerBhabhaModel::ComputeCrossSectionPerElectron(G4double kineticEnergy,
						    G4double cutEnergy,
						    G4double maxEnergy)
{
  //  if(!particle) { SetParticle(p); }

  G4double cross = 0.0;
  //  G4double tmax = MaxSecondaryEnergy(p, kineticEnergy);
  G4double tmax = MaxSecondaryEnergy(kineticEnergy);
  tmax = fmin(maxEnergy, tmax);

  if(cutEnergy < tmax) {

    G4double xmin  = cutEnergy/kineticEnergy;
    G4double xmax  = tmax/kineticEnergy;
    G4double tau   = kineticEnergy/electron_mass_c2;
    G4double gam   = tau + 1.0;
    G4double gamma2= gam*gam;
    G4double beta2 = tau*(tau + 2)/gamma2;

    //Moller (e-e-) scattering
    if (isElectron) {

      G4double gg = (2.0*gam - 1.0)/gamma2;
      cross = ((xmax - xmin)*(1.0 - gg + 1.0/(xmin*xmax)
			      + 1.0/((1.0-xmin)*(1.0 - xmax)))
            - gg*log( xmax*(1.0 - xmin)/(xmin*(1.0 - xmax)) ) ) / beta2;

    //Bhabha (e+e-) scattering
    } else {

      G4double y   = 1.0/(1.0 + gam);
      G4double y2  = y*y;
      G4double y12 = 1.0 - 2.0*y;
      G4double b1  = 2.0 - y2;
      G4double b2  = y12*(3.0 + y2);
      G4double y122= y12*y12;
      G4double b4  = y122*y12;
      G4double b3  = b4 + y122;

      cross = (xmax - xmin)*(1.0/(beta2*xmin*xmax) + b2
            - 0.5*b3*(xmin + xmax)
	    + b4*(xmin*xmin + xmin*xmax + xmax*xmax)/3.0)
            - b1*log(xmax/xmin);
    }

    cross *= twopi_mc2_rcl2/kineticEnergy;
  }
  return cross;
}

FQUALIFIER
G4double GPMollerBhabhaModel::ComputeCrossSectionPerAtom(G4double kineticEnergy,
							 G4double Z,
							 G4double cutEnergy,
							 G4double maxEnergy)
{
  return Z*ComputeCrossSectionPerElectron(kineticEnergy,cutEnergy,maxEnergy);
}

FQUALIFIER
G4double GPMollerBhabhaModel::CrossSectionPerVolume(GPMaterial* material,
						    G4double kinEnergy,
						    G4double cutEnergy,
						    G4double maxEnergy)
{
  G4double eDensity = GPMaterial_GetElectronDensity(material);
  return eDensity*ComputeCrossSectionPerElectron(kinEnergy,cutEnergy,maxEnergy);
}

FQUALIFIER
G4double GPMollerBhabhaModel::ComputeDEDXPerVolume(GPMaterial* material,
						   G4double kineticEnergy,
						   G4double cut)
{
  GPIonisParamMat aIonisParamMat;
  GPIonisParamMat_Constructor(&aIonisParamMat,material);

  //  if(!particle) { SetParticle(p); }

  // calculate the dE/dx due to the ionization by Seltzer-Berger formula
  // checl low-energy limit
  G4double electronDensity = GPMaterial_GetElectronDensity(material);
  
  G4double Zeff  = 
    electronDensity/GPMaterial_GetTotNbOfAtomsPerVolume(material);
  G4double th    = 0.25*sqrt(Zeff)*keV;
  G4double tkin  = kineticEnergy;
  if (kineticEnergy < th) { tkin = th; }
 
  G4double tau   = tkin/electron_mass_c2;
  G4double gam   = tau + 1.0;
  G4double gamma2= gam*gam;
  G4double bg2   = tau*(tau + 2);
  G4double beta2 = bg2/gamma2;

  //  G4double eexc  = material->GetIonisation()->GetMeanExcitationEnergy();
  G4double eexc  = GPIonisParamMat_GetMeanExcitationEnergy(&aIonisParamMat);
  eexc          /= electron_mass_c2;
  G4double eexc2 = eexc*eexc; 

  G4double d = fmin(cut, MaxSecondaryEnergy(tkin))/electron_mass_c2;
  G4double dedx;

  // electron
  if (isElectron) {

    dedx = log(2.0*(tau + 2.0)/eexc2) - 1.0 - beta2
         + log((tau-d)*d) + tau/(tau-d)
         + (0.5*d*d + (2.0*tau + 1.)*log(1. - d/tau))/gamma2;
   
  //positron
  } else {

    G4double d2 = d*d*0.5;
    G4double d3 = d2*d/1.5;
    G4double d4 = d3*d*0.75;
    G4double y  = 1.0/(1.0 + gam);
    dedx = log(2.0*(tau + 2.0)/eexc2) + log(tau*d)
         - beta2*(tau + 2.0*d - y*(3.0*d2 
         + y*(d - d3 + y*(d2 - tau*d3 + d4))))/tau;
  } 

  //density correction 
  G4double x = log(bg2)/twoln10;
  //  dedx -= material->GetIonisation()->DensityCorrection(x); 
  dedx -= GPIonisParamMat_DensityCorrection(&aIonisParamMat,x); 

  // now you can compute the total ionization loss
  dedx *= twopi_mc2_rcl2*electronDensity/beta2;
  if (dedx < 0.0) { dedx = 0.0; }
 
  // lowenergy extrapolation
 
  if (kineticEnergy < th) {
    x = kineticEnergy/th;
    if(x > 0.25) { dedx /= sqrt(x); }
    else { dedx *= 1.4*sqrt(x)/(0.1 + x); }
  }
  return dedx;
}

FQUALIFIER
void GPMollerBhabhaModel::SampleSecondaries(GXTrack* track,
					    GPMaterial* material,
					    G4double cutEnergy,
					    G4double maxEnergy)
{
  G4double kineticEnergy = track->E; //dp->GetKineticEnergy();
  G4double tmax;
  G4double tmin = cutEnergy;  
  if(isElectron) { 
    tmax = 0.5*kineticEnergy; 
  } else {
    tmax = kineticEnergy; 
  }
  if(maxEnergy < tmax) { tmax = maxEnergy; }
  if(tmin >= tmax) { return; }

  G4double energy = kineticEnergy + electron_mass_c2;
  G4double totalMomentum = sqrt(kineticEnergy*(energy + electron_mass_c2));
  G4double xmin   = tmin/kineticEnergy;
  G4double xmax   = tmax/kineticEnergy;
  G4double gam    = energy/electron_mass_c2;
  G4double gamma2 = gam*gam;
  G4double beta2  = 1.0 - 1.0/gamma2;
  G4double x, z, q, grej;

  //  G4ThreeVector direction = dp->GetMomentumDirection();
  GPThreeVector direction = 
    GPThreeVector_unit(GPThreeVector_create(track->px,track->py,track->pz));

  //Moller (e-e-) scattering
  if (isElectron) {

    G4double gg = (2.0*gam - 1.0)/gamma2;
    G4double y = 1.0 - xmax;
    grej = 1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));

    do {
      q = GPUniformRand(fDevStates,fThreadId);
      x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      y = 1.0 - x;
      z = 1.0 - gg*x + x*x*(1.0 - gg + (1.0 - gg*y)/(y*y));
    } while(grej * GPUniformRand(fDevStates,fThreadId) > z);

  //Bhabha (e+e-) scattering
  } else {

    G4double y   = 1.0/(1.0 + gam);
    G4double y2  = y*y;
    G4double y12 = 1.0 - 2.0*y;
    G4double b1  = 2.0 - y2;
    G4double b2  = y12*(3.0 + y2);
    G4double y122= y12*y12;
    G4double b4  = y122*y12;
    G4double b3  = b4 + y122;

    y    = xmax*xmax;
    grej = 1.0 + (y*y*b4 - xmin*xmin*xmin*b3 + y*b2 - xmin*b1)*beta2; 
    do {
      q = GPUniformRand(fDevStates,fThreadId);
      x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
      y = x*x;
      z = 1.0 + (y*y*b4 - x*y*b3 + y*b2 - x*b1)*beta2; 
    } while(grej * GPUniformRand(fDevStates,fThreadId) > z);
  }

  G4double deltaKinEnergy = x * kineticEnergy;

  G4double deltaMomentum =
           sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
                                   (deltaMomentum * totalMomentum);
  G4double sint = (1.0 - cost)*(1. + cost);
  if(sint > 0.0) { sint = sqrt(sint); }
  else { sint = 0.0; }

  G4double phi = twopi * GPUniformRand(fDevStates,fThreadId);

  GPThreeVector deltaDirection = 
    GPThreeVector_create(sint*cos(phi),sint*sin(phi), cost) ;

  //  deltaDirection.rotateUz(direction);
  GPThreeVector_rotateUz(&deltaDirection,direction);

  // primary change
  kineticEnergy -= deltaKinEnergy;

  pParticleChange->SetProposedKineticEnergy(kineticEnergy);
  // G4ThreeVector dir = totalMomentum*direction - deltaMomentum*deltaDirection;
  // direction = dir.unit();

  GPThreeVector dir = 
    GPThreeVector_sub(GPThreeVector_mult(direction,totalMomentum),
		      GPThreeVector_mult(deltaDirection,deltaMomentum));
  direction = GPThreeVector_unit(dir);

  pParticleChange->SetProposedMomentumDirection(direction);

  // create G4DynamicParticle object for delta ray
  //  G4DynamicParticle* delta = new G4DynamicParticle(theElectron,
  //				     deltaDirection,deltaKinEnergy);
  //  vdp->push_back(delta);
  //@@@G4FWP fill secondary (one electron)
  FillSecondary(track,deltaDirection,deltaKinEnergy,-1);

}

//---------------------------------------------------------------------------
//
// G4VEmModel
//
//---------------------------------------------------------------------------

FQUALIFIER 
void GPMollerBhabhaModel::SetParticleChange(GPVParticleChange* p)
//					    GPUniversalFluctuation* f)
{
  if(p && pParticleChange != p) { pParticleChange = p; }
  //  flucModel = f;
}

FQUALIFIER 
G4double GPMollerBhabhaModel::HighEnergyLimit()
{
  return highLimit;
}

FQUALIFIER 
G4double GPMollerBhabhaModel::LowEnergyLimit()
{
  return lowLimit;
}

FQUALIFIER 
G4double GPMollerBhabhaModel::HighEnergyActivationLimit()
{
  return eMaxActive;
}

FQUALIFIER 
G4double GPMollerBhabhaModel::LowEnergyActivationLimit()
{
  return eMinActive;
}

FQUALIFIER
void GPMollerBhabhaModel::SetHighEnergyLimit(G4double val)
{
  highLimit = val;
}

FQUALIFIER
void GPMollerBhabhaModel::SetLowEnergyLimit(G4double val)
{
  lowLimit = val;
}

FQUALIFIER 
void GPMollerBhabhaModel::SetActivationHighEnergyLimit(G4double val)
{
  eMaxActive = val;
}

FQUALIFIER 
void GPMollerBhabhaModel::SetActivationLowEnergyLimit(G4double val)
{
  eMinActive = val;
}

FQUALIFIER 
G4bool GPMollerBhabhaModel::IsActive(G4double kinEnergy)
{
  return (kinEnergy >= eMinActive && kinEnergy <= eMaxActive);
}

//---------------------------------------------------------------------------
//
// Utility
//
//---------------------------------------------------------------------------
FQUALIFIER 
GXTrack& GPMollerBhabhaModel::GetSecondary()
{
  return theDeltaRay;
}

FQUALIFIER 
void GPMollerBhabhaModel::FillSecondary(GXTrack* track, 
					GPThreeVector direction,
					G4double energy,
					G4double charge)
{
  G4double p = sqrt(energy*(energy+2.0*electron_mass_c2*charge*charge));
  theDeltaRay.x  = track->x;
  theDeltaRay.y  = track->y;
  theDeltaRay.z  = track->z;
  theDeltaRay.s  = 0;
  theDeltaRay.px = p*direction.x;
  theDeltaRay.py = p*direction.y;
  theDeltaRay.pz = p*direction.z;
  theDeltaRay.E  = energy;
  theDeltaRay.q  = charge;
  
}
