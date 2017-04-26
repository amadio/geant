#include "GUKleinNishinaComptonModel.h"
#include "GUConstants.h"
#include "ThreeVector.h"

#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"

#include "MaterialCuts.h"

#include "Electron.h"
#include "Gamma.h"

#include "LightTrack.h"
#include "PhysicsData.h"

#include <iostream>
#include <cmath>

namespace geantphysics {

GUKleinNishinaComptonModel::GUKleinNishinaComptonModel(const std::string &modelname)
  : EMModel(modelname),
    fMinPrimEnergy(1.0*geant::keV),
    fMaxPrimEnergy(1.0*geant::TeV)
{
  SetLowEnergyUsageLimit(fMinPrimEnergy);
  SetHighEnergyUsageLimit(fMaxPrimEnergy);

  //link to vecphys
  fVectorModel = new vecphys::ComptonKleinNishina(0,-1);
}

GUKleinNishinaComptonModel::~GUKleinNishinaComptonModel() {
  delete fVectorModel;
}

void GUKleinNishinaComptonModel::Initialize() {
  EMModel::Initialize();
  Initialise();
}

void GUKleinNishinaComptonModel::Initialise() {
  //initialization for vecphys::ComptonKleinNishina
  std::cout << "  ----> GUKleinNishinaComptonModel [" << GetLowEnergyUsageLimit() 
            << "," << GetHighEnergyUsageLimit() << "] in [geant::GeV]" << std::endl;

  fVectorModel->SetLowEnergyLimit(fMinPrimEnergy*vecphys::EScaleToGeant4);
  fVectorModel->SetHighEnergyLimit(fMaxPrimEnergy*vecphys::EScaleToGeant4);
  fVectorModel->Initialization();
}

double GUKleinNishinaComptonModel::ComputeMacroscopicXSection(const MaterialCuts *matcut,
                                                              double kinenergy, 
                                                              const Particle*) 
{
  double xsec = 0.0;
  if (kinenergy < GetLowEnergyUsageLimit() || kinenergy > GetHighEnergyUsageLimit()) {
    return xsec;
  }
  const Material *mat = matcut->GetMaterial();
  const double *cuts  = matcut->GetProductionCutsInEnergy();
  double gammacut     = cuts[0];

  xsec = ComputeXSectionPerVolume(mat, gammacut, kinenergy);

  return xsec;
}

double GUKleinNishinaComptonModel::ComputeXSectionPerVolume(const Material *mat,
                                                            double /* productionCutEnergy */, 
                                                            double particleekin) 
{
  double xsec = 0.;

  int nelm = mat->GetNumberOfElements();
  auto elementVec = mat->GetElementVector();
  const double *atomNumDensity = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect(); 

  for (int i = 0; i < nelm; ++i) {
    xsec += atomNumDensity[i] * ComputeXSectionPerAtom(elementVec[i],0,particleekin,0);
  }

  return xsec;
}

double GUKleinNishinaComptonModel::ComputeXSectionPerAtom(const Element *elem,
                                                          const MaterialCuts * /*matcut*/,
                                                          double kinEnergy,
                                                          const Particle * /*particle*/)
{
  //interface to vecphys and the unit conversion
  double xsec = fVectorModel->G4CrossSectionPerAtom(elem->GetZ(), kinEnergy*vecphys::EScaleToGeant4);
  return vecphys::invXsecScaleToGeant4*xsec;
}

double
GUKleinNishinaComptonModel::MinimumPrimaryEnergy(const MaterialCuts * /*matcut*/,
                                                 const Particle * /*part*/) const
{ return 1.e-6; }  //  1.0 * geant::keV;

int GUKleinNishinaComptonModel::SampleSecondaries(LightTrack &track,
                                                  std::vector<LightTrack> & /*sectracks*/,
                                                  Geant::GeantTaskData *td)
{
  int numSecondaries = 0;

  // conversion for vecphys
  double energyIn = track.GetKinE()*vecphys::EScaleToGeant4;

  // do nothing if the primary gamma is outside the valid energy range 
  // @syjun probably this is a redundant check and can be omitted in the vector mode
  if( energyIn < fVectorModel->GetLowEnergyLimit() ||
      energyIn > fVectorModel->GetHighEnergyLimit()) { 
    return numSecondaries; 
  };

  double energyOut = 0;
  double sinTheta = 0;
  const int targetElement = track.GetTargetZ();

  //sample a scattered photon
  fVectorModel-> template InteractKernel<vecphys::ScalarBackend>(energyIn, targetElement, energyOut, sinTheta);

  // update the primary track (photon) 
  double phi       = geant::kTwoPi*vecphys::UniformRandom<double>(0,-1);
  double cosTheta  = vecCore::math::Sqrt((1.-sinTheta)*(1+sinTheta));

  //rotate the sampled diection with respect to the line of flight
  vecphys::ThreeVector<double> gamDirection0(track.GetDirX(),track.GetDirY(),track.GetDirZ());
  vecphys::ThreeVector<double> gamDirection1(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta);
  gamDirection1.RotateUz(gamDirection0);

  // update kinematic for the scattered photon 
  double gamEnergy1 = energyOut/vecphys::EScaleToGeant4;
  if(gamEnergy1 > fMinPrimEnergy) {
    //@syjun or above the primary/secondary thresold
    track.SetKinE(gamEnergy1);
    track.SetDirX(gamDirection1.x());
    track.SetDirY(gamDirection1.y());
    track.SetDirZ(gamDirection1.z());
  }
  else {
    track.SetTrackStatus(LTrackStatus::kKill);
    track.SetKinE(0.0);
    track.SetEnergyDeposit(gamEnergy1);
  }

  // create the secondary partcile i.e. the scattered e-
  double deltaKinEnergy = (energyIn - energyOut)/vecphys::EScaleToGeant4;

  //@syjun if deltaKinEnergy is below the production thresold, 
  //       then deposit 'deltaKinEnergy' and return with numSecondaries = 0

  vecphys::ThreeVector<double> electronDir = energyIn*gamDirection0 - energyOut*gamDirection1;

  /*
  // Put a verbosity flag 
  // Check E/p relation for outgoing electron
  {
     double eKin = (energyIn - energyOut);
     double momentumMag2 = electronDir.Mag2();
     double eMass= geant::kElectronMassC2;
     
     double eKinFromP =  momentumMag2 /
        ( std::sqrt( momentumMag2 + eMass * eMass ) + eMass );
     double diffEkin= eKinFromP - eKin ;
     const double thrDiff= 1.0e-3;  // Threshold level for relative difference
     if( std::fabs(diffEkin) > thrDiff * eKin ){
        std::cerr << " WARNING> Difference in electron eKin = "
                  << diffEkin << " beyond threshold level= " << thrDiff
                  << "  eKin= "  << eKin << "  and from-P = " << eKinFromP
                  << std::endl;
     }
  }
  */
  
  electronDir = electronDir.Unit();

  //put one electron into the stack of secondaries
  numSecondaries = 1;

  int curSizeOfSecList = td->fPhysicsData->GetSizeListOfSecondaries();
  int curNumUsedSecs   = td->fPhysicsData->GetNumUsedSecondaries();

  if (curSizeOfSecList-curNumUsedSecs<numSecondaries) {
    td->fPhysicsData->SetSizeListOfSecondaries(2*curSizeOfSecList);
  }

  int secIndx = curNumUsedSecs;
  curNumUsedSecs +=numSecondaries;
  td->fPhysicsData->SetNumUsedSecondaries(curNumUsedSecs);
  std::vector<LightTrack>& sectracks = td->fPhysicsData->GetListOfSecondaries();

  //fill electron information and kinematic
  sectracks[secIndx].SetGVcode(Electron::Definition()->GetInternalCode()); 
  sectracks[secIndx].SetMass(geant::kElectronMassC2);
  sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); 

  sectracks[secIndx].SetKinE(deltaKinEnergy);
  sectracks[secIndx].SetDirX(electronDir.x());
  sectracks[secIndx].SetDirY(electronDir.y());
  sectracks[secIndx].SetDirZ(electronDir.z());

  secIndx++;

  return numSecondaries;
}

} // namespace geantphysics
