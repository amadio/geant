#include "GUSauterGavrilaModel.h"
#include "GUConstants.h"
#include "ThreeVector.h"

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

GUSauterGavrilaModel::GUSauterGavrilaModel(bool iselectron, const std::string &modelname)
  : EMModel(modelname),
    fIsElectron(iselectron),
    fMinPrimEnergy(1.0*geant::keV),
    fMaxPrimEnergy(1.0*geant::TeV)
{
  SetLowEnergyUsageLimit(fMinPrimEnergy);
  SetHighEnergyUsageLimit(fMaxPrimEnergy);

  //link to vecphys
  fVectorModel = new vecphys::PhotoElectronSauterGavrila(0,-1);
}

GUSauterGavrilaModel::~GUSauterGavrilaModel() {
  delete fVectorModel;
}

void GUSauterGavrilaModel::Initialize() {
  EMModel::Initialize();
  Initialise();
}

void GUSauterGavrilaModel::Initialise() {
  //initialization for vecphys::PhotoElectronSauterGavrila
  fVectorModel->SetLowEnergyLimit(fMinPrimEnergy*vecphys::EScaleToGeant4);
  fVectorModel->SetHighEnergyLimit(fMaxPrimEnergy*vecphys::EScaleToGeant4);
  fVectorModel->Initialization();
}

double GUSauterGavrilaModel::ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy,
                                                        const Particle*)
{
  double xsec = 0.0;
  if (kinenergy < GetLowEnergyUsageLimit() || kinenergy > GetHighEnergyUsageLimit()) {
    return xsec;
  }
  const Material *mat =  matcut->GetMaterial();
  const double *cuts  =  matcut->GetProductionCutsInEnergy();
  double gammacut     =  cuts[0];

  xsec = ComputeXSectionPerVolume(mat, gammacut, kinenergy);

  return xsec;
}

double GUSauterGavrilaModel::ComputeXSectionPerVolume(const Material *mat, double /* prodcutenergy */,
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

double GUSauterGavrilaModel::ComputeXSectionPerAtom(const Element *elem,
                                                    const MaterialCuts * /*matcut*/,
                                                    double kinEnergy,
                                                    const Particle * /*particle*/ )
{
  //interface to vecphys and the unit conversion
  double xsec = fVectorModel->G4CrossSectionPerAtom(elem->GetZ(), kinEnergy*vecphys::EScaleToGeant4);
  return vecphys::invXsecScaleToGeant4*xsec;
}

double GUSauterGavrilaModel::MinimumPrimaryEnergy(const MaterialCuts * /*matcut*/,
                                                  const Particle * /*part*/)  const
{
  return 1.0e-6; // 1.0 * geant::keV;
}

int GUSauterGavrilaModel::SampleSecondaries(LightTrack &track,  Geant::GeantTaskData *td) {

  int    numSecondaries      = 0;

  // conversion
  double energyIn = track.GetKinE()*vecphys::EScaleToGeant4;

  // do nothing if the primary gamma is outside the valid energy range
  // @syjun probably this is a redundant check and can be omitted in the vector mode

  double energyOut = 0;
  double sinTheta = 0;
  const int targetElement = track.GetTargetZ();

  //@syjun atomic deexcitation is not considerred)

  //sample a photo-electron

  fVectorModel-> template InteractKernel<vecphys::ScalarBackend>(energyIn, targetElement, energyOut, sinTheta);

  // update the primary track (photon) and the secondary (electron)
  double phi       = geant::kTwoPi*vecphys::UniformRandom<double>(0,-1);
  double cosTheta  = vecCore::math::Sqrt((1.-sinTheta)*(1+sinTheta));

  //rotate the sampled photo-electron diection with respect to the primary gamma
  vecphys::ThreeVector<double> gamDirection(track.GetDirX(),track.GetDirY(),track.GetDirZ());
  vecphys::ThreeVector<double> eleDirection(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta);
  eleDirection.RotateUz(gamDirection);

  // kill the primary photon
  track.SetKinE(0.0);
  track.SetTrackStatus(LTrackStatus::kKill);

  // create the secondary partcile i.e. the photo e-
  double deltaKinEnergy = (energyIn - energyOut)/vecphys::EScaleToGeant4;

  //@syjun if deltaKinEnergy is below the production thresold,
  //deposite deltaKinEnergy and return with numSecondaries = 0

  vecphys::ThreeVector<double> electronDir = eleDirection.Unit();

  //put photo-electron into the stack of secondaries
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

  //fill photo-electron information and kinematic
  sectracks[secIndx].SetGVcode(Electron::Definition()->GetInternalCode());
  sectracks[secIndx].SetMass(geant::kElectronMassC2);
  sectracks[secIndx].SetTrackIndex(track.GetTrackIndex());

  sectracks[secIndx].SetKinE(deltaKinEnergy);
  sectracks[secIndx].SetDirX(electronDir.x());
  sectracks[secIndx].SetDirY(electronDir.y());
  sectracks[secIndx].SetDirZ(electronDir.z());

  ++secIndx;

  return numSecondaries;
}

} // namespace geantphysics
