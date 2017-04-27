#include "GUBetheHeitlerConversionModel.h"
#include "GUConstants.h"
#include "ThreeVector.h"

#include "PhysicalConstants.h"
#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"
#include "MaterialCuts.h"
#include "MaterialHandler.h"

#include "Positron.h"
#include "Electron.h"
#include "Gamma.h"

#include "LightTrack.h"
#include "PhysicsData.h"

#include <iostream>
#include <cmath>

namespace geantphysics {

GUBetheHeitlerConversionModel::GUBetheHeitlerConversionModel(const std::string &modelname)
  : EMModel(modelname), 
    fMinPrimEnergy(2.0*geant::kElectronMassC2),
    fMaxPrimEnergy(1.0*geant::TeV)
{
  SetLowEnergyUsageLimit(fMinPrimEnergy);
  SetHighEnergyUsageLimit(fMaxPrimEnergy);

  //link to vecphys
  fVectorModel = new vecphys::ConversionBetheHeitler(0,-1);
}

GUBetheHeitlerConversionModel::~GUBetheHeitlerConversionModel() {
  delete fVectorModel;
}

void GUBetheHeitlerConversionModel::Initialize() {
  EMModel::Initialize();
  Initialise();

  //for the atomic dependent model
  EMModel::InitialiseElementSelectors(this,Gamma::Definition(),true);
}

void GUBetheHeitlerConversionModel::Initialise() {
  //initialization for vecphys::ConversionBetheHeitler
  fVectorModel->SetLowEnergyLimit(fMinPrimEnergy*vecphys::EScaleToGeant4);
  fVectorModel->SetHighEnergyLimit(fMaxPrimEnergy*vecphys::EScaleToGeant4);

  //pass material information to MaterialHandler
  Vector_t<Material*> mtable =  Material::GetTheMaterialTable();
  std::vector<int> etable;
  for(unsigned int im = 0 ; im < mtable.size() ; ++im) {
    Vector_t<Element*> ev = mtable[im]->GetElementVector();
    for(unsigned int ie = 0 ; ie < ev.size() ; ++ie) { 
      int z = static_cast<int>(ev[ie]->GetZ());
      if (std::find(etable.begin(), etable.end(), z) != etable.end()) continue;
      etable.push_back(z);
    }
  }
  vecphys::MaterialHandler::Instance()->BuildElementTable(etable);

  fVectorModel->Initialization();
}

double GUBetheHeitlerConversionModel::ComputeMacroscopicXSection(const MaterialCuts *matcut,
                                                                 double kinenergy, 
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

double GUBetheHeitlerConversionModel::ComputeXSectionPerVolume(const Material *mat, double /* prodcutenergy */, 
                                                                double particleekin) 
{
  double xsec = 0.;

  int nelm = mat->GetNumberOfElements();
  auto elementVec = mat->GetElementVector();
  const double *atomNumDensity = mat->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect(); 

  for (int i = 0; i < nelm; ++i) {
    xsec += atomNumDensity[i]*ComputeXSectionPerAtom(elementVec[i],0,particleekin,0);
  }

  return xsec;
}

double GUBetheHeitlerConversionModel::ComputeXSectionPerAtom(const Element *elem,
                                                             const MaterialCuts * /*matcut*/,
                                                             double kinEnergy,
                                                             const Particle * /*particle*/ )
{
  //interface to vecphys and the unit conversion
  double xsec = fVectorModel->G4CrossSectionPerAtom(elem->GetZ(), kinEnergy*vecphys::EScaleToGeant4);
  return vecphys::invXsecScaleToGeant4*xsec;
}

double
GUBetheHeitlerConversionModel::MinimumPrimaryEnergy(const MaterialCuts * /*matcut*/,
                                                    const Particle * /*part*/) const
{
  return 2.0*geant::kElectronMassC2;
}
   
int GUBetheHeitlerConversionModel::SampleSecondaries(LightTrack &track, std::vector<LightTrack> & /*sectracks*/,
                                                     Geant::GeantTaskData *td)
{
  int    numSecondaries      = 0;

  // conversion for vecphys
  double energyIn = track.GetKinE()*vecphys::EScaleToGeant4;

  // threshold
  if(track.GetKinE() < 2.0*geant::kElectronMassC2 ) return numSecondaries;

  double energyOut = 0;
  double sinTheta = 0;

  //select an element for atomic dependent models - this step has a serious vecotrization issue 
  const int cutIndex = track.GetMaterialCutCoupleIndex();
  const MaterialCuts *matcut = (MaterialCuts::GetTheMaterialCutsTable())[cutIndex];
  int elementIndex = EMModel::SampleTargetElementIndex(matcut, track.GetKinE(),vecphys::UniformRandom<double>(0,-1));
  int targetElement = ((matcut->GetMaterial()-> GetElementVector())[elementIndex])->GetZ();

  //sample a conversion electron
  fVectorModel-> template InteractKernel<vecphys::ScalarBackend>(energyIn, targetElement, energyOut, sinTheta);

  // update the primary track (photon) - i.e.,  kill the primary
  track.SetKinE(0.0);
  track.SetTrackStatus(LTrackStatus::kKill);

  // create secondary partciles i.e. a pair of e- and e+
  double electronKinEnergy = energyOut/vecphys::EScaleToGeant4 - geant::kElectronMassC2;
  double positronKinEnergy = (energyIn-energyOut)/vecphys::EScaleToGeant4 - geant::kElectronMassC2;

  double phi       = geant::kTwoPi*vecphys::UniformRandom<double>(0,-1);
  double cosTheta  = vecCore::math::Sqrt((1.-sinTheta)*(1+sinTheta));

  vecphys::ThreeVector<double> gamDirection(track.GetDirX(),track.GetDirY(),track.GetDirZ());
  vecphys::ThreeVector<double> eleDirection(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta);
  eleDirection.RotateUz(gamDirection);
  vecphys::ThreeVector<double> electronDir = eleDirection.Unit();

  double sinTheta1 = vecCore::math::ASin(sinTheta)*energyOut/(energyIn-energyOut);
  double cosTheta1  = vecCore::math::Sqrt((1.-sinTheta1)*(1+sinTheta1));

  vecphys::ThreeVector<double> posDirection(-sinTheta1*cos(phi), -sinTheta1*sin(phi), cosTheta1);
  posDirection.RotateUz(gamDirection);
  vecphys::ThreeVector<double> positronDir = posDirection.Unit();

  //put one electron into the stack
  numSecondaries = 2;

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

  sectracks[secIndx].SetKinE(electronKinEnergy);
  sectracks[secIndx].SetDirX(electronDir.x());
  sectracks[secIndx].SetDirY(electronDir.y());
  sectracks[secIndx].SetDirZ(electronDir.z());

  ++secIndx;

  //fill positron information and and kinematic
  sectracks[secIndx].SetGVcode(Positron::Definition()->GetInternalCode()); 
  sectracks[secIndx].SetMass(geant::kElectronMassC2);
  sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); 

  sectracks[secIndx].SetKinE(positronKinEnergy);
  sectracks[secIndx].SetDirX(positronDir.x());
  sectracks[secIndx].SetDirY(positronDir.y());
  sectracks[secIndx].SetDirZ(positronDir.z());

  ++secIndx;

  return numSecondaries;
}

} // namespace geantphysics
