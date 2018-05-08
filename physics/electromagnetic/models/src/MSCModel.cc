
#include <Geant/PhysicsData.h>
#include "Geant/MSCModel.h"

// from geantV
#include "Geant/TaskData.h"
#include "Geant/Track.h"
#include "Geant/NavigationInterface.h"

namespace geantphysics {

using geant::Double_v;
using geant::IndexD_v;
using geant::kVecLenD;
using geant::MaskD_v;
using vecCore::Get;
using vecCore::Set;
using vecCore::AssignMaskLane;
using vecCore::MaskFull;
using vecCore::MaskEmpty;

MSCModel::MSCModel(const std::string &name)
    : EMModel(name), fRangeFactor(0.06), fSafetyFactor(0.6), fGeomFactor(2.5), fSkin(3.),
      fMSCSteppingAlgorithm(MSCSteppingAlgorithm::kUseSaftey)
{
  fIsMSCModel = true;
}

MSCModel::~MSCModel()
{
}

void MSCModel::Initialize()
{
  // call the base class Initialize method
  EMModel::Initialize();
  const geant::TrackToken *mscdata = geant::TrackDataMgr::GetInstance()->GetToken("MSCdata");
  assert(mscdata);
  fMSCdata = *mscdata;
}

void MSCModel::AlongStepDoIt(geant::Track *gtrack, geant::TaskData *td)
{
  MSCdata &mscdata      = fMSCdata.Data<MSCdata>(gtrack);
  double truePathLength = mscdata.fTheTrueStepLenght;
  bool hasNewDir        = SampleScattering(gtrack, td);
  // compute displacement vector length
  double dl = mscdata.fTheDisplacementVectorX * mscdata.fTheDisplacementVectorX +
              mscdata.fTheDisplacementVectorY * mscdata.fTheDisplacementVectorY +
              mscdata.fTheDisplacementVectorZ * mscdata.fTheDisplacementVectorZ;
  if (dl > fGeomMinLimit2 && !gtrack->Boundary()) {
    // displace the post-step point
    dl            = std::sqrt(dl);
    double dir[3] = {mscdata.fTheDisplacementVectorX / dl, mscdata.fTheDisplacementVectorY / dl,
                     mscdata.fTheDisplacementVectorZ / dl};
    ScalarNavInterface::DisplaceTrack(*gtrack, dir, dl, GetGeomMinLimit());
  }
  // apply msc agular deflection
  //      if (!gtrack->Boundary() && hasNewDir) {
  if (hasNewDir) gtrack->SetDirection(mscdata.fTheNewDirectionX, mscdata.fTheNewDirectionY, mscdata.fTheNewDirectionZ);
  // update step length to store the true step length
  gtrack->SetStep(truePathLength);
}

void MSCModel::AlongStepDoIt(std::vector<geant::Track *> &gtracks, geant::TaskData *td)
{
  std::vector<bool> &hasNewDir         = td->fPhysicsData->fPhysicsScratchpad.hasNewDir;
  std::vector<double> &truePathLengths = td->fPhysicsData->fPhysicsScratchpad.truePathLengths;
  hasNewDir.resize(gtracks.size());
  truePathLengths.clear();

  for (auto gtrack : gtracks) {
    MSCdata &mscdata = fMSCdata.Data<MSCdata>(gtrack);
    // sample scattering: might have been done during the step limit phase
    // NOTE: in G4 the SampleScattering method won't use the possible shrinked truePathLength!!! but make it correct
    truePathLengths.push_back(mscdata.fTheTrueStepLenght);
  }

  SampleScattering(gtracks, hasNewDir, td);

  for (size_t i = 0; i < gtracks.size(); i += kVecLenD) {
    MaskD_v displMask = MaskD_v(false);
    Double_v displX, displY, displZ;

    for (int l = 0; l < kVecLenD; ++l) {
      auto gtrack      = gtracks[i + l];
      MSCdata &mscdata = fMSCdata.Data<MSCdata>(gtrack);
      Set(displMask, l, !gtrack->Boundary());
      Set(displX, l, mscdata.fTheDisplacementVectorX);
      Set(displY, l, mscdata.fTheDisplacementVectorY);
      Set(displZ, l, mscdata.fTheDisplacementVectorZ);
      if (hasNewDir[i + l])
        gtrack->SetDirection(mscdata.fTheNewDirectionX, mscdata.fTheNewDirectionY, mscdata.fTheNewDirectionZ);
      gtrack->SetStep(truePathLengths[i + l]);
    }

    Double_v dl = displX * displX + displY * displY + displZ * displZ;
    displMask   = displMask && dl > fGeomMinLimit2;
    if (!vecCore::MaskEmpty(displMask)) {
      Double_v dlSqrt = Math::Sqrt(dl);
      displX /= dlSqrt;
      displY /= dlSqrt;
      displZ /= dlSqrt;

      for (int l = 0; l < kVecLenD; ++l) {
        auto gtrack   = gtracks[i + l];
        double dir[3] = {Get(displX, l), Get(displY, l), Get(displZ, l)};
        if (Get(displMask, l)) ScalarNavInterface::DisplaceTrack(*gtrack, dir, Get(dlSqrt, l), GetGeomMinLimit());
      }
    }
  }
}

bool MSCModel::SamplingNeeded(geant::Track *gtrack, geant::TaskData *td)
{
  double geometricStepLength = gtrack->GetStep();
  double truePathLength      = geometricStepLength;
  MSCdata &mscdata           = fMSCdata.Data<MSCdata>(gtrack);
  if (mscdata.fTheTrueStepLenght <= GetGeomMinLimit()) {
    return false;
  }

  truePathLength = mscdata.fTheTrueStepLenght;
  // mscdata.fTheTrueStepLenght be updated during the conversion
  ConvertGeometricToTrueLength(gtrack, td);
  // protection againts wrong true-geometic-true gonversion
  truePathLength = std::min(truePathLength, mscdata.fTheTrueStepLenght);
  // optimization: scattring is not sampled if the particle is reanged out in this step or short step

  if (mscdata.fRange <= truePathLength || truePathLength <= GetGeomMinLimit()) {
    gtrack->SetStep(truePathLength);
    return false;
  }

  mscdata.fTheTrueStepLenght = truePathLength;

  return true;
}

} // namespace geantphysics
