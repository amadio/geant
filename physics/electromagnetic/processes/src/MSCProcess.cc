
#include "MSCProcess.h"

#include "SystemOfUnits.h"

#include "Electron.h"
#include "Positron.h"

#include "Region.h"

#include "MSCModel.h"
#include "EMModelManager.h"

// from geantV
#include "GeantTaskData.h"
#include "GeantTrack.h"

#include "Geant/NavigationInterface.h"

namespace geantphysics {

MSCProcess::MSCProcess(const std::string &name) : EMPhysicsProcess(name), fGeomMinLimit(0.05*geant::nm) {
  fGeomMinLimit2 = fGeomMinLimit*fGeomMinLimit;
  // process type is kElectromagnetic in the base EMPhysicsProcess calss so set it to kMSC
  SetType(ProcessType::kMSC);
  // set to be a pure continuous process
  SetIsContinuous(true);
  SetIsDiscrete(true); // just for the framework: msc has no discrete part (at the moment)
  // fill the list of particles that this process can be used to i.e. gamma particle
  AddToListParticlesAlloedToAssigned(Electron::Definition());
  AddToListParticlesAlloedToAssigned(Positron::Definition());
}

void MSCProcess::Initialize()
{
  // Call initialization via EMPhysicsProcess, then register MSCdata
  EMPhysicsProcess::Initialize();
  const Geant::TrackToken* mscdata = Geant::TrackDataMgr::GetInstance()->GetToken("MSCdata");
  assert(mscdata);
  fMSCdata = *mscdata;
}

MSCProcess::~MSCProcess() {}

// called at the PrePropagationStage(in the Handler)
double MSCProcess::AlongStepLimitationLength(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) const {
  // init all lengths to the current minimum physics step length (that is the true length)
  //
  MSCdata &mscdata = ((Geant::TrackToken)fMSCdata).Data<MSCdata>(gtrack);

  bool isOnBoundaryPostStp = gtrack->Boundary();
  gtrack->SetBoundary(gtrack->IsOnBoundaryPreStp());
  // get the phyics step limit due to all (other than msc) physics processes (that was used in the geometry stage)
  double minPhysicsStepLength   = gtrack->GetPstep();
  mscdata.fTheTrueStepLenght    = minPhysicsStepLength;
  mscdata.fTheTransportDistance = minPhysicsStepLength;
  mscdata.fTheZPathLenght       = minPhysicsStepLength;
  mscdata.SetDisplacement(0.,0.,0.);
  mscdata.SetNewDirectionMsc(0.,0.,1.);
  // select msc model
  double ekin        = gtrack->T();
  int    regIndx     = const_cast<vecgeom::LogicalVolume*>(gtrack->GetVolume())->GetRegion()->GetIndex();
  MSCModel *mscModel = static_cast<MSCModel*>(GetModelManager()->SelectModel(ekin,regIndx));
  // check if:
  // - current true step length is > limit
  // - current kinetic energy is above the minimum
  if (mscModel && minPhysicsStepLength>GetGeomMinLimit()) {
    // will check min/max usage limits, possible additional step limit, update the
    // mscdata.fTheTrueStepLenght and will convert the true->to->geometic length (that will be in mscdata.fTheZPathLenght)
    mscModel->StepLimit(gtrack, td);
    // check if msc limited the step: set that continuous process was the winer
    if (mscdata.fTheTrueStepLenght<minPhysicsStepLength) {
      gtrack->SetEindex(-1); // indicate that continuous process was the winer
      gtrack->SetProcess(GetGlobalIndex());       // set global indx of limiting process
    }
  }
  // update gtrack->fPstep to be the geometric step length (this is what propagation needs):
  // straight line discrete along the original direction
  // protection: geometric legth must always be <= than the true physics step length
  mscdata.fTheZPathLenght = std::min(mscdata.fTheZPathLenght, minPhysicsStepLength);
  gtrack->SetPstep(mscdata.fTheZPathLenght);
  // check if the geometrical physics step (true is change now by MSC to geometric) become shorter than snext and limit
  // snext to this geometrical physics step (track will be propagated to snext distance) and set the post-step point
  // boundary flag to false (we pulled back)
  if (gtrack->GetSnext() > gtrack->GetPstep()) {
    gtrack->SetSnext(gtrack->GetPstep());
    gtrack->SetBoundary(false);
  } else { // write back the original post-step point boundary flag and do nothing
    gtrack->SetBoundary(isOnBoundaryPostStp);
  }
  return 0.;// not used at all because the step limit is written directly into the GeantTrack (just the interface)
}

// called at the PostPropagationStage(in the Handler)
void MSCProcess::AlongStepDoIt(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) const {
  // get the step length (the geometric one)
  double geometricStepLength = gtrack->GetStep();
  double truePathLength      = geometricStepLength;
  MSCdata &mscdata = ((Geant::TrackToken)fMSCdata).Data<MSCdata>(gtrack);
  // select msc model
  double ekin        = gtrack->T();
  int    regIndx     = const_cast<vecgeom::LogicalVolume*>(gtrack->GetVolume())->GetRegion()->GetIndex();
  MSCModel *mscModel = static_cast<MSCModel*>(GetModelManager()->SelectModel(ekin,regIndx));
  // check if:
  // - current true step length is > limit
  // - current kinetic energy is above the minimum usea
  if (mscModel && mscdata.fTheTrueStepLenght>GetGeomMinLimit()) {
    truePathLength = mscdata.fTheTrueStepLenght;
    // mscdata.fTheTrueStepLenght be updated during the conversion
    mscModel->ConvertGeometricToTrueLength(gtrack, td);
    // protection againts wrong true-geometic-true gonversion
    truePathLength = std::min(truePathLength,mscdata.fTheTrueStepLenght);
    // optimization: scattring is not sampled if the particle is reanged out in this step or short step
    if (mscdata.fRange>truePathLength && truePathLength>GetGeomMinLimit()) {
      // sample scattering: might have been done during the step limit phase
      // NOTE: in G4 the SampleScattering method won't use the possible shrinked truePathLength!!! but make it correct
      mscdata.fTheTrueStepLenght = truePathLength;
      bool hasNewDir = mscModel->SampleScattering(gtrack, td);
      // compute displacement vector length
      double dl =   mscdata.fTheDisplacementVectorX*mscdata.fTheDisplacementVectorX
                  + mscdata.fTheDisplacementVectorY*mscdata.fTheDisplacementVectorY
                  + mscdata.fTheDisplacementVectorZ*mscdata.fTheDisplacementVectorZ;
      if (dl>fGeomMinLimit2 && !gtrack->Boundary()) {
        // displace the post-step point
        dl = std::sqrt(dl);
        double dir[3]={mscdata.fTheDisplacementVectorX/dl, mscdata.fTheDisplacementVectorY/dl, mscdata.fTheDisplacementVectorZ/dl};
        ScalarNavInterface::DisplaceTrack(*gtrack,dir,dl,GetGeomMinLimit());
      }
      // apply msc agular deflection
//      if (!gtrack->Boundary() && hasNewDir) {
      if (hasNewDir)
        gtrack->SetDirection(mscdata.fTheNewDirectionX, mscdata.fTheNewDirectionY, mscdata.fTheNewDirectionZ);
    }
  }
  // update step length to store the true step length
  gtrack->SetStep(truePathLength);
}



}  // namespace geantphysics
