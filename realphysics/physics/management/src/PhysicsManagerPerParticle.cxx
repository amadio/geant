#include "PhysicsManagerPerParticle.h"

#include "PhysicsParameters.h"
#include "Material.h"
#include "MaterialCuts.h"

#include "PhysicsProcess.h"
#include "Particle.h"

#include "SystemOfUnits.h"

#include "LightTrack.h"

#include <cmath>
#include <iostream>


namespace geantphysics {

PhysicsManagerPerParticle::PhysicsManagerPerParticle(const Particle *partcile, const PhysicsParameters *physpars)
: fParticle(partcile), fPhysicsParameters(physpars) {
  fIsHasMSCProcess            = false;
  fIsHasDecayProcess          = false;
  fIsHasElossProcess          = false;
}

// All process object pointers are stored in the PhysicsProcess and all process object is deleted by calling
// PhysicsProcess::ClearAllProcess() static method.
// So here we need to delete only local things
PhysicsManagerPerParticle::~PhysicsManagerPerParticle() {}


void PhysicsManagerPerParticle::Initialize() {
  // Initialize all processes assigned to the particle
  // Note: processes will initialise their models as well if any (especially in case of EM-processes)
  std::cerr<<"==================  Particle = "<<fParticle->GetName()<< "  has "<<fProcessVec.size()<<" process."<<std::endl;
  for (unsigned long ip=0; ip<fProcessVec.size(); ++ip) {
    fProcessVec[ip]->SetPhysicsParameters(fPhysicsParameters);
    fProcessVec[ip]->SetParticle(fParticle);
    fProcessVec[ip]->Initialize();
    if (fProcessVec[ip]->IsLambdaTablerequested()) {
//      std::cerr<< "  ===== Building Lambda table for process = " << fProcessVec[ip]->GetName() << std::endl;
      fProcessVec[ip]->BuildLambdaTables();
    }
  }
  // order processes and set some flags
  CheckProcesses();
}


void PhysicsManagerPerParticle::AddProcess( PhysicsProcess *process ) {
  if (process == nullptr) {
    std::cerr<<"   **** ERROR: PhysicsManagerPerParticle::AddProcess()\n"
             <<"    The process pointer is nullptr! \n";
    exit(-1);
  }
  if (process->GetForcedCondition()==ForcedCondition::kInActivated)
    return;
  // the process pointer will be added to this physics manager per particle
  process->SetIndex(fProcessVec.size());
  fProcessVec.push_back(process);
  // push also to sub categories
  bool isForced = (process->GetForcedCondition()!=ForcedCondition::kNotForced);
  if (process->GetIsContinuous()) {
    fAlongStepProcessVec.push_back( process );
  }
  if (process->GetIsDiscrete()) {
    fPostStepCandidateProcessVec.push_back(process);
    if (isForced) {
      fPostStepForcedProcessVec.push_back(process);
    }
  }
  if (process->GetIsAtRest())  {
    fAtRestCandidateProcessVec.push_back(process);
    if (isForced) {
      fAtRestForcedProcessVec.push_back(process);
    }
  }
}


const PhysicsProcess* PhysicsManagerPerParticle::GetMSCProcess() const {
  PhysicsProcess *mscProc = nullptr;
  if (fIsHasMSCProcess) {
    mscProc = fPostStepCandidateProcessVec[fPostStepCandidateProcessVec.size()-1];
  }
  return mscProc;
}


void PhysicsManagerPerParticle::CheckProcesses() {
  std::vector<int> indxVect;
  int indxMSC   = -1;
  for (unsigned long i=0; i<fPostStepCandidateProcessVec.size(); ++i) {
    PhysicsProcess *proc = fPostStepCandidateProcessVec[i];
    if (proc->GetType()==ProcessType::kMSC) { // msc will be pure continuous !!!
      fIsHasMSCProcess = true;
      indxMSC          = i;
    }
    if (proc->GetType()==ProcessType::kEnergyLoss) {
      fIsHasElossProcess         = true;
    }
    if (proc->GetType()!=ProcessType::kMSC) {
      indxVect.push_back(i);
    }
  }
  std::vector<PhysicsProcess*> theCopy(fPostStepCandidateProcessVec.size());
  for (unsigned long i=0; i<fPostStepCandidateProcessVec.size(); ++i) {
    theCopy[i] = fPostStepCandidateProcessVec[i];
  }
  for (unsigned long i=0; i<indxVect.size(); ++i) {
    fPostStepCandidateProcessVec[i] = theCopy[indxVect[i]];
  }
  if (fIsHasMSCProcess) {
    fPostStepCandidateProcessVec[fPostStepCandidateProcessVec.size()-1] = theCopy[indxMSC];

  }
}


void PhysicsManagerPerParticle::PrepareForRun() {
  // keep only one kEnergyLoss process pointer in the fAlongStepProcessVec if there are more than 1
  if (fIsHasElossProcess) {
    bool isElossProcessFound = false;
    for (unsigned long i=0; i<fAlongStepProcessVec.size(); ++i) {
      if (fAlongStepProcessVec[i]->GetType()==ProcessType::kEnergyLoss) {
        if (!isElossProcessFound) {
          isElossProcessFound = true;
        } else { // remove the process pointer if other kEnergyLoss process had been found before
          fAlongStepProcessVec.erase(fAlongStepProcessVec.begin()+i);
        }
      }
    }
  }
}


void PhysicsManagerPerParticle::ComputeIntLen(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) {
  // mfp (gtrack->fPhysicsInteractLength[]) is initialized to 1.0 for all possible processes to handle properly
  // particles that has no interaction when number-of-interaction length-left is updated
  //
  double limit      = PhysicsProcess::GetAVeryLargeValue();
  // 1. First the continuous step limit if any
  // if it has MSC than go to last but one
  size_t maxContIndx = fAlongStepProcessVec.size();
  if (fIsHasMSCProcess) {
    --maxContIndx;
  }
  for (size_t i=0; i<maxContIndx; ++i) {
    //if (!fAlongStepProcessVec[i]->IsApplicable(gtrack))
    //  continue;
    double cStepLimit = fAlongStepProcessVec[i]->AlongStepLimitationLength(gtrack, td);
    if (cStepLimit<limit) {
      limit = cStepLimit;
      gtrack->SetPhysicsProcessIndex(fAlongStepProcessVec[i]->GetIndex());
      gtrack->SetProcess(fAlongStepProcessVec[i]->GetGlobalIndex());       // set global indx of limiting process
    }
  }
  // set fEindex to -1.0 to indicate that continuous step limit happened; will be updated below if not
  gtrack->SetEindex(-1.0);
  //
  // 2. Discrete step limit if any
  // if it has MSC than go to last but one (UNTILL MSC is continuous-discrete)
  size_t maxDiscrtIndx = fPostStepCandidateProcessVec.size();
  if (fIsHasMSCProcess) {
    --maxDiscrtIndx;
  }
  // loop over the discrete processes and get their step limit
  bool haseloss = HasEnergyLossProcess();
  for (size_t i=0; i<maxDiscrtIndx; ++i) {
    //if (!fPostStepCandidateProcessVec[i]->IsApplicable(gtrack))
    //  continue;
    double dStepLimit = fPostStepCandidateProcessVec[i]->PostStepLimitationLength(gtrack, td, haseloss);
    if (dStepLimit<limit) {
      limit = dStepLimit;
      gtrack->SetPhysicsProcessIndex(fPostStepCandidateProcessVec[i]->GetIndex());
      gtrack->SetProcess(fPostStepCandidateProcessVec[i]->GetGlobalIndex());       // set global indx of limiting process
      // flag to indicate that discrete and NOT continous step limit
      gtrack->SetEindex(1000);
    }
  }
  // set the physics step limit (true one)
  gtrack->SetPstep(limit);
}


int PhysicsManagerPerParticle::AlongStepAction(LightTrack &track, Geant::GeantTaskData *td) {
  int numSecondaries = 0;
  for (unsigned long i=0; i<fAlongStepProcessVec.size(); ++i) {
    numSecondaries += fAlongStepProcessVec[i]->AlongStepDoIt(track, td);
    // check if tarck is still alive
    if (track.GetKinE()<=0.0) { // stopped
      // set track status
      if (fAtRestCandidateProcessVec.size()>0) {
        // invoke AtRestAction;
      } else {
        track.SetTrackStatus(LTrackStatus::kKill);
      }
      return numSecondaries;
    }
  }
  return numSecondaries;
}



int PhysicsManagerPerParticle::PostStepAction(LightTrack &track, Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) {
  int numSecondaries = 0;
  // reset number of interaction length left for the dicsrete process that just happened to indicate that it will need
  // to be resample
  size_t physicsProcessIndx = gtrack->GetPhysicsProcessIndex();
  // assert that physicsProcessIndx should never be < 0 at this point
  gtrack->SetPhysicsNumOfInteractLengthLeft(physicsProcessIndx, -1.0);
  // get material-cuts, kinetic energy and pre-step mfp of the selected process
  const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(track.GetMaterialCutCoupleIndex());
  double ekin          = track.GetKinE();
  double mass          = track.GetMass();
  double preStepLambda = gtrack->GetPhysicsInteractLength(physicsProcessIndx); // pre-step mfp
  PhysicsProcess *proc = fProcessVec[physicsProcessIndx];
  if (HasEnergyLossProcess()) {
    // get the current i.e. for the post-step kinetic energy 1/lambda and compare to the pre-step i.e. overestimated
    // 1/lambda to see if it is just a delta interaction and return immediately with null process if yes
    double curMacrXsec = proc->GetMacroscopicXSection(matCut, ekin, mass); // current 1/lambda
    // preStepLambda is lambda and not 1/lambda
    if (curMacrXsec<=0.0 || td->fRndm->uniform()>curMacrXsec*preStepLambda) { // delta interaction
      proc = nullptr;
    }
  }
  // if proc!=nullptr the selected process (discrete part is invoked)
  // should be done like this
//  if (!proc || !proc->IsApplicable(track)) {
//    return numSecondaries;
//  }
  if (!proc) {
    return numSecondaries;
  }
  // invoke the post step action of the selected discrete process
  numSecondaries = proc->PostStepDoIt(track, td);
  return numSecondaries;
}


}  // namespace geantphysics
