
#include "PhysicsProcess.h"

#include "LightTrack.h"
#include "Particle.h"

#include "LambdaTable.h"

#include <iostream>

namespace geantphysics {

std::vector<PhysicsProcess*> PhysicsProcess::gThePhysicsProcessTable;

const double PhysicsProcess::gAVeryLargeValue = 1.0e20;

PhysicsProcess::PhysicsProcess(const std::string &aName)
: fIndex(-1), fGlobalIndex(-1), fIsDiscrete(false), fIsContinuous(false), fIsAtRest(false), fIsLambdaTableRequested(false),
  fForcedCondition(ForcedCondition::kNotForced), fType(ProcessType::kNotDefined),
  fName(aName), fPhysicsParameters(nullptr), fParticle(nullptr), fLambdaTable(nullptr) {
  fGlobalIndex = gThePhysicsProcessTable.size();
  gThePhysicsProcessTable.push_back(this);
}


PhysicsProcess::PhysicsProcess(const bool aIsDiscrete, const bool aIsContinuous,
                               const bool aIsAtRest, const ForcedCondition aForcedCondition,
                               const ProcessType aType, const std::string &aName)
: fIndex(-1), fGlobalIndex(-1), fIsDiscrete(aIsDiscrete), fIsContinuous(aIsContinuous), fIsAtRest(aIsAtRest), fIsLambdaTableRequested(false),
  fForcedCondition(aForcedCondition), fType(aType), fName(aName), fPhysicsParameters(nullptr), fParticle(nullptr),
  fLambdaTable(nullptr) {
  fGlobalIndex = gThePhysicsProcessTable.size();
  gThePhysicsProcessTable.push_back(this);
}


PhysicsProcess::~PhysicsProcess() {
  // delete lambda tables if any
  if (fLambdaTable) {
    delete fLambdaTable;
  }
}

// the PhysicsParameters pointer is set by the PhysicsManagerPerParticle before calling this Initialize method
void PhysicsProcess::Initialize() {
  // check if the process is assigned only to allowed particles
  std::cerr<<"  ----> PhysicsProcess   Name = " << GetName() << "  is under initialization! "<< std::endl;
  for (unsigned long i=0; i<fListParticlesAssignedTo.size(); ++i) {
    Particle* part = fListParticlesAssignedTo[i];
    bool isok = false;
    for (unsigned long j=0; j<fListParticlesAlloedToAssigned.size(); ++j) {
      if (part==fListParticlesAlloedToAssigned[j]) {
        isok = true;
      }
    }
    if (!isok) {
      std::cerr << " *** ERROR: PhysicsProcess::Initialise()\n"
                << "   Process with Name = " << GetName() << "\n"
                << "   Is assigned to particle with Name = " << part->GetName() << "\n"
                << "   that is not in the allowed particle list of the process!"
                << std::endl;
      exit(-1);
    }
  }
}


double PhysicsProcess::AlongStepLimitationLength(Geant::GeantTrack * /*track*/, Geant::GeantTaskData * /*td*/) const {
  return gAVeryLargeValue;
}


double PhysicsProcess::PostStepLimitationLength(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td, bool haseloss) {
  double stepLimit = GetAVeryLargeValue();
  // get the material-cuts and kinetic energy
  const MaterialCuts *matCut = static_cast<const MaterialCuts*>((const_cast<vecgeom::LogicalVolume*>(gtrack->GetVolume())->GetMaterialCutsPtr()));
  double ekin                = gtrack->T();
  double mass                = gtrack->Mass(); // dynamic mass of the particle
  // get/compute the mean free path by (1)getting/(2)computing the macroscopic scross section by Accounting Possible
  // Energy Losses along the step:
  // - (1) from lambda table requested to built by the process
  // - (2) or by calling the ComputeMacroscopicXSection interface method directly if there is no lambda-table
  double macrXsec = GetMacroscopicXSectionForStepping(matCut, ekin, mass, haseloss);
  double mfp = GetAVeryLargeValue();
  if (macrXsec>0.) {
    mfp = 1./macrXsec;
  }
  // check if we need to sample new num.-of-int.-length-left: it is updated either after propagation(particles that
  // doesn't have MSC) or after the post-propagation (particles that has MSC)
  if (gtrack->GetPhysicsNumOfInteractLengthLeft(GetIndex())<=0.0) {
    double rndm = td->fRndm->uniform(); // use vecgeom RNG to get uniform random number
    gtrack->SetPhysicsNumOfInteractLengthLeft(GetIndex(), -std::log(rndm));
  }
  // save the mfp: to be used for the update num.-of-int.-length-left
  gtrack->SetPhysicsInteractLength(GetIndex(), mfp);
  //update the step length => length = lambda * -1. * log(rndm) = lambda * number of interaction leght left;
  stepLimit = mfp*gtrack->GetPhysicsNumOfInteractLengthLeft(GetIndex());
  return stepLimit;
}


double PhysicsProcess::AverageLifetime(const LightTrack & /*track*/) const {
  return gAVeryLargeValue;
}


void PhysicsProcess::RequestLambdaTables(bool ispermaterial) {
  fIsLambdaTableRequested = true;
  if (fLambdaTable) {
    delete fLambdaTable;
  }
  fLambdaTable = new LambdaTable(this, ispermaterial);
  // will be built by the ProcessManagerPerParticle after the process is initialized.
}


void PhysicsProcess::SetSpecialLambdaTableBinNum(int val) {
  if (!fLambdaTable) {
    RequestLambdaTables(); // by default it will be per-material !
  }
  fLambdaTable->SetSpecialLambdaTableBinNum(val);
}


void PhysicsProcess::BuildLambdaTables() {
  fLambdaTable->BuildLambdaTables();
}


double PhysicsProcess::GetMacroscopicXSectionMaximumEnergy(const MaterialCuts *matcut) {
  double maxOfMacXsecE = gAVeryLargeValue;
  if (fLambdaTable) {
    maxOfMacXsecE = fLambdaTable->GetMacroscopicXSectionMaximumEnergy(matcut);
  } else {
    maxOfMacXsecE = MacroscopicXSectionMaximumEnergy(matcut);
  }
  return maxOfMacXsecE;
}

// will be called only if GetMacroscopicXSectionMaximumEnergy < gAVeryLargeValue
double PhysicsProcess::GetMacroscopicXSectionMaximum(const MaterialCuts *matcut) {
  double maxOfMacXsec = 0.;
  if (fLambdaTable) {
    maxOfMacXsec = fLambdaTable->GetMacroscopicXSectionMaximum(matcut);
  } else {
    maxOfMacXsec = MacroscopicXSectionMaximum(matcut);
  }
  return maxOfMacXsec;
}


double PhysicsProcess::GetMacroscopicXSection(const MaterialCuts *matcut, double ekin, double mass) {
  double macrXsec = 0.;
  // Get the macroscopic cross section form the lambda table if it was requested to be built by the process or
  // call the ComputeMacroscopicXSection interface method to compute it on-the-fly
  if (fLambdaTable) {
    macrXsec = fLambdaTable->GetMacroscopicXSection(matcut, ekin);
  } else {
    macrXsec = ComputeMacroscopicXSection(matcut, ekin, fParticle, mass);
  }
  if (macrXsec<0.) {
    macrXsec =0.;
  }
  return macrXsec;
}


// called only at the pre-step point
double PhysicsProcess::GetMacroscopicXSectionForStepping(const MaterialCuts *matcut, double ekin, double mass, bool haseloss) {
  double macrXsec = 0.;
  // account possible energy loss along the step if the particle has energy loss process(es)
  if (haseloss) {
    // get the kinetic energy at which the macroscopic cross section has its maximum
    double maxOfMacXsecE = GetMacroscopicXSectionMaximumEnergy(matcut);
    // if the current kinetic energy is already on the left side of this maximum we provide 1/lambda for that just
    // because we assume that 1/lambda is already decreasing on this side with decareasing energy so we provide an
    // overestimate of 1/lambda
    // if the current kinetic energy is on the right side of this maximum point: more work to give an overestimate:
    if (ekin>maxOfMacXsecE) {
      // compute reduced energy: we assume that 1/lambda is higher at lower energy so we provide an overestimate
      double ekinReduced = 0.8*ekin;
      // check if it is still on the right side of the maximum point i.e. if our assumption is fine
      // if not: the reduced energy got to the left side of the maximum so we jumped the maximum so set the macroscopic
      // cross section to its maximum value: note that GetMacroscopicXSectionMaximumEnergy will return (by default) with
      // a very large value (=>we never get here) or with a proper value if lambda table was requested to build
      if (ekinReduced<maxOfMacXsecE) {
        macrXsec = GetMacroscopicXSectionMaximum(matcut);
        if (macrXsec<0.) {
          macrXsec = 0.;
        }
        return macrXsec;
      } else {
        // otherwise we are still on the right side of the maximum so provide 1/lambda at this reduced energy
        ekin = ekinReduced;
      }
    }
    // if we did not return earlier then we need to provide 1/lambda at ekin that has been set properly above to ensure
    // that the 1/lambda value at ekin energy will be higher than any 1/lambda values along the step i.e. between the
    // current, pre-step and post-step point energy:
  }
  //
  macrXsec = GetMacroscopicXSection(matcut, ekin, mass);
  return macrXsec;
}


void PhysicsProcess::ClearAllProcess() {
  for (unsigned long i=0; i<gThePhysicsProcessTable.size(); ++i) {
    std::cerr<<"  ********************* deleting proc = "<<gThePhysicsProcessTable[i]->GetName()<<std::endl;
    delete gThePhysicsProcessTable[i];
  }
  gThePhysicsProcessTable.clear();
}

} // namespace geantphysics
