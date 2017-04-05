
#include "PhysicsProcess.h"

#include "LightTrack.h"
#include "Particle.h"

#include <iostream>

namespace geantphysics {

std::vector<PhysicsProcess*> PhysicsProcess::gThePhysicsProcessTable;

const double PhysicsProcess::gAVeryLargeValue = 1.0e20;

PhysicsProcess::PhysicsProcess(const std::string &aName)
: fIndex(-1), fIsDiscrete(false), fIsContinuous(false), fIsAtRest(false),
  fForcedCondition(ForcedCondition::kNotForced), fType(ProcessType::kNotDefined),
  fName(aName) {
  fIndex = gThePhysicsProcessTable.size();
  gThePhysicsProcessTable.push_back(this);
}


PhysicsProcess::PhysicsProcess(const bool aIsDiscrete, const bool aIsContinuous,
                               const bool aIsAtRest, const ForcedCondition aForcedCondition,
                               const ProcessType aType, const std::string &aName)
: fIndex(-1), fIsDiscrete(aIsDiscrete), fIsContinuous(aIsContinuous), fIsAtRest(aIsAtRest),
  fForcedCondition(aForcedCondition), fType(aType), fName(aName) {
  fIndex = gThePhysicsProcessTable.size();
  gThePhysicsProcessTable.push_back(this);
}


/*
PhysicsProcess::PhysicsProcess(const PhysicsProcess &other) :
  fIsDiscrete(other.fIsDiscrete), fIsContinuous(other.fIsContinuous),
  fIsAtRest(other.fIsAtRest),
  fForcedCondition(other.fForcedCondition), fType(other.fType),
  fName(other.fName), fListActiveRegions(other.fListActiveRegions)
{}


PhysicsProcess& PhysicsProcess::operator=(const PhysicsProcess &other) {
  if (this != &other) {
    fIsDiscrete = other.fIsDiscrete;
    fIsContinuous = other.fIsContinuous;
    fIsAtRest = other.fIsAtRest;
    fForcedCondition = other.fForcedCondition;
    fType = other.fType;
    fName = other.fName;
    fListActiveRegions = other.fListActiveRegions;
  }
  return *this;
}
*/

PhysicsProcess::~PhysicsProcess() {}


void PhysicsProcess::Initialize() {
  // check if the process is assigned only to allowed particles
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

/*
double PhysicsProcess::GetAtomicCrossSection(const LightTrack &track) const {
  int particleCode = track.GetGVcode();
  double particleKinE = track.GetKinE();
  int targetZ = track.GetTargetZ();
  int targetN = track.GetTargetN();
  return GetAtomicCrossSection(particleCode, particleKinE, targetZ, targetN);
}
*/

/*
double PhysicsProcess::InverseLambda(const LightTrack &track) const {
  return 0.0;
}
*/

double PhysicsProcess::AlongStepLimitationLength(const LightTrack & /*track*/) const {
  return gAVeryLargeValue;
}

double PhysicsProcess::AverageLifetime(const LightTrack & /*track*/) const {
  return gAVeryLargeValue;
}


void PhysicsProcess::ClearAllProcess() {
  for (unsigned long i=0; i<gThePhysicsProcessTable.size(); ++i) {
    std::cerr<<"  ********************* deleting proc = "<<gThePhysicsProcessTable[i]->GetName()<<std::endl;
    delete gThePhysicsProcessTable[i];
  }
  gThePhysicsProcessTable.clear();
}

} // namespace geantphysics
