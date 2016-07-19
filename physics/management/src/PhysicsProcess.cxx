#include "NewPhysicsProcess.h"
#include "LightTrack.h"

using namespace geant;

//-----------------------------------------
//--- PhysicsProcess non-inline methods ---
//-----------------------------------------

const double aVeryLargeValue = 1.0e20;

 
PhysicsProcess::PhysicsProcess() : 
  fIsDiscrete( false ), fIsContinuous( false ), fIsAtRest( false ),
  fForcedCondition( ForcedCondition::kNotForced ), fType( ProcessType::kNotDefined ),
  fName( "" ) 
{}


PhysicsProcess::PhysicsProcess( const bool aIsDiscrete, const bool aIsContinuous, 
                                const bool aIsAtRest, const ForcedCondition aForcedCondition,
                                const ProcessType aType, const std::string &aName ) :
  fIsDiscrete( aIsDiscrete ), fIsContinuous( aIsContinuous ), fIsAtRest( aIsAtRest ),
  fForcedCondition( aForcedCondition ), fType( aType ), fName( aName ) 
{}


PhysicsProcess::PhysicsProcess( const PhysicsProcess &other ) : 
  fIsDiscrete( other.fIsDiscrete ), fIsContinuous( other.fIsContinuous ), 
  fIsAtRest( other.fIsAtRest ), 
  fForcedCondition( other.fForcedCondition ), fType( other.fType ),
  fName( other.fName )
{}


PhysicsProcess& PhysicsProcess::operator=( const PhysicsProcess &other ) {
  if ( this != &other ) {
    fIsDiscrete = other.fIsDiscrete;
    fIsContinuous = other.fIsContinuous;
    fIsAtRest = other.fIsAtRest;
    fForcedCondition = other.fForcedCondition;
    fType = other.fType;
    fName = other.fName;
  }
  return *this;
}


PhysicsProcess::~PhysicsProcess() {}


double PhysicsProcess::GetAtomicCrossSection( const LightTrack &track ) const {
  int particleCode = track.GetGVcode();
  double particleKinE = track.GetKinE();
  int targetZ = track.GetTargetZ();
  int targetN = track.GetTargetN();
  return GetAtomicCrossSection( particleCode, particleKinE, targetZ, targetN );
}

 
double PhysicsProcess::InverseLambda( const LightTrack &/*track*/ ) const {
  return 0.0;
}

double PhysicsProcess::AlongStepLimitationLength( const LightTrack &/*track*/ ) const {
  return aVeryLargeValue;
}
 
double PhysicsProcess::AverageLifetime( const LightTrack &/*track*/ ) const {
  return aVeryLargeValue;
}


void PhysicsProcess::SampleTarget( LightTrack &/*track*/ ) const {
  //
  // IMPLEMENTATION IS NEEDED HERE: IT SHOULD BE THE SAME FOR ALL PROCESSES!
  //  
}

