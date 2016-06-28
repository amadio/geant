#include "PhysicsManagerPerParticle.h"
#include "NewPhysicsProcess.h"

using namespace geant;

//----------------------------------------------------
//--- PhysicsManagerPerParticle non-inline methods ---
//----------------------------------------------------


PhysicsManagerPerParticle::PhysicsManagerPerParticle() : 
  fParticleGVcode( 0 )
{}


PhysicsManagerPerParticle::PhysicsManagerPerParticle( const int aParticleGVcode ) :
  fParticleGVcode( aParticleGVcode )
{}


PhysicsManagerPerParticle::~PhysicsManagerPerParticle() {
  // NOT CLEAR WHO SHOULD DELETE THE PROCESSES: FOR THE TIME BEING, 
  // NOTHING IS DONE IN THE DESTRUCTOR.
}


void PhysicsManagerPerParticle::AddProcess( PhysicsProcess *process ) {
  if ( process == nullptr || 
       ! process->IsApplicable( fParticleGVcode ) ||
       process->GetForcedCondition() == ForcedCondition::kInActivated ) return;

  bool isForced = ( process->GetForcedCondition() != ForcedCondition::kNotForced );

  if ( process->GetIsContinuous() ) fAlongStepProcessVec.push_back( process );

  if ( process->GetIsDiscrete() ) {
    fPostStepCandidateProcessVec.push_back( process );
    if ( isForced ) fPostStepForcedProcessVec.push_back( process );
  }

  if ( process->GetIsAtRest() ) {
    fAtRestCandidateProcessVec.push_back( process );
    if ( isForced ) fAtRestForcedProcessVec.push_back( process );
  }
}


void PhysicsManagerPerParticle::BuildTotalLambdaTables( /* Not defined yet */ ) {
  // SIGNATURE (INPUT PARAMETERS AND RETURN TYPE) TO BE REDEFINED LATER...
  // IMPLEMENTATION TO BE WRITTEN...
}


void PhysicsManagerPerParticle::GetTotalLambdaTable( /* Not defined yet */ ) const {
  // SIGNATURE (INPUT PARAMETERS AND RETURN TYPE) TO BE REDEFINED LATER...
  // IMPLEMENTATION TO BE WRITTEN...
}

