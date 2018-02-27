
#include "Geant/PhysicsList.h"

#include "Geant/PhysicsProcess.h"
#include "Geant/Particle.h"
#include "Geant/PhysicsParameters.h"

namespace geantphysics {

PhysicsList::PhysicsList(const std::string &name) : fName(name), fPhysicsParameters(nullptr) {
  fPhysicsParameters = new PhysicsParameters();
}

PhysicsList::~PhysicsList() {}


// user must implement!
void PhysicsList::Initialize( /* Not defined yet */ ) {
  // SIGNATURE (INPUT PARAMETERS AND RETURN TYPE) TO BE REDEFINED LATER...
}

void PhysicsList::AddProcessToParticle(Particle *particle, PhysicsProcess *process) {
  (particle->GetPhysicsProcessVector()).push_back(process);
}


}  // namespace geantphysics
