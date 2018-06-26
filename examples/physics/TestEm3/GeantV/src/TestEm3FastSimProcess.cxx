#include "TestEm3FastSimProcess.h"
#include "Geant/LightTrack.h"
#include "Geant/Isotope.h"
#include "Geant/Material.h"
#include "Geant/MaterialCuts.h"
#include "Geant/MaterialProperties.h"
#include "Geant/Particle.h"

#include "Geant/Electron.h"
#include "Geant/Gamma.h"

#include <iostream>

namespace userapplication {

//-----------------------------------
// TestEm3FastSimProcess non-inline methods
//-----------------------------------

TestEm3FastSimProcess::TestEm3FastSimProcess() : geantphysics::FastSimProcess("TestEm3FastSim")
{
  AddToListParticlesAlloedToAssigned(geantphysics::Electron::Definition());
  AddToListParticlesAlloedToAssigned(geantphysics::Gamma::Definition());
}

TestEm3FastSimProcess::TestEm3FastSimProcess(const std::vector<int> &particlecodevec)
    : geantphysics::FastSimProcess("TestEm3FastSim")
{
  SetParticleCodeVec(particlecodevec);
}

TestEm3FastSimProcess::~TestEm3FastSimProcess() {}

bool TestEm3FastSimProcess::IsApplicable(geant::Track *track) const
{

  bool isOK = false;

  if (track->E() > 5.0) isOK = true;

  return isOK;
}

int TestEm3FastSimProcess::FastSimDoIt(geantphysics::LightTrack &track, geant::TaskData *td)
{
  std::cerr << "****** TestEm3FastSimDoIt called for "
            << geantphysics::Particle::GetParticleByInternalCode(track.GetGVcode())->GetName() << " with energy "
            << track.GetKinE() << " GeV" << std::endl;
  return 0;
}
} // namespace userapplication
