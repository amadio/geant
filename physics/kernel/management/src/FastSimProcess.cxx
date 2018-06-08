#include "Geant/FastSimProcess.h"
#include "Geant/LightTrack.h"
#include "Geant/Isotope.h"
#include "Geant/Material.h"
#include "Geant/MaterialCuts.h"
#include "Geant/MaterialProperties.h"
#include "Geant/Particle.h"

using namespace geantphysics;

//-----------------------------------
// FastSimProcess non-inline methods
//-----------------------------------

FastSimProcess::FastSimProcess()
    : PhysicsProcess("")
{
  SetFastSim(true);
}

FastSimProcess::FastSimProcess(const std::string &name) : PhysicsProcess(name)
{
  SetFastSim(true);
}

FastSimProcess::FastSimProcess(const std::string &name, const std::vector<int> &particlecodevec)
    : PhysicsProcess(false, false, false, ForcedCondition::kNotForced, ProcessType::kFastSim, name)
{
  SetFastSim(true);  
  SetParticleCodeVec(particlecodevec);
}

FastSimProcess::~FastSimProcess()
{
}

bool FastSimProcess::IsApplicable(geant::Track *track) const
{
  return true;
}

int FastSimProcess::FastSimDoIt(LightTrack &track, geant::TaskData *td)
{
  return 0;
}
