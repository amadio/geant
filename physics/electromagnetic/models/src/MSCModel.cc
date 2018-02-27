
#include "MSCModel.h"

// from geantV
#include "TaskData.h"
#include "Geant/Track.h"

namespace geantphysics {


MSCModel::MSCModel(const std::string& name) : EMModel(name), fRangeFactor(0.06), fSafetyFactor(0.6), fGeomFactor(2.5)
 , fSkin(3.), fMSCSteppingAlgorithm(MSCSteppingAlgorithm::kUseSaftey) {}

MSCModel::~MSCModel() {}


void MSCModel::Initialize() {
  // call the base class Initialize method
  EMModel::Initialize();
}


} // namespace geantphysics
