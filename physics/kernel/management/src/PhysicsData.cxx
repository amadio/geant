
#include "PhysicsData.h"

#include "LightTrack.h"

namespace geantphysics {

std::vector<PhysicsData*> PhysicsData::gThePhysicsDataTable;

PhysicsData::PhysicsData() : fNumUsedSecondaries(0) {
  const int initSizeListOfSecondaries = 2;
  fListOfSecondaries.resize(initSizeListOfSecondaries);
  gThePhysicsDataTable.push_back(this);
}

void PhysicsData::ClearAll() {
  for (unsigned int i=0; i<gThePhysicsDataTable.size(); ++i) {
    delete gThePhysicsDataTable[i];
  }
  gThePhysicsDataTable.clear();
}

} // namespace geantphysics
