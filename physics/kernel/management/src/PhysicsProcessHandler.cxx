
#include "Geant/PhysicsProcessHandler.h"

#include "Geant/SystemOfUnits.h"

// realphysics material
#include "Geant/Material.h"
#include "Geant/Element.h"
#include "Geant/Region.h"
#include "Geant/MaterialCuts.h"

#include "Geant/PhysicsListManager.h"
#include "Geant/PhysicsProcess.h"
#include "Geant/PhysicsList.h"
#include "Geant/PhysicsList1.h"

#include "Geant/PhysicsManagerPerParticle.h"
#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

#include "Geant/Propagator.h"
#include "Geant/TaskData.h"

#include "Geant/PhysicsParameters.h"

#include <iostream>

namespace geantphysics {

PhysicsProcessHandler::~PhysicsProcessHandler()
{
  PhysicsListManager::Instance().ClearAll();
  Material::ClearAllMaterials(); // will delete all Element and Isotope as well
  MaterialCuts::ClearAll();      // delete all MaterialCuts objects
  PhysicsData::ClearAll();       // delete all PhysicsData objects
}

void PhysicsProcessHandler::AttachUserData(TaskData *td)
{
  //
  // attach physics data to tasks. This is called after the initialization of
  // physics and after the task data structures get created.
  td->fPhysicsData = new PhysicsData();
}

void PhysicsProcessHandler::Initialize()
{
  //
  // create all MaterialCuts
  MaterialCuts::CreateAll();
  //
  // print out the MaterialCuts-table
  size_t numMatCuts = MaterialCuts::GetTheMaterialCutsTable().size();
  std::cerr << " ====  Number of material-cuts = " << numMatCuts << std::endl;
  // std::cout<<MaterialCuts::GetTheMaterialCutsTable();
  //
  // print out the Material-table
  size_t numMats = Material::GetTheMaterialTable().size();
  std::cerr << " ====  Number of materials = " << numMats << std::endl;
  // std::cout<<Material::GetTheMaterialTable();

  //
  // set number of regions in the PhysicsListManager before we execute the user RegisterPhysicsList method(s)
  PhysicsListManager::Instance().SetNumberOfRegions(vecgeom::Region::GetNumberOfRegions());
  //
  //  Construct one user physics list and register it in the PhysicsListManager: this should be done in the application.
  //
  //  THIS IS VERY SIMILAR To Geant4 STYLE:
  //  We have only one physics list and the active region list vector is not provided: this only one physics list will
  //  be set to used(be active) in all regions automatically.
  //
  // this is what the user will need to do in their own application
  if (!PhysicsListManager::Instance().GetNumberOfRegisteredPhysicsLists()) {
    PhysicsList *physList1 = new PhysicsList1("Physics-List-1");
    PhysicsListManager::Instance().RegisterPhysicsList(physList1);
  }
  //
  // after the user has created and registered their PhysicsList(s): build them and initialize the whole physics
  PhysicsListManager::Instance().BuildPhysicsLists();
  //
  // print out the PhysicsParameters obejct: we ha only one physics list so we have only one physics parameter object
  std::cout << PhysicsParameters::GetThePhysicsParametersTable()[0];
  //
  // The physics data per task is connected when the framework calls AttachUserData later
}

} // namespace geantphysics
