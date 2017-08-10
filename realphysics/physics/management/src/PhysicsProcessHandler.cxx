
#include "PhysicsProcessHandler.h"

#include "SystemOfUnits.h"

// realphysics material
#include "Material.h"
#include "Element.h"
#include "Region.h"
#include "MaterialCuts.h"

#include "PhysicsListManager.h"
#include "PhysicsProcess.h"
#include "PhysicsList.h"
#include "PhysicsList1.h"

#ifdef USE_VECPHYS
#include "VecPhysList.h"
#endif

#include "PhysicsManagerPerParticle.h"
#include "LightTrack.h"
#include "PhysicsData.h"

#include "GeantPropagator.h"
#include "GeantTaskData.h"

#include "PhysicsParameters.h"

#include <iostream>

namespace geantphysics {

PhysicsProcessHandler::PhysicsProcessHandler() {
 // nothing to do so far
}


PhysicsProcessHandler::~PhysicsProcessHandler() {
  PhysicsListManager::Instance().ClearAll();
  Material::ClearAllMaterials(); // will delete all Element and Isotope as well
  MaterialCuts::ClearAll();  // delete all MaterialCuts objects
  PhysicsData::ClearAll();   // delete all PhysicsData objects
}

void PhysicsProcessHandler::AttachUserData(GeantTaskData *td) {
  //
  // attach physics data to tasks. This is called after the initialization of
  // physics and after the task data structures get created.
  td->fPhysicsData = new PhysicsData();
}

void PhysicsProcessHandler::Initialize() {
  //
  // create all MaterialCuts
  MaterialCuts::CreateAll();
  //
  // print out the MaterialCuts-table
  std::cout<<MaterialCuts::GetTheMaterialCutsTable();
  //
  // print out the MaterialCuts-table
  std::cout<<Material::GetTheMaterialTable();
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
#ifdef USE_VECPHYS
      PhysicsList *vecphysList = new VecPhysList("VecPhysList");
      PhysicsListManager::Instance().RegisterPhysicsList(vecphysList);
#else
      if (!PhysicsListManager::Instance().GetNumberOfRegisteredPhysicsLists()) {
        PhysicsList *physList1 = new PhysicsList1("Physics-List-1");
        PhysicsListManager::Instance().RegisterPhysicsList(physList1);
      }
#endif
  //
  // after the user has created and registered their PhysicsList(s): build them and initialize the whole physics
  PhysicsListManager::Instance().BuildPhysicsLists();
  //
  // print out the PhysicsParameters obejct: we ha only one physics list so we have only one physics parameter object
  std::cout<<PhysicsParameters::GetThePhysicsParametersTable()[0];
  //
  // The physics data per task is connected when the framework calls AttachUserData later
}


} // namespace geantphysics
