
#include "PhysicsProcessHandler.h"

#include "SystemOfUnits.h"

// vecgeom materials
//#ifdef USE_VECGEOM_NAVIGATOR
#include "materials/Material.h"
//#endif

// realphysics material
#include "Material.h"
#include "Element.h"
#include "Region.h"
#include "MaterialCuts.h"

#include "PhysicsListManager.h"
#include "PhysicsList.h"
#include "PhysicsList1.h"

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


void PhysicsProcessHandler::Initialize() {
  //
  // see the comments at its implementation
  BuildMaterials();
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
      PhysicsList *physList1 = new PhysicsList1("Physics-List-1");
      PhysicsListManager::Instance().RegisterPhysicsList(physList1);
  //
  // after the user has created and registered their PhysicsList(s): build them and initialize the whole physics
  PhysicsListManager::Instance().BuildPhysicsLists();
  //
  // print out the PhysicsParameters obejct: we ha only one physics list so we have only one physics parameter object
  std::cout<<PhysicsParameters::GetThePhysicsParametersTable()[0];
  //
  // THIS is the place where we should construct as many PhysicsData objects as working threads and set their pointers
  // in the corresponding GeantTaskData objects. HOWEVER, unfortunately those GeantTaskData objects are constructedhere
  // at the moment only after the call to this Initialize() method (so we need to do it in the Application::Initialize()
  // method at the moment.
  //
  // GeantTaskData are constructed later
  //
//  for (int i=0; i<GeantPropagator::Instance()->fNthreads; ++i) {
//    GeantPropagator::Instance()->fThreadData[i]->fPhysicsData = new geantphysics::PhysicsData();
//  }
}


//
//  This is a helper method to create proper geantphysics::Materials and one geantphysics::Region that are
//  esential for having a real physics simulation.
//
void PhysicsProcessHandler::BuildMaterials() {
  // Load elements from geometry, however in most cases it should already be done
 // vecgeom materials To geantphysics:;Material and creating 1 region and adding all material to that
 #ifdef USE_VECGEOM_NAVIGATOR
   std::vector<vecgeom::Material *> matlist = vecgeom::Material::GetMaterials();
   if (matlist.size()==0) {
     std::cerr<<"  ***** ERROR PhysicsProcessHandler::BuildMaterials(): \n"
              <<"   No material was found! (vecgeom::Material::GetMaterials() returns with 0 size list)\n"
              <<std::endl;
     exit(-1);
   }
   for (unsigned long i=0; i<matlist.size(); ++i) {
     std::cout<<"     -->  Creating Material with index = " << i <<" and vecgeom::name: "<<matlist[i]->GetName();
     vecgeom::Material *vgMat = matlist[i];
     int    numElem    = vgMat->GetNelements();
     double density    = vgMat->GetDensity()*geant::g/geant::cm3; // in g/cm3
     const std::string  name = vgMat->GetName();
     // check if it is a G4 NIST material
     std::string postName = "";
     bool isNistMaterial = false;
     if (name.substr(0,3)=="G4_") {
       postName = name.substr(3);
       isNistMaterial = true;
     }
     Material *matX = nullptr;
     if (isNistMaterial) {
       std::string nistName = "NIST_MAT_"+postName;
       matX = Material::NISTMaterial(nistName);
     } else {
       // create material
       matX = new Material(name, density, numElem);
       for (int j=0; j<numElem; ++j) {
         double va;
         double vz;
         double vw;
         vgMat->GetElementProp(va, vz, vw, j);
         // create NIST element
         Element *elX = Element::NISTElement(vz);
         // add to the Material
         matX->AddElement(elX, vw);
      }
    }
    std::cout<< "  geantphysics::name = " << matX->GetName() << std::endl;
   }
   std::cout<<" ================================================================ \n";
 #else
  std::cerr<<"  ***** ERROR PhysicsProcessHandler::BuildMaterials(): \n"
           <<"   RealPhysics support only VecGeom so build with -DUSE_VECGEOM_NAVIGATOR=ON\n"
           <<std::endl;
  exit(-1);
 #endif
}


void PhysicsProcessHandler::ComputeIntLen(Material_t *mat, int ntracks, GeantTrack_v &tracks, double * /*lengths*/,
                                          GeantTaskData *td) {
  for (int i=0; i<ntracks; ++i) {
    // here we will get the MaterialCuts from the LogicalVolume later
    int   matIndx = tracks.GetMaterial(i)->GetIndex();
    int   regIndx = const_cast<vecgeom::LogicalVolume*>(tracks.GetVolume(i))->GetRegion()->GetIndex();
    const MaterialCuts *matCut =  MaterialCuts::GetMaterialCut(regIndx,matIndx);
  /*
    if (mat) {
      matCut = MaterialCuts::GetMaterialCut(mat->GetIndex());
    } else {
      matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterial(i)->GetIndex());
    }
  */
    int   particleCode       = tracks.fGVcodeV[i];
    const Particle *particle = Particle::GetParticleByInteralCode(particleCode);
    // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
    if (!pManager) {
      tracks.fPstepV[i]  = PhysicsProcess::GetAVeryLargeValue();
      tracks.fIntLenV[i] = 1.0;
      continue; // no physics limit
    }
    // check if the partcile has anything along step and if yes get the minimum of their along-step limit
    // create a LightTrack
    if (pManager->GetListAlongStepProcesses().size()+pManager->GetListPostStepCandidateProcesses().size()>0) {
       LightTrack primaryLT;
       //LightTrack *lt = &(td->fPhysicsData->GetListOfSecondaries()[0]);
       // we will use members:
       //  fNintLen      <==>  fNintLen  // number of interaction left
       //  fTargetZ      <==>  fEindex   // will be set to flag if disc. or cont. step limit won
       //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
       //  fKinE         <==>  fE-fMass  // kinetic energy
       //  fGVcode       <==>  fGVcode   // internal particle code
       //  fIntLen       <==>  fIntLen   // will be set to store the current inverse total lambda
       //  fStepLength   <==>  fPstep    // will be set to store the physics step limit
       primaryLT.SetNumOfInteractionLegthLeft(tracks.fNintLenV[i]);
       primaryLT.SetKinE(tracks.fEV[i]-tracks.fMassV[i]);
       primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
       primaryLT.SetGVcode(tracks.fGVcodeV[i]);
       pManager->ComputeIntLen(primaryLT,td);
       // update GeantTrack
       tracks.fPstepV[i]   = primaryLT.GetStepLength();
       tracks.fIntLenV[i]  = primaryLT.GetTotalMFP();
       tracks.fEindexV[i]  = primaryLT.GetTargetZ();
       if (tracks.fNintLenV[i]<=0.0) { // was resampled
         tracks.fNintLenV[i] = primaryLT.GetNumOfInteractionLegthLeft();
       }
    }
  }
}


void PhysicsProcessHandler::AlongStepAction(Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout,
                                            GeantTaskData * /*td*/) {
  int numSecondaries = 0;
  for (int i=0; i<ntracks; ++i) {
    // here we will get the MaterialCuts from the LogicalVolume later
    int   matIndx = tracks.GetMaterial(i)->GetIndex();
    int   regIndx = const_cast<vecgeom::LogicalVolume*>(tracks.GetVolume(i))->GetRegion()->GetIndex();
    const MaterialCuts *matCut =  MaterialCuts::GetMaterialCut(regIndx,matIndx);
  /*
    const MaterialCuts *matCut = nullptr;
    if (mat) {
      matCut = MaterialCuts::GetMaterialCut(mat->GetIndex());
    } else {
      matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterial(i)->GetIndex());
    }
  */
    int particleCode         = tracks.fGVcodeV[i];
    const Particle *particle = Particle::GetParticleByInteralCode(particleCode);
    // check if the partcile has anything along step
    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
    if (!pManager) {
      continue; // no physics of this particle in this region
    }
    if (pManager->GetListAlongStepProcesses().size()>0) {
       LightTrack primaryLT;
       // we will use members:
       //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
       //  fKinE         <==>  fE-fMass  // kinetic energy; will be set to the new kinetic energy
       //  fGVcode       <==>  fGVcode   // internal particle code
       //  fStepLength   <==>  fStep     // current step length
       //  fEdep         <==>  fEdep     // init to 0.0; will be set to the current energy deposit
       primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
       primaryLT.SetKinE(tracks.fEV[i]-tracks.fMassV[i]);
       primaryLT.SetGVcode(tracks.fGVcodeV[i]);
       primaryLT.SetStepLength(tracks.fStepV[i]);
       primaryLT.SetEnergyDeposit(0.0);
       //primaryLT.SetTrackIndex(i);
       std::vector<LightTrack> secLt; // just dummy: no along step secondary production at the moment
       int nSecParticles = pManager->AlongStepAction(primaryLT,secLt);
       // update GeantTrack
       double newEkin = primaryLT.GetKinE();
       tracks.fEV[i]  = newEkin+tracks.fMassV[i];
       tracks.fPV[i]  = std::sqrt(newEkin*(newEkin+2.0*tracks.fMassV[i]));
       tracks.fEdepV[i] += primaryLT.GetEnergyDeposit();
       if (primaryLT.GetTrackStatus()==LTrackStatus::kKill) {
         tracks.fStatusV[i] = Geant::kKilled;
       }
       numSecondaries += nSecParticles;
       if (nSecParticles) {
         //insert secondary tracks
       }
    }
  }
  nout = numSecondaries; // number of inserted GeantTracks
}


void PhysicsProcessHandler::PostStepAction(Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout,
                                           GeantTaskData *td) {
  int numSecondaries = 0;
  for (int i=0; i<ntracks; ++i) {
    // here we will get the MaterialCuts from the LogicalVolume later
    int   matIndx = tracks.GetMaterial(i)->GetIndex();
    int   regIndx = const_cast<vecgeom::LogicalVolume*>(tracks.GetVolume(i))->GetRegion()->GetIndex();
    const MaterialCuts *matCut =  MaterialCuts::GetMaterialCut(regIndx,matIndx);
  /*
    const MaterialCuts *matCut = nullptr;
    if (mat) {
     matCut = MaterialCuts::GetMaterialCut(mat->GetIndex());
    } else {
     matCut = MaterialCuts::GetMaterialCut(tracks.GetMaterial(i)->GetIndex());
    }
  */
    int particleCode         = tracks.fGVcodeV[i];
    const Particle *particle = Particle::GetParticleByInteralCode(particleCode);
    // check if the partcile has anything along step
    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
    if (!pManager) {
      continue; // no physics of this particle in this reagion
    }
    if (pManager->GetListPostStepCandidateProcesses().size()>0) {
      LightTrack primaryLT;
      // we will use members:
      //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
      //  fKinE         <==>  fE-fMass  // kinetic energy;  will be set to the new kinetic energy
      //  fGVcode       <==>  fGVcode   // internal particle code
      //  fIntLen       <==>  fIntLen   // pre-step lambda for accounting energy loss i.e. to see if it is a delta inter.
      //  fXdir         <==>  fXdir     // direction vector x comp. will be set to the new direction x comp.
      //  fYdir         <==>  fYdir     // direction vector y comp. will be set to the new direction y comp.
      //  fZdir         <==>  fZdir     // direction vector z comp. will be set to the new direction z comp.
      primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
      primaryLT.SetKinE(tracks.fEV[i]-tracks.fMassV[i]);
      primaryLT.SetGVcode(tracks.fGVcodeV[i]);
      primaryLT.SetTrackIndex(i);
      primaryLT.SetDirX(tracks.fXdirV[i]);
      primaryLT.SetDirY(tracks.fYdirV[i]);
      primaryLT.SetDirZ(tracks.fZdirV[i]);
      primaryLT.SetTotalMFP(tracks.fIntLenV[i]);
      std::vector<LightTrack> secLt;  // dummy because we fill secondaries into GeantTaskData::PhysicsData
      //
      // clean the number of secondary tracks used (in PhysicsData)
      td->fPhysicsData->SetNumUsedSecondaries(0);
      //
      // invoke the PostStepAction of this particle PhysicsManagerPerParticle
      int nSecParticles = pManager->PostStepAction(primaryLT, secLt, td);
      //
      // update GeantTrack
      double newEkin    = primaryLT.GetKinE();
      tracks.fEV[i]     = newEkin+tracks.fMassV[i];
      tracks.fPV[i]     = std::sqrt(newEkin*(newEkin+2.0*tracks.fMassV[i]));
      tracks.fXdirV[i]  = primaryLT.GetDirX();
      tracks.fYdirV[i]  = primaryLT.GetDirY();
      tracks.fZdirV[i]  = primaryLT.GetDirZ();
      tracks.fEdepV[i] += primaryLT.GetEnergyDeposit();
      if (primaryLT.GetTrackStatus()==LTrackStatus::kKill) {
        tracks.fStatusV[i] = Geant::kKilled;
      }
      //
      // create secondary tracks if there are any
      if (nSecParticles) {
        for (int isec=0; isec<nSecParticles; ++isec) {
          // get the list of secondary tracks
          secLt = td->fPhysicsData->GetListOfSecondaries();
          int   secGVcode = secLt[isec].GetGVcode(); // GV index of this secondary particle
          const Particle *secParticle = Particle::GetParticleByInteralCode(secGVcode);
          // get a GeantTrack geantTrack;
          Geant::GeantTrack &geantTrack = td->GetTrack();
          // set the new track properties
          int t = secLt[isec].GetTrackIndex();          // parent GeantTrack index in the input GeantTrack_v
          geantTrack.fEvent    = tracks.fEventV[t];
          geantTrack.fEvslot   = tracks.fEvslotV[t];
          geantTrack.fGVcode   = secGVcode;
          geantTrack.fEindex   = 0;
          geantTrack.fCharge   = secParticle->GetPDGCharge();
          geantTrack.fProcess  = 0;
          geantTrack.fNsteps   = 0;
          geantTrack.fStatus   = Geant::kNew;           // this secondary is a new track
          geantTrack.fMass     = secParticle->GetPDGMass();
          geantTrack.fXpos     = tracks.fXposV[t]; // rx of this particle (same as parent)
          geantTrack.fYpos     = tracks.fYposV[t]; // ry of this particle (same as parent)
          geantTrack.fZpos     = tracks.fZposV[t]; // rz of this particle (same as parent)
          geantTrack.fXdir     = secLt[isec].GetDirX();     // dirx of this particle (before transform.)
          geantTrack.fYdir     = secLt[isec].GetDirY();     // diry of this particle before transform.)
          geantTrack.fZdir     = secLt[isec].GetDirZ();     // dirz of this particle before transform.)
          double secEkin       = secLt[isec].GetKinE();
          geantTrack.fP        = std::sqrt(secEkin*(secEkin-2.0*geantTrack.fMass)); // momentum of this secondadry particle
          geantTrack.fE        = secEkin+geantTrack.fMass;                          // total E of this secondary particle
          geantTrack.fTime     = tracks.fTimeV[t]; // global time
          geantTrack.fEdep     = 0.;
          geantTrack.fPstep    = 0.;
          geantTrack.fStep     = 0.;
          geantTrack.fSnext    = 0.;
          geantTrack.fSafety   = tracks.fSafetyV[t];
          geantTrack.fBoundary = tracks.fBoundaryV[t];
          geantTrack.fPending  = false;
          *geantTrack.fPath    = *tracks.fPathV[t];
          *geantTrack.fNextpath= *tracks.fPathV[t];
          geantTrack.fMother   = tracks.fParticleV[t];
          // add GeantTrack
          td->fPropagator->AddTrack(geantTrack);
          tracks.AddTrack(geantTrack);
          // increase the number of secondaries inserted
          ++numSecondaries;
        }
      }
    }
  }
  // set the number of secondaries insterted
  nout = numSecondaries;
}


void PhysicsProcessHandler::AtRestAction(Material_t * /*mat*/, int /*ntracks*/, GeantTrack_v & /*tracks*/,
                                         int & /*nout*/, GeantTaskData * /*td*/) {
  //
  // IMPLEMENTATION NEEDED!
  //
}

} // namespace geantphysics
