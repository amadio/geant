
#include "MyRunAction.hh"

#ifdef G4MULTITHREADED
  #include "G4MTRunManager.hh"
#else
  #include "G4RunManager.hh"
#endif
#include "G4Timer.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "MyRun.hh"
#include "MyPrimaryGeneratorAction.hh"

#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

MyRunAction::MyRunAction() : G4UserRunAction(), fIsPerformance(false), fRun(nullptr), fTimer(nullptr) {}


MyRunAction::~MyRunAction() {}


G4Run* MyRunAction::GenerateRun() {
  // don't Generate our Run in perfomance mode but return with nullptr:
  //   the RunManager::RunInitialization will create a base G4Run object if this method gives back nullptr
  if (!fIsPerformance) {
    fRun = new MyRun();
    return fRun;
  }
  return nullptr;
}


void MyRunAction::BeginOfRunAction(const G4Run* /*aRun*/) {
  if (isMaster) {
      std::vector<G4Region*>* regionVect =  G4RegionStore::GetInstance();
      int numRegions = regionVect->size();
      int sumNumMC = 0;
      G4ProductionCutsTable* pTable = G4ProductionCutsTable::GetProductionCutsTable();
      size_t numMC = pTable->GetTableSize();
      std::vector<int> mcRVect(23,0);
      for (size_t imc=0; imc<numMC; ++imc) {
          const G4MaterialCutsCouple* mc =pTable->GetMaterialCutsCouple(imc);
          if (!mc->IsUsed()) continue;
          G4Material* mat = const_cast<G4Material*>(mc->GetMaterial());
          int k=0;
          for (; k<numRegions && (!(*regionVect)[k]->FindCouple(mat)); ++k){}
              if (k<numRegions) {
                ++mcRVect[k];
              }
      }
      for (int ir=0; ir<numRegions; ++ir) {
          std::cout<< " ir = " << ir << "  region Name = " << (*regionVect)[ir]->GetName() <<" #mc = " << mcRVect[ir] << std::endl;
          sumNumMC += mcRVect[ir];
      }
      std::cout<< " === Total number of MC = " << sumNumMC << " vs " << numMC << " #regions = " << numRegions << std::endl;
      
      
#ifdef G4MULTITHREADED
    G4int nThreads = G4MTRunManager::GetMasterRunManager()->GetNumberOfThreads();
    G4cout << "\n  =======================================================================================\n"
           << "   Run started in MT mode with " << nThreads << "  threads    \n"
           << "  =======================================================================================  \n"
           << G4endl;
#else
    G4cout << "\n  =======================================================================================\n"
           << "   Run started in sequential mode (i.e. with 1 thread)        \n"
           << "  =======================================================================================  \n"
           << G4endl;
#endif
    fTimer = new G4Timer();
    fTimer->Start();
  }
}


void MyRunAction::EndOfRunAction(const G4Run*) {
  if (isMaster) {
    fTimer->Stop();
    // get number of events: even in case of perfomance mode when MyRun-s are not generated in GenerateRun()
#ifdef G4MULTITHREADED
    const G4Run* run = G4MTRunManager::GetMasterRunManager()->GetCurrentRun();
#else
    const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
#endif
    G4int nEvt   = run->GetNumberOfEvent();
    G4cout << "\n\n" << G4endl;
    G4cout << "  =======================================================================================  " << G4endl;
    G4cout << "   Run terminated:                                                                         " << G4endl;
    G4cout << "     Number of events transported    = " << nEvt                                             << G4endl;
    G4cout << "     Time:  "  << *fTimer                                                                   << G4endl;
    G4cout << "  =======================================================================================  "<< G4endl;
    delete fTimer;
    // print primary gun properties (not available at the begining of the run)
    MyPrimaryGeneratorAction::Print();
    if (!fIsPerformance) { // otherwise we do not even create any MyRun objects so fRun is nullptr
      fRun->EndOfRun();
    }
  }
}
