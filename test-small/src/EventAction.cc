//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "RunAction.hh"
#include "EventActionMessenger.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include "TabulatedDataManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction():
   runAct((RunAction*)G4RunManager::GetRunManager()->GetUserRunAction()),
   printModulo(1),
   eventMessenger(new EventActionMessenger(this))
#ifdef MAKESTAT
   ,EnergyAbs(0), EnergyGap(0),
   TrackLAbs(0), TrackLGap(0),
   fNSteps(0)
   startTime(0), endTime(0)
#endif

{
#ifdef MAKESTAT
   memset(fEdepGap,0,kNlayers*sizeof(G4double));
   memset(fLengthGap,0,kNlayers*sizeof(G4double));
   memset(fEdepAbs,0,kNlayers*sizeof(G4double));
   memset(fLengthAbs,0,kNlayers*sizeof(G4double));
   memset(fProcStat,0,kNProc*sizeof(unsigned long)); 
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*evt*/)
{  

   /*
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    CLHEP::HepRandom::showEngineStatus();
  }
   */
 
/*
  // initialisation per event
  EnergyAbs = EnergyGap = 0.;
  TrackLAbs = TrackLGap = 0.;
*/
  unsigned long sevent = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
  unsigned long cevent = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEvent()+1; 
  G4float ratio = cevent/(G4float)sevent;
//  std::cout<< cevent << "  " << sevent << std::endl;
  if( (cevent % (sevent/100+1) == 0) || (cevent==sevent) ) {
    G4int cc = ratio*55;
    std::cerr << std::setw(3) << (G4int)(ratio*100) << "% [";
    
    for(G4int x=0; x<cc; x++)
       std::cerr << "=";
    for(G4int x=cc; x<55; x++)
       std::cerr << " ";
    std::cerr << "]\r" << std::flush;
  }

#ifdef MAKESTAT
  // init at the beginning of each new event; for our historams  
  memset(fEdepGap,   0, kNlayers*sizeof(G4double));
  memset(fLengthGap, 0, kNlayers*sizeof(G4double));
  memset(fEdepAbs,   0, kNlayers*sizeof(G4double));
  memset(fLengthAbs, 0, kNlayers*sizeof(G4double));

  fNSteps = 0;

  memset(fProcStat,  0, kNProc*sizeof(unsigned long));

  TabulatedDataManager::killedTracks = 0;

  startTime = clock();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* /*evt*/)
{
/*
  //accumulates statistic
  //
  runAct->fillPerEvent(EnergyAbs, EnergyGap, TrackLAbs, TrackLGap);
  
  //print per event (modulo n)
  //

  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    G4cout << "---> End of event: " << evtNb << G4endl;        

    G4cout
       << "   Absorber: total energy: " << std::setw(7)
                                        << G4BestUnit(EnergyAbs,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(TrackLAbs,"Length")
       << G4endl
       << "        Gap: total energy: " << std::setw(7)
                                        << G4BestUnit(EnergyGap,"Energy")
       << "       total track length: " << std::setw(7)
                                        << G4BestUnit(TrackLGap,"Length")
       << G4endl;
  }
*/
#ifdef MAKESTAT
  endTime = clock();
  G4double time = ((G4double)(endTime-startTime)/CLOCKS_PER_SEC);  
  runAct->fillPerEvent(fEdepGap, fLengthGap, fEdepAbs, fLengthAbs, fNSteps, 
                       fProcStat, time);  
#endif
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::FillPerStep
#ifdef MAKESTAT
                             (G4int isgap, G4int layer, G4double edepo, 
                              G4double steplength,  G4int procIndex){
    if(isgap) {
        fEdepGap  [layer] += edepo;
        fLengthGap[layer] += steplength;
    } else {
        fEdepAbs  [layer] += edepo;
        fLengthAbs[layer] += steplength;
    }

    ++fProcStat[procIndex];
#else
                             (G4int, G4int, G4double,G4double,G4int){
#endif
}




