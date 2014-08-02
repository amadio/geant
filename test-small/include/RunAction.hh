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

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

class RunAction : public G4UserRunAction
{
public:
  RunAction();
  virtual ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void   EndOfRunAction(const G4Run*);
    
  void fillPerEvent(G4double, G4double, G4double, G4double); 
  void fillPerEvent(G4double *EdepGap, G4double *LengthGap, G4double *EdepAbs, 
                    G4double *LengthAbs, G4long numsteps, G4long *ProcStat,
                    G4double time);

  static const G4int kNlayers = 10;
  static const G4int kNProc  = 19; // 18+userCuts   

  static G4bool isStatistics;  // do you want stat. output ?
  static G4bool isTabPhys;     // running with TABPHYS ?
 

private:
  G4double sumEAbs, sum2EAbs;
  G4double sumEGap, sum2EGap;
    
  G4double sumLAbs, sum2LAbs;
  G4double sumLGap, sum2LGap;    

  // for the whole run
  G4double  fSumEdepGap[kNlayers];   // Energy deposition per layer
  G4double  fSumLengthGap[kNlayers]; // step length per layer
  G4double  fSumEdepAbs[kNlayers];   // Energy deposition per layer
  G4double  fSumLengthAbs[kNlayers]; // step length per layer
  G4double  fSum2EdepGap[kNlayers];   // Energy deposition per layer
  G4double  fSum2LengthGap[kNlayers]; // step length per layer
  G4double  fSum2EdepAbs[kNlayers];   // Energy deposition per layer
  G4double  fSum2LengthAbs[kNlayers]; // step length per layer
  G4double  sum2LAbs0, sum2LGap0;            
  G4double  sum2EAbs0, sum2EGap0;     

  G4long    fSumNSteps, fSum2NSteps;      

  G4long    fSumProcStat[kNProc];
  G4long    fSum2ProcStat[kNProc];

  G4long    fSumKilledSecs, fSum2KilledSecs;
  
  G4double  fSumTime, fSum2Time; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

