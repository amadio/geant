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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "TabulatedDataManager.hh"

#include "TPartIndex.h"

G4bool RunAction::isStatistics = FALSE;  // do you want stat. output ?
G4bool RunAction::isTabPhys = FALSE;     // running with TABPHYS ?


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4int sevent = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
  std::cerr << "<----- Run " << aRun->GetRunID() << " has started and will simulate " <<
            sevent << " evenets! ----->\n";

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
  //initialize cumulative quantities
  //
  sumEAbs = sum2EAbs =sumEGap = sum2EGap = 0.;
  sumLAbs = sum2LAbs =sumLGap = sum2LGap = 0.; 


  if(!RunAction::isStatistics)
    return;

  // init at the beginning of each new run; for our historams  
  memset(fSumEdepGap,   0, kNlayers*sizeof(G4double));
  memset(fSumLengthGap, 0, kNlayers*sizeof(G4double));
  memset(fSumEdepAbs,   0, kNlayers*sizeof(G4double));
  memset(fSumLengthAbs, 0, kNlayers*sizeof(G4double));
  memset(fSum2EdepGap,   0, kNlayers*sizeof(G4double));
  memset(fSum2LengthGap, 0, kNlayers*sizeof(G4double));
  memset(fSum2EdepAbs,   0, kNlayers*sizeof(G4double));
  memset(fSum2LengthAbs, 0, kNlayers*sizeof(G4double));

  sum2LAbs0 = sum2LGap0 = 0.0;          
  sum2EAbs0 = sum2EGap0 = 0.0;          

  fSumNSteps = fSum2NSteps = 0;
         
  memset(fSumProcStat,   0, kNProc*sizeof(G4long));
  memset(fSum2ProcStat,  0, kNProc*sizeof(G4long));

  TabulatedDataManager::killedSecs = 0; //secondaries that are below Elimit
  fSumKilledSecs = fSum2KilledSecs = 0;

  fSumTime = fSum2Time = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G4double EAbs, G4double EGap, G4double LAbs, G4double LGap)
{
  //accumulate statistic
  //
  sumEAbs += EAbs;  sum2EAbs += EAbs*EAbs;
  sumEGap += EGap;  sum2EGap += EGap*EGap;
  
  sumLAbs += LAbs;  sum2LAbs += LAbs*LAbs;
  sumLGap += LGap;  sum2LGap += LGap*LGap;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G4double *EdepGap, G4double *LengthGap, 
                             G4double *EdepAbs, G4double *LengthAbs, G4long numsteps,
                             G4long *ProcStat, G4double time){
     G4double s1,s2,s3,s4;
     s1 = s2 = s3 = s4 = 0.0; 
     for(G4int i=0; i<kNlayers; ++i) {
        fSumEdepGap[i]+=EdepGap[i];  s1+=EdepGap[i];  fSum2EdepGap[i]+=EdepGap[i]*EdepGap[i];
        fSumEdepAbs[i]+=EdepAbs[i];  s2+=EdepAbs[i];  fSum2EdepAbs[i]+=EdepAbs[i]*EdepAbs[i];
        fSumLengthGap[i]+=LengthGap[i]; s3+=LengthGap[i];  fSum2LengthGap[i]+=LengthGap[i]*LengthGap[i];
        fSumLengthAbs[i]+=LengthAbs[i]; s4+=LengthAbs[i]; fSum2LengthAbs[i]+=LengthAbs[i]*LengthAbs[i];
     }

     sum2EGap0+=s1*s1; sum2EAbs0+=s2*s2; sum2LGap0+=s3*s3; sum2LAbs0+=s4*s4;

     fSumNSteps+=numsteps; fSum2NSteps+=numsteps*numsteps;

     for(G4int i=0; i<kNProc; ++i) {
        fSumProcStat[i] +=ProcStat[i];
        fSum2ProcStat[i]+=ProcStat[i]*ProcStat[i];
     }

     fSumKilledSecs += TabulatedDataManager::killedSecs;
     fSum2KilledSecs+= TabulatedDataManager::killedSecs*TabulatedDataManager::killedSecs;

     fSumTime+=time; fSum2Time+=time*time;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  std::cerr<<"\n";
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
/*  
  //compute statistics: mean and rms
  //
  sumEAbs /= NbOfEvents; sum2EAbs /= NbOfEvents;
  G4double rmsEAbs = sum2EAbs - sumEAbs*sumEAbs;
  if (rmsEAbs >0.) rmsEAbs = std::sqrt(rmsEAbs); else rmsEAbs = 0.;
  
  sumEGap /= NbOfEvents; sum2EGap /= NbOfEvents;
  G4double rmsEGap = sum2EGap - sumEGap*sumEGap;
  if (rmsEGap >0.) rmsEGap = std::sqrt(rmsEGap); else rmsEGap = 0.;
  
  sumLAbs /= NbOfEvents; sum2LAbs /= NbOfEvents;
  G4double rmsLAbs = sum2LAbs - sumLAbs*sumLAbs;
  if (rmsLAbs >0.) rmsLAbs = std::sqrt(rmsLAbs); else rmsLAbs = 0.;
  
  sumLGap /= NbOfEvents; sum2LGap /= NbOfEvents;
  G4double rmsLGap = sum2LGap - sumLGap*sumLGap;
  if (rmsLGap >0.) rmsLGap = std::sqrt(rmsLGap); else rmsLGap = 0.;
  
  //print
  //
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << "\n mean Energy in Absorber : " << G4BestUnit(sumEAbs,"Energy")
     << " +- "                          << G4BestUnit(rmsEAbs,"Energy")  
     << "\n mean Energy in Gap      : " << G4BestUnit(sumEGap,"Energy")
     << " +- "                          << G4BestUnit(rmsEGap,"Energy")
     << G4endl;
     
  G4cout 
     << "\n---------  TRACK LENGTH FOR CHARGED PARTICLES ONLY ! ---------"
     << "\n mean trackLength in Absorber : " << G4BestUnit(sumLAbs,"Length")
     << " +- "                               << G4BestUnit(rmsLAbs,"Length")  
     << "\n mean trackLength in Gap      : " << G4BestUnit(sumLGap,"Length")
     << " +- "                               << G4BestUnit(rmsLGap,"Length")
     << "\n------------------------------------------------------------\n"
     << G4endl;
*/

  // for our histograms
  if(RunAction::isStatistics)
  {  

   std::cout<<"\n======================= RUN STAT ==============================\n"
            <<"---------- Average Energy Depo. in ABS and GAP [MeV] ----------\n"
            <<std::endl;
   G4double totalEdepAbs   = 0.0;
   G4double totalTracklAbs = 0.0;
   G4double totalEdepGap   = 0.0;
   G4double totalTracklGap = 0.0;
   for(G4int i=0; i<kNlayers; ++i) {
     G4double a = fSumEdepGap[i]/NbOfEvents;
       G4double rmsgap = fSum2EdepGap[i]/NbOfEvents-a*a;
       if(rmsgap > 0.0) rmsgap = std::sqrt(rmsgap);
       else             rmsgap = 0.0;
     G4double b = fSumEdepAbs[i]/NbOfEvents;
       G4double rmsabs = fSum2EdepAbs[i]/NbOfEvents-b*b;
       if(rmsabs > 0.0) rmsabs = std::sqrt(rmsabs);
       else             rmsabs = 0.0;
     totalEdepAbs += fSumEdepAbs[i];
     totalEdepGap += fSumEdepGap[i]; 
     printf("Layer %d: Eabs=%f sig.=%f  Egap=%f sig.=%f\n", i, b, rmsabs, a, rmsgap );
   }
   totalEdepAbs/=NbOfEvents;
   totalEdepGap/=NbOfEvents; 
   sum2EAbs0/=NbOfEvents;
   sum2EGap0/=NbOfEvents;
   G4double rmsabs0 = sum2EAbs0 - totalEdepAbs*totalEdepAbs; 
   if(rmsabs0 > 0.0) rmsabs0 = std::sqrt(rmsabs0);
   else              rmsabs0 = 0.0;
   G4double rmsgap0 = sum2EGap0 - totalEdepGap*totalEdepGap; 
   if(rmsgap0 > 0.0) rmsgap0 = std::sqrt(rmsgap0);
   else              rmsgap0 = 0.0;
   
   std::cout<<"______________________________________________________________\n"
            <<"TOTAL X: Eabs="<< totalEdepAbs << "  sig.="<<rmsabs0<<  "   Egap=" 
            << totalEdepGap << "  sig.="<< rmsgap0  <<"\n" 
            <<"--------------------------------------------------------------\n"
            << std::endl;

   std::cout<<"---------- Average Track length. in ABS and GAP [cm] ----------\n"
            <<std::endl;
   for(G4int i=0; i<kNlayers; ++i) {
     G4double a = fSumLengthGap[i]/NbOfEvents;
       G4double rmsgap = fSum2LengthGap[i]/NbOfEvents-a*a;
       if(rmsgap > 0.0) rmsgap = std::sqrt(rmsgap);
       else             rmsgap = 0.0;
     G4double b = fSumLengthAbs[i]/NbOfEvents;
       G4double rmsabs = fSum2LengthAbs[i]/NbOfEvents-b*b;
       if(rmsabs > 0.0) rmsabs = std::sqrt(rmsabs);
       else             rmsabs = 0.0;
     totalTracklAbs += fSumLengthAbs[i];
     totalTracklGap += fSumLengthGap[i]; 
     printf("Layer %d: Labs=%f sig.=%f  Lgap=%f sig.=%f\n", i, b, rmsabs, a, rmsgap );
   }
   totalTracklAbs/=NbOfEvents;
   totalTracklGap/=NbOfEvents; 
   sum2LAbs0/=NbOfEvents;
   sum2LGap0/=NbOfEvents;
   rmsabs0 = sum2LAbs0 - totalTracklAbs*totalTracklAbs; 
   if(rmsabs0 > 0.0) rmsabs0 = std::sqrt(rmsabs0);
   else              rmsabs0 = 0.0;
   rmsgap0 = sum2LGap0 - totalTracklGap*totalTracklGap; 
   if(rmsgap0 > 0.0) rmsgap0 = std::sqrt(rmsgap0);
   else              rmsgap0 = 0.0;
   
   std::cout<<"______________________________________________________________\n"
            <<"TOTAL X: Labs="<< totalTracklAbs << "  sig.="<<rmsabs0 <<"   Egap=" 
            << totalTracklGap << "  sig.="<< rmsgap0 << "\n" 
            <<"--------------------------------------------------------------\n"
            << std::endl;


   G4double avrgNSteps = fSumNSteps/((G4double)NbOfEvents);
   G4double rms = fSum2NSteps/((G4double)NbOfEvents) - avrgNSteps*avrgNSteps;
   if(rms>0.0) rms = std::sqrt(rms);
   else        rms = 0.0;
   std::cout<<"------------------- Number of steps --------------------------\n"
            <<" Average # steps = "<< avrgNSteps << " sig.= "<< rms 
            <<"  Total # steps = " << fSumNSteps <<"\n" 
            <<"--------------------------------------------------------------\n"
            << std::endl;

   G4double a = fSumTime/((G4double)NbOfEvents);
   G4double s = fSum2Time/((G4double)NbOfEvents);
   rms = s-a*a;
   if(rms>0.0) rms = std::sqrt(rms);
   else        rms = 0.0;
   std::cout<<"-------------------- Run time in [s] --------------------------\n"
            <<" Average time/event = "<< a << " sig.= "<< rms 
            <<"  Total = " << fSumTime <<"\n" 
            <<"--------------------------------------------------------------\n"
            << std::endl;

/*   
   std::cout<<"------------------ Process Number Stat. -----------------------\n"
            <<" [G5 Proc. Name]     [Total # ]      [Average #]       [Sigma] \n"  
            << std::endl;
   for(G4int i=0; i<kNProc; ++i) {
      G4double a = fSumProcStat[i]/((G4double)NbOfEvents);
      G4double s = fSum2ProcStat[i]/((G4double)NbOfEvents);
      rms = s-a*a;
      if(rms>0.0) rms = std::sqrt(rms);
      else        rms = 0.0;
      if(i!=kNProc-1) // userCut 
        printf(" %-13s: %14d %14.3f %14.3f\n",TPartIndex::I()->ProcName(i),fSumProcStat[i], a,rms);
      else if(!isTabPhys) {
             printf(" %-13s: %14d %14.3f %14.3f\n","UserCut", fSumProcStat[i], a, rms);
           } else {
             a = fSumKilledSecs/((G4double)NbOfEvents);
             s = fSum2KilledSecs/((G4double)NbOfEvents); 
             rms = s-a*a;
             if(rms>0.0) rms = std::sqrt(rms);
             else        rms = 0.0;
             printf(" %-13s: %14d %14.3f %14.3f\n","UserCut", fSumKilledSecs, a, rms);
           }
   }
   std::cout<<"--------------------------------------------------------------\n"
            <<"==============================================================="
            << std::endl;
*/

   std::cout<<"------------------ Process Number Stat. -----------------------\n"
            <<" [G5 Proc. Name]     [Total # ]      [Average #]       [Sigma] \n"  
            << std::endl;
   for(G4int i=0; i<kNProc; ++i) {
      a = fSumProcStat[i]/((G4double)NbOfEvents);
      s = fSum2ProcStat[i]/((G4double)NbOfEvents);
      rms = s-a*a;
      if(rms>0.0) rms = std::sqrt(rms);
      else        rms = 0.0;
      if(i!=kNProc-1) // userCut 
        printf(" %-13s: %14d %14.3f %14.3f\n",TPartIndex::I()->ProcName(i),fSumProcStat[i], a,rms);
      else 
        printf(" %-13s: %14d %14.3f %14.3f\n","UserCut", fSumProcStat[i], a, rms);
   }
   std::cout<<"---------------------------------------------------------------\n"
            << std::endl;
   if(isTabPhys) {
   std::cout<<"---Step and UserCut # corr. (when compared to G4) (estimate)---\n"
            <<"==============================================================="
            << std::endl;
     a = fSumKilledSecs/((G4double)NbOfEvents);
     s = fSum2KilledSecs/((G4double)NbOfEvents); 
     rms = s-a*a;
     if(rms>0.0) rms = std::sqrt(rms);
     else        rms = 0.0;
     std::cout<<" Average #/event = "<< a << " sig.= "<< rms 
              <<"  Total # = " << fSumKilledSecs <<"\n" 
              <<"---------------------------------------------------------------\n"
              << std::endl;
   }


   std::cout<<"==============================================================="
            << std::endl;



 }
 G4int sevent = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
 std::cerr << "<----- Run " << aRun->GetRunID() << " has finished the simulation of " <<
              sevent << " evenets! ---->\n";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
