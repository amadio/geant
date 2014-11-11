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

#include <time.h>

G4bool RunAction::isTabPhys = FALSE;     // running with TABPHYS ?


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction():
   fRunTime(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 

#ifdef MAKESTAT
std::cerr<< "\n ---MAKESTAT = ON  : detailed statistics will be collected.---\n\n";
#else
std::cerr<< "\n ---MAKESTAT = OFF : detailed statistics won't be collected.---\n\n";
#endif  

  G4int sevent = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
  std::cerr << "<----- Run " << aRun->GetRunID() << " has started and will simulate " <<
            sevent << " events ! ----->\n";

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    

#ifdef MAKESTAT
  //initialize cumulative quantities
  //
  sumEAbs = sum2EAbs =sumEGap = sum2EGap = 0.;
  sumLAbs = sum2LAbs =sumLGap = sum2LGap = 0.; 

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
         
  memset(fSumProcStat,   0, kNProc*sizeof(unsigned long));
  memset(fSum2ProcStat,  0, kNProc*sizeof(unsigned long));

  TabulatedDataManager::killedTracks = 0; //secondaries that are below Elimit
  fSumKilledSecs = fSum2KilledSecs = 0;

  fSumTime = fSum2Time = 0.0;
#endif
  
  fRunTime=clock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G4double EAbs, G4double EGap, G4double LAbs, G4double LGap)
{
#ifdef MAKESTAT
  //accumulate statistic
  //
  sumEAbs += EAbs;  sum2EAbs += EAbs*EAbs;
  sumEGap += EGap;  sum2EGap += EGap*EGap;
  
  sumLAbs += LAbs;  sum2LAbs += LAbs*LAbs;
  sumLGap += LGap;  sum2LGap += LGap*LGap;  
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G4double *EdepGap, G4double *LengthGap, 
                             G4double *EdepAbs, G4double *LengthAbs, unsigned long numsteps,
                             unsigned long *ProcStat, G4double time){
#ifdef MAKESTAT
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

     fSumKilledSecs += TabulatedDataManager::killedTracks;
     fSum2KilledSecs+= TabulatedDataManager::killedTracks*TabulatedDataManager::killedTracks;

     fSumTime+=time; fSum2Time+=time*time;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  std::cerr<<"\n";
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
 
  fRunTime=-1.0*(fRunTime-clock());
  std::cout<<"\n-------------------- Run time in [s] --------------------------\n"
            <<"  Total run time= " << ((G4double)(fRunTime)/CLOCKS_PER_SEC) 
            <<"\n--------------------------------------------------------------\n"
            << std::endl;

#ifdef MAKESTAT
   std::cout<<"\n======================= RUN STAT ==============================\n"
            <<"---------- Average Energy Depo. in ABS and GAP [MeV] ----------\n"
            <<std::endl;
   G4double totalEdepAbs   = 0.0;
   G4double totalTracklAbs = 0.0;
   G4double totalEdepGap   = 0.0;
   G4double totalTracklGap = 0.0;
   printf("%-10s %-15s %-15s %-15s %-15s\n","[#Layer]", "[Eabs]","[sig.]","[Egap]","[sig.]");
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
     printf("%4d\t %.6f\t %.6f\t %.6f\t %.6f\n", i, b, rmsabs, a, rmsgap );
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
   printf("%-10s %-15s %-15s %-15s %-15s\n","[#Layer]", "[Labs]","[sig.]","[Lgap]","[sig.]");
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
     printf("%4d\t %.6f\t %.6f\t %.6f\t %.6f\n", i, b, rmsabs, a, rmsgap );
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
            <<"TOTAL X: Labs="<< totalTracklAbs << "  sig.="<<rmsabs0 <<"   Lgap=" 
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
        printf(" %-13s: %14lu %14.3f %14.3f\n",TPartIndex::I()->ProcName(i),fSumProcStat[i], a,rms);
      else 
        printf(" %-13s: %14lu %14.3f %14.3f\n","UserCut", fSumProcStat[i], a, rms);
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
#endif

 G4int sevent = G4RunManager::GetRunManager()->GetCurrentRun()->GetNumberOfEventToBeProcessed();
 std::cerr << "<----- Run " << aRun->GetRunID() << " has finished the simulation of " <<
              sevent << " events! ---->\n";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
