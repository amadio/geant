
#include "GVectorPhysicsProcess.h"

#include <iostream>

#include "TPartIndex.h"

#include "globals.h"
#include "GeantTrack.h"
#include "TGeoMaterial.h"


#include "GVComptonProcess.h"

ClassImp(GVectorPhysicsProcess)

GVectorPhysicsProcess::GVectorPhysicsProcess()
      : PhysicsProcess(),
        fVComptonProcess(0), fEnergyLimit(0)
{
   TObject::SetBit(kDiscrete);
}

GVectorPhysicsProcess::GVectorPhysicsProcess(double energyLimit)
      : PhysicsProcess(),
        fVComptonProcess(0), fEnergyLimit(energyLimit)
{
   TObject::SetBit(kDiscrete);
}


GVectorPhysicsProcess::~GVectorPhysicsProcess()
{
  if(fVComptonProcess) delete fVComptonProcess;
}

void GVectorPhysicsProcess::Initialize()
{
// Initialize physics.
   if (fVComptonProcess) return;
   fVComptonProcess = 
     new GVComptonProcess( TPartIndex::I()->ProcIndex("Compton"), fEnergyLimit);
}

void GVectorPhysicsProcess::PostStepFinalStateSampling( TGeoMaterial *mat,
                                                        Int_t ntracks, 
                                                        GeantTrack_v &tracks,
                                                        Int_t &nout, 
                                                        Int_t tid)
{
    if(!fVComptonProcess) {
       std::cerr<<"************************************************************"
                <<" GVectorPhysicsProcess has not been initialized properly:   "
                <<"            fVComptonProcess member is still NULL           " 
                <<"************************************************************"
                <<std::endl;
      exit(EXIT_FAILURE);       
    }
/*
    int numCompton = 0;
    int numGammas  = 0; 
    for(int i=0; i<ntracks; ++i) {
     if(tracks.fPDGV[i] ==11 )   ++numGammas;
     if(tracks.fProcessV[i]==2) ++numCompton;
    }

    std::cout<< "-- numCompton = "<< numCompton 
             << " -- out of "     << ntracks
             << " -- from which "<< numGammas << " is gamma " 
             << " -- for t-ID "<< tid 
             <<std::endl;
*/     

//    fVComptonProcess->DoIWork();
}



