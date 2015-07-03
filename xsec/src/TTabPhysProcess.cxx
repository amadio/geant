// Toy physics processes for our propagator prototype. Currently including:
// - single scattering as a discrete process
// - energy loss as continuous process
// - generic interaction as discrete process, producing secondaries

#include "TTabPhysProcess.h"

#include "TTabPhysMgr.h"
#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "WorkloadManager.h"
#include "TMath.h"
#include "globals.h"
#include "GeantTrack.h"
#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
#include "management/GeoManager.h"
typedef vecgeom::GeoManager TGeoManager;
#define gGeoManager &vecgeom::GeoManager::Instance()
#else
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoBranchArray.h"
#endif
#include "TGeoMatrix.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"

ClassImp(TTabPhysProcess)

//______________________________________________________________________________
TTabPhysProcess::TTabPhysProcess()
        :PhysicsProcess(),
         fMgr(0),
         fXsecFileName(),
         fFinalSFileName()

{
// I/O ctor
   TObject::SetBit(kDiscrete);
}  

//______________________________________________________________________________
TTabPhysProcess::TTabPhysProcess(const char *name, const char *fxsec, const char *ffstate)
        :PhysicsProcess(name),
         fMgr(0),
         fXsecFileName(fxsec),
         fFinalSFileName(ffstate)
{
// Normal ctor
   TObject::SetBit(kDiscrete);
}  

//______________________________________________________________________________
void TTabPhysProcess::Initialize()
{
// Initialize physics.
   if (fMgr) return;
   if (!gGeoManager) {
      Fatal("Initialize", "Geometry not loaded.");
      return;
   }      
   fMgr = TTabPhysMgr::Instance(gGeoManager, fXsecFileName, fFinalSFileName);
   if (!fMgr) return;   
}



//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
void TTabPhysProcess::ApplyMsc(TGeoMaterial */*mat*/, 
                               Int_t /*ntracks*/, 
                               GeantTrack_v &/*tracks*/, 
                               GeantTaskData */*td*/)
{
   //Apply multiple scattering 
//   fMgr->ApplyMsc(mat, ntracks, tracks, td);
   
} 


//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
void TTabPhysProcess::Eloss(TGeoMaterial *mat,
                            Int_t ntracks, 
                            GeantTrack_v &tracks, 
                            Int_t &nout, 
                            GeantTaskData *td)
{
// Fill energy loss for the tracks according their fStepV

   nout = fMgr->Eloss(mat, ntracks, tracks, td);
}

//______________________________________________________________________________
void TTabPhysProcess::ComputeIntLen(TGeoMaterial *mat, 
                                    Int_t ntracks,
                                    GeantTrack_v &tracks,
                                    Double_t */*lengths*/,
                                    GeantTaskData *td)
{
// Tabulated cross section generic process computation of interaction length.

   fMgr->ProposeStep(mat, ntracks, tracks, td);
}                                      


//______________________________________________________________________________
void TTabPhysProcess::PostStepTypeOfIntrActSampling(TGeoMaterial *mat,
                                                    Int_t ntracks,
                                                    GeantTrack_v &tracks,
                                                    GeantTaskData *td)
{
   // # smapling: target atom and type of the interaction for each primary tracks
   //             all inf. regarding output of sampling is stored in the tracks  
   Int_t imat = -1;
   if (mat)
     imat = mat->GetIndex();
   fMgr->SampleTypeOfInteractions(imat, ntracks, tracks, td);
}


//______________________________________________________________________________
void TTabPhysProcess::PostStepFinalStateSampling(TGeoMaterial *mat,
                                                 Int_t ntracks,
                                                 GeantTrack_v &tracks,
                                                 Int_t &nout,
                                                 GeantTaskData *td)
{
   // # sampling final states for each primary tracks based on target atom and
   //    interaction type sampled in SampleTypeOfInteractionsInt;
   // # upadting primary track properties and inserting secondary tracks;
   // # return: number of inserted secondary tracks  
   Int_t imat = -1;
   if (mat)
     imat = mat->GetIndex();
   nout = fMgr->SampleFinalStates(imat, ntracks, tracks, td);
}


//______________________________________________________________________________
void TTabPhysProcess::AtRest(Int_t /*ntracks*/,
                             GeantTrack_v &/*tracks*/,
                             Int_t &/*nout*/,
                             GeantTaskData */*td*/)
{
// Do at rest actions on particle after generic tabxsec process. 
// Daughter tracks copied in trackout.
}
