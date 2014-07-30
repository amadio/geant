// Toy physics processes for our propagator prototype. Currently including:
// - single scattering as a discrete process
// - energy loss as continuous process
// - generic interaction as discrete process, producing secondaries

#include "TTabPhysProcess.h"

#include "TTabPhysMgr.h"
#include "GeantThreadData.h"
#include "GeantVApplication.h"
#include "WorkloadManager.h"
#include "TMath.h"
#include "globals.h"
#include "GeantTrack.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoBranchArray.h"
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
void TTabPhysProcess::ApplyMsc(TGeoMaterial *mat, 
                               Int_t ntracks, 
                               GeantTrack_v &tracks, 
                               Int_t tid)
{
// Temporary switch off MSC !!!
   //Apply multiple scattering 
//   Int_t imat = mat->GetIndex();
//   fMgr->ApplyMsc(imat, ntracks, tracks, tid);
   
} 


//______________________________________________________________________________
void TTabPhysProcess::Eloss(TGeoMaterial *mat,
                            Int_t ntracks, 
                            GeantTrack_v &tracks, 
                            Int_t &nout, 
                            Int_t tid)
{
// Fill energy loss for the tracks according their fStepV
   Int_t imat = mat->GetIndex();
   nout = fMgr->Eloss(imat, ntracks, tracks, tid);
}

//______________________________________________________________________________
void TTabPhysProcess::ComputeIntLen(TGeoMaterial *mat, 
                                      Int_t ntracks, 
                                      GeantTrack_v &tracks,
                                      Double_t */*lengths*/, 
                                      Int_t tid)
{
// Tabulated cross section generic process computation of interaction length.
   Int_t imat = mat->GetIndex();
   fMgr->ProposeStep(imat, ntracks, tracks, tid);
}                                      

//______________________________________________________________________________
void TTabPhysProcess::PostStep(TGeoMaterial *mat,
                                 Int_t ntracks,
                                 GeantTrack_v &tracks, 
                                 Int_t &nout, 
                                 Int_t tid)
{
// Do post-step actions on particle after generic tabxsec process. 
// Surviving tracks copied in trackout.
   Int_t imat = mat->GetIndex();
   nout = fMgr->SampleInt(imat, ntracks, tracks, tid);
}

//______________________________________________________________________________
void TTabPhysProcess::AtRest(Int_t /*ntracks*/,
                                 GeantTrack_v &/*tracks*/, 
                                 Int_t &/*nout*/, 
                                 Int_t /*tid*/)
{
// Do at rest actions on particle after generic tabxsec process. 
// Daughter tracks copied in trackout.
}
