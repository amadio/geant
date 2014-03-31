#include "TTabPhysMgr.h"

#include "TGeoMaterial.h"
#include "GeantTrack.h"
#include "TBits.h"
#include "TStopwatch.h"
#include "TError.h"
#include "TSystem.h"
#include "TPartIndex.h"
#include "TEXsec.h"
#include "TMXsec.h"
#include "TEFstate.h"

ClassImp(TTabPhysMgr)

TTabPhysMgr* TTabPhysMgr::fgInstance = 0;

//______________________________________________________________________________
TTabPhysMgr* TTabPhysMgr::Instance(TGeoManager* geom, const char* xsecfilename, 
                                   const char* finalsfilename) 
{
// Access to instance of TTabPhysMgr
   if(fgInstance) return fgInstance;
	if(!(geom && xsecfilename && finalsfilename)) {
      ::Error("TTabPhysMgr::Instance", "Create TTabPhysMgr instance providing geometry and xsec files");
      return 0;
   }   
   fgInstance = new TTabPhysMgr(geom, xsecfilename, finalsfilename);			
   return fgInstance;
}

//______________________________________________________________________________
TTabPhysMgr::~TTabPhysMgr()
{
// Destructor
	delete [] fMatXsec;
   delete [] fElemXsec;
   delete [] fElemFstate;
   fgInstance = 0;
}

//______________________________________________________________________________
TTabPhysMgr::TTabPhysMgr():
             TObject(),
             fNelements(0),
             fNmaterials(0),
             fElemXsec(0),
             fElemFstate(0),
             fMatXsec(0),
             fGeom(0)
{
// Dummy ctor.
   fgInstance = this;
}

//______________________________________________________________________________
TTabPhysMgr::TTabPhysMgr(TGeoManager* geom, const char* xsecfilename, 
                         const char* finalsfilename):
             TObject(),
             fNelements(0),
             fNmaterials(0),
             fElemXsec(0),
             fElemFstate(0),
             fMatXsec(0),
             fGeom(geom)
{

 fgInstance = this;
 TStopwatch timer;
 timer.Start();
 //Load elements from geometry
 TList *matlist = (TList*) geom->GetListOfMaterials();

 TIter next(matlist);
 TGeoMaterial *mat=0;

 //Open xsec_FTFP_BERT.root file (or other phys.lists)
 TFile *fxsec = TFile::Open(xsecfilename);
 if (!fxsec) {
    Fatal("TTabPhysMgr", "Cannot open %s", xsecfilename);
 }   
 fxsec->Get("PartIndex");
 // Open the fstate_FTFP_BERT.root file (or other phys.lists)
 TFile *fstate = TFile::Open(finalsfilename);
 if (!fstate) {
    Fatal("TTabPhysMgr", "Cannot open %s", finalsfilename);
 }   

 // Setting the energy grid in our current application (might be different than
 // the one that we used to sample the x-sections from G4)
 TPartIndex::I()->SetEnergyGrid(1e-3,1e3,100); // should be outside

 //INFO: print number of materials in the current TGeoManager
 printf("#materials:= %d \n",matlist->GetSize());

   // First loop on all materials to mark used elements
   TBits elements(NELEM);
   while((mat = (TGeoMaterial*) next())) {
      if(!mat->IsUsed() || mat->GetZ()<1.) continue;
      fNmaterials++;
      Int_t nelem = mat->GetNelements();
      // Check if we are on the safe side; should exit otherwise	
      if(nelem>MAXNELEMENTS){
         Fatal("TTabPhysMgr","Number of elements in %s is %d > TTabPhysMgr::MAXNELEMENTS=%d\n",
               mat->GetName(),nelem,MAXNELEMENTS);
      } 
      for(Int_t iel=0; iel<nelem; ++iel) {
         Double_t ad;
         Double_t zd;
         Double_t wd;
         mat->GetElementProp(ad,zd,wd,iel);
         if (zd<1 || zd>NELEM) {
            Fatal("TTabPhysMgr","In material %s found element with z=%d > NELEM=%d",
                  mat->GetName(), (Int_t)zd, NELEM);
         }
         elements.SetBitNumber(zd);
      }
   }
   fNelements = elements.CountBits();
   fElemXsec = new TEXsec*[NELEM];
   fElemFstate = new TEFstate *[NELEM];
   fMatXsec = new TMXsec*[fNmaterials];
   printf("Reading xsec data and final states for %d elements in %d materials\n",
          fNelements, fNmaterials);
   // Loop elements and load corresponding xsec and final states
   Int_t zel = elements.FirstSetBit();
   Int_t nbits = elements.GetNbits();
   TEXsec *exsec;
   TEFstate *estate;
   // Load elements xsec data in memory
   ProcInfo_t  procInfo1, procInfo2;
   gSystem->GetProcInfo(&procInfo1);
   while (zel<nbits) {
      exsec = TEXsec::GetElement(zel,0,fxsec);
      fElemXsec[zel] = exsec;
      estate = TEFstate::GetElement(zel,0,fstate);
      fElemFstate[zel] = estate;
      printf("   loaded xsec data and states for: %s\n", TPartIndex::I()->EleSymb(zel));
      zel = elements.FirstSetBit(zel+1);
   }
   gSystem->GetProcInfo(&procInfo2);
   Long_t mem = (procInfo2.fMemResident - procInfo1.fMemResident)/1024;
   fxsec->Close();
   fstate->Close();
   // xsec and states now in memory   
   // Go through all materials in the geometry and form the associated TMXsec 
   // objects. 
   Int_t *z = new Int_t[MAXNELEMENTS];
   Int_t *a = new Int_t[MAXNELEMENTS];
   Float_t *w = new Float_t[MAXNELEMENTS];
   fNmaterials = 0;
   next.Reset();
   while((mat = (TGeoMaterial*) next())) {
      if(!mat->IsUsed()) continue;
      Int_t nelem = mat->GetNelements();
      // loop over the elements of the current material in order to obtain the
      // z, a, w, arrays of the elements of this material
      Double_t ad;
      Double_t zd;
      Double_t wd;
      for(Int_t iel=0; iel<nelem; ++iel) {
	      mat->GetElementProp(ad,zd,wd,iel);
         a[iel]=ad;
         z[iel]=zd;
         w[iel]=wd;
      }
      //Construct the TMXsec object that corresponds to the current material
      fMatXsec[fNmaterials++] = new TMXsec(mat->GetName(),mat->GetTitle(),
		       z,a,w,nelem,mat->GetDensity(),kTRUE);
   }// End of while
   delete [] z;
   delete [] a;
   delete [] w;
 

 // After setting up all the necessary TMXsec objects we have the arra of the
 // loaded elemental TEXsec object pointers in: static TEXsec::TEXsec *fElements[NELEM]
 // Since the static TEXsec *fElements[NELEM] is private in TEXsec, I added a getter:
 // IN TEXsec.h:
 // static TEXsec** GetElements(){ return fElements; } 
 // that I will use here to set TTabPhysMgr::fElemXsec ) 
// fElemXsec = TEXsec::GetElements();
  Int_t nelements = TEXsec::NLdElems();
  if (nelements != fNelements) Error("TTabPhysMgr", "Number of elements not matching");
 

 // INFO: print some info for checking	
 printf("number of materials in fMatXsec[]:= %d\n", fNmaterials);
  for(Int_t i=0; i<fNmaterials; ++i)
     printf("   fMatXsec[%d]: %s\n",i,fMatXsec[i]->GetName());
  timer.Stop();   
  printf("Memory taken by xsec and states: %ld [MB] loaded in: %g [sec]\n", mem, timer.CpuTime());
}

//______________________________________________________________________________
void TTabPhysMgr::TransformLF(Int_t indref, GeantTrack_v &tracks, 
                              Int_t nproducts, Int_t indprod, GeantTrack_v &output)
{
// Transform tracks taken from the final state from the local frame to the lab 
// frame (LF). Not clear what parameters to add yet.
// Input: reference track (mother) described as vector container + index of ref track
// Input: number of tracks in the final state, start index and vector container
// Output: roto-boosted tracks in the output vector
}

//______________________________________________________________________________
void TTabPhysMgr::ApplyMsc(Int_t imat, Int_t ntracks, GeantTrack_v &tracks)
{
// Compute MSC angle at the beginning of the step and apply it to the vector
// of tracks.
// Input: material index, number of tracks in the tracks vector to be used
// Output: fXdirV, fYdirV, fZdirV modified in the track container for ntracks
   TGeoMaterial *mat = (TGeoMaterial*)fGeom->GetListOfMaterials()->At(imat);
}

//______________________________________________________________________________
void TTabPhysMgr::Eloss(Int_t imat, Int_t ntracks, GeantTrack_v &tracks)
{
// Apply energy loss for the input material for ntracks in the vector of 
// tracks. Output: modified tracks.fEV array
   TGeoMaterial *mat = (TGeoMaterial*)fGeom->GetListOfMaterials()->At(imat);
}

//______________________________________________________________________________
void TTabPhysMgr::ProposeStep(Int_t imat, Int_t ntracks, GeantTrack_v &tracks)
{
// Sample element in the mixture (still to find where to store it), sample
// physics process based on xsec and store it in tracks.fProcessV. Fill 
// tracks.fPstepV
   TGeoMaterial *mat = (TGeoMaterial*)fGeom->GetListOfMaterials()->At(imat);
}

//______________________________________________________________________________
Int_t TTabPhysMgr::SampleDecay(Int_t ntracks, GeantTrack_v &tracksin, GeantTrack_v &tracksout)
{
// Sample decay for the tracks in the input vector and push the resulting tracks in 
// the output vector. Change status of decayed tracks. Returns number of new tracks.
}

//______________________________________________________________________________
Int_t TTabPhysMgr::SampleInt(Int_t ntracks, GeantTrack_v &tracksin, GeantTrack_v &tracksout)
{
// Sample interaction using the tracksin.fProcessV which is already selected
// Store new tracks in the output vector. Returns number of new tracks.
}
