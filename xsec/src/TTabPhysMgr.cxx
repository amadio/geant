#include "TTabPhysMgr.h"

#include "TGeoMaterial.h"
#include "TGeoExtension.h"
#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/navigationstate.h"
#else
#include "TGeoBranchArray.h"
#endif
#include "GeantTrack.h"
#include "globals.h"
#include "GeantPropagator.h"
#include "GeantThreadData.h"
#include "TRandom.h"
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
             fGeom(0),
             fIsRestProcOn(kTRUE) 
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
             fGeom(geom),
             fIsRestProcOn(kTRUE) 
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
      fElemXsec[zel]-> SetIndex(zel); //for quick access to the corresponding fstate 
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
      TMXsec *mxs = new TMXsec(mat->GetName(),mat->GetTitle(),
		       z,a,w,nelem,mat->GetDensity(),kTRUE);
      fMatXsec[fNmaterials++] = mxs;       
      // Connect to TGeoMaterial 
      mat->SetFWExtension(new TGeoRCExtension(mxs));      
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
void TTabPhysMgr::TransformLF(Int_t /*indref*/, GeantTrack_v &/*tracks*/, 
                              Int_t /*nproducts*/, Int_t /*indprod*/, GeantTrack_v &/*output*/)
{
// Transform tracks taken from the final state from the local frame to the lab 
// frame (LF). Not clear what parameters to add yet.
// Input: reference track (mother) described as vector container + index of ref track
// Input: number of tracks in the final state, start index and vector container
// Output: roto-boosted tracks in the output vector
}

//______________________________________________________________________________
void TTabPhysMgr::ApplyMsc(Int_t imat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid)
{
// Compute MSC angle at the beginning of the step and apply it to the vector
// of tracks.
// Input: material index, number of tracks in the tracks vector to be used
// Output: fXdirV, fYdirV, fZdirV modified in the track container for ntracks
   TGeoMaterial *mat = (TGeoMaterial*)fGeom->GetListOfMaterials()->At(imat);
   TMXsec *mxs = ((TMXsec*)((TGeoRCExtension*)mat->GetFWExtension())->GetUserObject());
    
   Float_t msTheta;
   Float_t msPhi;

   Double_t *rndArray = GeantPropagator::Instance()->fThreadData[tid]->fDblArray;
   GeantPropagator::Instance()->fThreadData[tid]->fRndm->RndmArray(ntracks, rndArray);

   for(Int_t i = 0; i < ntracks; ++i){
      msTheta = mxs->MS(tracks.fG5codeV[i], tracks.fEV[i]-tracks.fMassV[i]);     
      msPhi = 2.*TMath::Pi()*rndArray[i];
      RotateTrack(tracks, i, msTheta, msPhi);
   }

}

//______________________________________________________________________________
void TTabPhysMgr::Eloss(Int_t imat, Int_t ntracks, GeantTrack_v &tracks)
{
// Apply energy loss for the input material for ntracks in the vector of 
// tracks. Output: modified tracks.fEV array
   TGeoMaterial *mat = (TGeoMaterial*)fGeom->GetListOfMaterials()->At(imat);
   TMXsec *mxs = ((TMXsec*)((TGeoRCExtension*)mat->GetFWExtension())->GetUserObject());
   mxs->Eloss(ntracks, tracks);
}

//______________________________________________________________________________
void TTabPhysMgr::ProposeStep(Int_t imat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid)
{
// Sample free flight/proposed step for the firts ntracks tracks and store them 
// in tracks.fPstepV  
   TGeoMaterial *mat = (TGeoMaterial*)fGeom->GetListOfMaterials()->At(imat);
   TMXsec *mxs = ((TMXsec*)((TGeoRCExtension*)mat->GetFWExtension())->GetUserObject());
   mxs->ProposeStep(ntracks, tracks, tid);
}

//______________________________________________________________________________
Int_t TTabPhysMgr::SampleDecay(Int_t /*ntracks*/, GeantTrack_v &/*tracksin*/, GeantTrack_v &/*tracksout*/)
{
// Sample decay for the tracks in the input vector and push the resulting tracks in 
// the output vector. Change status of decayed tracks. Returns number of new tracks.
   return 0;
}

//______________________________________________________________________________
Int_t TTabPhysMgr::SampleInt(Int_t imat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid)
{
//-should be called only for particles with reaction (first fNPartReac particle 
// in TPartIndex::fPDG[]); the case ipartindex>=fNPartReac is handled now everywhere
//
// 0. ntracks contains particles with status of Alive 
// 1.Sampling the element of the material for interaction based on the relative 
// total X-secs of the elements; Sampling the type of the interaction (on the 
// sampled element) based on the realtive total X-secs of the interactions ;
// OUT:-indices of the TEXsec* in fElemXsec, that correspond to the sampled 
//      elements, will be in GeantTrack_v::fEindexV array; GeantTrack_v::fEindexV[i] 
//      will be -1 if no reaction for i-th particle
//     -the G5 reaction indices will be in GeantTrack_v::fProcessV array; 
//      GeantTrack_v::fEindexV[i] will be -1 if no reaction for i-th particle     
// 2.Sampling the finale states for the selected interaction and store the secondary
// tracks in tracks; only those traks go into tracks that can pass the energyLimit,
// Rest process final states will be sampled in case of those secondaries that fail to 
// pass the energyLimit (recursion!),
// so the status of tracks can be only kAlive in tracks  
// 4.number of secondary trecks will be returned and original track status will 
// be updated (if they have been killed or still alive) 

   Double_t energyLimit = gPropagator->fEmin;

   TGeoMaterial *mat = (TGeoMaterial*)fGeom->GetListOfMaterials()->At(imat);
   TMXsec *mxs = ((TMXsec*)((TGeoRCExtension*)mat->GetFWExtension())->GetUserObject());

   //1.	
   mxs->SampleInt(ntracks, tracks, tid);

//   for(Int_t i = 0; i<ntracks;++i)
//	printf("[%d]-th Fstate element name:= %s index:= %d\n",i ,fElemXsec[fstateindx[i]]->GetTitle(), fstateindx[i]);



   //2. at Rest story makes this part a bit complicated! (only 'trackable' secondaries can go to tracks)
   // tid-based rng
   Double_t *rndArray = GeantPropagator::Instance()->fThreadData[tid]->fDblArray;
   GeantPropagator::Instance()->fThreadData[tid]->fRndm->RndmArray(2*ntracks, rndArray);

   Int_t nSecPart     = 0;  //number of secondary particles per reaction
   Int_t nTotSecPart  = 0;  //total number of secondary particles in tracks
   const Int_t *pid   = 0;  //GeantV particle codes [nSecPart]
   const Float_t *mom = 0;  //momentum vectors the secondaries [3*nSecPart]
   Float_t  ener      = 0;  //energy at the fstate (Ekin of primary after the interc.)
   Float_t  kerma     = 0;  //released energy
   Float_t  weight    = 0;  //weight of the fstate (just a dummy parameter now)
   Char_t   isSurv    = 0;  //is the primary survived the interaction 	
   
   for(Int_t t = 0; t < ntracks; ++t){
    if(tracks.fProcessV[t] < 0) //if reaction was not selected for this partilce
     continue;	

     isSurv = fElemFstate[tracks.fEindexV[t]]->SampleReac(tracks.fG5codeV[t], 
		tracks.fProcessV[t], tracks.fEV[t]-tracks.fMassV[t], 
		nSecPart, weight, kerma, ener, pid, mom, rndArray[2*t], rndArray[2*t+1]);

     tracks.fEdepV[t] += kerma;    //add the deposited energy (in the interaction) to the energy depositon	

    //if we have secondaries from the current interaction
    if(nSecPart){   
      Double_t oldXdir  = tracks.fXdirV[t];       //old X direction of the primary
      Double_t oldYdir  = tracks.fYdirV[t];       //old Y direction of the primary
      Double_t oldZdir  = tracks.fZdirV[t];       //old Z direction of the primary
      Int_t j = 0;

      if(isSurv && ener >= energyLimit){ //primary particle survived -> it is in the list of secondaries [0]: 
//        tracks.fStatusV[t] = kAlive;   //1. ener=Ekin > energyLimit -> Alive (need to update in tracks) 
 
        //update primary in tracks
        Double_t secPtot2 = mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2];//total P^2 [GeV^2]
        Double_t secPtot  = TMath::Sqrt(secPtot2);     //total P [GeV]
        Double_t secEtot  = ener + tracks.fMassV[t]; //total energy in [GeV]

        tracks.fPV[t]   = secPtot;		//momentum of this particle 
        tracks.fEV[t]   = secEtot;		//total E of this particle 
        tracks.fXdirV[t] = mom[0]/secPtot;	//dirx of this particle (before transform.)
        tracks.fYdirV[t] = mom[1]/secPtot;	//diry of this particle (before transform.)
        tracks.fZdirV[t] = mom[2]/secPtot;	//dirz of this particle (before transform.)

        //Rotate parent track in tracks to original parent track's frame 
        //(a boost will be here as well; before the rotation)		           
        RotateNewTrack(oldXdir, oldYdir, oldZdir, tracks, t);

        j = 1; 
      } else { //2. (primary survived but Ekin < energyLimit) || (primary hasn't survived) -> kill in tracks and invoke at rest
        tracks.fStatusV[t] = kKilled;   //set status of primary in tracks to kKilled;

//      j = 0;
//      tracks.fEdepV[t]  += ener;    //add the after interaction Ekin of the primary to the energy depositon	
        //note: even if the primary in the list of secondaries, its energy (ener)
        //      is lower than the energyLimit so it will be handled properly as 
        //      the other secondaries
      } 

      //loop over the secondaries and put them into tracks:
      // j=0 -> including stopped primary as well if isSurv = kTRUE; 
      // j=1 -> skipp the primary in the list of secondaries (was already updated in tracks above) 
      for(Int_t i = j; i < nSecPart; ++i) {
        Int_t secPDG = TPartIndex::I()->PDG(pid[i]); //Geant V particle code -> particle PGD code
        TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
        Double_t secMass  = secPartPDG->Mass();
        Double_t secPtot2 = mom[3*i]*mom[3*i]+mom[3*i+1]*mom[3*i+1]+mom[3*i+2]*mom[3*i+2];//total P^2 [GeV^2]
        Double_t secPtot  = TMath::Sqrt(secPtot2);//total P [GeV]
        Double_t secEtot  = TMath::Sqrt(secPtot2+ secMass*secMass); //total energy in [GeV]
        Double_t secEkin  = secEtot - secMass; //kinetic energy in [GeV]
        // Ekin of the i-th secondary is higher than the threshold
        if(secEkin >= energyLimit) { //insert secondary into OUT tracks_v and rotate 
//          GeantTrack &gTrack = GeantPropagator::Instance()->GetTempTrack(tid);
          GeantTrack gTrack;
          gTrack.fPath = new TGeoBranchArray();
          gTrack.fNextpath = new TGeoBranchArray();
          

          //set the new track properties
          gTrack.fEvent    = tracks.fEventV[t];
          gTrack.fEvslot   = tracks.fEvslotV[t];
//          gTrack.fParticle = nTotSecPart;          //index of this particle
          gTrack.fPDG      = secPDG;               //PDG code of this particle
          gTrack.fG5code   = pid[i];               //G5 index of this particle
          gTrack.fEindex   = 0;
          gTrack.fCharge   = secPartPDG->Charge()/3.; //charge of this particle
          gTrack.fProcess  = 0;
          gTrack.fIzero    = 0;
          gTrack.fNsteps   = 0;
//          gTrack.fSpecies  = 0;
          gTrack.fStatus   = kNew;                 //status of this particle
          gTrack.fMass     = secMass;              //mass of this particle
          gTrack.fXpos     = tracks.fXposV[t];     //rx of this particle (same as parent)
          gTrack.fYpos     = tracks.fYposV[t];     //ry of this particle (same as parent)
          gTrack.fZpos     = tracks.fZposV[t];     //rz of this particle (same as parent)
          gTrack.fXdir     = mom[3*i]/secPtot;     //dirx of this particle (before transform.)
          gTrack.fYdir     = mom[3*i+1]/secPtot;   //diry of this particle before transform.)
          gTrack.fZdir     = mom[3*i+2]/secPtot;   //dirz of this particle before transform.)
          gTrack.fP        = secPtot;              //momentum of this particle 
          gTrack.fE        = secEtot;              //total E of this particle 
          gTrack.fEdep     = 0.;
          gTrack.fPstep    = 0.;
          gTrack.fStep     = 0.;
          gTrack.fSnext    = 0.;
          gTrack.fSafety   = tracks.fSafetyV[t];
          gTrack.fFrombdr  = tracks.fFrombdrV[t];
          gTrack.fPending  = kFALSE;
          *gTrack.fPath    = *tracks.fPathV[t];
          *gTrack.fNextpath = *tracks.fPathV[t];

          gPropagator->AddTrack(gTrack);
          tracks.AddTrack(gTrack);

          //Rotate new track to parent track's frame (a boost will be here as well; before the rotation)		           
          RotateNewTrack(oldXdir, oldYdir, oldZdir, tracks, nTotSecPart);
           		 
          ++nTotSecPart;
        } else { // {secondary Ekin < energyLimit} -> kill this secondary and call GetRestFinalSates
          tracks.fEdepV[t]   += secEkin;    //add the Ekin of this secondary to the energy depositon	
          GetRestFinSates(pid[i], fElemFstate[tracks.fEindexV[t]], energyLimit, tracks, t, nTotSecPart, tid);  
        } 
      } //end loop over the secondaries
     } else { //nSecPart = 0 i.e. there is no any secondaries -> primary was killed as well
      tracks.fStatusV[t] = kKilled; //set status of primary in tracks to kKilled;
     }

   } //end loop over tracks   

  return nTotSecPart;
}


//______________________________________________________________________________
//will be called recursively; only CaptureAtRest at the moment
void TTabPhysMgr::GetRestFinSates(Int_t partindex, TEFstate *elemfstate, 
        Double_t energyLimit, GeantTrack_v &tracks, Int_t iintrack, 
        Int_t &nTotSecPart, Int_t tid)
{

//   if(!fIsRestProcOn)//Secondaries from rest proc. can be turned off
     return;
   
   Double_t randn = GeantPropagator::Instance()->fThreadData[tid]->fRndm->Rndm();

   Int_t nSecPart     = 0;   //number of secondary particles per reaction
   const Int_t   *pid = 0;  //GeantV particle codes [nSecPart]
   const Float_t *mom = 0;  //momentum vectors the secondaries [3*nSecPart]
   Float_t  ener      = 0;	    //energy at the fstate
   Float_t  kerma     = 0;	    //released energy
   Float_t  weight    = 0;     //weight of the fstate (just a dummy parameter now)
   Char_t   isSurv    = 0;	    //is the primary survived the interaction 	
   
   //sample RestCapture final states for this particle; it is the only at rest proc. at the moment
   isSurv = elemfstate->SampleRestCaptFstate(partindex, nSecPart, weight, kerma, ener, pid, mom, randn);
   //note: parent was already stopped because an at Rest process happend; -> primary is not in the list of secondaries

   //isSurv should always be kFALSE here because primary was stopped -> just a check
   if(isSurv) printf("A stopped particle survived its rest process in TTabPhysMgr::GetRestFinSates!\n"); 

   //if tehere was any energy deposit add it to parent track doposited energy
   tracks.fEdepV[iintrack] += kerma;   

   //loop over the secondaries
   for(Int_t i = 0; i< nSecPart; ++i){
     Int_t secPDG = TPartIndex::I()->PDG(pid[i]); //Geant V particle code -> particle PGD code
     TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
     Double_t secMass  = secPartPDG->Mass();
     Double_t secPtot2 = mom[3*i]*mom[3*i]+mom[3*i+1]*mom[3*i+1]+mom[3*i+2]*mom[3*i+2];//total P^2 [GeV^2]
     Double_t secPtot  = TMath::Sqrt(secPtot2);//total P [GeV]
     Double_t secEtot  = TMath::Sqrt(secPtot2+ secMass*secMass); //total energy in [GeV]
     Double_t secEkin  = secEtot - secMass; //kinetic energy in [GeV]
     // Ekin of the i-th secondary is higher than the threshold
     if(secEkin >= energyLimit) { //insert secondary into tracks_v 
//       GeantTrack &gTrack = GeantPropagator::Instance()->GetTempTrack(tid);
       GeantTrack gTrack;
       gTrack.fPath = new TGeoBranchArray();
       gTrack.fNextpath = new TGeoBranchArray();

       //set the new track properties
       gTrack.fEvent    = tracks.fEventV[iintrack];
       gTrack.fEvslot   = tracks.fEvslotV[iintrack];
//       gTrack.fParticle = nTotSecPart;          //index of this particle
       gTrack.fPDG      = secPDG;               //PDG code of this particle
       gTrack.fG5code   = pid[i];               //G5 index of this particle
       gTrack.fEindex   = 0;
       gTrack.fCharge   = secPartPDG->Charge()/3.; //charge of this particle
       gTrack.fProcess  = 0;
       gTrack.fIzero    = 0;
       gTrack.fNsteps   = 0;
//       gTrack.fSpecies  = 0;
       gTrack.fStatus   = kNew;                 //status of this particle
       gTrack.fMass     = secMass;              //mass of this particle
       gTrack.fXpos     = tracks.fXposV[iintrack];     //rx of this particle (same as parent)
       gTrack.fYpos     = tracks.fYposV[iintrack];     //ry of this particle (same as parent)
       gTrack.fZpos     = tracks.fZposV[iintrack];     //rz of this particle (same as parent)
       gTrack.fXdir     = mom[3*i]/secPtot;     //dirx of this particle (before transform.)
       gTrack.fYdir     = mom[3*i+1]/secPtot;   //diry of this particle before transform.)
       gTrack.fZdir     = mom[3*i+2]/secPtot;   //dirz of this particle before transform.)
       gTrack.fP        = secPtot;              //momentum of this particle 
       gTrack.fE        = secEtot;              //total E of this particle 
       gTrack.fEdep     = 0.;
       gTrack.fPstep    = 0.;
       gTrack.fStep     = 0.;
       gTrack.fSnext    = 0.;
       gTrack.fSafety   = tracks.fSafetyV[iintrack];
       gTrack.fFrombdr  = tracks.fFrombdrV[iintrack];
       gTrack.fPending  = kFALSE;
       *gTrack.fPath    = *tracks.fPathV[iintrack];
       *gTrack.fNextpath = *tracks.fPathV[iintrack];

       gPropagator->AddTrack(gTrack);
       tracks.AddTrack(gTrack);

       //These are secodaries from Rest process -> were sampled at rest -> no rotation and boost
      		 
       ++nTotSecPart; //increase # of secondaries in tracks_v 
       // end if {secondary Ekin >= energyLimit}
     } else { // {secondary Ekin < energyLimit} -> kill this secondary and call GetRestFinalSates
       tracks.fEdepV[iintrack] += secEkin;    //add the Ekin of this secondary to the energy depositon	
       GetRestFinSates(pid[i], elemfstate, energyLimit, tracks, iintrack, nTotSecPart, tid);  //RECURSION
     } 
   }//end loop over the secondaries
}


//_____________________________________________________________________________
// (oldXdir, oldYdir, oldZdir) is the direction vector of parent track in lab. 
// frame; direction vector of the current track, measured from local Z is 
// already updated in GeantTrack_v tracks; here we rotate it to lab. frame
void TTabPhysMgr::RotateNewTrack(Double_t oldXdir, Double_t oldYdir, Double_t oldZdir,
     GeantTrack_v &tracks, Int_t itrack)
{
     const Double_t one  = 1.0f;  
     const Double_t zero = 0.0f; 
     const Double_t amin = 1.0e-10f; 
//     const Double_t one5 = 1.5f; 
//     const Double_t half = 0.5f; 

     Double_t cosTheta0 = oldZdir; 
     Double_t sinTheta0 = TMath::Sqrt(oldXdir*oldXdir + oldYdir*oldYdir);
     Double_t cosPhi0;
     Double_t sinPhi0;

     if(sinTheta0 > amin) {
       cosPhi0 = oldXdir/sinTheta0;
       sinPhi0 = oldYdir/sinTheta0;                     
     } else {
       cosPhi0 = one;
       sinPhi0 = zero;                     
     }
    
     Double_t h0 = tracks.fXdirV[itrack];
     Double_t h1 = sinTheta0*tracks.fZdirV[itrack] + cosTheta0*h0;
     Double_t h2 = tracks.fYdirV[itrack];
 
     tracks.fXdirV[itrack] = h1*cosPhi0 - h2*sinPhi0;
     tracks.fYdirV[itrack] = h1*sinPhi0 + h2*cosPhi0;
     tracks.fZdirV[itrack] = tracks.fZdirV[itrack]*cosTheta0 - h0*sinTheta0;

/* off now; not realy necessary now because (fXdir,fYdir,fZdir) always 
            normalized at entry   
     //renormalization: 
     Double_t delta = one5-half*( tracks.fXdirV[itrack]*tracks.fXdirV[itrack] + 
				  tracks.fYdirV[itrack]*tracks.fYdirV[itrack] +
				  tracks.fZdirV[itrack]*tracks.fZdirV[itrack] );
     tracks.fXdirV[itrack]*=delta;
     tracks.fYdirV[itrack]*=delta;
     tracks.fZdirV[itrack]*=delta;
*/
}


void TTabPhysMgr::RotateTrack(GeantTrack_v &tracks, Int_t itrack, Float_t theta, 
        Float_t phi)
{
     const Double_t one  = 1.0f;  
     const Double_t zero = 0.0f; 
     const Double_t amin = 1.0e-10f; 
     const Double_t one5 = 1.5f; 
     const Double_t half = 0.5f; 

     Double_t cosTheta0 = tracks.fZdirV[itrack]; 
     Double_t sinTheta0 = TMath::Sqrt(tracks.fXdirV[itrack]*tracks.fXdirV[itrack] +
                                      tracks.fYdirV[itrack]*tracks.fYdirV[itrack]
                                     );
     Double_t cosPhi0;
     Double_t sinPhi0;
     Double_t cosTheta = TMath::Cos(theta);
     Double_t sinTheta = TMath::Sin(theta);
//     Double_t cosPhi = TMath::Cos(phi);
//     Double_t sinPhi = TMath::Sin(phi);


     if(sinTheta0 > amin) {
       cosPhi0 = tracks.fXdirV[itrack]/sinTheta0;
       sinPhi0 = tracks.fYdirV[itrack]/sinTheta0;                     
     } else {
       cosPhi0 = one;
       sinPhi0 = zero;                     
     }

     Double_t h0 = sinTheta*TMath::Cos(phi);
     Double_t h1 = sinTheta0*cosTheta + cosTheta0*h0;
     Double_t h2 = sinTheta*TMath::Sin(phi);
 
     tracks.fXdirV[itrack] = h1*cosPhi0 - h2*sinPhi0;
     tracks.fYdirV[itrack] = h1*sinPhi0 + h2*cosPhi0;
     tracks.fZdirV[itrack] = cosTheta*cosTheta0 - h0*sinTheta0;

    //renormalization: ensure normality to avoid accumulated numerical errors
    //                 due to sequential calls of rotation 
     Double_t delta = one5-half*( tracks.fXdirV[itrack]*tracks.fXdirV[itrack] + 
				  tracks.fYdirV[itrack]*tracks.fYdirV[itrack] +
				  tracks.fZdirV[itrack]*tracks.fZdirV[itrack] );
     tracks.fXdirV[itrack]*=delta;
     tracks.fYdirV[itrack]*=delta;
     tracks.fZdirV[itrack]*=delta;

}


