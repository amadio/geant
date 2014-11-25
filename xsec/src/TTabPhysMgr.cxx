#include "TTabPhysMgr.h"

#include "TGeoMaterial.h"
#include "TGeoExtension.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
#include "management/GeoManager.h"
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
#include "TFile.h"
#include "TList.h"
#include "TSystem.h"

#include "TPartIndex.h"
#include "TEXsec.h"
#include "TMXsec.h"
#include "TEFstate.h"
#include "TPDecay.h"

#include <iostream>

#include "base/RNG.h"
#include "Geant/Math.h"

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
   delete fDecay;
   delete fHasNCaptureAtRest; 
   fgInstance = 0;
}

//______________________________________________________________________________
TTabPhysMgr::TTabPhysMgr():
             fNelements(0),
             fNmaterials(0),
             fElemXsec(0),
             fElemFstate(0),
             fMatXsec(0),
             fDecay(0),
             fGeom(0),
             fHasNCaptureAtRest(0)
{
// Dummy ctor.
   fgInstance = this;
}

//______________________________________________________________________________
TTabPhysMgr::TTabPhysMgr(TGeoManager* geom, const char* xsecfilename, 
                         const char* finalsfilename):
   //             TObject(),
             fNelements(0),
             fNmaterials(0),
             fElemXsec(0),
             fElemFstate(0),
             fMatXsec(0),
             fDecay(0),
             fGeom(geom),
             fHasNCaptureAtRest(0) 
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

 // check version of the data files
 if(fgVersion != TPartIndex::I()->Version()) {
    std::cerr
     <<"\n\n*************************************************************\n"
     <<"  ---------------------------ERROR-----------------------------\n" 
     <<"    Your xsec_*.root and fstate_*.root data files at           \n"
     <<"    -> "<< xsecfilename                                     <<"\n"
     <<"    -> "<< finalsfilename                                   <<"\n"                                             
     <<"    Version is       : "<< TPartIndex::I()->VersionMajor()  <<"."
                                << TPartIndex::I()->VersionMinor()  <<"."
                                << TPartIndex::I()->VersionSub()     <<"\n"
     <<"    Required version : " << GetVersion()                     <<"\n"
     <<"    Update your xsec_*.root and fstate_*.root data files !     "
     <<"\n*************************************************************\n\n";                    
     exit(EXIT_FAILURE);  
 }

 // get the decay table from the final state file 
 fDecay = (TPDecay*)fstate->Get("DecayTable");

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
      // init : does the particle have nuclear cpature at rest? array
      if(!fHasNCaptureAtRest) {
        Int_t numParticles = TPartIndex::I()->NPart(); 
        fHasNCaptureAtRest = new Bool_t[numParticles];
        for(Int_t ip=0; ip<numParticles; ++ip)
           fHasNCaptureAtRest[ip] = estate->HasRestCapture(ip);
      }
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
      if (nelem==0) {
         mat->Dump();
         Fatal("TTabPhysMgr","The material (%s) seems to have no elements",mat->GetName());
      }
      //Construct the TMXsec object that corresponds to the current material
      TMXsec *mxs = new TMXsec(mat->GetName(),mat->GetTitle(),
                               z,a,w,nelem,mat->GetDensity(),kTRUE,fDecay);
      fMatXsec[fNmaterials++] = mxs;       
      // Connect to TGeoMaterial
      mat->SetFWExtension(new TGeoRCExtension(new TOMXsec(mxs)));
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
                              Int_t /*nproducts*/, Int_t /*indprod*/, GeantTrack_v &/*output*/){
// Transform tracks taken from the final state from the local frame to the lab 
// frame (LF). Not clear what parameters to add yet.
// Input: reference track (mother) described as vector container + index of ref track
// Input: number of tracks in the final state, start index and vector container
// Output: roto-boosted tracks in the output vector
}

// NOT ACTIVE NOW
//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
void TTabPhysMgr::ApplyMsc(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid){
// Compute MSC angle at the beginning of the step and apply it to the vector
// of tracks.
// Input: material index, number of tracks in the tracks vector to be used
// Output: fXdirV, fYdirV, fZdirV modified in the track container for ntracks
 
#ifndef GEANT_CUDA_DEVICE_BUILD
   TMXsec *mxs = ((TOMXsec*)((TGeoRCExtension*)mat->GetFWExtension())->GetUserObject())->MXsec();
#else
   TMXsec *mxs = 0; // NOTE: we need to get it from somewhere ....
   assert(mxs!=0);
#endif
//   static Int_t icnt=0;
   Double_t msTheta;
   Double_t msPhi;

#ifndef GEANT_CUDA_DEVICE_BUILD
   Double_t *rndArray = GeantPropagator::Instance()->fThreadData[tid]->fDblArray;
   GeantPropagator::Instance()->fThreadData[tid]->fRndm->RndmArray(ntracks, rndArray);
#else
   Double_t *rndArray = 0; // NOTE: we need to get it from somewhere ....
   VECGEOM_NAMESPACE::RNG::Instance().uniform_array(ntracks,rndArray,0.,1.);
#endif

//   Double_t dir[3] = {0.,0.,0.};
   for(Int_t i = 0; i < ntracks; ++i){
      msTheta = mxs->MS(tracks.fG5codeV[i], tracks.fEV[i]-tracks.fMassV[i]);
      msPhi = 2.*Math::Pi()*rndArray[i];
/*
      if (icnt<100 && mat->GetZ()>10) {
         Printf("theta=%g  phi=%g", msTheta*TMath::RadToDeg(), msPhi*TMath::RadToDeg());
         dir[0] = tracks.fXdirV[i];
         dir[1] = tracks.fYdirV[i];
         dir[2] = tracks.fZdirV[i];
      }   
*/
      RotateTrack(tracks, i, msTheta, msPhi);
/*
      if (icnt<100 && mat->GetZ()>10) {
         icnt++;
         Double_t dot = dir[0]*tracks.fXdirV[i] + dir[1]*tracks.fYdirV[i] +dir[2]*tracks.fZdirV[i];
         Double_t angle = TMath::ACos(dot)*TMath::RadToDeg();
         Printf("new angle=%g   delta=%g", angle, TMath::Abs(angle-msTheta*TMath::RadToDeg()));
      }   
*/      
   }
}

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
Int_t TTabPhysMgr::Eloss(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid){
// Apply energy loss for the input material for ntracks in the vector of 
// tracks. Output: modified tracks.fEV array

#ifndef GEANT_CUDA_DEVICE_BUILD
   TMXsec *mxs = ((TOMXsec*)((TGeoRCExtension*)mat->GetFWExtension())->GetUserObject())->MXsec();
#else
   TMXsec *mxs = 0; // NOTE: we need to get it from somewhere ....
   assert(mxs!=0);
#endif
   mxs->Eloss(ntracks, tracks);

   //call atRest sampling for tracks that have been stopped by Eloss and has at-rest
   Int_t nTotSecPart  = 0;  //total number of new tracks
   Double_t energyLimit = gPropagator->fEmin;    
   for(Int_t i = 0; i < ntracks; ++i)
     if( tracks.fProcessV[i] == -2 && HasRestProcess(tracks.fG5codeV[i]) )
       GetRestFinStates(tracks.fG5codeV[i], mxs, energyLimit, tracks, i, 
                        nTotSecPart, tid);  

   return nTotSecPart;
}

//______________________________________________________________________________
void TTabPhysMgr::ProposeStep(TGeoMaterial *mat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid)
{
// Sample free flight/proposed step for the firts ntracks tracks and store them 
// in tracks.fPstepV  

   TMXsec *mxs = ((TOMXsec*)((TGeoRCExtension*)mat->GetFWExtension())->GetUserObject())->MXsec();
   mxs->ProposeStep(ntracks, tracks, tid);
}


// Implemented in a different way
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
// Rest process final states will be sampled in case of those secondaries that 
// stopped. So the status of tracks can be only kAlive in tracks  
// 4.number of secondary tracks will be returned and original track status will 
// be updated (if they have been killed) 

   GeantPropagator *propagator = GeantPropagator::Instance();
   Double_t energyLimit = propagator->fEmin;

   TGeoMaterial *mat = (TGeoMaterial*)fGeom->GetListOfMaterials()->At(imat);
   TMXsec *mxs = ((TOMXsec*)((TGeoRCExtension*)mat->GetFWExtension())->GetUserObject())->MXsec();

   //1. sampling: a. decay or something else
   //             b. if else then what on what target?
   //   output of sampling is stored in the tracks 	
   mxs->SampleInt(ntracks, tracks, tid);


   //2.
 
   // tid-based rng
   Double_t *rndArray = propagator->fThreadData[tid]->fDblArray;
   GeantPropagator::Instance()->fThreadData[tid]->fRndm->RndmArray(2*ntracks, rndArray);

   Int_t nTotSecPart  = 0;  //total number of secondary particles in tracks
   
   for(Int_t t = 0; t < ntracks; ++t) {
    // if no interaction was selected for this track (because it doesn't have any)
    if(tracks.fProcessV[t] < 0)
      continue;

    Int_t nSecPart     = 0;  //number of secondary particles per reaction
    const Int_t *pid   = 0;  //GeantV particle codes [nSecPart]
    const Float_t *mom = 0;  //momentum vectors the secondaries [3*nSecPart]
    Float_t  ener      = 0;  //energy at the fstate (Ekin of primary after the interc.)
    Float_t  kerma     = 0;  //released energy
    Float_t  weight    = 0;  //weight of the fstate (just a dummy parameter now)
    Char_t   isSurv    = 0;  //is the primary survived the interaction 	
    Int_t    ebinindx  = -1; //energy bin index of the selected final state  

    // firts check the results of interaction sampling:
    if(tracks.fProcessV[t] == 3) {
      // decay : in-flight decay was selected 
      // kill the primary tarck
      tracks.fStatusV[t] = kKilled;
      // sample in-flight decay final state 
      SampleDecayInFlight(tracks.fG5codeV[t], mxs, energyLimit, tracks, t, 
                          nTotSecPart, tid);            
      continue;
    }

    // not decay but something else was selected
    Double_t curPrimEkin = tracks.fEV[t]-tracks.fMassV[t];
    isSurv = fElemFstate[tracks.fEindexV[t]]->SampleReac(tracks.fG5codeV[t], 
		tracks.fProcessV[t], curPrimEkin, nSecPart, weight, kerma, ener, 
                pid, mom, ebinindx, rndArray[2*t], rndArray[2*t+1]);

    // it is the case of: pre-step energy sigma is not zero of this interaction
    //                    but post step is zero-> we don't have final state for 
    //                    this interaction at the postStep energy bin-> do nothing
    // This can happen if interaction type is selected based on the pre-step energy 
    // and there is some energy loss along the step (or something like this)
    if(isSurv && ener<0) // let it go further as it is 
      continue;    

    // we should correct the kerma as well but we don't have enough information
    tracks.fEdepV[t] += kerma;    	

    //if we have secondaries from the current interaction
    if(nSecPart){   
      Double_t oldXdir  = tracks.fXdirV[t];       //old X direction of the primary
      Double_t oldYdir  = tracks.fYdirV[t];       //old Y direction of the primary
      Double_t oldZdir  = tracks.fZdirV[t];       //old Z direction of the primary
      Int_t j = 0;

      // setting the final state correction factor (we scale only the 3-momentums)
      //-get mass of the primary   
      Double_t primMass  = tracks.fMassV[t]; // mass [GeV]
      //-compute corFactor = P_current/P_original = Pz_current/Pz_original 
      // (normaly a check would be good but not necessary: if(ebinindx<0 -> ...) 
      Double_t orgPrimEkin = (TPartIndex::I()->EGrid())[ebinindx]; 
      Double_t corFactor   = Math::Sqrt( curPrimEkin*(curPrimEkin+2.0*primMass) /
                                          (orgPrimEkin*(orgPrimEkin+2.0*primMass)) );
      //-if corFactor is set here to 1.0 --> no correction of the final states
      // corFactor = 1.0;
    
      // check if we need to correct the post-interaction Ekin of the primary:
      // if the primary is survived and has non-zero Ekin --> compute its corrected Ekin
      Double_t postEkinOfParimary = ener;
      if(isSurv && (postEkinOfParimary > 0.0)) { //survived and not stopped 
        // get corrected 3-momentum of the post-interaction primary
        Double_t px = mom[0];
        Double_t py = mom[1];
        Double_t pz = mom[2];
        px *= corFactor;
        py *= corFactor;
        pz *= corFactor;
        // compute corrected P^2 in [GeV^2]
        Double_t postPrimP2 = px*px+py*py+pz*pz;
        // recompute post-interaction Ekin of the primary with corrected 3-momentum
        postEkinOfParimary = Math::Sqrt(postPrimP2 + primMass*primMass) - primMass;
      }
 
      if(postEkinOfParimary > energyLimit) { // survived even after the correction and the E-limit.
        // keep alive
//        tracks.fStatusV[t] = kAlive; 
        Double_t px = mom[0];
        Double_t py = mom[1];
        Double_t pz = mom[2];
        px *= corFactor;
        py *= corFactor;
        pz *= corFactor;
        // compute corrected P^2 in [GeV^2]
        Double_t postPrimP2 = px*px+py*py+pz*pz;
        // recompute post-interaction Ekin of the primary with corrected 3-momentum
        postEkinOfParimary = Math::Sqrt(postPrimP2 + primMass*primMass) - primMass;

        //update primary in tracks
        Double_t secPtot  = Math::Sqrt(postPrimP2);               //total P [GeV]
        Double_t secEtot  = postEkinOfParimary + tracks.fMassV[t]; //total energy in [GeV]
        tracks.fPV[t]    = secPtot;		//momentum of this particle 
        tracks.fEV[t]    = secEtot;		//total E of this particle 
        tracks.fXdirV[t] = px/secPtot;	//dirx of this particle (before transform.)
        tracks.fYdirV[t] = py/secPtot;	//diry of this particle (before transform.)
        tracks.fZdirV[t] = pz/secPtot;	//dirz of this particle (before transform.)

        //Rotate parent track in tracks to original parent track's frame 
        RotateNewTrack(oldXdir, oldYdir, oldZdir, tracks, t);          
        // primary track is updated
      } else {
        // Primary particle energy is below tracking limit  
        //-set status of primary in tracks to kKilled;
        tracks.fStatusV[t] = kKilled;
        tracks.fEdepV[t] += postEkinOfParimary;
        // if the primary is stopped i.e. Ekin <= 0 then call at-rest if it has 
        if( isSurv && postEkinOfParimary<=0.0 && HasRestProcess(tracks.fG5codeV[t]) ) 
          GetRestFinStates(tracks.fG5codeV[t], mxs, energyLimit, tracks, t, 
                           nTotSecPart, tid);  
      }
      
      if(isSurv) j=1;
      //loop over the secondaries and put them into tracks if they good to track:
      // j=0 -> including stopped primary as well if isSurv = kTRUE; 
      // j=1 -> skipp the primary in the list of secondaries (was already updated in tracks above) 
      for(Int_t i = j; i < nSecPart; ++i) {
        if(pid[i]>=TPartIndex::I()->NPart()) { // fragment: put its Ekin to energy deposit
          Int_t idummy     = pid[i] - 1000000000;
          Int_t Z          = idummy/10000.;
          Int_t A          = (idummy - Z*10000)/10.;
          Double_t secMass = TPartIndex::I()->GetAprxNuclearMass(Z, A);
          Double_t px = mom[3*i];
          Double_t py = mom[3*i+1];
          Double_t pz = mom[3*i+2];
          px *= corFactor;
          py *= corFactor;
          pz *= corFactor;
          Double_t secPtot2 = px*px+py*py+pz*pz;  //total P^2 [GeV^2]
          tracks.fEdepV[t]+= Math::Sqrt( secPtot2 + secMass*secMass) - secMass;
          continue;
        }
        Int_t secPDG = TPartIndex::I()->PDG(pid[i]); //Geant V particle code -> particle PGD code
        TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
        Double_t secMass  = secPartPDG->Mass();
        Double_t px = mom[3*i];
        Double_t py = mom[3*i+1];
        Double_t pz = mom[3*i+2];
        px *= corFactor;
        py *= corFactor;
        pz *= corFactor;
        Double_t secPtot2 = px*px+py*py+pz*pz;  //total P^2 [GeV^2]
        Double_t secPtot  = Math::Sqrt(secPtot2);//total P [GeV]
        Double_t secEtot  = Math::Sqrt(secPtot2+ secMass*secMass); //total energy in [GeV]
        Double_t secEkin  = secEtot - secMass; //kinetic energy in [GeV]
        // Ekin of the i-th secondary is higher than the threshold
        if(secEkin >= energyLimit) { //insert secondary into OUT tracks_v and rotate 
          GeantTrack &gTrack = propagator->GetTempTrack(tid);
//          GeantTrack gTrack;
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
          gTrack.fXdir     = px/secPtot;     //dirx of this particle (before transform.)
          gTrack.fYdir     = py/secPtot;     //diry of this particle before transform.)
          gTrack.fZdir     = pz/secPtot;     //dirz of this particle before transform.)
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

          // Rotate new track to parent track's frame      
          RotateNewTrack(oldXdir, oldYdir, oldZdir, gTrack);

          propagator->AddTrack(gTrack);
          tracks.AddTrack(gTrack);

          ++nTotSecPart;
        } else { // {secondary Ekin < energyLimit} -> kill this secondary
          tracks.fEdepV[t]   += secEkin;    //add the Ekin of this secondary to the energy depositon	
          // is secEkin <=0 then call at-rest process if the sec. particle has any     
          if( secEkin<=0.0 && HasRestProcess(pid[i]) )
            GetRestFinStates(pid[i], mxs, energyLimit, tracks, t, nTotSecPart, tid); 
        } 
      } //end loop over the secondaries
     } else { //nSecPart = 0 i.e. there is no any secondaries -> primary was killed as well
      tracks.fStatusV[t] = kKilled; //set status of primary in tracks to kKilled;
     }
   } //end loop over tracks   

  return nTotSecPart;
}

// Will be called only if the particle has decay or/and nuclear capture at-rest
//______________________________________________________________________________
//will be called recursively if necessary
GEANT_CUDA_DEVICE_CODE
void TTabPhysMgr::GetRestFinStates(Int_t partindex, TMXsec *mxs,
        Double_t energyLimit, GeantTrack_v &tracks, Int_t iintrack,
        Int_t &nTotSecPart, Int_t tid) {
   // current track should have already been killed before calling  
   const Double_t mecc = 0.00051099906; //e- mass c2 in [GeV] 

   GeantPropagator *propagator =  GeantPropagator::Instance();
   
   Double_t rndArray[3];
#ifndef GEANT_CUDA_DEVICE_BUILD
   propagator->fThreadData[tid]->fRndm->RndmArray(3, rndArray);
#else
   VECGEOM_NAMESPACE::RNG::Instance().uniform_array(3,rndArray,0.,1.);
#endif

   Int_t nSecPart     = 0;  //number of secondary particles per reaction
   const Int_t   *pid = 0;  //GeantV particle codes [nSecPart]
   const Float_t *mom = 0;  //momentum vectors the secondaries [3*nSecPart]
   Float_t  ener      = 0;  //energy at the fstate
   Float_t  kerma     = 0;  //released energy
   Float_t  weight    = 0;  //weight of the fstate (just a dummy parameter now)
   Char_t   isSurv    = 0;  //is the primary survived the interaction 	

   // check if particle is e+ : e+ annihilation at rest if $E_{limit}< m_{e}c^{2}$ 
   if(partindex == TPartIndex::I()->GetSpecGVIndex(1) ) { 
     if(energyLimit < mecc) {
       Double_t randDirZ = 1.0-2.0*rndArray[0];
       Double_t randSinTheta = Math::Sqrt(1.0-randDirZ*randDirZ);
       Double_t randPhi      = 2.0*rndArray[1]*Math::Pi();
       Double_t randDirX = randSinTheta*Math::Cos(randPhi);
       Double_t randDirY = randSinTheta*Math::Sin(randPhi);

    // need to do it one-by-one 
    // 1. gamma   
       GeantTrack &gTrack1 = propagator->GetTempTrack(tid);
       //set the new track properties: 2 gamma with m_{e}*c*c
       gTrack1.fEvent   = tracks.fEventV[iintrack];
       gTrack1.fEvslot  = tracks.fEvslotV[iintrack];
//       gTrack.fParticle = nTotSecPart;          //index of this particle
       gTrack1.fPDG     = 22;  //gamma PDG code
       gTrack1.fG5code  = TPartIndex::I()->GetSpecGVIndex(2); //gamma G5 index
       gTrack1.fEindex  = 0;
       gTrack1.fCharge  = 0.; // charge
       gTrack1.fProcess = 0;
       gTrack1.fIzero   = 0;
       gTrack1.fNsteps  = 0;
//       gTrack.fSpecies  = 0;
       gTrack1.fStatus  = kNew;   //status of this particle
       gTrack1.fMass    = 0.;              //mass of this particle
       gTrack1.fXpos    = tracks.fXposV[iintrack];     //rx of this particle (same as parent)
       gTrack1.fYpos    = tracks.fYposV[iintrack];     //ry of this particle (same as parent)
       gTrack1.fZpos    = tracks.fZposV[iintrack];     //rz of this particle (same as parent)
       gTrack1.fXdir    = randDirX;
       gTrack1.fYdir    = randDirY;
       gTrack1.fZdir    = randDirZ;
       gTrack1.fP       = mecc;            //momentum of this particle 
       gTrack1.fE       = mecc;           //total E of this particle 
       gTrack1.fEdep    = 0.;
       gTrack1.fPstep   = 0.;
       gTrack1.fStep    = 0.;
       gTrack1.fSnext   = 0.;
       gTrack1.fSafety  = tracks.fSafetyV[iintrack];
       gTrack1.fFrombdr = tracks.fFrombdrV[iintrack];
       gTrack1.fPending = kFALSE;
       *gTrack1.fPath   = *tracks.fPathV[iintrack];
       *gTrack1.fNextpath = *tracks.fPathV[iintrack];

       gPropagator->AddTrack(gTrack1);
       tracks.AddTrack(gTrack1);

    // 2. gamma : everything is the same but the direction
       gTrack1.fXdir    = -1.*randDirX;     
       gTrack1.fYdir    = -1.*randDirY;     
       gTrack1.fZdir    = -1.*randDirZ;     

       gPropagator->AddTrack(gTrack1);
       tracks.AddTrack(gTrack1);
  
       nTotSecPart+=2;
       return;
     } else {
       return; 
     }
   }
   // If the stopped particle doesn't have nuclear capture at-rest then decay it
   if(!fHasNCaptureAtRest[partindex]) {
     // Decay at-rest
     // sample final state for decay
     isSurv = fDecay->SampleDecay(partindex, nSecPart, pid, mom);      
   } else {
     // It has nuclear capture at rest so invoke that
     // sample one element of the material
     TEFstate *elemfstate = fElemFstate[mxs->SampleElement(tid)];
     // stample final state for nuclear capture at-rest
     isSurv = elemfstate->SampleRestCaptFstate(partindex, nSecPart, weight, kerma,
                                               ener, pid, mom, rndArray[0]);
   } 

   Double_t randDirX=0;
   Double_t randDirY=0;
   Double_t randDirZ=1;
   Double_t randSinTheta; 
   Double_t randPhi; 

   //note: parent was already stopped because an at Rest process happend; 
   //      -> primary is not in the list of secondaries

   //isSurv should always be kFALSE here because primary was stopped -> just a check
   if(isSurv) 
     printf("A stopped particle survived its rest process in TTabPhysMgr::GetRestFinSates!\n"); 

   // for a random rotation   
   if(nSecPart) {
       randDirZ = 1.0-2.0*rndArray[1];
       randSinTheta = Math::Sqrt((1.0-randDirZ)*(1.0+randDirZ));
       randPhi      = Math::TwoPi()*rndArray[2];
       randDirX = randSinTheta*Math::Cos(randPhi);
       randDirY = randSinTheta*Math::Sin(randPhi);
   }

   //if tehere was any energy deposit add it to parent track doposited energy
   tracks.fEdepV[iintrack] += kerma;   

   //loop over the secondaries
   for(Int_t i = 0; i< nSecPart; ++i){
     if(pid[i]>=TPartIndex::I()->NPart()) { // fragment: put its Ekin to energy deposit
       Int_t idummy            = pid[i] - 1000000000;
       Int_t Z                 = idummy/10000.;
       Int_t A                 = (idummy - Z*10000)/10.;
       Double_t secMass        = TPartIndex::I()->GetAprxNuclearMass(Z, A);
       Double_t px = mom[3*i];
       Double_t py = mom[3*i+1];
       Double_t pz = mom[3*i+2];
       Double_t secPtot2 = px*px+py*py+pz*pz;  //total P^2 [GeV^2]
       tracks.fEdepV[iintrack]+= Math::Sqrt( secPtot2 + secMass*secMass) - secMass;
       continue;
     }

     Int_t secPDG = TPartIndex::I()->PDG(pid[i]); //Geant V particle code -> particle PGD code
     TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
     Double_t secMass  = secPartPDG->Mass();
     Double_t px = mom[3*i];
     Double_t py = mom[3*i+1];
     Double_t pz = mom[3*i+2];
     Double_t secPtot2 = px*px+py*py+pz*pz;  //total P^2 [GeV^2]
     Double_t secPtot  = Math::Sqrt(secPtot2);//total P [GeV]
     Double_t secEtot  = Math::Sqrt(secPtot2+ secMass*secMass); //total energy in [GeV]
     Double_t secEkin  = secEtot - secMass; //kinetic energy in [GeV]
     // Ekin of the i-th secondary is higher than the threshold
     if(secEkin > energyLimit) { //insert secondary into tracks_v 
       GeantTrack &gTrack = GeantPropagator::Instance()->GetTempTrack(tid);
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
       gTrack.fXdir     = px/secPtot;     //dirx of this particle (before transform.)
       gTrack.fYdir     = py/secPtot;     //diry of this particle before transform.)
       gTrack.fZdir     = pz/secPtot;     //dirz of this particle before transform.)
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

       // rotate at-rest secondary by a common random theta and random phi 
       RotateNewTrack(randDirX, randDirY, randDirZ, gTrack);

       gPropagator->AddTrack(gTrack);
       tracks.AddTrack(gTrack);
		 
       ++nTotSecPart; //increase # of secondaries in tracks_v 
     } else {
       // add the Ekin of this secondary to the energy depositon 
       tracks.fEdepV[iintrack] += secEkin;
       // check if it is a stopped particle and call at-rest sampling if necessary
       if( secEkin<=0.0 && HasRestProcess(pid[i]) )
         GetRestFinStates(pid[i], mxs, energyLimit, tracks, iintrack, nTotSecPart, tid);  //RECURSION
     } 
   }//end loop over the secondaries
}

//______________________________________________________________________________
void TTabPhysMgr::SampleDecayInFlight(Int_t partindex, TMXsec *mxs, 
        Double_t energyLimit, GeantTrack_v &tracks, Int_t iintrack, 
        Int_t &nTotSecPart, Int_t tid ) {
    Int_t nSecPart     = 0;  //number of secondary particles per reaction
    const Int_t *pid   = 0;  //GeantV particle codes [nSecPart]
    const Float_t *mom = 0;  //momentum vectors the secondaries [3*nSecPart]
    Char_t   isSurv    = 0;  //is the primary survived the interaction

    isSurv = fDecay->SampleDecay(partindex, nSecPart, pid, mom);
    // isSurv should always be FALSE here because primary was stopped
    if( isSurv ) 
      std::cout<< "\n---       A particle survived its decay!!!       ---\n"
               << "----    In TTabPhysMgr::SampleFinalStateAtRest     ---\n"
               << std::endl;

   if(nSecPart) { 
      // Go for the secondaries
     Double_t beta = tracks.fPV[iintrack]/tracks.fEV[iintrack];
     Double_t bx   = tracks.fXdirV[iintrack]*beta;
     Double_t by   = tracks.fYdirV[iintrack]*beta; 
     Double_t bz   = tracks.fZdirV[iintrack]*beta;
     Double_t b2   = bx*bx + by*by + bz*bz; //it is beta*beta
     Double_t gam  = 1.0 / Math::Sqrt(1.0 - b2);
     Double_t gam2 = b2>0.0 ? (gam - 1.0)/b2 : 0.0;

     for(Int_t isec=0; isec<nSecPart; ++isec) {
       if(pid[isec]>=TPartIndex::I()->NPart()) { // fragment: put its Ekin to energy deposit
         Int_t idummy      = pid[isec] - 1000000000;
         Int_t Z           = idummy/10000.;
         Int_t A           = (idummy - Z*10000)/10.;
         Double_t secMass  = TPartIndex::I()->GetAprxNuclearMass(Z, A);
         Double_t px       = mom[3*isec];
         Double_t py       = mom[3*isec+1];
         Double_t pz       = mom[3*isec+2];
         Double_t secPtot2 = px*px+py*py+pz*pz;  //total P^2 [GeV^2]
         tracks.fEdepV[iintrack] += Math::Sqrt( secPtot2 + secMass*secMass) - secMass;
         continue;
       }

       Int_t secPDG = TPartIndex::I()->PDG(pid[isec]); //GV part.code -> PGD code
       TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
       Double_t secMass  = secPartPDG->Mass(); // mass [GeV]
       Double_t px = mom[3*isec];
       Double_t py = mom[3*isec+1];
       Double_t pz = mom[3*isec+2];
       Double_t secP2 = px*px+py*py+pz*pz;  //total P^2 [GeV^2]
       Double_t secEtot  = Math::Sqrt(secP2 + secMass*secMass); //total E [GeV]
        //Double_t secEkin  = secEtot - secMass; //kinetic energy in [GeV]

       Double_t bp = bx*px + by*py + bz*pz;
       px      = px + gam2*bp*bx +gam*bx*secEtot;
       py      = py + gam2*bp*by +gam*by*secEtot;
       pz      = pz + gam2*bp*bz +gam*bz*secEtot;
       secEtot = gam*(secEtot+bp);

       Double_t secPtot = Math::Sqrt((secEtot-secMass)*(secEtot+secMass));
       Double_t secEkin = secEtot-secMass;
       if(secEkin > energyLimit) { //insert secondary into tracks_v 
         GeantTrack &gTrack = GeantPropagator::Instance()->GetTempTrack(tid);
         //set the new track properties
         gTrack.fEvent    = tracks.fEventV[iintrack];
         gTrack.fEvslot   = tracks.fEvslotV[iintrack];
//         gTrack.fParticle = nTotSecPart;          //index of this particle
         gTrack.fPDG      = secPDG;                 //PDG code of this particle
         gTrack.fG5code   = pid[isec];              //G5 index of this particle
         gTrack.fEindex   = 0;
         gTrack.fCharge   = secPartPDG->Charge()/3.; //charge of this particle
         gTrack.fProcess  = -1;
         gTrack.fIzero    = 0;
         gTrack.fNsteps   = 0;
//         gTrack.fSpecies  = 0;
         gTrack.fStatus   = kNew;                 //status of this particle
         gTrack.fMass     = secMass;              //mass of this particle
         gTrack.fXpos     = tracks.fXposV[iintrack];     //rx of this particle (same as parent)
         gTrack.fYpos     = tracks.fYposV[iintrack];     //ry of this particle (same as parent)
         gTrack.fZpos     = tracks.fZposV[iintrack];     //rz of this particle (same as parent)
         gTrack.fXdir     = px/secPtot;     //dirx of this particle (before transform.)
         gTrack.fYdir     = py/secPtot;     //diry of this particle before transform.)
         gTrack.fZdir     = pz/secPtot;     //dirz of this particle before transform.)
         gTrack.fP        = secPtot;              //momentum of this particle 
         gTrack.fE        = secEtot;              //total E of this particle 
         gTrack.fEdep     = 0.;
         gTrack.fPstep    = 0.;
         gTrack.fStep     = 0.;
         gTrack.fSnext    = 0.;
         gTrack.fSafety   = tracks.fSafetyV[iintrack];
         gTrack.fFrombdr  = tracks.fFrombdrV[iintrack];
         gTrack.fPending  = kFALSE;
        *gTrack.fPath     = *tracks.fPathV[iintrack];
        *gTrack.fNextpath = *tracks.fPathV[iintrack];

         gPropagator->AddTrack(gTrack);
         tracks.AddTrack(gTrack);
		 
         ++nTotSecPart; //increase # of secondaries in tracks_v 
       } else {
         // add the Ekin of this secondary to the energy depositon 
         tracks.fEdepV[iintrack] += secEkin;
         // check if it is a stopped particle and call at-rest sampling if necessary
         if( secEkin<=0.0 && HasRestProcess(pid[isec]) )
           GetRestFinStates(pid[isec], mxs, energyLimit, tracks, iintrack, nTotSecPart, tid); 
       }
     }// end loop over secondaries
   }// end if has secondaries                   
}


//_____________________________________________________________________________
// FOR A SINGLE GeantTrack 
// (oldXdir, oldYdir, oldZdir) is the direction vector of parent track in lab. 
// frame; direction vector of the current track, measured from local Z is 
// already updated in GeantTrack track; here we rotate it to lab. frame
void TTabPhysMgr::RotateNewTrack(Double_t oldXdir, Double_t oldYdir, Double_t oldZdir,
                                 GeantTrack &track){
     const Double_t one  = 1.0;  
     const Double_t zero = 0.0; 
     const Double_t amin = 1.0e-10; 
     const Double_t one5 = 1.5; 
     const Double_t half = 0.5; 

     Double_t cosTheta0 = oldZdir; 
     Double_t sinTheta0 = Math::Sqrt(oldXdir*oldXdir + oldYdir*oldYdir);
     Double_t cosPhi0;
     Double_t sinPhi0;

     if(sinTheta0 > amin) {
       cosPhi0 = oldXdir/sinTheta0;
       sinPhi0 = oldYdir/sinTheta0;                     
     } else {
       cosPhi0 = one;
       sinPhi0 = zero;                     
     }
    
     Double_t h0 = track.fXdir;
     Double_t h1 = sinTheta0*track.fZdir + cosTheta0*h0;
     Double_t h2 = track.fYdir;
 
     track.fXdir = h1*cosPhi0 - h2*sinPhi0;
     track.fYdir = h1*sinPhi0 + h2*cosPhi0;
     track.fZdir = track.fZdir*cosTheta0 - h0*sinTheta0;

     //renormalization: -use 1-th order Taylor aprx. of 1/sqrt(x) around 1.0
     // that should be almost exact since the vector almost normalized!  
     Double_t delta = one5-half*( track.fXdir*track.fXdir + 
				  track.fYdir*track.fYdir +
				  track.fZdir*track.fZdir );
     track.fXdir*=delta;
     track.fYdir*=delta;
     track.fZdir*=delta;

}

//______________________________________________________________________________
// FOR THE itrack-th element of a GeantTrack_v   
// (oldXdir, oldYdir, oldZdir) is the direction vector of parent track in lab. 
// frame; direction vector of the current track, measured from local Z is 
// already updated in GeantTrack track; here we rotate it to lab. frame
GEANT_CUDA_DEVICE_CODE
void TTabPhysMgr::RotateNewTrack(Double_t oldXdir, Double_t oldYdir, Double_t oldZdir,
                                 GeantTrack_v &tracks, Int_t itrack){
     const Double_t one  = 1.0;  
     const Double_t zero = 0.0; 
     const Double_t amin = 1.0e-10; 
     const Double_t one5 = 1.5; 
     const Double_t half = 0.5; 

     Double_t cosTheta0 = oldZdir; 
     Double_t sinTheta0 = Math::Sqrt(oldXdir*oldXdir + oldYdir*oldYdir);
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

     //renormalization: -use 1-th order Taylor aprx. of 1/sqrt(x) around 1.0
     // that should be almost exact since the vector almost normalized!  
     Double_t delta = one5-half*( tracks.fXdirV[itrack]*tracks.fXdirV[itrack] + 
				  tracks.fYdirV[itrack]*tracks.fYdirV[itrack] +
				  tracks.fZdirV[itrack]*tracks.fZdirV[itrack] );
     tracks.fXdirV[itrack]*=delta;
     tracks.fYdirV[itrack]*=delta;
     tracks.fZdirV[itrack]*=delta;

}


//______________________________________________________________________________
// FOR A SINGLE GeantTrack 
// GeantTrack track contains the original direction in lab frame; theta and 
// phi are the scattering angles measured form the particle local Z
void TTabPhysMgr::RotateTrack(GeantTrack &track, Double_t theta, Double_t phi){
     const Double_t one  = 1.0;  
     const Double_t zero = 0.0; 
     const Double_t amin = 1.0e-10; 
     const Double_t one5 = 1.5; 
     const Double_t half = 0.5; 

     Double_t cosTheta0 = track.fZdir; 
     Double_t sinTheta0 = Math::Sqrt(track.fXdir*track.fXdir +track.fYdir*track.fYdir);
     Double_t cosPhi0;
     Double_t sinPhi0;
     Double_t cosTheta = Math::Cos(theta);
     Double_t sinTheta = Math::Sin(theta);


     if(sinTheta0 > amin) {
       cosPhi0 = track.fXdir/sinTheta0;
       sinPhi0 = track.fYdir/sinTheta0;                     
     } else {
       cosPhi0 = one;
       sinPhi0 = zero;                     
     }

     Double_t h0 = sinTheta*Math::Cos(phi);
     Double_t h1 = sinTheta0*cosTheta + cosTheta0*h0;
     Double_t h2 = sinTheta*Math::Sin(phi);
 
     track.fXdir = h1*cosPhi0 - h2*sinPhi0;
     track.fYdir = h1*sinPhi0 + h2*cosPhi0;
     track.fZdir = cosTheta*cosTheta0 - h0*sinTheta0;

    //renormalization: -ensure normality to avoid accumulated numerical errors
    //    due to sequential calls of rotation; avoid 1/sqrt(x) computation by
    //    using the 1-th order Taylor aprx. around 1.0 that should be almost 
    //    exact since the vector almost normalized!
     Double_t delta = one5-half*( track.fXdir*track.fXdir + 
                                  track.fYdir*track.fYdir +
				  track.fZdir*track.fZdir );
     track.fXdir*=delta;
     track.fYdir*=delta;
     track.fZdir*=delta;
}

//______________________________________________________________________________
// FOR THE itrack-th element of a GeantTrack_v   
// GeantTrack_v contains the original direction in lab frame; theta and 
// phi are the scattering angles measured form the particle local Z
GEANT_CUDA_DEVICE_CODE
void TTabPhysMgr::RotateTrack(GeantTrack_v &tracks, Int_t itrack, Double_t theta, 
                              Double_t phi){
     const Double_t one  = 1.0;  
     const Double_t zero = 0.0; 
     const Double_t amin = 1.0e-10; 
     const Double_t one5 = 1.5; 
     const Double_t half = 0.5; 

     Double_t cosTheta0 = tracks.fZdirV[itrack]; 
     Double_t sinTheta0 = Math::Sqrt(tracks.fXdirV[itrack]*tracks.fXdirV[itrack] +
                                     tracks.fYdirV[itrack]*tracks.fYdirV[itrack]
                                     );
     Double_t cosPhi0;
     Double_t sinPhi0;
     Double_t cosTheta = Math::Cos(theta);
     Double_t sinTheta = Math::Sin(theta);


     if(sinTheta0 > amin) {
       cosPhi0 = tracks.fXdirV[itrack]/sinTheta0;
       sinPhi0 = tracks.fYdirV[itrack]/sinTheta0;                     
     } else {
       cosPhi0 = one;
       sinPhi0 = zero;                     
     }

     Double_t h0 = sinTheta*Math::Cos(phi);
     Double_t h1 = sinTheta0*cosTheta + cosTheta0*h0;
     Double_t h2 = sinTheta*Math::Sin(phi);
 
     tracks.fXdirV[itrack] = h1*cosPhi0 - h2*sinPhi0;
     tracks.fYdirV[itrack] = h1*sinPhi0 + h2*cosPhi0;
     tracks.fZdirV[itrack] = cosTheta*cosTheta0 - h0*sinTheta0;

    //renormalization: -ensure normality to avoid accumulated numerical errors
    //    due to sequential calls of rotation; avoid 1/sqrt(x) computation by
    //    using the 1-th order Taylor aprx. around 1.0 that should be almost 
    //    exact since the vector almost normalized!
     Double_t delta = one5-half*( tracks.fXdirV[itrack]*tracks.fXdirV[itrack] + 
				  tracks.fYdirV[itrack]*tracks.fYdirV[itrack] +
				  tracks.fZdirV[itrack]*tracks.fZdirV[itrack] );
     tracks.fXdirV[itrack]*=delta;
     tracks.fYdirV[itrack]*=delta;
     tracks.fZdirV[itrack]*=delta;
}

//______________________________________________________________________________
char* TTabPhysMgr::GetVersion(){
    char *ver = new char[512];
    sprintf(ver,"%d.%d.%d",VersionMajor(),VersionMinor(),VersionSub());
    return ver;
}

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
Bool_t TTabPhysMgr::HasRestProcess(Int_t gvindex){
    return fDecay->HasDecay(gvindex) || fHasNCaptureAtRest[gvindex] ||
           (gvindex == TPartIndex::I()->GetSpecGVIndex(1));
} 

