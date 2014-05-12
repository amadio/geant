//
// S.Y. Jun & J. Apostolakis, April 2014
//
#include "G4Track.hh"
#include "TabulatedDataManager.hh"

#include "TSystem.h"
#include "TFile.h"
#include "TError.h"
#include "TBits.h"
#include "TMath.h"

#include "TPartIndex.h"
#include "TEXsec.h"
#include "TMXsec.h"
#include "TEFstate.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoExtension.h"
#include "G4ParticleChange.hh"

TabulatedDataManager* TabulatedDataManager::fgInstance = 0;
TGeoManager *TabulatedDataManager::fgGeom= 0 ;         // Pointer to the geometry manager   

TabulatedDataManager* TabulatedDataManager::Instance()
{
  if (fgInstance == 0) {
    // char* gdmlFileName = getenv("VP_GEOM_GDML");
    char* xsecFileName = getenv("VP_DATA_XSEC");
    char* fstaFileName = getenv("VP_DATA_FSTA");
    
    if( !fgGeom ) { 
          // if(gdmlFileName && xsecFileName && fstaFileName) {    
          // fGeom = TGeoManager::Import(gdmlFileName);

          ::Error("TabulatedDataManager::Instance", "Missing pointer to TGeomManager");
          return 0;
    }
    else{
       if(xsecFileName && fstaFileName) {
          fgInstance = new TabulatedDataManager(fgGeom,xsecFileName,fstaFileName);
       }
       else {
          ::Error("TabulatedDataManager::Instance",
                  "Missing VP_DATA_XSEC VP_DATA_FSTA");
          exit(1);
          return 0;
       }
    }
  }
  return fgInstance;
}

TabulatedDataManager::~TabulatedDataManager() {
  fgInstance = 0;
  delete [] fMatXsec;
  delete [] fElemXsec;
  delete [] fElemFstate;
}

TabulatedDataManager::TabulatedDataManager() :
  fNelements(0),
  fNmaterials(0),
  fElemXsec(0),
  fElemFstate(0),
  fMatXsec(0)
  // fgGeom(0)
{
}

TabulatedDataManager::TabulatedDataManager(TGeoManager* geom,
					   const char* xsecfilename,
					   const char* finalsfilename) :
  fNelements(0),
  fNmaterials(0),
  fElemXsec(0),
  fElemFstate(0),
  fMatXsec(0)
  //, fgGeom(geom)
{
  //this is clone of TTabPhysMgr::TTabPhysMgr(TGeoManager* geom, 
  //                 const char* xsecfilename, const char* finalsfilename): 

  std::cout << "TabulatedDataManager - constructor called." << std::endl;

  if( fgGeom != geom ) { 
     Fatal("TabulateDataManager", "Conflicting pointers to TGeoManager"); 
  }

  //Open xsec_FTFP_BERT.root file and fstate_FTFP_BERT.root
  TFile *fxsec = TFile::Open(xsecfilename);
  if (!fxsec) {
    Fatal("TabulatedDataManager", "Cannot open %s", xsecfilename);
  }   
  fxsec->Get("PartIndex");
  
  TFile *fstate = TFile::Open(finalsfilename);
  if (!fstate) {
    Fatal("TabulatedDataManager", "Cannot open %s", finalsfilename);
  }   
  
  //Load elements from geometry
  TList *matlist = (TList*) geom->GetListOfMaterials();
  
  TIter next(matlist);
  TGeoMaterial *mat=0;
  
  // Setting the energy grid in our current application (might be different than
  // the one that we used to sample the x-sections from G4)
  TPartIndex::I()->SetEnergyGrid(1e-3,1e3,100); // should be outside
  
  //INFO: print number of materials in the current TGeoManager
  printf("#materials:= %d \n",matlist->GetSize());
  
  // First loop on all materials to mark used elements
  TBits elements(NELEM);
  while((mat = (TGeoMaterial*) next())) {
    std::cout << "TabulatedDataManager> Checking material " << mat->GetName() << std::endl;
    if(!mat->IsUsed() || mat->GetZ()<1.) continue;
    fNmaterials++;
    Int_t nelem = mat->GetNelements();
    // Check if we are on the safe side; should exit otherwise        
    if(nelem>MAXNELEMENTS){
      Fatal("TabulatedDataManager",
	    "Number of elements in %s is %d > MAXNELEMENTS=%d\n",
	    mat->GetName(),nelem,MAXNELEMENTS);
    } 
    for(Int_t iel=0; iel<nelem; ++iel) {
      Double_t ad;
      Double_t zd;
      Double_t wd;
      mat->GetElementProp(ad,zd,wd,iel);
      if (zd<1 || zd>NELEM) {
	Fatal("TabulatedDataManager",
	      "In material %s found element with z=%d > NELEM=%d",
	      mat->GetName(), (Int_t)zd, NELEM);
      }
      elements.SetBitNumber(zd);
    }
  }
  fNelements = elements.CountBits();
  fElemXsec = new TEXsec*[NELEM];
  fElemFstate = new TEFstate *[NELEM];
  fMatXsec = new TMXsec*[fNmaterials];
  printf("Reading xsec and final states for %d elements in %d materials\n",
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
    fElemXsec[zel]-> SetIndex(zel); //quick access to the corresponding fstate 
    estate = TEFstate::GetElement(zel,0,fstate);
    fElemFstate[zel] = estate;
    printf("   loaded xsec data and states for: %s\n", 
	   TPartIndex::I()->EleSymb(zel));
    zel = elements.FirstSetBit(zel+1);
  }
  
  gSystem->GetProcInfo(&procInfo2);
  // Long_t mem = (procInfo2.fMemResident - procInfo1.fMemResident)/1024;
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
  
  Int_t nelements = TEXsec::NLdElems();
  if (nelements != fNelements) Error("TabulatedDataManager",
				     "Number of elements not matching");
  
  // INFO: print some info for checking  
  printf("number of materials in fMatXsec[]:= %d\n", fNmaterials);
  for(Int_t i=0; i<fNmaterials; ++i)
    printf("   fMatXsec[%d]: %s\n",i,fMatXsec[i]->GetName());
}

G4double TabulatedDataManager::GetInteractionLength(G4int imat, 
						    const G4Track& track)
{
  G4double x = DBL_MAX;

  G4int inpid = track.GetDynamicParticle()->GetPDGcode();
  G4int ipart = TPartIndex::I()->PartIndex(inpid);
  G4double en = track.GetKineticEnergy()/CLHEP::GeV; //E is [GeV] in the tab.data

  if (imat >= 0 ||imat < fNmaterials) x = fMatXsec[imat]->Xlength(ipart,en);

  return x*cm; // length is in [cm] in the tab.data
}

void TabulatedDataManager::SampleSecondaries(std::vector<GXTrack*>* vdp, 
					     G4int imat, 
					     const G4Track* atrack,
					     G4int ireac)
{
  //select random atom
  TGeoMaterial *mat = (TGeoMaterial*)fgGeom->GetListOfMaterials()->At(imat);
  TMXsec *mxs = 
    ((TMXsec*)((TGeoRCExtension*)mat->GetFWExtension())->GetUserObject());

  G4int inpid = atrack->GetDynamicParticle()->GetPDGcode();
  G4int ipart = TPartIndex::I()->PartIndex(inpid);
  G4double kineticEnergy = atrack->GetKineticEnergy();

  G4int indexElem = mxs->SelectElement(ipart,ireac,kineticEnergy);
  
  if(indexElem<0) return;
  
  //sample secondaries and store to vdp
  G4int nSecPart     = 0;  //number of secondary particles per reaction
  // G4int nTotSecPart  = 0;  //total number of secondary particles in tracks
  const G4int *pid   = 0;  //GeantV particle codes [nSecPart]
  const G4float *mom = 0;  //momentum vectors the secondaries [3*nSecPart]
  Float_t  ener      = 0;  //energy at the fstate (Ekin of primary after the interc.)
  G4float  kerma     = 0;  //released energy
  G4float  weight    = 0;  //weight of the fstate (just a dummy parameter now)
  char     isSurv    = 0;  //is the primary survived the interaction   

  isSurv = fElemFstate[indexElem]->SampleReac(ipart, ireac, kineticEnergy, 
				   nSecPart, weight, kerma, ener, pid, mom);
  if(nSecPart>0) {

    for(G4int is = 0 ; is < nSecPart ; ++is) {
      G4int secPDG = TPartIndex::I()->PDG(pid[is]); //Geant V particle code

      TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
      G4double secMass  = secPartPDG->Mass();
      G4double secPtot2 = mom[3*is]*mom[3*is]+mom[3*is+1]*mom[3*is+1]
	                + mom[3*is+2]*mom[3*is+2]; //total P^2 [GeV^2]
      // G4double secPtot  = TMath::Sqrt(secPtot2); //total P [GeV]
      G4double secEtot  = TMath::Sqrt(secPtot2+ secMass*secMass); //total energy in [GeV]
      G4double secEkin  = secEtot - secMass; //kinetic energy in [GeV]

      GXTrack* secTrack = (GXTrack *) malloc(sizeof(GXTrack));
      secTrack->id = secPDG;
      secTrack->px = mom[3*is];
      secTrack->py = mom[3*is+1];
      secTrack->pz = mom[3*is+2];
      secTrack->E  = secEkin;

      vdp->push_back(secTrack);
    }
  }
}

//sampling element for interaction and type of interaction on that element
Int_t TabulatedDataManager::SampleInteraction(  const G4int gvmaterialindex, // root material index
                                                const G4Track &atrack,
                                                Int_t &reactionid) {    
  G4int     partIndex;   // GV particle index 
  G4double  kinEnergy;   // kinetic energy of the particle in GeV
 
  partIndex = TPartIndex::I()->PartIndex( atrack.GetParticleDefinition()->GetPDGEncoding() );
  kinEnergy = atrack.GetKineticEnergy()/CLHEP::GeV;
 
  // sampling element for intercation based on element-wise relative tot-xsecs
  // and sampling the interaction itself based on the relative total xsections 
  // on the selected element
  TEXsec *elemXsec = fMatXsec[gvmaterialindex]->SampleInt( partIndex, kinEnergy,
                                                           reactionid);

  std::cout<< "\n=========SAMPLING INTERACTION====USING=GV-TABPHYS=================\n" 
           << "***  Particle is       = " << TPartIndex::I()->PartName(partIndex) << " \n"  
           << "***  Particle E kin.   = " << kinEnergy << " [GeV]"                << " \n"
           << "***  Selected element  = " << elemXsec->GetName()                  << " \n" 
           << "***  Selected reaction = index: " << reactionid << " which is "
           <<                      TPartIndex::I()->ProcName(reactionid)
           << "\n==================================================================\n"
           << std::endl;

  return elemXsec->Index();
}

//sampling final state, put 2ndaries into the particle change, update primary, 
//do proper transformations if necessary
void TabulatedDataManager::SampleFinalState(const Int_t elementindex, 
                       const Int_t reactionid, const G4Track &atrack, 
                       G4ParticleChange *particlechange){

  Int_t nSecPart     = 0;  //number of secondary particles per reaction
  const Int_t *pid   = 0;  //GeantV particle codes [nSecPart]
  const Float_t *mom = 0;  //momentum vectors the secondaries [3*nSecPart]
  Float_t  energyFst = 0;  //energy at the fstate (Ekin of primary after the interc.)
  Float_t  kerma     = 0;  //released energy
  Float_t  weightFst = 0;  //weight of the fstate (just a dummy parameter now)
  Char_t   isSurv    = 0;  //is the primary survived the interaction
  
  Int_t  partindex    = TPartIndex::I()->PartIndex(atrack.GetParticleDefinition()->GetPDGEncoding());
  Double_t  kinEnergy = atrack.GetDynamicParticle()->GetKineticEnergy()/CLHEP::GeV;
    std::cout << "EKINE " << kinEnergy << std::endl; 
 
  isSurv = fElemFstate[elementindex]->SampleReac(partindex,
                                                 reactionid,
                                                 kinEnergy,
                                                 nSecPart,
                                                 weightFst,
                                                 kerma,
                                                 energyFst,
                                                 pid,
                                                 mom);
  printf("NUM SECONDARIES:= %d  kerma %f  finalEkin %f \n", nSecPart, kerma, energyFst);

  Double_t oldXdir = atrack.GetMomentumDirection().x();
  Double_t oldYdir = atrack.GetMomentumDirection().y();
  Double_t oldZdir = atrack.GetMomentumDirection().z();
  
  // Fill the particle change after transformation
//!! THIS PART IS EXPERIMENTAL !!! JUST FOR TESTING AT THE MOMENT !!
  // 1. update primaty
  if(isSurv && (energyFst > 1.e-5) ) { //survived 
    particlechange->ProposeTrackStatus(fAlive);
    particlechange->ProposeLocalEnergyDeposit(kerma*GeV);
    particlechange->ProposeNonIonizingEnergyDeposit(0.);

    particlechange->ProposeEnergy(energyFst*GeV);
    // rotate direction of primary particle
        //Double_t mass     = track.GetDynamicParticle()->GetMass() 
        Double_t secPtot2 = mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2];//total P^2 [GeV^2]
        Double_t secPtot  = std::sqrt(secPtot2);     //total P [GeV]
        //Double_t secEtot  = ener + tracks.fMassV[t]; //total energy in [GeV]

        //tracks.fPV[t]   = secPtot;		//momentum of this particle 
        //tracks.fEV[t]   = secEtot;		//total E of this particle 
        G4ThreeVector newDir(mom[0]/secPtot,mom[1]/secPtot,mom[2]/secPtot);
        RotateNewTrack(oldXdir, oldYdir, oldZdir, newDir);
        //tracks.fXdirV[t] = mom[0]/secPtot;	//dirx of this particle (before transform.)
        //tracks.fYdirV[t] = mom[1]/secPtot;	//diry of this particle (before transform.)
        //tracks.fZdirV[t] = mom[2]/secPtot;	//dirz of this particle (before transform.)
        particlechange->ProposeMomentumDirection(newDir);

  


       particlechange->ProposeProperTime(atrack.GetProperTime());
       particlechange->ProposePosition(atrack.GetPosition());
       particlechange->ProposeGlobalTime(atrack.GetGlobalTime());
       particlechange->ProposeEnergy(energyFst*GeV);

  }
  else {
     particlechange->ProposeTrackStatus(fStopAndKill);
     particlechange->ProposeLocalEnergyDeposit(kinEnergy*GeV);
     particlechange->ProposeEnergy(0.0);  
  }

    particlechange->SetNumberOfSecondaries(0);

  // 2. fill secondaries
}   
 

//_____________________________________________________________________________
// (oldXdir, oldYdir, oldZdir) is the direction vector of parent track in lab. 
// frame; direction vector of the current track, measured from local Z is 
// already updated in GeantTrack_v tracks; here we rotate it to lab. frame
void TabulatedDataManager::RotateNewTrack(Double_t oldXdir, Double_t oldYdir, 
       Double_t oldZdir, G4ThreeVector &newDir){

     const Double_t one  = 1.0;  
     const Double_t zero = 0.0; 
     const Double_t amin = 1.0e-10; 
     const Double_t one5 = 1.5; 
     const Double_t half = 0.5; 

     Double_t cosTheta0 = oldZdir; 
     Double_t sinTheta0 = std::sqrt(oldXdir*oldXdir + oldYdir*oldYdir);
     Double_t cosPhi0;
     Double_t sinPhi0;

     if(sinTheta0 > amin) {
       cosPhi0 = oldXdir/sinTheta0;
       sinPhi0 = oldYdir/sinTheta0;                     
     } else {
       cosPhi0 = one;
       sinPhi0 = zero;                     
     }
    
     Double_t h0 = newDir.x(); //tracks.fXdirV[itrack];
     Double_t h1 = sinTheta0*newDir.z() /*tracks.fZdirV[itrack]*/ + cosTheta0*h0;
     Double_t h2 = newDir.y(); //tracks.fYdirV[itrack];
 
     newDir.setX( h1*cosPhi0 - h2*sinPhi0 );
     newDir.setY( h1*sinPhi0 + h2*cosPhi0 );
     newDir.setZ( newDir.z()*cosTheta0 - h0*sinTheta0 );

//  not realy necessary now because (fXdir,fYdir,fZdir) always 
//  normalized at entry   
     //renormalization: 
/*   
     Double_t delta = std::sqrt(newDir.x()*newDir.x() + newDir.y()*newDir.y() +
                                  + newDir.z()*newDir.z() );  
     newDir.setX( newDir.x()/delta );
     newDir.setY( newDir.y()/delta );
     newDir.setZ( newDir.z()/delta );
*/
     Double_t delta = one5-half*(newDir.x()*newDir.x() + newDir.y()*newDir.y() +
                                  + newDir.z()*newDir.z() );  
     newDir.setX( newDir.x()*delta );
     newDir.setY( newDir.y()*delta );
     newDir.setZ( newDir.z()*delta );

}

/*
void TTabPhysMgr::RotateTrack(GeantTrack_v &tracks, Int_t itrack, Double_t theta, 
        Double_t phi)
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
*/
 
