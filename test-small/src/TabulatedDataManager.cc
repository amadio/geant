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
#include "G4ParticleTable.hh"

using CLHEP::GeV;
using CLHEP::cm;

TabulatedDataManager* TabulatedDataManager::fgInstance = 0;
TGeoManager *TabulatedDataManager::fgGeom= 0 ; // Pointer to the geometry manager   

const char* TabulatedDataManager::tStatus[] = {"fAlive", "fStopButAlive", 
                                               "fStopAndKill"};
G4int TabulatedDataManager::fgVerboseLevel = 0;

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
  
  // the common energy grid must be exactly the same as used in tabxsec for the 
  // tabulated data (x-sections and final states) extraction because fstates 
  // cannot be interpolated ! So get the correct grid from the xsec file.
  fxsec->Get("PartIndex"); 
  
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

//_____________________________________________________________________________
void TabulatedDataManager::EnergyLoss(G4int imat, const G4Track &atrack,  
                     const G4Step &astep, G4ParticleChange *particlechange, 
                     G4double energylimit){
  G4int    partIndex;   // GV particle index 
  G4double kinEnergy;   // kinetic energy of the particle in GeV
 
  partIndex = TPartIndex::I()->PartIndex( atrack.GetParticleDefinition()->
                                          GetPDGEncoding() );
  kinEnergy = atrack.GetKineticEnergy()/CLHEP::GeV; // from MeV->GeV

  if (imat < 0 || imat >= fNmaterials) {
    std::cout<< "\n!!*******SAMPLING ENERY LOSS*****USING*GV-TABPHYS***************!!\n" 
             << "***  Particle is       = " << 
                                   TPartIndex::I()->PartName(partIndex) << " \n"  
             << "***  Particle E kin.   = " << kinEnergy << " [GeV]"    << " \n"
             << "***  Material index:   = " << imat 
                    << " OUT OF RANGE: [ 0 ," << fNmaterials  <<  " ]!" << " \n"
             << "!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
             << std::endl;
    Fatal("TabulatedDataManager::EnergyLoss", "Material index is out of range!");
   }

   // Get dE/dx and go for energy loss computation:
   // We need to know one of the elements (index) to check if partice has AtRest 
   G4int    firstElemIndx; 
   G4double dedx  = fMatXsec[imat]->DEdx(partIndex, kinEnergy, firstElemIndx);
   G4double edepo = dedx*(astep.GetStepLength()/CLHEP::cm); // from mm->cm

   if( fgVerboseLevel >=2 )
     std::cout<< "\n=====COMPUTING ALONG-STEP ENERY LOSS==USING=GV-TABPHYS=(dE/dx)====\n" 
            << "***  Particle is       = " << 
                                   TPartIndex::I()->PartName(partIndex) << " \n"  
            << "***  Particle E kin.   = " << kinEnergy << " [GeV]"     << " \n"
            << "***  Selected reaction = index: " << 
                               TPartIndex::I()->ProcIndex("Ionisation") << 
                                               " which is Ionisation. " << " \n"
            << "***  Material name     = " << fMatXsec[imat]->GetName() << " \n"
            << "***  Step lenght       = " << 
                            astep.GetStepLength()/CLHEP::cm << "[cm]  " <<
                                     astep.GetStepLength()/cm << "[cm]" << " \n"
            << "***  dE/dx             = " << dedx << "[GeV/cm] Eloss:= " << 
                                                       edepo << " [GeV]" << "\n"  
            << "\n==================================================================\n"
            << std::endl;

   // If Ekin-EnergyLoss is above the cut then update the current track.
   // Otherwise: particle is stopped, Ekin goes to energy deposit, status is set 
   // to StopButAlive or StopAndKill depending on if the particle does/doesn't 
   // have NuclearCaptureAtRest process.  
   if( (kinEnergy-edepo) > energylimit) {
     particlechange->ProposeEnergy((kinEnergy-edepo)*GeV); //from GeV->MeV
     particlechange->ProposeTrackStatus(fAlive);
     particlechange->ProposeLocalEnergyDeposit(edepo*GeV); //from GeV->MeV
   } else {
     particlechange->ProposeEnergy(0.0);
     particlechange->ProposeLocalEnergyDeposit(kinEnergy*GeV);
     if(fElemFstate[firstElemIndx]->HasRestCapture(partIndex) )
       particlechange->ProposeTrackStatus(fStopButAlive);  
     else
       particlechange->ProposeTrackStatus(fStopAndKill); 
   }    
}

//_____________________________________________________________________________
G4double TabulatedDataManager::GetInteractionLength(G4int imat, 
                                                    const G4Track &atrack){
  G4double x = DBL_MAX;
  G4int    partIndex;   // GV particle index 
  G4double kinEnergy;   // kinetic energy of the particle in GeV
 
  partIndex = TPartIndex::I()->PartIndex( atrack.GetParticleDefinition()->
                                          GetPDGEncoding() );
  kinEnergy = atrack.GetKineticEnergy()/CLHEP::GeV; // from MeV->GeV
 
  if (imat < 0 || imat >= fNmaterials) {
    std::cout<< "\n!!********SAMPLING-INTERACTION-LENGTH***USING*GV-TABPHYS*********!!\n" 
             << "***  Particle is       = " << 
                                   TPartIndex::I()->PartName(partIndex) << " \n"  
             << "***  Particle E kin.   = " << kinEnergy << " [GeV]"    << " \n"
             << "***  Material index:   = " << imat 
                    << " OUT OF RANGE: [ 0 ," << fNmaterials  <<  " ]!" << " \n"
             << "!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
             << std::endl;
    Fatal("TabulatedDataManager::GetInteractionLength", "Material index is out of range!");
   }

  // Get the total mean free path 
  x = fMatXsec[imat]->Xlength(partIndex,kinEnergy);
  return x*cm; // from cm->mm
}

//____________________________________________________________________________
//sampling element for interaction and type of interaction on that element
Int_t TabulatedDataManager::SampleInteraction(  const G4int imat, 
                                                const G4Track &atrack,
                                                Int_t &reactionid) {    
  G4int     partIndex;   // GV particle index 
  G4double  kinEnergy;   // kinetic energy of the particle in GeV
 
  partIndex = TPartIndex::I()->PartIndex( atrack.GetParticleDefinition()->
                                          GetPDGEncoding() );
  kinEnergy = atrack.GetKineticEnergy()/CLHEP::GeV; // from MeV->GeV
 
  // sampling element for intercation based on element-wise relative tot-xsecs
  // and sampling the interaction itself based on the relative total xsections 
  // on the selected element
  TEXsec *elemXsec = fMatXsec[imat]->SampleInt( partIndex, kinEnergy, 
                                                reactionid);

  if( elemXsec ) {                                                           
    if( fgVerboseLevel >= 2)
      std::cout<< "\n=========SAMPLING INTERACTION====USING=GV-TABPHYS=================\n" 
               << "***  Particle is       = " << 
                                   TPartIndex::I()->PartName(partIndex) << " \n"  
               << "***  Particle E kin.   = " << kinEnergy << " [GeV]"  << " \n"
               << "***  Selected element  = " << elemXsec->GetName()    << " \n" 
               << "***  Selected reaction = index: " << reactionid 
                        << " which is " 
                        <<  TPartIndex::I()->ProcName(reactionid)       << " \n"
               << "==================================================================\n"
               << std::endl;

    
    return elemXsec->Index();
  } else {
      std::cout<< "\n!!*******SAMPLING INTERACTION****USING*GV-TABPHYS***************!!\n" 
               << "***  Particle is       = " << 
                                   TPartIndex::I()->PartName(partIndex) << " \n"  
               << "***  Particle E kin.   = " << kinEnergy << " [GeV]"  << " \n"
               << "***  Selected element  = " << elemXsec->GetName()    << " \n" 
               << "***  Selected reaction = index: " << reactionid 
                        << " which is "
               <<       TPartIndex::I()->ProcName(reactionid)           << " \n"
               << "!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
               << std::endl;
      Fatal("TabulatedDataManager::SampleInteraction", "No element selected!");
    //return -1;
  }
}

//_____________________________________________________________________________
//sampling final state, put 2ndaries into the particle change, update primary, 
//do proper transformations if necessary
void TabulatedDataManager::SampleFinalState(const Int_t elementindex, 
                       const Int_t reactionid, const G4Track &atrack, 
                       G4ParticleChange *particlechange, Double_t energylimit){

  Double_t totEdepo  = 0.0;
  Int_t nSecPart     = 0;  //number of secondary particles per reaction
  const Int_t *pid   = 0;  //GeantV particle codes [nSecPart]
  const Float_t *mom = 0;  //momentum vectors the secondaries [3*nSecPart]
  Float_t  energyFst = 0;  //Ekin of primary after the interaction
  Float_t  kerma     = 0;  //released energy
  Float_t  weightFst = 0;  //weight of the fstate (just a dummy parameter now)
  Char_t   isSurv    = 0;  //is the primary survived the interaction
  
  Int_t  partindex    = TPartIndex::I()->PartIndex(
                              atrack.GetParticleDefinition()->GetPDGEncoding());
  Double_t  kinEnergy = atrack.GetDynamicParticle()->GetKineticEnergy()/CLHEP::GeV;
 
  isSurv = fElemFstate[elementindex]->SampleReac(partindex,
                                                 reactionid,
                                                 kinEnergy,
                                                 nSecPart,
                                                 weightFst,
                                                 kerma,
                                                 energyFst,
                                                 pid,
                                                 mom);

  totEdepo = kerma;
  // store original (pre-interaction) direction of the primary particle 
  Double_t oldXdir = atrack.GetMomentumDirection().x();
  Double_t oldYdir = atrack.GetMomentumDirection().y();
  Double_t oldZdir = atrack.GetMomentumDirection().z();
  
  if(isSurv && (energyFst > energylimit) ) { //survived 
    particlechange->ProposeTrackStatus(fAlive);

    particlechange->ProposeEnergy(energyFst*GeV); // from GeV->MeV
    // rotate direction of primary particle
        Double_t secPtot2 = mom[0]*mom[0]+mom[1]*mom[1]+
                            mom[2]*mom[2];  //total P^2 [GeV^2]
        Double_t secPtot  = std::sqrt(secPtot2);     //total P [GeV]

        G4ThreeVector newDir(mom[0]/secPtot,mom[1]/secPtot,mom[2]/secPtot);
        RotateNewTrack(oldXdir, oldYdir, oldZdir, newDir);
        particlechange->ProposeMomentumDirection(newDir);
  }
  else {
     // Particle is stopped, Ekin goes to energy deposit, status is set 
     // to StopButAlive or StopAndKill depending on if the particle does/doesn't 
     // have NuclearCaptureAtRest process.  
     totEdepo += kinEnergy;
     particlechange->ProposeEnergy(0.0);  
     particlechange->ProposeMomentumDirection(0.0,0.0,1.0); // not since <-Ekin =0.
     if(fElemFstate[elementindex]->HasRestCapture(partindex) ) 
       particlechange->ProposeTrackStatus(fStopButAlive);  
     else
       particlechange->ProposeTrackStatus(fStopAndKill); 
  }

  Int_t isec, j = 0;
  if( isSurv ) ++j;  // skipp the first that is the post-interaction primary

  if( fgVerboseLevel >= 2)
    std::cout<<"============SECONDERS=FROM=TAB.=PHYS========================\n" 
           << "***  Primary is a  : = " << 
                                  TPartIndex::I()->PartName(partindex)  << " \n"
           << "***  Kinetic Enery : = " << kinEnergy << " [GeV] "       << " \n"
           << "***  Element is    : = " << 
                                  fElemXsec[elementindex]->GetName()    << " \n"
           << "***  Interaction is: = " << 
                                  TPartIndex::I()->ProcName(reactionid) << " \n"
           << "***  LIST OF SECONDARIES: \n"    
           << "***  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                                                                        << " \n" 
           << "***  Primary status: = "  << 
                           tStatus[particlechange->GetTrackStatus()]    << " \n"
           << "***  Kinetic Enery : = " << 
                           particlechange->GetEnergy()/GeV << " [GeV]"  << " \n"
           << "***  Momentum dir. : = [ "<< 
                        *(particlechange->GetMomentumDirection())<< " ]"<< " \n"
           << "***  ______________________________________________________"
           << std::endl;

  std::vector<G4Track*> secTracks;
  G4int totalNumSec = 0;
  G4double time = atrack.GetGlobalTime();
  G4double weight = particlechange->GetParentWeight();

  for(isec=j; isec < nSecPart; ++isec){
     if(pid[isec] > TPartIndex::I()->NPart())
       continue; // will need to handle fragments and ions later

     Int_t secPDG = TPartIndex::I()->PDG(pid[isec]); //GV part.code -> PGD code
     TParticlePDG *secPartPDG = TDatabasePDG::Instance()->GetParticle(secPDG);
     Double_t secMass  = secPartPDG->Mass(); // mass [GeV]
     Double_t secPtot2 = mom[3*isec]*mom[3*isec]+mom[3*isec+1]*mom[3*isec+1]+
                         mom[3*isec+2]*mom[3*isec+2]; // total P^2 [GeV^2]
     Double_t secPtot  = TMath::Sqrt(secPtot2);       // total P [GeV]
     Double_t secEtot  = TMath::Sqrt(secPtot2+ secMass*secMass); //total E [GeV]
     Double_t secEkin  = secEtot - secMass; //kinetic energy in [GeV]
  
     // Check if Ekin of this secondary is above our energy limit and add it to 
     // the list of secondary tracks. 
     // Otherwise: Ekin of this secondary goes to energy deposit and check if 
     // this secondary has NuclearCaptureAtRest process:
     // If it does has: insert into the list of secondary tracks with a status  
     // of StopButAlive and Ekin = 0.0 (momentum direction is unimoprtant)).
     // If it doesn't have: there is nothing else to do! (drink a coffee!)
     if(secEkin > energylimit) {       
       ++totalNumSec;
      
       G4ParticleDefinition *particleDef = 
                      G4ParticleTable::GetParticleTable()->FindParticle(secPDG);

       G4ThreeVector newDir(mom[3*isec]/secPtot,mom[3*isec+1]/secPtot,
                            mom[3*isec+2]/secPtot);
       RotateNewTrack(oldXdir, oldYdir, oldZdir, newDir);

       G4DynamicParticle* dynamicParticle = 
                        new G4DynamicParticle(particleDef, newDir, secEkin*GeV);        

       G4Track* secTrack = 
                       new G4Track(dynamicParticle, time, atrack.GetPosition());
       secTrack->SetKineticEnergy(secEkin*GeV);
       secTrack->SetTrackStatus(fAlive);
       secTrack->SetTouchableHandle(atrack.GetTouchableHandle());
       secTrack->SetWeight(weight); 

       secTracks.push_back(secTrack);
       if( fgVerboseLevel >= 2)
         std::cout<< "###  "<<totalNumSec<< "-th SECONDARY:"            << " \n" 
                  << "***  Seconder is a : = " << 
                            TPartIndex::I()->PartName(pid[isec])        << " \n"
                  << "***  Momentum dir. : = [ "<< newDir[0] << ", " 
                           << newDir[1] << ", " << newDir[2] << " ]"    << " \n"
                  << "***  Kinetic Energy: = " << secEkin << " [GeV]"   << " \n"
                  << "***  ______________________________________________________"
                  << std::endl;

     } else { 
       totEdepo+=secEkin;
       if(fElemFstate[elementindex]->HasRestCapture(partindex) ) { 
         ++totalNumSec;
      
         G4ParticleDefinition *particleDef = 
                      G4ParticleTable::GetParticleTable()->FindParticle(secPDG);

         G4ThreeVector newDir(0,0,1); // not important since Ekin = 0.;
         G4DynamicParticle* dynamicParticle = new G4DynamicParticle(particleDef,
                                                                    newDir, 0.0);        

         G4Track* secTrack = new G4Track(dynamicParticle, time, 
                                         atrack.GetPosition());
         secTrack->SetKineticEnergy(0.0);
         secTrack->SetTrackStatus(fStopButAlive);
         secTrack->SetTouchableHandle(atrack.GetTouchableHandle());
         secTrack->SetWeight(weight); 
         
         secTracks.push_back(secTrack);
       } // else: do nothing just add secEkin to Edepo that is already done
     }
  }

  // all this verbosity story could be done in a nicer way using the partcile
  // change object after the for-loop below or our G4Track* vector here but I 
  // don't care now (can be polished later) 
  if( fgVerboseLevel >= 2)
    std::cout<<"============================================================"<<std::endl;

  particlechange->SetNumberOfSecondaries(totalNumSec);
  for(G4int i=0; i< totalNumSec; ++i)
     particlechange->AddSecondary(secTracks[i]);

  // Set the overall energy deposit  
  particlechange->ProposeLocalEnergyDeposit(totEdepo*GeV); // from GeV->MeV
}   
 

//_____________________________________________________________________________
// (oldXdir, oldYdir, oldZdir) are the direction vector of parent track in lab. 
// frame; direction vector of the current track, measured from local Z is in the
// G4ThreeVector &newDir; here we rotate it to lab. frame
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
    
     Double_t h0 = newDir.x(); 
     Double_t h1 = sinTheta0*newDir.z() + cosTheta0*h0;
     Double_t h2 = newDir.y(); 
 
     newDir.setX( h1*cosPhi0 - h2*sinPhi0 );
     newDir.setY( h1*sinPhi0 + h2*cosPhi0 );
     newDir.setZ( newDir.z()*cosTheta0 - h0*sinTheta0 );

     Double_t delta = one5-half*(newDir.x()*newDir.x() + newDir.y()*newDir.y() +
                                  + newDir.z()*newDir.z() );  
     newDir.setX( newDir.x()*delta );
     newDir.setY( newDir.y()*delta );
     newDir.setZ( newDir.z()*delta );

}



// This is not for MISI. I don't know what is it for so I will keep it.
//_____________________________________________________________________________
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

