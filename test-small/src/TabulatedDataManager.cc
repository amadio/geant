//
// S.Y. Jun & J. Apostolakis, April 2014
//
#include "G4Track.hh"
#include "TabulatedDataManager.hh"

#include "TSystem.h"
#include "TFile.h"
#include "TError.h"
#include "TBits.h"

#include "TPartIndex.h"
#include "TEXsec.h"
#include "TMXsec.h"
#include "TEFstate.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoExtension.h"

TabulatedDataManager* TabulatedDataManager::theInstance = 0;

TabulatedDataManager* TabulatedDataManager::Instance()
{
  if (theInstance == 0) {
    char* gdmlFileName = getenv("VP_GEOM_GDML");
    char* xsecFileName = getenv("VP_DATA_XSEC");
    char* fstaFileName = getenv("VP_DATA_FSTA");
    
    if(gdmlFileName && xsecFileName && fstaFileName) {    
      TGeoManager* geom = TGeoManager::Import(gdmlFileName);
      theInstance = new TabulatedDataManager(geom,xsecFileName,fstaFileName);
    }
    else {
      ::Error("TabulatedDataManager::Instance",
	      "Missing VP_GEOM_GDML VP_DATA_XSEC VP_DATA_FSTA");
      return 0;
    }
  }
  return theInstance;
}

TabulatedDataManager::~TabulatedDataManager() {
  theInstance = 0;
  delete [] fMatXsec;
  delete [] fElemXsec;
  delete [] fElemFstate;
}

TabulatedDataManager::TabulatedDataManager() :
  fNelements(0),
  fNmaterials(0),
  fElemXsec(0),
  fElemFstate(0),
  fMatXsec(0),
  fGeom(0)
{
}

TabulatedDataManager::TabulatedDataManager(TGeoManager* geom,
					   const char* xsecfilename,
					   const char* finalsfilename) :
  fNelements(0),
  fNmaterials(0),
  fElemXsec(0),
  fElemFstate(0),
  fMatXsec(0),
  fGeom(geom)
{
  //this is clone of TTabPhysMgr::TTabPhysMgr(TGeoManager* geom, 
  //                 const char* xsecfilename, const char* finalsfilename): 

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
  
  Int_t nelements = TEXsec::NLdElems();
  if (nelements != fNelements) Error("TabulatedDataManager",
				     "Number of elements not matching");
  
  // INFO: print some info for checking  
  printf("number of materials in fMatXsec[]:= %d\n", fNmaterials);
  for(Int_t i=0; i<fNmaterials; ++i)
    printf("   fMatXsec[%d]: %s\n",i,fMatXsec[i]->GetName());
}

G4double TabulatedDataManager::GetInteractionLength(G4int imat, 
						    G4int part, 
						    G4double en) 
{
  G4double x = DBL_MAX;
  if (imat >= 0 ||imat< fNmaterials) x = fMatXsec[imat]->Xlength(part,en);
  return x;
}

void TabulatedDataManager::SampleSecondaries(std::vector<G4Track*>* vdp, 
					     G4int imat, 
					     G4Track* atrack)
{
  //select random atom
  G4int indexElement = 0;

  //sample secondaries and store to vdp

}
