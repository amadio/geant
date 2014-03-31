
#include "TTabPhysMgr.h"

ClassImp(TTabPhysMgr)

TTabPhysMgr* TTabPhysMgr::fgInstance = 0;

//______________________________________________________________________________
TTabPhysMgr* TTabPhysMgr::Instance(TGeoManager* geom, const char* xsecfilename, 
                                   const char* finalsfilename) 
{
// Access to instance of TTabPhysMgr
   if(fgInstance) return fgInstance;
	if(!(geom && xsecfilename && finalsfilename)) {
      Error("Instance", "Create TTabPhysMgr instance providing geometry and xsec files");
      return 0;
   }   
   fgInstance = new TTabPhysMgr(geom, xsecfilename, finalsfilename);			
   return fgInstance;
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
             fGeom(0)
{

 //Load elements from geometry
 TList *matlist = (TList*) geom->GetListOfMaterials();
 fMatXsec = new TMXsec*[matlist->GetSize()];

 TIter next(matlist);
 TGeoMaterial *mat=0;

 //Open xsec_FTFP_BERT.root file (or other phys.lists)
 TFile *f = new TFile(xsecfilename);
 f->Get("PartIndex");

 // Setting the energy grid in our current application (might be different than
 // the one that we used to sample the x-sections from G4)
 TPartIndex::I()->SetEnergyGrid(1e-3,1e3,100);

 //INFO: print number of materials in the current TGeoManager
 printf("#materials:= %d \n",matlist->GetSize());

 // Go through all materials in the geometry and form the associated TMXsec 
 // objects. The necessary elemental TEXsec objects will be loaded (and 
 // interpollated if necessary) as well. (TEXsec::GetElement(...))
 Int_t *z = new Int_t[MAXNELEMENTS];
 Int_t *a = new Int_t[MAXNELEMENTS];
 Float_t *w = new Float_t[MAXNELEMENTS];
 fNmaterials = 0;
 while((mat = (TGeoMaterial*) next())) {
     if(!mat->IsUsed()) continue;
     Int_t nelem = mat->GetNelements();

     // Check if we are on the safe side; should exit otherwise	
     if(nelem>MAXNELEMENTS){
	printf("ERROR: number of elements in %s is %d > TTabPhysMgr::MAXNELEMENTS=%d\n",
		mat->GetName(),nelem,MAXNELEMENTS);
     } 


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
 
 // We can close the xsec_FTFP_BERT.root that we have used so far
 f->Close();

 // After setting up all the necessary TMXsec objects we have the arra of the
 // loaded elemental TEXsec object pointers in: static TEXsec::TEXsec *fElements[NELEM]
 // Since the static TEXsec *fElements[NELEM] is private in TEXsec, I added a getter:
 // IN TEXsec.h:
 // static TEXsec** GetElements(){ return fElements; } 
 // that I will use here to set TTabPhysMgr::fElemXsec ) 
 fElemXsec = TEXsec::GetElements();
 fNelements = TEXsec::NLdElems();
 

 // INFO: print some info for checking	
 printf("number of materials in fMatXsec[]:= %d\n", fNmaterials);
 for(Int_t i=0; i<fNmaterials; ++i)
  printf("Name of %d-th material is fMatXsec[] is %s .\n",i,fMatXsec[i]->GetName());
 // INFO: print some info for checking
 for(Int_t i=0; i<fNelements; ++i)
  printf("%d-th loaded TEXsec is: %s\n",i, TPartIndex::I()->EleName(fElemXsec[i]->Ele()/10000)); 


 // GO FOR FINAL STATES: we know what elemental final states we need based on
 // the elemental TXEsec-s ponters that are already in fElemXsec. So just loop
 // over it end load the corresponding elemental finale states i.e. TEFstate-s 
 
 // Here we could do it in a few different ways. I'm going to be consistent now. 

 // First open the fstate_FTFP_BERT.root file (or other phys.lists)
 f = new TFile(finalsfilename);
 // Than go for the necessary finale states
 for(Int_t i=0; i<fNelements; ++i)
   TEFstate::GetElement(fElemXsec[i]->Ele()/10000,0,f); 

 // At this point, all the necessray TEFstate objects are loaded and their 
 // pointers are stored in: static TEFstate *fElements[NELEM]. So let's get it.
 // Since the static TEFstate *fElements[NELEM] is private in TEFstate, I added 
 // a getter:
 // IN TEFstate.h:
 // static TEFstate** GetElements(){ return fElements; } 
 // that I will use here to set TTabPhysMgr::fElemFstate ) 
 fElemFstate = TEFstate::GetElements();

 // We can close the fstate_FTFP_BERT.root file (or other phys.lists)
 f->Close();

 // INFO: print some info regarding the loaded finale states.
 printf("# TEFstate-s loaded is: %d\n", TEFstate::NLdElems());

 //WE HAVE DONE. EVERYTHING THAT WE NEED IS IN THE MEMORY (IN THEORY)
	
}


   	
//______________________________________________________________________________
TTabPhysMgr::~TTabPhysMgr(){
	delete[] fMatXsec;
}




