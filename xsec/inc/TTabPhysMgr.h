#ifndef ROOT_TTabPhysMgr
#define ROOT_TTabPhysMgr

#include <TFile.h>
#include <TGeoExtension.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TList.h>

#define MAXNELEMENTS 20 //max number of elements in one material(TMXsec)

// Singleton that handles tabulated physics data loaded for the 
// materials present in the geometry. 

class TEXsec;
class TMXsec;
class TEFstate;
class GeantTrack_v;
class TGeoMaterial;

class TTabPhysMgr : public TObject
{
private:
   Int_t            fNelements;     // Total number of elements in the geometry
   Int_t            fNmaterials;    // Total number of materials in the geometry	
   TEXsec         **fElemXsec;      // Array of x-section pointers per element
   TEFstate       **fElemFstate;    // Array of final state pointers per element
   TMXsec 	  **fMatXsec;	    // Array of x-section pointers per material	
   TGeoManager     *fGeom;	    // Pointer to the geometry manager   
   Bool_t          fIsRestProcOn;   // Use at rest process   
	
   static TTabPhysMgr *fgInstance;	    // Singleton instance

public:   
   TTabPhysMgr();
   TTabPhysMgr(TGeoManager* geom, const char* xsecfilename, 
	       const char* finalsfilename);
   ~TTabPhysMgr();
   static TTabPhysMgr* Instance(TGeoManager* geom=0, const char* xsecfilename=0, 
                                   const char* finalsfilename=0);
   // Rotation+boost utility
   void TransformLF(Int_t indref, GeantTrack_v &tracks, Int_t nproducts, Int_t 
            indprod, GeantTrack_v &output);//not. imp. but done
   // API used by particle transport
   void  ApplyMsc(Int_t imat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid);
   Int_t Eloss(Int_t imat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid);
   void  ProposeStep(Int_t imat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid);
   Int_t SampleDecay(Int_t ntracks, GeantTrack_v &tracksin, GeantTrack_v &tracksout);//not. imp.
   Int_t SampleInt(Int_t imat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid);
   void  GetRestFinSates(Int_t partindex, TEFstate *elemfstate, Double_t energyLimit,
            GeantTrack_v &tracks, Int_t iintrack, Int_t &nTotSecPart, Int_t tid);
   void  RotateNewTrack(Double_t oldXdir, Double_t oldYdir, Double_t oldZdir,
            GeantTrack_v &tracks, Int_t itrack);
   void  RotateTrack(GeantTrack_v &tracks, Int_t itrack, Double_t theta, Double_t phi);

   void SetIsRestProcOn(Bool_t boolval){fIsRestProcOn = boolval;}

private:
   TTabPhysMgr(const TTabPhysMgr &);//no imp.	
   TTabPhysMgr& operator=(const TTabPhysMgr &);//no imp.

   ClassDef(TTabPhysMgr,1)
};

#endif
