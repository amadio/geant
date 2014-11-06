#ifndef ROOT_TTabPhysMgr
#define ROOT_TTabPhysMgr

#include <TGeoExtension.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>

#define MAXNELEMENTS 20 //max number of elements in one material(TMXsec)

// Singleton that handles tabulated physics data loaded for the 
// materials present in the geometry. 

class TEXsec;
class TMXsec;
class TEFstate;
class TPDecay;
class GeantTrack_v;
class GeantTrack;
class TGeoMaterial;

class TTabPhysMgr
{
private:
   Int_t            fNelements;     // Total number of elements in the geometry
   Int_t            fNmaterials;    // Total number of materials in the geometry	
   TEXsec         **fElemXsec;      // Array of x-section pointers per element
   TEFstate       **fElemFstate;    // Array of final state pointers per element
   TMXsec 	  **fMatXsec;	    // Array of x-section pointers per material	
   TPDecay         *fDecay;         // Decay tables for each particles
   TGeoManager     *fGeom;	    // Pointer to the geometry manager   
   Bool_t          *fHasNCaptureAtRest; // do the particle have nCapture at rest?

	
   static TTabPhysMgr *fgInstance;	    // Singleton instance

public:   
   TTabPhysMgr();
   TTabPhysMgr(TGeoManager* geom, const char* xsecfilename, 
	       const char* finalsfilename);
   virtual ~TTabPhysMgr();
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
   void  GetRestFinStates(Int_t partindex, TMXsec *mxs, Double_t energyLimit,
            GeantTrack_v &tracks, Int_t iintrack, Int_t &nTotSecPart, Int_t tid);
   void  SampleDecayInFlight(Int_t partindex, TMXsec *mxs, Double_t energyLimit,
            GeantTrack_v &tracks, Int_t iintrack, Int_t &nTotSecPart, Int_t tid );

   Bool_t HasRestProcess(Int_t gvindex);

   void  RotateNewTrack(Double_t oldXdir, Double_t oldYdir, Double_t oldZdir,
            GeantTrack &track);
   void  RotateNewTrack(Double_t oldXdir, Double_t oldYdir, Double_t oldZdir,
            GeantTrack_v &tracks, Int_t itrack);
   void  RotateTrack(GeantTrack &track, Double_t theta, Double_t phi);
   void  RotateTrack(GeantTrack_v &tracks, Int_t itrack, Double_t theta, Double_t phi);


   // get current version number
   Int_t VersionMajor() const {return fgVersion/1000/1000;}
   Int_t VersionMinor() const {return fgVersion/1000-VersionMajor()*1000;}
   Int_t VersionSub()   const {return fgVersion-VersionMajor()*1000000-VersionMinor()*1000;}
   char* GetVersion();

private:
   TTabPhysMgr(const TTabPhysMgr &);//no imp.	
   TTabPhysMgr& operator=(const TTabPhysMgr &);//no imp.

   // current version number
   static const Int_t fgVersion=1000001;

   ClassDef(TTabPhysMgr,1)
};

#endif
