#ifndef ROOT_TTabPhysMgr
#define ROOT_TTabPhysMgr

#include <TFile.h>
#include <TGeoExtension.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TList.h>

#define MAXNELEMENTS 20 //max number of elements in one material(TMXsec)

// Singleton-like class that handles tabulated physics data loaded for the 
// materials present in the geometry. 

class TEXsec;
class TMXsec;
class TEFstate;

class TTabPhysMgr : public TObject
{
private:
   Int_t            fNelements;     // Total number of elements in the geometry
   Int_t            fNmaterials;    // Total number of materials in the geometry	
   TEXsec         **fElemXsec;      // Array of x-section pointers per element
   TEFstate       **fElemFstate;    // Array of final state pointers per element
   TMXsec 	      **fMatXsec;	      // Array of x-section pointers per material	
   TGeoManager     *fGeom;	         // Pointer to the geometry manager   
	
   static TTabPhysMgr *fgInstance;	    // Singleton instance

public:   
   TTabPhysMgr();
   TTabPhysMgr(TGeoManager* geom, const char* xsecfilename, 
	       const char* finalsfilename);
   ~TTabPhysMgr();
   static TTabPhysMgr* Instance(TGeoManager* geom=0, const char* xsecfilename=0, 
                                   const char* finalsfilename=0);
private:
   TTabPhysMgr(const TTabPhysMgr &);//no imp.	
   TTabPhysMgr& operator=(const TTabPhysMgr &);//no imp.

   ClassDef(TTabPhysMgr,1)
};

#endif
