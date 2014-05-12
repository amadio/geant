//
// S.Y. Jun & J. Apostolaskis, April 2014
//
#ifndef TabulatedDataManager_HH
#define TabulatedDataManager_HH 1

#include "globals.hh"
#include <TGeoManager.h>
#include "GXTrack.hh"
#include "G4ParticleChange.hh"



#define MAXNELEMENTS 20 //max number of elements in one material(TMXsec)

class G4Track;
class G4DynamicParticle;

class TEXsec;
class TMXsec;
class TEFstate;
class TGeoMaterial;

class TabulatedDataManager
{
public:

  static TabulatedDataManager* Instance();

  TabulatedDataManager();
  ~TabulatedDataManager();	
   TabulatedDataManager(TGeoManager* geom, 
                       const char* xsecfilename, 
		       const char* finalsfilename);	
  
  G4double GetInteractionLength(G4int imat, const G4Track& atrack);
  void SampleSecondaries(std::vector<GXTrack*>* vdp, G4int imat, 
                         const G4Track* atrack, G4int ireac);

  //sampling element for intercation and type of intercation based on the tab. data
  //element index is the return and type will be in reactionid after termination   
  Int_t SampleInteraction( const G4int gvmatindex, const G4Track &atrack, 
                           Int_t &reactionid);
  //sampling a tabulated final state, insert secondaries into the particle change 
  //after proper transforamtion, update properties of primary in track 
  void SampleFinalState(const Int_t elementindex, const Int_t reactionid, 
                const G4Track &atrack, G4ParticleChange *particlechange);   
  void RotateNewTrack(Double_t oldXdir, Double_t oldYdir, Double_t oldZdir, 
                G4ThreeVector &newDir);

 
  TEFstate **   GetElemFstate()   { return fElemFstate; }

  static void   SetTGeomManager(TGeoManager* tgm) { fgGeom= tgm; } 
  TGeoManager*  GetTGeomManager() { return fgGeom; } 

private:

  static TabulatedDataManager* fgInstance;

  G4int            fNelements;    // Total number of elements in the geometry
  G4int            fNmaterials;   // Total number of materials in the geometry

  TEXsec         **fElemXsec;     // Array of x-section pointers per element
  TEFstate       **fElemFstate;   // Array of final state pointers per element
  TMXsec         **fMatXsec;      // Array of x-section pointers per material 
  static TGeoManager *fgGeom;         // Pointer to the geometry manager   
};

#endif
