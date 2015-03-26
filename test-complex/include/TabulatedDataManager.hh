#ifndef TabulatedDataManager_HH
#define TabulatedDataManager_HH 1

#include "Geant/Config.h"

#include "globals.hh"
#include "TGeoManager.h"
#include "G4ParticleChange.hh"

#define MAXNELEMENTS 20 //max number of elements in one material(TMXsec)

class G4Track;
class G4DynamicParticle;

class TEXsec;
class TMXsec;
class TEFstate;
class TPDecay;
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

  //sampling element for intercation and type of intercation based on the tab. data
  //element index is the return and type will be in reactionid after termination   
  Int_t SampleInteraction(const G4int imat, const G4Track &atrack, 
                          Int_t &reactionid);
  //Get range of the particle
  G4double GetRange(const G4int imat, const G4Track &atrack);

  //sampling a tabulated final state, insert secondaries into the particle change 
  //after proper transforamtion, update properties of primary in track 
  void SampleFinalState(const Int_t elementindex, const Int_t reactionid,
                        const G4Track &atrack, G4ParticleChange *particlechange, 
                        Double_t energylimit);   
  //compute energy loss for charged particle from dE/dx and update particle change
  void EnergyLoss(G4int imat, const G4Track &atrack, const G4Step &astep, 
                  G4ParticleChange *particlechange, G4double energylimit);
  void SampleFinalStateAtRest(const Int_t imat, const G4Track &atrack,
                        G4ParticleChange *particlechange, Double_t energylimit);   
  void SampleDecayInFlight(const Int_t partindex, const G4Track &atrack, 
                       G4ParticleChange *particlechange, Double_t energylimit );
  void ApplyMsc(G4int imat, const G4Track &atrack);
  void RotateNewTrack(Double_t oldXdir, Double_t oldYdir, Double_t oldZdir, 
                G4ThreeVector &newDir);
  void RotateTrack(G4ThreeVector &newdir, Double_t costheta, Double_t phi);

  G4bool HasRestProcess(Int_t gvindex);


  TEFstate **   GetElemFstate()   { return fElemFstate; }

  static void   SetTGeomManager(TGeoManager* tgm) { fgGeom= tgm; } 
  TGeoManager*  GetTGeomManager() { return fgGeom; } 

  static void SetVerboseLevel(G4int level) {fgVerboseLevel = level;}
//#ifdef MAKESTAT
//  static unsigned long killedTracks;
//#endif
  // use or not ranges, inv. ranges and along step limit
  // will also change if the discrete interaction is selected by using the pre 
  // or the post step energy:
  //                -isUseRange = false post step Ekin is sused (def.)
  //                -isUseRange = true  pre step Ekin is sused (when comp. to G4)   
  static G4bool fgIsUseRange; 

  // get current version number
  Int_t VersionMajor() const {return fgVersion/1000/1000;}
  Int_t VersionMinor() const {return fgVersion/1000-VersionMajor()*1000;}
  Int_t VersionSub()   const {return fgVersion-VersionMajor()*1000000-VersionMinor()*1000;}
  char* GetVersion();

private:
  TabulatedDataManager(const TabulatedDataManager&);
  TabulatedDataManager& operator=(const TabulatedDataManager&);

  static TabulatedDataManager* fgInstance;

  // current version number
  static const Int_t fgVersion=1000002;

  G4int            fNelements;    // Total number of elements in the geometry
  G4int            fNmaterials;   // Total number of materials in the geometry

  TEXsec         **fElemXsec;     // Array of x-section pointers per element
  TEFstate       **fElemFstate;   // Array of final state pointers per element
  TMXsec         **fMatXsec;      // Array of x-section pointers per material 
  TPDecay         *fDecay;        // Decay tables for each particles 
  static TGeoManager *fgGeom;     // Pointer to the geometry manager   

  static const char* tStatus[];   // Possible track status for printing info. 
  static G4int fgVerboseLevel;    // Verbosity level in tabulated physics
  G4bool      *fHasNCaptureAtRest;// [#particles] 

};

#endif
