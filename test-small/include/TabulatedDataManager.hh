//
// S.Y. Jun & J. Apostolaskis, April 2014
//
#ifndef TabulatedDataManager_HH
#define TabulatedDataManager_HH 1

#include "globals.hh"
#include <TGeoManager.h>

#define MAXNELEMENTS 20 //max number of elements in one material(TMXsec)

class G4Track;

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
  TabulatedDataManager(TGeoManager* geom, const char* xsecfilename, 
		       const char* finalsfilename);	
  
  G4double GetInteractionLength(G4int imat, G4int part, G4double en);
  void SampleSecondaries(std::vector<G4Track*>* vdp, G4int imat, 
			 G4Track* atrack);

private:

  static TabulatedDataManager* theInstance;

  G4int            fNelements;    // Total number of elements in the geometry
  G4int            fNmaterials;   // Total number of materials in the geometry

  TEXsec         **fElemXsec;     // Array of x-section pointers per element
  TEFstate       **fElemFstate;   // Array of final state pointers per element
  TMXsec         **fMatXsec;      // Array of x-section pointers per material 
  TGeoManager     *fGeom;         // Pointer to the geometry manager   
};

#endif
