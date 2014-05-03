//
// S.Y. Jun & J. Apostolakis, April 2014
//

//  Create Root Materials for all materials in the Geant4 Material Table
//  Maintain an index for the correspondance between G4 -> Root material
//    and the reverse ( Root -> G4 material )
//
#ifndef MaterialConverter_hh
#define MaterialConverter_hh 1

// #include "globals.hh"
// typedef int G4int;
#include <vector>
class TGeoManager;

class MaterialConverter
{
 public:
   static MaterialConverter* Instance();
  
   void Initialize(){ CreateRootMaterials(); }
  
   void CreateRootMaterials();
    // Build the list of materials in Root, given the materials in Geant4

   // Access methods
   int GetG4Material( int iRT)
         { iRT=std::min(iRT, fMaxRootMaterial); iRT=std::max(0,iRT); return fG4MatIndices[iRT]; } 
   int GetRootMaterial( int iG4)
         { iG4=std::min(iG4, fMaxG4Material); iG4=std::max(0,iG4); return fG4MatIndices[iG4]; }
  
   static TGeoManager*  GetGeomMgr() { return fTGeomMgr;}

private:
   MaterialConverter(); 
   ~MaterialConverter(); 

  void ExpandRtIndices(int idG4);
  void ExpandG4Indices(int idRt);
  
 private:
   static TGeoManager *fTGeomMgr;    // Should be singleton 
  
   std::vector<int> fRootMatIndices; // [ key =  G4  Mat index ]
   std::vector<int> fG4MatIndices;   // [ key = Root Mat index ]
   int                fMaxG4Material;
   int                fMaxRootMaterial;
};

#endif // MaterialConverter_hh 1
