//
// S.Y. Jun & J. Apostolakis, April 2014
//

//  Create Root Materials for all materials in the Geant4 Material Table
//    Maintain an index for the correspondance between G4 -> Root material
//    and the reverse ( Root -> G4 material )
//
#include "MaterialConverter.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "TGeoManager.h"
#include "TGeoMaterial.h"

// #include "assert.h"
#include <cstdio>

TGeoManager* MaterialConverter::fTGeomMgr= 0;

MaterialConverter::MaterialConverter():
   fRootMatIndices(),
   fG4MatIndices(),
   fMaxG4Material(0),
   fMaxRootMaterial(0),
   fUsedExistingGeomMgr(true),
   fInitialized(false)
{ 
   int defSize= 4; // Must be at least 1
   fRootMatIndices.reserve(defSize); 
   fG4MatIndices.reserve(defSize); 
   fRootMatIndices[0]= -1; 
   fG4MatIndices[0]= -1;
  
   if( fTGeomMgr == 0){
      fTGeomMgr = new TGeoManager(); // ::Instance();
      fUsedExistingGeomMgr= false;
   } else {
      // fTGeomMgr= existingGeomMgr;  // Must be already set - no arguments
      fUsedExistingGeomMgr= true;
   }
}

MaterialConverter::~MaterialConverter() { }

MaterialConverter* MaterialConverter::Instance()
{
   static MaterialConverter MatConvObj;
   return &MatConvObj;
}


void MaterialConverter::CreateRootMaterials()
{
  // Build the list of materials in Root, given the materials in Geant4

  const G4MaterialTable *theG4MaterialTable = G4Material::GetMaterialTable();

  int nmaterials= theG4MaterialTable->size();
  // G4MatIndex.reserve(numRootMaterials); 
  ExpandG4Indices(nmaterials); // Ensure that G4MatIndex has the necessary space (known)
  // fRootMatIndex.reserve(nmaterials); 
  ExpandRtIndices(nmaterials); // Ensure that RTIndex    has the predicted space (true max tbd)

  for(G4int imatG4=0; imatG4<nmaterials; ++imatG4) {
     G4Material *g4mat = (*theG4MaterialTable)[imatG4];

     // Check the G4 index
     size_t g4matIndex= g4mat->GetIndex();
     // assert( g4matIndex == imatG4 );
     if( g4matIndex != (unsigned long)imatG4 )
        std::cerr << "ERROR in G4 material index - "
                 << " Expected " << imatG4 << " and found " << g4matIndex
                 << std::endl;
    
     // Treat elemental materials separately
     // G4Element *g4ele= 
     size_t numElements= g4mat->GetNumberOfElements(); 

     // Create corresponding TGeoMaterial
     TGeoMaterial *tgeoMaterial;

     if( numElements == 1 ) {
        tgeoMaterial =
            new TGeoMaterial(g4mat->GetName(), g4mat->GetA(), g4mat->GetZ(), 
                             g4mat->GetDensity(),   // Units => Root Units ?
                             g4mat->GetRadlen(), 
                             g4mat->GetNuclearInterLength() );
     } else {
        // G4ElementVector* elementVec= g4mat->GetElementVector();
        const G4double*  g4elemFractions= g4mat->GetFractionVector();

        TGeoMixture *tgeoMixture = new TGeoMixture(g4mat->GetName(), numElements, g4mat->GetDensity() );      
        for( int ielem = 0; ielem < (G4int) numElements ; ielem++ )
        {
           const G4Element* g4elem= g4mat->GetElement(ielem);
           tgeoMixture->AddElement(g4elem->GetA(), g4elem->GetZ(), g4elemFractions[ielem]);
        }
        tgeoMaterial = tgeoMixture; 
     }

    
     // Create index for correspondence between G4 and TGeo Materials
    
     // int AddMaterial(const TGeoMaterial *material);  // Creates and returns the index of new material
     int rtMatIdx= fTGeomMgr->AddMaterial(tgeoMaterial);
     fRootMatIndices.push_back(rtMatIdx);
     ExpandG4Indices(rtMatIdx); // In case Root creates larger indices for some reason
     fG4MatIndices[rtMatIdx]= imatG4;

  }
  IdentifyUsedMaterials();

  fInitialized= true;
}

void MaterialConverter::ConnectG4andRootMaterials()
{
  // Check whether the materials in Geant4 already exist in Root
  // If an element already exists, reuse it
  // If a material with the same name exists, use it (do not check contents)
  // Connect all G4 materials with the corresponding Root material
  // - Reverse connects are not guaranteed.
  
  const G4MaterialTable *theG4MaterialTable = G4Material::GetMaterialTable();
  
  int nmaterials= theG4MaterialTable->size();
  // G4MatIndex.reserve(numRootMaterials);
  ExpandG4Indices(nmaterials); // Ensure that G4MatIndex has the necessary space (known)
                               // fRootMatIndex.reserve(nmaterials);
  ExpandRtIndices(nmaterials); // Ensure that RTIndex    has the predicted space (true max tbd)

  std::cout<< "\n==================MATERIAL=CONVERTER=START================="<<std::endl;  
  for(G4int imatG4=0; imatG4<nmaterials; ++imatG4) {
    std::cout<< "\n==========================================================="<<std::endl;
             
    G4Material *g4mat = (*theG4MaterialTable)[imatG4];
    
    // Check the G4 index
    size_t g4matIndex= g4mat->GetIndex();
    // assert( g4matIndex == imatG4 );
    if( g4matIndex != (unsigned long) imatG4 )
      std::cerr << "ERROR in G4 material index - "
      << " Expected " << imatG4 << " and found " << g4matIndex
      << std::endl;
    
    // Treat elemental materials separately
    // G4Element *g4ele=
    size_t numElements= g4mat->GetNumberOfElements();
    
    // Find or Create corresponding TGeoMaterial
    TGeoMaterial *tgeoMaterial;
    int rtMatIdx= -1;
      
    // TGeoMaterial     *GetMaterial(const char *matname) const;
    //tgeoMaterial= fTGeomMgr->GetMaterial(g4mat->GetName()) ;
    const G4String rootMatName = g4mat->GetName();//+'0'; 
    tgeoMaterial= fTGeomMgr->GetMaterial(rootMatName) ;
   
    if( tgeoMaterial != 0 ){
      // Find the Root index for this material ?
      rtMatIdx= tgeoMaterial->GetIndex();

      std::cout << "MaterialConverter::" << std::endl;
      std::cout << "*** FOUND Material with Name= "
                << g4mat->GetName() << "  in TGeoManager." << std::endl;
      std::cout << "*** " << "Index of " << g4mat->GetName() 
                << " in [ Geant4 , TGeoMgr ]= [ " << imatG4 << " , " 
                << rtMatIdx << " ]"  << std::endl;

      // Could check here whether the material is correct, and
      // return an error if otherwise
      
      // Could even use the following method to remove - for rectification
      // void      TGeoManager::RemoveMaterial(int index);
    }
    else
    {
      std::cout << "MaterialConverter::" << std::endl;
      std::cout << "*** Did NOT find Material with  Name= " 
                << g4mat->GetName() << " in TGeoManager. " << std::endl;
      std::cout << "*** " << g4mat->GetName() << " is composed from "<< numElements 
                << " elements." << std::endl;

      if( numElements == 1 ) {
        tgeoMaterial =
          new TGeoMaterial(rootMatName, g4mat->GetA(), g4mat->GetZ(),
                           g4mat->GetDensity(),   // Units => Root Units ?
                           g4mat->GetRadlen(),
                           g4mat->GetNuclearInterLength() );
      } else {
        // G4ElementVector* elementVec= g4mat->GetElementVector();
        const G4double*  g4elemFractions= g4mat->GetFractionVector();
        TGeoMixture *tgeoMixture = new TGeoMixture(rootMatName, numElements, 
                                                   g4mat->GetDensity() );
        for( int ielem = 0; ielem < (G4int) numElements ; ielem++ )
        {
          const G4Element* g4elem= g4mat->GetElement(ielem);
          tgeoMixture->AddElement(g4elem->GetA(), g4elem->GetZ(), g4elemFractions[ielem]);
        }
        //tgeoMaterial = tgeoMixture;
      }
      //rtMatIdx= fTGeomMgr->AddMaterial(tgeoMaterial); already added by the CTRs
      rtMatIdx = fTGeomMgr->GetMaterialIndex(rootMatName);
      std::cout << "*** " << g4mat->GetName() << " has been added to TGeoMgr!"<< std::endl;
      std::cout << "*** " << "Index of " << g4mat->GetName() << " in [ Geant4 , TGeoMgr ]= [ " 
                << imatG4 << " , " << rtMatIdx << " ]"  << std::endl;

    }
    
    // Create index for correspondence between G4 and TGeo Materials
    
    // int AddMaterial(const TGeoMaterial *material);  // Creates and returns the index of new material
    fRootMatIndices[imatG4]= rtMatIdx;
    ExpandG4Indices(rtMatIdx); // In case Root creates larger indices for some reason
    fG4MatIndices[rtMatIdx]= imatG4;
  }
  std::cout<< "\n===================MATERIAL=CONVERTER=END==================\n\n"<<std::endl;  

  // IdentifyUsedMaterials();

  fInitialized= true;
}

#include "G4LogicalVolumeStore.hh"
#include <TList.h>

//  TGeoMaterial method
//  void SetUsed(bool flag=true)
void
MaterialConverter::IdentifyUsedMaterials()
{
  G4LogicalVolumeStore* logVolStore= G4LogicalVolumeStore::GetInstance();
  // TList* MaterialList = fTGeomMgr->GetListOfMaterials();
  
  G4cout << "MaterialConverter::IdentifyUsedMaterials called." << G4endl;
  
  const int numVolumes= logVolStore->size();
  Size_t iVol;
  
  //  int matRtMaxIndex= fRootMatIndices.size();
  
  for( iVol=0;  iVol<numVolumes; iVol++ )
  {
    G4LogicalVolume* logVol= logVolStore->at(iVol);
    G4Material*      matG4= logVol->GetMaterial();
    Size_t           imatG4= matG4->GetIndex();
    int              imatR= fRootMatIndices.at(imatG4);
    TGeoMaterial*    matRt=  (TGeoMaterial*)fTGeomMgr->GetListOfMaterials()->At(imatR);
    // (TGeoMaterial*) MaterialList->At(imat);
    
    G4cout << " Material " << matG4->GetName() << " is used in " << logVol->GetName() << G4endl;
    matRt->SetUsed(true);
  }
  G4cout << "MaterialConverter::IdentifyUsedMaterials ended." << G4endl;
}

void
MaterialConverter::ExpandG4Indices(int idRt)
{
   if( idRt > fMaxRootMaterial){
      fMaxRootMaterial= idRt;
      fG4MatIndices.reserve(fMaxRootMaterial);// The G4Mat 'array' is indexed by Root material
   }
}

void
MaterialConverter::ExpandRtIndices(int idG4)
{
   if( idG4 > fMaxG4Material){
      fMaxG4Material= idG4;
      fRootMatIndices.reserve(fMaxG4Material);  // The RtMat 'array' is indexed by G4 material
   }
}

#include <TList.h>
#include <TCollection.h>
#define NELEM 118         // Total number of elements

void MaterialConverter::DumpListOfMaterials(bool /*onlyUsed*/)
{
  int noUsedMaterials=0;
  
  //Load elements from geometry
  TList *matlist = (TList*) fTGeomMgr->GetListOfMaterials();
  
  TIter next(matlist);
  TGeoMaterial *mat=0;
    
  //INFO: print number of materials in the current TGeoManager
  std::cout << "#materials:= " << matlist->GetSize() << std::endl;
  
  // First loop on all materials to mark used elements
  while((mat = (TGeoMaterial*) next())) {
    if(!mat->IsUsed() || mat->GetZ()<1.) {
      // std::cout << "  Z= " << mat->GetZ() << "  Z= " << mat->GetZ()
      printf(" Index=%3d Z=%6.1f  Name=%15s - Unused\n", mat->GetIndex(), mat->GetZ(), mat->GetName());
      // if( onlyUsed )
      continue;
    } else {
      noUsedMaterials++;
    }
    printf(" Index=%3d Z=%6.1f  Name=%15s -   Used\n", mat->GetIndex(), mat->GetZ(), mat->GetName());

    // Check the elements
    int nelem = mat->GetNelements();
 
    for(int iel=0; iel<nelem; ++iel) {
      double ad;
      double zd;
      double wd;
      mat->GetElementProp(ad,zd,wd,iel);
      if (zd<1 || zd>NELEM) {
        std::cerr << " Fatal in MaterialConverter::DumpListOfMaterials"; 
	printf( "In material %s found element with z=%d > NELEM=%d",
	      mat->GetName(), (int)zd, NELEM);
      }
    }
  }
  // int numElements = elements.CountBits();
  
  // std::cout << " Total num of elements used= " << numElements << std::endl;
  std::cout << "END of DumpListOfMaterials" << std::endl;
}
