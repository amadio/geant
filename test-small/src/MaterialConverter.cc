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
   fMaxG4Material(0),
   fMaxRootMaterial(0)
{ 
   int defSize= 4; // Must be at least 1
   fRootMatIndices.reserve(defSize); 
   fG4MatIndices.reserve(defSize); 
   fRootMatIndices[0]= -1; 
   fG4MatIndices[0]= -1;
  
   if( fTGeomMgr == 0)
      fTGeomMgr = new TGeoManager(); // ::Instance();
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
     if( g4matIndex != imatG4 )
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
     }else{
        // G4ElementVector* elementVec= g4mat->GetElementVector();
        const G4double*  g4elemFractions= g4mat->GetFractionVector();

        TGeoMixture *tgeoMixture = new TGeoMixture(g4mat->GetName(), numElements, g4mat->GetDensity() );      
        for( int ielem = 0; ielem < numElements ; ielem++ )
        {
           const G4Element* g4elem= g4mat->GetElement(ielem);
           tgeoMixture->AddElement(g4elem->GetA(), g4elem->GetZ(), g4elemFractions[ielem]);
        }
        tgeoMaterial = tgeoMixture; 
     }

    
     // Create index for correspondence between G4 and TGeo Materials
    
     // Int_t AddMaterial(const TGeoMaterial *material);  // Creates and returns the index of new material
     Int_t rtMatIdx= fTGeomMgr->AddMaterial(tgeoMaterial);
     fRootMatIndices[imatG4]= rtMatIdx;
     ExpandG4Indices(rtMatIdx); // In case Root creates larger indices for some reason
     fG4MatIndices[rtMatIdx]= imatG4;

  }
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


#if 0
PotentialCodeFromTabXsec()
{
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
      Fatal("TabulatedDataConverter",
	    "Number of elements in %s is %d > MAXNELEMENTS=%d\n",
	    mat->GetName(),nelem,MAXNELEMENTS);
    } 
    for(Int_t iel=0; iel<nelem; ++iel) {
      Double_t ad;
      Double_t zd;
      Double_t wd;
      mat->GetElementProp(ad,zd,wd,iel);
      if (zd<1 || zd>NELEM) {
	Fatal("TabulatedDataCon2verter",
	      "In material %s found element with z=%d > NELEM=%d",
	      mat->GetName(), (Int_t)zd, NELEM);
      }
      elements.SetBitNumber(zd);
    }
  }
  fNelements = elements.CountBits();
}
#endif
