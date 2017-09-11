///////////////////////////////////////////////////////////////////////////////////////////
//
//			CaloDetectorConstruction.cxx
//			Created: 20 June 2017
//			Author: Ryan Schmitz
//
// Description: A (linear) calorimeter implemented using VecGeom libraries, this detector construction
// 		is fully customizable and can easily pass data about its structure to other
// 		classes.
//
// 		Customizable features include:
// 		-Number of calorimeter layers
// 		-Number of absorbers per layer
//		-Production cuts (on energy, length, or particle-specific for gammas, electrons, or positrons)
//		-YZ calorimeter cross-section
//
//		Other features like incident particle energy and particle type may be set in the
//		macro which executes this application, caloAppRP.cc, as well as in its accompanying
//		macro.
//
///////////////////////////////////////////////////////////////////////////////////////////

#include "CaloDetectorConstruction.h"
#include "GeantVDetectorConstruction.h"
#include <iostream>
#include <vector>

//Material includes
#include "Isotope.h"
#include "Element.h"
#include "Material.h"
#include "MaterialProperties.h"
#include "NISTElementData.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

//Region and productionCut includes
#include "Region.h"
#include "PhysicsParameters.h"
#include "MaterialCuts.h"
/////////////////////////////////
//VecGeom includes
#include "management/GeoManager.h"
#include "volumes/Box.h"
#include "volumes/LogicalVolume.h"
/////////////////////////////////

namespace userapplication {

CaloDetectorConstruction::~CaloDetectorConstruction()
{
}

void CaloDetectorConstruction::CreateMaterials() {

  //unit definitions
  using geant::g;
  using geant::mg;
  using geant::mole;
  using geant::cm3;
//  using geant::bar;
  using geant::pascal;
  using geant::kelvin;
  using geant::atmosphere;
  using geant::kUniverseMeanDensity;
//  using geant::perCent;
//  using geant::kSTPTemperature;

  //useful variable declarations
  std::string name   = "";
  std::string symbol = "";
  double a       = 1.*g/mole,   // molar mass in internal [weight/mole] unit
         z       = 1.,          // mean numnber of protons
         density = 1.*g/cm3;    // material density in internal [weight/length^3] unit
//         abundance;             // relative abundance of the i-th isotope
//  int isoZ, isoN;               // number of protons/nucleons in an isotope;
//  int numcomponents,            // number of components the naterial is built up
//      numatoms;                 // number of i-th atoms in the molecule (for Way 1) some useful variable declaration

  double pressure    = 1.*atmosphere;  // pressure
  double temperature = 273.15*kelvin; // temperature

  //Isotope from geantphysics namespace
  using geantphysics::Isotope;
  // Element from geantphysics namespace
  using geantphysics::Element;
  // Material from the geantphysics namespace
  using geantphysics::Material;
  // MaterialState from geantphysics namespace
  using geantphysics::MaterialState;
  using geantphysics::NISTElementData;
  using geantphysics::MaterialProperties;

    //Define vacuum
    density     = kUniverseMeanDensity;
    pressure    = 3.e-18*pascal;
    temperature = 2.73*kelvin;
    Material *Galactic = new Material(name="Galactic", z=1., a=1.01*g/mole, density,
                 MaterialState::kStateGas,temperature,pressure);

    //Explicit Pb definition
 //   density = 11.35*g/cm3;
//    a = 207.19*g/mole;
//    Material* Pb = new Material(name="Lead"     , z=82., a, density);

    //Declare calo materials
    fWorldMaterialName = "Galactic";
    fWorldMaterial     =  Galactic ;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void CaloDetectorConstruction::SetProductionCutsByEnergy(double energy){

  prodCutByLength=false;
  fGammaCut=energy*geant::eV;
  fElectronCut=energy*geant::eV;
  fPositronCut=energy*geant::eV;
}

void CaloDetectorConstruction::SetProductionCutsByLength(double length){
  fGammaCut=length*geant::mm;
  fElectronCut=length*geant::mm;
  fPositronCut=length*geant::mm;
}

//all these settings are in length
void CaloDetectorConstruction::SetDetectorGammaProductionCut(double length){  fGammaCut=length*geant::mm;}
void CaloDetectorConstruction::SetDetectorElectronProductionCut(double length){  fElectronCut=length*geant::mm;}
void CaloDetectorConstruction::SetDetectorPositronProductionCut(double length){  fPositronCut=length*geant::mm;}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void CaloDetectorConstruction::SetAbsorberMaterialName(int absNum, std::string matName){
	if (absNum<=numAbsorbers){
		userMaterial[absNum]=true;
		fAbsMaterialName[absNum]=matName;
	} else {std::cerr << "ERROR: There are too few absorbers! Increase the number of absorbers or set a different absorber's material.\n";}
}
std::string CaloDetectorConstruction::GetAbsorberMaterialName(int absNum){return fAbsMaterialName[absNum];}
geantphysics::Material* CaloDetectorConstruction::GetAbsorberMaterial(int absNum){
	if (absNum<=numAbsorbers){
		return fAbsMaterial[absNum];
	} else {std::cerr << "ERROR: There are too few absorbers! Increase the number of absorbers or set a different absorber's material.\n";
	return nullptr;
	}
}
void CaloDetectorConstruction::SetDetectorMaterials(){ //private, set each absorber to its material name
	for (int i=1; i<=numAbsorbers; i++){
		fAbsMaterial[i]=geantphysics::Material::NISTMaterial(fAbsMaterialName[i]);
	}
	fWorldMaterial=geantphysics::Material::NISTMaterial(fWorldMaterialName);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void CaloDetectorConstruction::SetNumLayers(int nLayers){
	userLayerNum=true;
	numLayers=nLayers;
}
void CaloDetectorConstruction::SetNumAbsorbers(int nAbsorbers){
	userAbsorberNum=true;
	numAbsorbers=nAbsorbers;
}
int CaloDetectorConstruction::GetNumLayers() {return numLayers;}
int CaloDetectorConstruction::GetNumAbsorbers() {return numAbsorbers;}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void CaloDetectorConstruction::SetAbsorberThickness(int absNum, double thickness){
	if (absNum<=numAbsorbers){
		userThickness[absNum]=true;
		fAbsThickness[absNum]=thickness;
	} else {std::cerr << "ERROR: There are too few absorbers! Increase the number of absorbers or set a different absorber's material.\n"; }
}
double CaloDetectorConstruction::GetAbsorberThickness(int absNum){
	if (absNum<=numAbsorbers){
		return fAbsThickness[absNum];
	} else {std::cerr << "ERROR: There are too few absorbers! Increase the number of absorbers or set a different absorber's material.\n";
	return 0;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void CaloDetectorConstruction::SetDetectorYZ(double yz){
	userCaloYZ=true;
	fCaloSizeYZ=yz;
}

double CaloDetectorConstruction::GetDetectorX() {return fLayerThickness*numLayers;}
double CaloDetectorConstruction::GetDetectorYZ() {return fCaloSizeYZ;}
double CaloDetectorConstruction::GetWorldX() {return fWorldSizeX;}
double CaloDetectorConstruction::GetWorldYZ() {return fWorldSizeYZ;}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

int CaloDetectorConstruction::GetAbsorberLogicalVolumeID(int absorber){ return fAbsLogicVolumeID[absorber]; }
int CaloDetectorConstruction::GetDetectorRegionIndex() {return fDetectorRegionIndex; }

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void CaloDetectorConstruction::CreateGeometry() {

  using geant::mm;
  using geant::cm;

  //create default detector parameters if the user hasn't defined them
  if (!userLayerNum) numLayers = 10;
  if (!userAbsorberNum) numAbsorbers = 2;
  //Abs1 = liquidArgon by default
  if (!userThickness[1]) fAbsThickness[1] = 5*mm;
  if (!userMaterial[1]) fAbsMaterialName[1]="NIST_MAT_lAr";
  //Abs2 = Lead by default
  if (!userThickness[2]) fAbsThickness[2] = 10*mm;
  if (!userMaterial[2]) fAbsMaterialName[2]="NIST_MAT_Pb";
  if (numAbsorbers>=3){
  	for(int i=3;i<=numAbsorbers;i++){
		if (userThickness[i] && userMaterial[i]) {}//good, the thickness and material are defined for this non-default absorber
		else { std::cerr << "ERROR: Invalid absorber parameters. Try using SetAbsorberThickness and SetAbsorberMaterialName for absorbers 3+\n";}
	}
  }

  std::cout << "There are "<< numLayers << " layers\n";
  std::cout << "There are "<< numAbsorbers << " absorbers\n";
  for (int p=1;p<=numAbsorbers;p++)
	std::cout << "Absorber number " << p << " has material " << fAbsMaterialName[p] << " and thickness " << fAbsThickness[p] << "cm"<< std::endl;


  if (!userCaloYZ) fCaloSizeYZ = 10*cm;
  std::cout << "The calo YZ cross-section has side length "<< fCaloSizeYZ << " cm\n";
  fLayerThickness=0;
  for (int i=1;i<=numAbsorbers;i++){fLayerThickness+=fAbsThickness[i];}
  fWorldSizeX=1.2*fLayerThickness*numLayers;
  fWorldSizeYZ=1.2*fCaloSizeYZ;
  SetDetectorMaterials();

//define regions
//for now both regions have the same cuts; can change this if a user wants, but it's not very relevant for this example
  vecgeom::Region *worldRegion = new vecgeom::Region("WorldRegion",prodCutByLength,fGammaCut,fElectronCut,fPositronCut);
  vecgeom::Region *calRegion = new vecgeom::Region("CalRegion",prodCutByLength,fGammaCut,fElectronCut,fPositronCut);
  fDetectorRegionIndex = calRegion->GetIndex();

  //define world
  vecgeom::UnplacedBox *worldUnplaced      = new vecgeom::UnplacedBox(fWorldSizeX/2, fWorldSizeYZ/2, fWorldSizeYZ/2);
  vecgeom::LogicalVolume *world = new vecgeom::LogicalVolume("world", worldUnplaced);
  world->SetRegion(worldRegion);
  world->SetMaterialPtr(fWorldMaterial);

  //define unplaced solids and unfilled logical volumes for absorbers
  vecgeom::UnplacedBox *fAbsBoxes[numAbsorbers+1];
  for (int i=1; i<=numAbsorbers; i++) {
  	fAbsBoxes[i] = new vecgeom::UnplacedBox("Absorber", fAbsThickness[i]/2,fCaloSizeYZ/2,fCaloSizeYZ/2);
  }
  vecgeom::LogicalVolume *fAbsLogic[numAbsorbers+1];

  //define layers
  vecgeom::UnplacedBox *fLayerBox = new vecgeom::UnplacedBox("Layer", fLayerThickness/2,fCaloSizeYZ/2,fCaloSizeYZ/2);
  vecgeom::LogicalVolume *fLayerLogic = new vecgeom::LogicalVolume("Layer",fLayerBox);
  fLayerLogic->SetRegion(calRegion);
  fLayerLogic->SetMaterialPtr(fWorldMaterial);

  //initialize variables for placement of absorbers
  char *volName = new char;
  double xfront=-fLayerThickness/2;
  double xcenter=0;

  //create and place absorbers
  for (int k=1; k<=numAbsorbers; k++) {
   	sprintf(volName,"abs%d",k);
        fAbsLogic[k] = new vecgeom::LogicalVolume(volName,fAbsBoxes[k]);
    	fAbsLogic[k]->SetMaterialPtr(fAbsMaterial[k]);
    	fAbsLogic[k]->SetRegion(calRegion);
	fAbsLogicVolumeID[k]=fAbsLogic[k]->id();
	std::cout << "logic volume id for " << k << " is " << fAbsLogicVolumeID[k] << std::endl;
    	xcenter = xfront+0.5*fAbsThickness[k];
    	xfront += fAbsThickness[k];
    	vecgeom::Transformation3D *calPlace = new vecgeom::Transformation3D(xcenter,0,0, 0,0,0);
    	fLayerLogic->PlaceDaughter(volName,fAbsLogic[k],calPlace);
   }

  double xfront_l=-0.5*numLayers*fLayerThickness;
  double xcenter_l=0;

  for (int j=1; j<=numLayers; j++) {
	sprintf(volName,"Layer%d",j);
        xcenter_l = xfront_l+fLayerThickness/2;
	xfront_l += fLayerThickness;
    	vecgeom::Transformation3D *layerPlace = new vecgeom::Transformation3D(xcenter_l,0,0, 0,0,0);
	world->PlaceDaughter(volName,fLayerLogic,layerPlace);
  }

  //place world volume, close geometry
  vecgeom::VPlacedVolume *w = world->Place();
  vecgeom::GeoManager::Instance().SetWorld(w);
  vecgeom::GeoManager::Instance().CloseGeometry();
}

}
