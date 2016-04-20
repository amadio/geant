//
// Simple demonstrator: how to create isotopes, elements, materials.
//
// compile
// g++ -std="c++11" -I../../management/inc -I../inc -O2 -o TestMaterial01 TestMaterial01.cc ../src/*.cc

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include "PhysicalConstants.h"

#include "Isotope.h"
#include "Element.h"
#include "Material.h"

//#include "MaterialProperties.h"


int main(){

  // Units and physical constants used in the example ()
  using geant::g;
  using geant::mg;
  using geant::mole;
  using geant::cm3;
  using geant::bar;
  using geant::pascal;
  using geant::kelvin;
  using geant::atmosphere;
  using geant::perCent;
  using geant::kUniverseMeanDensity;
  using geant::kSTPTemperature;

  // Isotope from geant namespace
  using geant::Isotope;
  // Element from geant namespace
  using geant::Element;
  // Material from the geant namespace
  using geant::Material;
  // MaterialState from geant namespace
  using geant::MaterialState;
  // MaterialProperties from geant namespace
  using geant::MaterialProperties;


  // some useful variable declaration with default values
  std::string name   = "",
              symbol = "";
  double a       = 1.*g/mole,   // molar mass in internal [weight/mole] unit
         z       = 1.,          // atomic number
         density = 1.*g/cm3,    // material density in internal [weight/length^3] unit
         abundance;             // relative abundance of the i-th isotope
  int isoZ, isoN;               // number of protons/nucleons in an isotope;
  int numcomponents,            // number of components the naterial is built up
      numatoms;                 // number of i-th atoms in the molecule (for Way 1)
                                // material construction


  // ===========================================================================
  // Examples of creating(/obtaining) isotopes from the internal data base and
  // building up an elemnt from them
  // ===========================================================================
  // creating the special isotopes
  Isotope *U5 = Isotope::GetIsotope(isoZ=92, isoN=235);
  Isotope *U8 = Isotope::GetIsotope(isoZ=92, isoN=238);
  // build up an element from the special isotopes by subsequently adding them
  // to the created element object together with their realtive abundance
  Element *elU = new Element(name="enriched Uranium",symbol="U",numcomponents=2);
  elU->AddIsotope(U5, abundance= 90.*perCent);
  elU->AddIsotope(U8, abundance= 10.*perCent);
  // build up a material from this special element by adding the (single) element
  // to the created material object
  Material *matEU = new Material(name="EnrichedUranium", 19.1*g/cm3, 1);
  matEU->AddElement(elU, 1.0);
  //____________________________________________________________________________


  // ===========================================================================
  // Examples of creating elements
  // ===========================================================================
  // # crete Hydrogen
  // natural isotopes will be created and added automatically
  a = 1.01*g/mole;
  Element* elH  = new Element(name="Hydrogen",symbol="H" , z= 1., a);


  // # create element Carbon
  // here we do not give the effective atomic mass so it will be computed
  // automatically based on the natual isotope properties i.e. mean atomic mass
  // of the isotopes weighted by their natural abundance.
  // NOTE: in this case this element will be exactly the same as we were building
  // with the Element::NISTElement() method !!!! so it can lead to duplication:
  //  a = 12.01*g/mole;
  Element* elC  = new Element(name="Carbon"  ,symbol="C" , z= 6.);

  // # create element Nitrogen
  a = 14.01*g/mole;
  Element* elN  = new Element(name="Nitrogen",symbol="N" , z= 7., a);

  // # create element Oxygen
  a = 16.00*g/mole;
  Element* elO  = new Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  // # create element Silicon
  a = 28.09*g/mole;
  Element* elSi = new Element(name="Silicon",symbol="Si" , z= 14., a);

  // # create element Iron
  a = 55.85*g/mole;
  Element* elFe = new Element(name="Iron"    ,symbol="Fe", z=26., a);
  //____________________________________________________________________________



  // ===========================================================================
  // Examples of creating simple(having only one component) materials
  // ===========================================================================
  // # create material Aluminim
  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  Material* Al = new Material(name="Aluminium"   , z=13., a, density);

  // # create material Silicon
  density = 2.3300*g/cm3;
  a = 28.085*g/mole;
  Material* Si = new Material(name="Silicon"     , z=14., a, density);

  // # create material liquidArgon
  density = 1.390*g/cm3;
  a = 39.95*g/mole;
  Material* lAr = new Material(name="liquidArgon", z=18., a, density);

  // # create material Copper
  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  Material* Cu = new Material(name="Copper"      , z=29., a, density);

  // # create material Lead
  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  Material* Pb = new Material(name="Lead"        , z=82., a, density);
  //____________________________________________________________________________



  // ===========================================================================
  // Examples of creating materials from elements by subsequently adding elements
  // together with their 'number of atoms in the molecule'.(Way 1 material
  // building i.e. based on atom count)
  // ===========================================================================
  density = 1.000*g/cm3;
  Material* H2O = new Material(name="Water", density, numcomponents=2);
  H2O->AddElement(elH, numatoms=2);
  H2O->AddElement(elO, numatoms=1);

  density = 1.032*g/cm3;
  Material* Sci = new Material(name="Scintillator", density, numcomponents=2);
  Sci->AddElement(elC, numatoms=9);
  Sci->AddElement(elH, numatoms=10);

  density = 2.200*g/cm3;
  Material* SiO2 = new Material(name="quartz", density, numcomponents=2);
  SiO2->AddElement(elSi, numatoms=1);
  SiO2->AddElement(elO , numatoms=2);
  //____________________________________________________________________________


  // ===========================================================================
  // Examples of creatig materials from elements by subsequently adding elements
  // together with their fractional mass (i.e. their ratio by weight).  (Way 2/1:
  // mixture of elements)
  // ===========================================================================
  //
  // useful variable declaration
  double massfraction = 0.;  // mass fraction of the i-th element
  density = 1.290*mg/cm3;
  Material* Air = new Material(name="Air  "  , density, numcomponents=2);
  Air->AddElement(elN, massfraction=0.7);
  Air->AddElement(elO, massfraction=0.3);
  //____________________________________________________________________________


  // ===========================================================================
  // Examples of creatig materials from elements and other matrials by
  // subsequently adding the elements and materials together with their fractio-
  // nal mass (i.e. their ratio by weight).  (Way 2/2: mixture of mixtures)
  // ===========================================================================
  density = 0.200*g/cm3;
  Material* Aerog = new Material(name="Aerogel", density, numcomponents=3);
  Aerog->AddMaterial(SiO2, massfraction=0.625);
  Aerog->AddMaterial(H2O , massfraction=0.374);
  Aerog->AddElement (elC , massfraction=0.001);
  //______________________________________________________________________________



  // ===========================================================================
  // Examples of gas in non Standard Temperature and Pressure (STP) conditions
  // ===========================================================================
  //
  // useful variable declaration
  double pressure    = 1.*atmosphere;  // pressure
  double temperature = 273.15*kelvin;  // temperature

  // Create materials in Way 1 under non STP consdition
  density     = 27.*mg/cm3;
  pressure    = 50.*atmosphere;
  temperature = 325.*kelvin;
  Material* CO2 = new Material(name="Carbonic gas", density, numcomponents=2,
                               MaterialState::kStateGas,temperature,pressure);
  CO2->AddElement(elC, numatoms=1);
  CO2->AddElement(elO, numatoms=2);

  density     = 2.67*mg/cm3;
  pressure    = 1.*atmosphere;
  temperature = 273.15*kelvin;
  Material* C4H10 = new Material(name="isobutane", density, numcomponents=2,
                                 MaterialState::kStateGas,temperature,pressure);
  C4H10->AddElement(elC, numatoms=4);
  C4H10->AddElement(elH, numatoms=10);

  // Create materials in Way 2 under non STP consdition
  density     = 0.3*mg/cm3;
  pressure    = 2.*atmosphere;
  temperature = 500.*kelvin;
  Material* steam = new Material(name="Water steam ", density, numcomponents=1,
                                 MaterialState::kStateGas,temperature,pressure);
  steam->AddMaterial(H2O, massfraction=1.);
  //______________________________________________________________________________



  // ===========================================================================
  // Examples of creating vacuum
  // ===========================================================================
  density     = kUniverseMeanDensity;
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  Material* galactic =  new Material(name="Galactic", z=1., a=1.01*g/mole, density,
                                     MaterialState::kStateGas,temperature,pressure);

  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = kSTPTemperature;
  Material* beam = new Material(name="Beam ", density, numcomponents=1,
                                MaterialState::kStateGas,temperature,pressure);
  beam->AddMaterial(Air, massfraction=1.);
  //____________________________________________________________________________



  // ===========================================================================
  // Examples of creating NIST Materials
  // ===========================================================================
  //
  // Material Aluminium has already been creted earlier so a pointer to that will
  // be returned and NIST_MAT_Al is not created. However, NIST_MAT_Al can be
  // created by changing the parameter name from "Aluminim" to "NIST_MAT_Al"
  Material *nist_mat_Al    = Material::NISTMaterial("Aluminium");
  // Material NIST_MAT_Si has not been creted earlier and we have data for this name
  // in the NISTMaterialData so that NIST material will be created and a pointer to
  // that will be returned
  Material *nist_mat_Si    = Material::NISTMaterial("NIST_MAT_Si");
  // NIST amaterial water
  Material *nist_mat_Water = Material::NISTMaterial("NIST_MAT_WATER");
  // NIST amaterial Concrete
  Material *nist_mat_Concr = Material::NISTMaterial("NIST_MAT_CONCRETE");

  // Material NIST_MAT_XX has not been creted earlier and we do NOT have data for this
  // name in the NISTMaterialData so that NIST material CANNOT be created. ERROR
  //Material *nist_X = Material::NISTMaterial("NIST_MAT_XX");


  // ===========================================================================
  // printing out the isotope table (i.e. all isotopes that were created)
  //
  std::cout<< Isotope::GetTheIsotopeTable();
  // ===========================================================================

  // ===========================================================================
  // printing out the element table (i.e. all elements that were created)
  //
  std::cout<< Element::GetTheElementTable();
  // ===========================================================================

  // ===========================================================================
  // printing out the material table (i.e. all materials that were created)
  //
  std::cout<< Material::GetTheMaterialTable();
  // ===========================================================================

  return 0;
}
