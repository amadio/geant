// compile
//  g++ -std="c++11" -I../../../base/inc -I../../inc -o genNISTElementDoc genNISTElementDataDoc.cc ../../src/Isotope.cc ../../src/NISTElementData.cc ../../src/NISTElementData1.cc ../../src/Element.cc ../../src/NISTMaterialData.cc ../../src/NISTMaterialData1.cc ../../src/Material.cc ../../src/DensityEffectData.cc ../../src/DensityEffectData1.cc ../../src/MaterialProperties.cc
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include "Isotope.h"
#include "Element.h"
#include "Material.h"

#include "MaterialProperties.h"

#include "NISTElementData.h"
//#include "NISTMaterialData.h"
//using geant::NISTMaterialData;

#include "PhysicalConstants.h"

int main(){
  // Unist and physical constants used in the example
  using global::g;
  using global::mg;
  using global::eV;
  using global::mole;
  using global::cm3;
  using global::bar;
  using global::pascal;
  using global::kelvin;
  using global::atmosphere;
  using global::perCent;
  using global::kUniverseMeanDensity;
  using global::kSTPTemperature;
  using global::kNTPTemperature;
  using global::kSTPPressure;

  // Isotope from geant namespace
  using geant::Isotope;
  // Element from geant namespace
  using geant::Element;
  // Material from the geant namespace
  using geant::Material;
  // MaterialState from geant namespace
  using geant::MaterialState;
  //
  using geant::NISTElementData;

  //
  using geant::MaterialProperties;

 std::cout
 <<"/**\n"
 <<"* @brief   Data used in the geant::NISTElementData database.\n"
 <<"* @page    NISTElementDataDoc\n"
 <<"* @author  M Novak, A Ribon\n"
 <<"* @date    december 2015\n"
 <<"*\n"
 <<"* Atomic weights and isotopic composition data are taken from NIST database \n"
 <<"* \\cite nist2015isotopcomp (data were taken 10 December 2015). Total binding energy of \n"
 <<"* the electrons are computed based on the parametrisation given by Lunney et. \n"
 <<"* al. \\cite lunney2003recent. Isotope mass is the mass of the bare nucleus. Some of the isotope \n"
 <<"* have been changed such that the corresponding isotope masses are equal to \n"
 <<"* the rest mass of Deuteron, Triton, He3 and Alpha in the internal database \n"
 <<"* i.e. consistency.\\n\n"
 <<"*\n"
 <<"*\\n\n"
 <<"<table>"<<std::endl
 <<"<caption id=\"NISTElementDataDoc\">NIST Atomic Weights and Isotopic Compositions"
 <<" data used in geant::NISTElementData</caption>"<<std::endl
 //
 <<"<tr><th rowspan=\"2\">  Elem symbol <th rowspan=\"2\"> Atomic number (Z)"
 <<"<th rowspan=\"2\"> Number of isotopes <th colspan=\"5\"> Composition"
 //
 <<"<tr><th> Number of nucleons (N) <th> Natural abundance <th> Atomic mass [u] <th> Isotope mass [GeV] <th> Total e- b.e. [keV]"<<std::endl;

 for (int i=0; i<NISTElementData::Instance().GetNumberOfNISTElements();++i) {
 int Z = i+1;
 int numcomps = NISTElementData::Instance().GetNumberOfIsotopes(Z);
 std::cout
 <<"<tr>"
 <<"<td align=\"center\" rowspan=\""<<numcomps+1<<"\">\\c "<<NISTElementData::Instance().GetElementSymbol(Z)
 <<"<td align=\"center\" rowspan=\""<<numcomps+1<<"\">"<<Z
 <<"<td align=\"center\" rowspan=\""<<numcomps+1<<"\">"<<numcomps;
 const int* Ns    = NISTElementData::Instance().GetIsotopeNucleonNums(Z);
 const double* Ws = NISTElementData::Instance().GetIsotopeNaturalAbundances(Z);
 const double* As = NISTElementData::Instance().GetIsotopeAtomicMasses(Z);

 for(int j=0;j<numcomps;++j){
   std::cout<<std::endl<<"<tr><td align=\"center\"> "<<Ns[j]<<" <td align=\"right\" > "<<Ws[j]
                       <<" <td align=\"right\" > "<<As[j]/global::kAtomicMassUnit
                       <<" <td align=\"right\" > "<<NISTElementData::Instance().GetIsotopeMass(Z,Ns[j])/global::GeV
                       <<" <td align=\"right\" > "<<NISTElementData::Instance().GetBindingEnergy(Z,Ns[j])/global::keV;

 }
 std::cout<<std::endl;
}

std::cout<<"*/\n";





  return 0;
}
