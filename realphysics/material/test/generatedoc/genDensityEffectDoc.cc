// compile
// g++ -std="c++11" -I../../../base/inc -I../../inc -o genDensityEffectDoc genDensityEffectDoc.cc ../../src/Isotope.cc ../../src/NISTElementData.cc ../../src/NISTElementData1.cc ../../src/Element.cc ../../src/NISTMaterialData.cc ../../src/NISTMaterialData1.cc ../../src/Material.cc ../../src/DensityEffectData.cc ../../src/DensityEffectData1.cc ../../src/MaterialProperties.cc
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include "Isotope.h"
#include "Element.h"
#include "Material.h"

#include "MaterialProperties.h"

//#include "NISTElementData.h"
//#include "NISTMaterialData.h"
#include "DensityEffectData.h"

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
  using geant::DensityEffectData;

  //
  using geant::MaterialProperties;

 std::cout
 <<"/**\n"
 <<"* @brief   Data used in the geant::DensityEffectData database.\n"
 <<"* @page    DensityEffectDataDoc\n"
 <<"* @author  M Novak, A Ribon\n"
 <<"* @date    January 2016\n"
 <<"*\n"
 <<"* The DensityEffectData internal database contains parameters to compute the density effect parameter.\n"
 <<"* These data are taken from \\cite sternheimer1984density . Mean excitation energies (i.e. I(mean) [eV]) and \n"
 <<"* material densities are available in the \\ref NISTMaterialDataDoc.\n"
 <<"*\n"
 <<"<table>"<<std::endl
 <<"<caption id=\"DensityEffectDataDoc\">Density effect parameters used in geant::MaterialProperties</caption>"<<std::endl
 //
 <<"<tr><th rowspan=\"1\">  Name <th rowspan=\"1\"> Plasma energy  [eV] <th rowspan=\"1\"> -C"
 <<"<th rowspan=\"1\"> X_0 <th rowspan=\"1\"> X_1 <th rowspan=\"1\"> a "
 <<"<th rowspan=\"1\"> m <th rowspan=\"1\"> Delta_0 ";
 //

 for (int i=0; i<278;++i) {
 std::cout<<"<tr><td> \\c "<<DensityEffectData::Instance().GetName(i)
          <<" <td> "<<DensityEffectData::Instance().GetPlasmaEnergy(i)/eV
          <<" <td> "<<DensityEffectData::Instance().GetParameterC(i)
          <<" <td> "<<DensityEffectData::Instance().GetParameterX0(i)
          <<" <td> "<<DensityEffectData::Instance().GetParameterX1(i)
          <<" <td> "<<DensityEffectData::Instance().GetParameterA(i)
          <<" <td> "<<DensityEffectData::Instance().GetParameterM(i)
          <<" <td> "<<DensityEffectData::Instance().GetParameterDelta0(i)
          <<std::endl;
}

std::cout<<"*/\n";

  return 0;
}
