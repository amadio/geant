// compile

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include "Isotope.h"
#include "Element.h"
#include "Material.h"

#include "MaterialProperties.h"

//#include "NISTElementData.h"
#include "NISTMaterialData.h"
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
  using geant::NISTMaterialData;

  //
  using geant::MaterialProperties;

 std::cout
 <<"/**\n"
 <<"* @brief   Data used in the geant::NISTMaterialData database.\n"
 <<"* @page    NISTMaterialDataDoc\n"
 <<"* @author  M Novak, A Ribon\n"
 <<"* @date    december 2015\n"
 <<"*\n"
 <<"* Data are taken from NIST database [http://physics.nist.gov/cgi-bin/Star/compos.pl] and some of them has\n"
 <<"* has already been corrected (see the comments in NISTMaterialData::BuildTable() method for more details).\n"
 <<"*\n"
 <<"<table>"<<std::endl
 <<"<caption id=\"NISTMaterialDataDoc\">NIST material data used in geant::NISTMaterialData</caption>"<<std::endl
 //
 <<"<tr><th rowspan=\"2\">  Name <th rowspan=\"2\"> Density [g/cm3] <th rowspan=\"2\"> I(mean) [eV]"
 <<"<th rowspan=\"2\"> State <th rowspan=\"2\"> Temperature [K]<th rowspan=\"2\"> Pressure [atm]"
 <<"<th rowspan=\"2\"> \\#Components <th rowspan=\"2\"> By mass fraction <br/> or <br/> \\#atoms  <th colspan=\"2\"> Composition"
 //
 <<"<tr><th> Z <th> Mass fraction <br/> OR <br/> \\#atoms"<<std::endl;

 for (int i=0; i<NISTMaterialData::Instance().GetNumberOfNISTMaterials();++i) {
 int numcomps = NISTMaterialData::Instance().GetNumberOfComponents(i);
 std::cout
 <<"<tr>"
 <<"<td rowspan=\""<<numcomps+1<<"\">\\c "<<NISTMaterialData::Instance().GetName(i);

 if(NISTMaterialData::Instance().GetDensity(i)==kUniverseMeanDensity) {
  std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<"\\c global::kUniverseMeanDensity";
 } else {
  std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<NISTMaterialData::Instance().GetDensity(i)/(g/cm3);
 }

 std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<NISTMaterialData::Instance().GetMeanExcitationEnergy(i)/eV;

 if(NISTMaterialData::Instance().GetMaterialState(i)==MaterialState::kStateSolid){
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<"\\c geant::MaterialState::kStateSolid";
 } else if(NISTMaterialData::Instance().GetMaterialState(i)==MaterialState::kStateGas){
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<"\\c geant::MaterialState::kStateGas";
 } else if(NISTMaterialData::Instance().GetMaterialState(i)==MaterialState::kStateLiquid){
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<"\\c geant::MaterialState::kStateLiquid";
 } else if(NISTMaterialData::Instance().GetMaterialState(i)==MaterialState::kStateUndefined){
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<"\\c geant::MaterialState::kStateUndefined";
 } else {std::cerr<<" ***** What the fuck is this material state in case of NISTMaterial index= "<<i<<std::endl;}

 if(NISTMaterialData::Instance().GetTemperature(i)==kSTPTemperature){
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<"\\c global::kSTPTemperature";
 } else if(NISTMaterialData::Instance().GetTemperature(i)==kNTPTemperature){
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<"\\c global::kNTPTemperature";
 } else {
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<NISTMaterialData::Instance().GetTemperature(i)/kelvin;
 }

 if(NISTMaterialData::Instance().GetPressure(i)==kSTPPressure){
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<"\\c global::kSTPPressure";
 } else {
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<NISTMaterialData::Instance().GetPressure(i)/atmosphere;
 }

 std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<numcomps;

 if(NISTMaterialData::Instance().IsToBuildByAtomCount(i)){
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<"\\c By \\c Atom \\c count";
 } else {
   std::cout<<"<td rowspan=\""<<numcomps+1<<"\">"<<"\\c By \\c Mass \\c fraction";
 }

 const double *w = NISTMaterialData::Instance().GetListOfElementFractions(i);
 const int    *z = NISTMaterialData::Instance().GetListOfElements(i);
 for(int j=0;j<numcomps;++j){
   std::cout<<std::endl<<"<tr><td> "<<z[j]<<" <td> "<<w[j];
 }
 std::cout<<std::endl;
}

std::cout<<"*/\n";

/*

   static const int ElemsZ_314[]        = {  1,   6,   7,   8,  15};
   static const double FractionsZ_314[] = {  9,   9,   2,   9,   1};
   fNISTMaterialDataTable[314].fName                 = "NIST_MAT_DNA_U"; X
   fNISTMaterialDataTable[314].fDensity              = 1 * (g/cm3);      X
   fNISTMaterialDataTable[314].fMeanExcitationEnergy = 72 * eV;          X
   fNISTMaterialDataTable[314].fTemperature          = kNTPTemperature;  X
   fNISTMaterialDataTable[314].fPressure             = kSTPPressure;     X
   fNISTMaterialDataTable[314].fNumComponents        = 5;
   fNISTMaterialDataTable[314].fState                = MaterialState::kStateSolid; X
   fNISTMaterialDataTable[314].fElementList          = ElemsZ_314;
   fNISTMaterialDataTable[314].fElementFraction      = FractionsZ_314;
   fNISTMaterialDataTable[314].fIsBuiltByAtomCount   = true;             X

 */




  return 0;
}
