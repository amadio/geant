// compile
// g++ -std="c++11" -I../../../base/inc -I../../inc -o tDensityEffectData tDensityEffectData.cc ../../src/DensityEffectData.cc ../../src/DensityEffectData1.cc

#include <iostream>
#include "DensityEffectData.h"

#include "MaterialState.h"

int main(){
  //use geant::DensityEffectData
  using geant::DensityEffectData;
  // using geant::MaterialState
  using geant::MaterialState;

 //
 // we will set the requested material state to kStateUndefined in all cases now
 //
  int Z = 1;
  int indx = DensityEffectData::Instance().GetElementalIndex(Z,MaterialState::kStateUndefined);
  if (indx>0)
    std::cout<<"  Z = 1 (Hydrogen; H)  => "<< DensityEffectData::Instance().GetName(indx)<<std::endl;
  else
    std::cout<<"**** NO data for Z ="<<Z<<"  ****"<<std::endl;

  Z = 14;
  indx = DensityEffectData::Instance().GetElementalIndex(Z,MaterialState::kStateUndefined);
  if (indx>0)
    std::cout<<"  Z = 14 (Silicon; Si)  => "<< DensityEffectData::Instance().GetName(indx)<<std::endl;
  else
    std::cout<<"**** NO data for Z ="<<Z<<"  ****"<<std::endl;

  Z = 84;
  indx = DensityEffectData::Instance().GetElementalIndex(Z,MaterialState::kStateUndefined);
  if (indx>0)
    std::cout<<"  Z = 84 (Polonium; Po)  => "<< DensityEffectData::Instance().GetName(indx)<<std::endl;
  else
    std::cout<<"**** NO data for Z ="<<Z<<"  ****"<<std::endl;

  Z = 85;
  indx = DensityEffectData::Instance().GetElementalIndex(Z,MaterialState::kStateUndefined);
  if (indx>0)
    std::cout<<"  Z = 85 (Astatine; At)  => "<< DensityEffectData::Instance().GetName(indx)<<std::endl;
  else
    std::cout<<"**** NO data for Z ="<<Z<<"  ****"<<std::endl;

  Z = 86;
  indx = DensityEffectData::Instance().GetElementalIndex(Z,MaterialState::kStateUndefined);
  if (indx>0)
    std::cout<<"  Z = 86 (Radon; Rn)  => "<< DensityEffectData::Instance().GetName(indx)<<std::endl;
  else
    std::cout<<"**** NO data for Z ="<<Z<<"  ****"<<std::endl;

  Z = 87;
  indx = DensityEffectData::Instance().GetElementalIndex(Z,MaterialState::kStateUndefined);
  if (indx>0)
    std::cout<<"  Z = 87 (Francium; Fr)  => "<< DensityEffectData::Instance().GetName(indx)<<std::endl;
  else
    std::cout<<"**** NO data for Z ="<<Z<<"  ****"<<std::endl;

  Z = 88;
  indx = DensityEffectData::Instance().GetElementalIndex(Z,MaterialState::kStateUndefined);
  if (indx>0)
    std::cout<<"  Z = 88 (Radium; Ra)  => "<< DensityEffectData::Instance().GetName(indx)<<std::endl;
  else
    std::cout<<"**** NO data for Z ="<<Z<<"  ****"<<std::endl;

  Z = 97;
  indx = DensityEffectData::Instance().GetElementalIndex(Z,MaterialState::kStateUndefined);
  if (indx>0)
    std::cout<<"  Z = 97 (Berkelium; Bk)  => "<< DensityEffectData::Instance().GetName(indx)<<std::endl;
  else
    std::cout<<"**** NO data for Z ="<<Z<<"  ****"<<std::endl;

  Z = 98;
  indx = DensityEffectData::Instance().GetElementalIndex(Z,MaterialState::kStateUndefined);
  if (indx>0)
    std::cout<<"  Z = 98 (Califormium; Cf)  => "<< DensityEffectData::Instance().GetName(indx)<<std::endl;
  else
    std::cout<<"**** NO data for Z ="<<Z<<"  ****"<<std::endl;

  return 0;
}
