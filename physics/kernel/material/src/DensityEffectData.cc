
#include "Geant/DensityEffectData.h"
#include "Geant/PhysicalConstants.h"
#include "Geant/MaterialState.h"

#include <iostream>

namespace geantphysics {


DensityEffectData& DensityEffectData::Instance() {
   static DensityEffectData instance;
   return instance;
}


DensityEffectData::DensityEffectData() {
   BuildTable();
   // Add all material name from the density effect data DB to the name -> index internal map
   for (int i=0; i<gNumberOfDensityEffectData; ++i) {
     fMapMaterialNameToDenistyEffectDataIndex[fDensityEffectDataTable[i].fName] = i;
   }
}


// We have data for simple, elemental materials for Z = 1-gMaxElementalZet and
// data for Z=85 and Z=87 are not available
int DensityEffectData::GetElementalIndex(int z, MaterialState state) const {
  int indx = -1;
  if (z>0 && z<=gMaxElementalZet) {
    if (z==85 || z==87)
      return indx;
    indx = z;
    if (z>85) --indx;
    if (z>87) --indx;
  }
  if (indx<0) // we don't have data for this elemental material in the DB
    return indx;
  // we have data for this elemental material in the DB so check if the material
  // state the same as requested:
  // - if the stored state = requested or the requested = undefined we are fine
  // - retrun -1 otherwise because the state are diferent
  if (state == MaterialState::kStateUndefined || state == fDensityEffectDataTable[indx].fState)
    return indx;
  return -1;
}


int DensityEffectData::FindDensityEffectDataIndex(const std::string &name) {
  int indx = -1;
  const Map_t<std::string,int>::iterator itr = fMapMaterialNameToDenistyEffectDataIndex.find(name);
  if (itr!=fMapMaterialNameToDenistyEffectDataIndex.end())
    indx = itr->second;
  return indx;
}


} // namespace geantphysics
