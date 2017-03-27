
#include "NISTMaterialData.h"
#include "PhysicalConstants.h"
#include "MaterialState.h"

namespace geantphysics {

NISTMaterialData& NISTMaterialData::Instance() {
   static NISTMaterialData instance;
   return instance;
}


//______________________________________________________________________________
NISTMaterialData::NISTMaterialData() {
   BuildTable();
   // Add all NIST materials to the NIST material name -> index internal map
   for (int i = 0; i < gNumberOfNISTMaterials; ++i)
     fMapNISTMaterialNameToIndex[fNISTMaterialDataTable[i].fName] = i;
}


//______________________________________________________________________________
int NISTMaterialData::FindNISTMaterialDataIndex(const std::string &name) {
  int indx = -1;
  const std::map<const std::string,int>::iterator itr = fMapNISTMaterialNameToIndex.find(name);
  if (itr != fMapNISTMaterialNameToIndex.end())
    indx = itr->second;
  return indx;
}

} //namespace geantphysics
