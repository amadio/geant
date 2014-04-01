#ifndef GXEMPhysicsUtils_H
#define GXEMPhysicsUtils_H

#include <cmath>
#include "GPPhysicsTable.h"
#include "GXPhysics2DVector.h"

void readTable(GPPhysicsTable* table, const char* fname) {
  std::ifstream fIn;
  fIn.open(fname,std::ios::in);
  if(!fIn){
    std::cout << fname << " was not found!!!" << std::endl;
    return;
  }
  fIn >> table->nPhysicsVector;
  for(int idx=0; idx<table->nPhysicsVector; idx++){
    int vType=0;
    fIn >> vType;
    table->physicsVectors[idx].type = vType;
    fIn >> table->physicsVectors[idx].edgeMin;
    fIn >> table->physicsVectors[idx].edgeMax;
    fIn >> table->physicsVectors[idx].numberOfNodes;
    int siz=0;
    fIn >> siz;
    for(int j=0; j<siz; j++){
      fIn >> table->physicsVectors[idx].binVector[j];
      fIn >> table->physicsVectors[idx].dataVector[j];
    }// j
    G4double theEmin = table->physicsVectors[idx].binVector[0];
    table->physicsVectors[idx].dBin = 
      std::log10(table->physicsVectors[idx].binVector[1]/theEmin);
    table->physicsVectors[idx].baseBin = 
      std::log10(theEmin)/table->physicsVectors[idx].dBin;
  }// idx

}

G4bool RetrieveSeltzerBergerData(std::ifstream& in, GXPhysics2DVector *vec2D)
{
  // binning
  G4int k;
  G4int dummyX; // 32 fixed up to Z = 92
  G4int dummyY; // 57 fixed up to Z = 92
  //  in >> k >> numberOfXNodes >> numberOfYNodes;
  in >> k >> dummyX >> dummyY;
  if (in.fail())  { return false; }

  // contents
  G4double valx, valy, val;
  for(size_t i = 0; i<numberOfXNodes; ++i) {
    in >> valx;
    if (in.fail())  { return false; }
    vec2D->PutX(i,valx);
   }
  for(size_t j = 0; j<numberOfYNodes; ++j) {
    in >> valy;
    if (in.fail())  { return false; }
    vec2D->PutY(j,valy);
   }
  for(size_t j = 0; j<numberOfYNodes; ++j) {
    for(size_t i = 0; i<numberOfXNodes; ++i) {
      in >> val;
      if (in.fail())  { return false; }
      vec2D->PutValue(i, j, val);
     }
  }
  in.close();
  return true;
}

#endif
