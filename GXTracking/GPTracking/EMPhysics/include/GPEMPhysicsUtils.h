#ifndef GPEMPhysicsUtils_H
#define GPEMPhysicsUtils_H

#include <cstdlib>
#include <cmath>
#include <ctime>
#include "GPPhysicsTable.h"
#include "GPPhysics2DVector.h"

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

void readTableAndSetSpline(GPPhysicsTable* table, const char* fname) {

  bool useSpline = true;
  readTable(table,fname);
  G4int nv = table->nPhysicsVector;
  for(int j=0; j < nv; j++){
    table->physicsVectors[j].SetSpline(useSpline);
  }
}

G4bool RetrieveSeltzerBergerData(std::ifstream& in, GPPhysics2DVector *vec2D)
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

double Rndm() 
{ //random (0,1]
#ifdef GPUNONRANDOM
  return 0.123456;
#else
  return ( static_cast<double> (rand()) )/(RAND_MAX);
#endif
}

void InitcurandState(unsigned seed=0)
{ //set the initial seed for random number generator
  if( seed==0 ) {
    time_t curtime;    
    time(&curtime);    
    srand(static_cast<unsigned> (curtime));
  }
  else {
    srand(seed);
  }
}

void CopyTrack(GXTrack *This, GXTrack *track)
{ 
  track->status = This->status;
  track->q      = This->q;
  track->x      = This->x;
  track->y      = This->y;
  track->z      = This->z;
  track->px     = This->px;
  track->py     = This->py;
  track->pz     = This->pz;
  track->E      = This->E;
  track->s      = This->s;
}

#endif
