#ifndef EMPHYSICS_GPU_H
#define EMPHYSICS_GPU_H

#include "GPTrack.h"
#include "GPPhysicsTable.h"
#include <cmath>

// EMPhysics_kernel wrapper
void EMPhysics_gpu(GPTrack *track, size_t nTrackSize,
		    GPPhysicsTable* eBrem_table, 
		    GPPhysicsTable* eIoni_table, 
		    GPPhysicsTable* msc_table, 
		    bool useIntegral, bool useLambdaTable,
		    int nBlocks,
		    int nThreads); 

void EMPhysics_cpu(GPTrack *trackIn, size_t nTrackSize,
		    GPPhysicsTable* eBrem_table,
		    GPPhysicsTable* eIoni_table,
		    GPPhysicsTable* msc_table,
		    bool useIntegral, bool useLambdaTable);


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
    table->physicsVectors[idx].dBin = std::log10(table->physicsVectors[idx].binVector[1]/theEmin);
    table->physicsVectors[idx].baseBin = std::log10(theEmin)/table->physicsVectors[idx].dBin;
  }// idx

}

#endif
