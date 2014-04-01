#ifndef GXTRACK_H
#define GXTRACK_H 

#include "GPTypeDef.h"

struct
#ifdef __CUDA_
__align__(16)
#endif
GXTrack
{
  G4int status;
  G4int id;
  G4int proc;
  G4double x; 
  G4double y;
  G4double z;
  G4double px;
  G4double py;
  G4double pz;
  G4double E;
  G4double q;
  G4double s;
} ;

struct
#ifdef __CUDA_
__align__(16)
#endif
GXTrack7
{
  G4double x; 
  G4double y;
  G4double z;
  G4double px;
  G4double py;
  G4double pz;
  G4double s;
};

#endif
