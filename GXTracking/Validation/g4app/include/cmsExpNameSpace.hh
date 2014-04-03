#ifndef cmsExpNameSpace_H
#define cmsExpNameSpace_H 1

#include "globals.hh"

namespace cmsExp {

  //tracker volume for magnetic field: unit in [mm]

  const G4double trackerRmax = 3000.0;
  const G4double trackerZmax = 3000.0;

  //detector volume for magnetic field: unit in [mm]

  const G4double minR =      0.0;
  const G4double maxR =   9000.0;
  const G4double minZ = -16000.0;
  const G4double maxZ =  16000.0;

  //field map array
  const G4int nbinR = 900+1;
  const G4int nbinZ = 2*1600+1;
  const G4int noffZ = 1600;

  //cms prameteric field constants

  const G4double Bp_3_8T[9] = {4.24326,15.0201,3.81492,0.0178712,0.000656527,
			       2.45818,0.00778695,2.12500,1.77436};
  //Magnetic field types

  enum BFieldType {
    kNull = -1,         // Not defined
    kVolumebase,        // Volumebase grid (default for all regions)
    kParametric,        // Parameteric magnetic field inside the tracker
    kUniform,           // Uniform mangnetic field along +z
    kNumberBFieldType
  };

  struct cmsFieldMapData { int   iz; 
                           int   ir;
                           float Bz; 
                           float Br; };

  //utility functions
  void ffunkti(G4double u, G4double *ff);
}

#endif
