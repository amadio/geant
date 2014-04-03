#ifndef GPConstants_hh
#define GPConstants_hh 1

#include "GPTypeDef.h"
#include "GPPhysicalConstants.h"

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157e+308
#endif

#ifndef DBL_MIN 
#define DBL_MIN 2.2250738585072014e-308
#endif

//tracker volume for magnetic field: unit in [mm]

CONSTTYPE const G4double trackerRmax = 3000.0;
CONSTTYPE const G4double trackerZmax = 3000.0;

//parameters of a simple EM calorimeter (a.k.a the CMS ECAL barrel)

//dimension:  unit in [mm]
CONSTTYPE const G4double ecalZmin = -3000.0;
CONSTTYPE const G4double ecalZmax =  3000.0;
CONSTTYPE const G4double ecalRmim =  1290.0;
CONSTTYPE const G4double ecalRmax =  1520.0;

//number of crystals / 100 = 4*number of modules
CONSTTYPE const int ecalNphi  = 360/10;
CONSTTYPE const int ecalNz    = 2*85/10;

//detector volume for magnetic field: unit in [mm]

CONSTTYPE const G4double minR =     00.0;
CONSTTYPE const G4double maxR =   9000.0;
CONSTTYPE const G4double minZ = -16000.0;
CONSTTYPE const G4double maxZ =  16000.0;

//field map array
CONSTTYPE const int nbinR = 900+1;
CONSTTYPE const int nbinZ = 2*1600+1;
CONSTTYPE const int noffZ = 1600;

//EM Physics
CONSTTYPE const int nElements = 3;
CONSTTYPE const int maximumNumberOfSecondaries = 2;
CONSTTYPE const int nPhotonProcessMax = 3;

//G4LEDATA for G4SeltzerBergerModel (Z: 1-92)
CONSTTYPE const int maxElements = 92;
CONSTTYPE const unsigned int numberOfXNodes = 32;
CONSTTYPE const unsigned int numberOfYNodes = 57;
CONSTTYPE const unsigned int numberOfSBData = 101;

//Material
CONSTTYPE const int kElementMax = 4;

CONSTTYPE const int GXFieldTrack_ncompSVEC = 6;
CONSTTYPE const int maxPhysicsVector = 2;
//const int maxPVDataVector = 71;
CONSTTYPE const int maxPVDataVector = 78;

#endif
