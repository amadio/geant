#ifndef GPConstants_hh
#define GPConstants_hh 1

#include "CLHEP/Units/PhysicalConstants.h" // From CLHEP

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157e+308
#endif

#ifndef DBL_MIN 
#define DBL_MIN 2.2250738585072014e-308
#endif

// const double kInfinity = 9.0E99;

//tracker volume for magnetic field: unit in [mm]

const double trackerRmax = 3000.0;
const double trackerZmax = 3000.0;

//detector volume for magnetic field: unit in [mm]

const double minR =     00.0;
const double maxR =   9000.0;
const double minZ = -16000.0;
const double maxZ =  16000.0;

//field map array
const int nbinR = 900+1;
const int nbinZ = 2*1600+1;
const int noffZ = 1600;


/*
const double millimeter  = 1.;
const double centimeter  = 10.*millimeter;
const double meter  = 1000.*millimeter;
const double micrometer = 1.e-6 *meter;

const double c_light   = 2.99792458e+8 * 1000./1.e+9; //m/s
const double eplus = 1. ;// positron charge
const double perCent     = 0.01 ;
const double perThousand = 0.001;
const double perMillion  = 0.000001;

const double keV  = 0.001;
const double MeV  = 1.;
const double GeV  = 1000.*MeV;
const double TeV  = 1000.*GeV;

const double tesla = 0.001;


const double mole = 1.;
const double Avogadro = 6.02214179e+23/mole;
*/

#endif
