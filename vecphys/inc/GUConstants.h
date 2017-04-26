#ifndef GUCONSTANTS_H
#define GUCONSTANTS_H 1

namespace vecphys {

// unit conversion from CLHEP (MeV) to GeantV (GeV) for Energy
constexpr double EScaleToGeant4 = 1000.;  
// unit conversion from CLHEP (mm^2) to GeantV (cm^2) for Energy
constexpr double XsecScaleToGeant4 = 100.;
constexpr double invXsecScaleToGeant4 = 1.0 / XsecScaleToGeant4;

// maximum of the atomic number
constexpr int maximumZ                = 92;

// Physics2DVector (SeltzerBerger data) - valid up to maximumZ = 92
// with following number of XY nodes - may extend data up to 100
// with another set of [numberOfXNodes,numberOfYNodes] = (14,31)
constexpr unsigned int numberOfXNodes = 32;
constexpr unsigned int numberOfYNodes = 57;

// ConversionBetheHeitler::CrossSectionKernel
constexpr double a0= 8.7842e+2;
constexpr double a1=-1.9625e+3;
constexpr double a2= 1.2949e+3;
constexpr double a3=-2.0028e+2;
constexpr double a4= 1.2575e+1;
constexpr double a5=-2.8333e-1;

constexpr double b0=-1.0342e+1;
constexpr double b1= 1.7692e+1;
constexpr double b2=-8.2381   ;
constexpr double b3= 1.3063   ;
constexpr double b4=-9.0815e-2;
constexpr double b5= 2.3586e-3;

constexpr double c0=-4.5263e+2;
constexpr double c1= 1.1161e+3;
constexpr double c2=-8.6749e+2;
constexpr double c3= 2.1773e+2;
constexpr double c4=-2.0467e+1;
constexpr double c5= 6.5372e-1;

} // end namespace vecphys

#endif
