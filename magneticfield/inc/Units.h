
#ifndef FIELDUNITS_H
#define FIELDUNITS_H 1

// Temporary solution, until there is SystemOfUnits.h

#include "Constants.h"

namespace fieldUnits
{
   using fieldConstants::pi;
   
   static constexpr double gigaElectronVolt = 1.0; // native unit of energy
   static constexpr double centimeter       = 1.0; // native unit of length
   static constexpr double second           = 1.0; //

   static constexpr double eplus   = 1.0 ;         // TBC
  
   static constexpr double electronVolt = 1.0e-9 * gigaElectronVolt;

   
   static constexpr double megaElectronVolt = 0.001 * gigaElectronVolt;

   static constexpr double GeV        = gigaElectronVolt;
   static constexpr double MeV        = megaElectronVolt;

   static constexpr double nanosecond   = 1.e-9 *second;
   
   static constexpr double meter      = 100.0 * centimeter; 
   static constexpr double millimeter = 0.1 * centimeter;

   static constexpr double meter2     = meter * meter; 

   static constexpr double volt       = electronVolt/eplus;
   static constexpr double tesla      = volt*second/meter2;// tesla =0.001*megavolt*ns/mm2

  // Angle
   static const double radian      = 1.;                  
   static const double milliradian = 1.e-3*radian;
   static const double degree = (pi/180.0)*radian;
   
   // To be moved into Constants.h etc
   static constexpr double c_light  = 2.99792458e+8 * meter/second;
   static constexpr double c_light2 = c_light * c_light;
}
#endif
