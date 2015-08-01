// -*- C++ -*-
// $Id:$
// ----------------------------------------------------------------------
// HEP coherent system of Units
//
// This file has been provided to CLHEP by Geant4 (simulation toolkit for HEP).
//
// The basic units are :
// millimeter              (millimeter)
// nanosecond              (nanosecond)
// Mega electron Volt      (MeV)
// positron charge         (eplus)
// degree Kelvin           (kelvin)
// the amount of substance (mole)
// luminous intensity      (candela)
// radian                  (radian)
// steradian               (steradian)
//
// Below is a non exhaustive list of derived and pratical units
// (i.e. mostly the SI units).
// You can add your own units.
//
// The SI numerical value of the positron charge is defined here,
// as it is needed for conversion factor : positron charge = e_SI (coulomb)
//
// The others physical constants are defined in the header file :
//PhysicalConstants.h
//
// Authors: M.Maire, S.Giani
//
// History:
//
// 06.02.96   Created.
// 28.03.96   Added miscellaneous constants.
// 05.12.97   E.Tcherniaev: Redefined pascal (to avoid warnings on WinNT)
// 20.05.98   names: meter, second, gram, radian, degree
//            (from Brian.Lasiuk@yale.edu (STAR)). Added luminous units.
// 05.08.98   angstrom, picobarn, microsecond, picosecond, petaelectronvolt
// 01.03.01   parsec    
// 31.01.06   kilogray, milligray, microgray    
// 29.04.08   use PDG 2006 value of e_SI
// 03.11.08   use PDG 2008 value of e_SI

#ifndef HEP_SYSTEM_OF_UNITS_H
#define HEP_SYSTEM_OF_UNITS_H

#include "base/Global.h"

namespace CLHEP {

  //
  //
  //
  const double     pi  = 3.14159265358979323846;
  const double  twopi  = 2*pi;
  const double halfpi  = pi/2;
  const double     pi2 = pi*pi;

  // 
  // Length [L]
  //
  const double millimeter  = 1.;                        
  const double millimeter2 = millimeter*millimeter;
  const double millimeter3 = millimeter*millimeter*millimeter;

  const double centimeter  = 10.*millimeter;   
  const double centimeter2 = centimeter*centimeter;
  const double centimeter3 = centimeter*centimeter*centimeter;

  const double meter  = 1000.*millimeter;                  
  const double meter2 = meter*meter;
  const double meter3 = meter*meter*meter;

  const double kilometer = 1000.*meter;                   
  const double kilometer2 = kilometer*kilometer;
  const double kilometer3 = kilometer*kilometer*kilometer;

  const double parsec = 3.0856775807e+16*meter;

  const double micrometer = 1.e-6 *meter;             
  const double  nanometer = 1.e-9 *meter;
  const double  angstrom  = 1.e-10*meter;
  const double  fermi     = 1.e-15*meter;

  const double      barn = 1.e-28*meter2;
  const double millibarn = 1.e-3 *barn;
  const double microbarn = 1.e-6 *barn;
  const double  nanobarn = 1.e-9 *barn;
  const double  picobarn = 1.e-12*barn;

  // symbols
  const double nm  = nanometer;                        
  const double um  = micrometer;                        

  const double mm  = millimeter;                        
  const double mm2 = millimeter2;
  const double mm3 = millimeter3;

  const double cm  = centimeter;   
  const double cm2 = centimeter2;
  const double cm3 = centimeter3;

  const double m  = meter;                  
  const double m2 = meter2;
  const double m3 = meter3;

  const double km  = kilometer;                   
  const double km2 = kilometer2;
  const double km3 = kilometer3;

  const double pc = parsec;

  //
  // Angle
  //
  const double radian      = 1.;                  
  const double milliradian = 1.e-3*radian;
  const double degree = (pi/180.0)*radian;

  const double   steradian = 1.;
  
  // symbols
  const double rad  = radian;
  const double mrad = milliradian;
  const double sr   = steradian;
  const double deg  = degree;

  //
  // Time [T]
  //
  const double nanosecond  = 1.;
  const double second      = 1.e+9 *nanosecond;
  const double millisecond = 1.e-3 *second;
  const double microsecond = 1.e-6 *second;
  const double  picosecond = 1.e-12*second;

  const double hertz = 1./second;
  const double kilohertz = 1.e+3*hertz;
  const double megahertz = 1.e+6*hertz;

  // symbols
  const double ns = nanosecond;
  const double  s = second;
  const double ms = millisecond;

  //
  // Electric charge [Q]
  //
  const double eplus = 1. ;// positron charge
  const double e_SI  = 1.602176487e-19;// positron charge in coulomb
  const double coulomb = eplus/e_SI;// coulomb = 6.24150 e+18 * eplus

  //
  // Energy [E]
  //
  const double megaelectronvolt = 1. ;
  const double     electronvolt = 1.e-6*megaelectronvolt;
  const double kiloelectronvolt = 1.e-3*megaelectronvolt;
  const double gigaelectronvolt = 1.e+3*megaelectronvolt;
  const double teraelectronvolt = 1.e+6*megaelectronvolt;
  const double petaelectronvolt = 1.e+9*megaelectronvolt;

  const double joule = electronvolt/e_SI;// joule = 6.24150 e+12 * MeV

  // symbols
  const double MeV = megaelectronvolt;
  const double  eV = electronvolt;
  const double keV = kiloelectronvolt;
  const double GeV = gigaelectronvolt;
  const double TeV = teraelectronvolt;
  const double PeV = petaelectronvolt;

  //
  // Mass [E][T^2][L^-2]
  //
  const double  kilogram = joule*second*second/(meter*meter);   
  const double      gram = 1.e-3*kilogram;
  const double milligram = 1.e-3*gram;

  // symbols
  const double  kg = kilogram;
  const double   g = gram;
  const double  mg = milligram;

  //
  // Power [E][T^-1]
  //
  const double watt = joule/second;// watt = 6.24150 e+3 * MeV/ns

  //
  // Force [E][L^-1]
  //
  const double newton = joule/meter;// newton = 6.24150 e+9 * MeV/mm

  //
  // Pressure [E][L^-3]
  //
#define pascal hep_pascal                          // a trick to avoid warnings 
  const double hep_pascal = newton/m2;   // pascal = 6.24150 e+3 * MeV/mm3
  const double bar        = 100000*pascal; // bar    = 6.24150 e+8 * MeV/mm3
  const double atmosphere = 101325*pascal; // atm    = 6.32420 e+8 * MeV/mm3

  //
  // Electric current [Q][T^-1]
  //
  const double      ampere = coulomb/second; // ampere = 6.24150 e+9 * eplus/ns
  const double milliampere = 1.e-3*ampere;
  const double microampere = 1.e-6*ampere;
  const double  nanoampere = 1.e-9*ampere;

  //
  // Electric potential [E][Q^-1]
  //
  const double megavolt = megaelectronvolt/eplus;
  const double kilovolt = 1.e-3*megavolt;
  const double     volt = 1.e-6*megavolt;

  //
  // Electric resistance [E][T][Q^-2]
  //
  const double ohm = volt/ampere;// ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

  //
  // Electric capacitance [Q^2][E^-1]
  //
  const double farad = coulomb/volt;// farad = 6.24150e+24 * eplus/Megavolt
  const double millifarad = 1.e-3*farad;
  const double microfarad = 1.e-6*farad;
  const double  nanofarad = 1.e-9*farad;
  const double  picofarad = 1.e-12*farad;

  //
  // Magnetic Flux [T][E][Q^-1]
  //
  const double weber = volt*second;// weber = 1000*megavolt*ns

  //
  // Magnetic Field [T][E][Q^-1][L^-2]
  //
  const double tesla     = volt*second/meter2;// tesla =0.001*megavolt*ns/mm2

  const double gauss     = 1.e-4*tesla;
  const double kilogauss = 1.e-1*tesla;

  //
  // Inductance [T^2][E][Q^-2]
  //
  const double henry = weber/ampere;// henry = 1.60217e-7*MeV*(ns/eplus)**2

  //
  // Temperature
  //
  const double kelvin = 1.;

  //
  // Amount of substance
  //
  const double mole = 1.;

  //
  // Activity [T^-1]
  //
  const double becquerel = 1./second ;
  const double curie = 3.7e+10 * becquerel;
  const double kilobecquerel = 1.e+3*becquerel;
  const double megabecquerel = 1.e+6*becquerel;
  const double gigabecquerel = 1.e+9*becquerel;
  const double millicurie = 1.e-3*curie;
  const double microcurie = 1.e-6*curie;
  const double Bq = becquerel;
  const double kBq = kilobecquerel;
  const double MBq = megabecquerel;
  const double GBq = gigabecquerel;
  const double Ci = curie;
  const double mCi = millicurie;
  const double uCi = microcurie;

  //
  // Absorbed dose [L^2][T^-2]
  //
  const double      gray = joule/kilogram ;
  const double  kilogray = 1.e+3*gray;
  const double milligray = 1.e-3*gray;
  const double microgray = 1.e-6*gray;

  //
  // Luminous intensity [I]
  //
  const double candela = 1.;

  //
  // Luminous flux [I]
  //
  const double lumen = candela*steradian;

  //
  // Illuminance [I][L^-2]
  //
  const double lux = lumen/meter2;

  //
  // Miscellaneous
  //
  const double perCent     = 0.01 ;
  const double perThousand = 0.001;
  const double perMillion  = 0.000001;

}  // namespace CLHEP

#endif /* HEP_SYSTEM_OF_UNITS_H */
