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
  VECPHYS_GLOBAL Precision     pi  = 3.14159265358979323846;
  VECPHYS_GLOBAL Precision  twopi  = 2*pi;
  VECPHYS_GLOBAL Precision halfpi  = pi/2;
  VECPHYS_GLOBAL Precision     pi2 = pi*pi;

  // 
  // Length [L]
  //
  VECPHYS_GLOBAL Precision millimeter  = 1.;                        
  VECPHYS_GLOBAL Precision millimeter2 = millimeter*millimeter;
  VECPHYS_GLOBAL Precision millimeter3 = millimeter*millimeter*millimeter;

  VECPHYS_GLOBAL Precision centimeter  = 10.*millimeter;   
  VECPHYS_GLOBAL Precision centimeter2 = centimeter*centimeter;
  VECPHYS_GLOBAL Precision centimeter3 = centimeter*centimeter*centimeter;

  VECPHYS_GLOBAL Precision meter  = 1000.*millimeter;                  
  VECPHYS_GLOBAL Precision meter2 = meter*meter;
  VECPHYS_GLOBAL Precision meter3 = meter*meter*meter;

  VECPHYS_GLOBAL Precision kilometer = 1000.*meter;                   
  VECPHYS_GLOBAL Precision kilometer2 = kilometer*kilometer;
  VECPHYS_GLOBAL Precision kilometer3 = kilometer*kilometer*kilometer;

  VECPHYS_GLOBAL Precision parsec = 3.0856775807e+16*meter;

  VECPHYS_GLOBAL Precision micrometer = 1.e-6 *meter;             
  VECPHYS_GLOBAL Precision  nanometer = 1.e-9 *meter;
  VECPHYS_GLOBAL Precision  angstrom  = 1.e-10*meter;
  VECPHYS_GLOBAL Precision  fermi     = 1.e-15*meter;

  VECPHYS_GLOBAL Precision      barn = 1.e-28*meter2;
  VECPHYS_GLOBAL Precision millibarn = 1.e-3 *barn;
  VECPHYS_GLOBAL Precision microbarn = 1.e-6 *barn;
  VECPHYS_GLOBAL Precision  nanobarn = 1.e-9 *barn;
  VECPHYS_GLOBAL Precision  picobarn = 1.e-12*barn;

  // symbols
  VECPHYS_GLOBAL Precision nm  = nanometer;                        
  VECPHYS_GLOBAL Precision um  = micrometer;                        

  VECPHYS_GLOBAL Precision mm  = millimeter;                        
  VECPHYS_GLOBAL Precision mm2 = millimeter2;
  VECPHYS_GLOBAL Precision mm3 = millimeter3;

  VECPHYS_GLOBAL Precision cm  = centimeter;   
  VECPHYS_GLOBAL Precision cm2 = centimeter2;
  VECPHYS_GLOBAL Precision cm3 = centimeter3;

  VECPHYS_GLOBAL Precision m  = meter;                  
  VECPHYS_GLOBAL Precision m2 = meter2;
  VECPHYS_GLOBAL Precision m3 = meter3;

  VECPHYS_GLOBAL Precision km  = kilometer;                   
  VECPHYS_GLOBAL Precision km2 = kilometer2;
  VECPHYS_GLOBAL Precision km3 = kilometer3;

  VECPHYS_GLOBAL Precision pc = parsec;

  //
  // Angle
  //
  VECPHYS_GLOBAL Precision radian      = 1.;                  
  VECPHYS_GLOBAL Precision milliradian = 1.e-3*radian;
  VECPHYS_GLOBAL Precision degree = (pi/180.0)*radian;

  VECPHYS_GLOBAL Precision   steradian = 1.;
  
  // symbols
  VECPHYS_GLOBAL Precision rad  = radian;
  VECPHYS_GLOBAL Precision mrad = milliradian;
  VECPHYS_GLOBAL Precision sr   = steradian;
  VECPHYS_GLOBAL Precision deg  = degree;

  //
  // Time [T]
  //
  VECPHYS_GLOBAL Precision nanosecond  = 1.;
  VECPHYS_GLOBAL Precision second      = 1.e+9 *nanosecond;
  VECPHYS_GLOBAL Precision millisecond = 1.e-3 *second;
  VECPHYS_GLOBAL Precision microsecond = 1.e-6 *second;
  VECPHYS_GLOBAL Precision  picosecond = 1.e-12*second;

  VECPHYS_GLOBAL Precision hertz = 1./second;
  VECPHYS_GLOBAL Precision kilohertz = 1.e+3*hertz;
  VECPHYS_GLOBAL Precision megahertz = 1.e+6*hertz;

  // symbols
  VECPHYS_GLOBAL Precision ns = nanosecond;
  VECPHYS_GLOBAL Precision  s = second;
  VECPHYS_GLOBAL Precision ms = millisecond;

  //
  // Electric charge [Q]
  //
  VECPHYS_GLOBAL Precision eplus = 1. ;// positron charge
  VECPHYS_GLOBAL Precision e_SI  = 1.602176487e-19;// positron charge in coulomb
  VECPHYS_GLOBAL Precision coulomb = eplus/e_SI;// coulomb = 6.24150 e+18 * eplus

  //
  // Energy [E]
  //
  VECPHYS_GLOBAL Precision megaelectronvolt = 1. ;
  VECPHYS_GLOBAL Precision     electronvolt = 1.e-6*megaelectronvolt;
  VECPHYS_GLOBAL Precision kiloelectronvolt = 1.e-3*megaelectronvolt;
  VECPHYS_GLOBAL Precision gigaelectronvolt = 1.e+3*megaelectronvolt;
  VECPHYS_GLOBAL Precision teraelectronvolt = 1.e+6*megaelectronvolt;
  VECPHYS_GLOBAL Precision petaelectronvolt = 1.e+9*megaelectronvolt;

  VECPHYS_GLOBAL Precision joule = electronvolt/e_SI;// joule = 6.24150 e+12 * MeV

  // symbols
  VECPHYS_GLOBAL Precision MeV = megaelectronvolt;
  VECPHYS_GLOBAL Precision  eV = electronvolt;
  VECPHYS_GLOBAL Precision keV = kiloelectronvolt;
  VECPHYS_GLOBAL Precision GeV = gigaelectronvolt;
  VECPHYS_GLOBAL Precision TeV = teraelectronvolt;
  VECPHYS_GLOBAL Precision PeV = petaelectronvolt;

  //
  // Mass [E][T^2][L^-2]
  //
  VECPHYS_GLOBAL Precision  kilogram = joule*second*second/(meter*meter);   
  VECPHYS_GLOBAL Precision      gram = 1.e-3*kilogram;
  VECPHYS_GLOBAL Precision milligram = 1.e-3*gram;

  // symbols
  VECPHYS_GLOBAL Precision  kg = kilogram;
  VECPHYS_GLOBAL Precision   g = gram;
  VECPHYS_GLOBAL Precision  mg = milligram;

  //
  // Power [E][T^-1]
  //
  VECPHYS_GLOBAL Precision watt = joule/second;// watt = 6.24150 e+3 * MeV/ns

  //
  // Force [E][L^-1]
  //
  VECPHYS_GLOBAL Precision newton = joule/meter;// newton = 6.24150 e+9 * MeV/mm

  //
  // Pressure [E][L^-3]
  //
#define pascal hep_pascal                          // a trick to avoid warnings 
  VECPHYS_GLOBAL Precision hep_pascal = newton/m2;   // pascal = 6.24150 e+3 * MeV/mm3
  VECPHYS_GLOBAL Precision bar        = 100000*pascal; // bar    = 6.24150 e+8 * MeV/mm3
  VECPHYS_GLOBAL Precision atmosphere = 101325*pascal; // atm    = 6.32420 e+8 * MeV/mm3

  //
  // Electric current [Q][T^-1]
  //
  VECPHYS_GLOBAL Precision      ampere = coulomb/second; // ampere = 6.24150 e+9 * eplus/ns
  VECPHYS_GLOBAL Precision milliampere = 1.e-3*ampere;
  VECPHYS_GLOBAL Precision microampere = 1.e-6*ampere;
  VECPHYS_GLOBAL Precision  nanoampere = 1.e-9*ampere;

  //
  // Electric potential [E][Q^-1]
  //
  VECPHYS_GLOBAL Precision megavolt = megaelectronvolt/eplus;
  VECPHYS_GLOBAL Precision kilovolt = 1.e-3*megavolt;
  VECPHYS_GLOBAL Precision     volt = 1.e-6*megavolt;

  //
  // Electric resistance [E][T][Q^-2]
  //
  VECPHYS_GLOBAL Precision ohm = volt/ampere;// ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

  //
  // Electric capacitance [Q^2][E^-1]
  //
  VECPHYS_GLOBAL Precision farad = coulomb/volt;// farad = 6.24150e+24 * eplus/Megavolt
  VECPHYS_GLOBAL Precision millifarad = 1.e-3*farad;
  VECPHYS_GLOBAL Precision microfarad = 1.e-6*farad;
  VECPHYS_GLOBAL Precision  nanofarad = 1.e-9*farad;
  VECPHYS_GLOBAL Precision  picofarad = 1.e-12*farad;

  //
  // Magnetic Flux [T][E][Q^-1]
  //
  VECPHYS_GLOBAL Precision weber = volt*second;// weber = 1000*megavolt*ns

  //
  // Magnetic Field [T][E][Q^-1][L^-2]
  //
  VECPHYS_GLOBAL Precision tesla     = volt*second/meter2;// tesla =0.001*megavolt*ns/mm2

  VECPHYS_GLOBAL Precision gauss     = 1.e-4*tesla;
  VECPHYS_GLOBAL Precision kilogauss = 1.e-1*tesla;

  //
  // Inductance [T^2][E][Q^-2]
  //
  VECPHYS_GLOBAL Precision henry = weber/ampere;// henry = 1.60217e-7*MeV*(ns/eplus)**2

  //
  // Temperature
  //
  VECPHYS_GLOBAL Precision kelvin = 1.;

  //
  // Amount of substance
  //
  VECPHYS_GLOBAL Precision mole = 1.;

  //
  // Activity [T^-1]
  //
  VECPHYS_GLOBAL Precision becquerel = 1./second ;
  VECPHYS_GLOBAL Precision curie = 3.7e+10 * becquerel;
  VECPHYS_GLOBAL Precision kilobecquerel = 1.e+3*becquerel;
  VECPHYS_GLOBAL Precision megabecquerel = 1.e+6*becquerel;
  VECPHYS_GLOBAL Precision gigabecquerel = 1.e+9*becquerel;
  VECPHYS_GLOBAL Precision millicurie = 1.e-3*curie;
  VECPHYS_GLOBAL Precision microcurie = 1.e-6*curie;
  VECPHYS_GLOBAL Precision Bq = becquerel;
  VECPHYS_GLOBAL Precision kBq = kilobecquerel;
  VECPHYS_GLOBAL Precision MBq = megabecquerel;
  VECPHYS_GLOBAL Precision GBq = gigabecquerel;
  VECPHYS_GLOBAL Precision Ci = curie;
  VECPHYS_GLOBAL Precision mCi = millicurie;
  VECPHYS_GLOBAL Precision uCi = microcurie;

  //
  // Absorbed dose [L^2][T^-2]
  //
  VECPHYS_GLOBAL Precision      gray = joule/kilogram ;
  VECPHYS_GLOBAL Precision  kilogray = 1.e+3*gray;
  VECPHYS_GLOBAL Precision milligray = 1.e-3*gray;
  VECPHYS_GLOBAL Precision microgray = 1.e-6*gray;

  //
  // Luminous intensity [I]
  //
  VECPHYS_GLOBAL Precision candela = 1.;

  //
  // Luminous flux [I]
  //
  VECPHYS_GLOBAL Precision lumen = candela*steradian;

  //
  // Illuminance [I][L^-2]
  //
  VECPHYS_GLOBAL Precision lux = lumen/meter2;

  //
  // Miscellaneous
  //
  VECPHYS_GLOBAL Precision perCent     = 0.01 ;
  VECPHYS_GLOBAL Precision perThousand = 0.001;
  VECPHYS_GLOBAL Precision perMillion  = 0.000001;

}  // namespace CLHEP

#endif /* HEP_SYSTEM_OF_UNITS_H */
