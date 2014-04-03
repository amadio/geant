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

#include "GPTypeDef.h"

  // 
  // Length [L]
  //
  CONSTTYPE const G4double millimeter  = 1.;                        
//  CONSTTYPE const G4double millimeter2 = millimeter*millimeter;
//  CONSTTYPE const G4double millimeter3 = millimeter*millimeter*millimeter;

  CONSTTYPE const G4double centimeter  = 10.*millimeter;   
  CONSTTYPE const G4double centimeter2 = centimeter*centimeter;
  CONSTTYPE const G4double centimeter3 = centimeter*centimeter*centimeter;

  CONSTTYPE const G4double meter  = 1000.*millimeter;                  
//  CONSTTYPE const G4double meter2 = meter*meter;
//  CONSTTYPE const G4double meter3 = meter*meter*meter;

  CONSTTYPE const G4double kilometer = 1000.*meter;                   
//  CONSTTYPE const G4double kilometer2 = kilometer*kilometer;
//  CONSTTYPE const G4double kilometer3 = kilometer*kilometer*kilometer;

//  CONSTTYPE const G4double parsec = 3.0856775807e+16*meter;

  CONSTTYPE const G4double micrometer = 1.e-6 *meter;             
  CONSTTYPE const G4double  nanometer = 1.e-9 *meter;
  CONSTTYPE const G4double  angstrom  = 1.e-10*meter;
  CONSTTYPE const G4double  fermi     = 1.e-15*meter;

  CONSTTYPE const G4double      barn = 1.e-28*meter*meter;
  CONSTTYPE const G4double millibarn = 1.e-3 *barn;
  CONSTTYPE const G4double microbarn = 1.e-6 *barn;
  CONSTTYPE const G4double  nanobarn = 1.e-9 *barn;
  CONSTTYPE const G4double  picobarn = 1.e-12*barn;

  // symbols
  CONSTTYPE const G4double nm  = nanometer;                        
  CONSTTYPE const G4double um  = micrometer;                        

  CONSTTYPE const G4double mm  = millimeter;                        
//  CONSTTYPE const G4double mm2 = millimeter2;
//  CONSTTYPE const G4double mm3 = millimeter3;

  CONSTTYPE const G4double cm  = centimeter;   
  CONSTTYPE const G4double cm2 = centimeter2;
  CONSTTYPE const G4double cm3 = centimeter3;

  CONSTTYPE const G4double m  = meter;                  
//  CONSTTYPE const G4double m2 = meter2;
//  CONSTTYPE const G4double m3 = meter3;

  CONSTTYPE const G4double km  = kilometer;                   
//  CONSTTYPE const G4double km2 = kilometer2;
//  CONSTTYPE const G4double km3 = kilometer3;

//  CONSTTYPE const G4double pc = parsec;

  //
  // Angle
  //
  CONSTTYPE const G4double radian      = 1.;                  
  CONSTTYPE const G4double milliradian = 1.e-3*radian;
  CONSTTYPE const G4double degree = (3.14159265358979323846/180.0)*radian;

//  CONSTTYPE const G4double   steradian = 1.;
  
  // symbols
  CONSTTYPE const G4double rad  = radian;
  CONSTTYPE const G4double mrad = milliradian;
//  CONSTTYPE const G4double sr   = steradian;
  CONSTTYPE const G4double deg  = degree;

  //
  // Time [T]
  //
  CONSTTYPE const G4double nanosecond  = 1.;
  CONSTTYPE const G4double second      = 1.e+9 *nanosecond;
  CONSTTYPE const G4double millisecond = 1.e-3 *second;
  CONSTTYPE const G4double microsecond = 1.e-6 *second;
  CONSTTYPE const G4double  picosecond = 1.e-12*second;

  CONSTTYPE const G4double hertz = 1./second;
  CONSTTYPE const G4double kilohertz = 1.e+3*hertz;
  CONSTTYPE const G4double megahertz = 1.e+6*hertz;

  // symbols
  CONSTTYPE const G4double ns = nanosecond;
  CONSTTYPE const G4double  s = second;
  CONSTTYPE const G4double ms = millisecond;

  //
  // Electric charge [Q]
  //
  CONSTTYPE const G4double eplus = 1. ;// positron charge
  CONSTTYPE const G4double e_SI  = 1.602176487e-19;// positron charge in coulomb
  CONSTTYPE const G4double coulomb = eplus/e_SI;// coulomb = 6.24150 e+18 * eplus

  //
  // Energy [E]
  //
  CONSTTYPE const G4double megaelectronvolt = 1. ;
  CONSTTYPE const G4double     electronvolt = 1.e-6*megaelectronvolt;
  CONSTTYPE const G4double kiloelectronvolt = 1.e-3*megaelectronvolt;
  CONSTTYPE const G4double gigaelectronvolt = 1.e+3*megaelectronvolt;
  CONSTTYPE const G4double teraelectronvolt = 1.e+6*megaelectronvolt;
//  CONSTTYPE const G4double petaelectronvolt = 1.e+9*megaelectronvolt;

  CONSTTYPE const G4double joule = electronvolt/e_SI;// joule = 6.24150 e+12 * MeV

  // symbols
  CONSTTYPE const G4double MeV = megaelectronvolt;
  CONSTTYPE const G4double  eV = electronvolt;
  CONSTTYPE const G4double keV = kiloelectronvolt;
  CONSTTYPE const G4double GeV = gigaelectronvolt;
  CONSTTYPE const G4double TeV = teraelectronvolt;
//  CONSTTYPE const G4double PeV = petaelectronvolt;

  //
  // Mass [E][T^2][L^-2]
  //
  CONSTTYPE const G4double  kilogram = joule*second*second/(meter*meter);   
  CONSTTYPE const G4double      gram = 1.e-3*kilogram;
  CONSTTYPE const G4double milligram = 1.e-3*gram;

  // symbols
  CONSTTYPE const G4double  kg = kilogram;
  CONSTTYPE const G4double   g = gram;
  CONSTTYPE const G4double  mg = milligram;

  //
  // Power [E][T^-1]
  //
  CONSTTYPE const G4double watt = joule/second;// watt = 6.24150 e+3 * MeV/ns

  //
  // Force [E][L^-1]
  //
  CONSTTYPE const G4double newton = joule/meter;// newton = 6.24150 e+9 * MeV/mm

  //
  // Pressure [E][L^-3]
  //
#define pascal hep_pascal                          // a trick to avoid warnings 
  CONSTTYPE const G4double hep_pascal = newton/(meter*meter);   // pascal = 6.24150 e+3 * MeV/mm3
  CONSTTYPE const G4double bar        = 100000*pascal; // bar    = 6.24150 e+8 * MeV/mm3
  CONSTTYPE const G4double atmosphere = 101325*pascal; // atm    = 6.32420 e+8 * MeV/mm3

  //
  // Electric current [Q][T^-1]
  //
  CONSTTYPE const G4double      ampere = coulomb/second; // ampere = 6.24150 e+9 * eplus/ns
  CONSTTYPE const G4double milliampere = 1.e-3*ampere;
  CONSTTYPE const G4double microampere = 1.e-6*ampere;
  CONSTTYPE const G4double  nanoampere = 1.e-9*ampere;

  //
  // Electric potential [E][Q^-1]
  //
  CONSTTYPE const G4double megavolt = megaelectronvolt/eplus;
  CONSTTYPE const G4double kilovolt = 1.e-3*megavolt;
  CONSTTYPE const G4double     volt = 1.e-6*megavolt;

  //
  // Electric resistance [E][T][Q^-2]
  //
  CONSTTYPE const G4double ohm = volt/ampere;// ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

  //
  // Electric capacitance [Q^2][E^-1]
  //
  CONSTTYPE const G4double farad = coulomb/volt;// farad = 6.24150e+24 * eplus/Megavolt
  CONSTTYPE const G4double millifarad = 1.e-3*farad;
  CONSTTYPE const G4double microfarad = 1.e-6*farad;
  CONSTTYPE const G4double  nanofarad = 1.e-9*farad;
  CONSTTYPE const G4double  picofarad = 1.e-12*farad;

  //
  // Magnetic Flux [T][E][Q^-1]
  //
  CONSTTYPE const G4double weber = volt*second;// weber = 1000*megavolt*ns

  //
  // Magnetic Field [T][E][Q^-1][L^-2]
  //
  CONSTTYPE const G4double tesla     = volt*second/(meter*meter);// tesla =0.001*megavolt*ns/mm2

  CONSTTYPE const G4double gauss     = 1.e-4*tesla;
  CONSTTYPE const G4double kilogauss = 1.e-1*tesla;

  //
  // Inductance [T^2][E][Q^-2]
  //
  CONSTTYPE const G4double henry = weber/ampere;// henry = 1.60217e-7*MeV*(ns/eplus)**2

  //
  // Temperature
  //
  CONSTTYPE const G4double kelvin = 1.;

  //
  // Amount of substance
  //
  CONSTTYPE const G4double mole = 1.;

  //
  // Activity [T^-1]
  //
//  CONSTTYPE const G4double becquerel = 1./second ;
//  CONSTTYPE const G4double curie = 3.7e+10 * becquerel;

  //
  // Absorbed dose [L^2][T^-2]
  //
//  CONSTTYPE const G4double      gray = joule/kilogram ;
//  CONSTTYPE const G4double  kilogray = 1.e+3*gray;
//  CONSTTYPE const G4double milligray = 1.e-3*gray;
//  CONSTTYPE const G4double microgray = 1.e-6*gray;

  //
  // Luminous intensity [I]
  //
//  CONSTTYPE const G4double candela = 1.;

  //
  // Luminous flux [I]
  //
//  CONSTTYPE const G4double lumen = candela*steradian;

  //
  // Illuminance [I][L^-2]
  //
//  CONSTTYPE const G4double lux = lumen/meter2;

  //
  // Miscellaneous
  //
  CONSTTYPE const G4double perCent     = 0.01 ;
  CONSTTYPE const G4double perThousand = 0.001;
  CONSTTYPE const G4double perMillion  = 0.000001;


#endif /* HEP_SYSTEM_OF_UNITS_H */
