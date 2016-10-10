#include "NISTElementData.h"
#include "PhysicalConstants.h"

#include <iostream>
#include <iomanip>

namespace geantphysics {

NISTElementData& NISTElementData::Instance(){
   static NISTElementData instance;
   return instance;
}

NISTElementData::NISTElementData() {
   BuildTable();
   // init list of NIST element indices(indices in the global element table) that has
   // been already built to -1
   for (int i = 0; i < gNumberOfNISTElements; ++i)
     fIndicesOfBuiltNISTElements[i] = -1;
}

// we compute from the isotope mass to conserve consistency with some particle masses
double NISTElementData::GetAtomicMass(int Z, int N) {
   using geant::kElectronMassC2;
   using geant::kAvogadro;
   using geant::kCLightSquare;
   constexpr double unitconv = kAvogadro/kCLightSquare;

   double theMass = 0.0;
   if (Z>0 && Z<=gNumberOfNISTElements) {
     int numisos = fNISTElementDataTable[Z-1].fNumOfIsotopes;
     int indxN   = N - fNISTElementDataTable[Z-1].fNIsos[0];
     if (indxN>=0 && indxN<numisos)
       theMass = fNISTElementDataTable[Z-1].fMassIsos[indxN] + Z*kElectronMassC2
                 - fBindingEnergies[Z-1];
       theMass *= unitconv; // convert energy to [weight/mole]
   }
  return theMass;
}

double NISTElementData::GetIsotopeMass(int Z, int N) {
   double theMass = 0.0;
   if (Z>0 && Z<=gNumberOfNISTElements) {
     int numisos = fNISTElementDataTable[Z-1].fNumOfIsotopes;
     int indxN   = N - fNISTElementDataTable[Z-1].fNIsos[0];
     if (indxN>=0 && indxN<numisos)
       theMass = fNISTElementDataTable[Z-1].fMassIsos[indxN];
   }
  return theMass;
}

// total electron binind energy in internal [energy] unit
double NISTElementData::GetBindingEnergy(int Z, int N) {
   double theBE = 0.0;
   if (Z>0 && Z<=gNumberOfNISTElements) {
     int numisos = fNISTElementDataTable[Z-1].fNumOfIsotopes;
     int indxN   = N - fNISTElementDataTable[Z-1].fNIsos[0];
     if (indxN>=0 && indxN<numisos)
       theBE = fBindingEnergies[Z-1];
   }
  return theBE;
}

void NISTElementData::PrintData(int Z) {
   using geant::g;
   using geant::mole;
   using geant::kAvogadro;
   using geant::perCent;
   using geant::GeV;

   std::cout<<"   *** NIST element data for "<<GetElementSymbol(Z)<<" Z = "<<Z<<" :"<<std::endl;
   std::cout<<"   Mean Atomic mass = "<<GetMeanAtomicMass(Z)/(g/mole) << " [g/mole]" <<std::endl;
   int numisos = GetNumberOfIsotopes(Z);
   std::cout<<"   Number of known isotopes = "<<GetNumberOfIsotopes(Z) <<std::endl;
   const std::string *symbols = GetIsotopeSymbols(Z);
   const int         *N       = GetIsotopeNucleonNums(Z);
   const double      *A       = GetIsotopeAtomicMasses(Z);
   const double      *W       = GetIsotopeNaturalAbundances(Z);
   const double      *M       = GetIsotopeMasses(Z);
   for ( int i = 0; i < numisos; ++i ) {
     std::cout<<"    "   <<std::setw(6)<<symbols[i]
              << "  N = "<<std::setw(4)<< N[i]
              << "  A = "<<std::setw(12)<<std::setprecision(8)<< A[i]*kAvogadro/(g/mole)      << " [g/mole]"
              << " natural abundance = " <<std::setw(12)<<std::setprecision(8)<< W[i]/perCent << " [%]"
              << " isotope mass = "     <<std::setw(12)<<std::setprecision(8)<< M[i]/GeV      << " [GeV]"
              <<std::endl;
   }
}

} // namespace geantphysics
