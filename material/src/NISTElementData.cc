#include "NISTElementData.h"
#include "PhysicalConstants.h"

#include <iostream>
#include <iomanip>

namespace geantphysics {

NISTElementData& NISTElementData::Instance() {
   static NISTElementData instance;
   return instance;
}

NISTElementData::NISTElementData() {
   BuildTable();
   // init list of NIST element indices(indices in the global element table) that has
   // been already built to -1
   for (int i=0; i<gNumberOfNISTElements; ++i)
     fIndicesOfBuiltNISTElements[i] = -1;
}

// we compute from the isotope mass to conserve consistency with some particle masses
double NISTElementData::GetAtomicMass(int z, int n) {
   using geant::kElectronMassC2;
   using geant::kAvogadro;
   using geant::kCLightSquare;
   constexpr double unitconv = kAvogadro/kCLightSquare;

   double theMass = -1.0;
   if (z>0 && z<=gNumberOfNISTElements) {
     int numisos = fNISTElementDataTable[z-1].fNumOfIsotopes;
     int indxN   = n - fNISTElementDataTable[z-1].fNIsos[0];
     if (indxN>=0 && indxN<numisos)
       theMass = fNISTElementDataTable[z-1].fMassIsos[indxN] + z*kElectronMassC2
                 - fBindingEnergies[z-1];
       theMass *= unitconv; // convert energy to [weight/mole]
   }
   if (theMass<0.0) {
     std::cerr << " *** ERROR NISTElementData::GetAtomicMass \n "
               << "   unknown isotope: atomic number = " << z << " nucleon number = " << n
               << std::endl;
     exit(1);
   }
  return theMass;
}

double NISTElementData::GetIsotopeMass(int z, int n) {
   double theMass = -1.0;
   if (z>0 && z<=gNumberOfNISTElements) {
     int numisos = fNISTElementDataTable[z-1].fNumOfIsotopes;
     int indxN   = n - fNISTElementDataTable[z-1].fNIsos[0];
     if (indxN>=0 && indxN<numisos)
       theMass = fNISTElementDataTable[z-1].fMassIsos[indxN];
   }
   if (theMass<0.0) {
     std::cerr << " *** ERROR NISTElementData::GetIsotopeMass \n "
               << "   unknown isotope: atomic number = " << z << " nucleon number = " << n
               << std::endl;
     exit(1);
   }
  return theMass;
}

// total electron binind energy in internal [energy] unit
double NISTElementData::GetBindingEnergy(int z, int n) {
   double theBE = -1.0;
   if (z>0 && z<=gNumberOfNISTElements) {
     int numisos = fNISTElementDataTable[z-1].fNumOfIsotopes;
     int indxN   = n - fNISTElementDataTable[z-1].fNIsos[0];
     if (indxN>=0 && indxN<numisos)
       theBE = fBindingEnergies[z-1];
   }
   if (theBE<0.0) {
     std::cerr << " *** ERROR NISTElementData::GetBindingEnergy \n "
               << "   unknown isotope: atomic number = " << z << " nucleon number = " << n
               << std::endl;
     exit(1);
   }
  return theBE;
}

void NISTElementData::PrintData(int z) {
   using geant::g;
   using geant::mole;
   using geant::kAvogadro;
   using geant::perCent;
   using geant::GeV;

   std::cout<<"   *** NIST element data for "<<GetElementSymbol(z)<<" Z = "<<z<<" :"<<std::endl;
   std::cout<<"   Mean Atomic mass = "<<GetMeanAtomicMass(z)/(g/mole) << " [g/mole]" <<std::endl;
   int numisos = GetNumberOfIsotopes(z);
   std::cout<<"   Number of known isotopes = "<<GetNumberOfIsotopes(z) <<std::endl;
   const std::string *symbols = GetIsotopeSymbols(z);
   const int         *N       = GetIsotopeNucleonNums(z);
   const double      *A       = GetIsotopeAtomicMasses(z);
   const double      *W       = GetIsotopeNaturalAbundances(z);
   const double      *M       = GetIsotopeMasses(z);
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
