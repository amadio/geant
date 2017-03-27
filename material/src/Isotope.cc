
// for printouts
#include <iostream>
#include <iomanip>

//
#include "Isotope.h"
#include "NISTElementData.h"
#include "PhysicalConstants.h"

namespace geantphysics {

// static data member init
std::vector<Isotope*> Isotope::gTheIsotopeTable;  // global isotope table
std::map<int,int>     Isotope::gTheIsotopeMap;    // map to get iso. the index (used only internaly)


/**
 *  If a(atomic mass) and/or isomass(mass of the bare nucleus) are not given
 *  they will be taken from the internal NIST element database (NISTElementData).
 *  Similarly, if the name/symbol of the isotope is not given it will be set
 *  automatically.
 */
Isotope* Isotope::GetIsotope(int z, int n, double a, int isol, double isomass,
                             const std::string &name) {
   Isotope *theIsotope;
   // check if this isotope has already been created
   int isoindx = GetIsotopeIndex(z,n,isol);
   if (isoindx>-1) { // i.e. has already been created so return with the existing one
     theIsotope = Isotope::GetTheIsotopeTable()[isoindx];
   } else {             // i.e. has not been created yet so create it
     if (z<1) {
       std::cerr << "Wrong Isotope in Isotope::Isotope() " << name << " Z= "
                 << z << std::endl;
       exit(1);
     }

     if (n<z) {
       std::cerr << "Wrong Isotope in Isotope::Isotope()" << name << " Z= "
                 << z << " > N= " << n << std::endl;
       exit(1);
     }

     if (a<=0.) { // consistency of isotope mass with some particles guaranted only in this case !!!
       a = NISTElementData::Instance().GetAtomicMass(z,n);
       if (isomass<=0.) {
         isomass = NISTElementData::Instance().GetIsotopeMass(z,n);
       }
     } else if (isomass<=0.) {
       isomass =  a/geant::kAvogadro * geant::kCLightSquare
                - z*geant::kElectronMassC2
                + NISTElementData::Instance().GetBindingEnergy(z,n);
     }

     if (name!="") {
       theIsotope = new Isotope(name, z, n, a, isomass, isol);
     } else {
       theIsotope = new Isotope(NISTElementData::Instance().GetElementSymbol(z)
                                +std::to_string(n), z, n, a, isomass, isol);
     }
   }
  return theIsotope;
}

//
//ctr
//______________________________________________________________________________
Isotope::Isotope(const std::string &name, int z, int n, double a, double isomass, int isol)
  : fName(name), fZ(z), fN(n), fA(a), fIsoMass(isomass), fIsoL(isol) {
  fIndex = gTheIsotopeTable.size();
  gTheIsotopeTable.push_back(this);
  // add the index to the internal map
  gTheIsotopeMap[GetKey(z,n,isol)] = fIndex;
}

//
// dtr
//______________________________________________________________________________
Isotope::~Isotope(){
  const std::map<int,int>::iterator itr = gTheIsotopeMap.find(GetKey(fZ,fN,fIsoL));
  gTheIsotopeMap.erase(itr);
  gTheIsotopeTable[fIndex] = nullptr;
}

//______________________________________________________________________________
void Isotope::ClearAllIsotopes() {
  for (unsigned int i=0; i<gTheIsotopeTable.size(); ++i) {
    delete gTheIsotopeTable[i];
  }
  gTheIsotopeTable.clear();
}

//
// get index of the specified isotope in the global isotope table or -1 if the
// isotope has not been created yet
//______________________________________________________________________________
int Isotope::GetIsotopeIndex(int z, int n, int isol) {
  int indx = -1;
  int key = GetKey(z,n,isol);
  const std::map<int,int>::iterator itr = gTheIsotopeMap.find(key);
  if (itr != gTheIsotopeMap.end())
    indx = itr->second;
  return indx;
}


//
// Printouts
//
//______________________________________________________________________________
std::ostream& operator<<(std::ostream& flux, const Isotope* isotope){
  std::ios::fmtflags mode = flux.flags();
  flux.setf(std::ios::fixed,std::ios::floatfield);
  long prec = flux.precision(3);

  flux
    << " Isotope: " << std::setw(5) << isotope->fName
    << "   Z = " << std::setw(2)    << isotope->fZ
    << "   N = " << std::setw(3)    << isotope->fN
    << "   A = " << std::setw(6)    << std::setprecision(6)
    << (isotope->fA)/(geant::g/geant::mole) << " [g/mole]"
    << "   mass of nucleus = "  << isotope->fIsoMass/(geant::GeV) << " [GeV]";

  flux.precision(prec);
  flux.setf(mode,std::ios::floatfield);
  return flux;
}

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& flux, const Isotope& isotope)
{
  flux << &isotope;
  return flux;
}

//______________________________________________________________________________
std::ostream& operator<<(std::ostream& flux, std::vector<Isotope*> isotable)
{
 //dump info for all known isotopes
   flux
     << "\n***** Table : Nb of isotopes = " << isotable.size()
     << " *****\n" << std::endl;

   for (int i=0; i<isotable.size(); ++i)
     flux << isotable[i] << std::endl;

   return flux;
}

} // namespace geantphysics
