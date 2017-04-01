
#include "Element.h"

#include "Isotope.h"
#include "NISTElementData.h"
#include "PhysicalConstants.h"

#include <iostream>
#include <iomanip>
#include <cmath>

namespace geantphysics {

//static data member init
VectorHelper<Element*>::Vector_t Element::gTheElementTable;  // the global element table


/**
 * If element symbol is an empty string it will be set automatically based on the
 * provided atomic number.
 * \warning
 * If effective atomic mass is not specified it will be set automatically as the
 * natural abundance weighted average isotopic atomic mass. Therefore, the created
 * element is exactly the same as the corresponding NIST element i.e. if it were
 * constructed by the NISTElement() method except its name. So :
 *  - the corresponding NIST element will be created with a name staring with
 *    NIST_ELM_ if the corresponding NIST element has not been created yet.
 *  - the element as specified by the user will be created if the corresponding
 *    NIST element has already been created but a warning will be produced to
 *    suggest to obtain this element through the NISTElement() method.
 * \warning
 * In order to reduce the possibility of having duplicated elements we recommend
 * to use the NISTElement() method to obtain the element if the atomic mass is not
 * specified (instead of using this constructor directly). Using the NISTElement()
 *  method guarantees uniqueness of NIST elements.
 *
 * \warning
 * When the isotopes of the element are not specified, the natural isotopes (i.e.
 * isotopes with non-zero natural isotope abundance) will be taken from NISTElementData
 * isotopic composition database. In case of elements, that does not have any
 * stable isotopes, the most stable isotope will be inserted with abundance = 1.0 and
 * the user will be informed by a warning. It is always possible to build an element
 * by specifying its isotopes with their abundances. In the case of elements, that
 * has no stable isotope we recommend to use this opportunity to build the lement.
 */
//
//ctr for direct build
Element::Element(const std::string &name, const std::string &symbol, double  zeff, double  aeff)
  : fName(name), fSymbol(symbol){
     using geant::perMillion;

     int iz = std::lrint(zeff);
     if (iz<1) {
       std::cerr << "Fail to create Element in Element::Element():" << name
                 << " Z= " << zeff << " < 1 !";
       exit(-1);
     }
     if (std::abs(zeff - iz)>perMillion) {
       std::cerr << "Element Warning in Element::Element():  " << name
                 << " Z= " << zeff
                 << std::endl;
     }
     //
     InitialiseMembers();
     fZeff   = zeff;
     AddNaturalIsotopes();
     fAeff   = aeff;

     // if effective atomic mass is not given then:
     //  - effective atomic mass will be taken from the NIST database i.e. natural
     //    abundance weighted mean isotopic atomic mass
     //  - therefore, this element is exactly the same as the corresponding NISTElement
     //    with the same Z(atomic number) except its name.
     //  - So if the corresponding NIST element has already been created:  inform
     //    the user that this element is a duplicate of that NIST element we recommend
     //    to use that one in order to avoid duplication.
     //  - If the corresponding NIST element has NOT been created yet: create the
     //    corresponding NIST element and inform the user that NIST element has been
     //    created and the name of the element has been changed.
     if (aeff<=0.0) {
       // NIST element data must be used so check if NISTElementData has data
       // for this z and error if not.
       int maxZet = NISTElementData::Instance().GetNumberOfNISTElements();
       if (iz>maxZet) {
         std::cerr << "Failed to create NIST element with Z= " << iz
                   << " - available atomic numbers in the NISTElementData 1 <= Z <= " << maxZet
                   << std::endl;
         exit(-1);
       }
       // check if the corresponding NIST element has already been created
       int indx = NISTElementData::Instance().GetNISTElementIndex(iz);
       if (indx>-1) { // has already been created
       std::cerr << "Element Warning : in case of Element :"
                 << "  Name = "<< fName << " with Z = " << zeff << " : " << std::endl
                 << "  Since natural isotope composition is taken and effective atomic mass is \n"
                 << "  not given this Element is exactly the same as the already existing NIST \n"
                 << "  element with the same atomic number (Z). We recommend to use that NIST  \n"
                 << "  element through the Element::NISTElement(Z) method instead of creating  \n"
                 << "  this element as an exact duplicate (except its name)."
                 << std::endl;
         fAeff = NISTElementData::Instance().GetMeanAtomicMass(iz);
       } else {       // has NOT been created yet so create it and inform the user
         // create the corresponding NIST element
         fName   = "NIST_ELM_"+NISTElementData::Instance().GetElementSymbol(iz);
         fSymbol = NISTElementData::Instance().GetElementSymbol(iz);
         fAeff   = NISTElementData::Instance().GetMeanAtomicMass(iz);
         // inform the user
         std::cerr << "Element Warning : in case of Element :"
                   << "  Name = "<< name << " with Z = " << zeff << " : " << std::endl
                   << "  Since natural isotope composition is taken and effective atomic mass is    \n"
                   << "  not given, this Element is exactly the same as if it were created as a NIST \n"
                   << "  element by means of Element::NISTElement(Z) except its name. In order to   \n"
                   << "  to avoide duplication, the corresponding NIST element has been created with\n"
                   << "  the Element Name = "<< fName
                   << std::endl;
       }
     }

/*     {

       fAeff = 0.0;
       for (int i = 0; i < fCurNumOfIsotopes; ++i)
         fAeff +=  fRelativeIsotopeAbundanceVector[i]*(fIsotopeVector[i])->GetA();
     }
*/

     fIndex = gTheElementTable.size();
     gTheElementTable.push_back(this);
}


//
// ctr to build an element from isotopes by adding the isotopes one by one via
// AddIsotope method
Element::Element(const std::string &name, const std::string &symbol, int nisotopes)
  : fName(name),fSymbol(symbol) {
  InitialiseMembers();
  int n = nisotopes;
  if (0>=nisotopes) {
    std::cerr << "Fail to create Element in Element::Element()" << name
              << " <" << symbol << "> with " << nisotopes
              << " isotopes";
    exit(-1);
  } else {
    fIsotopeVector.resize(n,0);
    fRelativeIsotopeAbundanceVector = new double[nisotopes];
  }
}


//
// Add an isotope to the element
void Element::AddIsotope(Isotope* isotope, double abundance) {
  if (fIsotopeVector.size()==0) {
    std::cerr << "Fail to add Isotope to Element in Element::AddIsotope()"
              << fName << " because number of isotopes was set to <= 0 !";
    exit(-1);
  }
  int iz = isotope->GetZ();
  // filling ...
  if (fCurNumOfIsotopes<(int)fIsotopeVector.size()) {
    // check same Z
    if (fCurNumOfIsotopes==0) {
      fZeff = double(iz);
    } else if (double(iz)!=fZeff) {
      std::cerr << "Fail to add Isotope Z= " << iz << " to Element " << fName
                << " with different Z= " << fZeff /*<< fNeff */<< " in G4Element::AddIsotope()";
      exit(-1);
    }
    //Z ok
    fRelativeIsotopeAbundanceVector[fCurNumOfIsotopes] = abundance;
    fIsotopeVector[fCurNumOfIsotopes]                  = isotope;
    ++fCurNumOfIsotopes;
  } else {
    std::cerr << "Fail to add Isotope Z= " << iz << " to Element " << fName
             << " - more isotopes than declaired ";
    exit(-1);
  }
  // filled.
  if (fCurNumOfIsotopes==(int)fIsotopeVector.size()) {
    double wtSum=0.0;
    fAeff = 0.0;
    for (int i=0; i<fCurNumOfIsotopes; ++i) {
      fAeff +=  fRelativeIsotopeAbundanceVector[i]*(fIsotopeVector[i])->GetA();
      wtSum +=  fRelativeIsotopeAbundanceVector[i];
    }
    if (wtSum>0.0) {
      fAeff /= wtSum;
    }
    // renormalise abundances
    if (wtSum!=1.0) {
      wtSum = 1./wtSum;
      for (int i=0; i<fCurNumOfIsotopes; ++i) {
	      fRelativeIsotopeAbundanceVector[i] *= wtSum;
      }
    }
    // Set element symbol if it was not given
    if (""==fSymbol) {
      fSymbol = NISTElementData::Instance().GetElementSymbol(iz);
    }
    // set the index of this element and add to the global element table
    fIndex = gTheElementTable.size();
    gTheElementTable.push_back(this);
  }
}


//
// member initialiser
void Element::InitialiseMembers() {
  fZeff      = 0.;
  fAeff      = 0.;

  fCurNumOfIsotopes               = 0;
  fRelativeIsotopeAbundanceVector = nullptr;
  fIsotopeVector.clear();

  fIndex               = -1;
  fIsNaturalAbundance  = false;
}


//
// dtr
Element::~Element(){
  if (fRelativeIsotopeAbundanceVector) {
    delete [] fRelativeIsotopeAbundanceVector;
  }
  //remove this element from theElementTable
  gTheElementTable[fIndex] = nullptr;
}


void Element::ClearAllElements() {
  for (size_t i=0; i<gTheElementTable.size(); ++i) {
    delete gTheElementTable[i];
  }
  gTheElementTable.clear();
  // clear all isotopes as well
  Isotope::ClearAllIsotopes();
}


//
// Get (find or build) a NIST element i.e. based on NIST natural isotope composition.
Element* Element::NISTElement(double zeff) {
  int iz = std::lrint(zeff);
  //check if NISTElementData has data for this z and error if not.
  int maxZet = NISTElementData::Instance().GetNumberOfNISTElements();
  if (iz<1 || iz>maxZet) {
    std::cerr << "Failed to create NIST element with Z= " << iz
              << " - available atomic number in the NISTElementData 1 <= Z <= " << maxZet
              << std::endl;
    exit(-1);
  }
  Element *elem = nullptr;
  // check if NIST element with this atomic number has already been created
  int indx = NISTElementData::Instance().GetNISTElementIndex(iz);
  if (indx>-1) {   // already built so just return with its pointer from the global element table
    elem = gTheElementTable[indx];
  } else {         // has not been built yet so build it now and set its index
    double effAmass = NISTElementData::Instance().GetMeanAtomicMass(iz);
    elem = new Element("NIST_ELM_"+NISTElementData::Instance().GetElementSymbol(iz),"",(double)iz,effAmass);
    //set its index
    NISTElementData::Instance().SetNISTElementIndex(iz, elem->GetIndex());
  }
  return elem;
}


//
// add natuarl isotopes (taken from NIST database)
void Element::AddNaturalIsotopes() {
  int Z = std::lrint(fZeff);
  // check if we have data for this Z in the NIST element data table
  if (Z>NISTElementData::Instance().GetNumberOfNISTElements())
    return;

  // Get number of isotopes
  int numisos = NISTElementData::Instance().GetNumberOfIsotopes(Z);
  // Set element symbol if it was not given
  if (""==fSymbol)
    fSymbol = NISTElementData::Instance().GetElementSymbol(Z);

  // count number of natural isotopes i.e. with natural abundance > 0.
  const double *natabs = NISTElementData::Instance().GetIsotopeNaturalAbundances(Z);
  const int    *arrayN = NISTElementData::Instance().GetIsotopeNucleonNums(Z);
  fCurNumOfIsotopes = 0;
  for (int i=0; i<numisos; ++i) {
    if (natabs[i]>0.0) {
      ++fCurNumOfIsotopes;
    }
  }
  int     idx = 0;
  double xsum = 1.0;
  // if the element has no natural isotope then we take the most stable one but we give a warning
  if (fCurNumOfIsotopes==0) {
    fCurNumOfIsotopes = 1;
    idx = 1;
    // set container size
    fIsotopeVector.resize(fCurNumOfIsotopes, 0);
    fRelativeIsotopeAbundanceVector = new double[fCurNumOfIsotopes];
    // get the index of the most stable isotope for this Z
    int indxN = NISTElementData::Instance().GetIndexOfTheMostStableIsotope(Z);
    // get or create the isotope
    fIsotopeVector[0] = Isotope::GetIsotope(Z, arrayN[indxN]);
    fRelativeIsotopeAbundanceVector[0] = 1.0;
    std::cerr << "Element Warning : in case of AddNaturalIsotopes :\n"
              << "Element Z = "<< Z << " has no stable isotope i.e. all natural abundances are zero.\n"
              << "The most stable isotope Z = " << Z << "  N = "<< arrayN[indxN]<< " has been used with \n"
              << "\"natural\" abudance = 1.0 in order to avoid having an element without any isotopes.\n"
              << "Note, that elements can be built up by specifying their isotopes. We strongly recommend\n"
              << "to use this later possibility in case of element Z = "<<Z<<" !"<<std::endl;

  } else {
    // set container size
    fIsotopeVector.resize(fCurNumOfIsotopes, 0);
    fRelativeIsotopeAbundanceVector = new double[fCurNumOfIsotopes];
    idx = 0;
    xsum = 0.0;
    for (int i=0; i<numisos; ++i) {
      if (natabs[i]>0.0) {
        // get or create the isotope
        fIsotopeVector[idx] = Isotope::GetIsotope(Z, arrayN[i]);
        fRelativeIsotopeAbundanceVector[idx] = natabs[i];
        xsum += natabs[i];
        ++idx;
      }
    }
  }

  // make sure that relative isotope abundances are properly normalized
  if (xsum!=0.0 && xsum!=1.0) {
    xsum = 1./xsum;
    for (int i=0; i<idx; ++i) {
      fRelativeIsotopeAbundanceVector[i] *= xsum;
    }
  }
  // the the flag that natural isotope abundances are used for this element
  fIsNaturalAbundance = true;
}


//
// Printouts
std::ostream& operator<<(std::ostream& flux, const Element* element) {
  using geant::g;
  using geant::mole;
  using geant::perCent;

  std::ios::fmtflags mode = flux.flags();
  flux.setf(std::ios::fixed,std::ios::floatfield);
  long prec = flux.precision(3);

  flux
    << " Element: " << element->fName   << " (" << element->fSymbol << ")"
    << "   Z = "    << std::setw(4) << std::setprecision(1) <<  element->fZeff
//    << "   N = "    << std::setw(4) << std::setprecision(1) <<  std::lrint(element->fNeff)
    << "   A = "    << std::setw(4) << std::setprecision(3)
                    << (element->fAeff)/(g/mole) << " [g/mole]";

  for (int i=0; i<element->fCurNumOfIsotopes; ++i) {
    flux
      << "\n         ---> " << (element->fIsotopeVector)[i]
      << "   abundance: " << std::setw(6) << std::setprecision(3)
      << (element->fRelativeIsotopeAbundanceVector[i])/perCent << " [%]";
  }
  flux.precision(prec);
  flux.setf(mode,std::ios::floatfield);
  return flux;
}


 std::ostream& operator<<(std::ostream& flux, const Element& element) {
  flux << &element;
  return flux;
}


std::ostream& operator<<(std::ostream& flux, VectorHelper<Element*>::Vector_t elementtable) {
 //Dump info for all known elements
   flux << "\n***** Table : Nb of elements = " << elementtable.size()
        << " *****\n" << std::endl;
   for (size_t i=0; i<elementtable.size(); i++) {
     flux << elementtable[i] << std::endl << std::endl;
   }
   return flux;
}

} // namespace geantphysics
