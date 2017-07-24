
#include "Material.h"

#include "Element.h"
#include "NISTElementData.h"
#include "NISTMaterialData.h"
#include "MaterialProperties.h"

#include <iostream>
#include <iomanip>
#include <cmath>

namespace geantphysics {
inline namespace GEANT_IMPL_NAMESPACE {

Vector_t<Material*>    Material::gTheMaterialTable;
Map_t<std::string,int> Material::gMapMaterialNameToIndex;

/**
 * If state of the material is MaterialState::kStateUndefined at construction the
 * state will be determined automatically based on the relation of the given
 * density to the kGasThreshold density constant and set internally .
 */
Material::Material(const std::string &name, double z, double a, double density, MaterialState state, double temp,
                   double pressure) : fName(name) {
  using geant::kUniverseMeanDensity;
  using geant::g;
  using geant::cm3;
  using geant::kGasThreshold;

  InitialiseMembers();

  if (density<kUniverseMeanDensity) {
    std::cerr << " Material WARNING in Material::Material() :"
              << " define a material with density=0 is not allowed. \n"
              << " The material " << name << " will be constructed with the"
              << " default minimal density: " << kUniverseMeanDensity/(g/cm3)
              << " [g/cm3]"
              << std::endl;
    density = kUniverseMeanDensity;
  }
  fIsUsed      = false;
  fDensity     = density;
  fState       = state;
  fTemperature = temp;
  fPressure    = pressure;
  // This material is made of one component and the corresponding element is
  // currently added. So set total number of components, current number of
  // components and number of elements to 1.
  fNumOfComponents  = fCurNumOfComponents = fNumOfElements = 1;
  // Create the corresponding element
  //
  // get the element name from the internal NIST database if we have data for this Z
  std::string enam, snam;
  int iz = std::lrint(z);
  if (iz<NISTElementData::Instance().GetNumberOfNISTElements()) {
    snam = NISTElementData::Instance().GetElementSymbol(iz);
    enam = snam;
  } else {
    enam = "ELM_" + name;
    snam = name;
  }
  // create the element and store its pointer in the element vector of this material
  fElementVector.push_back(new Element(enam, snam, z, a));
  // create and set mass fraction and realtive number of atoms per volume vectors
  fMassFractionVector     = new double[1];
  fMassFractionVector[0]  = 1.;
  fRelNumOfAtomsPerVol    = new double[1];
  fRelNumOfAtomsPerVol[0] = 1.;
  // set state (if state was not given)
  if (fState == MaterialState::kStateUndefined) {
    if (fDensity > kGasThreshold) {
      fState  = MaterialState::kStateSolid;
    } else {
      fState  = MaterialState::kStateGas;
    }
  }
  // create the associated MaterialProperties object
  fMaterialProperties = new MaterialProperties(this);
}


// Ctr to create a material from a combination of elements and/or materials
// subsequently added via AddElement and/or AddMaterial
Material::Material(const std::string &name, double density,int numcomponents, MaterialState state, double temp,
                   double pressure) : fName(name) {
  using geant::kUniverseMeanDensity;
  using geant::g;
  using geant::cm3;
  using geant::kGasThreshold;

  InitialiseMembers();

  if (density<kUniverseMeanDensity) {
    std::cerr << " Material WARNING in Material::Material() :"
	              << " define a material with density=0 is not allowed. \n"
	              << " The material " << name << " will be constructed with the"
	              << " default minimal density: " << kUniverseMeanDensity/(g/cm3)
	              << " [g/cm3]"
              << std::endl;
    density = kUniverseMeanDensity;
  }
  fIsUsed      = false;
  fDensity     = density;
  fState       = state;
  fTemperature = temp;
  fPressure    = pressure;
  // This material will be made of 'numcomponents' component.
  fNumOfComponents    = numcomponents;
  fCurNumOfComponents = 0; // zero components has been added so far
  fNumOfElements      = 0; // zero elements has been added so far
  // Prepare space for fCurArrayLenght = numcomponents components
  // These sizes are ok if we compose the material by adding elements by their
  // atom count i.e. AddElement(Element*,int) or by adding elements by their
  // mass fraction.
  // These sizes might need to be extended later if we compose the material by
  // adding (1)elements+materials, (2)materials+materials and at least one of the
  // materials is composed from more than one element.
  fCurArrayLenght = numcomponents;
  fElementVector.reserve(fCurArrayLenght);
  fMassFractionVector     = new double[fCurArrayLenght];
  fRelNumOfAtomsPerVol    = new double[fCurArrayLenght];
  // set material state if it was not defined
  // set state (if state was not given)
  if (fState==MaterialState::kStateUndefined) {
    if (fDensity>kGasThreshold) {
      fState = MaterialState::kStateSolid;
    } else {
      fState = MaterialState::kStateGas;
    }
  }
  fMaterialProperties = nullptr;
}


Material::~Material() {
  if (fMassFractionVector)    delete [] fMassFractionVector;
  if (fRelNumOfAtomsPerVol)   delete [] fRelNumOfAtomsPerVol;
  if (fMaterialProperties)    delete    fMaterialProperties;
  gTheMaterialTable[fIndex]  = nullptr;
}


void Material::ClearAllMaterials() {
  for (size_t i=0; i<gTheMaterialTable.size(); ++i) {
    delete gTheMaterialTable[i];
  }
  gTheMaterialTable.clear();
  // clear all elements as well
  Element::ClearAllElements();
}


// Add an element by atom count (not mixable with any fraction based Add.. )
void Material::AddElement(Element* element, int numatoms) {
  // filling ...
  if (fNumOfElements<fNumOfComponents) {
    fElementVector.push_back(element);
    // store number of atoms from this element in the fRelNumOfAtomsPerVol vect.
    fRelNumOfAtomsPerVol[fCurNumOfComponents] = numatoms;
    fCurNumOfComponents = ++fNumOfElements;
  } else {
    std::cerr << "Material::AddElement ERROR : Attempt to add more than the "
              << "declared number of elements for material: "<< std::endl
              << " Name = " << fName
              << " # components(elements) declared = " <<  fNumOfElements
              << std::endl;
    exit(-1);
  }
  // filled ..
  if (fCurNumOfComponents==fNumOfComponents) {
    double sumnumatoms = 0;
    double molarMass = 0.;
    for (int i=0; i<fNumOfElements; ++i) {
      // get mass of a mole (in internal [weight/mole] unit) of the i-th element
      // and multiply by the number of i-th atom in the molecule
      double w   = fRelNumOfAtomsPerVol[i]*fElementVector[i]->GetA();
      molarMass += w;
      fMassFractionVector[i] = w;
      // compute normalisation factor for realtive number of atoms
      sumnumatoms += fRelNumOfAtomsPerVol[i];
    }
    // normalise the masses to get the mass fraction
    // normalise the realtive number of atoms for each element
    for (int i=0; i<fNumOfElements; ++i) {
      fMassFractionVector[i]  /= molarMass;     // becomes ratio by weight
      fRelNumOfAtomsPerVol[i] /= sumnumatoms;   // becomes ratio by number of atoms
    }
    // create the associated MaterialProperties object
    if (fMaterialProperties) delete fMaterialProperties;
    fMaterialProperties = new MaterialProperties(this);
  }
}


// Add an element by its mass fraction (can be mixed with AddMaterial)
void Material::AddElement(Element *element, double massfraction) {
   // check if input data is correct
  if (massfraction<0.0 || massfraction>1.0) {
    std::cerr << "Material::AddElement ERROR : Attempt to add element with "
              << " wrong mass fraction in case of material: " << std::endl
              << " Material name = " << fName
              << " Element name  = " <<  element->GetName()
              << " Element mass fraction = " << massfraction
              << std::endl;
    exit(-1);
  }
  // filling ...
  if (fCurNumOfComponents<fNumOfComponents) {
    int el = 0;
    // check if the element has already been added to this material earlier
    // IF yes use that index ELSE use the current max elemnt index and insert element
    while ((el<fNumOfElements) && (element!=fElementVector[el])) ++el;
    // check if this element has already been found in the element vector of this material
    if (el<fNumOfElements) {  // already there
      fMassFractionVector[el] += massfraction;
    } else {                      // not there -> insert
      fElementVector.push_back(element);
      fMassFractionVector[el] = massfraction;
      ++fNumOfElements;
    }
    // one component has been added so increase the already added components counter
    ++fCurNumOfComponents;
  } else {
    std::cerr << "Material::AddElement ERROR : Attempt to add more than the "
              << "declared number of components for material: "<<std::endl
              << " Name = " << fName
              << " # components declared = " <<  fNumOfComponents
              << std::endl;
    exit(-1);
  }
  // filled.
  if (fCurNumOfComponents==fNumOfComponents) {
    // - check if sum of the mass fractions is unity ?
    // - prepare realtive number of atoms vector for the elements
    // - compute normalisation factor for relative number of atoms
    double normFactor = 0.;
    double sumWeight  = 0.;
    for (int i=0; i<fNumOfElements; ++i) {
      sumWeight += fMassFractionVector[i];
      fRelNumOfAtomsPerVol[i] = fMassFractionVector[i]/fElementVector[i]->GetA();
      normFactor += fRelNumOfAtomsPerVol[i];
    }
    if (std::fabs(1.-sumWeight)>geant::perThousand) {
      std::cerr << "WARNING !! Material::AddElement by mass fraction :"
                << std::endl <<" Sum of fractional masses = " << sumWeight
                << " is not = 1. This can lead to wrong results! "
                << std::endl << "  Material name = "<< fName << std::endl;
    }
    // compute relative number of atoms ( only a normalisation is left )
    normFactor = 1./normFactor;
    for (int i=0; i<fNumOfElements; ++i) {
      fRelNumOfAtomsPerVol[i] *= normFactor;
    }
    // create the associated MaterialProperties object
    if (fMaterialProperties) delete fMaterialProperties;
    fMaterialProperties = new MaterialProperties(this);
  }
}


// Add a Material by its mass fraction (can be mixed with AddElement by mass fraction)
void Material::AddMaterial(Material *material, double massfraction) {
   // check if input data is correct
  if (massfraction<0.0 || massfraction>1.0) {
    std::cerr << "Material::AddMaterial ERROR : Attempt to add material with "
              << " wrong mass fraction in case of material: " << std::endl
              << " Material name = " << fName
              << " Material name(try to add)  = " <<  material->GetName()
              << " Material mass fraction(try to add) = " << massfraction
              << std::endl;
    exit(-1);
  }
  // -get the element vector of the material that we try to add
  // -check how many new elements (that has not been added yet to this material)
  //  are among them
  // -if this is > 1 then we need to extend the space we allocated in the CTR
  //  for some arrays because we will add more than one elements at a component
  const Vector_t<Element*> elemVector = material->GetElementVector();
  int numNewElems = 0;
  for (auto elementptr : elemVector) {
    int el = 0;
    while ((el<fNumOfElements) && (elementptr!=fElementVector[el])) ++el;
    // if the element has not been found increase new element counter
    if (el==fNumOfElements)
      ++numNewElems;
  }
  if (numNewElems>1) {
    int nold         = fCurArrayLenght;
    fCurArrayLenght += numNewElems-1;
    // allocate new space for the mass fraction and rel. number of atoms vectors
    double *newMFR   = new double[fCurArrayLenght];
    double *newRNA   = new double[fCurArrayLenght];
    // copy old data
    for (int i=0; i<nold; ++i) {
      newMFR[i] = fMassFractionVector[i];
      newRNA[i] = fRelNumOfAtomsPerVol[i];
    }
    // delete old arrays
    delete [] fRelNumOfAtomsPerVol;
    delete [] fMassFractionVector;
    // set new pointers
    fMassFractionVector  = newMFR;
    fRelNumOfAtomsPerVol = newRNA;
  }
  // filling ...
  if (fCurNumOfComponents<fNumOfComponents) {
      // loop over the elements of the material that we try to add
      // check if the element has already been added to this material earlier
      // IF yes use that index ELSE use the current max elemnt index and insert element
      for (size_t ielem=0; ielem<elemVector.size(); ++ielem) {
        int el = 0;
        while ((el<fNumOfElements) && (elemVector[ielem]!=fElementVector[el])) ++el;
        // check if this element has already been found in the element vector of this material
        if (el<fNumOfElements) { // already there
          fMassFractionVector[el] += massfraction*(material->GetMassFractionVector())[ielem];
        } else {                     // not there -> insert
          fElementVector.push_back(elemVector[ielem]);
          fMassFractionVector[el] += massfraction*(material->GetMassFractionVector())[ielem];
          ++fNumOfElements;
        }
      }
      // one component has been added so increase the already added components counter
      ++fCurNumOfComponents;
  } else {
    std::cerr << "Material::AddMaterial ERROR : Attempt to add more than the "
              << "declared number of components for material: "<<std::endl
              << " Name = " << fName
              << " # components declared = " <<  fNumOfComponents
              << std::endl;
    exit(-1);
  }

  // filled. (same as in the case of AddElement by mass fraction )
  if (fCurNumOfComponents==fNumOfComponents) {
    // - check if sum of the mass fractions is unity ?
    // - prepare realtive number of atoms vector for the elements
    // - compute normalisation factor for relative number of atoms
    double normFactor = 0.;
    double sumWeight  = 0.;
    for (int i=0; i<fNumOfElements; ++i) {
      sumWeight += fMassFractionVector[i];
      fRelNumOfAtomsPerVol[i] = fMassFractionVector[i]/fElementVector[i]->GetA();
      normFactor += fRelNumOfAtomsPerVol[i];
    }
    if (std::fabs(1.-sumWeight)>geant::perThousand) {
      std::cerr << "WARNING !! Material::AddMaterial by mass fraction :"
                << std::endl <<" Sum of fractional masses = " << sumWeight
                << " is not = 1. This can lead to wrong results! "
                << std::endl << "  Material name = "<< fName << std::endl;
    }
    // compute relative number of atoms ( only a normalisation is left )
    normFactor = 1./normFactor;
    for (int i=0; i<fNumOfElements; ++i) {
      fRelNumOfAtomsPerVol[i] *= normFactor;
    }
    // create the associated MaterialProperties object
    if (fMaterialProperties) delete fMaterialProperties;
    fMaterialProperties = new MaterialProperties(this);
  }
}


void Material::InitialiseMembers() {
  //fElementVector.resize(0);
  fElementVector.clear();

  fMassFractionVector      = nullptr;
  fRelNumOfAtomsPerVol     = nullptr;
  fMaterialProperties      = nullptr;

  // initilized data members
  fIsUsed              = false;
  fDensity             = 0.0;
  fTemperature         = 0.0;
  fPressure            = 0.0;
  fNumOfComponents     = 0;
  fCurNumOfComponents  = 0;
  fNumOfElements       = 0;
  fCurArrayLenght      = 0;

  fState               = MaterialState::kStateUndefined;

  // Store in the static Table of Materials if there has not been created other
  // material before with this name. Error if the material name is not unique.
  if (FindMaterialIndex(fName)>-1) {
    std::cerr << "Material ERROR: duplicated name of material = "
              << fName
              << std::endl;
    exit(-1);
  }
  fIndex = gTheMaterialTable.size();
  // add this material to the global material table
  gTheMaterialTable.push_back(this);
  // add this index to the name -> index map
  gMapMaterialNameToIndex[fName] = fIndex;
}


int Material::FindMaterialIndex(const std::string &name) {
  int indx = -1;
  const Map_t<std::string,int>::iterator itr = gMapMaterialNameToIndex.find(name);
  if (itr!=gMapMaterialNameToIndex.end()) {
    indx = itr->second;
  }
  return indx;
}


Material* Material::GetMaterial(const std::string &name) {
  Material *mat = nullptr;
  int indx = FindMaterialIndex(name);
  if (indx>-1) {
    mat = gTheMaterialTable[indx];
  }
  return mat;
}


Material* Material::NISTMaterial(const std::string &name) {
  Material *mat = nullptr;
  // first check if material (NIST or other) has already been created with the
  // given name before and return with a ponter to that if yes.
  int indx = FindMaterialIndex(name);
  if (indx>-1) {
    mat = gTheMaterialTable[indx];
    return mat;
  }
  // Since no material was found (i.e. material with the given name has not been
  // created yet) we need to see if we can create the requested material as a NIST
  // material. So try to get the index of the NISTMaterialData structure.
  indx = NISTMaterialData::Instance().FindNISTMaterialDataIndex(name);
  if (indx>-1) {  // we have NISTMaterialData for this material name so buil the material
    // get NISTMaterialData from the corresponding data structure
    const std::string name = NISTMaterialData::Instance().GetName(indx);
    double        density         = NISTMaterialData::Instance().GetDensity(indx);
// will be set in the corresponding MaterialProperties
//    double        meanExcEnergy   = NISTMaterialData::Instance().GetMeanExcitationEnergy(indx);
    double        temperature     = NISTMaterialData::Instance().GetTemperature(indx);
    double        pressure        = NISTMaterialData::Instance().GetPressure(indx);
    MaterialState state           = NISTMaterialData::Instance().GetMaterialState(indx);
    int           numComp         = NISTMaterialData::Instance().GetNumberOfComponents(indx);
    const int    *listZ           = NISTMaterialData::Instance().GetListOfElements(indx);
    const double *listW           = NISTMaterialData::Instance().GetListOfElementFractions(indx);
    bool          isByAtomCount   = NISTMaterialData::Instance().IsToBuildByAtomCount(indx);

    // 1. create the Material
    mat = new Material(name, density, numComp, state, temperature, pressure);
    // 2. add numComp NIST elements to the material
    for (int i=0; i<numComp; ++i) {
      Element *elem = Element::NISTElement(listZ[i]);
      if (isByAtomCount) { // by number of atoms in the molecule
        mat->AddElement(elem, (int)(std::lrint(listW[i])));
      } else {               // by mass fraction
        mat->AddElement(elem, listW[i]);
      }
    }
  } else { // material was not found (has not been created before) and we don't have
           // NIST material data for the given material name ==> we cannot create the material
    std::cerr << "  ERROR in Material::GetNISTMaterial() : \n"
              << "   - material Name = "<<name << " has not been created earlier and\n"
              <<"      NISTMaterialData cannot be found for this material name."
              << std::endl;
    exit(-1);
  }
  return mat;
}


std::ostream& operator<<(std::ostream& flux, const Material* material) {
  std::ios::fmtflags mode = flux.flags();
  flux.setf(std::ios::fixed,std::ios::floatfield);
  long prec = flux.precision(3);

  using geant::g;
  using geant::cm3;
  using geant::atmosphere;
  using geant::kelvin;
  using geant::perCent;

  flux
    << " Material: "          << std::setw(8) <<  material->fName
    << "  with density = "    << std::setw(6) << std::setprecision(5)
    << material->fDensity / (g/cm3) << " [g/cm3]"
    << "\n";

  if (material->fMaterialProperties)
    flux << material->fMaterialProperties;
/*
  if (material->fMaterialProperties) {
    MaterialProperties *matp = material->fMaterialProperties;
    flux
    << std::setfill(' ') << std::setw(10) << " "
      << " Material properties:  "
      << " Imean =  " << matp->GetMeanExcitationEnergy()/geant::eV << " [eV]"
      << "\n";
  }
*/
  if (material->fState==MaterialState::kStateGas) {
    flux
      << std::setfill(' ') << std::setw(10) << " "
      << " Material state is kStateGas :"
      <<"  temperature = " << std::setw(6) << std::setprecision(2)
      << material->fTemperature / kelvin << " [K]"
      << "  pressure = "    << std::setw(6) << std::setprecision(2)
      << material->fPressure / atmosphere << " [atm]";
  }
  flux << "\n";
  for (int i=0; i<material->fNumOfElements; ++i) {
    flux
      << "\n   ---> " << (material->fElementVector)[i]
      << "\n          ElmMassFraction: "
      << std::setw(6)<< std::setprecision(2)
      << (material->fMassFractionVector)[i] / perCent << " [%]"
      << "  ElmAbundance "     << std::setw(6)<< std::setprecision(2)
      << 100*((material->fRelNumOfAtomsPerVol)[i])
      << " [%] \n";
  }
  flux.precision(prec);
  flux.setf(mode,std::ios::floatfield);
  return flux;
}


std::ostream& operator<<(std::ostream& flux, const Material& material) {
  flux << &material;
  return flux;
}


std::ostream& operator<<(std::ostream& flux, Vector_t<Material*> MaterialTable) {
  //Dump info for all known materials
  flux << "\n***** Table : Nb of materials = " << MaterialTable.size()
       << " *****\n\n";
  for (size_t i=0; i<MaterialTable.size(); ++i) {
    flux << MaterialTable[i] << std::endl << std::endl;
  }
  return flux;
}

} // inline namespace GEANT_IMPL_NAMESPACE
} // namespace geantphysics
