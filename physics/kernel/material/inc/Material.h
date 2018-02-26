

#ifndef MATERIAL_H
#define MATERIAL_H

#include "Types.h"

#include "Geant/PhysicalConstants.h"
#include "MaterialState.h"
#include "Element.h"

#include <string>

GEANT_DEVICE_DECLARE_CONV(geantphysics, class, Material);

namespace geantphysics {
inline namespace GEANT_IMPL_NAMESPACE {

class MaterialProperties;
class Element;

/**
 * @brief   Class to describe a material.
 * @class   Material
 * @author  M Novak, A Ribon
 * @date    december 2015
 *
 * Name of the materials must be unique i.e. there cannot be more than one
 * material with the same material name.\n
 * Materials can be created in different ways:
 * 1. Simple materials, that consists of only one element can be created directly.
 * 2. Complex materials, that consists of more than one elements can be build up by:
 *   -  First creating a material object by declaring the number of components the
 *      material will be built up.
 *      1. Then adding elements to the material one by one together with the
 *         (integer) number of atoms in the molecule using the AddElement() method.
 *      2. Or adding elements and/or materials to the material one by one together
 *         with their mass fraction using the AddElement() and/or AddMaterial()
 *         methods.
 * 3. Several pre-defined (NIST) materials can be constructed by their names. See
 *    the NISTMaterial() mathod for more details.
 *
 * Each material stores a list of pointers to the elements from which the material
 * is built up. This can be accessed by the user through the GetElementVector() method.
 *
 * Each material pointer is stored in a global material table when created.
 * This table can be accessed by the user through the GetTheMaterialTable() method.
 *
 * \warning It's strongly not recommended to delete any materials! Each material
 * will be deleted automatically at the end of the run.
 *
 */


class Material {
 public:
/**
* @name Constructors, destructor and construction realted public methods:
*/
//@{
   /**
    * @brief Constructor to create simple material from single element
    *
    * @param[in] name Name of the material.
    * @param[in] z Atomic number of the single element that this material is built up.
    * @param[in] a Atomic mass of the single element that this material is built up
    *              in internal [weight/amount of substance] units.
    * @param[in] density Density of the material in internal [weight/length\f$^3\f$] units.
    * @param[in] state State of the material (see MaterialState).
    * @param[in] temp  Temperature of the material in internal [temperature] units.
    * @param[in] pressure  Pressure of the material in internal [pressure] units.
    */
    //
    // Ctr to create a material from single element
   Material(const std::string &name, double z, double a, double density,
            MaterialState state  = MaterialState::kStateUndefined,
            double temp          = geant::units::kNTPTemperature,
            double pressure      = geant::units::kSTPPressure);

   /**
    * @brief Constructor to create material from elements and/or materials.
    *
    * Materials can be created as combination of elements and/or materials added
    * subsequently by using the AddElement() and/or AddMaterial() methods.
    *
    * @param[in] name Name of the material.
    * @param[in] density Density of the material in internal [weight/length\f$^3\f$] units.
    * @param[in] numcomponents Number of components(elements and/or materials) this material will be bult up.
    * @param[in] state State of the material (see MaterialState).
    * @param[in] temp  Temperature of the material in internal [temperature] units.
    * @param[in] pressure  Pressure of the material in internal [pressure] units.
    */
    // Ctr to create a material from a combination of elements and/or materials
    // subsequently added via AddElement and/or AddMaterial
   Material(const std::string &name, double density,int numcomponents,
             MaterialState state  = MaterialState::kStateUndefined,
             double temp          = geant::units::kNTPTemperature,
             double pressure      = geant::units::kSTPPressure);

   /**
    * @brief Method to add an element to the material based on molecular composition.
    *
    * @param[in] element Pointer to the currently added element object.
    * @param[in] numatoms Number of the currently added atoms in the molecule.
    */
   // Way 1 of building up complex material with the above CTR
   // Add an element, giving the number of this atoms in the molecule
   void AddElement(Element *element, int numatoms);

   // Way 2 of building up complex material with the above CTR
   // Add an element or material, giving fraction of mass
   /**
    * @brief Method to add an element to the material based on its mass fraction.
    *
    * @param[in] element Pointer to the currently added Element object.
    * @param[in] massfraction Mass fraction of the currently added element in the material.
    */
   void AddElement (Element  *element , double   massfraction);
   /**
    * @brief Method to add a material to the material based on its mass fraction.
    *
    * @param[in] material Pointer to the currently added Material object.
    * @param[in] massfraction Mass fraction of the currently added Material in the material.
    */
   void AddMaterial(Material *material, double   massfraction);

   /** @brief Destructor.
    *
    *  \warning It's strongly not recommended to delete any materials! Each
    *  material will be deleted automatically at the end of the run.
    */
  ~Material();
//@}


/**
 * @name Field related getters:
 */
//@{
  /** @brief Public method to get the name of this material.
   *  @return Name of this material. */
  const std::string& GetName() const {return fName;}

  /** @brief Public method to get the flag that indicates whether the material is used or not in the geometry.
   *  @return  that indicates whether the material is used or not in the geometry */
  bool   IsUsed() const {return fIsUsed;}

  /** @brief Public method to set the flag that indicates whether the material is used or not in the geometry.*/
  void   SetIsUsed(bool val) {fIsUsed = val;}

  /** @brief Public method to get the density of this material.
   *  @return Density of this material in internal [weight/ length/f$^3/f$] units. */
  double  GetDensity() const {return fDensity;}

  /** @brief Public method to get the temperature of this material.
   *  @return Temperature of this material in internal [temperature] units. */
  double  GetTemperature() const {return fTemperature;}

  /** @brief Public method to get the pressure of this material.
   *  @return Pressure of this material in internal [pressure] units. */
  double  GetPressure() const {return fPressure;}

  /** @brief Public method to get the index of this material in the global material table.
   *  @return index of this material in the global material table. */
  int     GetIndex() const {return fIndex;}

  /** @brief Public method to get the state of this material.
   *  @return State of this material. */
  MaterialState  GetMaterialState() const {return fState;}

  /** @brief Public method to access the material properties object associated
   *         to this material.
   *  @return Material properties object of this material. */
  MaterialProperties* GetMaterialProperties() const {return fMaterialProperties;}
//@}



/**
 * @name Methods to access data regarding the elements constituting this Material:
 */
//@{
  /** @brief Public method to get the number of elements this material is built up.
   *  @return Number of elements this material is bult up. */
  int     GetNumberOfElements() const {return fNumOfElements;}

  /** @brief Public method to get the list of elements this material is built up.
   *  @return Vector of elements this material is bult up. */
  const Vector_t<Element*>& GetElementVector() const {return fElementVector;}

  /** @brief Public method to get the mas fraction for each elements this material is built up.
   *  @return Mass fraction for each elements this material is built up [GetNumberOfElements()]. */
  const double*  GetMassFractionVector() const {return fMassFractionVector;}

  /** @brief Public method to get the relative number of atoms for each elements this material is built up.
   *  @return Relative number of atoms for each elements this material is built up [GetNumberOfElements()]. */
  const double*  GetRelativeNumberOfAtomsPerVolume() const {return fRelNumOfAtomsPerVol;}
//@}


/**
 * @name NIST material related public mathod:
 */
//@{
  /**
   * @brief Public method to retrieve/create a (NIST) material based on material
   *        name.
   *
   * Several pre-defined (NIST) materials can be constructed through this method
   * by specifying the requested material with its name. The NISTMaterialData
   * internal database stores all the data that are necessary to build the
   * corresponding material. These pre-defined materials are identified by their
   * names that always starts with the NIST_MAT_ prefix. See \ref NISTMaterialDataDoc
   * or the implementation of NISTMaterialData::BuildTable() method for
   * the list of available pre-defined material names.
   *
   * No new material is created if a material(NIST or other) with the given name
   * has already been constructed earlier i.e. can be found in the global material
   * table. A pointer to the already existing material will be returned in this
   * case.\n
   * If material, with the given name, has not been created before:
   *  - a NIST material will be created if the requested material can be found
   *    in the internal NISTMaterialData database (by a perfect match between
   *    the given name of the material and one of those material data stored in
   *    the internal NISTMaterialData database).
   *  - error otherwise
   *
   * @param[in] name Name of the requested NIST material.
   * @return    Pointer to the requested (NIST or other) material if:
   *              - a material has already been created earlier with the given name.
   *              - the given material name corresponds to a NIST material name in
   *                the internal NISTMaterialData internal database.
   *
   */
   static Material* NISTMaterial(const std::string &name);
//@}

/**
 * @name Methods to get data regarding the global material table:
 */
//@{
   /**
    * @brief Public method to get the global material table.
    * @return Vector of pointers to all the materials that has been created so far.
    */
  static const Vector_t<Material*>& GetTheMaterialTable() { return gTheMaterialTable;}
  /**
   * @brief Public method to delete all Material objects that has been created.
   *
   * This method is called by the PhysicsProcessHandler when it is deleted to clean all Material objects.
   * Users should never call this method!
   */
  static void ClearAllMaterials();
  /**
   * @brief Public method to get the number of materials in the global material table.
   * @return Number of materials in the global material table i.e. number of materials
   *         have been created so far.
   */
  static int GetNumberOMaterias() {return gTheMaterialTable.size();}

  /**
   * @brief Public method to retrieve a material from the global material table
   *        by material name.
   *
   * @param[in] name Name of the requested material.
   * @return    - Pointer to the requested material if a material has already been
   *              created earlier with the given name.
   *            - nullptr otherwise.
   */
   static Material* GetMaterial(const std::string &name);
//@}

/**
 * @name Printouts:
 */
//@{
  friend std::ostream& operator<<(std::ostream&, const Material*);
  friend std::ostream& operator<<(std::ostream&, const Material&);
  friend std::ostream& operator<<(std::ostream&, Vector_t<Material*>);
//@}


// THESE ARE ONLY FOR THE TABULATED PHYSICS AND WILL BE REMOVED SOON
  void  SetXsecPtr(void* xsec) { fTabXsecPtr = xsec; }
  void* GetXsecPtr() const { return fTabXsecPtr; }
//

 private:
   /** @brief Copy constructor. Not implemented yet */
   Material(const Element&);                // not implemneted

   /** @brief Equal operator. Not implemented yet */
   Material& operator=(const Material&);    // not implemented

   /** @brief Helper method to initialise some mebers */
   void InitialiseMembers();

   /** @brief Helper method to find the index of a material given by its name.
    *  @param[in]  Name of the material.
    * @return      - Index of the material in the global material table if there
    *                is a material created earlier with the given name.
    *              - -1 otherwise.
    */
   static int FindMaterialIndex(const std::string &name);

 //
 // data members
 private:
   /** @brief Name of the material. */
   std::string     fName;                 // material name
   /** @brief Does the material used in the current geometry? */
   bool            fIsUsed;
   /** @brief Density of the material in internal [weight/length/f$ ^3 /f$] units. */
   double          fDensity;              // material density
   /** @brief Temperature of the material in internal [temperature] units. */
   double          fTemperature;          // temperature (defaults: NTP)
   /** @brief Pressure in internal [pressure] units. */
   double          fPressure;             // pressure    (defaults: STP)
   /** @brief Number of components this material is built up (internal only). */
   int             fNumOfComponents;      // total number of components in the material (used only internally, go getter)
   /** @brief Number of components that has been added to so far (internal only). */
   int             fCurNumOfComponents;   // number of components added so far(used only internally, no getter)
   /** @brief Number of elements this material is built up. */
   int             fNumOfElements;        // number of elements this material is composed from
   /** @brief Index of this material in the global material table. */
   int             fIndex;                // the position of this material in the global material table
   /** @brief Current length of some arrays at construction (internal only). */
   int             fCurArrayLenght;       // just a helper at ctr: lenght of currently pre-allocated space(no getter)
   /** @brief Sate of the material. */
   MaterialState   fState;                // state of the material (default kStateUndefined)

   /** @brief Ratio by mass for each element this material is bult up. */
   double       *fMassFractionVector;     // ratio by mass for each element this material is built up
   /** @brief Relative number of atoms per volume for each element this material is buit up. */
   double       *fRelNumOfAtomsPerVol;    // relative number of atoms per volume for each element this
                                          // material is built up
   /** @brief List of elements this material built up. */
   Vector_t<Element*> fElementVector;  // vector of element pointers this material is built up

   /** @brief The global material table. */
   static Vector_t<Material*> gTheMaterialTable; // the global material table

   /** @brief Internal map to store the already created material indices in the global material
    *         table with a key = material name. */
   static Map_t<std::string,int> gMapMaterialNameToIndex;

   /** @brief Object to store additional properties realted to this material (the class owns the object)*/
   MaterialProperties *fMaterialProperties;

   /** @brief THIS MEMBER IS ONLY FOR THE TABULATED PHSYICS AND WILL BE REMOVED SOON */
   void     *fTabXsecPtr;
};

} // GEANT_IMPL_NAMESPACE
} // namespace geantphysics

#endif // MATERIAL_H
