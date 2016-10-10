#include <iostream>

#ifndef NISTELEMENTDATA_H
#define NISTELEMENTDATA_H

#include <string>
#include <vector>

namespace geantphysics {
  /**
   * @brief   Internal(singletone) database to store NIST element data.
   * @class   NISTElementData
   * @author  M Novak, A Ribon
   * @date    december 2015
   *
   * Internal database to store Atomic Weights and Isotopic Compositions for All
   * Data are taken from NIST database \cite nist2015isotopcomp (data were taken
   * 10 December 2015). Data used in the database can be found at \ref NISTElementDataDoc.
   * Total binding energy of the electrons are computed based on the parametrisation
   * given by Lunney et. al. \cite lunney2003recent. Isotope mass is the mass of
   * the bare nucleus. Some of the isotope atomic masses have been changed such
   * that the corresponding isotope masses are equal to the rest mass of Deuteron,
   * Triton, He3 and Alpha in the internal database  (for full consistency).
   *
   * Isotopic composition of the elements are used to add the natural isotopes
   * to the elements if the isotopic composition is not given by the user. In case
   * of elements having no stable (natural) isotope, the most stable isotope will
   * be used (and the user will be informed by a warning). This information is
   * also stored in this databes.
   *
   * \todo In the private NISTElementData::BuildTable() we need to make sure
   * that we use masses for Deuteron, Triton, He3 and Alpha that consistent with
   * the masses of these particles in their definition.
   */

class NISTElementData{
 public:
   /** @brief Public method to access the singletone database instance. */
   static NISTElementData& Instance();

/**
 * @name Public getters:
 */
//@{
   /** @brief Public method to access the number of elements in the database. */
   int                GetNumberOfNISTElements()           const { return gNumberOfNISTElements; }

   /** @brief Public method to access element symbols.
    *
    * @param[in] Z Atomic number of the requested element.
    * @return    Symbol of the element specified by the given atomic number.
    */
   std::string        GetElementSymbol(const int Z)       const { return fNISTElementDataTable[Z-1].fElemSymbol; }

   /** @brief Public method to access number of known isotopes of a given element.
    *
    * @param[in] Z Atomic number of the requested element.
    * @return    Number of known isotopes of the element specified by the given atomic number.
    */
   int                GetNumberOfIsotopes(const int Z)    const { return fNISTElementDataTable[Z-1].fNumOfIsotopes; }

   /** @brief Public method to get the index of the most stable isotope of this element.
    *
    * If the element has no stable isotope we take the most stable isotope: when user build an element without
    * specifying the isotopes of the element, this database is used and isotopes of the given element with
    * non zero natural abundances are insterted into the element. However, there are elements that has no stable
    * isotope. In this case the this procedure would lead to an element without having any isotopes. IN order to avoid
    * this situtation, we will take the most stable isotope of the element and will inform the user.
    *
    * @param[in] Z Atomic number of the requested element.
    * @return    Index of the most stable isotope if the elemnet has no stable isotope. Will return with -1 if the
    *            element has stable isotope.
    */
   int                GetIndexOfTheMostStableIsotope(const int Z)    const { return fNISTElementDataTable[Z-1].fIndxOfMostStableIsotope; }

   /** @brief Public method to get the flag that indicates if this element has no stable isotope.
    *
    * See more at GetIndexOfTheMostStableIsotope().
    *
    * @param[in] Z Atomic number of the requested element.
    * @return    True if the element has no stable isotope. The index of the most stable isotope can be got through GetIndexOfTheMostStableIsotope().
    *            False if the element has stable isotope.
    */
   bool                IsHaveNoStableIsotope(const int Z)    const { return fNISTElementDataTable[Z-1].fIsNoStableIsotope; }

   /** @brief Public method to access symbols of all known isotopes of a given element.
    *
    * @param[in] Z Atomic number of the requested element.
    * @return    Symbols of all known isotopes of the element specified by the given atomic number.
    *            Size of the array is equal to GetNumberOfIsotopes().
    */
   const std::string* GetIsotopeSymbols(const int Z)      const { return fNISTElementDataTable[Z-1].fSymbols; }

   /** @brief Public method to access nucleon numbers of all known isotopes of a given element.
    *
    * @param[in] Z Atomic number of the requested element.
    * @return    Nucleon numbers of all known isotopes of the element specified by the given atomic number.
    *            Size of the array is equal to GetNumberOfIsotopes().
    */
   const int*         GetIsotopeNucleonNums(const int Z)  const { return fNISTElementDataTable[Z-1].fNIsos; }

   /** @brief Public method to access atomic masses of all known isotopes of a given element.
    *
    * Atomic masses include the mass of the electron shell as well. Atomic masses are stored in
    * internal [weigth] units and can be converted to internal [weight/amount of substance] units
    * multiplying by the Avogadro(kAvogadro) number OR to internal [energy] unit multiplying by
    * speed of light square (kCLightSquare).
    *
    * @param[in] Z Atomic number of the requested element.
    * @return    Atomic masses of all known isotopes of the element specified by the given atomic number.
    *            Size of the array is equal to GetNumberOfIsotopes().
    */
   const double*      GetIsotopeAtomicMasses(const int Z) const { return fNISTElementDataTable[Z-1].fAIsos; }

   /** @brief Public method to access isotope masses of all known isotopes of a given element.
    *
    * Isotope masses are the masses of the bare nucleus in internal [energy] units i.e. without
    * the electron shell and the corresponding total electron binding energy.
    *
    * @param[in] Z Atomic number of the requested element.
    * @return    Masses of all known isotopes of the element specified by the given atomic number.
    *            Size of the array is equal to GetNumberOfIsotopes().
    */
   const double*      GetIsotopeMasses(const int Z) const { return fNISTElementDataTable[Z-1].fMassIsos; }

   /** @brief Public method to access natural abundances of all known isotopes of a given element.
    *
    * @param[in] Z Atomic number of the requested element.
    * @return    Natural abundances of all known isotopes of the element specified by the given atomic number.
    *            Size of the array is equal to GetNumberOfIsotopes().
    */
   const double*      GetIsotopeNaturalAbundances(const int Z) const { return fNISTElementDataTable[Z-1].fWIsos; }

   /** @brief Public method to access the mean atomic mass of a given element.
    *
    * The mean atomic mass is the natural abundance weighted mean isotope atomic mass.
    * The mean atomic mass is stored in internal [weight/amount of substance] units.
    *
    * @param[in] Z Atomic number of the requested element.
    * @return    Mean atomic mass of the element specified by the given atomic number.
    */
   double             GetMeanAtomicMass(const int Z) const { return fNISTElementDataTable[Z-1].fMeanAtomicMass; }

   /** @brief Public method to access the atomic mass of a given isotope.
    *
    * The atomic mass is returned in internal [weight/amount of substance] units
    * and includes the electron shell.
    *
    * @param[in] Z Atomic number of the requested isotope.
    * @param[in] N Nucleon number of the requested isotope.
    * @return    Atomic mass of the isotope specified by the given atomic number
    *            and number of nucleons.
    */
   double GetAtomicMass(int Z, int N);

   /** @brief Public method to access the mass of a given isotope.
    *
    * It is the mass of the bare nucleus and returned in internal [energy] units.
    *
    * @param[in] Z Atomic number of the requested isotope.
    * @param[in] N Nucleon number of the requested isotope.
    * @return    Mass of the isotope specified by the given atomic number
    *            and number of nucleons.
    */
   double GetIsotopeMass(int Z, int N);

   /** @brief Public method to access the total electron binding energy of a given isotope.
    *
    * The total electron binding energy is returned in internal [energy] units.
    *
    * @param[in] Z Atomic number of the requested isotope.
    * @param[in] N Nucleon number of the requested isotope.
    * @return    Total electron binding energy of the isotope specified by the
    *            given atomic number and number of nucleons.
    */
   double GetBindingEnergy(int Z, int N);
//@}

/**
 * @name Already built NIST element realted getters/setters:
 */
//@{
  /** @brief Public method to access the already built NIST element indices in the global element table.
   *
   * When NIST elements are built through the Element::NISTElement() method, their indices
   * (in the global element table) will be stored here.
   *
   * @param[in] Z Atomic number of the requested NIST element.
   * @return
   *            - Index of the specified NIST element if that NIST element has already built.
   *            - -1 otherwise.
   */
  int    GetNISTElementIndex(int Z) const {return fIndicesOfBuiltNISTElements[Z-1];}

  /** @brief Public method to store the NIST element index in the global element table.
   *
   * When NIST elements are built through the Element::NISTElement() method, their indices
   * (in the global element table) will be stored here.
   *
   * @param[in] Z Atomic number of the NIST element.
   * @param[in] indx Index of the NIST element in the global element table.
   */
  void   SetNISTElementIndex(int Z, int indx) {fIndicesOfBuiltNISTElements[Z-1] = indx;}
//@}


/**
 * @name Printouts:
 */
//@{
   /** @brief Public method to print some element data related to its isotopic composition.
    *
    *
    * @param[in] Z Atomic number of the requested element.
    */
   void PrintData(int Z);
//@}

/**
 * @name Delted methods:
 */
//@{
  /** @brief Copy constructor (deleted to avoid implicit declaration and public
   *         for better error msgs).
   */
   NISTElementData(const NISTElementData&) = delete;
   /** @brief Assignment operator (deleted to avoid implicit declaration and public
    *         for better error msgs).
    */
   NISTElementData& operator=(const NISTElementData&) = delete;
//@}

 private:
   /** @brief Constructor. */
   NISTElementData();
   /** @brief Internal method to set up the database. */
   void BuildTable();

 private:
   /** @brief Number of elements in the database. */
   static const int gNumberOfNISTElements = 118;

   /** @brief Internal data structure to represent one element isotopic composition. */
   struct NISTElement{
     /** @brief Symbol of this element. */
     std::string        fElemSymbol;
     /** @brief Symbols for each isotopes [fNumOfIsotopes]. */
     const std::string  *fSymbols;
     /** @brief Atomic number of this element. */
     int                fZ;
     /** @brief Number of known isotopes of this element. */
     int                fNumOfIsotopes;
     /** @brief Index of the most stable isotope if this element has no stable isotope. */
     int                fIndxOfMostStableIsotope;
     /** @brief Flag to indicate that this element has no stable isotope. */
     bool               fIsNoStableIsotope;
     /** @brief Mean atomic mass of this element in internal [weight/mole] units. */
     double             fMeanAtomicMass;
     /** @brief Number of nucleons for each known isotopes of this element [fNumOfIsotopes]. */
     const int          *fNIsos;
     /** @brief Atomic masses in internal [weight] units for each known isotopes of this element [fNumOfIsotopes]. */
     const double       *fAIsos;
     /** @brief Natural abundances for each known isotopes of this element [fNumOfIsotopes]. */
     const double       *fWIsos;
     /** @brief Masses for each known isotopes of this element in internal [energy] units [fNumOfIsotopes]. */
     double             *fMassIsos;
   }/** @brief Table to store isotopic composition data for each element.*/
    fNISTElementDataTable[gNumberOfNISTElements];

   /** @brief Total electron binding energies for each elements in internal [energy] units. */
   double fBindingEnergies[gNumberOfNISTElements];

   /** @brief Indices of NIST elements that has already built by Element::NISTElement(). */
   int fIndicesOfBuiltNISTElements[gNumberOfNISTElements];
};

} // namespace geant

#endif // NISTELEMENTDATA_H
