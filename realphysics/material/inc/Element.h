
#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include <string>

namespace geantphysics {

class Isotope;

  /**
   * @brief   Class to describe an element.
   * @class   Element
   * @author  M Novak, A Ribon
   * @date    december 2015
   *
   * Elements can be created:
   * 1. Directly, without specifying its isotopes (natural isotopes of the
   *    element will be taken from the NIST database i.e. NISTElementData
   *    utomatically).
   * 2. By adding its isotopes one by one using the AddIsotope() method.
   * 3. Constructing NIST elements by using the NISTElement() method.
   *
   * Each element stores a list of pointers to the isotopes from which the element
   * is built up. This can be accessed by the user through the GetIsotopeVector()
   * method.
   *
   * Each element pointer is stored in a global element table when created.
   * This table can be accessed by the user through the GetTheElementTable() method.
   *
   * \warning It's strongly not recommended to delete any elements! Each element
   * will be deleted automatically at the end of the run.
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

class Element{
 public:
/**
 * @name Constructors, destructor and construction realted public methods:
 */
//@{
  /**
    * @brief Constructor to build an element directly i.e. no direct reference
    *        to isotopes.
    *
    * Natural isotopes of the element will be taken from the NIST database i.e.
    * NISTElementData utomatically.
    *
    * @param[in] name Name of this element.
    * @param[in] symbol Symbol of this element.
    * @param[in] zeff Effective atomic number of this element.
    * @param[in] aeff Effective atomic mass of this element in internal [weight/
    *            amount of substance] units if given.
    */
   Element(const std::string &name, const std::string &symbol, double  zeff,
           double  aeff = 0.0);

   /**
     * @brief Constructor to build an element from isotopes via AddIsotope() method.
     *
     * @param[in] name Name of this element.
     * @param[in] symbol Symbol of this element.
     * @param[in] numisotopes Number of isotopes this element will be built up.
     */
   Element(const std::string &name, const std::string &symbol, int numisotopes);

   /**
    * @brief Method to add the isotopes to the elements one by one.
    *
    * @param[in] isotope Pointer to the currently added isotope object.
    * @param[in] relabundance Relative abundance of this isotope.
    */
   void AddIsotope(Isotope* isotope, double relabundance);

    /** @brief Destructor.
     *
     *  \warning It's strongly not recommended to delete any elements! Each
     *  elements will be deleted automatically at the end of the run.
     */
   ~Element();
//@}


   /**
    * @brief Method to create/retrieve a NIST element.
    *
    * A NIST element has natural isotope composition defined in the NISTElementData,
    * effective atomic mass as the natural abundance weighted mean isotopic atomic
    * masses, normal element symbol and name starting with NIST_ELM_ + symbol.
    *
    * Uniqueness of NIST elements is guaranted:
    * when a NIST element is requested through this method, first it will be
    * checked if the requested element had already been created earlier:
    *  - new NIST element is not constructed if it has already been created earlier
    *    just a pointer to the already existing NIST element will be returned.
    *  - new NIST element will be created otherwise and its pointer will be
    *    returned.
    *
    * @param[in] zeff Effective atomic number of the requested NIST element.
    * @return A pointer to the requested NIST element.
    */
   static Element* NISTElement(double zeff);



/**
 * @name Field related getters/setters:
 */
 //@{
   /**
    * @brief Public method to get the name of this element.
    *
    * @return Name of this element.
    */
   const std::string& GetName()   const {return fName;}

   /**
    * @brief Public method to get the symbol of this element.
    *
    * @return Symbol of this element.
    */
   const std::string& GetSymbol() const {return fSymbol;}

   /**
    * @brief Public method to get the effective atomic number of this element.
    *
    * @return Effective atomic number of this element.
    */
   double GetZ()   const {return fZeff;}

   /**
    * @brief Public method to get the effective atomic mass of this element.
    *
    * @return Effective atomic mass of this element in internal [weight/amount of
    *         substance] units.
    */
   double GetA()   const {return fAeff;}

   /**
    * @brief Public method to get the index of this element in the global element table.
    *
    * @return Index of this element in the global element table.
    */
   int GetIndex() const {return fIndex;}

   /**
    * @brief Public method to check if natural isotope abundance was used to build
    *        this element.
    *
    * @return \li True: if natural isotope abundance was used to build this element.
    *         \li False: otherwise.
    */
   bool   GetNaturalAbundanceFlag()   const {return fIsNaturalAbundance;}

   /**
    * @brief Public method to set if natural isotope abundance was used to build
    *        this element.
    *
    * @param[in] val Boolean value indicating if natual isotope abundance was used
    *            to build this element.
    */
   void   SetNaturalAbundanceFlag(bool val) {fIsNaturalAbundance = val;}
//@}


/**
 * @name Methods to access data regarding the isotopes constituting this element:
 */
//@{
   /**
    * @brief Public method to get the number of isotopes constituting this element.
    *
    * It will return to the number of isotopes that has been added so far.
    *
    * @return Number of isotopes constituting this element.
    */
   int  GetNumberOfIsotopes() const {return fCurNumOfIsotopes;}
   // Get vector of relative abundance of each isotope in this element
   /**
    * @brief Public method to get the relative abundance of each isotope constituting
    *        this element.
    *
    * @return Relative abundance of each isotope constituting this element.
    *         (The size of the array is GetNumberOfIsotopes().)
    */
   const double* GetRelativeAbundanceVector() const {return fRelativeIsotopeAbundanceVector;}

   // Get vector of pointers to isotopes constituting this element
   /**
    * @brief Public method to get the list of isotopes that this element is built up.
    *
    * @return Vector of pointer to the isotope objects this element is built up.
    */
   const std::vector<Isotope*>& GetIsotopeVector() const {return fIsotopeVector;}
//@}


/**
 * @name Methods to get data regarding the global element table:
 */
//@{
   /**
    * @brief Public method to get the global element table.
    *
    * @return Vector of pointers to all the elements that has been created so far.
    */
   static const std::vector<Element*>& GetTheElementTable(){ return gTheElementTable;}
   /**
    * @brief Public method to delete all Element objects that has been created.
    *
    * This method is called by the PhysicsProcessHandler when it is deleted to clean all Element objects.
    * Users should never call this method!
    */
   static void ClearAllElements();
   // Get number of elements in the static element table
   /**
    * @brief Public method to get the number of elements in the global element table.
    *
    * @return Number of elements in the global element table i.e. number of elements
    *         have been created so far.
    */
   static int GetNumberOfElements() {return gTheElementTable.size();}
//@}


/**
 * @name Printouts:
 */
//@{
   friend std::ostream& operator<<(std::ostream&, const Element*);
   friend std::ostream& operator<<(std::ostream&, const Element&);
   friend std::ostream& operator<<(std::ostream&, std::vector<Element*>);
//@}


 private:
   /** @brief Copy constructor. Not implemented yet */
   Element(const Element&);               // not implemneted

   /** @brief Equal operator. Not implemented yet */
   Element& operator=(const Element&);    // not implemented

   /** @brief Helper method to initialise some mebers */
   void InitialiseMembers();

   /** @brief Method to add all natural isotopes to an element */
   void AddNaturalIsotopes();


 private:
   /** @brief Name of this element. */
   std::string fName;
   /** @brief Symbol of this element. */
   std::string fSymbol;
   /** @brief Effective atomic number of this element. */
   double      fZeff;
   /** @brief Effective atomic mass of this element in internal [weight/amount of substance] units. */
   double      fAeff;
   /** @brief Number of isotopes added to this element so far. */
   int         fCurNumOfIsotopes;
   /** @brief Index of this element in the global element table. */
   int         fIndex;
   /** @brief Flag to indicate is natural isotope abundances were used in case of this element. */
   bool        fIsNaturalAbundance;
   /** @brief Relative abundances for each isotope constituting this element [GetNumberOfIsotopes()]. */
   double*     fRelativeIsotopeAbundanceVector;
   /** @brief List of pointers to the isotope objects constituting this element. */
   std::vector<Isotope*> fIsotopeVector;

   /** @brief The global element table that contains all elements that has been created so far. */
   static std::vector<Element*> gTheElementTable; // the global element table
};

} // namespace geantphysics

#endif  // ELEMENT_H
