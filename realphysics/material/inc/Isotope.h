

#ifndef ISOTOPE_H
#define ISOTOPE_H


#include <vector>
#include <map>
#include <string>

namespace geantphysics {
  /**
   * @brief   Class to describe an isotope.
   * @class   Isotope
   * @author  M Novak, A Ribon
   * @date    december 2015
   *
   * Isotopes must be unique i.e. there cannot be more than one isotope with the
   * the same values of Z(atomic numer), N(number of nucleons) and isomer level.
   * In order to ensure this:\n
   * isotopes cannot be created directly by the user. Instead, isotopes can be
   * accessed by the user through the static GetIsotope() method by specifying the
   * properties of the requested isotope. This method will return with a pointer
   * to the requested isotope (the requested isotope will be created internaly
   * only if it has not been created yet).
   *
   * Each isotope pointer is stored in a global isotope table when created.
   * This table can be accessed by the user through the static GetTheIsotopeTable()
   * method.
   *
   * \warning It's strongly not recommended to delete any isotopes! Each isotope will
   * be deleted automatically at the end of the run.
   *
   * \todo \li We might make the destructor private. With this we can prevent
   *           the user to delete any of the isotopes during the run.
   *       \li Then we can have a static DeleteAll() method that will delete all
   *           isotopes at the end of the run. Even better if we have a material
   *           manager that is friend of the Isotope (and Element and Material)
   *           class and that manager will have the DeleteAllIsotopes() method.
   *       \li Then we might have a CreateAll() method (or in the manager friend)
   *           that can create all the isotopes that has not been created yet.So
   *           first the user creates some of the isotopes (e.g. with the prefered
   *           atomic mass) then all the remaining isotopes can be created through
   *           this method.
   *      \li This manager friend class can be called at physics initialisation and
   *          at the end of the run.
   */

class Isotope{
 public:
/**
 * @name Construction related method:
 */
//@{
    /**
     * @brief The only public method to get an isotope specified by its properties.
     *
     * Isotopes must be unique i.e. there cannot be more than one isotope with the
     * the same values of Z(atomic numer), N(number of nucleous) and isomer level.
     * in order to ensure this, isotopes are created internally(private ctr) when
     * a given isotope is requested first time by the user through this public
     * GetIsotope() method. Each isotope pointer is stored in a global
     * isotope table at construction. If an isotope, that has already been
     * created earlier (i.e. it can be found in the global isotope table),
     * is requested by the user then the pointer of the corresponding isotope is
     * returned.
     *
     * @param[in] z Atomic numer.
     * @param[in] n Number of nucleons.
     * @param[in] a Atomic mass in internal [weight/amount of substance] units.
     * @param[in] isol Isomer level.
     * @param[in] isomass Mass of the bare nucleus in internal [energy] units.
     * @param[in] name Name/symbol of the isotope.
     * @return    Pointer to the requested isotope specified by the given input
     *            parameter values.
     */
    static Isotope* GetIsotope(int z, int n, double a = 0., int isol = 0,
                               double isomass = 0., const std::string &name = "");
//@}


    /** @brief Destructor.
     *
     *  \warning It's strongly not recommended to delete any isotopes! Each
     *  isotope will be deleted automatically at the end of the run.
     */
   ~Isotope();



/**
 * @name Field related getters:
 */
//@{
    /** @brief  Public method to access the name/symbol of this isotope. */
    const std::string& GetName()  const {return fName;}

    /** @brief  Public method to access the atomic number of this isotope. */
    int    GetZ() const {return fZ;}

    /** @brief  Public method to access the umber of nucleons in this isotope. */
    int    GetN() const {return fN;}

    /** @brief Public method to access the atomic mass of this isotope.
     *
     *  The atomic mass will be returned in internal [weight/amount of substance]
     *  units and the mass of the electron shell as well as the corresponding
     *  total binding energy is included.
     */
    double GetA() const {return fA;}

    /** @brief Public method to access the isomer level of this isotope. */
    int    GetIsoLevel() const {return fIsoL;}

    /** @brief Public method to access the mass of this isotope.
     *
     * This method will return with the mass of the bare nucleus of this isotope
     * i.e. without the electron shell and the corresponding total binding
     * energy in internal [energy] units.
     */
    double GetIsoMass() const {return fIsoMass;}

    /** @brief Public method to access the index of this isotope.
     *
     * Each isotopes with unique values of Z(atomic numer), N(number of nucleous)
     * and isomer level are stored in a global isotope table (that can be
     * accessed through the GetTheIsotopeTable() method) when created. This
     * method returns with the index of this isotope in that global table.
     */
    int    GetIndex() const {return fIndex;}
//@}

/**
 * @name Methods to get data regarding the global isotope table:
 */
//@{
    /** @brief Public method to access the global isotope table. */
    static const std::vector<Isotope*>& GetTheIsotopeTable() {return gTheIsotopeTable;}

    /**
     * @brief Public method to delete all Element objects that has been created.
     *
     * This method is called by the PhysicsProcessHandler when it is deleted to clean all Element objects.
     * Users should never call this method!
     */
    static void ClearAllIsotopes();

    /** @brief Public method to access the number of isotopes in the global
      *        isotope table.
      * @return Number of isotopes in the global isotope table i.e. number of
      *         isotopes have been created so far.
      */
    static int   GetNumberOfIsotopes() {return gTheIsotopeTable.size();}
//@}

/**
 * @name Printouts:
 */
//@{
    friend std::ostream& operator<<(std::ostream&, const Isotope*);
    friend std::ostream& operator<<(std::ostream&, const Isotope&);
    friend std::ostream& operator<<(std::ostream&, std::vector<Isotope*>);
//@}

 private:
   /** @brief Isotope constructor(private)
     *
     * Isotopes are created internally when a given isotope is requested first
     * time by the user through the GetIsotope() public method.
     *
     * @param[in] name Name/symbol of the isotope.
     * @param[in] z Atomic numer.
     * @param[in] n Number of nucleons.
     * @param[in] a Atomic mass in internal [weight/amount of substance] units.
     * @param[in] isomass Mass of the bare nucleus in internal [energy] units.
     * @param[in] isol Isomer level.
     */
   Isotope(const std::string &name, int z, int n, double a = 0.,
           double isomass = 0., int isol = 0);

   /** @brief Copy constructor. Not implemented yet. */
   Isotope(const Isotope&);              // not implemented

   /** @brief Equal operator. Not implemented yet.*/
   Isotope& operator=(const Isotope&);   // not implemented

   /** @brief Private method to compute key value for a given isotope
    *
    *  Each isotopes stored in a global isotope table when created. The index
    *  of the isotope in this global table is stored in a map with a key value
    *  computed from z(atomic mass), n(nucleon number) and isol(isomer level).
    *
    * @param[in] z Atomic numer.
    * @param[in] n Number of nucleons.
    * @param[in] isol Isomer level.
    * @return    Key value of the specified isotope.
    */
   static int GetKey(int z, int n, int isol=0) { return 10000*z + 10*n + isol;}

   /** @brief Private method to get the index of the specified isotope.
    *
    * This method is used internally to find the index of the specified isotope
    * in the global isotope table. It will return -1 if the isotope has not been
    * created yet.
    *
    * @param[in] z Atomic numer.
    * @param[in] n Number of nucleons.
    * @param[in] isol Isomer level.
    * @return     \li If the isotope has already been created: index of the isotope
    *                in the global isotope table.
    *            \li Otherwise: -1
    */
   static int   GetIsotopeIndex(int z, int n, int isol =0);

 //
 // data members
 private:
    std::string   fName;     /** @brief Name/symbol of the isotope */
    int           fZ;        /** @brief Atomic number */
    int           fN;        /** @brief Number of nucleons */
    int           fIsoL;     /** @brief Isomer level */
    int           fIndex;    /** @brief Index in this isotope in the global table */
    double        fA;        /** @brief Atomic mass in internal [weight/mole] units (including electron shell) */
    double        fIsoMass;  /** @brief Mass of the bare nucleus in internal [energy] units */

    /** The global isotope table that stores the pointers of all isotopes that
     *  has already been created.
     */
    static std::vector<Isotope*>  gTheIsotopeTable;
    /** Global map to store the indices of the already created isotopes in the
     * global isotope table. The key values are computed from Z(atomic mass),
     * N(nucleon number) and isoL(isomer level).
     */
    static std::map<int,int>      gTheIsotopeMap;
};

} // namespace geantphysics

#endif // ISOTOPE_H
