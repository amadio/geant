
#ifndef NISTMATERIALDATA_H
#define NISTMATERIALDATA_H

#include <string>
// change to vecgeom::map
//#include <map>
#include "base/Map.h"


namespace geantphysics {
  /**
   * @brief   Internal(singletone) database to store NIST material data.
   * @class   NISTMaterialData
   * @author  M Novak, A Ribon
   * @date    december 2015
   *
   * Internal database to store material composition data for several materials
   * taken from NIST database [http://physics.nist.gov/cgi-bin/Star/compos.pl].
   * These materials can be constructed through the Material::NISTMaterial()
   * method by specifying their name. The list of material names available in
   * the database can be found  \ref NISTMaterialDataDoc.
   */

enum class MaterialState;

class NISTMaterialData{

 public:
   /** @brief Public method to access the singletone database instance. */
   static NISTMaterialData& Instance();

/**
* @name Public getters:
*/
//@{
   /** @brief Public method to get the number of materials in the database. */
   int  GetNumberOfNISTMaterials() const { return gNumberOfNISTMaterials; }

   /** @brief Public method to get the number of elemental materials in the database. */
   int  GetNumberOfElementalNISTMaterialData() const { return gNumberOfElementalNISTMaterials; }

   /** @brief Method to find the index of a NISTMaterialData structure given by material name.
    *  @param[in]  name Name of the NIST material.
    * @return
    *              - Index of the corresponding NIST material data structure in the internal database.
    *                if the database contains data with the given NIST material name.
    *              - -1 otherwise.
    */
   int FindNISTMaterialDataIndex(const std::string &name);

   /** @brief Public mathod to get the name of the NIST material.
    *  @param[in] indx Index of the NIST matrial in the database.
    *  @return Name of the NIST material.
    */
   const std::string& GetName(const int indx) const                {return fNISTMaterialDataTable[indx].fName;}

   /** @brief Public mathod to get the density of the NIST material.
    *  @param[in] indx Index of the NIST matrial in the database.
    *  @return Density of the NIST material in internal [weight/length/f$^3/f$] units.
    */
   double            GetDensity(const int indx) const              {return fNISTMaterialDataTable[indx].fDensity;}

   /** @brief Public mathod to get the mean excitation energy of the NIST material.
    *  @param[in] indx Index of the NIST matrial in the database.
    *  @return Mean excitation energy of the NIST material in internal [energy] units.
    */
   double            GetMeanExcitationEnergy(const int indx) const {return fNISTMaterialDataTable[indx].fMeanExcitationEnergy;}

   /** @brief Public mathod to get the temperature of the NIST material.
    *  @param[in] indx Index of the NIST matrial in the database.
    *  @return Temperature of the NIST material in internal [temperature] units.
    */
   double            GetTemperature(const int indx) const          {return fNISTMaterialDataTable[indx].fTemperature;}

   /** @brief Public mathod to get the pressure of the NIST material.
    *  @param[in] indx Index of the NIST matrial in the database.
    *  @return Pressure of the NIST material in internal [pressure] units.
    */
   double            GetPressure(const int indx) const             {return fNISTMaterialDataTable[indx].fPressure;}

   /** @brief Public mathod to get the number of components this NIST material is built up.
    *  @param[in] indx Index of the NIST matrial in the database.
    *  @return Number of components(elements) this NIST material is built up.
    */
   int               GetNumberOfComponents(const int indx) const   {return fNISTMaterialDataTable[indx].fNumComponents;}

   /** @brief Public mathod to get the state of this NIST material.
    *  @param[in] indx Index of the NIST matrial in the database.
    *  @return State of this NIST material.
    */
   MaterialState     GetMaterialState(const int indx) const        {return fNISTMaterialDataTable[indx].fState;}

   /** @brief Public mathod to get the list of atomic number of the (NIST)elements this NIST material is built up.
    *  @param[in] indx Index of the NIST matrial in the database.
    *  @return List of atomic number of the (NIST)elements this NIST material is built up [fNumComponents].
    */
   const int*        GetListOfElements(const int indx) const       {return fNISTMaterialDataTable[indx].fElementList;}

   /** @brief Public mathod to get the list of mass fractions or number of atoms (depending on the fIsBuiltByAtomCount)
    *         for each constituting elements.
    *  @param[in] indx Index of the NIST matrial in the database.
    *  @return List fraction/count for each element this material is built up [fNumComponents].
    */
   const double*     GetListOfElementFractions(const int indx) const {return fNISTMaterialDataTable[indx].fElementFraction;}

   /** @brief Public mathod to get if this NIST material is specified by atom counts or by fractional masses.
    *  @param[in] indx Index of the NIST matrial in the database.
    *  @return
    *           - True if this material is specified by atom counts of the elments.
    *           - False if this material is specified by mass fractions of the elements.
    */
   bool              IsToBuildByAtomCount(const int indx) const       {return fNISTMaterialDataTable[indx].fIsBuiltByAtomCount;}
//@}

/**
 * @name Delted methods:
 */
//@{
   /** @brief Copy constructor (deleted to avoid implicit declaration and public
    *         for better error msgs).
    */
   NISTMaterialData(const NISTMaterialData&) = delete;
   /** @brief Assignment operator (deleted to avoid implicit declaration and public
    *         for better error msgs).
    */
   NISTMaterialData& operator=(const NISTMaterialData&) = delete;
//@}

 private:
   /** @brief Constructor. */
   NISTMaterialData();
   /** @brief Internal method to set up the database. */
   void BuildTable();

 private:
   /** @brief Number of materials in the database. */
   static const int gNumberOfNISTMaterials = 316;

   /** @brief Number of elemental materials in the database (they are the first)*/
   static const int gNumberOfElementalNISTMaterials = 98;

   /** @brief Internal data structure to represent one NIST material composition and basic properties. */
   struct NISTMaterial{
   /** @brief Name of this NIST material. */
   std::string        fName;
   /** @brief Density of this NIST material in internal [weight/length/f$^3/f$] units. */
   double             fDensity;
   /** @brief Mean excitation energy of this NIST material in internal [energy] units. */
   double             fMeanExcitationEnergy;
   /** @brief Temperature of this NIST material. */
   double             fTemperature;
   /** @brief Pressure of this NIST material. */
   double             fPressure;
   /** @brief Number of components(elements) this NIST material is composed. */
   int                fNumComponents;
   /** @brief State of this NIST material. */
   MaterialState      fState;
   /** @brief List of atomic number of the (NIST)elements this NIST material is built up [fNumComponents]. */
   const int          *fElementList;
   /** @brief List of mass fractions or number of atoms (depending on the fIsBuiltByAtomCount) for each
    *         constituting elements [fNumComponents]. */
   const double       *fElementFraction;
   /** @brief Flag to indicate that in this case of NIST material atom counts are stored in fElementW. */
   bool               fIsBuiltByAtomCount;

  }/** @brief Table to store material composition data for some NIST materials. */
   fNISTMaterialDataTable[gNumberOfNISTMaterials];

   /** @brief Internal map to store NIST material indices in the internal NIST material
    *         database with a key = material name. */
  // change to vecgeom::map
  // std::map<const std::string,int> fMapNISTMaterialNameToIndex;
  vecgeom::map<const char*,int> fMapNISTMaterialNameToIndex;
};

} // namespace geant

#endif // NISTMATERIALDATA_H
