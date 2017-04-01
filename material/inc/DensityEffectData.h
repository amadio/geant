
#ifndef DENSITYEFFECTDATA_H
#define DENSITYEFFECTDATA_H

#include "Types.h"

#include <string>

namespace geantphysics {
  /**
   * @brief   Internal(singletone) database to store density effect parameters.
   * @class   DensityEffectData
   * @author  M Novak, A Ribon
   * @date    January 2016
   *
   * Internal database to store parameters that are needed to compute the density
   * effect correction parameter in the MaterialProperties. Data are
   * taken from \cite sternheimer1984density and
   * available for simple (elemental) materials [composed of single element with
   * Z=[1,97] with no data available for Z=85 and Z=87] and a set of compounds.
   * Each material can be identified by name (that starts with the NIST_MAT_
   * prefix) or the index of the elemental materials can be retrieved by atomic
   * number. See \ref DensityEffectDataDoc for more details regarding materials
   * parameter values available in this database.
   *
   *\warning (It is a warning here but I put it under todo to see this under Related Pages as well)
   *\todo
   * These parameters were determined in the reference (see below) through
   * a fitting procedure. The fitted equations contains some parameters that a
   * given fixed value was used for during the fit:
   * - \f$ I \f$ -mean excitation energy
   * - \f$ -C = 2\ln(I/(h\nu_{p})) +1 \f$) where \f$ h\nu_{p} \f$ is the plasma energy
   *   that depends on the density and molar mass of the material
   *\warning (It is a warning here but I put it under todo to see this under Related Pages as well)
   *\todo
   * The values of these parameters (mean excitation energy, density and molar mass
   * of the material) were fixed in the reference to those given in the NISTMaterialData.
   * Therefore, these parameters will lead to incorrect density effect parameter for
   * materials with mean excitation energy and/or density and/or molar mass different
   * from those given in the NISTMaterialData. Therefore, we use these parameters only
   * for NIST materials at the moment in order to avoid inconsistencies.
   *
   * \note (It is a note but I put it under todo to see this under Related Pages as well)
   * \todo Geant4 will use these values not only in the case of name match i.e.
   * not only for NIST materials but also for simple (elemental) materials irrespectively
   * if the above conditions are true or not. But this can cause inconsistency and lead
   * to inaccurate density effect parameter. We will always compute the parameters in case
   * of non NIST materials. Actually, even in the case of NIST materials we never take
   * parameter C from here to avoid possible inconsistencies coming from a re-setted
   * mean excitation energy that different from the original NIST material data value.
   *
   */

enum class MaterialState;

class DensityEffectData{

 public:
  /** @brief Public method to access the singletone database instance. */
  static DensityEffectData& Instance();

  /** @brief Public method to retrieve the index of an elemental material in the database.
   *  @param[in]  z Atomic number of the single element the material is consist of.
   *  @param[in]  state State of the requested material.
   *  @return
   *             - Index of the density effect data in the internal database if
   *               the statements are true below:
   *                - the database contains data for the given simple material
   *                - the requested state either match the stored state or MaterialState::kStateUndefined
   *             - -1 otherwise
   */
  int    GetElementalIndex(int z, MaterialState state) const;

  /** @brief Public method to retrieve the index of a density effect data in the database by material name.
   *  @param[in]  name Name of the requested material. (NIST_MAT_ prefix format).
   *  @return
   *              - Index of the density effect data in the internal database
   *                if the database contains data with the given material name.
   *              - -1 otherwise.
   */
  int    FindDensityEffectDataIndex(const std::string &name);


  /** @brief Get name of the material.*/
  const std::string& GetName(const int indx) const { return fDensityEffectDataTable[indx].fName; }

  /** @brief Plasma energy \f$ h\nu_{p} \f$ of the material in internal [energy] units.
   *
   *  Important to note, that the plasma energy of the electrons of the material considered as free
   *  electron gas:
   *  \f[
   *     h\nu_{p} = \sqrt{4\pi n_{el} r_e} \hbar c
   *  \f]
   * where \f$ n_{el}\f$ is the electron density of the material and \f$ r_e \f$ is the classical electron radius.
   * The electron density for simple (elemental) materials
   * \f[
   *    n_{el} = Z \frac{\mathcal{N}_{av}\rho}{A}
   * \f]
   * and compound materials
   * \f[
   *    n_{el} = \sum_{i}^{\#elements} Z_i \frac{\mathcal{N}_{av}w_i\rho}{A_i}
   * \f]
   * where \f$ \mathcal{N}_{av} \f$ is the Avogadro number, \f$ \rho \f$ is the material density,
   * \f$ w_i \f$ is the fractional mass of the \f$ i^{th} \f$ element of the material
   * with molar mass of \f$ A_i \f$ . Therefore, the plasma energy depends on the electron density
   * that depends on material density and molar mass. So this plasma energy will be correct only in the case
   * if the these data (density, molar mass) are exactly the same as in the corresponding NIST material.
   */
  double GetPlasmaEnergy(const int indx)      const { return fDensityEffectDataTable[indx].fPlasmaEnergy;}

  /** @brief Get the Sternheimer adjustment factor of the material. */
  double GetAdjustmentFactor(const int indx)  const { return fDensityEffectDataTable[indx].fSternheimerFactor;}

  /** @brief Get density parameter \f$ -C \f$ of the material.
   *
   *  Important to note, that density parameter \f$ -C=2\ln(I/(h\nu_{p})) +1 \f$ where \f$ I \f$ is the mean excitation energy
   *  of the material, \f$ h\nu_{p} \f$ is the plasma energy. The plasma energy depends on the material density and molar mass.
   *  Therefore, this density parameter will be corret only in the case if the mean excitation energy, density and molar mass
   *  of the material are exactly the same as in the corresponding NIST material.
   */
  double GetParameterC(const int indx)        const { return fDensityEffectDataTable[indx].fParameterC;}

  /** @brief Get the fitted density parameter \f$ X_0 \f$ .
   *
   *  This is one of the parameters in Eqs.(9,10,14) of the reference so its value was determined throught fit. Since
   *  parameter \f$ C \f$ is also involved in these equations all the comments made at GetParameterC() are true for
   *  for this parameter as well.
   */
  double GetParameterX0(const int indx)       const { return fDensityEffectDataTable[indx].fParameterFitX0;}

  /** @brief Get the fitted density parameter \f$ X_1 \f$ .
   *
   *  This is one of the parameters in Eqs.(9,10,14) of the reference so its value was determined throught fit. Since
   *  parameter \f$ C \f$ is also involved in these equations all the comments made at GetParameterC() are true for
   *  for this parameter as well.
   */
  double GetParameterX1(const int indx)       const { return fDensityEffectDataTable[indx].fParameterFitX1;}

  /** @brief Get the fitted density parameter \f$ a \f$ .
   *
   *  This is one of the parameters in Eqs.(9,10,14) of the reference so its value was determined throught fit. Since
   *  parameter \f$ C \f$ is also involved in these equations all the comments made at GetParameterC() are true for
   *  for this parameter as well.
   */
  double GetParameterA(const int indx)        const { return fDensityEffectDataTable[indx].fParameterFitA;}

  /** @brief Get the fitted density parameter \f$ m \f$ .
   *
   *  This is one of the parameters in Eqs.(9,10,14) of the reference so its value was determined throught fit. Since
   *  parameter \f$ C \f$ is also involved in these equations all the comments made at GetParameterC() are true for
   *  for this parameter as well.
   */
  double GetParameterM(const int indx)        const { return fDensityEffectDataTable[indx].fParameterFitM;}

  /** @brief Get the density parameter \f$ \delta_0 \f$ .
   *
   *  This is the \f$ \delta_0 \f$ parameter in Eqs.(14) of the reference.
   */
  double GetParameterDelta0(const int indx)   const { return fDensityEffectDataTable[indx].fParameterDelta0;}

  /** @brief Get the upper bound of the error of the density effect parameter \f$ \delta \f$ inherent in the fitting procedure. */
  double             GetDeltaErrorMax(const int indx)     const { return fDensityEffectDataTable[indx].fDeltaErrorMax;}

  /** @brief Get the state of this material. */
  MaterialState      GetMaterialState(const int indx)     const { return fDensityEffectDataTable[indx].fState;}



/**
 * @name Delted methods:
 */
//@{
   /** @brief Copy constructor (deleted to avoid implicit declaration and public
    *         for better error msgs).
    */
   DensityEffectData(const DensityEffectData&) = delete;
   /** @brief Assignment operator (deleted to avoid implicit declaration and public
    *         for better error msgs).
    */
   DensityEffectData& operator=(const DensityEffectData&) = delete;
//@}

 private:
   /** @brief Constructor. */
   DensityEffectData();
   /** @brief Internal method to set up the database. */
   void BuildTable();

 private:
   /** @brief Number of materials in the database that we have density effect data. */
   static const int gNumberOfDensityEffectData = 278;
   /** @brief Maximum atomic number of simple(elemental) material that we have data for.
    *
    *  Data are available from Z=1 to Z=97 but there is no data for Z=85 and Z=87.
    */
   static const int gMaxElementalZet           =  97;  // Z=85 and Z=87 are not available

   struct DensityEffectParameter{
   std::string fName;         // material name
   double fPlasmaEnergy;      // \f$ h\nu_{p} \f$ plasma energy in internal [energy] units
   double fSternheimerFactor; // Sternheimer adjustment factor
   double fParameterC;        // \f$ =2\ln(I/(h\nu_{p})) +1 \f$ where \f$ I \f$ is the mean excitation energy
   double fParameterFitX0;    // 'X0' parameter of the fit [in Eqs.(9,10) and (14) in ref.]
   double fParameterFitX1;    // 'X1' parameter of the fit [in Eqs.(9,10) and (14) in ref.]
   double fParameterFitA;     // 'a' parameter of the fit [in Eqs.(9,10) and (14) in ref.]
   double fParameterFitM;     // 'm' parameter of the fit [in Eqs.(9,10) and (14) in ref.]
   double fParameterDelta0;   // '\f$ \delta_0 \f$' parameter of the fit [in Eq.(14) in ref.]
   double fDeltaErrorMax;     // upper bound of the error inherent in the fitting procedure.
   MaterialState fState;      // state of the material (only for simple elemental materials fStateUndefined otherwise)
  } fDensityEffectDataTable[gNumberOfDensityEffectData];

  /** @brief Internal map to store density effect data indices in the internal DB with a key = material name. */
  MapHelper<std::string,int>::Map_t fMapMaterialNameToDenistyEffectDataIndex;
};


} // namespace geantphysics


#endif // DENSITYEFFECTDATA_H
