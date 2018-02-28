#ifndef MATERIALCUTS_H
#define MATERIALCUTS_H

// for inline namespace VECGEOM_IMPL_NAMESPACE
#include "base/TypeMap.h"
// for inlice namespace GEANT_IMPL_NAMESPACE
#include "Geant/Config.h"

#include <vector>
#include <iostream>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {
class Region;
}
}

namespace geantphysics {
inline namespace GEANT_IMPL_NAMESPACE {
class Material;
}
}

namespace geantphysics {

// class Material;
/**
 * @brief   Material - particle production cut object (dummy class at the moment).
 * @class   MaterialCuts
 * @author  M Novak, A Ribon
 * @date    january 2016
 *
 * Object to store partcile production cuts in material. Production cuts are currently available for gamma, electron and
 * positron particles and stored both in energy threshold and length. Each constructed material - cuts object pointer is
 * stored in a global table that can be obtained through the static MaterialCuts::GetTheMaterialCutsTable method. Each
 * matrial-cuts object stores its own index in the global table that can be obtained by MaterialCuts::GetIndex method.
 * Each matrial-cuts object stores a pointer to the material object (does not own the material object) that can be
 * obtained by MaterialCuts::GetMaterial.
 *
 * Production cuts in energy/lenght are stored in arrays index by particle type:
 * \begin{listings}
 *  \mathrm{index} = 0  & \mathrm{gamma particle} \\
 *  \mathrm{index} = 1  & \mathrm{electron}       \\
 *  \mathrm{index} = 2  & \mathrm{positron}
 * \end{listings}
 *
 * The energy/length production cuts arrays can be obtained by
 * MaterialCuts::GetProductionCutsInEnergy/MaterialCuts::GetProductionCutsInLength. Production cuts in energy/length
 * for a given particle can be set by MaterialCuts::SetProductionCutEnergy/MaterialCuts::SetProductionCutLength by
 * providing the particle index (see above) and the appropriate production cut in appropriate internal units.
 *
 *
 * \todo: this is a very basic implementation with the minimal functionality. It will be developed further together with
 *        the general physics framework and functionalities like converting production cuts from length to energy
 *        threshold or possibility to set production cuts either in energy threshold or length will be added.
 *        More detailed documentation will be provided later with the implementation of these new functionalities.
 */
class MaterialCuts {
public:
  // creates all MaterialCuts and converts length/energy to energy/lenght
  static void CreateAll();
  // deletes all MaterialCuts objects and set the table to default (zero) size
  static void ClearAll();

  /** @brief Public method to obtain the index of this material-cuts object in the global table. */
  int GetIndex() const { return fIndex; }
  int GetRegionIndex() const { return fRegionIndex; }

  const double *GetProductionCutsInLength() const { return fProductionCutsInLength; }
  const double *GetProductionCutsInEnergy() const { return fProductionCutsInEnergy; }
  bool IsProductionCutsGivenInLength() const { return fIsProductionCutsGivenInLength; }

  const Material *GetMaterial() const { return fMaterial; }
  // get a MaterialCuts object pointer by its index
  static const MaterialCuts *GetMaterialCut(int indx);
  // get a MaterialCuts object by specifying the Region index and the Material index
  // TODO: we should get rid of it by writing the MaterialCuts global index into the Reagion and requested the kernel
  //       to provide it at each physics call
  static const MaterialCuts *GetMaterialCut(int regionindx, int materialindx);

  // get the global mategrial-cuts table
  static const std::vector<MaterialCuts *> &GetTheMaterialCutsTable() { return gTheMaterialCutsTable; }

  /**
   * @name Printouts:
   */
  //@{
  friend std::ostream &operator<<(std::ostream &, const MaterialCuts *);
  friend std::ostream &operator<<(std::ostream &, const MaterialCuts &);
  friend std::ostream &operator<<(std::ostream &, std::vector<MaterialCuts *>);
  //@}

private:
  /**
    * @brief Dummy constructor for testing physics models. Will be removed.
    *
    * Production cuts are set to the provided values.
    *
    * @param[in] mat          Pointer to specify the material object part of this mategrial-cuts pair.
    * @param[in] gcutlength   Production cut for gamma particle in length.
    * @param[in] emcutlength  Production cut for electron in length.
    * @param[in] epcutlength  Production cut for positron in length.
    * @param[in] gcutenergy   Production cut for gamma particle in energy.
    * @param[in] emcutenergy  Production cut for electron in energy.
    * @param[in] epcutenergy  Production cut for positron in energy.
    */
  MaterialCuts(int regionindx, const Material *mat, bool iscutinlength, double gcut, double emcut, double epcut);
  /** @brief Destructor */
  ~MaterialCuts(){};

  // NOTE: they might be private if we create all matcuts automatically
  // if not   // TODO: add check of indices
  void SetProductionCutEnergy(int indx, double val) { fProductionCutsInEnergy[indx] = val; }
  void SetProductionCutLength(int indx, double val) { fProductionCutsInLength[indx] = val; }

  // checks if MaterialCuts has already been created for the given Material in this Region:
  // if not, it will create a new MaterialCuts and will return a pointer to the corresponding MaterialCuts object
  static MaterialCuts *CheckMaterialForRegion(const vecgeom::Region *region, const Material *mat);
  // convert gamma,e-,e+ cuts length to energy or energy to length
  static void ConvertAll();

private:
  static const int kNumProdCuts = 3;
  int fRegionIndex;                    // index of the region where this material-cuts is located
  int fIndex;                          // in the global table
  bool fIsProductionCutsGivenInLength; // is production cuts are given by the user in length
  double fProductionCutsInLength[kNumProdCuts];
  double fProductionCutsInEnergy[kNumProdCuts];

  const Material *fMaterial; // does not own this material onject

  // the global mategrial-cuts table; filled constatntly when MaterialCuts obejcts are created
  static std::vector<MaterialCuts *> gTheMaterialCutsTable;
  // the material-cuts table per Regions; size will be [#regions][#MaterialCuts in that Region]
  // filled constatntly when MaterialCuts obejcts are created
  static std::vector<std::vector<MaterialCuts *>> gTheMaterialCutsPerRegion;
  // the material-cuts table per Regions; size will be [#regions][#Matrials]
  // filled constatntly when MaterialCuts obejcts are created
  static std::vector<std::vector<MaterialCuts *>> gTheMaterialCutsPerRegionPerMaterial;
};

} // namespace geantphysics

#endif // MATERIALCUT_H
