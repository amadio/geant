#ifndef MATERIALCUTS_H
#define MATERIALCUTS_H

#include <vector>

namespace geant {

class Material;
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
 * Production cuts in energy/lenght are stored in arrays and indexed by particle type:
 * \f[
 *  \mathrm{index} = 
 * \begin{cases}
 *  0  & \mathrm{gamma}\; \mathrm{particle} \\
 *  1  & \mathrm{electron}       \\
 *  2  & \mathrm{positron}
 * \end{cases}
 * \f]
 *
 * The energy/length production cuts arrays can be obtained by
 * MaterialCuts::GetProductionCutsInEnergy / MaterialCuts::GetProductionCutsInLength. Production cuts in energy/length
 * for a given particle can be set by MaterialCuts::SetProductionCutEnergy / MaterialCuts::SetProductionCutLength by
 * providing the particle index (see above) and the appropriate production cut in appropriate internal units.
 *
 *
 * \todo: this is a very basic implementation with the minimal functionality. It will be developed further together with
 *        the general physics framework and functionalities like converting production cuts from length to energy
 *        threshold or possibility to set production cuts either in energy threshold or length will be added.
 *        More detailed documentation will be provided later with the implementation of these new functionalities.
 */
class MaterialCuts{
public:
  /**
    * @brief Constructor.
    *
    * Production cuts are set to default (1 mm) for all particles.
    *
    * @param[in] mat  Pointer to specify the material object part of this mategrial-cuts pair.
    *
    */
  MaterialCuts(Material *mat);
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
  MaterialCuts(Material *mat, double gcutlength, double emcutlenght, double epcutlength, 
               double gcutenergy  = -1.0, double emcutenergy = -1.0, double epcutenergy = -1.0);
  /** @brief Destructor */
 ~MaterialCuts(){};

  /** @brief Public method to obtain the index of this material-cuts object in the global table. */
  int             GetIndex() const { return fIndex; }
  void            SetProductionCutEnergy(int indx, double val) {fProductionCutsInEnergy[indx] = val;}
  void            SetProductionCutLength(int indx, double val) {fProductionCutsInLength[indx] = val;}

  const double*   GetProductionCutsInLength() const { return fProductionCutsInLength;}
  const double*   GetProductionCutsInEnergy() const { return fProductionCutsInEnergy;}

  const Material* GetMaterial() const { return fMaterial; }

 // get the global mategrial-cuts table
 static const std::vector<MaterialCuts*>& GetTheMaterialCutsTable() { return gTheMaterialCutsTable; }

private:
  static const int kNumProdCuts = 3;
  int    fIndex; // in the global table
  double fProductionCutsInLength[kNumProdCuts];
  double fProductionCutsInEnergy[kNumProdCuts];

  Material *fMaterial;

  // the global mategrial-cuts table
  static std::vector<MaterialCuts*> gTheMaterialCutsTable;
};

}  // namespace geant

#endif // MATERIALCUT_H
