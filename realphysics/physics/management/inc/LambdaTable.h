
#ifndef LAMBDATABLE_H
#define LAMBDATABLE_H

#include <vector>

namespace geantphysics {

/**
 * @brief   Optional per-process table to build, handle macroscopic cross sections pre-computed at initialisation time.
 * @class   LambdaTable
 * @author  M Novak
 * @date    July 2017
 */

class Spline;
class PhysicsProcess;
class Particle;
class MaterialCuts;

class LambdaTable {
public:

  /** @brief Constructor.
   *
   * @param[in] process       Pointer to the process objcet this table belongs to.
   * @param[in] ispermaterial Flag to indicate if the table is requested to build per-material (default). Tables will
   *                          be built per-material-cuts otherwise.
   */
  LambdaTable(const PhysicsProcess *process, bool ispermaterial=true);

  /** @brief Destructor. */
 ~LambdaTable();

  /** @brief Public method to provide possibility for some processes to set sepcial number of lambda table bins. Must be
    *        called before initialisation of the corresponding process.
    *
    * @param[in] val Value of the requested number of lambda table bins (will be adujted)
    */
  void SetSpecialLambdaTableBinNum(int val) {
    fIsSpecialLambdaTableBinNum = true;
    fNumSpecialLambdaTableBins  = val;
  }

  /** @brief Public method to provide (interpolated) macroscopic cross section for the given process in a given material
    *         or material-cuts, particle kinetic energy.
    *
    * @param[in]  matcut  Material-cuts in which the cross section is required.
    * @param[in]  ekin    Kinetic energy at which the cross section is required (in internal [energy] unit).
    * @return     Macroscopic cross section for the specified input of the process-particle couple which the table
    *             belongs to (in internal [1/length] unit).
    */
  double  GetMacroscopicXSection(const MaterialCuts *matcut, double ekin);

  /** @brief Public method to provide the kinetic energy at which the macroscopic cross section of the process-particle
    *        couple reach its maximum value in the given material/material-cuts.
    * @param [in] matcut  Pointer to the material-cuts object to specify the material or material-cuts.
    * @return     Kinetic energy at which the cross section has its maximum in the given material/material-cuts (in
    *             internal [energy] unit.)
    */
  double  GetMacroscopicXSectionMaximumEnergy(const MaterialCuts *matcut) const;

  /** @brief Public method to provide the maximum value of the macroscopic cross section of the process-particle couple
    *        in the given material/material-cuts.
    * @param [in] matcut  Pointer to the material-cuts object to specify the material or material-cuts.
    * @return     Maximum value of the macroscopic cross section in the given material/material-cuts (in internal
    *             [energy] unit.)
    */
  double  GetMacroscopicXSectionMaximum(const MaterialCuts *matcut) const;


  /** @brief Public method to build the macroscopic cross section tables (called automatically at the initialisation of
    *        the process this table talongs to)
    */
  void    BuildLambdaTables();

private:
  /** brief Private method to clear all tables built by calling BuildLambdaTables .*/
  void   ClearAllTables();
  // deleted methods
  LambdaTable() = delete;
  LambdaTable(const LambdaTable &other) = delete;
  LambdaTable& operator=(const LambdaTable &other) = delete;

private:
  /** @brief Pointer to the physics process object this table belongs to */
  const PhysicsProcess *fProcess;          // not owned by the object
  /** @brief Flag to indicate if this table is built per-material (per-material-cuts otherwise) */
  bool        fIsLambdaTablesPerMaterial;  // by def true and if false then per material-cuts
  /** @brief Flag to indicate if the process requested special number of kinetic energy bins */
  bool        fIsSpecialLambdaTableBinNum; // by def false
  //
  // these are used only if the lambda tables are requested to build per material because then the kinetic energy grid
  // is common for each material

/**
 * @name Members to describe one macroscopic cross section table if the tables were requested to build by the process
 *       per-material. NOTE:
 *       - these members are not used if the tables were requested to build per-material-cuts
 *       - the kinetic energy grid is common for all tables when the tables are requested to build per-materials
 */
//@{
  /** @brief Number of (kinetic energy and cross section) bins in the table. */
  int                       fNumLambdaTableBins;
  /** @brief Special number of bins in the table (also in case of tables per material-cuts but only requested). */
  int                       fNumSpecialLambdaTableBins;
  /** @brief Minimum kinetic energy of the table. */
  double                    fMinLambdaTableEnergy;
  /** @brief Maximum kinetic energy of the table. */
  double                    fMaxLambdaTableEnergy;
  /** @brief Logarithm of the minimum kinetic energy of the table. */
  double                    fLogMinLambdaTableEnergy;
  /** @brief Inverse of the logarithmic bin width of the kinetic energy grid. */
  double                    fEnergyILDelta;
  /** @brief The kinetic energy grid (size is fNumLambdaTableBins). */
  std::vector<double>       fEnergyGrid;
  /** @brief Data structure to describe one cross section table ie for one material */
  struct ALambdaTable {
    /** @brief Kinetic energy at which the macr. cross section reach its maximum. */
    double                    fLambdaMaxEnergy;
    /** @brief Maximum value of the macr. cross section (used to account possible energy loss along the step) */
    double                    fLambdaMax;
    /** @brief Macroscopic cross section values for one given material (size is fNumLambdaTableBins). */
    std::vector<double>       fOneLambdaTable;
    /** @brief Pointer to a spline interpolator object set up of this cross section table. */
    Spline                   *fSpline;
  };
//@}

/**
 ** @name Data structure to describe one macroscopic cross section table if the tables were requested to build by the
 *       process per-material-cuts. NOTE:
 *       - these members are not used if the tables were requested to build per-material
 *       - each table has its own kinetic energy grid when the tables are requested to build per-materials-cuts
 */
//@{
  struct LambdaTableForAMaterialCuts {
    /** @brief Number of bins in the table (accounts special bin numbers if requested). */
    int                       fNumLambdaTableBins;
    /** @brief Minimum kinetic energy of the table. */
    double                    fMinLambdaTableEnergy;
    /** @brief Maximum kinetic energy of the table. */
    double                    fMaxLambdaTableEnergy;
    /** @brief Logarithm of the minimum kinetic energy of the table. */
    double                    fLogMinLambdaTableEnergy;
    /** @brief Inverse of the logarithmic bin width of the kinetic energy grid. */
    double                    fEnergyILDelta;
    /** @brief Kinetic energy at which the macr. cross section reach its maximum. */
    double                    fLambdaMaxEnergy;
    /** @brief Maximum value of the macr. cross section (used to account possible energy loss along the step) */
    double                    fLambdaMax;
    /** @brief The kinetic energy grid (size is fNumLambdaTableBins). */
    std::vector<double>       fEnergyGrid;
    /** @brief Macroscopic cross section values for one given material (size is fNumLambdaTableBins). */
    std::vector<double>       fLambdaTable;
    /** @brief Pointer to a spline interpolator object set up of this cross section table. */
    Spline                   *fSpline;
  };
//@}

  /** @brief Container to store lambda tables per-material indexed by material index. Only those indices will be non-
   *        nullptr that correspond to material for which the model needs to respond at run-time. (used only if tables
   *        were requested to build by the process per-material.
   */
  std::vector<ALambdaTable*>                   fLambdaTablesPerMaterial;
  /** @brief Container to store lambda tables per-material-cuts indexed by material-cut index. Only those indices will be
   *        non-nullptr that correspond to material-cuts for which the model needs to respond at run-time. (used only
   *        if tables were requested to build by the process per-material-cuts.
   */
  std::vector<LambdaTableForAMaterialCuts*>    fLambdaTablesPerMaterialCuts;

private:
  /** @brief Private method to generate the kinetic energy grid of the tables */
  void GenerateEnergyGrid(const MaterialCuts *matcut, struct LambdaTableForAMaterialCuts *data=nullptr);

};

}        // LAMBDATABLE_H

#endif   // LAMBDATABLE_H
