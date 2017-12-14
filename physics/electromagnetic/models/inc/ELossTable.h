
#ifndef ELOSSTABLE_H
#define ELOSSTABLE_H

#include <vector>

namespace geantphysics {

/**
 * @brief   Class to build, store and provide energy loss relted data like dE/dx, range, inverse range.
 * @class   ELossTable
 * @author  M Novak, A Ribon
 * @date    august 2016
 *
 * Each ELossTable object handels energy loss related tables like restricted stopping power, restricted range, restric-
 * ted inverse range for a set of MaterialCuts for all partciles. There will be one ELossTable objcet created in the
 * ELossTableManager for each different PhysicsParameters(PhysicsList) object and the number of MaterialCuts size
 * ELossTable* vector elements of the ELossTableManager will be set to point to the one that the MaterialCuts with the
 * same index belongs to.
 * Energy loss related data for a given MaterialCuts, Particle and kinetic energy can be obtained at run-time by
 * interpolation through the GetRestrictedDEDX(), GetRestrictedRange(), GetEnergyForRestrictedRange() methods, that are
 * called by the corresponding methods of the ELossTableManager singletone.
 * Full CSDA range tables per material are also built and available at run-time if the fIsComputeCSDARange flag in the
 * PhysicsParameters object is set to true before initialization.
 */

// forward declarations
class PhysicsParameters;
class EMPhysicsProcess;
class MaterialCuts;
class Particle;
class GLIntegral;
class Spline;

class ELossTable {

public:
  /**
   * @brief CTR; As many ELossTable objects are created in the ELossTableManager as many different PhysicsParameters
   *        objects exists. Note, each PhysicsList has its own PhysicsParameters object.
   */
  ELossTable(PhysicsParameters *physpar);
 ~ELossTable();

  /** @brief  Run time method to obtain restricted stopping power for the given MaterialCuts index, particle and kinetic
   *          energy.
   *
   * The method is expected to be called by the ELossTableManager.
   * The restricted stopping power table is built at initialization for each particle that has at least one EnergyLoss
   * process (i.e. ionization or bremsstrahlung) for each MaterialCuts that belong to region where the corresponding
   * kEnergyLoss process/es is/are active over a kinetic energy grid defined by the PhysicsParameters object that
   * belongs to the corresponding region. The run-time value is obtained by spline interpolation from this table.
   *
   * @param[in] matcutindx  Index of MaterialCuts object in which the restricted stopping power is requested.
   * @param[in] partindx    Internal index of the particle for which the restricted stopping power is requested.
   * @param[in] kinenergy   Kinetic energy of the particle at which the restricted stopping power is requested.
   * @return    Restricted total stopping power in the specified MaterialCuts for the specified particle at the given
   *            kinetic energy if the corresponding stopping power table was built in internal [energy/length] units.
   *            Zero otherwise.
   */
  double GetRestrictedDEDX(int matcutindx, int partindx, double kinenergy);

  /** @brief  Run time method to obtain restricted range for the given MaterialCuts index, particle and kinetic energy.
   *
   * The method is expected to be called by the ELossTableManager.
   * The restricted range table is built at initialization, by integration of the restricted stopping power table, for
   * each particle that has at least one EnergyLoss process (i.e. ionization or bremsstrahlung) for each MaterialCuts
   * that belong to region where the corresponding kEnergyLoss process/es is/are active over a kinetic energy grid
   * defined by the PhysicsParameters object that belongs to the corresponding region. The run-time value is obtained
   * by spline interpolation from this table.
   *
   * @param[in] matcutindx  Index of MaterialCuts object in which the restricted range is requested.
   * @param[in] partindx    Internal index of the particle for which the restricted range is requested.
   * @param[in] kinenergy   Kinetic energy of the particle at which the restricted range is requested.
   * @return    Restricted range in the specified MaterialCuts for the specified particle at the given kinetic energy
   *            if the corresponding range table was built in internal [length] units. A big (1.0e+20) otherwise.
   */
  double GetRestrictedRange(int matcutindx, int partindx, double kinenergy);


  /** @brief  Run time method to obtain the kinetic energy that corresponds to a given restricted range in the
   *          MaterialCuts specified by its index.
   *
   * The method is expected to be called by the ELossTableManager.
   * The restricted inverse range table (i.e. kinetic energy v.s. range) is based on the the resctricted range table
   * (i.e. range v.s. kinetic energy) built at initialization, by integration of the restricted stopping power table,
   * for each particle that has at least one EnergyLoss process (i.e. ionization or bremsstrahlung) for each
   * MaterialCuts that belong to region where the corresponding kEnergyLoss process/es is/are active over a kinetic
   * energy grid defined by the PhysicsParameters object that belongs to the corresponding region. The run-time value
   * is obtained by spline interpolation set up at initialisation as kinetic energy v.s. restricted range.
   * It is used to compute the energy loss of partciles (having continuous energy losses) after traveling a given s path
   * E_final = What this method provides by passing range = Range(at the pre-step point energy)-s.
   *
   * @param[in] matcutindx  Index of MaterialCuts object in which the kinetic energy is requested.
   * @param[in] partindx    Internal index of the particle for which the kinetic energy is requested.
   * @param[in] range       Restricted range of the particle at which the corresponding kinetic energy is requested.
   * @return    Kinetic energy of the particle that corresponds to the provided range in the specified MaterialCuts
   *            internal [energy] units if the restricted range table was built. Zero otherwise (should indicate that
   *            there is no any energy losses i.e. should return with the current energy).
   */
  // TODO:: think about this! now each time there is a binary search to find the lower index of the range parameter
  //       - at the moment we try to avoid to use this such that we provide methods like GetMeanEnergyAfterAStep...
  //        that should be used instead of more fragmented computations
  // returns with the kinetic energy that corresponds to a given restricted range
  double GetEnergyForRestrictedRange(int matcutindx, int partindx, double range);

  /**
   * @brief  Method to obtain full(CSDA) range for the given Material index, particle and kinetic energy. Available only
   *         if the corresponding (full) range table was requested to be built by setting the fIsComputeCSDARange member
   *         of the corresponding PhysicsParameters object before initialisation.
   *
   * The method is expected to be called by the ELossTableManager.
   * The full (CSDA) range table is built at initialization if it was requested (false by default), by computing the
   * full stopping power (i.e. setting the upper limit of the integral to be very large or equal to the primary kinetic
   * energy) and from that the corresponding full CSDA range. These range tables are built for each particle that has
   * at least one kEnergyLoss process (i.e. ionization or bremsstrahlung) for each Material that belong to region where
   * the corresponding kEnergyLoss process/es is/are active. The run-time value is obtained by spline interpolation from
   * this table.
   *
   * @param[in] matindx     Index of the Material object in which the CSDA range is requested.
   * @param[in] partindx    Internal index of the particle for which the CSDA range is requested.
   * @param[in] kinenergy   Kinetic energy of the particle at which the rCSDA range is requested.
   * @return    Full CSDA range in the specified Material for the specified particle at the given kinetic energy
   *            if the corresponding range table was built in internal [length] units. A high value (1.0e+20) otherwise.
   */
  double GetRange(int matindx, int partindx, double kinenergy);

/*
//
//  THESE ARE ONLY FOR TESTING
//
 void PrintRestrictedDEDX(int matcutindx, int partindx) {
   ELossData *lossData = nullptr;
   if (fELossDataPerMaterialCutsPerParticle[matcutindx].size()>partindx &&
       (lossData = fELossDataPerMaterialCutsPerParticle[matcutindx][partindx])) {
     for (int i=0; i<lossData->fNumData; ++i) {
       std::cout<< std::setprecision(16)<<lossData->fEnergyGridData[i]/geant::MeV << "  "
                << std::setprecision(12)<< lossData->fRestrictedDEDXData[i]/(geant::MeV/geant::mm)
                <<std::endl;
     }
   } else {
     std::cerr << " ======= ELossTable: no ELossData for Particle = "
               << Particle::GetParticleByInteralCode(partindx)->GetName()
               << "  in MaterialCuts: \n";
     std::cerr << MaterialCuts::GetTheMaterialCutsTable()[matcutindx] <<std::endl;
   }
 }

 void PrintRestrictedRange(int matcutindx, int partindx) {
   ELossData *lossData = nullptr;
   if (fELossDataPerMaterialCutsPerParticle[matcutindx].size()>partindx &&
       (lossData = fELossDataPerMaterialCutsPerParticle[matcutindx][partindx])) {
     for (int i=0; i<lossData->fNumData; ++i) {
       std::cout<< std::setprecision(6)<< lossData->fEnergyGridData[i]/geant::MeV << "  "
                << std::setprecision(12)<< lossData->fRestrictedRangeData[i]/(geant::mm)
                <<std::endl;
     }
   } else {
     std::cerr << " ======= ELossTable: no ELossData for Particle = "
               << Particle::GetParticleByInteralCode(partindx)->GetName()
               << "  in MaterialCuts: \n";
     std::cerr << MaterialCuts::GetTheMaterialCutsTable()[matcutindx] <<std::endl;
   }
 }

 void PrintRange(int matindx, int partindx) {
   if (!fIsComputeCSDARange) {
     std::cerr << "\n ======= ELossTable::PrintRange : CSDA Rnage computation was not requested. \n"
               << " Make sure that CSDA range computation is set to requested in PhysicsParameters \n"
               << " because no Range data are avaialbe for Particle = "
               << Particle::GetParticleByInteralCode(partindx)->GetName()
               << "  in Material: \n";
               std::cerr << Material::GetTheMaterialTable()[matindx] <<std::endl;
     return;
   }
   ELossData *lossData = nullptr;
   if (fELossDataPerMaterialPerParticle[matindx].size()>partindx
       && (lossData = fELossDataPerMaterialPerParticle[matindx][partindx])
       && lossData->fRangeData
      ) {
     for (int i=0; i<lossData->fNumData; ++i) {
       std::cout<< std::setprecision(16)<< lossData->fEnergyGridData[i]/geant::MeV << " "
                << std::setprecision(16)<< lossData->fRangeData[i]*Material::GetTheMaterialTable()[matindx]->GetDensity()/(geant::g/(geant::cm*geant::cm))
                <<std::endl;
     }
   } else {
     std::cerr << " ======= ELossTable::PrintRange : no Range data are avaialbe for Particle = "
               << Particle::GetParticleByInteralCode(partindx)->GetName()
               << "  in Material: \n";
     std::cerr << Material::GetTheMaterialTable()[matindx] <<std::endl;
   }
 }
*/

  /**
   * @brief Method to build and set up all energy loss related tables for all particles that has kEnergyLoss process(es)
   *        for a set of MaterialCuts (that belongs to region(s) that this ELossTable handels). The method is called
   *        from the ELossTableManager.
   *
   * @param[in]  elosstablespermatcut Vector from ELossTableManager that stores the energy loss related tables for each
   *                                  MaterialCuts.
   */
  void BuildELossTable(std::vector<ELossTable*> &elosstablespermatcut);

  // data structure to store both restricted and total energy loss related data;
  // - (restricted)ELossData are created and stored per partcile per MaterialCuts:
  //     # restricted dedx
  //     # the corresponding restricted range and inverse range
  // - (total)ELossData are created and stored per partcile per Material if it was set to requested in the
  //   PhysicsParameters through the SetIsComputeCSDARange(false/true) method (by def. false):
  //     # total range of the particle in the given material that corresponds to the total dedx
  struct ELossData {
    int                  fNumData;             // size of data arrays determined by some members of the fPhysicsParameters
    const MaterialCuts  *fMaterialCuts;        // pointer to the MaterialCuts object this data belongs to; the calss do
                                               // NOT own the MaterialCuts object;
    const Particle      *fParticle;            // pointer to the particle object this data belongs to; the class do NOT
                                               // own the Particle object
    double              *fEnergyGridData;      // pointer to the common energy grid determined by some members of the
                                               // fPhysicsParameters; the struct do NOT own this data because there is
                                               // only one energy grid per ELossTable that is owned by the ELossTable
    // the struct do own the following data i.e. need to delete them.
    double              *fRestrictedDEDXData;    // restricted dedx for the given particle
                                                 // and MaterialCuts; the struct own the data;
    double              *fRestrictedRangeData;   // restricted range for the given particle
                                                 // and MaterialCuts; the struct own the data;
    double              *fRangeData;             // set only if total data was requested in BuildOneELossData i.e. if
                                                 // the parameter fIsComputeCSDARange was set to true in the
                                                 // corresponding PhysicsParameters object;
                                                 // total range of the particle in the given material that corresponds
                                                 // to the total  dedx; the struct do own the data;
    // spline interpolators
    Spline              *fSplineRestrictedDEDX;    // spline interpolator to obtain
                                                   // restricted dedx values at run time for the given partcile and
                                                   // MaterialCuts; the struct do own the data;
    Spline              *fSplineRestrictedRange;   // spline interpolator to obtain
                                                   // restricted range values at run time for the given partcile and
                                                   // MaterialCuts; the struct do own the data;
    Spline              *fSplineRestrictedInvRange;// spline interpolator to obtain
                                                   // restricted inverse range values at run time for the given partcile
                                                   // and MaterialCuts; the struct do own the data;
    Spline              *fSplineRange;             // set only if total data was requested in BuildOneELossData;
                                                   // spline interpolator to obtain total range values at run time for
                                                   // the given partcile and Material; the struct do own the data;
    std::vector<EMPhysicsProcess*>  fLossProcesses; // list of energy loss processes that are active in the set of regions
                                                    // handeled by this ELossTable for the given partcile; the class do NOT
                                                    // own this PhysicsProcess-es

  };

//
// private methods
//
private:

  /**
   * @brief Method to set the common energy grid for all energy loss related tables.
   *
   * The properties of the common energy grid like min, max values or number of bins are determined by some of the
   * PhysicsParameters class members that this ELossTable belongs to.
   */
  void InitializeEnergyGrid();


  /**
   * @brief Helper method to build and set up one ELossData structure i.e. for a given MaterialCuts and Particle.
   *
   * The method is called from the BuildELossTable where the pointed ELossData structure (pointer is an input parameter)
   * is created and its number of data, MaterialCuts, Particle, energy grid has been already set.
   *
   * @param[in] lossdata           Pointer to one ELossData structure (created and set up i.e. energy grid, grid size,
   *                               MatrialCuts and Particle pointers are set properly in the BuildELossTable() before).
   * @param[in] iscomputetotaldata Flag to indicate if the total data i.e. full CSDA range data must be computed for as well.
   */
  void BuildOneELossData(ELossData *lossdata, bool iscomputetotaldata);


  /**
   * @brief Helper method to compute restricted stopping power table for a given MaterialCuts and Particle.
   *
   * The method is expected to be called from the BuildOneELossData helper method to handle the dEdx computation part.
   * Computes dE/dx by summing up contributions from all kEnergyLoss EMPhysicsProcess(es) assigned to the given Particle.
   * The ComputeDEDX method of the registred kEnergyLoss EMPhysicsProcess-es are invoked and a spline interpollatior is
   * set up on the computed dEdx table for run-time interpolation.
   *
   * @param[in] lossdata  Pointer to the ELossData structure with energy grid, MaterialCuts, Particle and kEnergyLoss
   *                      process pointer list already set properly and this method fills the restricted dEdx table part
   *                      and then sets up a spline interpolator on the dEdx table.
   */
  void BuildRestrictedDEDXTable(ELossData *lossdata);


  /**
   * @brief Helper method to compute restricted range table for a given MaterialCuts and Particle.
   *
   * The method is expected to be called from the BuildOneELossData helper method to handle the range computation part.
   * Computes range table by integrating data taken from the dE/dx table using 16 point GL integral at each energy bin.
   * Then sets up a spline interpolator for the run-time restricted range interpolation (range v.s. kinetic energy) as
   * well as a spline interpolator for the run-time inverse range interpolation (kinetic energy v.s. range).
   *
   * @param[in] lossdata  Pointer to the ELossData structure in which the restricted dedx has already been computed i.e.
   *                      the BuildRestrictedDEDXTable method has already been called for this ELossData structure. This
   *                      method will fill the restricted range table part and sets up interpolator for both restricted
   *                      range and restricted inverse range interpolations.
   */
  void BuildRestrictedRangeTable(ELossData *lossdata);


  /**
   * @brief Helper method to compute full CSDA range table for a given MaterialCuts and Particle.
   *
   * The total (CSDA) range is computed only if fIsComputeCSDARange is true in the corresponding PhysicsParameters
   * object. The full CSDA range is NOT used during the simulation (restricted range is used). The method is expected to
   * be called from the BuildOneELossData helper method to handle the full range computation part:
   *  The total dedx i.e. that corresponds to production cut equal to the particle kinetic energy is computed by summing
   *  up contributions from all registered EnergyLoss EMPhysicsProcess; the ComputeDEDX method of the registred
   *  EnergyLoss EMPhysicsProcess-es are invoked with indicating that total dedx i.e. production cut equal to the
   *  partcile energy is requested.
   *  The inverse total dedx is integrated to obtain the corresponding total range.
   *  A spline interpollatior is also set on the total range data.
   *
   * @param[in] lossdata  Pointer to one ELossData structure (created and set up i.e. energy grid, grid size,
   *                      MatrialCuts and Particle pointers are set properly in the BuildELossTable() before).
   */
  void  BuildTotalRangeTable(ELossData *lossdata);

  /**
   * @brief Method to clean up all allocated memory, data structures and reset container sizes.
   */
  void Clear();

//
// private data members
//
private:
  PhysicsParameters   *fPhysicsParameters;  // PhysicsParameters that this ELossTable need to use; the class do NOT own
                                            // the PhysicsParameters object


  bool        fIsComputeCSDARange;
  int         fNGL;
  int         fNumLossTableBins;
  double      fMinLossTableEnergy;
  double      fMaxLossTableEnergy;
  double      fLogMinLossTableEnergy;
  double      fEnergyILDelta;
  double     *fEnergyGrid;     // common energy grid determined by some of the members of the PhysicsParameters; owned



  // Restricted data:
  // first dimension will be number of MaterialCuts in the global MaterialCuts table; element, with index that
  // correspond to a MatrialCuts that belongs to any of the regions where this ELossTable is active (active regions
  // can be obtained from the PhysicsParameters member) can store a vector of ELossData structure per partcile indexed
  // by particle internal code;
  // fELossDataPerMaterialCutsPerParticle.size()    :
  //                                     => number of MatrialCuts in the global MatrialCuts table
  // fELossDataPerMaterialCutsPerParticle[i].size() :
  //                                     => 0 if the MatrialCuts with index =i does not belong to region where this
  //                                          ELossTable is active; these MatrialCuts are handeled by other ELossTable
  //                                     =  maximum internal code of the particle that has any EnergyLoss
  //                                        EMPhysicsProcess registered in ELossTableRegister
  // fELossDataPerMaterialCutsPerParticle[i].[j]    :
  //                                     = ELossData structure for Partcile with internal partcile code j and for
  //                                       MatrialCuts with index =i if the Partcile had any EnergyLoss EMPhysicsProcess
  //                                       registered in the ELossTableRegister with active region set equal to set of
  //                                       regions where this ELossTable is active
  //                                     = nullptr otherwise, i.e. there is no any EnergyLoss EMPhysicsProcess-es
  //                                       registered for the Particle with internal particle code j in
  //                                       ELossTableRegister OR the registered EnergyLoss processes are active in
  //                                       different set of regions where this ELossTable is active
  std::vector<std::vector<ELossData*> > fELossDataPerMaterialCutsPerParticle;
  // Unrestricted data:
  // first dimension will be number of Material-s in the global Material table; element, with index that
  // correspond to a Matrial that belongs to any of the regions where this ELossTable is active (active regions
  // can be obtained from the PhysicsParameters member) can store a vector of ELossData structure per partcile indexed
  // by particle internal code;
  // fELossDataPerMaterialPerParticle.size()        :
  //                                     => number of Matrial in the global Matrial table
  // fELossDataPerMaterialPerParticle[i].size()
  //                                     => 0 if the Matrial with index =i does not belong to region where this
  //                                          ELossTable is active; these Matrial-s are handeled by other ELossTable
  //                                     =  maximum internal code of the particle that has any EnergyLoss
  //                                        EMPhysicsProcess registered in ELossTableRegister
  // fELossDataPerMaterialPerParticle[i].[j]
  //                                     = ELossData structure for Partcile with internal partcile code j and for
  //                                       Matrial with index =i if the Partcile had any EnergyLoss EMPhysicsProcess
  //                                       registered in the ELossTableRegister with active region set equal to set of
  //                                       regions where this ELossTable is active
  //                                     = nullptr otherwise, i.e. there is no any EnergyLoss EMPhysicsProcess-es
  //                                       registered for the Particle with internal particle code j in
  //                                       ELossTableRegister OR the registered EnergyLoss processes are active in
  //                                       different set of regions where this ELossTable is active
  // The stored ELossData* are shared with fELossDataPerMaterialCutsPerParticle and deleted properly when that is cleared.
  std::vector<std::vector<ELossData*> > fELossDataPerMaterialPerParticle;

  GLIntegral   *fGL;
};

}  // namespace geantphysics

#endif // ELOSSTABLE_H
