
#ifndef ELOSSTABLEMANAGER_H
#define ELOSSTABLEMANAGER_H

#include <vector>

namespace geantphysics {

class ELossTable;
class MaterialCuts;
class Particle;

/**
 * @brief   Singletone to hande energy ELossTable-s for all MaterialCuts. Each ELossTable will handle energy loss
 *          related tables for a set of MaterialCuts-s and in case of each MaterialCuts for all particles.
 * @class   ELossTableManager
 * @author  M Novak, A Ribon
 * @date    august 2016
 *
 * All energy loss related tables is built in the main ELossTableManager::BuildELossTables() method, that must be called
 * at initialisation from the PhysicsListManager after all processes(with their models) are initialised properly to
 * build all energy loss related tables.
 *
 * Each different physics list has its own PhysicsParameters so energy loss related data like min, max energy of the
 * energy loss table might be different. So as many ELossTable-s will be created as many different PhysicsParameters
 * objects and each of them will handle only those MaterialCuts and kEnergyLoss physics process(es) that belongs to
 * region in which the PhysicsParameters object is active. Energy loss data like dedx, range, inverse range ... can be
 * obtained at run-time through the GetRestrictedDEDX, GetRestrictedRange, GetEnergyForRestrictedRange methods by passing
 * pointers to the MaterialCuts and Particle object and specifying the kinetic energy of the particle.
 *
 * At the end of the run or in case of re-initialization one needs to call ELossTableManager::Clear() method to delete
 * all ELossTable objects and reset all container sizes! For makeing sure it happens, it is done in ELossTableManager::
 * BuildELossTables() before starting anything.
 */

class ELossTableManager {

public:
  /**
   * @brief Static method to obtain the singletone object.
   */
  static ELossTableManager& Instance();


  // copy CTR and assignment operators as deleted
  ELossTableManager(const ELossTableManager&) = delete;
  ELossTableManager& operator=(const ELossTableManager&) = delete;


  /**
   *  @brief The main method to build all energy loss tables for all regions (for all PhysicsParameters objects), for
   *         MaterialCuts, for all Particle(s) that have at leat one kEnergyLoss process registered in the
   *         ELossTableRegister.
   *
   * Each kEnergyLoss EMPhysicsProcess will be automatically registered in the ELossTableRegister singletone for the
   * particle that it is assigned to at initialisation of the EMPhysicsProcess. This method must be called at
   * initialisation from the PhysicsListManager after all processes(with their models) are initialised properly to
   * build all energy loss related tables.
   *
   * It will set up the number of MaterialCuts size vector of ELossTable*; then it reads the registered kEnergyLoss
   * EMPhysicsProcess-es per particle from the ELossTableRegister and set up the ELossTable objects and initilize each
   * of them.
   */
  void BuildELossTables();


  /**
   * @brief  Method to clear up all memory i.e. delete all ELossTable objects and reset container sizes. One needs to
   *         class this method at the end of the run or in case of re-initialization (before calling ELossTableManager::
   *         BuildELossTables()).
   */
  void Clear();


  /** @brief  Run time method to obtain restricted stopping power for the given MaterialCuts, Particle and kinetic
   *          energy.
   *
   * The restricted stopping power table is built at initialization for each MaterialCuts, for each particle that has at
   * least one kEnergyLoss process (i.e. ionization or bremsstrahlung) active in the reagion to which the given
   * MaterialCuts belongs to. The run-time value is obtained by spline interpolation from this table.
   *
   * @param[in] matcut    Pointer to the MaterialCuts object in which the restricted stopping power is requested.
   * @param[in] part      Pointer to the Particle for which the restricted stopping power is requested.
   * @param[in] kinenergy Kinetic energy of the particle at which the restricted stopping power is requested.
   * @return    Restricted total stopping power in the given MaterialCuts, Particle, kinetic energy(if the corresponding
   *            stopping power table was built) in internal [energy/length] units.
   *            Zero otherwise.
   */
  double GetRestrictedDEDX(const MaterialCuts *matcut, const Particle *part, double kinenergy);


  /** @brief  Run time method to obtain restricted range for the given MaterialCuts, Particle and kinetic energy.
   *
   * The restricted range table is built at initialization for each MaterialCuts, for each particle that has at least
   * one kEnergyLoss process (i.e. ionization or bremsstrahlung) active in the reagion to which the given MaterialCuts
   * belongs to. The run-time value is obtained by spline interpolation from this table.
   *
   * @param[in] matcut    Pointer to the MaterialCuts object in which the restricted range is requested.
   * @param[in] part      Pointer to the Particle for which the restricted range is requested.
   * @param[in] kinenergy Kinetic energy of the particle at which the restricted range is requested.
   * @return    Restricted range in the given MaterialCuts, Particle, kinetic energy(if the corresponding restricted
   *            range table was built) in internal [length] units.
   *            A high (1.0e+20) value otherwise.
   */
  double GetRestrictedRange(const MaterialCuts *matcut, const Particle *part, double kinenergy);


  /** @brief  Run time method to obtain the kinetic energy that corresponds to a given restricted range in the given
   *          MaterialCuts and Particle.
   *
   * The restricted range table is built at initialization for each MaterialCuts, for each particle that has at least
   * one kEnergyLoss process (i.e. ionization or bremsstrahlung) active in the reagion to which the given MaterialCuts
   * belongs to. The run-time value is obtained by spline interpolation set up as kinetic energy v.s. range.
   *
   * @param[in] matcut    Pointer to the MaterialCuts object in which the restricted range is requested.
   * @param[in] part      Pointer to the Particle for which the restricted range is requested.
   * @param[in] range     Restricted range of the particle at which the corresponding kinetic energy is requested.
   * @return    Kinetic energy of the particle that corresponds to the provided range in the MaterialCuts in internal
   *            [energy] units if the restricted range table was built. Zero otherwise (should indicate that there is no
   *            any energy losses i.e. should return with the current energy).
   */
  double GetEnergyForRestrictedRange(const MaterialCuts *matcut, const Particle *part, double range);


  /**
   * @brief   Method to obtain full (CSDA) range for the given Material(that specified by the MaterialCut), Particle and
   *          kinetic energy. Available only if the corresponding (full) range table was requested to be built by
   *          setting the fIsComputeCSDARange member of the corresponding PhysicsParameters object before initialisation.
   *
   * The full (CSDA) range table is built at initialization if it was requested (false by default), by computing the
   * full stopping power (i.e. setting the upper limit of the integral to be very large or equal to the primary kinetic
   * energy) and from that the corresponding full CSDA range. These range tables are built for each particle that has
   * at least one kEnergyLoss process (i.e. ionization or bremsstrahlung) for each Material that belong to region where
   * the corresponding kEnergyLoss process/es is/are active. The run-time value is obtained by spline interpolation from
   * this table.
   *
   * @param[in] matcut    Pointer to the MaterialCuts object in which the full(CSDA) range is requested.
   * @param[in] part      Pointer to the Particle for which the full(CSDA) range is requested.
   * @param[in] kinenergy Kinetic energy of the particle at which the full(CSDA) range is requested.
   * @return    The full(CSDA) range in the given MaterialCuts, Particle, kinetic energy(if the corresponding restricted
   *            range table was built) in internal [length] units.
   *            A high value (1.0e+20) value otherwise.
   */
  double GetRange(const MaterialCuts *matcut, const Particle *part, double kinenergy);


/*
  // just for testing
  void PrintRestrictedDEDX(const MaterialCuts *matcut, const Particle *part) {
    if (fElossTablePerMaterialCuts[matcut->GetIndex()]) {
      fElossTablePerMaterialCuts[matcut->GetIndex()]->PrintRestrictedDEDX(matcut->GetIndex(),part->GetInternalCode());
    } else { // there is no any ELossTable active in the region where the MatrialCut belongs to; should never happen
      std::cerr<<"  ====  ELossTableManager:  No ELossTable for MaterialCuts: \n";
      std::cerr<<matcut<<std::endl;
    }
  }
  // just for testing
  void PrintRestrictedRange(const MaterialCuts *matcut, const Particle *part) {
    if (fElossTablePerMaterialCuts[matcut->GetIndex()]) {
      fElossTablePerMaterialCuts[matcut->GetIndex()]->PrintRestrictedRange(matcut->GetIndex(),part->GetInternalCode());
    } else { // there is no any ELossTable active in the region where the MatrialCut belongs to; should never happen
      std::cerr<<"  ====  ELossTableManager:  No ELossTable for MaterialCuts: \n";
      std::cerr<<matcut<<std::endl;
    }
  }
  // just for testing
  void PrintRange(const MaterialCuts *matcut, const Particle *part) {
    if (fElossTablePerMaterialCuts[matcut->GetIndex()]) {
      fElossTablePerMaterialCuts[matcut->GetIndex()]->PrintRange(matcut->GetMaterial()->GetIndex(),part->GetInternalCode());
    } else { // there is no any ELossTable active in the region where the MatrialCut belongs to; should never happen
      std::cerr<<"  ====  ELossTableManager:  No ELossTable for MaterialCuts: \n";
      std::cerr<<matcut<<std::endl;
    }
  }
*/


private:
  // CTR private
  ELossTableManager() {}


private:
  std::vector<ELossTable*>      fElossTablePerMaterialCuts;  // size will be #MaterialCuts and each element will
                                                             // point to the ELossTable object that that MatrialCuts
                                                             // belongs to; this should not be deleted because stored
                                                             // below
  std::vector<ELossTable*>      fTheElossTables;  // as many ELossTables object as registered PhysicsParameters; the
                                                  // class owns these obejcts

};

}  // namespace geantphysics

#endif  // ELOSSTABLEMANAGER_H
