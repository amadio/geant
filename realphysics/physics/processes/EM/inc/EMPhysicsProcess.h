
#ifndef EMPHYSICSPROCESS_H
#define EMPHYSICSPROCESS_H

#include "PhysicsProcess.h"
#include <vector>

namespace geantphysics {

// forward declarations
class MaterialCuts;
class Partcile;
class EMModel;
class EMModelManager;

class LightTrack;

/**
 * @brief   Base class to provide interface for all ordinary electromagnetic physics processes.
 * @class   EMPhysicsProcess
 * @author  M Novak, A Ribon
 * @date    august 2016
 *
 * The interface provides common implementations for some methods: to initialize the EMPhysicsProcess, compute
 * macroscopic cross section, dE/dx (for kEnergyLoss processes) and to add EMModel-s. Each EMPhysicsProcess has an
 * EMModelManager member. EMModels can be added and accessed through this EMModelManager member of the EMPhysicsProcess.
 * The EMModelManager object owns all the models that have been added to the EMPhysicsProcess.
 *
 * TODO: implement method for AlongStepLimit, AlongStepDoIt, (post step limit is handled by the PhysicsManagerPerParticle)
 * PostStepDoIt.
 */

class EMPhysicsProcess : public PhysicsProcess {
public:
  /**
   * @brief CTR.
   */
  EMPhysicsProcess(const std::string &name);


  /**
   * @brief DTR.
   */
  virtual ~EMPhysicsProcess();


  /**
   * @brief  Virtual method from the PhysicsProcess base class to handle initialization of this EMPhysicsProcess.
   *
   * This method must be called for each electromagnetic physics processes: those derived processes that has their own
   * Initialize() method implementation needs to call this base class Initialize() method explicitely!
   *
   * Fisrt the PhysicsProcess base class Initialize method will be called that: it will be checked there if the process
   * is assigned only to particles that the process is allowed to. Then the Initialize method of the EMModelManager
   * member of this EMPhysicsProcess will be called: (1) will set for each EMModel-s that have been added to the
   * EMPhysicsProcess the list of active regions: (a) first the default active regions which are determined by the
   * active regions of the EMPhysicsProcess; (b) on the top of this, the user requested inactive regions are considered
   * in case of each individual EMModel-s; (2) the EMModelManager will prepare the list of active EMModel-s per region
   * collections. (3) After each EMModel-s active regions list is set, the EMModelManager will initilise all EMModel-s:
   * each EMModel-s will be initilized only for reagions where they are active.
   * After the initialization of the EMModelManager member, it is checked if this EMPhysicsProcess is a kEnergyLoss
   * process and if yes then it will be registered in the ELossTableRegister singletone for particles this kEnergyLoss
   * EMPhysicsProcess is assigned to. This information will be used later to build all energy loss related tables like
   * dEdx or range table by the ELossTable objects that are created and handled by the ELossTableManager singletone.
   */
  virtual void   Initialize();

  /**
   * @brief Common method for EMPhysicsProcess to compute stopping power for a given MaterialCuts, Partcile, kinetic
   *        energy.
   *
   * This method is expected to be called at initialization i.e. from the ELossTable object when the energy loss
   * related tables are built.
   *
   * EMPhysicsProcess-s that have enegy loss along the step must set to be kEnergyLoss process. The corresponding
   * EMModel-s must implement the base EMModel::ComputeDEDX() method. Each kEnergyLoss EMPhysicsProcess is registered
   * automatically in the ELossTableRegister for the Particle(s) the process is assigned to. This information is used
   * later when the ELossTable objects (to build, store and handle energy loss related tables) are set up by the
   * ELossTableManager singletone.
   *
   * @param[in] matcut      Pointer to the MaterialCuts object in which the dEdx must be computed.
   * @param[in] kinenergy   Kinetic energy of the Partcile at which the dEdx must be computed.
   * @param[in] particle    Pointer to the Partcile object for which the dEdx must be computed.
   * @param[in] istotal     Flag to indicate if full dEdx must be computed. False by default i.e. restricted dEdx is
   *                        required.
   * @return    Restricted or full stopping power computed by one of the EMModel models that had been added to this
   *            EMPhysicsProcess. One EMModel is selected based on the kinetic energy of the particle, the low/high
   *            energy usage limits of the EMModel-s (added to this EMPhysicsProcess) that are active in the region to
   *            which the given MaterialCuts beongs to. The restricted stopping power is provided in internal [energy/
   *            length] units.
   */
  double ComputeDEDX(const MaterialCuts *matcut, double kinenergy, const Particle *particle, bool istotal=false);


  /**
   * @brief Common implementation of the base PhysicsProcess::ComputeMacroscopicXSection() method for each
   *        EMPhysicsProcess.
   *
   * This method is expected to be called at initialization from the PhysicsManagerPerParticle object when the lambda
   * (i.e. macroscopic cross section tables are built) for those EMPhysicsProcess-es that have discrete part
   * (fIsDiscrete member was set to be true).
   *
   * Each EMModel, that describes an interaction that can happen at the post step point, must implement the
   * corresponding base EMModel::ComputeMacroscopicXSection() method.
   *
   * @param[in] matcut      Pointer to the MaterialCuts object in which the x-sec must be computed.
   * @param[in] kinenergy   Kinetic energy of the Partcile at which the x-sec must be computed.
   * @param[in] particle    Pointer to the Partcile object for which the x-sec must be computed.
   * @param[in] mass        Dynamic mas of the particle.
   * @return    Macroscopic cross section computed by one of the EMModel models that had been added to this
   *            EMPhysicsProcess. One EMModel is selected based on the kinetic energy of the particle, the low/high
   *            energy usage limits of the EMModel-s (added to this EMPhysicsProcess) that are active in the region to
   *            which the given MaterialCuts beongs to. The macroscopic cross section is provided in internal [1/length]
   *            units.
   */
  virtual double ComputeMacroscopicXSection(const MaterialCuts *matcut, double kinenergy, const Particle *particle, double mass) const;

  // used to determine the lowest energy of the lambda table when lambda table is requested to build by a given process
  virtual double GetMinimumLambdaTableKineticEnergy(const MaterialCuts *matcut, const Particle*) const;

  /**
   * @brief Common implementation of the continuous step limit method of the base PhysicsProcess class for ordinary
   *        electromagnetic processes.
   *
   * Only kEnergyLoss processes will limit the step. A step limit function is used with 2 parameters set in the
   * PhysicsParameters object. In case of particles that have any kEnergyLoss processes: this will be the cumulative
   * along step limit from all kEnergyLoss processes i.e.: should make sure that only one kEnergyLoss process stays in
   * the continuous process vector because we use the cumulative energy loss related tables.
   *
   * @return    Continuous step limit from this EMPhysicsProcess. In case of particles that have any kEnergyLoss
   *            processes: this will be the cumulative along step limit from all kEnergyLoss processes i.e.: should make
   *            sure that only one kEnergyLoss process stays in the continuous process vector because we use the
   *            cumulative energy loss related tables.
   */
   virtual double AlongStepLimitationLength(Geant::GeantTrack* /*gtrack*/, Geant::GeantTaskData* /*td*/) const;

  /**
   * @brief Common implementation of the AlongStepDoIt method of the base PhysicsProcess class for ordinary
   *        electromagnetic processes.
   *
   * Only kEnergyLoss processes will do something: along step energy loss will be computed.  This will be the cumulative
   * along step limit from all kEnergyLoss processes i.e.: should make sure that only one kEnergyLoss process stays in
   * the continuous process vector because we use the cumulative energy loss related tables.
   *
   * @param[in/out]  track     A LightTrack object that contains all the necessary primary track information.
   * @param[in/out]  sectracks List of secondary tracks created in this DoIt method.
   * @return    Number of secondary tracks created and stored in the sectracks vector.
   */
  virtual  int AlongStepDoIt(LightTrack &track, Geant::GeantTaskData *td);

  // for msc: no secondaries and acts directly on GeantTrack
  virtual  void AlongStepDoIt(Geant::GeantTrack* /*gtrack*/, Geant::GeantTaskData* /*td*/) const {}


  // Will be called only if disceret interaction was selected
  virtual  int PostStepDoIt(LightTrack &track, Geant::GeantTaskData *td);


  /**
   * @brief Method to add EMModel to the EMPhysicsProcess.
   *
   * Each EMPhysicsProcess has an EMModelManager member that handels the EMModel-s of the process. This member will be
   * used to store the EMModel-s.
   *
   * @param[in] mod  Pointer to an EMModel object to be added to the physics process.
   * @return    The index of the EMModel pointer in the EMModelManager container.
   */
  int AddModel(EMModel *model);


  /**
   * @brief Method to get a pointer to the EMModelManager of this EMPhysicsProcess.
   *
   * Each EMPhysicsProcess has an EMModelManager member that handels the EMModel-s of the process.
   * NOTE: this method has been added only for testing reagions and not used during a normal simulation.
   *
   * @return    Pointer to the EMModelManager of this EMPhysicsProcess.
   */
  EMModelManager* GetModelManager() const { return fModelManager; }


  /**
   * @brief Method to print out some info. Must be change later to a DumpProcess method. */
  friend std::ostream& operator<<(std::ostream& flux,  EMPhysicsProcess* proc);

//
// data members
//
private:
  EMModelManager    *fModelManager; /** each EMPhysicsProcess has an EMModelManager member to handle the EMModel-s of
                                        the process. */
  // continuous step limit related parameters; set only for kEnergyLoss processes at initialization
  double fFinalRange;
  double fDRoverRange;

  double fLowestKineticEnergy;   // this is a kind of tracking cut for charged partciles.
  double fLinearEnergyLossLimit; // fraction of initial kinetic energy allowed to be lost under linear approximation

};

} // namespace geantphysics

#endif // EMPHYSICSPROCESS_H
