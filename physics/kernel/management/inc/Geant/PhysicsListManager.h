
#ifndef PHYSICSLISTMANAGER_H
#define PHYSICSLISTMANAGER_H

#include <vector>

namespace geantphysics {
  class PhysicsList;
  class PhysicsManagerPerParticle;
  class PhysicsParameters;
/**
 * @brief  One of the main physics manager classes(singletone).
 * @class  PhysicsListManager
 * @author M Novak, A Ribon
 * @date   july 2016
 *
 * The main method of the class is the BuildPhysicsLists() that is responsible for doing the followings:
 *    - creates all particles
 *    - loops over all PhysicsList registered by the user in PhysicsListManager and calls one-by-one the
 *      Initialize() method of the registered PhysicsLists:
 *      -  each call to a PhysicsList Initialize() method will assigne processes to the temporary PhysicsProcess vector
 *         of the static particle object
 *      -  create one PhysicsManagerPerParticle object for each partcile that the current PhysicsList has added at
 *         least one PhysicsProcess
 *      -  add all processes from the current physics list to this PhysicsManagerPerParticle
 *         and push the particle to the particle list of the phsyics process that the physics process is assigned to
 *      -  also set the active region indices both in the PhysicsProcess-es(a) and in the PhysicsManagerPerParticle(b)
 *         objects
 *      -  store the created PhysicsManagerPerParticle object pointer in a table (will be used to delete and initilize
 *         only once them)
 *      -  set pointers to regional PhysicsManagerPerParticle in the static Particle definition to point to this
 *         object where the current physics list is active
 *    -  each created PhysicsManagerPerParticle is initilised by calling their Initialize() method:
 *      -  each processes are initialised by calling their Initialize() method:
 *         -  it is checked in the base PhysicsProcess Initialize() method if the process is assigned only to allowed
 *            particles
 *         -  in case of EMPhysicsProcess-es (after calling the base i.e. PhysicsProcess Initialize() method)
 *            the EMModelManager member of the EMPhysicsProcess is initialized that will initialize the EMModel-s
 *            as well togeter with setting reagions where they active. The default active regions for EMModel-s are
 *            determined by the active regions of the EMPhysicsProcess(PhysicsProcess) that they belong to. On the top
 *            of this, user requested inactive regions in case of the individual EMModel-s are considered by the
 *            EMModelManager when it sets the active regions of the individual EMModel-s. At the end of initialisation
 *            of the EMModelManager member of the EMPhysicsProcess: EMModel-s per region ponters are set and the
 *            corresponding EMModel-s are initialised and each EMModel knows its list of active regions.
 *         -  if the base EMModel class InitialiseElementSelectors() method was explicitly called from the derived
 *            emmodel class Initialize() method (at the end i.e. after the model is properly initialised), then
 *            target element selector is built for the given model and tagget elemnt can be sampled at run time by the
 *            int EMModel::SampleTargetElementIndex(const MaterialCuts* matcut, double ekin, double rndm) method.
 *         -  PhysicsProcess-es set to be kEnergyLoss processes are registered in the ELossTableRegister for the given
 *            particle.
 *      -   lambda tables (if any, see more at PhysicsManagerPerParticle::BuildLambdaTables()) i.e. total(a) and per-
 *          discrete process(b) macroscopic cross section tables are built over an kinetic energy grid and set up for
 *           -  run time total discrete step limit(a) (i.e. total inverse lambda)
 *              double PhysicsManagerPerParticle::GetInvTotalLambda(const MaterialCuts *matcut, double kinenergy)
 *           -  discrete process sampling(b)
 *              const PhysicsProcess* PhysicsManagerPerParticle:::SelectDiscreteInteraction(const MaterialCuts *matcut,
 *                                                                double kinenergy, double rndm)
 *    - the ELossTableManager is initialised i.e. for all partciles that at least one PhysicsProcess::kEnergyLoss
 *      process had been registered in the ELossTableRegister: dedx, range, inverse range tables are built over an
 *      energy grid per material-cuts that belongs to region where the given process is active and run time
 *      interpolation is set up and can be used at run time to get:
 *       - restricted stopping power:
 *            double ELossTableManager::GetRestrictedDEDX(const MaterialCuts *matcut, const Particle *part,
 *                                                        double kinenergy)
 *       - restricted range:
 *            double ELossTableManager::GetRestrictedRange(const MaterialCuts *matcut, const Particle *part,
 *                                                         double kinenergy)
 *       - energy that corresponds to a given restricted range i.e. from the inverse range table:
 *            double ELossTableManager::GetEnergyForRestrictedRange(const MaterialCuts *matcut, const Particle *part,
 *                                                                  xdouble range)
 */
class PhysicsListManager {
public:
  static PhysicsListManager& Instance();

  // copy CTR and assignment operators as deleted
  PhysicsListManager(const PhysicsListManager&) = delete;
  PhysicsListManager& operator=(const PhysicsListManager&) = delete;

  // must be set before we start to use i.e. before we call the BuildPhysicsLists methods!!!!
  void SetNumberOfRegions(const int numregions) { fNumOfRegions=numregions; }
  int  GetNumberOfRegions() const { return fNumOfRegions; }

  void RegisterPhysicsList(PhysicsList *physlist, std::vector<bool> activeregionlist);
  void RegisterPhysicsList(PhysicsList *physlist);

  int  GetNumberOfRegisteredPhysicsLists() const { return fPhysicsListVector.size(); } 

  void BuildPhysicsLists();

  // at the end we can clear all physics processes and clear the physics process vector and process manager per particle
  // vector of the static Particle properties
  // we also delete the physics lists added by the user and clear the local physics lists vector and active indices vector
  void ClearAll();

  // this is just for testing
  void PrintAll();


private:
  // CTR
  PhysicsListManager(){}
  // create all particles
  void CreateAllParticles();

private:
  int                                       fNumOfRegions;  // number of regions

  // stores the registered physics lists pointers; size will be at the end as many as physics lists registered by
  // the user; the object owns these PhysicsList obejcts and they are cleand in the ClearAll method
  std::vector<PhysicsList*>                  fPhysicsListVector;
  // for each registered physics list it stores a mask for the active regions; size is #regions times #physics lists
  std::vector<std::vector<bool> >            fActiveRegionMasks;
  // a unique store of pointers to all physics manager per particle objects that has been created; the object owns these
  // PhysicsManagerPerParticle objects and they are cleand in the ClearAll method
  std::vector<PhysicsManagerPerParticle*>    fPhysicsManagerPerParticleTable;
};

} // namespace geantphysics

#endif // PHYSICSLISTMANAGER_H
