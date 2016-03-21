//==========================================================================
// PhysicsList.h - Geant-V Prototype
/**
 * @file PhysicsList.h
 * @brief Physics List, i.e. the list of physics processes assigned to
 *        each particle in the Geant-V prototype
 *
 * This is the base class from which any physics list must derived from.
 * We assume that all the GeantV particles are present in the simulation,
 * and those which do not have any physics process assigned to them 
 * explicitly in the physics list will be transported without physics
 * interactions (like geantinos).
 * This class is responsible to do the following at initialization:
 * - construct all the physics processes;
 * - construct all the PhysicsManagerPerParticle objects (associated to each
 *   of the particle types we want to include in the simulation);
 * - create all needed physics tables. 
 *
 * @author Mihaly Novak & Alberto Ribon (Nov 2015)
 */
//==========================================================================

#ifndef PHYSICS_LIST
#define PHYSICS_LIST


namespace geant {

/**
 * @brief Class PhysicsList
 */
class PhysicsList {
private:
  // EMPTY FOR THE TIME BEING!
 
public:
  /** @brief PhysicsList default constructor */
  PhysicsList();

  /** @brief PhysicsList copy constructor */
  // Not clear if needed: if not it should be put in the private part of the class.
  PhysicsList( const PhysicsList &other );

  /** @brief Operator = */
  // Not clear if needed: if not it should be put in the private part of the class.
  PhysicsList& operator=( const PhysicsList &other );

  /** @brief PhysicsList destructor */
  virtual ~PhysicsList();

  /** @brief Main method, called at the initialization, for the physics list.
   *
   *  The method should do the following:
   *  - construct all the physics processes;
   *  - construct all the PhysicsManagerPerParticle objects (associated to each
   *    of the particle types we want to include in the simulation);
   *  - create all needed physics tables. 
   *  The method is called by GeantPhysicsInterface::Initialize
   *  (which is called by GeantPropagator::Initialize ,
   *   which is called by GeantPropagator::PropagatorGeom )
   */
  virtual void Initialize( /* Not defined yet */ );

};

}  // end of namespace geant

#endif
