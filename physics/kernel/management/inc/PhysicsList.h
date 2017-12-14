
#ifndef PHYSICS_LIST
#define PHYSICS_LIST

#include <string>

namespace geantphysics {

class Particle;
class PhysicsProcess;
class PhysicsParameters;
/**
 * @brief  Physics List, i.e. the list of physics processes assigned to each particle in the Geant-V prototype.
 * @class  PhysicsList
 * @author M Novak, A Ribon
 * @date   november 2015
 *
 * This is the base class from which any physics list must derived from.  We assume that all the GeantV particles are
 * present in the simulation, and those which do not have any physics process assigned to them explicitly in the physics
 * list will be transported without physics interactions (like geantinos). User can assigne PhysicsProcess-es to
 * particles by using this class. The PhysicsListManager will invoke all registered PhysicsList Initialize method.
 */
class PhysicsList {
public:
  /** @brief PhysicsList constructor.
   *
   * @param[in]  name  Name of the physics list.
   */
  PhysicsList(const std::string &name);

  /** @brief PhysicsList destructor */
  virtual ~PhysicsList();

  /** @brief Main method, called at the initialization, for the physics list by the PhysicsListManager.
   *
   *  The method should do the following:
   *  - construct all the physics processes for all particles and assigne the processes to partciles;
   *  - the PhysicsParameters object will be available in this method so the user can change any values of it
   */
  virtual void Initialize( /* Not defined yet */ ) = 0;

  /** @brief Get the name of this physics list.
   *
   *  @return Name of this phyics list.
   */
  const std::string& GetName() const {return fName;}

  /** @brief Get the PhysicsParameters member pointer of this physics list.
    *
    * The physics list owns the physics parameter object;
    *
    * @return Pointer to the physics parameter object of this physics list.
    */
  PhysicsParameters* GetPhysicsParameters() { return fPhysicsParameters; }

protected:
   /** @brief Insert a given process object pointer to the particle process list.
    *
    * @param[in] particle Pointer to the particle object to which the process should be assigned.
    * @param[in] process  Pointer to a process obejct that should be assigned to the particle.
    */
   void AddProcessToParticle(Particle *particle, PhysicsProcess *process);

//
// data members
//
private:
  std::string         fName;               /** Name of this physics list. */

protected:
  PhysicsParameters  *fPhysicsParameters;  /** the PhysicsParameters object of this physics list;
                                               will be created in the PhysicsList base ctr;
                                               the class do NOT own the object because they are owned by the
                                               PhysicsParameters itself and all PhysicsParameters can be removed by
                                               PhysicsParameters::Clear() */
};

}  // namespace geantphysics

#endif  // PHYSICS_LIST
