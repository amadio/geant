

#ifndef ELOSSTABLEREGISTER_H
#define ELOSSTABLEREGISTER_H

#include <vector>

namespace geantphysics {

class EMPhysicsProcess;
/**
 * @brief   Singletone used only at initialization to register kEnergyLoss EMPhysicsProcess for particles.
 * @class   ELossTableRegister
 * @author  M Novak, A Ribon
 * @date    august 2016
 *
 * Each kEnergyLoss EMPhysicsProcess will register itself at initialization in this ELossTableRegister for the
 * Particle(s) the given process is assigned to. This information will be used by the ELossTable object to build all
 * energy loss related data tables for a given Particle (in all MaterialCuts that the given ELossTable needs to handle).
 *
 * The ELossTableRegister::Clear() method, that cleans up the register, must be called at re-initialization.
 */

class ELossTableRegister {
public:
  /**
   * @brief Static method to obtain the singletone object.
   */
  static ELossTableRegister& Instance();


  // copy CTR and assignment operators as deleted
  ELossTableRegister(const ELossTableRegister&) = delete;
  ELossTableRegister& operator=(const ELossTableRegister&) = delete;

  /**
   * @brief Public method to register kEnergyLoss EMPhysicsProcess for a given Particle.
   *
   * At the initialization of each kEnergyLoss EMPhysicsProcess, the EMPhysicsProcess will be automatically registered
   * for the Particle it is assigned to by calling this method.
   *
   * @param[in] internalpartindx  Internal index of the Particle to which this kEnergyLoss EMPhysicsProcess is assigned to.
   * @param[in] elossprocess      Pointer to the kEnergyLoss EMPhysicsProcess that must be registered for the Particle.
   */
  void RegisterEnergyLossProcessForParticle(int internalpartindx, EMPhysicsProcess *elossprocess);


  /**
   * @brief Public method to obtain the list of kEnergyLoss EMPhysicsProcess-es registered for a given Particle.
   *
   * @param[in] internalpartindx Internal index of the Particle for which the list of kEnergyLoss EMPhysicsProcess is
   *                             requested
   * @return    List of pointers to kEnergyLoss EMPhysicsProcess objects that has been registered (assigned) to the
   *            given Particle.
   */
  const std::vector<EMPhysicsProcess*>& GetListEnergyLossProcessesPerParticle(int internalpartindx) const {
    return fELossProcessesPerParticle[internalpartindx];
  }


  /**
   * @brief Public method to obtain the list of kEnergyLoss EMPhysicsProcess-es registered for all Particles.
   *
   * @return    Vector of lists per internal Particle indices of pointers to kEnergyLoss EMPhysicsProcess objects that
   *            has been registered (assigned) to each Particles. Size of the vector is the maximum internal index + 1
   *            of the Particle-s that has any kEnergyLoss EMPhysicsProcess registered.
   */
  const std::vector<std::vector<EMPhysicsProcess*> >& GetListEnergyLossProcesses() const {
    return fELossProcessesPerParticle;
  }

  /**
   * @brief Method to clean up the register. Must be called before at each re-initialization before the PhysicsProcess-es
   *        are initialised to clean the register.
   */
  void Clear();

//
// private methods
//
private:
  // CTR is private.
  ELossTableRegister() {}

//
// private data members
//
private:
  // the i-th element contains a vector of registered kEnergyLoss EMPhysicsProcess-es for
  // particle with internal code of i; the size will be set to the maximum of the particles internal index + 1 that has
  // any kEnergyLoss EMPhysicsProcess registered; the class do NOT own the EMPhysicsProcess-es
  std::vector<std::vector<EMPhysicsProcess*> >  fELossProcessesPerParticle;
};


}  // namespace geantphysics

#endif // ELOSSTABLEREGISTER_H
