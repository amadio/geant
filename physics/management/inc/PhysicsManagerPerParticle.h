//==========================================================================
// PhysicsList.h - Geant-V Prototype
/**
 * @file PhysicsManagerPerParticle.h
 * @brief This class keeps the list of physics processes associated to one
 *        type of particle in the Geant-V prototype
 *
 * For each type of particle (electron, positron, gamma, proton, etc.) there
 * should be one, unique object of this class that keeps the list of physics
 * processes associated to that type of particle.
 * It is foreseen that the physics list instantiates at initialization
 * all the PhysicsManagerPerParticle objects, corresponding to each
 * particle type that has at least one physics process, and fills them
 * with the corresponding physics processes defined in the physics list.
 *
 * @author Mihaly Novak & Alberto Ribon (Dec 2015)
 */
//==========================================================================

#ifndef PHYSICS_MANAGER_PER_PARTICLE
#define PHYSICS_MANAGER_PER_PARTICLE


#include <string>
#include <vector>


namespace geant {

// Forward declarations
class PhysicsProcess;


/**
 * @brief Class PhysicsManagerPerParticle
 */
class PhysicsManagerPerParticle {
private:
  int fParticleGVcode;  /** GV code of the particle associated to this physics manager */

  /** The same process can be inserted in none, one or more of the following lists: */
  std::vector< PhysicsProcess* > fAlongStepProcessVec;          /** List of along-step processes */
  std::vector< PhysicsProcess* > fPostStepCandidateProcessVec;  /** List of post-step candidate processes */
  std::vector< PhysicsProcess* > fPostStepForcedProcessVec;     /** List of post-step forced processes */
  std::vector< PhysicsProcess* > fAtRestCandidateProcessVec;    /** List of at-rest candidate processes */
  std::vector< PhysicsProcess* > fAtRestForcedProcessVec;       /** List of at-rest forced processes */

  // POINTER TO THE TOTAL LAMBDA TABLE (FOR EACH MATERIAL-CUT COUPLE) IS NEEDED HERE!

  /** @brief PhysicsManagerPerParticle copy constructor is not defined */
  PhysicsManagerPerParticle( const PhysicsManagerPerParticle &other );

  /** @brief Operator = is not defined */
  PhysicsManagerPerParticle& operator=( const PhysicsManagerPerParticle &other );

public:
  /** @brief PhysicsManagerPerParticle constructors 
   *  
   *  It is expected that one single object of this class is instantiated
   *  for each particle type that has at least one active physics process
   *  associated to it.
   *  It is the physics list, at initialization, that is expected to call
   *  the constructor(s) of this class.
   */
  PhysicsManagerPerParticle();
  explicit PhysicsManagerPerParticle( const int aParticleGVcode );

  /** @brief PhysicsManagerPerParticle destructor */
  ~PhysicsManagerPerParticle();

  /** @brief Add a physics process to the list of processes handled by this manager 
   *  @param aPhysicsProcess a pointer to a physics process object
   *
   *  The method is responsible for checking that the process is defined for
   *  the same particle type expected for this manager, and if it is active:
   *  only an active process associated to the correct particle will be inserted.
   *  The lists, if any, in which the physics process is added to depend on the
   *  properties of the process.
   *  Note that a "forced" process is always considered also as a "candidate" one.
   *
   *  This method should be called only at initialization by the physics list.
   *
   */
  void AddProcess( PhysicsProcess *process );

  /** @brief Method that is called once only at initialization. 
   *  SIGNATURE (INPUT PARAMETERS AND RETURN TYPE) TO BE REDEFINED LATER...
   *
   *  THE TOTAL LAMBDA TABLES (FOR EACH MATERIAL-CUT COUPLE) ARE OBTAINED
   *  AS THE INVERSE OF THE SUM OF THE INVERSE OF THE LAMBDA TABLES FOR
   *  EACH DISCRETE PROCESS ASSOCIATED TO THE PARTICLE.  
   *
   *  This method is expected to be called at initialization by the physics list.
   */
  void BuildTotalLambdaTables( /* Not defined yet */ );
 
  /** @brief Method that returns the total lambda table.
   *  RETURN TYPE AND EVENTUAL INPUT PARAMETERS TO BE REDEFINED LATER...
   */
  void GetTotalLambdaTable( /* Not defined yet */ ) const;

  /** @brief Method that returns the internal GeantV code of the particle type
             associated to this physics manager */
  int GetParticleType() const;

  /** @brief Method that sets the particle type associated to this physics manager 
   *  @param aParticleGVcode internal GeantV code of the particle
   */
  void SetParticleType( const int aParticleGVcode ) {
    fParticleGVcode = aParticleGVcode;
  }

  /** @brief Method that return the list of along-step processes */
  const std::vector< PhysicsProcess* >& GetListAlongStepProcesses() const {
    return fAlongStepProcessVec;
  }

  /** @brief Method that return the list of candidate discrete processes */
  const std::vector< PhysicsProcess* >& GetListPostStepCandidateProcesses() const {
    return fPostStepCandidateProcessVec;
  }

  /** @brief Method that return the list of forced discrete processes */  
  const std::vector< PhysicsProcess* >& GetListPostStepForcedProcesses() const {
    return fPostStepForcedProcessVec;
  }

  /** @brief Method that return the list of candidate at-rest processes */
  const std::vector< PhysicsProcess* >& GetListAtRestCandidateProcesses() const {
    return fAtRestCandidateProcessVec;
  }

  /** @brief Method that return the list of forced at-rest processes */
  const std::vector< PhysicsProcess* >& GetListAtRestForcedProcesses() const {
    return fAtRestForcedProcessVec;
  }

};

}  // end of namespace geant

#endif
