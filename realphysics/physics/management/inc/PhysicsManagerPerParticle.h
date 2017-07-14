
#ifndef PHYSICS_MANAGER_PER_PARTICLE
#define PHYSICS_MANAGER_PER_PARTICLE

// for inlice namespace GEANT_IMPL_NAMESPACE
#include "Geant/Config.h"


#include <string>
#include <vector>

#include "GeantTaskData.h"
#include "GeantTrack.h"

namespace geantphysics {
  inline namespace GEANT_IMPL_NAMESPACE {
    class Material;
  }
}


namespace geantphysics {

 // forward declarations
class PhysicsParameters;
//class Material;
class MaterialCuts;
class PhysicsProcess;
class Particle;
class LightTrack;

/**
 * @brief  This class keeps the list of physics processes associated to one type of particle.
 * @class  PhysicsManagerPerParticle
 * @author M Novak, A Ribon
 * @date   december 2015
 *
 * For each type of particle (electron, positron, gamma, proton, etc.) there should be one, unique object of this class
 * in each set of active regions (i.e. per physics list) that keeps the list of physics processes associated to that
 * type of particle. The PhysicsListManager instantiates at initialization all the PhysicsManagerPerParticle objects,
 * corresponding to each particle type that has at least one physics process, and fills them with the corresponding
 * physics processes defined in the physics list.
 */
class PhysicsManagerPerParticle {
public:

  /** @brief PhysicsManagerPerParticle constructors
   *
   *  It is expected that one single object of this class is instantiated
   *  for each particle type that has at least one active physics process
   *  associated to it.
   *  It is the physics list manager, at initialization, that is expected to call
   *  the constructor of this class.
   *
   *  @param[in] particle  Pointer to the particle object this manager belongs to.
   *  @param[in] physpars  Pointer to the PhysicsParameters object that belongs to the same set of regions that this
   *                       manager.
   */
  PhysicsManagerPerParticle(const Particle *particle, const PhysicsParameters *physpars);

  /** @brief PhysicsManagerPerParticle destructor */
  ~PhysicsManagerPerParticle();

  /**
   *   Initialize all processes assigned to the particle that this manager belongs to.
   *
   */
  void Initialize();

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
  void AddProcess(PhysicsProcess *process);

  /** @brief Method that return the list of all processes */
  const std::vector<PhysicsProcess*>& GetListProcesses() const {
    return fProcessVec;
  }

  /** @brief Method that return the list of along-step processes */
  const std::vector<PhysicsProcess*>& GetListAlongStepProcesses() const {
    return fAlongStepProcessVec;
  }

  /** @brief Method that return the list of candidate discrete processes */
  const std::vector<PhysicsProcess*>& GetListPostStepCandidateProcesses() const {
    return fPostStepCandidateProcessVec;
  }

  /** @brief Method that return the list of forced discrete processes */
  const std::vector<PhysicsProcess*>& GetListPostStepForcedProcesses() const {
    return fPostStepForcedProcessVec;
  }

  /** @brief Method that return the list of candidate at-rest processes */
  const std::vector<PhysicsProcess*>& GetListAtRestCandidateProcesses() const {
    return fAtRestCandidateProcessVec;
  }

  /** @brief Method that return the list of forced at-rest processes */
  const std::vector<PhysicsProcess*>& GetListAtRestForcedProcesses() const {
    return fAtRestForcedProcessVec;
  }

  /** Methods that are needed for the new concept of physics-per-region:
      is this physics manager per particle active in the i-th region? */
  std::vector< bool >& GetListActiveRegions() { return fListActiveRegions; }
  bool IsActiveRegion( const int regionindx ) const { return fListActiveRegions[ regionindx ]; }


  void PrepareForRun();

  void ComputeIntLen(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td);
  int  AlongStepAction(LightTrack &track, Geant::GeantTaskData *td);
  int  PostStepAction(LightTrack &track, Geant::GeantTrack *gtrack, Geant::GeantTaskData *td);

  bool  HasEnergyLossProcess() const { return fIsHasElossProcess; }
  bool  HasMSCProcess() const { return fIsHasMSCProcess; }
  const PhysicsProcess* GetMSCProcess() const;

// private methods
//
private:
  /** @brief PhysicsManagerPerParticle copy constructor is not defined */
  PhysicsManagerPerParticle(const PhysicsManagerPerParticle &other);
  /** @brief Operator = is not defined */
  PhysicsManagerPerParticle& operator=(const PhysicsManagerPerParticle &other);
  // checks the processes assigned to the particle and sets some flags
  void CheckProcesses();



private:
  /** Pointer to the partcile object to which this physics manager belongs to.
   *  The class doesn't own the partcile object.
   */
 const Particle            *fParticle;           // not owned;
 const PhysicsParameters   *fPhysicsParameters;  // not owned;

 // some flags to indicate ...
 bool        fIsHasMSCProcess;
 bool        fIsHasDecayProcess;
 bool        fIsHasElossProcess;

 /** The same process can be inserted in none, one or more of the following lists:
  *  The class doesn't own any PhysicsProcess objects! Each PhysicsProcess object pointer that are created stored in
  *  global table and the corresponding objects are deleted by the main manager calling PhysicsProcess::ClearAllProcess
  */
 std::vector<PhysicsProcess*> fProcessVec;                   /** List of processes, contains all process pointers that has been added */
 std::vector<PhysicsProcess*> fAlongStepProcessVec;          /** List of along-step processes */
 std::vector<PhysicsProcess*> fPostStepCandidateProcessVec;  /** List of post-step candidate processes */
 std::vector<PhysicsProcess*> fPostStepForcedProcessVec;     /** List of post-step forced processes */
 std::vector<PhysicsProcess*> fAtRestCandidateProcessVec;    /** List of at-rest candidate processes */
 std::vector<PhysicsProcess*> fAtRestForcedProcessVec;       /** List of at-rest forced processes */

 std::vector<bool> fListActiveRegions;  /** is this physics manager per particle active in the i-th region? */

};

}  // end of namespace geantphysics

#endif
