//==========================================================================
// PhysicsList.h - Geant-V Prototype
/**
 * @file PhysicsProcess.h
 * @brief Base class for all single physics processes in the Geant-V prototype.
 *
 * Similarly to the Geant4 G4VProcess concept, this class specifies
 * one single physics process for one particular particle type.
 * This is the base class from which all single physics processes, either
 * in-flight and/or at-rest, must derived from.
 * As in Geant4, it has two main types of methods, "GetLength" (called
 * here differently but with the same property) and "DoIt", and (at most)
 * three main "components": along-the-step (continuous), post-step (discrete)
 * and at-rest.
 *
 * @author Mihaly Novak & Alberto Ribon (Nov 2015)
 */
//==========================================================================

#ifndef PHYSICS_PROCESS
#define PHYSICS_PROCESS


// REMINDER-TO-BE-REMOVED: THE NAME OF THIS CLASS, PhysicsProcess, REQUIRES
// ----------------------  THAT THE CURRENT KERNEL CLASS WITH THE SAME NAME
//                         BE RENAMED AS GeantPhysicsInterface .
//                         IF NOT, THEN THE CLASS IN THIS FILE SHOULD BE
//                         RENAMED AS SinglePhysicsProcess .


#include <string>


namespace geant {

// Forward declarations
class LightTrack;
class LightTrack_v;


enum class ForcedCondition { 
             kInActivated,       /** for inactivated processes */
             kNotForced,         /** for normal, common processes: default */
             kForced,            /** for special processes, which should be always invoked */
             kExclusivelyForced  /** for fast-sim processes, which not only should be always
                                   * invoked, but should prevent other processes to be called */
           };


enum class ProcessType {
             kNotDefined,         /** for a process that does not fit in any other type: default */
             kDecay,              /** for a decay process */
             kElectromagnetic,    /** for an electromagnetic process */
             kOptical,            /** for a process involving optical photons */
             kHadronic,           /** for a hadronic process */
             kPhotoLeptonHadron,  /** for a gamma- or lepto-nuclear process */
             kFastSim             /** for a fast-sim (e.g. parameterization) process */
           };


/**
 * @brief Class PhysicsProcess
 */
class PhysicsProcess {
private:
  bool fIsDiscrete;    /** "true" if the process has a discrete part; "false" otherwise */
  bool fIsContinuous;  /** "true" if the process has a continuous part; "false" otherwise */
  bool fIsAtRest;      /** "true" if the process has an at-rest part; "false" otherwise */
  ForcedCondition fForcedCondition;  /** type of ForcedCondition for this process */
  ProcessType fType;                 /** type of this process (MAYBE USEFUL IN THE FUTURE) */
  std::string fName;                 /** name of the process (useful for debugging) */

  // LIKELY POINTERS TO PHYSICS TABLES ARE NEEDED HERE:
  // -  ATOMIC CROSS SECTION TABLES (ONE FOR EACH ELEMENT)
  // -  INVERSE LAMBDA TABLES (ONE FOR EACH MATERIAL-CUT COUPLE)
 
public:
  /** @brief PhysicsProcess default constructor */
  PhysicsProcess();

  /** @brief PhysicsProcess complete constructor */
  PhysicsProcess( const bool aIsDiscrete, const bool aIsContinuous, const bool aIsAtRest,
                  const ForcedCondition aForcedCondition, const ProcessType aType,
                  const std::string &aName );

  /** @brief PhysicsProcess copy constructor */
  PhysicsProcess( const PhysicsProcess &other );

  /** @brief Operator = */
  PhysicsProcess& operator=( const PhysicsProcess &other );

  /** @brief PhysicsProcess destructor */
  virtual ~PhysicsProcess();

public:

  /** @brief Method that is called once only at initialization. 
   *  SIGNATURE (INPUT PARAMETERS AND RETURN TYPE) TO BE REDEFINED LATER...
   *
   *  ATOMIC CROSS SECTION TABLES (ONE FOR EACH ELEMENT), AND
   *  INVERSE LAMBDA TABLE (ONE FOR EACH MATERIAL-CUT COUPLE, ALTHOUGH MOST PROCESSES
   *                        E.G. HADRONIC, HAVE ONLY MATERIAL (NO CUT) DEPENDENCY)
   *  SHOULD BE BUILT.
   *  
   *  This method is expected to be called at initialization by the physics list.
   */
  virtual void BuildPhysicsTables( /* Not yet defined */ ) {}

  /** @brief Method that returns "true" is the particle can have this process,
   *         "false" otherwise.
   *  @param particleGVcode internal code of the particle
   *
   *  Note: this method is meant to be called only at initialization,
   *        and not in the simulation event loop.
   */
  virtual bool IsApplicable( const int particleGVcode ) const {return true;}

  /** @brief Methods that return the atomic (i.e. microscopic) cross section 
   *         (unit: 1/length^2) of the discrete part of this process.
   * 
   *  These methods are meant to be called directly from user code 
   *  (e.g. for testing), and not in the simulation event loop.
   *
   *  The first method corresponds to the minimal signature:
   *  @param projectileCode is the GV particle code of the projectile
   *  @param projectileKineticEnergy is the projectile kinetic energy in MeV
   *  @param targetZ is the atomic number (Z) of the target atom
   *  @param targetN is the number of the nucleons of the target atom
   *         (Note: the default, 0, value means that the average natural abundance
   *                of the element should be considered)
   *
   *  The second method has a LightTrack object as input, and what it does is
   *  simply to extra the information needed by the first method and then calls it.
   */
  virtual double GetAtomicCrossSection( const int projectileCode,
                                        const double projectileKineticEnergy,
                                        const int targetZ, const int targetN = 0 ) const {return 0;}

  double GetAtomicCrossSection( const LightTrack &track ) const;

  /** @brief Method that returns the "macroscopic cross section"
   *         (unit: of 1/length).
   *  
   *  The "macroscopic cross section" is defined as the number of atoms
   *  per unit of volume multiplied by the atomic cross section of the
   *  discrete part of this process.
   *  This method is meant both at initialization - to build the total
   *  lambda table for the particle - and in the simulation event loop -
   *  for the sampling of the type of interaction.
   */
  virtual double InverseLambda( const LightTrack &track ) const;

  /** @brief Method that returns the along-the-step limitation length
   *         (unit: length)
   *
   *  This applies only for the continuous part of this process.
   *  If the process does not have a continuous part, or if it has one
   *  that does not limit the step, then the method returns an arbitrary,
   *  very large value.
   *  The method is meant for the sampling of the type of interactions.
   */
  virtual double AlongStepLimitationLength( const LightTrack &track ) const;
 
  /** @brief Method that returns the average lifetime of this process
   *         (unit: time)
   *
   *  This applies only for the at-rest part of this process.
   *  If the process does not have an at-rest part, or if it has one
   *  that does not limit the time, then the method returns an arbitrary,
   *  very large value.
   *  The method is meant for the sampling of the type of interactions.
   */
  virtual double AverageLifetime( const LightTrack &track ) const;

  /** @brief Method that samples the target atom (Z, N) where the discrete
   *         part of the in-flight process happens.
   *   
   *  The method uses the input parameter, a LightTrack object, also
   *  to store the output of the method, i.e. setting the Z and N
   *  fields of that object.
   */ 
  void SampleTarget( LightTrack &track ) const;


  // REMINDER-TO-BE-REMOVED: THE NEXT 3 METHODS HAVE LightTrack_v AS OUTPUT.
  // ----------------------  THIS COULD BE CHANGED TO:  LightTrack_v*
  //                         OR TO: std::vector< LightTrack > 


  /** @brief Method that does the along-the-step, i.e. the continuous part,
             action of the in-flight process.
   *
   *  The input parameter, a LightTrack object, can be modified: typically
   *  the track status, the local energy deposit and the non-ionizing
   *  local energy deposit can be changed in this method.
   *  The output of the method is a LightTrack_v object, which corresponds
   *  to the new particles created along-the-step (this is rare, but possible:
   *  e.g. for sub-cutoff production; or from atomic de-excitation).
   */ 
  //  virtual LightTrack_v AlongStepDoIt( LightTrack &track ) const;

  /** @brief Method that does the post-step, i.e. the discrete part, action
   *         of the in-flight process.
   *
   *  The input parameter, a LightTrack object, can be modified: typically
   *  the track status is changed in this method.
   *  The output of the method is a LightTrack_v object, which corresponds
   *  to the new particles created by the discrete part of the process.
   */ 
  //  virtual LightTrack_v PostStepDoIt( LightTrack &track ) const;

  /** @brief Method that does the at-rest action of the process.
   *
   *  The input parameter, a LightTrack object, can be modified: typically
   *  the track status is changed in this method.
   *  The output of the method is a LightTrack_v object, which corresponds
   *  to the new particles created by the at-rest part of the process.
   *  Note: this method also includes the sampling of the target atom (Z, N)
   *        where the at-rest process happens.
   */ 
  //  virtual LightTrack_v AtRestDoIt( LightTrack &track ); 

  //--- Getters ---

  /** Method that returns whether this process has a discrete part or not */
  bool GetIsDiscrete() const { return fIsDiscrete; }

  /** Method that returns whether this process has a continuous part or not */
  bool GetIsContinuous() const { return fIsContinuous; }

  /** Method that returns whether this process has an at-rest part or not */
  bool GetIsAtRest() const { return fIsAtRest; }

  /** Method that returns the ForcedCondition type of this process */
  ForcedCondition GetForcedCondition() const { return fForcedCondition; }

  /** Method that returns the type of this process */
  ProcessType GetType() const { return fType; }

  /** Method that returns the name of this process */
  std::string GetName() const { return fName; }

  //--- Setters ---

  /**
   * @brief Method that sets whether this process has a discrete part or not
   * @param aIsDiscrete has the process a discrete part?
   */
  void SetIsDiscrete( const bool aIsDiscrete ) { fIsDiscrete = aIsDiscrete; }

  /**
   * @brief Method that sets whether this process has a continuous part or not
   * @param aIsContinuous has the process a continuous part?
   */
  void SetIsContinuous( const bool aIsContinuous ) { fIsContinuous = aIsContinuous; }

  /**
   * @brief Method that sets whether this process has an at-rest part or not
   * @param aIsAtRest has the process an at-rest part?
   */
  void SetIsAtRest( const bool aIsAtRest ) { fIsAtRest = aIsAtRest; }

  /**
   * @brief Method that sets the ForcedCondition type of this process
   * @param aForcedCondition is the ForcedCondition type to set
   */
  void SetForcedCondition( const ForcedCondition aForcedCondition ) { 
    fForcedCondition = aForcedCondition;
  }

  /**
   * @brief Method that sets the type of this process
   * @param aProcessType is the process type to set
   */
  void SetType( const ProcessType aType ) { 
    fType = aType;
  }

  /**
   * @brief Method that sets the name of this process
   * @param aName is the name of the process to set
   */
  void SetName( const std::string &aName ) { fName = aName; }

};

}  // end of namespace geant

#endif
