
#ifndef PHYSICS_PROCESS
#define PHYSICS_PROCESS

#include <string>
#include <vector>

#include "GeantTrack.h"
#include "GeantTaskData.h"

namespace geantphysics {

// Forward declarations
class LightTrack;
class Particle;
class LambdaTable;
class PhysicsParameters;

class MaterialCuts;
class Particle;
/**
 * @brief  Base class for all single physics processes.
 * @class  PhysicsProcess
 * @author M Novak, A Ribon
 * @date   november 2015
 *
 * Similarly to the Geant4 G4VProcess concept, this class specifies one single physics process for one particular
 * particle type. This is the base class from which all single physics processes, either in-flight and/or at-rest, must
 * derived from. As in Geant4, it has two main types of methods, "GetLength" (called here differently but with the same
 * property) and "DoIt", and (at most) three main "components": along-the-step (continuous), post-step (discrete) and
 * at-rest.
 */
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
             kMSC,                /** for a multiple scattering process */
             kEnergyLoss,         /** for an energy loss process i.e. ioni and brem */
             kOptical,            /** for a process involving optical photons */
             kHadronic,           /** for a hadronic process */
             kPhotoLeptonHadron,  /** for a gamma- or lepto-nuclear process */
             kFastSim             /** for a fast-sim (e.g. parameterization) process */
           };


/**
 * @brief Class PhysicsProcess
 */
class PhysicsProcess {
public:
  /** @brief PhysicsProcess default constructor */
  PhysicsProcess(const std::string &name);

  /** @brief PhysicsProcess complete constructor */
  PhysicsProcess(const bool aIsDiscrete, const bool aIsContinuous, const bool aIsAtRest,
                 const ForcedCondition aForcedCondition, const ProcessType aType, const std::string &aName);

  /** @brief PhysicsProcess destructor */
  virtual ~PhysicsProcess();

  //
  virtual void Initialize();

// NOTE: set to always applicable now; each process needs to have a list of particles that the process can be assigned
//       to, this list is checked when a process is assigned to a given particle at initialisation time
  /** @brief Method that returns "true" if the process is applicable.
   *  @param A LightTrack containing the current particle and its kinematics.
   *
   */
  virtual bool IsApplicable(const LightTrack &/*track*/) const { return true; }


  /** @brief Method that returns the macroscopic cross section in internal [1/length] unit.
   *
   *  Each process that has discrete part (except some special processes like multiple scattering) needs to implement
   *  this method. The macroscopic cross section for a given material, interaction is defined in terms of the
   *  corresponding elementwise atomic cross sections as their weighted sum. The weights are the number of atoms of the
   *  given element per unit volume of the material.
   *  This method is used:
   *   - either only at initialisation by the PhysicsManagerPerParticle to build both the lambda table per process and
   *     the total lambda table for the partcile. The lambda table per process is used to sample the type of discrete
   *     intercation while the total lambda table is used to sample the discrete physics step limit.
   *   - or only at run time for particles that has only one process for which the macroscopic cross section is easily
   *     computable (e.g. particles that has only decay).
   *
   *  @param[in] matcut     Material production cuts couple object in which the macroscopic cross section is requested.
   *                        Note, that the material of a MaterialCuts object can be obtained from the MaterialCuts
   *                        object.
   *  @param[in] kinenergy  Kinetic energy of the particle in internal [energy] unit.
   *  @param[in] particle   Pointer to the partcile object (NOTE:we pass this for processes that can be used for more than one particle).
   *  @param[in] mass       Dynamic mass of the particle in internal [energy] units.
   *  @return               Computed macroscopic cross section in internal [1/length] unit.
   */
  virtual double ComputeMacroscopicXSection(const MaterialCuts * /*matcut*/, double /*kinenergy*/,
                                            const Particle * /*particle*/, double /*mass*/) const {return 0.;}


  // the minimum kinetic energy for (optional) lambda table (macroscopic cross section) table: will be used only if
  // the lambda table was requested per-material-cuts and even in that case the max of this and the
  // PhysicsParameters->GetMinLambdaTableEnergy() will be taken. Therefore, only production cut dependent processes
  // needs to be reimplement this method.
  virtual double GetMinimumLambdaTableKineticEnergy(const MaterialCuts * /*matcut*/,
                                                    const Particle * /*particle*/) const { return 0.;}


  /** @brief Method that returns the along-the-step limitation length
   *         (unit: length)
   *
   *  This applies only for the continuous part of this process.
   *  If the process does not have a continuous part, or if it has one
   *  that does not limit the step, then the method returns an arbitrary,
   *  very large value.
   */
   virtual double AlongStepLimitationLength(Geant::GeantTrack * /*track*/, Geant::GeantTaskData * /*td*/) const;


   /** @brief Method that returns the post-step limitation length
    *         (unit: length)
    *
    *  This applies only for the cdiscrete part of this process.
    *  If the process does not have a discrete part this method is not called.
    */
   virtual double PostStepLimitationLength(Geant::GeantTrack * /*track*/, Geant::GeantTaskData * /*td*/, bool haseloss=false);


  /** @brief Method that returns the average lifetime of this process
   *         (unit: time)
   *
   *  This applies only for the at-rest part of this process.
   *  If the process does not have an at-rest part, or if it has one
   *  that does not limit the time, then the method returns an arbitrary,
   *  very large value.
   *  The method is meant for the sampling of the type of interactions.
   */
  virtual double AverageLifetime(const LightTrack &track) const;


  // REMINDER-TO-BE-REMOVED: THE NEXT 3 METHODS HAVE LightTrack_v AS OUTPUT.
  // ----------------------  THIS COULD BE CHANGED TO:  LightTrack_v*
  //                         OR TO: std::vector< LightTrack >


  /** @brief Method that does the along-the-step, i.e. the continuous part,
             action of the in-flight process.
   *
   *  The input parameter, a LightTrack object, can be modified: typically
   *  the track status, the local energy deposit and the non-ionizing
   *  local energy deposit can be changed in this method.
   *  New particles, created along-the-step (this is rare, but possible:
   *  e.g. for sub-cutoff production; or from atomic de-excitation) are stored
   *  in the GeantTaskData::PhysicsDada object. Will return
   *  with the number of created secondary tracks that are stored in the sectracks
   *  vector.
   */
  virtual  int AlongStepDoIt(LightTrack & /*track*/, Geant::GeantTaskData * /*td*/) {return 0;}

  /** @brief Method that does the post-step, i.e. the discrete part, action.
   *
   *  The input parameter, a LightTrack object, can be modified: typically
   *  the track status is changed in this method.
   *  The new particles created by the discrete part of the process are stored
   *  in the GeantTaskData::PhysicsDada object.
   */
  virtual int PostStepDoIt(LightTrack & /*track*/, Geant::GeantTaskData * /*td*/)  {return 0;}

  /** @brief Method that does the at-rest action of the process.
   *
   *  The input parameter, a LightTrack object, can be modified: typically
   *  the track status is changed in this method.
   *  The new particles created by the at-rest part of the process are stored
   *  in the GeantTaskData::PhysicsDada object.
   *  Note: this method also includes the sampling of the target atom (Z, N)
   *        where the at-rest process happens.
   */
  virtual void AtRestDoIt(LightTrack & /*track*/, Geant::GeantTaskData * /*td*/) {return;}


  virtual double MacroscopicXSectionMaximumEnergy(const MaterialCuts * /*matcut*/) { return gAVeryLargeValue; }
  // will be called only if GetMacroscopicXSectionMaximumEnergy < gAVeryLargeValue
  virtual double MacroscopicXSectionMaximum(const MaterialCuts * /*matcut*/) { return 0.; }

  double GetMacroscopicXSection(const MaterialCuts *matcut, double ekin, double mass);
  double GetMacroscopicXSectionForStepping(const MaterialCuts *matcut, double ekin, double mass, bool haseloss=false);


  //--- Getters ---
  /** Method to get the index of this process in the per particle process manager */
  size_t  GetIndex() const { return fIndex; }

  /** Method to get the index of this process in the global process table */
  size_t  GetGlobalIndex() const { return fGlobalIndex; }

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

  bool  IsLambdaTablerequested() const { return fIsLambdaTableRequested; }

  const PhysicsParameters* GetPhysicsParameters() const { return fPhysicsParameters; }

  const Particle* GetParticle() const { return fParticle; }

  /** Methods that are needed for the new concept of physics-per-region:
      is this physics process active in the i-th region? */
  std::vector< bool >& GetListActiveRegions() { return fListActiveRegions; }
  bool IsActiveRegion(const int regionindx) const { return fListActiveRegions[regionindx]; }

  //--- Setters ---

  void SetIndex(size_t val) { fIndex = val; }

  /**
   * @brief Method that sets whether this process has a discrete part or not
   * @param aIsDiscrete has the process a discrete part?
   */
  void SetIsDiscrete(const bool aIsDiscrete) { fIsDiscrete = aIsDiscrete; }

  /**
   * @brief Method that sets whether this process has a continuous part or not
   * @param aIsContinuous has the process a continuous part?
   */
  void SetIsContinuous(const bool aIsContinuous) { fIsContinuous = aIsContinuous; }

  /**
   * @brief Method that sets whether this process has an at-rest part or not
   * @param aIsAtRest has the process an at-rest part?
   */
  void SetIsAtRest(const bool aIsAtRest) { fIsAtRest = aIsAtRest; }

  /**
   * @brief Method that sets the ForcedCondition type of this process
   * @param aForcedCondition is the ForcedCondition type to set
   */
  void SetForcedCondition(const ForcedCondition aForcedCondition) {
    fForcedCondition = aForcedCondition;
  }

  /**
   * @brief Method that sets the type of this process
   * @param aProcessType is the process type to set
   */
  void SetType(const ProcessType aType) {
    fType = aType;
  }

  void SetPhysicsParameters(const PhysicsParameters *physpars) { fPhysicsParameters = physpars; }

  void SetParticle(const Particle *particle) { fParticle = particle; }

  void RequestLambdaTables(bool ispermaterial=true);

  void SetSpecialLambdaTableBinNum(int val);

  void BuildLambdaTables();

  // kinetic energy of the maximum of macroscopic cross section: might be used in case of particles that has energy-loss p
  // very large value indicates that actualy it won't play any role (and the GetMacroscopicXSectionMaximum won't be used)
  double GetMacroscopicXSectionMaximumEnergy(const MaterialCuts * /*matcut*/);
  double GetMacroscopicXSectionMaximum(const MaterialCuts * /*matcut*/);

  void AddToListParticlesAssignedTo(Particle* part) { fListParticlesAssignedTo.push_back(part); }
  const std::vector<Particle*>& GetListParticlesAssignedTo() const { return fListParticlesAssignedTo; }

  void AddToListParticlesAlloedToAssigned(Particle* part) { fListParticlesAlloedToAssigned.push_back(part); }
  const std::vector<Particle*>& GetListParticlesAlloedToAssigned () const { return fListParticlesAlloedToAssigned; }

  // to delete all created process object; must be called only by the main manager
  static void ClearAllProcess();


  static double GetAVeryLargeValue() {return gAVeryLargeValue;}

  static const PhysicsProcess*  GetProcessByGlobalIndex(size_t gprocindx) {
    PhysicsProcess *proc = nullptr;
    if (gprocindx<gThePhysicsProcessTable.size())
      proc = gThePhysicsProcessTable[gprocindx];
    return proc;
  }

private:
  /** @brief PhysicsProcess copy constructor is not defined */
  PhysicsProcess(const PhysicsProcess &other);
  /** @brief Operator = is not defined*/
  PhysicsProcess& operator=(const PhysicsProcess &other);



private:
  size_t  fIndex;         /** Index of this process in the per particle process manager */
  size_t  fGlobalIndex;   /** Index of this process in the global process table */
  bool    fIsDiscrete;    /** "true" if the process has a discrete part; "false" otherwise */
  bool    fIsContinuous;  /** "true" if the process has a continuous part; "false" otherwise */
  bool    fIsAtRest;      /** "true" if the process has an at-rest part; "false" otherwise */
  bool    fIsLambdaTableRequested; /** false to store if lambda table was requested by the process (flase by default) */

  ForcedCondition fForcedCondition;  /** type of ForcedCondition for this process */
  ProcessType fType;                 /** type of this process (MAYBE USEFUL IN THE FUTURE) */
  std::string fName;                 /** name of the process (useful for debugging) */
  std::vector< bool > fListActiveRegions;  /** is this process active in the i-th region?;
                                               will be set by the PhysicsListManager */
  std::vector<Particle*> fListParticlesAssignedTo; /** list of particles this process is assigned to;
                                                       this list is determined by the physics list; do NOT own the
                                                       Particle objects
                                                   */
  std::vector<Particle*> fListParticlesAlloedToAssigned; /** list of particles that this process can be used;
                                                             this list is determined by the developer, do NOT own the
                                                             Particle objects */

  static const double    gAVeryLargeValue; // for the along step limitation of those that do not limit the step
  // unique collection of process object pointers that has been created so far; will be used to delete all
  // processes; each process stores its index in this table in its fGlobalIndex member and processes can be obtained
  // by this index using the GetProcessByGlobalIndex(global_index) static metho
  static std::vector<PhysicsProcess*> gThePhysicsProcessTable;

  // pointer to the PhysicsParameters object active in region(s) where this process is active (will be set by the
  // PhysicsManagerPerParticle )
  const PhysicsParameters  *fPhysicsParameters;  // not owned by the process

  // particle to which the process is assinged (will be set by the PhysicsManagerPerParticle)
  const Particle           *fParticle;

  // optional lambda table: table of macroscopic cross sections per materia/material-cuts
  LambdaTable              *fLambdaTable;   // owned by the process”
};

}  // end of namespace geantphysics

#endif
