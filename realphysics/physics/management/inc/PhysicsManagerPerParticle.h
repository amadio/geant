
#ifndef PHYSICS_MANAGER_PER_PARTICLE
#define PHYSICS_MANAGER_PER_PARTICLE

#include <string>
#include <vector>

#include "GeantTaskData.h"

namespace geantphysics {

 // forward declarations
class PhysicsParameters;
class Material;
class MaterialCuts;
class PhysicsProcess;
class Particle;
class Spline;
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
   *   Initialize all processes assigned to the particle that this manager belongs to and build all lambda tables.
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

  /**
   * @brief Run time method to obtain total macroscopic cross section for a given MaterialCuts/Material and kinetic
   *         energy.
   *
   *  This method is used at run time to obtain the restricted/full macroscopic cross section that is the inverse of the
   *  mean path length (mean free path) that the particle travels in the given material with the given kinetic energy
   *  between two discrete interactions. It is used to sample the path, that the particle travels till the next discrete
   *  interaction at run time. The total macroscopic scross section is obtained by spline interpolation based on the
   *  total lambda table built at initialisation by computation. The MaterialCuts object provided as input parameter
   *  contains pointer to the corresponding Material object. The PhysicsManagerPerParticle object knows if the lambda
   *  tables are per MaterialCuts(if the particle has any EnergyLoss process i.e. ionization or bremsstrahlung) or per
   *  Material and the appropriate index is used to obtain the lambda table.
   *
   *  @param[in]  matcut    Pointer to the MaterialCuts object where the total macroscopic cross section is requested.
   *  @param[in]  kinenergy Kinetic energy of the particle at which the total macroscopic cross section is requested in
   *                        in internal [energy] units.
   *  @return     Total macroscopic cross section in the given MaterialCuts/Material, at the given kinetic energy in
   *              internal [1/length] units.
   */
  double GetInvTotalLambda(const MaterialCuts *matcut, double kinenergy);
  // when possible energy loss along the step is accounted during the discrete step limit
  double ComputeInvTotalLambdaForStepLimit(const MaterialCuts *matcut, double ekin);


  /**
   * @brief Run time method to select one of the discrete interactions assigned to the particle for post-step interaction.
   *
   * The discrete interaction, for the provided MaterialCuts/Material and particle kinetic energy, is selected based on
   * the per-process partial macroscopic cross section table (if there are more than one discrete processes assigned to
   * the particle). The partial macroscopic cross sections are obtained by interpolation based on the approriate lambda
   * table built at initialisation by computation.
   *
   *  @param[in]  matcut    Pointer to the MaterialCuts object where the total macroscopic cross section is requested.
   *  @param[in]  kinenergy Kinetic energy of the particle at which the total macroscopic cross section is requested in
   *                        in internal [energy] units.
   *  @param[in]  presteplambda  Total mean free path for discrete interaction that was set at the step limit.
   *  @param[in]  td         Pointer to the GeantTaskData.
   *  @return     Pointer to the selected discrete physics process. nullptr if delta interaction.
   */
  PhysicsProcess* SelectDiscreteInteraction(const MaterialCuts *matcut, double kinenergy, double presteplambda,
                                            Geant::GeantTaskData *td);


  /**
   * @brief Method to get macroscopic cross section for a given discrete process (used only for testing!)
   *
   * This method will use the lambda table(s), built at initialisation to interpolate the per-process macroscopic
   * corss section. This method is used only for testing and not used during a normal simulation.
   *
   *  @param[in]  matcut    Pointer to the MaterialCuts object where the total macroscopic cross section is requested.
   *  @param[in]  kinenergy Kinetic energy of the particle at which the total macroscopic cross section is requested in
   *                        in internal [energy] units.
   *  @param[in]  procindex Index of the (discrete) process in the fPostStepCandidateProcessVec.
   *  @return     (Interpolated) macroscopic cross section for the given process, MaterialCuts(Material) at the given
   *              kinetic energy in internal [1/length] units.
   */
  double GetMacroscopicXSectionForProcess(const MaterialCuts *matcut, double kinenergy, size_t procindex);


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

  void ComputeIntLen(LightTrack &track, Geant::GeantTaskData *td);
  int  AlongStepAction(LightTrack &track, std::vector<LightTrack> &sectracks);
  int  PostStepAction(LightTrack &track, std::vector<LightTrack> &sectracks, Geant::GeantTaskData *td);

/*
// only for testing; should be removed later
  void PrintLambda(const MaterialCuts *matcut){
    LambdaTable *lambTab = nullptr;
    if (!fIsHasTotalLambdaTable) {
      std::cerr<<"   ****** No total LambdaTable for Particle = "<<fParticle->GetName()<<std::endl;
      return;
    }
    if (fIsLambdaTablesPerMaterial) {
      lambTab = fLambdaTables[matcut->GetMaterial()->GetIndex()];
    } else {
      lambTab = fLambdaTables[matcut->GetIndex()];
    }
    if (lambTab && lambTab->fInvTotalLambda) {
      for (int i=0; i<lambTab->fNumData; ++i) {
        std::cout<<fEnergyGrid[i]/geant::MeV<< "  " <<lambTab->fInvTotalLambda[i]/(1.0/geant::cm);
        if (fIsHasPerProcessLambdaTable) {
          for (int j=0; j<lambTab->fInvLambdaPerProcess.size(); ++j)
            std::cout<<"  "<<lambTab->fInvLambdaPerProcess[j][i]/(1.0/geant::cm);
        }
        std::cout<<std::endl;
      }
    }
  }
*/

// Inverse lambda tables i.e. macroscopic cross sections per material/material-cut: total and per (discrete)processes:
struct LambdaTable {
  int                    fNumData;              // number of data points i.e. size of the arrays
  const  MaterialCuts   *fMaterialCut;          // not owned
  const  Material       *fMaterial;             // not owned
  double                *fInvTotalLambda;       // owned
  double                 fLambdaMax;            // maximum of the total macroscopic cross section
  double                 fLambdaMaxEnergy;      // energy where the total macroscopic cross section has its maximum
  Spline                *fSplineInvTotalLambda; // owned
  std::vector<double*>   fInvLambdaPerProcess;  // owned
  std::vector<Spline*>   fSplinesInvLambdaPerProcess; // owned
};


//
// private methods
//
private:
  /** @brief PhysicsManagerPerParticle copy constructor is not defined */
  PhysicsManagerPerParticle(const PhysicsManagerPerParticle &other);
  /** @brief Operator = is not defined */
  PhysicsManagerPerParticle& operator=(const PhysicsManagerPerParticle &other);


  /** @brief Method to build lambda tables at initialization.
   *
   *  Lambda tables (if any) i.e. total(a) and per-discrete process(b) macroscopic cross section tables are built over a
   *  kinetic energy grid defined in the corresponding PhysicsParametes object and set up for
   *  - run time total discrete step limit(a) (i.e. total inverse lambda)
   *    double PhysicsManagerPerParticle::GetInvTotalLambda(const MaterialCuts *matcut, double kinenergy)
   *  - discrete process sampling(b)
   *    const PhysicsProcess* PhysicsManagerPerParticle:::SelectDiscreteInteraction(const MaterialCuts *matcut,
   *                                                                                double kinenergy, double rndm)
   *
   *  This method is expected to be called at initialization.
   */
  void BuildLambdaTables();

  // check the dicrete processes assigned to the particle
  // 1. check if it has any process that has discrete part that are neither MSC nor decay:
  //            - MSC do not have contribution to the macroscopic cross section (i.e. inv-lambda) table
  //            - if the only discrete process (beyond possible MSC) is deacy, we do not build any (neither total nor
  //              per process) lambda tables since lambda from decay can be computed on the fly easily
  // 2. if it has at least one discrete process that is not MSC or decay:
  //            -  we will  build the total lambda table:
  //                  a. if it has decay as well we handle decay as normal discrete process i.e. the particle has 2
  //                     discrete processes so both total lambda and per process lambda tables will be built with the
  //                     second discrete process being the decay
  //                  b. if it doesn't have decay, then the partcile has only one discrete process so we build only
  //                     the total lambda table, but not the per process table -> if discrete interaction happens it's
  //                     always the only one that the particle has
  // 3. if it has more than one discrete processes that is not MSC or decay:
  //            -  we will build both the total lambda table and the lambda tables per processes:
  //                   # if it has decay as well we handle decay as normal discrete process i.e. decay will be one of
  //                     the discrete processes so its lambda contribution will go both to the total lambda and per
  //                     process lambda tables as well
  // 4. discrete processes will also be reordered:
  //     - all the discrete processes (including decay it the particle has, and starting with that) will be at the
  //       beginning of the discrete process vector and MSC (if the particle has) will be the last
  //
  void CheckForLambdaTable();

  // utility method to build lambda table for one MaterialCuts/Material
  void BuildOneLambdaTable(const MaterialCuts *matcut, int indx);

  // set up the common energy grid that is determined by some members of the fPhysicsParameters
  void InitializeEnergyGrid();

  // delete all lambda tables
  void ClearLambdaTables();



private:
   /** Pointer to the partcile object to which this physics manager belongs to.
    *  The class doesn't own the partcile object.
    */
  const Particle            *fParticle;           // not owned;
  const PhysicsParameters   *fPhysicsParameters;  // not owned;

  // some flags to indicate ...
  bool        fIsLambdaTablesPerMaterial;  // by def false
  bool        fIsHasTotalLambdaTable;
  bool        fIsHasPerProcessLambdaTable;
  bool        fIsHasMSCProcess;
  bool        fIsHasDecayProcess;
  bool        fIsHasElossProcess;
  // the energy grid
  int         fNumLambdaTableBins;
  double      fMinLambdaTableEnergy;
  double      fMaxLambdaTableEnergy;
  double      fLogMinLambdaTableEnergy;
  double      fEnergyILDelta;
  double     *fEnergyGrid;     // common energy grid determined by some of the members of the PhysicsParameters; owned


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


  std::vector<LambdaTable*>  fLambdaTables; // lambda tables per material/material-cut; size will #of material/materilaCuts
                                            // all lambda tables are owned by the class

  std::vector<bool> fListActiveRegions;  /** is this physics manager per particle active in the i-th region? */

};

}  // end of namespace geantphysics

#endif
