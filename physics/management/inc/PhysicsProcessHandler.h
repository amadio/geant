//==========================================================================
// PhysicsList.h - Geant-V Prototype
/**
 * @file PhysicsProcessHandler.h
 * @brief Class that handles the physics simulation in one step
 *
 * This class decides and calls the appropriate methods to simulate the
 * physics of any given particle happening in one step.
 * This final class derives from the  GeantPhysicsInterface class
 *
 * @author Mihaly Novak & Alberto Ribon (Dec 2015)
 */
//==========================================================================

#ifndef PHYSICS_PROCESS_HANDLER
#define PHYSICS_PROCESS_HANDLER


namespace geant {

class GeantPhysicsInterface {
  // THIS IS A FAKE CLASS, JUST TO COMPILE, TO TAKE THE PLACE OF THE CURRENT
  // GEANTV CLASS  PhysicsProcess  WHICH SHOULD BE RENAMED  GeantPhysicsInterface
};


/**
 * @brief Class PhysicsProcessHandler
 */
class PhysicsProcessHandler : public GeantPhysicsInterface {
private:
  // EMPTY FOR THE TIME BEING!


  /** @brief PhysicsProcessHandler copy constructor is not defined */
  PhysicsProcessHandler( const PhysicsProcessHandler &other );

  /** @brief Operator = is not defined */
  PhysicsProcessHandler& operator=( const PhysicsProcessHandler &other );
 
public:
  /** @brief PhysicsList default constructor */
  PhysicsProcessHandler();

  /** @brief PhysicsProcessHandler destructor */
  virtual ~PhysicsProcessHandler();

  /** @brief Initialize the physics
   *
   *  This method, which of course is called at initialization, fetches the
   *  production and tracking cuts, the material-cut couple table, and the
   *  number of working threads; then calls  PhysicsList::Initialize .
   */
  virtual void Initialize();

  /** @brief The method proposes the step-length from the physics 
   *  
   *  @param mat Material_t material
   *  @param ntracks Number of tracks
   *  @param tracks Vector of tracks_v
   *  @param lengths Partial process lengths
   *  @param td Thread data
   *
   *  This length is the minimum of the following lengths:
   *  - the sampled value from an exponential distribution whose parameter
   *    is given by the total lambda table obtained from PhysicsManagerPerParticle
   *  - the returned value of the method PhysicsProcess::AlongStepLimitationLength
   *    for each of the continuous processes associated to the particle.
   *  Note: if the winner is a continuous process, then the process index
   *        should be set negative.
   *  Note: this method should never be called for a particle at rest.
   */
  virtual void ComputeIntLen( /* Material_t *mat, int ntracks, GeantTrack_v &tracks, 
                                 double *lengths, GeantTaskData *td */ );

  /** @brief The method invokes the "AlongStepDoIt" method of each continuous process
   *         associated to the particle
   *
   *  @param mat Material_t material
   *  @param ntracks Number of tracks
   *  @param tracks Vector of tracks_v
   *  @param nout Number of tracks in the output
   *  @param td Thread data
   *
   *  Note that continuous processes do not compete but collaborate with each other.
   */
  virtual void AlongStepAction( /* Material_t *mat, int ntracks, GeantTrack_v &tracks, 
                                   int &nout, GeantTaskData *td */ );
  /** @brief The method selects the winner discrete process, then the target,
   *         and finally produces a final-state.
   *
   *  @param mat Material_t material
   *  @param ntracks Number of tracks
   *  @param tracks Vector of tracks_v
   *  @param nout Number of tracks in the output
   *  @param td Thread data
   *
   *  This method does the following:
   *  1. Selects the winner process between all the discrete processes
   *     associated to the particle by drawing a random number and using
   *     their  PhysicsProcess::InverseLambda  values.
   *  2. Selects the target (Z, N) by calling PhysicsProcess::SampleTarget .
   *  3. Calls, for the winner discrete process (and eventually the "Forced"
   *     discrete processes), the method PhysicsProcess::PostStepDoIt .
   *  Note that discrete processes compete against each other. 
   *  
   *  NOTE-TO-BE-DELETED: NOT YET CLEAR AT THE MOMENT WHETHER WE NEED
   *                      KEEP THE TWO METHODS:
   *                      PostStepTypeOfIntrActSampling AND PostStepFinalStateSampling
   *                      WHICH ARE USED BY THE TABULATED PHYSICS.
   */
  virtual void PostStepAction( /* Material_t *mat, int ntracks, GeantTrack_v &tracks, 
                                  int &nout, GeantTaskData *td */ );

  /** @brief The method selects the winner at-rest process, then the target
   *         if needed (not for decays) and finally produces a final-state.
   *
   *  @param mat Material_t material
   *  @param ntracks Number of tracks
   *  @param tracks Vector of tracks_v
   *  @param nout Number of tracks in the output
   *  @param td Thread data
   *
   *  This method does the following:
   *  1. Selects the winner process between all the at-rest processes
   *     associated to the particle by drawing a random number and using
   *     their  PhysicsProcess::AverageLifetime  values.
   *  2. Calls, for the winner process (and eventually the "Forced" at-rest
   *     processes) the method  PhysicsProcess::AtRestDoIt , which, if needed
   *     (e.g. nuclear capture, but not for decays), selects also the target (Z, N).
   *  Note that at-rest processes compete against each other. 
   *
   *  NOTE-TO-BE-DELETED: THE SIGNATURE IS SLIGHTLY DIFFERENT FROM THE
   *                      ORIGINAL: NEED TO ADD THE MATERIAL POINTER, AS
   *                      FOR AlongStepAction AND AtRestAction.
   */
  virtual void AtRestAction( /* Material_t *mat, int ntracks, GeantTrack_v &tracks, 
                                int &nout, GeantTaskData *td */ );

};

}  // end of namespace geant

#endif
