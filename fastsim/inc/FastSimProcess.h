//==========================================================================
// FastSimProcess.h - Geant-V Prototype
/**
 * @file FastSimProcess.h
 * @brief Fast simulation example for the Geant-V prototype
 *
 * This example is based on the Geant4 fast simulation extended example:
 *   examples/extended/parameterisations/Par02/
 * which, in turn, is based on a simplified version of a Geant4-based
 * fast simulation application written by Anna Zaborowska for 
 * Future Circular Collider (FCC) studies.
 *
 * This example shows how to do "track and energy smearing", in order
 * to have a very fast simulation based on assumed detector resolutions.
 * Any charged particle is smeared in the tracker: its momentum (the module
 * only, not the direction) is smeared according to a gaussian, at the
 * exit of the tracker.
 * Any electron, or positron, or gamma is smeared in the electromagnetic
 * calorimeter: the particle is killed at the entrance of the ECAL, with
 * a deposited energy obtained by gaussian smearing of its kinetic energy.
 * Any hadron is smeared in the hadronic calorimeter: the particle is killed
 * at the entrance of the HCAL, with a deposited energy obtained by gaussian
 * smearing of its kinetic energy.
 *
 * This class interfaces with the GeantV kernel, while the actual 
 * fast simulation work is delegated to the class Smearer. 
 *
 * @author Mihaly Novak & Alberto Ribon (Apr 2016)
 */
//==========================================================================

#ifndef GEANT_FASTSIMPROCESS
#define GEANT_FASTSIMPROCESS

#include "Geant/Config.h"

#ifndef GEANT_PHYSICSPROCESS
#include "PhysicsProcess.h"
#endif

#include "base/Global.h"
#include "Geant/Typedefs.h"

#include "GeantFwd.h"


class Smearer;


/**
 * @brief Class FastSimProcess
 */
class FastSimProcess : public PhysicsProcess {

  public:

    using GeantTrack = Geant::GeantTrack;
    using GeantTrack_v = Geant::GeantTrack_v;

    /** @brief Default constructor */
    FastSimProcess();

    /** @brief Destructor */
    virtual ~FastSimProcess();

    /** @brief Initialization: an object of the Smearer class is created */
    virtual void Initialize();

    /** @brief Compute the proposed step-lengths for the vector of tracks.
     *
     *  The proposed step-lengths are assigned to  tracks.fPstepV[i] .
     *  We also set  tracks.fEindexV[i]  to the value "1000" in order to treat
     *  the fast simulation as a continuous & discrete process.
     */
    virtual void ComputeIntLen( Material_t *mat, int ntracks, GeantTrack_v &tracks, 
                                double *lengths, GeantTaskData *td );

    /** @brief This method does the "continuous" part of the fast simulation
     *
     *  This method is used only for the tracker parameterisation, in order to apply
     *  the smearing of the module of the momentum of each charged particle track
     *  at the exit of the tracker (while the transportation of the track inside
     *  the tracker is done normally).
     */
    virtual void Eloss( Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout,
                        GeantTaskData *td );

    /** @brief Dummy method */
    virtual void PostStep( Material_t * /* mat */, int /* ntracks */, 
                           GeantTrack_v & /* tracks */, int & /* nout */, 
                           GeantTaskData * /* td */ ) {}

    /** @brief First part of PostStep: not used for fast simulation. */
    virtual void PostStepTypeOfIntrActSampling( Material_t *mat, int ntracks, 
                                                GeantTrack_v &tracks, GeantTaskData *td );

    /** @brief Second part of PostStep: discrete part of the fast simulation
     *
     *  Currently this method is used for both electromagnetic and hadronic
     *  calorimeter parameterisations.
     *  When the fast simulation is applied (e.g. for electrons, positrons and
     *  gammas for the ECAL, and for hadrons in the HCAL), the particle is killed
     *  at the entrance of the calorimeter, depositing an energy obtained by
     *  gaussian smearing of the kinetic energy of the particle.
     */
    virtual void PostStepFinalStateSampling( Material_t *mat, int ntracks, GeantTrack_v &tracks,
                                             int &nout, GeantTaskData *td );

    /** @brief FastSim at rest: not used, assumed only discrete & continuous. */ 
    virtual void AtRest( int /* ntracks */, GeantTrack_v & /*tracks*/ , int & /* nout */, 
                         GeantTaskData * /* td */ ) {};

    /** @brief FastSim does nothing for multiple scattering. */
    virtual void ApplyMsc( Material_t * /* mat */, int /* ntracks */, 
                           GeantTrack_v & /* tracks */, GeantTaskData * /* td */ ) {};

  private:

    /** @brief Copy constructor (not allowed) */
    FastSimProcess( const FastSimProcess& );

    /** @brief Operator= (not allowed) */
    FastSimProcess& operator=( const FastSimProcess& );

    Smearer* fSmearer;  /** Pointer to the object that does the fast simulation */
#ifdef USE_ROOT
    ClassDef( FastSimProcess, 1 )
#endif
};

#endif
