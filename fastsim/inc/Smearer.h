//==========================================================================
// Smearer.h - Geant-V Prototype
/**
 * @file Smearer.h
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
 * The geometry which is considered is a simplified collider detector set-up,
 * inspired by ALEPH/ATLAS/CMS detectors. Although it is much simpler than
 * a realistic detector, it is anyhow fairly complex and therefore build up
 * from a GDML file, Par02FullDetector.root .
 * In this example:
 * - Particles with transverse momentum less than 1 MeV or pseudorapidity
 *   larger (in module) than 5.5 are neglected (i.e. the corresponding
 *   track is killed).
 * - Any charged particle is smeared in the tracker as follows:
 *   its momentum (the module only, not the direction) is smeared
 *   according to a gaussian, with mean equal to 1.0 and sigma taken from
 *   the momentum resolution of the CMS tracker (with ALEPH or ATLAS tracker
 *   as a possible alternative), and then placed at the end of the tracker,
 *   at the position that it would reach if normally transported (i.e.
 *   without smearing).
 * - Any electron, or positron, or gamma is smeared in the electromagnetic
 *   calorimeter as follows: 
 *   it is killed at the entrance of the electromagnetic calorimeter, with
 *   a deposited energy equal to the gaussian smearing (with mean equal to
 *   1.0 and sigma taken from the energy resolution of the CMS electromagnetic
 *   calorimeter - with ALEPH or ATLAS electromagnetic calorimeter as a
 *   possible alternative) of its kinetic energy (at the entrance of the
 *   electromagnetic calorimeter).
 * - Any hadron is smeared in the hadronic calorimeter as follows:
 *   it is killed at the entrance of the hadronic calorimeter, with a
 *   deposited energy equal to the gaussian smearing (with mean equal to 1.0
 *   and sigma taken from the energy resolution of the CMS hadronic 
 *   calorimeter - with ALEPH or ATLAS hadronic calorimeter as a possible
 *   alternative) of its kinetic energy (at the entrance of the hadronic
 *   calorimeter).
 *
 * @author Mihaly Novak & Alberto Ribon (Apr 2016)
 */
//==========================================================================

#ifndef GEANT_Smearer
#define GEANT_Smearer

#include "Geant/Config.h"
#include "base/Global.h"
#include "Geant/Typedefs.h"
#include "Geant/Typedefs.h"
#include "GeantFwd.h"


class TGeoManager;
class GeantTrack;
class GeantTrack_v;
class GeantTaskData;


/** @brief Class Smearer */
class Smearer {

  public:

    using GeantTrack = Geant::GeantTrack;
    using GeantTrack_v = Geant::GeantTrack_v;
    using GeantTaskData = Geant::GeantTaskData;

    /** @brief Smearer default constructor
     *  
     *  Currently in GeantV there is no association between volumes and eventual
     *  fast simulation parameterisations.
     *  Therefore, for each type of fast simulation parameterisation, we use a
     *  vector of integers for bookkeeping whether a volume has such fast-sim
     *  model defined (1) or not (0); we use the unique number identifier of
     *  each volume as the index of the vector.
     */
    Smearer();

    /** @brief Smearer destructor */
    ~Smearer() {}

    enum ActivatedParamType { 
      NONE,           /** No fast-sim parameterization */
      TRACKER_PARAM,  /** Activated tracker type of fast-sim parameterization */
      ECAL_PARAM,     /** Activated EM Calorimeter type of fast-sim parameterization */
      HCAL_PARAM,     /** Activated HAD Calorimeter type of fast-sim parameterization */
      MUON_PARAM      /** Activated Muon detector type of fast-sim parameterization */
    };  

    enum ParameterisationType { 
      eCMS,           /** CMS-like detector resolutions */
      eATLAS,         /** ATLAS-like detector resolutions */
      eALEPH          /** ALEPH-like detector resolutions */
    };    

    enum DetectorType { 
      eTRACKER,       /** Tracker type of detector */
      eEMCAL,         /** EM Calorimeter type of detector */
      eHCAL           /** HAD Calorimeter type of detector */
    };


    /** @brief Method that returns the vector of proposed step-lengths for the
     *         vector of tracks by the fast simulation.
     *
     *  For a given track, the returned proposed step-length is a very large value
     *  in the case that there is no fast simulation applicable for that track
     *  in its current position.
     *  If there is a fast simulation parameterisation applicable to a given track,
     *  then the returned value is either 0.0 - if the parameterisation does 
     *  something as soon as the track enters a volume - or a very large value -
     *  if the parameterisation does something only at the end of a volume, i.e.
     *  at the level of "Eloss".
     *  In our case, we use the first approach (i.e. returning 0.0) for 
     *  calorimeter (ECAL and HCAL) parameterisations, whereas the second approach
     *  (i.e. returning a very large value) is used for the tracker parameterisation.
     */
    std::vector< double > StepLengthProposedByParameterisation( int ntracks, 
                                                                GeantTrack_v& tracks,
                                                                GeantTaskData& aTaskData );

    /** @brief Method that applies the fast simulation whenever it is appropriate
     *         for the vector of tracks
     *
     *  The tracker parameterisation is applied to all charged particles at the
     *  exit of the tracker.
     *  The ECAL (electromagnetic calorimeter) parameterisation is applied to
     *  all electrons, positrons, and gammas as soon as they enter the ECAL.
     *  The HCAL (hadronic calorimeter) parameterisation is applied to all
     *  types of hadrons as soon as they enter the HCAL.
     *  Note: We use the arbitrarily chosen number 999 as process number identifier
     *        for fast simulation. This is used only in the method: 
     *          FastSimApplication::StepManagerFast
     *        to decide when to fill the histograms.
     * 
     *  @param isElossCalling is "true" if the method is called from Eloss, 
     *         "false" otherwise (i.e. called from PostStepFinalStateSampling)
     */
    void ApplyParameterisation( int ntracks, GeantTrack_v& tracks, GeantTaskData& aTaskData,
                                bool isElossCalling = false );

  private:

    /** @brief Copy constructor (not allowed) */
    Smearer( const Smearer& );

    /** @brief Operator= (not allowed) */
    Smearer& operator=( const Smearer& );

    /** @brief Methods that return the step-length proposed by the different types
     *         of fast simulation parameterisations */
    double StepLengthProposedByTrackerParameterisation( GeantTrack& aTrack );
    double StepLengthProposedByEcalParameterisation( GeantTrack& aTrack );
    double StepLengthProposedByHcalParameterisation( GeantTrack& aTrack );
    double StepLengthProposedByMuonParameterisation( GeantTrack& aTrack );

    /** @brief Methods that tell whether a track has a given type of parameterisation */
    bool IsTrackerParameterisationApplicable( GeantTrack& aTrack );
    bool IsEcalParameterisationApplicable( GeantTrack& aTrack );
    bool IsHcalParameterisationApplicable( GeantTrack& aTrack );
    bool IsMuonParameterisationApplicable( GeantTrack& aTrack );

    /** @brief Methods that apply the different types of fast simulation
     *         parameterisations whenever it is appropriate for the vector of tracks */
    void ApplyTrackerParameterisation( GeantTrack_v& tracks, GeantTaskData& aTaskData, int index );
    void ApplyEcalParameterisation( GeantTrack_v& tracks, GeantTaskData& aTaskData, int index );
    void ApplyHcalParameterisation( GeantTrack_v& tracks, GeantTaskData& aTaskData, int index );
    void ApplyMuonParameterisation( GeantTrack_v& tracks, GeantTaskData& aTaskData, int index );

    /** @brief Methods that perform the smearing
     *
     *  The smearing of the energy depends on the kinetic energy of the
     *  particle and the energy resolution of the calorimeter (EM or HAD).
     *  The smearing of the momentum - only the module, while the direction
     *  is kept unchanged - depends on the momentum of the particle and
     *  the momentum resolution of the tracker.
     */
    double EnergySmearing( const double ekin, const double resolution );
    double MomentumSmearing( const double p, const double resolution );

    /** @brief Method that returns the resolution
     *
     *  The resolution depends on the detector type, type of parameterisation,
     *  and momentum of the particle.
     */
    double GetResolution( DetectorType aDetector, ParameterisationType aParameterisation,
                          double aMomentum );

    /** @brief Method that returns the efficiency (between 0 and 1)
     *
     *  The efficiency depends on the detector type, type of parameterisation,
     *  and momentum of the particle.
     *  Currently this method is useless, i.e. it returns always 1.
     */    
    double GetEfficiency( DetectorType aDetector, ParameterisationType aParameterisation,
                          double aMomentum );

    /** The following vectors of integer tell whether a volume, whose unique identifier
     *  number is the index of the vector, has a given type of fast-sim parameterization
     *  or not, e.g. hasVolumeHcalParameterisation[ 7 ]  is either "1" - if the volume
     *  with number 7 has a HCAL type of fast-sim parameterisation, or "0" otherwise.
     *  These vectors are set at initialization and then remain unchanged.
     */
    std::vector< int > hasVolumeTrackerParameterisation;
    std::vector< int > hasVolumeEcalParameterisation;
    std::vector< int > hasVolumeHcalParameterisation;
    std::vector< int > hasVolumeMuonParameterisation;
    std::vector< int > hasVolumeParameterisation;
};


#endif
