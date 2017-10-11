//===--- FastSimApplication.h - Geant-V ------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file FastSimApplication.h
 * @brief Implementation of FastSim example for Geant-V prototype 
 * @author Mihaly Novak & Alberto Ribon (Apr 2016)
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_FastSimApplication
#define GEANT_FastSimApplication

#ifndef GEANT_VAPPLICATION
#include "GeantVApplication.h"
#endif

#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif

#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif

#include "Geant/Typedefs.h"

#include "GeantFwd.h"

#include "base/Vector.h"

class TH1F;


/** @brief FastSimApplication class */
class FastSimApplication : public Geant::GeantVApplication {

  using GeantRunManager = Geant::GeantRunManager;
  using GeantEvent = Geant::GeantEvent;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;

  private:
    bool fInitialized;  /** Initialized flag */

    std::vector< int > isTrackerVolume;  /** is the volume of type Tracker? */
    std::vector< int > isEcalVolume;     /** is the volume of type ECAL ? */
    std::vector< int > isHcalVolume;     /** is the volume of type HACL ? */

    std::mutex fMHist;  /** Mutex for concurrent histogram filling */

    TH1F* fRatioMomentumInTracker;  /** Histogram of p_smeared/p_true in the Tracker */
    TH1F* fRatioEnergyInEcal;       /** Histogram of Ekin_smeared/Ekin_true in the ECAL */
    TH1F* fRatioEnergyInHcal;       /** Histogram of Ekin_smeared/Ekin_true in the HCAL */

    /** @brief Copy constructor (not allowed) */
    FastSimApplication( const FastSimApplication & );

    /** @brief Operator= (not allowed) */
    FastSimApplication &operator=( const FastSimApplication & );

  public:
    /** @brief Default constructor */
    FastSimApplication(GeantRunManager *runmgr);

    /** @brief Destructor */
    virtual ~FastSimApplication() {}

    /** @brief Initialization */
    virtual bool Initialize();

    /**
     * @brief Provides access to the information of the step
     *        (similar to the Geant4 stepping user action)
     */
    virtual void StepManager( int npart, const GeantTrack_v &tracks, GeantTaskData *td );

    /**
     * @brief Digitization
     * 
     * @param event Event that should be digitized
     */
    virtual void FinishEvent(int /*evt*/, int /*slot*/) {}

    /** @brief User FinishRun function */
    virtual void FinishRun();
};

#endif
