
#ifndef LIGHT_TRACK
#define LIGHT_TRACK

#include <atomic>

namespace geantphysics {

/**
 * @brief  Light track for the physics of Geant-V prototype
 * @class  LightTrack
 * @author M Novak, A Ribon
 * @date   november 2015
 *
 * As the name suggests, this class is supposed to be a smaller, lighter version of GeantTrack, to be used by the
 * physics processes. It is foreseen that different processes, in particular electromagnetic and hadronic ones, need
 * different information: this should be placed in a new object, whose class should be inherited from the (now empty)
 * ExtraInfo class, and the pointer to this object is kept in the LightTrack.
 */

/**
 * @brief Class ExtraInfo : an empty class to inherit from in order
 *                          to put extra information.
 */
class ExtraInfo {};  // Base class for any "ExtraInfo" class

enum class LTrackStatus {
  kAlive,         /** track is alive i.e. has enough kinetic energy to track further */
  kStopButAlive,  /** stopped but need to check if it has any at-rest interaction */
  kKill,          /** track reached the end of its life */
  kNew            /** new track */
};

/**
 * @brief Class LightTrack
 */
class LightTrack {
private:
  LTrackStatus fTrackStatus;    /** Status of the track */
  int fGVcode;                  /** GV particle code */
  int fGTrackIndex;             /** Index of the track in the GeantTrack_v vector */
  int fMaterialCutCoupleIndex;  /** Index of the material-cut couple where the track is */
  int fProcessIndex;            /** Index of the selected process (in the list of
                                    active discrete or at-rest processes kept by the
                                    PhysicsManagerPerParticle; negative if the winner
                                    is a continuous process) */
  int fTargetZ;                 /** Atomic number (Z) of the target atom where the
                                    particle has a physical interaction */
  int fTargetN;                 /** Number of nucleons of the target atom where the
                                    particle has a physical interaction */

  double fXdir;        /** X direction (normalized, adimensional) of the particle */
  double fYdir;        /** Y direction (normalized, adimensional) of the particle */
  double fZdir;        /** Z direction (normalized, adimensional) of the particle */
  double fKinE;        /** Kinetic energy (in GeV) of the particle */
  double fMass;        /** Dynamic mass (in GeV) of the particle */
  double fTime;        /** Time (global, in sec) of the particle */
  double fWeight;      /** Weight (adimensional) of the particle */
  double fStepLength;  /** True, physical length (in cm) of the last step of the particle */
  double fEdep;        /** Energy deposit (in GeV) in the last step of the particle */
                       /** Note: in the future we could add the non-ionizing energy deposit
                                 in the last step of the particle */
  double fNintLen;     /** Number of discrete interaction left */
  double fIntLen;      /** Total mean free path i.e. macroscopic cross section */

  ExtraInfo* fExtraInfo;  /** Pointer to an arbitrary struct to keep extra information */

public:
  /** @brief LightTrack default constructor */
  LightTrack();

  /** @brief LightTrack complete constructor */
  LightTrack( const LTrackStatus aTrackStatus, const int aGVcode, const int aGTrackIndex, const int aMaterialCutCoupleIndex,
              const int aProcessIndex, const int aTargetZ,
              const int aTargetN, const double aXdir, const double aYdir, const double aZdir,
              const double aKinE, const double aMass, const double aTime, const double aWeight,
              const double aStepLength, const double aEdep, const double aNintLen, const double aIntLen,
              ExtraInfo* aExtraInfo = 0 );

  /** @brief LightTrack copy constructor */
  LightTrack( const LightTrack &other );

  /** @brief Operator = */
  LightTrack& operator=( const LightTrack &other );

  /** @brief LightTrack destructor */
  ~LightTrack();

  //--- Getters ---

  /** @brief Method that returns the track status */
  LTrackStatus GetTrackStatus() const { return fTrackStatus; }

  /** @brief Method that returns the GV particle code */
  int GetGVcode() const { return fGVcode; }

  /** @brief Method that returns the track index */
  int GetTrackIndex() const { return fGTrackIndex; }

  /** @brief Method that returns the material-cut couple index */
  int GetMaterialCutCoupleIndex() const { return fMaterialCutCoupleIndex; }

  /** @brief Method that returns the process index */
  int GetProcessIndex() const { return fProcessIndex; }

  /** @brief Method that returns the target atomic number, Z */
  int GetTargetZ() const { return fTargetZ; }

  /** @brief Method that returns the target number of nucleons */
  int GetTargetN() const { return fTargetN; }

  /** @brief Method that returns the X direction value (normalized, adimensional) */
  double GetDirX() const { return fXdir; }

  /** @brief Method that returns the Y direction value (normalized, adimensional) */
  double GetDirY() const { return fYdir; }

  /** @brief Method that returns the Z direction value (normalized, adimensional) */
  double GetDirZ() const { return fZdir; }

  /** @brief Method that returns the kinetic energy (unit: energy) */
  double GetKinE() const { return fKinE; }

  /** @brief Method that returns the dynamic mass (unit: energy) */
  double GetMass() const { return fMass; }

  /** Method that returns the global time (unit: time) */
  double GetTime() const { return fTime; }

  /** Method that returns the weight (adimensional) */
  double GetWeight() const { return fWeight; }

  /** Method that returns the true, physical length of the last step (unit: length) */
  double GetStepLength() const { return fStepLength; }

  /** Method that returns the energy deposit in the last step (unit: energy) */
  double GetEnergyDeposit() const { return fEdep; }

  /** Method to get the number of interaction length left. */
  double GetNumOfInteractionLegthLeft() const { return fNintLen; }

  /** Method to get the total mean free path i.e. the macroscopic cross section . */
  double GetTotalMFP() const { return fIntLen; }

  /** Method that returns the pointer to extra information */
  ExtraInfo* GetExtraInfo() const { return fExtraInfo; }

  //--- Setters ---

  /**
   * @brief Method that sets the GV particle code
   * @param aGVcode GV particle code
   */
  void SetGVcode( const int aGVcode ) { fGVcode = aGVcode; }

  /**
   * @brief Method that sets the track index
   * @param aTrackIndex track index
   */
  void SetTrackIndex( const int aTrackIndex ) { fGTrackIndex = aTrackIndex; }

  /**
   * @brief Method that sets the material-cut couple index
   * @param aMaterialCutCoupleIndex material-cut couple index
   */
  void SetMaterialCutCoupleIndex( const int aMaterialCutCoupleIndex ) {
    fMaterialCutCoupleIndex = aMaterialCutCoupleIndex;
  }

  /**
   * @brief Method that sets the track status
   * @param aTrackStatus track status
   */
  void SetTrackStatus( const LTrackStatus aTrackStatus ) { fTrackStatus = aTrackStatus; }

  /**
   * @brief Method that sets the process index
   * @param aProcessIndex process index
   */
  void SetProcessIndex( const int aProcessIndex ) { fProcessIndex = aProcessIndex; }

  /**
   * @brief Method that sets the target atomic number, Z
   * @param aTargetZ target Z
   */
  void SetTargetZ( const int aTargetZ ) { fTargetZ = aTargetZ; }

  /**
   * @brief Method that sets the target number of nucleons
   * @param aTargetN target number of nucleons
   */
  void SetTargetN( const int aTargetN ) { fTargetN = aTargetN; }

  /**
   * @brief Method that sets direction
   * @param aXdir, aYdir, aZdir direction components (normalized, adimensional)
   */
  void SetDirection( const double aXdir, const double aYdir, const double aZdir)
  {
    fXdir = aXdir;
    fYdir = aYdir;
    fZdir = aZdir;
  }

  /**
   * @brief Method that sets the X direction value
   * @param aXdir X direction (normalized, adimensional)
   */
  void SetDirX( const double aXdir ) { fXdir = aXdir; }

  /** @brief Method that sets the Y direction value
   * @param aYdir Y direction (normalized, adimensional)
   */
  void SetDirY( const double aYdir ) { fYdir = aYdir; }

  /** @brief Method that sets the Z direction value
   * @param aZdir Z direction (normalized, adimensional)
   */
  void SetDirZ( const double aZdir ) { fZdir = aZdir; }

  /**
   * @brief Method that sets the kinetic energy
   * @param aKinE kinetic energy (unit: energy)
   */
  void SetKinE( const double aKinE ) { fKinE = aKinE; }

  /**
   * @brief Method that sets the dynamic mass
   * @param aMass dynamic mass (unit: energy)
   */
  void SetMass( const double aMass ) { fMass = aMass; }

  /**
   * @brief Method that sets the global time
   * @param aTime global time (unit: time)
   */
  void SetTime( const double aTime ) { fTime = aTime; }

  /**
   * @brief Method that sets the weight
   * @param aWeight weight (adimensional)
   */
  void SetWeight( const double aWeight ) { fWeight = aWeight; }

  /**
   * @brief Method that sets the true, physical length of the last step
   * @param aStepLength (unit: length)
   */
  void SetStepLength( const double aStepLength ) { fStepLength = aStepLength; }

  /**
   * @brief Method that sets the energy deposit in the last step
   * @param aEdep (unit: energy)
   */
  void SetEnergyDeposit( const double aEdep ) { fEdep = aEdep; }

  /**
   * @brief Method that sets the number of discrete interaction left.
   * @param[in] val Number of discrete interactions.
   */
   void SetN(const double val)  { fIntLen = val; }

   /**
    * @brief Method that sets the number of interaction length left .
    * @param[in] val Number of interaction length left.
    */
   void SetNumOfInteractionLegthLeft(const double val) { fNintLen = val; }

  /**
   * @brief Method that sets the total mean free path i.e. the macroscopic cross section .
   * @param[in] val Value of the total mean free path in internal [1/length] unit.
   */
   void SetTotalMFP(const double val)  { fIntLen = val; }

  /**
   * @brief Method that sets the pointer to extra information
   * @param aExtraInfo pointer to extra information
   */
  void SetExtraInfo( ExtraInfo* aExtraInfo ) {
    // If the pointer is already set, and it is different from the new one,
    // then delete the old extra information before pointing to the new one.
    // Note: we assume that the light track has the ownership of the extra information.
    if ( fExtraInfo != nullptr  &&  fExtraInfo != aExtraInfo ) delete fExtraInfo;
    fExtraInfo = aExtraInfo;
  }

};


/**
 * @brief SOA (Structure Of Arrays) for LightTrack
 *
 */
class LightTrack_v {
private:
  std::atomic_int fNtracks;  /** number of tracks contained */

  LTrackStatus *fTrackStatusV;    /** Status of the tracks */
  int *fGVcodeV;                  /** GV particle codes */
  int *fGTrackIndexV;             /** Indices of the tracks in the GeantTrack_v vector */
  int *fMaterialCutCoupleIndexV;  /** Indices of the material-cut couples where the tracks are */
  int *fProcessIndexV;            /** Indices of the selected processes (in the list of
                                      active discrete or at-rest processes kept by the
                                      PhysicsManagerPerParticle; negative if the winner
                                      is a continuous process) */
  int *fTargetZV;                 /** Atomic numbers (Z) of the target atoms where the
                                      particles have a physical interaction */
  int *fTargetNV;                 /** Numbers of nucleons of the target atoms where the
                                      particles have a physical interaction */

  double *fXdirV;        /** X directions (normalized, adimensional) of the particles */
  double *fYdirV;        /** Y directions (normalized, adimensional) of the particles */
  double *fZdirV;        /** Z directions (normalized, adimensional) of the particles */
  double *fKinEV;        /** Kinetic energies (in GeV) of the particles */
  double *fMassV;        /** Dynamic masses (in GeV) of the particles */
  double *fTimeV;        /** Times (global, in sec) of the particles */
  double *fWeightV;      /** Weights (adimensional) of the particles */
  double *fStepLengthV;  /** True, physical lengths (in cm) of the last step of the particles */
  double *fEdepV;        /** Energy deposits (in GeV) in the last step of the particles */
  double *fNintLenV;     /** Number of discrete interaction left */
  double *fIntLenV;      /** Total mean free path i.e. macroscopic scross section */

  ExtraInfo* *fExtraInfoV;  /** Pointers to arbitrary structs to keep extra information */

private:
  /** @brief Copy constructor (not allowed) */
  LightTrack_v( const LightTrack_v &track_v );

  /** @brief Operator= (not allowed) */
  LightTrack_v& operator=( const LightTrack_v &track_v );

public:
  /** @brief LightTrack_v constructor */
  LightTrack_v();

  /** @brief LightTrack_v destructor */
  ~LightTrack_v();

  /** @brief Method that returns the number of tracks contained */
  int GetNtracks() const { return fNtracks; }

  /** @brief Method that set number of tracks contained */
  void SetNtracks( const int ntracks ) { fNtracks = ntracks; }

  /**
   * @brief Method that returns a track
   * @param i Bit number 'i'
   * @param aLightTrack track to be returned
   */
  void GetTrack( const int i, LightTrack &aLightTrack ) const;

  /**
   * @brief Method that add a track
   * @param aLightTrack track that should be added
   */
  void AddTrack( LightTrack &aLightTrack );

};

}  // end of namespace geantphysics

#endif
