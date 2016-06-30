#ifndef GEANT_LHCbApp
#define GEANT_LHCbApp

#include <map>

#include "globals.hh"

class G4Step;
class TH1F;
class G4LogicalVolume;

/** @brief LHCbApp class */
class LHCbApp {
  static const G4int kNvolumes     = 4900;
  static const G4int kNECALModules = 30;
  static const G4int kNHCALModules = 30;

private:
  G4bool  fInitialized;                           /** Initialized flag */
  std::map<G4LogicalVolume*,G4int> fLogicVolToIndex;       /** Map G4LogicalVolume* -> its index in the store */
  G4int necal;
  G4int nhcal;

  G4bool  fSensFlags[kNvolumes];                  /** Array marking sensitive volumes */
  G4double fEdepECAL[kNECALModules];              /** Energy deposition in ECAL */
  G4double fEdepHCAL[kNHCALModules];              /** Energy deposition in HCAL */
  G4int   fECALid[kNECALModules];                 /** ECAL volume id's */
  G4int   fHCALid[kNHCALModules];                 /** HCAL volume id's */
  G4double fCubicVolumes[kNvolumes];              /** Precomputed volume capacity */ 
  std::map<int,int> fECALMap;                     /** Map of ECAL modules */
  std::map<int,int> fHCALMap;                     /** Map of ECAL modules */
  TH1F   *fFluxElec;                              /** Flux histogram for electrons */
  TH1F   *fFluxGamma;                             /** Flux histogram for gammas */
  TH1F   *fFluxP;                                 /** Flux histogram for protons */
  TH1F   *fFluxPi;                                /** Flux histogram for pions */
  TH1F   *fFluxK;                                 /** Flux histogram for kaons */
  TH1F   *fEdepElec;                              /** Edep histogram for electrons */
  TH1F   *fEdepGamma;                             /** Edep histogram for gammas */
  TH1F   *fEdepP;                                 /** Edep histogram for protons */
  TH1F   *fEdepPi;                                /** Edep histogram for pions */
  TH1F   *fEdepK;                                 /** Edep histogram for kaons */ 

  /**
   * @brief Copy constructor LHCbApp
   * * @todo Still not implemented
   */
  LHCbApp(const LHCbApp &);

  /**
   * @brief Operator=
   * @todo Still not implemented
   */
  LHCbApp &operator=(const LHCbApp &);
public:

  /** @brief Constructor LHCbApp */
  LHCbApp();

  /** @brief Destructor LHCbApp */
  ~LHCbApp() {}

  /**
   * @brief Function of initialization
   */
  G4bool Initialize();

  /**
   * @brief Step by step to do
   */
  void SteppingAction(const G4Step *);

  /**
   * @brief End of run to do
   */
  void EndOfRunAction(G4int);

public:
  static G4bool fgIsScoreActive; 

};
#endif

