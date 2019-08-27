
#ifndef MyRunDataPerPrimary_h
#define MyRunDataPerPrimary_h 1

#include "globals.hh"

class MyRunDataPerPrimary {
public:
  MyRunDataPerPrimary();
  ~MyRunDataPerPrimary();

  void Clear();

  MyRunDataPerPrimary& operator+=(const MyRunDataPerPrimary& other);

  G4double fNumPrimaries;   // total number of primary particles
  G4double fEdep;           // sum of energy deposit (per event)
  G4double fEdep2;          // sum of energy deposit square
  G4double fTrackLCh;       // sum of charged step length (per event)
  G4double fTrackLCh2;      // sum of charged step length square
  G4double fTrackLNe;       // sum of neutral step length (per event)
  G4double fTrackLNe2;      // sum of neutral step length square

  G4double fTrackTCh;       // sum of charged track (global) time (per event)
  G4double fTrackTCh2;      // sum of charged track (global) time square
  G4double fTrackTNe;       // sum of neutral track (global) time (per event)
  G4double fTrackTNe2;      // sum of neutral track (global) time square

  G4double fChargedStep;    // sum of number of charged steps (per event)
  G4double fChargedStep2;   // sum of number of charged steps square
  G4double fNeutralStep;    // sum of number of neutral steps (per event)
  G4double fNeutralStep2;   // sum of number of neutral steps  square

  G4double fNGamma;         // sum of number of secondary gamma (per event)
  G4double fNGamma2;        // sum of number of secondary gamma square
  G4double fNElec;          // sum of number of secondary e- (per event)
  G4double fNElec2;         // sum of number of secondary e- square
  G4double fNPosit;         // sum of number of secondary e+ (per event)
  G4double fNPosit2;        // sum of number of secondary e+ square
};

#endif
