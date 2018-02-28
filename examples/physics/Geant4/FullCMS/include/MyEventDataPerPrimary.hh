
#ifndef MyEventDataPerPrimary_h
#define MyEventDataPerPrimary_h 1

#include "globals.hh"
#include <iostream>

class MyEventDataPerPrimary {

public:
  MyEventDataPerPrimary();
  ~MyEventDataPerPrimary();

  void Clear();

  friend std::ostream &operator<<(std::ostream &, const MyEventDataPerPrimary &);

  G4double fEdep;        // sum of energy deposit
  G4double fTrackLCh;    // sum of charged step length
  G4double fTrackLNe;    // sum of neutral step length
  G4double fChargedStep; // sum of number of charged steps
  G4double fNeutralStep; // sum of number of neutral steps

  G4double fNGamma; // sum of number of secondary gamma
  G4double fNElec;  // sum of number of secondary e-
  G4double fNPosit; // sum of number of secondary e+
};

#endif
