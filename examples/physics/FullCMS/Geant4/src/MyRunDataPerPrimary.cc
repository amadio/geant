
#include "MyRunDataPerPrimary.hh"


MyRunDataPerPrimary::MyRunDataPerPrimary() {
  Clear();
}


MyRunDataPerPrimary::~MyRunDataPerPrimary() {}


void MyRunDataPerPrimary::Clear() {
  fNumPrimaries = 0.;
  fEdep         = 0.;
  fEdep2        = 0.;
  fTrackLCh     = 0.;
  fTrackLCh2    = 0.;
  fTrackLNe     = 0.;
  fTrackLNe2    = 0.;

  fTrackTCh     = 0.;
  fTrackTCh2    = 0.;
  fTrackTNe     = 0.;
  fTrackTNe2    = 0.;

  fChargedStep  = 0.;
  fChargedStep2 = 0.;
  fNeutralStep  = 0.;
  fNeutralStep2 = 0.;

  fNGamma       = 0.;
  fNGamma2      = 0.;
  fNElec        = 0.;
  fNElec2       = 0.;
  fNPosit       = 0.;
  fNPosit2      = 0.;
}


MyRunDataPerPrimary& MyRunDataPerPrimary::operator+=(const MyRunDataPerPrimary& right) {
  fNumPrimaries += right.fNumPrimaries;
  fEdep         += right.fEdep;
  fEdep2        += right.fEdep2;
  fTrackLCh     += right.fTrackLCh;
  fTrackLCh2    += right.fTrackLCh2;
  fTrackLNe     += right.fTrackLNe;
  fTrackLNe2    += right.fTrackLNe2;

  fTrackTCh     += right.fTrackTCh;
  fTrackTCh2    += right.fTrackTCh2;
  fTrackTNe     += right.fTrackTNe;
  fTrackTNe2    += right.fTrackTNe2;

  fChargedStep  += right.fChargedStep;
  fChargedStep2 += right.fChargedStep2;
  fNeutralStep  += right.fNeutralStep;
  fNeutralStep2 += right.fNeutralStep2;

  fNGamma       += right.fNGamma;
  fNGamma2      += right.fNGamma2;
  fNElec        += right.fNElec;
  fNElec2       += right.fNElec2;
  fNPosit       += right.fNPosit;
  fNPosit2      += right.fNPosit2;

  return *this;
}
