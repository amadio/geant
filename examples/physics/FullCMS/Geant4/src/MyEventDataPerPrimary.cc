
#include "MyEventDataPerPrimary.hh"

#include "G4SystemOfUnits.hh"

MyEventDataPerPrimary::MyEventDataPerPrimary() {
  fEdep         = 0.;
  fTrackLCh     = 0.;
  fTrackLNe     = 0.;
  fChargedStep  = 0;
  fNeutralStep  = 0;
  fNGamma       = 0.;
  fNElec        = 0.;
  fNPosit       = 0.;
}


MyEventDataPerPrimary::~MyEventDataPerPrimary() {}


void MyEventDataPerPrimary::Clear() {
  fEdep         = 0.;
  fTrackLCh     = 0.;
  fTrackLNe     = 0.;
  fChargedStep  = 0 ;
  fNeutralStep  = 0 ;
  fNGamma       = 0.;
  fNElec        = 0.;
  fNPosit       = 0.;
}


std::ostream& operator<<(std::ostream& flux, const MyEventDataPerPrimary& evtdata) {
  std::ios::fmtflags mode = flux.flags();
  flux.setf(std::ios::fixed,std::ios::floatfield);
  long prec = flux.precision(3);
  flux << "    .............................................................  \n"
       << "    Edep                   = " << evtdata.fEdep/GeV     << " [GeV] \n"
       << "    Track lenght (charged) = " << evtdata.fTrackLCh/cm  << "  [cm] \n"
       << "    Track lenght (neutral) = " << evtdata.fTrackLNe/cm  << "  [cm] \n"
       << "    Steps (charged)        = " << evtdata.fChargedStep  << "       \n"
       << "    Steps (neutral)        = " << evtdata.fNeutralStep  << "       \n"
       << "    Secondary Gammas       = " << evtdata.fNGamma       << "       \n"
       << "    Secondary Electrons    = " << evtdata.fNElec        << "       \n"
       << "    Secondary Positrons    = " << evtdata.fNPosit       << "       \n";
  flux.precision(prec);
  flux.setf(mode,std::ios::floatfield);
  return flux;
}
