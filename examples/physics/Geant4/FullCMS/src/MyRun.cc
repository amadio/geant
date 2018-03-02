
#include "MyRun.hh"

#include "MyEventDataPerPrimary.hh"
#include "MyPrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>


MyRun::MyRun() : G4Run() {
  fNumPrimaryTypes = MyPrimaryGeneratorAction::GetNumberOfPrimaryTypes();
  fRunDataPerPrimary.resize(fNumPrimaryTypes);
  for (G4int ipt=0; ipt<fNumPrimaryTypes; ++ipt) {
    fRunDataPerPrimary[ipt].Clear();
  }
}


MyRun::~MyRun() {
  fRunDataPerPrimary.clear();
}


void MyRun::FillPerEvent(const MyEventDataPerPrimary& data, G4int primtypeindx) {
  MyRunDataPerPrimary& runData = fRunDataPerPrimary[primtypeindx];
  runData.fEdep        += data.fEdep;         runData.fEdep2        += data.fEdep*data.fEdep;
  runData.fTrackLCh    += data.fTrackLCh;     runData.fTrackLCh2    += data.fTrackLCh*data.fTrackLCh;
  runData.fTrackLNe    += data.fTrackLNe;     runData.fTrackLNe2    += data.fTrackLNe*data.fTrackLNe;
  runData.fChargedStep += data.fChargedStep;  runData.fChargedStep2 += data.fChargedStep*data.fChargedStep;
  runData.fNeutralStep += data.fNeutralStep;  runData.fNeutralStep2 += data.fNeutralStep*data.fNeutralStep;
  runData.fNGamma      += data.fNGamma;       runData.fNGamma2      += data.fNGamma*data.fNGamma;
  runData.fNElec       += data.fNElec;        runData.fNElec2       += data.fNElec*data.fNElec;
  runData.fNPosit      += data.fNPosit;       runData.fNPosit2      += data.fNPosit*data.fNPosit;
}


void MyRun::Merge(const G4Run* run) {
  const MyRun* localRun = static_cast<const MyRun*>(run);
  if (localRun) {
    for (G4int ipt=0; ipt<fNumPrimaryTypes; ++ipt) {
      fRunDataPerPrimary[ipt] += localRun->GetRunDataPerPrimary(ipt);
    }
  }
  G4Run::Merge(run);
}


void MyRun::EndOfRun() {
  G4int  numEvents    = GetNumberOfEvent();
  G4int  numPrimaries = 0;
  for (G4int ipt=0; ipt<fNumPrimaryTypes; ++ipt) {
    numPrimaries += fRunDataPerPrimary[ipt].fNumPrimaries;
  }
  std::ios::fmtflags mode = G4cout.flags();
  G4int  prec = G4cout.precision(2);
  G4cout<< " \n ==================================   Run summary   ===================================== \n" << G4endl;
  G4cout<< std::setprecision(4);
  G4cout<< "    Number of events        = " << numEvents                                                     << G4endl;
  G4cout<< "    Total number of primary = " << numPrimaries                                                  << G4endl;
  G4cout<< " \n ---------------------------------------------------------------------------------------- \n" << G4endl;
  // compute and print run statistics per primary type per primary
  for (G4int ipt=0; ipt<fNumPrimaryTypes; ++ipt) {
    MyRunDataPerPrimary& runData = fRunDataPerPrimary[ipt];
    G4int     nPrimaries = runData.fNumPrimaries;
    G4double  norm = static_cast<G4double>(nPrimaries);
    if (norm>0.) {
      norm = 1./norm;
    } else {
      continue;
    }
    //compute and print statistic
    //
    G4String primName   = MyPrimaryGeneratorAction::GetPrimaryName(ipt);
    G4double meanEdep   = runData.fEdep*norm;
    G4double rmsEdep    = std::sqrt(std::abs(runData.fEdep2*norm-meanEdep*meanEdep));
    G4double meanLCh    = runData.fTrackLCh*norm;
    G4double rmsLCh     = std::sqrt(std::abs(runData.fTrackLCh2*norm-meanLCh*meanLCh));
    G4double meanLNe    = runData.fTrackLNe*norm;
    G4double rmsLNe     = std::sqrt(std::abs(runData.fTrackLNe2*norm-meanLNe*meanLNe));
    G4double meanStpCh  = runData.fChargedStep*norm;
    G4double rmsStpCh   = std::sqrt(std::abs(runData.fChargedStep2*norm-meanStpCh*meanStpCh));
    G4double meanStpNe  = runData.fNeutralStep*norm;
    G4double rmsStpNe   = std::sqrt(std::abs(runData.fNeutralStep2*norm-meanStpNe*meanStpNe));
    G4double meanNGam   = runData.fNGamma*norm;
    G4double rmsNGam    = std::sqrt(std::abs(runData.fNGamma2*norm-meanNGam*meanNGam));
    G4double meanNElec  = runData.fNElec*norm;
    G4double rmsNElec   = std::sqrt(std::abs(runData.fNElec2*norm-meanNElec*meanNElec));
    G4double meanNPos   = runData.fNPosit*norm;
    G4double rmsNPos    = std::sqrt(std::abs(runData.fNPosit2*norm-meanNPos*meanNPos));

    G4cout<< "  Number of primaries        = " << nPrimaries  << "  " << primName                             <<G4endl;
    G4cout<< "  Total energy deposit per primary = " << meanEdep/GeV <<  " +- " << rmsEdep/GeV << " [GeV]"    <<G4endl;
    G4cout<< G4endl;
    G4cout<< "  Total track length (charged) per primary = " << meanLCh/cm << " +- " << rmsLCh/cm <<  " [cm]" <<G4endl;
    G4cout<< "  Total track length (neutral) per primary = " << meanLNe/cm << " +- " << rmsLNe/cm <<  " [cm]" <<G4endl;
    G4cout<< G4endl;
    G4cout<< "  Number of steps (charged) per primary = " << meanStpCh << " +- " << rmsStpCh << G4endl;
    G4cout<< "  Number of steps (neutral) per primary = " << meanStpNe << " +- " << rmsStpNe << G4endl;
    G4cout<< G4endl;
    G4cout<< "  Number of secondaries per primary : " << G4endl
          << "     Gammas    =  " << meanNGam      <<  " +- " << rmsNGam  << G4endl
          << "     Electrons =  " << meanNElec     <<  " +- " << rmsNElec << G4endl
          << "     Positrons =  " << meanNPos      <<  " +- " << rmsNPos  << G4endl;
    G4cout<< " ......................................................................................... \n" << G4endl;
  }
  G4cout<< " \n ======================================================================================== \n" << G4endl;

  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
}
