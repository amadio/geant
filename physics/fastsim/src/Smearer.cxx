#include "Smearer.h"

#include "GeantScheduler.h"
#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "WorkloadManager.h"
#include "globals.h"
#include "GeantTrackVec.h"

#ifdef USE_VECGEOM_NAVIGATOR
  #include "ParticleOld.h"
  #include "management/GeoManager.h"
  #include "navigation/NavigationState.h"
  #include "base/Global.h"
  #include "base/RNG.h"
#else
  #include "TGeoManager.h"
  #include "TGeoMaterial.h"
  #include "TGeoBranchArray.h"
  #include "TGeoExtension.h"
  #include "TList.h"
#endif

#ifdef USE_ROOT
  #include "Rtypes.h"
  #include "TRandom.h"
  #include "TFile.h"
  #include "TBits.h"
  #include "TError.h"
  #include "TSystem.h"
  #include "TDatabasePDG.h"
#endif
#include <strings.h>


Smearer::Smearer(GeantTaskData */*td*/) {

  std::cout << "Smearer::Smearer : Start" << std::endl;  // Debug

  #ifdef USE_VECGEOM_NAVIGATOR
//    td->fPropagator->LoadVecGeomGeometry();
    std::vector< vecgeom::LogicalVolume* > vecgeomVolumes;
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes( vecgeomVolumes );
    int numberOfVolumes = vecgeomVolumes.size();
  #else
    TObjArray* allVolumes = gGeoManager->GetListOfVolumes();
    int numberOfVolumes = allVolumes->GetEntries();
  #endif
  std::cout << "\t numberOfVolumes=" << numberOfVolumes << std::endl;  // Debug

  hasVolumeTrackerParameterisation.reserve( numberOfVolumes+1 );
  hasVolumeEcalParameterisation.reserve( numberOfVolumes+1 );
  hasVolumeHcalParameterisation.reserve( numberOfVolumes+1 );
  hasVolumeMuonParameterisation.reserve( numberOfVolumes+1 );
  hasVolumeParameterisation.reserve( numberOfVolumes+1 );
  
  for ( auto ivol = 0; ivol < numberOfVolumes; ivol++ ) {

    int idVol = 0;
    #ifdef USE_VECGEOM_NAVIGATOR
      vecgeom::LogicalVolume* aVolume = vecgeomVolumes[ ivol ];
      if ( aVolume ) idVol = aVolume->id();
    #else
      TGeoVolume* aVolume = (TGeoVolume*) allVolumes->At( ivol );
      if ( aVolume ) idVol = aVolume->GetNumber();
    #endif
      const char* nameVolume = ( aVolume ) ? aVolume->GetName() : "";

    int flagTracker = 0;
    int flagEcal = 0;
    int flagHcal = 0;
    int flagMuon = 0;
    if (strcasecmp(nameVolume, "tracker")>0 ) {
      flagTracker = 1;
    } else if (strcasecmp(nameVolume, "ecal")>0 ) {
      flagEcal = 1;
    } else if (strcasecmp(nameVolume, "hcal")>0 ) {
      flagHcal = 1;
    } else if (strcasecmp(nameVolume, "muon")>0 ) {
    //flagMuon = 1;  // No muon parameterisation
    }
    hasVolumeTrackerParameterisation[ idVol ] = flagTracker;
    hasVolumeEcalParameterisation[ idVol ] = flagEcal;
    hasVolumeHcalParameterisation[ idVol ] = flagHcal;
    hasVolumeMuonParameterisation[ idVol ] = flagMuon;
    hasVolumeParameterisation[ idVol ] = 0;
    if ( flagTracker || flagEcal || flagHcal || flagMuon ) {
      hasVolumeParameterisation[ idVol ] = 1;
    }

    std::cout << "\t name=" << nameVolume << " ; id=" << idVol
              << " ; (my-counter=" << ivol << ")"
              << " ; flagTracker=" << flagTracker << " ; flagEcal=" << flagEcal 
              << " ; flagHcal=" << flagHcal << " ; flagMuon=" << flagMuon 
              << "  -> " << hasVolumeParameterisation[ idVol ] << std::endl;  // Debug
  }

  std::cout << "Smearer::Smearer : --- End ---" << std::endl;  // Debug
}


std::vector< double > Smearer::StepLengthProposedByParameterisation( int ntracks, GeantTrack_v& tracks, GeantTaskData& aTaskData ) {
  //std::cout << "Smearer::StepLengthProposedByParameterisation : Start" << std::endl;  // Debug
  const double aLargeValue = 9999999.9;
  std::vector< double > stepLengths;
  GeantTrack& aTrack = aTaskData.GetNewTrack();
  for ( int i = 0; i < ntracks; i++ ) {
    double value = aLargeValue;
     tracks.GetTrack(i, aTrack);
    // Check in which part of the detector (Tracker, ECAL, HCAL, Muon system)
    // the track is, and then check if a parameterisation should be applied
    // for that track in that volume.
    #ifdef USE_VECGEOM_NAVIGATOR
      int idVol = aTrack.GetVolume()->id();
    #else
      int idVol = aTrack.GetVolume()->GetNumber();
    #endif
    //std::cout << " Smearer::StepLengthProposedByParameterisation : pdg=" << aTrack.PDG() 
    //          << " ; volume=" << aTrack.GetVolume()->GetName() << std::endl;  // Debug
    if ( hasVolumeParameterisation[ idVol ] ) {
      if ( hasVolumeTrackerParameterisation[ idVol ] ) {
        value = StepLengthProposedByTrackerParameterisation( aTrack );
      } else if ( hasVolumeEcalParameterisation[ idVol ] ) { 
        value = StepLengthProposedByEcalParameterisation( aTrack );
      } else if ( hasVolumeHcalParameterisation[ idVol ] ) { 
        value = StepLengthProposedByHcalParameterisation( aTrack );
      } else if ( hasVolumeMuonParameterisation[ idVol ] ) { 
        value = StepLengthProposedByMuonParameterisation( aTrack );
      }
    }
    stepLengths.push_back( value );
  }
  aTaskData.ReleaseTrack(aTrack);
  //std::cout << "Smearer::StepLengthProposedByParameterisation : --- End ---" << std::endl;  // Debug
  return stepLengths;
}


double Smearer::StepLengthProposedByTrackerParameterisation( GeantTrack& aTrack ) {
  // The tracker parameterisation should be applied for all primary charged particles.
  // Given that currently there is no way to check whether a track is a primary,
  // for the time being we apply it to all charged particles.
  // The proposed step length which is returned is a large value, because the
  // parameterisation must be applied at the end (exit) of the tracker.
  //std::cout << "Smearer::StepLengthProposedByTrackerParameterisation : Start" << std::endl;  // Debug
  const double aLargeValue = 9999999.9;
  double proposedLength = aLargeValue;
  if ( IsTrackerParameterisationApplicable( aTrack ) ) {
    proposedLength = aLargeValue;
  }
  //std::cout << " Smearer::StepLengthProposedByTrackerParameterisation=" << proposedLength << std::endl
  //          << "Smearer::StepLengthProposedByTrackerParameterisation : --- End ---" << std::endl;  // Debug
  return proposedLength;
}


double Smearer::StepLengthProposedByEcalParameterisation( GeantTrack& aTrack ) {
  // The ECAL parameterisation should be applied for all primary electrons, 
  // positrons and gammas.
  // Given that currently there is no way to check whether a track is a primary,
  // for the time being we apply it to all electrons, positrons and gammas. 
  // The proposed step length which is returned is 0.0, because the
  // parameterisation must be applied at the entrance of the EM calorimeter.
  //std::cout << "Smearer::StepLengthProposedByEcalParameterisation : Start" << std::endl;  // Debug
  const double aLargeValue = 9999999.9;
  double proposedLength = aLargeValue;
  if ( IsEcalParameterisationApplicable( aTrack ) ) {
    proposedLength = 0.0;
  }
  //std::cout << " Smearer::StepLengthProposedByEcalParameterisation=" << proposedLength << std::endl
  //          << "Smearer::StepLengthProposedByEcalParameterisation : --- End ---" << std::endl;  // Debug
  return proposedLength;
}


double Smearer::StepLengthProposedByHcalParameterisation( GeantTrack& aTrack ) {
  // The HCAL parameterisation should be applied for all primary hadrons.
  // Given that currently there is no way to check whether a track is a primary,
  // for the time being we apply it to all hadrons. 
  // The proposed step length which is returned is 0.0, because the
  // parameterisation must be applied at the entrance of the HAD calorimeter.
  //std::cout << "Smearer::StepLengthProposedByHcalParameterisation : Start" << std::endl;  // Debug
  const double aLargeValue = 9999999.9;
  double proposedLength = aLargeValue;
  if ( IsHcalParameterisationApplicable( aTrack ) ) {
    proposedLength = 0.0;
  }
  //std::cout << " Smearer::StepLengthProposedByHcalParameterisation=" << proposedLength << std::endl
  //          << "Smearer::StepLengthProposedByHcalParameterisation : --- End ---" << std::endl;  // Debug
  return proposedLength;
}


double Smearer::StepLengthProposedByMuonParameterisation( GeantTrack& aTrack ) {
  // The muon parameterisation (in principle) should be applied for all primary muons
  // (but we don't have any parameterisation for muons).
  // Given that currently there is no way to check whether a track is a primary,
  // for the time being we apply it (in principle) to all muons. 
  // The proposed step length which is returned should be either 0.0, if the
  // parameterisation must be applied at the entrance of the muon sub-detector,
  // or a large value, if it should be applied at the exit.
  //std::cout << "Smearer::StepLengthProposedByMuonParameterisation : Start" << std::endl;  // Debug
  const double aLargeValue = 9999999.9;
  double proposedLength = aLargeValue;
  if ( IsMuonParameterisationApplicable( aTrack ) ) {
    //proposedLength = 0.0;
  }
  //std::cout << " Smearer::StepLengthProposedByMuonParameterisation=" << proposedLength << std::endl
  //          << "Smearer::StepLengthProposedByMuonParameterisation : --- End ---" << std::endl;  // Debug
  return proposedLength;
}


bool Smearer::IsTrackerParameterisationApplicable( GeantTrack& aTrack ) {
  // The tracker parameterisation should be applied for all primary charged particles.
  // Given that currently there is no way to check whether a track is a primary,
  // for the time being we apply it to all charged particles.
  //std::cout << "Smearer::IsTrackerParameterisationApplicable : Start" << std::endl;  // Debug
  bool isApplicable = false;
  #ifdef USE_VECGEOM_NAVIGATOR
    int idVol = aTrack.GetVolume()->id();
  #else
    int idVol = aTrack.GetVolume()->GetNumber();
  #endif
  if ( hasVolumeTrackerParameterisation[ idVol ] ) {
    if ( aTrack.Charge() ) {  // any charged particle
      isApplicable = true;
    }
  }
  //std::cout << " Smearer::IsTrackerParameterisationApplicable=" << isApplicable << std::endl
  //          << "Smearer::IsTrackerParameterisationApplicable : --- End ---" << std::endl;  // Debug
  return isApplicable;
}


bool Smearer::IsEcalParameterisationApplicable( GeantTrack& aTrack ) {
  // The ECAL parameterisation should be applied for all primary electrons, 
  // positrons and gammas.
  // Given that currently there is no way to check whether a track is a primary,
  // for the time being we apply it to all electrons, positrons and gammas. 
  //std::cout << "Smearer::IsEcalParameterisationApplicable : Start" << std::endl;  // Debug
  bool isApplicable = false;
  #ifdef USE_VECGEOM_NAVIGATOR
    int idVol = aTrack.GetVolume()->id();
  #else
    int idVol = aTrack.GetVolume()->GetNumber();
  #endif
  if ( hasVolumeEcalParameterisation[ idVol ] ) {
    int absPDG = std::abs( aTrack.PDG() );
    if ( absPDG == 11 || absPDG == 22 ) {  // electron, positron and gamma
      isApplicable = true;
    }
  }
  //std::cout << " Smearer::IsEcalParameterisationApplicable=" << isApplicable << std::endl
  //          << "Smearer::IsEcalParameterisationApplicable : --- End ---" << std::endl;  // Debug
  return isApplicable;
}


bool Smearer::IsHcalParameterisationApplicable( GeantTrack& aTrack ) {
  // The HCAL parameterisation should be applied for all primary hadrons.
  // Given that currently there is no way to check whether a track is a primary,
  // for the time being we apply it to all hadrons. 
  //std::cout << "Smearer::IsHcalParameterisationApplicable : Start" << std::endl;  // Debug
  bool isApplicable = false;
  #ifdef USE_VECGEOM_NAVIGATOR
    int idVol = aTrack.GetVolume()->id();
  #else
    int idVol = aTrack.GetVolume()->GetNumber();
  #endif
  if ( hasVolumeHcalParameterisation[ idVol ] ) {
    int absPDG = std::abs( aTrack.PDG() );
    if ( 
         absPDG == 111  ||  absPDG == 211  ||                     // pions
         absPDG == 113  ||  absPDG == 213  ||                     // rho
         absPDG == 221  ||  absPDG == 331  ||                     // eta, eta_prime
         absPDG == 223  ||  absPDG == 333  ||                     // omega and phi
         absPDG == 130  ||  absPDG == 310  ||  absPDG == 321  ||  // kaon 
         absPDG == 411  ||  absPDG == 421  ||                     // D-meson
         absPDG == 511  ||  absPDG == 521  ||  absPDG == 531  ||  // B-meson
         absPDG == 2212 ||  absPDG == 2112 ||                     // nucleon
         ( absPDG >= 3122  &&  absPDG <= 5554 )                   // hyperon, charm- and bottom-baryons 
       ) {
      isApplicable = true;     
    }
  }
  //std::cout << " Smearer::IsHcalParameterisationApplicable=" << isApplicable << std::endl
  //          << "Smearer::IsHcalParameterisationApplicable : --- End ---" << std::endl;  // Debug
  return isApplicable;
}


bool Smearer::IsMuonParameterisationApplicable( GeantTrack& aTrack ) {
  // The muon parameterisation (in principle) should be applied for all primary muons
  // (but we don't have any parameterisation for muons).
  // Given that currently there is no way to check whether a track is a primary,
  // for the time being we apply it (in principle) to all muons. 
  //std::cout << "Smearer::IsMuonParameterisationApplicable : Start" << std::endl;  // Debug
  bool isApplicable = false;
  #ifdef USE_VECGEOM_NAVIGATOR
    int idVol = aTrack.GetVolume()->id();
  #else
    int idVol = aTrack.GetVolume()->GetNumber();
  #endif
  if ( hasVolumeMuonParameterisation[ idVol ] ) {
    int absPDG = std::abs( aTrack.PDG() );
    if ( absPDG == 13 ) {  // muon
      isApplicable = true;
    }
  }
  //std::cout << " Smearer::IsMuonParameterisationApplied=" << isApplicable << std::endl
  //          << "Smearer::IsMuonParameterisationApplicable : --- End ---" << std::endl;  // Debug
  return isApplicable;
}


void Smearer::ApplyParameterisation( int ntracks, GeantTrack_v& tracks, GeantTaskData& aTaskData, bool isElossCalling ) {
  // Apply the parameterisation, if it is the case, by calling the proper 
  // parameterisation; else do nothing.
  //std::cout << "Smearer::ApplyParameterisation : Start" << std::endl;  // Debug
  const double aSmallValue = 1.0e-9;
  GeantTrack& aTrack = aTaskData.GetNewTrack();
  for ( int i = 0; i < ntracks; i++ ) {
    //std::cout << " Smearer::ApplyParameterisation : pdg=" << tracks.fPDGV[i] 
    //          << " ; volume=" << tracks.GetVolume(i)->GetName() << std::endl;  // Debug
    // Kill particles with low Pt (< 1 MeV) or large eta (> 5.5 in module).
    // We put here this check because this method is called for all particles.
    double p = tracks.fPV[i];
    double pz = p*tracks.fZdirV[i];
    double pt = std::sqrt( p*p - pz*pz );
    double eta = 0.5 * std::log( (p + pz)/(p - pz) );
    //std::cout << "\t \t p=" << p << " ; pz=" << pz << " ; pt=" << pt << " GeV; eta=" << eta << std::endl;  // Debug 
    if ( pt < 0.001 || std::abs( eta ) > 5.5 ) {
      tracks.fStatusV[i] = Geant::kKilled;
      //std::cout << " Smearer::ApplyParameterisation : KILLED track with  pt=" << pt 
      //          << " ; eta=" << eta << std::endl;  // Debug
      continue;
    }
    // Skip it if it is called from "Eloss" but the proposed step-length is zero.
    if ( isElossCalling && tracks.fPstepV[i] < aSmallValue ) {
      //std::cout << " Smearer::ApplyParameterisation : SKIP, because isElossCalling="
      //          << isElossCalling << " && tracks.fPstepV[i]=" << tracks.fPstepV[i] 
      //          << " < " << aSmallValue << std::endl;  // Debug
      continue;
    }
    tracks.fProcessV[i] = 999;  // Arbitrary number to signal fast-simulation process
    //std::cout << " Smearer::ApplyParameterisation : isElossCalling=" << isElossCalling
    //          << " ; tracks.fPstepV[i]=" << tracks.fPstepV[i] << std::endl;  // Debug
    tracks.GetTrack(i, aTrack);
    // Check in which part of the detector (Tracker, ECAL, HCAL, Muon system)
    // the track is, and then check if a parameterisation should be applied
    // for that track in that volume.
    #ifdef USE_VECGEOM_NAVIGATOR
      int idVol = aTrack.GetVolume()->id();
    #else
      int idVol = aTrack.GetVolume()->GetNumber();
    #endif
    if ( hasVolumeParameterisation[ idVol ] ) {
      if ( hasVolumeTrackerParameterisation[ idVol ] ) {
        if ( IsTrackerParameterisationApplicable( aTrack ) ) {
          ApplyTrackerParameterisation( tracks, aTaskData, i );
        }
      } else if ( hasVolumeEcalParameterisation[ idVol ] ) { 
        if ( IsEcalParameterisationApplicable( aTrack ) ) {
          ApplyEcalParameterisation( tracks, aTaskData, i );
        }
      } else if ( hasVolumeHcalParameterisation[ idVol ] ) { 
        if ( IsHcalParameterisationApplicable( aTrack ) ) {
          ApplyHcalParameterisation( tracks, aTaskData, i );
        }
      } else if ( hasVolumeMuonParameterisation[ idVol ] ) { 
        if ( IsMuonParameterisationApplicable( aTrack ) ) {
          ApplyMuonParameterisation( tracks, aTaskData, i );
        }
      }
    }
  }
  aTaskData.ReleaseTrack(aTrack);
  //std::cout << "Smearer::ApplyParameterisation : --- End ---" << std::endl;  // Debug
}


void Smearer::ApplyTrackerParameterisation( GeantTrack_v& tracks, GeantTaskData& aTaskData, int index ) {
  //std::cout << "Smearer::ApplyTrackerParameterisation : Start" << std::endl;  // Debug
  if ( tracks.fEindexV[index] != 1000 ) {  // Consider only discrete & continuous interactions.
    //std::cout << "Smearer::ApplyTrackerParameterisation : --- End --- EARLY-RETURN : tracks.fEindexV[index]=" 
    //          << tracks.fEindexV[index] << " != 1000" << std::endl;
    return;
  }

  // The charged particle is transported normally to the end of the Tracker,
  // and then the module of the momentum is smeared (while the direction is unchanged).
  double p = tracks.fPV[index];            // In GeV
  double res = GetResolution( DetectorType::eTRACKER, ParameterisationType::eCMS, p );
  //double eff = GetEfficiency( DetectorType::eTRACKER, ParameterisationType::eCMS, p );
  // Smeared the momentum of the particle according to a Gaussian distribution
  // with mean 1.0 and sigma equal to the resolution.
  // Make sure that the smeared momentum is non negative.       
  double smeared_p = p;
  if ( res > 0.0 ) {
    do {
      smeared_p = p * aTaskData.fRndm->Gaus( 1.0, res );
    } while ( smeared_p < 0.0 );
  }
  tracks.fPV[index] = smeared_p;  // Set the smeared momentum

  //std::cout << " Smearer::ApplyTrackerParameterisation : p=" << p 
  //          << " ; smeared_p=" << smeared_p << std::endl
  //          << "Smearer::ApplyTrackerParameterisation : --- End ---" << std::endl;  // Debug
}


void Smearer::ApplyEcalParameterisation( GeantTrack_v& tracks, GeantTaskData& aTaskData, int index ) {
  //std::cout << "Smearer::ApplyEcalParameterisation : Start" << std::endl;  // Debug
  // Consider only discrete & continuous interactions.
  if ( tracks.fEindexV[index] != 1000 ) {
    //std::cout << "Smearer::ApplyEcalParameterisation : --- End --- EARLY-RETURN : tracks.fEindexV[index]=" 
    //          << tracks.fEindexV[index] << " != 1000" << std::endl;
    return;
  }
  // As soon as the electron or positron or gamma enters the electromagnetic calorimeter,
  // its kinetic energy is smeared and deposited and then the particle is killed.
  #ifdef USE_VECGEOM_NAVIGATOR
    const geant::ParticleOld *const & partPDG = &geant::ParticleOld::GetParticle( tracks.fPDGV[index] );
  #else
    TParticlePDG* partPDG = TDatabasePDG::Instance()->GetParticle( tracks.fPDGV[index] );
  #endif

  double mass = partPDG->Mass();   // In GeV
  double eKin = tracks.fEV[index] - mass;  // In GeV
  double p = tracks.fPV[index];            // In GeV
  double res = GetResolution( DetectorType::eEMCAL, ParameterisationType::eCMS, p );
  //double eff = GetEfficiency( DetectorType::eEMCAL, ParameterisationType::eCMS, p );
  // Smeared the original kinetic energy of the particle according to a Gaussian
  // distribution with mean 1.0 and sigma equal to the resolution.
  // Make sure that the smeared kinetic energy is non negative.
  double smeared_eKin = eKin;
  if ( res > 0.0 ) {
    do {
      smeared_eKin = eKin * aTaskData.fRndm->Gaus( 1.0, res );
    } while ( smeared_eKin < 0.0 );
  }
  tracks.fEdepV[index] = smeared_eKin;  // Set the smeared kinetic energy
  tracks.fStatusV[index] = Geant::kKilled;     // Kill the track

  //std::cout << " Smearer::ApplyEcalParameterisation : p=" << p 
  //          << " ; eKin=" << eKin << " ; smeared_eKin=" << smeared_eKin << std::endl
  //          << "Smearer::ApplyEcalParameterisation : --- End ---" << std::endl;  // Debug
}


void Smearer::ApplyHcalParameterisation( GeantTrack_v& tracks, GeantTaskData& aTaskData, int index ) {
  //std::cout << "Smearer::ApplyHcalParameterisation : Start" << std::endl;  // Debug
  // Consider only discrete interactions.
  if ( tracks.fEindexV[index] != 1000 ) {
    //std::cout << "Smearer::ApplyHcalParameterisation : --- End --- EARLY-RETURN : tracks.fEindexV[index]=" 
    //          << tracks.fEindexV[index] << " != 1000" << std::endl;
    return;
  }
  // As soon as the hadron enters the hadronic calorimeter, its kinetic energy
  // is smeared and deposited and then the particle is killed.
  #ifdef USE_VECGEOM_NAVIGATOR
    const geant::ParticleOld *const & partPDG = &geant::ParticleOld::GetParticle( tracks.fPDGV[index] );
  #else
    TParticlePDG* partPDG = TDatabasePDG::Instance()->GetParticle( tracks.fPDGV[index] );
  #endif

  double mass = partPDG->Mass();   // In GeV
  double eKin = tracks.fEV[index] - mass;  // In GeV
  double p = tracks.fPV[index];            // In GeV
  double res = GetResolution( DetectorType::eHCAL, ParameterisationType::eCMS, p );
  //double eff = GetEfficiency( DetectorType::eHCAL, ParameterisationType::eCMS, p );
  // Smeared the original kinetic energy of the particle according to a Gaussian
  // distribution with mean 1.0 and sigma equal to the resolution.
  // Make sure that the smeared kinetic energy is non negative.
  double smeared_eKin = eKin;
  if ( res > 0.0 ) {
    do {
      smeared_eKin = eKin * aTaskData.fRndm->Gaus( 1.0, res );
    } while ( smeared_eKin < 0.0 );
  }
  tracks.fEdepV[index] = smeared_eKin;  // Set the smeared kinetic energy
  tracks.fStatusV[index] = Geant::kKilled;     // Kill the track
  //std::cout << " Smearer::ApplyHcalParameterisation : p=" << p 
  //          << " ; eKin=" << eKin << " ; smeared_eKin=" << smeared_eKin << std::endl
  //          << "Smearer::ApplyHcalParameterisation : --- End ---" << std::endl;  // Debug
}


void Smearer::ApplyMuonParameterisation( GeantTrack_v& /*tracks*/, GeantTaskData& /*aTaskData*/, int /*index*/ ) {
  // Nothing for the time being.
  //std::cout << "Smearer::ApplyMuonParameterisation : Start & End" << std::endl;  // Debug
}


double Smearer::GetResolution( DetectorType aDetector, ParameterisationType aParam, 
                               double aMomentum ) {
  //std::cout << "Smearer::GetResolution : Start" << std::endl;  // Debug
  // Note: the momentum's unit is in GeV
  double res = 1.0;
  if ( aParam == eCMS ) {
    switch ( aDetector ) {
      case Smearer::eTRACKER :
        res = 0.013;
        break;
      case Smearer::eEMCAL :
        res = std::sqrt(   std::pow( 0.03 / std::sqrt( aMomentum ), 2 )  // stochastic
                         + std::pow( 0.12 / aMomentum, 2 )               // noise
                         + std::pow( 0.003, 2 ) );                       // constant
        break;
      case Smearer::eHCAL :
        res = std::sqrt(   std::pow( 1.1 / std::sqrt( aMomentum ), 2 )   // stochastic
                         + std::pow( 0.09, 2 ) );                        // constant
        break;
    }
  } else if ( aParam == eATLAS ) {
    switch ( aDetector ) {
      case Smearer::eTRACKER :
        res = 0.01;
        break;
      case Smearer::eEMCAL :
        res = std::sqrt(   std::pow( 0.1 / std::sqrt( aMomentum ), 2 )   // stochastic
                         + std::pow( 0.0017, 2 ) );                      // constant
        break;
      case Smearer::eHCAL :
        res = std::sqrt(   std::pow( 0.55 / std::sqrt( aMomentum ), 2 )  // stochastic
                         + std::pow( 0.06, 2 ) );                        // constant
        break;
    }
  } else if ( aParam == eALEPH ) {
    switch ( aDetector ) {
      case Smearer::eTRACKER :
        res = 0.01;
        break;
      case Smearer::eEMCAL :
        res = std::sqrt(   std::pow( 0.18 / std::sqrt( aMomentum ), 2 )  // stochastic
                         + std::pow( 0.009, 2 ) );                       // constant
        break;
      case Smearer::eHCAL :
        res = 0.85 / std::sqrt( aMomentum );                             // stochastic
        break;
    }
  }

  //std::cout << " Smearer::GetResolution : p=" << aMomentum 
  //          << " ; aDetector=" << aDetector << " ; aParam=" << aParam << std::endl
  //          << "Smearer::GetResolution : --- End ---" << std::endl;  // Debug
  return res;
}


double Smearer::GetEfficiency( DetectorType aDetector, ParameterisationType /*aParam*/,
                               double /*aMomentum*/ ) {
  //std::cout << "Smearer::GetEfficiency : Start" << std::endl;  // Debug
  // For the time being, we set the efficiency to 1.0
  double eff = 1.0;
  switch ( aDetector ) {
    case Smearer::eTRACKER :
      eff = 1.0;
      break;
    case Smearer::eEMCAL :
      eff = 1.0;
      break;
    case Smearer::eHCAL :
      eff = 1.0;
      break;
  }
  //std::cout << "Smearer::GetResolution : --- End ---" << std::endl;  // Debug
  return eff;
}

