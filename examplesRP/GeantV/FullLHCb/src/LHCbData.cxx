
#include "LHCbData.h"

#include <iostream>

namespace lhcbapp {

//--------------------------------------------------------------------------------------------------------------------//
// LHCbDataPerPrimary
LHCbDataPerPrimary::LHCbDataPerPrimary() { Clear(); }

void LHCbDataPerPrimary::Clear() {
  fNumChargedSteps = fNumNeutralSteps = fChargedTrackL  = fNeutralTrackL = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons   = 0.;
  fEdep            = 0.;
}

LHCbDataPerPrimary& LHCbDataPerPrimary::operator+=(const LHCbDataPerPrimary& other) {
  fNumChargedSteps += other.fNumChargedSteps;
  fNumNeutralSteps += other.fNumNeutralSteps;
  fChargedTrackL   += other.fChargedTrackL;
  fNeutralTrackL   += other.fNeutralTrackL;
  fNumGammas       += other.fNumGammas;
  fNumElectrons    += other.fNumElectrons;
  fNumPositrons    += other.fNumPositrons;
  fEdep            += other.fEdep;
  return *this;
}


void LHCbDataPerPrimary::Print() {
  std::cout << "    .............................................................  \n"
            << "    Edep                   = " << fEdep                 << " [GeV] \n"
            << "    Track lenght (charged) = " << fChargedTrackL        << "  [cm] \n"
            << "    Track lenght (neutral) = " << fNeutralTrackL        << "  [cm] \n"
            << "    Steps (charged)        = " << fNumChargedSteps      << "       \n"
            << "    Steps (neutral)        = " << fNumNeutralSteps      << "       \n"
            << "    Secondary Gammas       = " << fNumGammas            << "       \n"
            << "    Secondary Electrons    = " << fNumElectrons         << "       \n"
            << "    Secondary Positrons    = " << fNumPositrons         << "       \n"
            << std::endl;
}


//--------------------------------------------------------------------------------------------------------------------//
// LHCbDataPerPrimaryType
LHCbDataPerPrimaryType::LHCbDataPerPrimaryType() { Clear(); }

void LHCbDataPerPrimaryType::Clear() {
  fNumPrimaries    = 0.;
  fNumChargedSteps = fNumNeutralSteps = fNumChargedSteps2 = fNumNeutralSteps2 = 0.;
  fChargedTrackL   = fNeutralTrackL   = fChargedTrackL2   = fNeutralTrackL2   = 0.;
  fNumGammas       = fNumGammas2      = fNumElectrons     = fNumElectrons2    = 0.;
  fNumPositrons    = fNumPositrons2   = fEdep             = fEdep2            = 0.;
}

void LHCbDataPerPrimaryType::AddDataPerPrimary(LHCbDataPerPrimary& data) {
  AddOnePrimary();
  AddChargedSteps(data.GetChargedSteps());
  AddNeutralSteps(data.GetNeutralSteps());
  AddChargedTrackL(data.GetChargedTrackL());
  AddNeutralTrackL(data.GetNeutralTrackL());
  AddGammas   (data.GetGammas()   );
  AddElectrons(data.GetElectrons());
  AddPositrons(data.GetPositrons());
  AddEdep(data.GetEdep());
}


//--------------------------------------------------------------------------------------------------------------------//
// LHCbData
LHCbData::LHCbData(int nprimtypes) : fNumPrimaryTypes(nprimtypes) {
  fDataPerPrimaryType.resize(fNumPrimaryTypes);
  Clear();
}

LHCbData::~LHCbData() {
  fDataPerPrimaryType.clear();
}

void LHCbData::Clear() {
  for (int ipt=0; ipt<fNumPrimaryTypes; ++ipt) {
    fDataPerPrimaryType[ipt].Clear();
  }
}

void LHCbData::AddDataPerPrimaryType(LHCbDataPerPrimary& data, int ptypeindx) {
  fDataPerPrimaryType[ptypeindx].AddDataPerPrimary(data);
}



//--------------------------------------------------------------------------------------------------------------------//
// LHCbDataPerEvent
LHCbDataPerEvent::LHCbDataPerEvent(int nprimperevent) : fNumPrimaryPerEvent(nprimperevent) {
  fPerPrimaryData.resize(fNumPrimaryPerEvent);
  Clear();
}

void LHCbDataPerEvent::Clear() {
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i].Clear();
  }
}

LHCbDataPerEvent& LHCbDataPerEvent::operator+=(const LHCbDataPerEvent &other) {
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i] += other.fPerPrimaryData[i];
  }
  return *this;
}


//--------------------------------------------------------------------------------------------------------------------//
// LHCbThreadDataEvents
  LHCbThreadDataEvents::LHCbThreadDataEvents(int nevtbuffered, int nprimperevent) : fHitsTree(0), fHitsBlock(0),
										    fNumBufferedEvents(nevtbuffered) {
    fPerEventData.reserve(fNumBufferedEvents);
    for (int i=0; i<fNumBufferedEvents; ++i) {
      fPerEventData.push_back(LHCbDataPerEvent(nprimperevent));
    }
  }
  
bool  LHCbThreadDataEvents::Merge(int evtslotindx, const LHCbThreadDataEvents& other) {
  fPerEventData[evtslotindx] += other.GetDataPerEvent(evtslotindx);
  return true;
}



} // namespace lhcbapp
