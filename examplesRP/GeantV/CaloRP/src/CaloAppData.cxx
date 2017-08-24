

#ifdef USE_ROOT
  #include "TH1F.h"
#else
  #include "Hist.h"
#endif

#include "CaloAppData.h"


namespace userapplication {

//
// CaloAppDataPerPrimary
CaloAppDataPerPrimary::CaloAppDataPerPrimary() { Clear(); }

void CaloAppDataPerPrimary::Clear() {

  for (int k=1;k<maxAbsorbers; k++){
  	fChargedTrackL[k]  = fNeutralTrackL[k] = fEdepInAbsorber[k] = 0.;
  }

  fNumChargedSteps = fNumNeutralSteps = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons   = 0.;
  fELeakPrimary    = fELeakSecondary = 0.;
}

CaloAppDataPerPrimary& CaloAppDataPerPrimary::operator+=(const CaloAppDataPerPrimary& other) {
  for(int k=1;k<maxAbsorbers;k++) {
  	fChargedTrackL[k]   += other.fChargedTrackL[k];
  	fNeutralTrackL[k]   += other.fNeutralTrackL[k];
  	fEdepInAbsorber[k]    += other.fEdepInAbsorber[k];
  }
  fNumChargedSteps += other.fNumChargedSteps;
  fNumNeutralSteps += other.fNumNeutralSteps;
  fNumGammas       += other.fNumGammas;
  fNumElectrons    += other.fNumElectrons;
  fNumPositrons    += other.fNumPositrons;
  fELeakPrimary    += other.fELeakPrimary;
  fELeakSecondary  += other.fELeakSecondary;
  return *this;
}




//
// CaloAppData
CaloAppData::CaloAppData() { Clear(); }

void CaloAppData::Clear() {
  for (int k=1; k<maxAbsorbers;k++){
  	fChargedTrackL[k]   = fNeutralTrackL[k]   = fChargedTrackL2[k]   = fNeutralTrackL2[k]   = 0.;
  	fEdepInAbsorber[k]    = fEdepInAbsorber2[k]   = 0.;
  }

  fNumChargedSteps = fNumNeutralSteps = fNumChargedSteps2 = fNumNeutralSteps2 = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons     = 0.;
  fELeakPrimary    = fELeakSecondary  = fELeakPrimary2    = fELeakSecondary2  = 0.;
}

void CaloAppData::AddDataPerPrimary(CaloAppDataPerPrimary& data) {
  AddChargedSteps(data.GetChargedSteps());
  AddNeutralSteps(data.GetNeutralSteps());

  for (int k=1; k<maxAbsorbers; k++){
  	AddChargedTrackL(data.GetChargedTrackL(k),k);
  	AddNeutralTrackL(data.GetNeutralTrackL(k),k);
  	AddEdepInAbsorber(data.GetEdepInAbsorber(k),k);
  }
  AddGammas   (data.GetGammas()   );
  AddElectrons(data.GetElectrons());
  AddPositrons(data.GetPositrons());


  AddELeakPrimary(data.GetELeakPrimary());
  AddELeakSecondary(data.GetELeakSecondary());
}

//
// CaloAppDataPerEvent
CaloAppDataPerEvent::CaloAppDataPerEvent(int nprimperevent) : fNumPrimaryPerEvent(nprimperevent) {
  fPerPrimaryData.reserve(fNumPrimaryPerEvent);
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData.push_back(CaloAppDataPerPrimary());
  }
}

void CaloAppDataPerEvent::Clear() {
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i].Clear();
  }
}

CaloAppDataPerEvent& CaloAppDataPerEvent::operator+=(const CaloAppDataPerEvent &other) {
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i] += other.fPerPrimaryData[i];
  }
  return *this;
}




//
// CaloAppDataEvents
CaloAppThreadDataEvents::CaloAppThreadDataEvents(int nevtbuffered, int nprimperevent) : fNumBufferedEvents(nevtbuffered) {
  fPerEventData.reserve(fNumBufferedEvents);
  for (int i=0; i<fNumBufferedEvents; ++i) {
    fPerEventData.push_back(CaloAppDataPerEvent(nprimperevent));
  }
}

bool CaloAppThreadDataEvents::Merge(int evtslotindx, const CaloAppThreadDataEvents& other) {
  fPerEventData[evtslotindx] += other.GetDataPerEvent(evtslotindx);
  return true;
}




//
// CaloAppThreadDataRun
CaloAppThreadDataRun::CaloAppThreadDataRun() : fHisto1(nullptr) {}

CaloAppThreadDataRun::~CaloAppThreadDataRun() {
  if (fHisto1) {
    delete fHisto1;
  }
  fHisto1 = nullptr;
}

void CaloAppThreadDataRun::CreateHisto1(int nbins, double min, double max) {
  if (fHisto1) {
    delete fHisto1;
  }
#ifdef USE_ROOT
  fHisto1= new TH1F("HistName", "Hist Title", nbins, min, max);
#else
  fHisto1= new Hist(min, max, nbins);
#endif
}

bool CaloAppThreadDataRun::Merge(int /*evtslotindx*/, const CaloAppThreadDataRun& other) {
#ifdef USE_ROOT
  fHisto1->Add(other.GetHisto1(),1);
#else
  (*fHisto1) += *(other.GetHisto1());
#endif
  return true;
}




} // namespace userapplication
