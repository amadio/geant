

#ifdef USE_ROOT
  #include "TH1F.h"
#else
  #include "Hist.h"
#endif

#include "CaloAppData.h"


namespace userapplication {

//
// CaloAppDataPerPrimary
CaloAppDataPerPrimary::CaloAppDataPerPrimary(int numabs) : fNumAbsorbers(numabs) {
  fEdepInAbsorber.resize(fNumAbsorbers,0.);
  fChargedTrackL.resize(fNumAbsorbers ,0.);
  fNeutralTrackL.resize(fNumAbsorbers ,0.);
  Clear();
}

CaloAppDataPerPrimary::~CaloAppDataPerPrimary() {
  fEdepInAbsorber.clear();
  fChargedTrackL.clear();
  fNeutralTrackL.clear();
}

void CaloAppDataPerPrimary::Clear() {
  for (int k=0;k<fNumAbsorbers; k++){
  	fChargedTrackL[k]  = fNeutralTrackL[k] = fEdepInAbsorber[k] = 0.;
  }
  fNumChargedSteps = fNumNeutralSteps = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons   = 0.;
}

CaloAppDataPerPrimary& CaloAppDataPerPrimary::operator+=(const CaloAppDataPerPrimary& other) {
  for(int k=0;k<fNumAbsorbers;k++) {
  	fChargedTrackL[k]   += other.fChargedTrackL[k];
  	fNeutralTrackL[k]   += other.fNeutralTrackL[k];
  	fEdepInAbsorber[k]  += other.fEdepInAbsorber[k];
  }
  fNumChargedSteps += other.fNumChargedSteps;
  fNumNeutralSteps += other.fNumNeutralSteps;
  fNumGammas       += other.fNumGammas;
  fNumElectrons    += other.fNumElectrons;
  fNumPositrons    += other.fNumPositrons;
  return *this;
}




//
// CaloAppData
CaloAppData::CaloAppData(int numabs) : fNumAbsorbers(numabs) {
  fEdepInAbsorber.resize(fNumAbsorbers,0.);
  fEdepInAbsorber2.resize(fNumAbsorbers,0.);
  fChargedTrackL.resize(fNumAbsorbers,0.);
  fChargedTrackL2.resize(fNumAbsorbers,0.);
  fNeutralTrackL.resize(fNumAbsorbers,0.);
  fNeutralTrackL2.resize(fNumAbsorbers,0.);
  Clear();
}

CaloAppData::~CaloAppData(){
  fEdepInAbsorber.clear();
  fEdepInAbsorber2.clear();
  fChargedTrackL.clear();
  fChargedTrackL2.clear();
  fNeutralTrackL.clear();
  fNeutralTrackL2.clear();
  Clear();
}

void CaloAppData::Clear() {
  for (int k=0; k<fNumAbsorbers;k++) {
  	fChargedTrackL[k]   = fNeutralTrackL[k]     = fChargedTrackL2[k]   = fNeutralTrackL2[k]   = 0.;
  	fEdepInAbsorber[k]  = fEdepInAbsorber2[k]   = 0.;
  }
  fNumChargedSteps = fNumNeutralSteps = fNumChargedSteps2 = fNumNeutralSteps2 = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons     = 0.;
}

void CaloAppData::AddDataPerPrimary(CaloAppDataPerPrimary& data) {
  AddChargedSteps(data.GetChargedSteps());
  AddNeutralSteps(data.GetNeutralSteps());
  for (int k=0; k<fNumAbsorbers; k++) {
  	AddChargedTrackL(data.GetChargedTrackL(k),k);
  	AddNeutralTrackL(data.GetNeutralTrackL(k),k);
  	AddEdepInAbsorber(data.GetEdepInAbsorber(k),k);
  }
  AddGammas   (data.GetGammas()   );
  AddElectrons(data.GetElectrons());
  AddPositrons(data.GetPositrons());
}

//
// CaloAppDataPerEvent
CaloAppDataPerEvent::CaloAppDataPerEvent(int nprimperevent, int numabs) : fNumPrimaryPerEvent(nprimperevent) {
  fPerPrimaryData.reserve(fNumPrimaryPerEvent);
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData.push_back(CaloAppDataPerPrimary(numabs));
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
CaloAppThreadDataEvents::CaloAppThreadDataEvents(int nevtbuffered, int nprimperevent, int numabs) : fNumBufferedEvents(nevtbuffered) {
  fPerEventData.reserve(fNumBufferedEvents);
  for (int i=0; i<fNumBufferedEvents; ++i) {
    fPerEventData.push_back(CaloAppDataPerEvent(nprimperevent,numabs));
  }
}

bool CaloAppThreadDataEvents::Merge(int evtslotindx, const CaloAppThreadDataEvents& other) {
  fPerEventData[evtslotindx] += other.GetDataPerEvent(evtslotindx);
  return true;
}



/*
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

bool CaloAppThreadDataRun::Merge(int evtslotindx, const CaloAppThreadDataRun& other) {
(void)evtslotindx;
#ifdef USE_ROOT
  fHisto1->Add(other.GetHisto1(),1);
#else
  (*fHisto1) += *(other.GetHisto1());
#endif
  return true;
}
*/



} // namespace userapplication
