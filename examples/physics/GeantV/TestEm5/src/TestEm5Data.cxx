
#include "TestEm5Data.h"

#include "Hist.h"


namespace userapplication {

//
// TestEm5DataPerPrimary
TestEm5DataPerPrimary::TestEm5DataPerPrimary() { Clear(); }

void TestEm5DataPerPrimary::Clear() {
  fNumChargedSteps = fNumNeutralSteps = fChargedTrackL  = fNeutralTrackL = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons   = 0.;
  fNumPrimaryTrans = fNumPrimaryRefl  = fNumOneTrans    = fNumOneRefl = 0.;
  fEdepInTarget    = fELeakPrimary    = fELeakSecondary = 0.;
}

TestEm5DataPerPrimary& TestEm5DataPerPrimary::operator+=(const TestEm5DataPerPrimary& other) {
  fNumChargedSteps += other.fNumChargedSteps;
  fNumNeutralSteps += other.fNumNeutralSteps;
  fChargedTrackL   += other.fChargedTrackL;
  fNeutralTrackL   += other.fNeutralTrackL;
  fNumGammas       += other.fNumGammas;
  fNumElectrons    += other.fNumElectrons;
  fNumPositrons    += other.fNumPositrons;
  fNumPrimaryTrans += other.fNumPrimaryTrans;
  fNumPrimaryRefl  += other.fNumPrimaryRefl;
  fNumOneTrans     += other.fNumOneTrans;
  fNumOneRefl      += other.fNumOneRefl;
  fEdepInTarget    += other.fEdepInTarget;
  fELeakPrimary    += other.fELeakPrimary;
  fELeakSecondary  += other.fELeakSecondary;
  return *this;
}




//
// TestEm5Data
TestEm5Data::TestEm5Data() { Clear(); }

void TestEm5Data::Clear() {
  fNumChargedSteps = fNumNeutralSteps = fNumChargedSteps2 = fNumNeutralSteps2 = 0.;
  fChargedTrackL   = fNeutralTrackL   = fChargedTrackL2   = fNeutralTrackL2   = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons     = 0.;
  fNumPrimaryTrans = fNumPrimaryRefl  = fNumOneTrans      = fNumOneRefl       = 0.;
  fEdepInTarget    = fEdepInTarget2   = 0.;
  fELeakPrimary    = fELeakSecondary  = fELeakPrimary2    = fELeakSecondary2  = 0.;
}

void TestEm5Data::AddDataPerPrimary(TestEm5DataPerPrimary& data) {
  AddChargedSteps(data.GetChargedSteps());
  AddNeutralSteps(data.GetNeutralSteps());
  AddChargedTrackL(data.GetChargedTrackL());
  AddNeutralTrackL(data.GetNeutralTrackL());
  AddGammas   (data.GetGammas()   );
  AddElectrons(data.GetElectrons());
  AddPositrons(data.GetPositrons());
  AddPrimaryTransmitted(data.GetPrimaryTransmitted());
  AddPrimaryReflected(data.GetPrimaryReflected());
  AddOneTransmitted(data.GetOneTransmitted());
  AddOneReflected(data.GetOneReflected());
  AddEdepInTarget(data.GetEdepInTarget());
  AddELeakPrimary(data.GetELeakPrimary());
  AddELeakSecondary(data.GetELeakSecondary());
}




//
// TestEm5DataPerEvent
TestEm5DataPerEvent::TestEm5DataPerEvent(int nprimperevent) : fNumPrimaryPerEvent(nprimperevent) {
  fPerPrimaryData.reserve(fNumPrimaryPerEvent);
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData.push_back(TestEm5DataPerPrimary());
  }
}

void TestEm5DataPerEvent::Clear() {
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i].Clear();
  }
}

TestEm5DataPerEvent& TestEm5DataPerEvent::operator+=(const TestEm5DataPerEvent &other) {
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i] += other.fPerPrimaryData[i];
  }
  return *this;
}




//
// TestEm5ThreadDataEvents
TestEm5ThreadDataEvents::TestEm5ThreadDataEvents(int nevtbuffered, int nprimperevent) : fNumBufferedEvents(nevtbuffered) {
  fPerEventData.reserve(fNumBufferedEvents);
  for (int i=0; i<fNumBufferedEvents; ++i) {
    fPerEventData.push_back(TestEm5DataPerEvent(nprimperevent));
  }
}

bool  TestEm5ThreadDataEvents::Merge(int evtslotindx, const TestEm5ThreadDataEvents& other) {
  fPerEventData[evtslotindx] += other.GetDataPerEvent(evtslotindx);
  return true;
}




//
// TestEm5ThreadDataRun
TestEm5ThreadDataRun::TestEm5ThreadDataRun() : fHisto1(nullptr) {}

TestEm5ThreadDataRun::~TestEm5ThreadDataRun() {
  if (fHisto1) {
    delete fHisto1;
  }
  fHisto1 = nullptr;
}

void TestEm5ThreadDataRun::CreateHisto1(int nbins, double min, double max) {
  if (fHisto1) {
    delete fHisto1;
  }
  fHisto1= new Hist(min, max, nbins);
}

bool TestEm5ThreadDataRun::Merge(int /*evtslotindx*/, const TestEm5ThreadDataRun& other) {
  (*fHisto1) += *(other.GetHisto1());
  return true;
}




} // namespace userapplication
