#include "TestEm3Data.h"

#include <stddef.h>

namespace userapplication {

//
// TestEm3DataPerPrimary
TestEm3DataPerPrimary::TestEm3DataPerPrimary(int numabs) : fNumAbsorbers(numabs)
{
  fEdepInAbsorber.resize(fNumAbsorbers, 0.);
  fChargedTrackL.resize(fNumAbsorbers, 0.);
  fNeutralTrackL.resize(fNumAbsorbers, 0.);
  Clear();
}

TestEm3DataPerPrimary::~TestEm3DataPerPrimary()
{
  fEdepInAbsorber.clear();
  fChargedTrackL.clear();
  fNeutralTrackL.clear();
}

void TestEm3DataPerPrimary::Clear()
{
  for (int k = 0; k < fNumAbsorbers; k++) {
    fChargedTrackL[k] = fNeutralTrackL[k] = fEdepInAbsorber[k] = 0.;
  }
  fNumChargedSteps = fNumNeutralSteps = 0.;
  fNumGammas = fNumElectrons = fNumPositrons = 0.;
}

TestEm3DataPerPrimary &TestEm3DataPerPrimary::operator+=(const TestEm3DataPerPrimary &other)
{
  for (int k = 0; k < fNumAbsorbers; k++) {
    fChargedTrackL[k] += other.fChargedTrackL[k];
    fNeutralTrackL[k] += other.fNeutralTrackL[k];
    fEdepInAbsorber[k] += other.fEdepInAbsorber[k];
  }
  fNumChargedSteps += other.fNumChargedSteps;
  fNumNeutralSteps += other.fNumNeutralSteps;
  fNumGammas += other.fNumGammas;
  fNumElectrons += other.fNumElectrons;
  fNumPositrons += other.fNumPositrons;
  return *this;
}

//
// TestEm3Data
TestEm3Data::TestEm3Data(int numabs) : fNumAbsorbers(numabs)
{
  fEdepInAbsorber.resize(fNumAbsorbers, 0.);
  fEdepInAbsorber2.resize(fNumAbsorbers, 0.);
  fChargedTrackL.resize(fNumAbsorbers, 0.);
  fChargedTrackL2.resize(fNumAbsorbers, 0.);
  fNeutralTrackL.resize(fNumAbsorbers, 0.);
  fNeutralTrackL2.resize(fNumAbsorbers, 0.);
  Clear();
}

TestEm3Data::~TestEm3Data()
{
  fEdepInAbsorber.clear();
  fEdepInAbsorber2.clear();
  fChargedTrackL.clear();
  fChargedTrackL2.clear();
  fNeutralTrackL.clear();
  fNeutralTrackL2.clear();
  Clear();
}

void TestEm3Data::Clear()
{
  for (int k = 0; k < fNumAbsorbers; k++) {
    fChargedTrackL[k] = fNeutralTrackL[k] = fChargedTrackL2[k] = fNeutralTrackL2[k] = 0.;
    fEdepInAbsorber[k] = fEdepInAbsorber2[k] = 0.;
  }
  fNumChargedSteps = fNumNeutralSteps = fNumChargedSteps2 = fNumNeutralSteps2 = 0.;
  fNumGammas = fNumElectrons = fNumPositrons = 0.;
}

void TestEm3Data::AddDataPerPrimary(TestEm3DataPerPrimary &data)
{
  AddChargedSteps(data.GetChargedSteps());
  AddNeutralSteps(data.GetNeutralSteps());
  for (int k = 0; k < fNumAbsorbers; k++) {
    AddChargedTrackL(data.GetChargedTrackL(k), k);
    AddNeutralTrackL(data.GetNeutralTrackL(k), k);
    AddEdepInAbsorber(data.GetEdepInAbsorber(k), k);
  }
  AddGammas(data.GetGammas());
  AddElectrons(data.GetElectrons());
  AddPositrons(data.GetPositrons());
}

//
// TestEm3DataPerEvent
TestEm3DataPerEvent::TestEm3DataPerEvent(int nprimperevent, int numabs) : fNumPrimaryPerEvent(nprimperevent)
{
  fPerPrimaryData.reserve(fNumPrimaryPerEvent);
  for (int i = 0; i < fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData.push_back(TestEm3DataPerPrimary(numabs));
  }
}

void TestEm3DataPerEvent::Clear()
{
  for (int i = 0; i < fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i].Clear();
  }
}

TestEm3DataPerEvent &TestEm3DataPerEvent::operator+=(const TestEm3DataPerEvent &other)
{
  for (int i = 0; i < fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i] += other.fPerPrimaryData[i];
  }
  return *this;
}

//
// TestEm3DataEvents
TestEm3ThreadDataEvents::TestEm3ThreadDataEvents(int nevtbuffered, int nprimperevent, int numabs)
    : fNumBufferedEvents(nevtbuffered)
{
  fPerEventData.reserve(fNumBufferedEvents);
  for (int i = 0; i < fNumBufferedEvents; ++i) {
    fPerEventData.push_back(TestEm3DataPerEvent(nprimperevent, numabs));
  }
}

bool TestEm3ThreadDataEvents::Merge(int evtslotindx, const TestEm3ThreadDataEvents &other)
{
  fPerEventData[evtslotindx] += other.GetDataPerEvent(evtslotindx);
  return true;
}

//
// TestEm3ThreadDataRun
bool TestEm3ThreadDataRun::Merge(int /*evtslotindx*/, const TestEm3ThreadDataRun &other)
{
  const std::vector<double> &chv = other.GetCHTrackLPerLayer();
  const std::vector<double> &edv = other.GetEDepPerLayer();
  size_t nLayers                 = chv.size();
  for (size_t i = 0; i < nLayers; ++i) {
    fChargedTrackLPerLayer[i] += chv[i];
    fEdepPerLayer[i] += edv[i];
  }
  return true;
}

} // namespace userapplication
