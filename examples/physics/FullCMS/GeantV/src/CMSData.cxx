
#include "CMSData.h"

#include "Geant/SystemOfUnits.h"

#include <iostream>
#include <iomanip>

namespace cmsapp {

//--------------------------------------------------------------------------------------------------------------------//
// CMSDataPerPrimary
CMSDataPerPrimary::CMSDataPerPrimary()
{
  Clear();
}

void CMSDataPerPrimary::Clear()
{
  fNumChargedSteps = fNumNeutralSteps = fChargedTrackL = fNeutralTrackL = 0.;
  fChargedTrackT = fNeutralTrackT = 0.;
  fNumGammas = fNumElectrons = fNumPositrons = 0.;
  fEdep                                      = 0.;
  unsigned long numRegions                   = vecgeom::Region::GetTheRegionTable().size();
  fNumStepPerRegion.resize(numRegions);
  for (size_t i = 0; i < numRegions; i++)
    fNumStepPerRegion[i] = 0;
}

CMSDataPerPrimary &CMSDataPerPrimary::operator+=(const CMSDataPerPrimary &other)
{
  fNumChargedSteps += other.fNumChargedSteps;
  fNumNeutralSteps += other.fNumNeutralSteps;
  fChargedTrackL += other.fChargedTrackL;
  fNeutralTrackL += other.fNeutralTrackL;
  fChargedTrackT += other.fChargedTrackT;
  fNeutralTrackT += other.fNeutralTrackT;
  fNumGammas += other.fNumGammas;
  fNumElectrons += other.fNumElectrons;
  fNumPositrons += other.fNumPositrons;
  fEdep += other.fEdep;
  for (size_t i = 0; i < fNumStepPerRegion.size(); ++i)
    fNumStepPerRegion[i] += other.fNumStepPerRegion[i];
  return *this;
}

void CMSDataPerPrimary::Print()
{
  std::cout << "    .............................................................  \n"
            << "    Edep                    = " << std::setw(9) << std::right << fEdep << " [GeV] \n"
            << "    Track lenght (charged)  = " << std::setw(9) << std::right << fChargedTrackL << "  [cm] \n"
            << "    Track lenght (neutral)  = " << std::setw(9) << std::right << fNeutralTrackL << "  [cm] \n"
            << "    Time (global-charged)   = " << std::setw(9) << std::right << fChargedTrackT/geant::units::ns << "  [ns] \n"
            << "    Time (global-neutral)   = " << std::setw(9) << std::right << fNeutralTrackT/geant::units::ns << "  [ns] \n"
            << "    Steps (charged)         = " << std::setw(9) << std::right << fNumChargedSteps << "       \n"
            << "    Steps (neutral)         = " << std::setw(9) << std::right << fNumNeutralSteps << "       \n"
            
            << "    Secondary Gammas        = " << std::setw(9) << std::right  << fNumGammas << "       \n"
            << "    Secondary Electrons     = " << std::setw(9) << std::right  << fNumElectrons << "       \n"
            << "    Secondary Positrons     = " << std::setw(9) << std::right  << fNumPositrons << "       \n";
  for (size_t i = 0; i < fNumStepPerRegion.size(); ++i)
    if( fNumStepPerRegion[i] > 0 ) 
       std::cout << "    Steps per reg. " << std::setw(4) << i << "      = " << std::setw(13) << std::right
                 << fNumStepPerRegion[i] << "       \n";
}

//--------------------------------------------------------------------------------------------------------------------//
// CMSDataPerPrimaryType
CMSDataPerPrimaryType::CMSDataPerPrimaryType()
{
  Clear();
}

void CMSDataPerPrimaryType::Clear()
{
  fNumPrimaries    = 0.;
  fNumChargedSteps = fNumNeutralSteps = fNumChargedSteps2 = fNumNeutralSteps2 = 0.;
  fChargedTrackL = fNeutralTrackL = fChargedTrackL2 = fNeutralTrackL2 = 0.;
  fChargedTrackT = fNeutralTrackT = fChargedTrackT2 = fNeutralTrackT2 = 0.;
  fNumGammas = fNumGammas2 = fNumElectrons = fNumElectrons2 = 0.;
  fNumPositrons = fNumPositrons2 = fEdep = fEdep2 = 0.;
  unsigned long numRegions                        = vecgeom::Region::GetTheRegionTable().size();
  fNumStepPerRegion.resize(numRegions);
  for (size_t i = 0; i < numRegions; i++)
    fNumStepPerRegion[i] = 0;
}

void CMSDataPerPrimaryType::AddDataPerPrimary(CMSDataPerPrimary &data)
{
  AddOnePrimary();
  AddChargedSteps(data.GetChargedSteps());
  AddNeutralSteps(data.GetNeutralSteps());
  AddChargedTrackL(data.GetChargedTrackL());
  AddNeutralTrackL(data.GetNeutralTrackL());
  AddChargedTrackT(data.GetChargedTrackT());
  AddNeutralTrackT(data.GetNeutralTrackT());
  AddGammas(data.GetGammas());
  AddElectrons(data.GetElectrons());
  AddPositrons(data.GetPositrons());
  AddEdep(data.GetEdep());
  for (size_t i = 0; i < fNumStepPerRegion.size(); ++i) {
    AddStepPerRegion(i, data.GetStepsPerRegion(i));
  }
}

//--------------------------------------------------------------------------------------------------------------------//
// CMSData
CMSData::CMSData(int nprimtypes) : fNumPrimaryTypes(nprimtypes)
{
  fDataPerPrimaryType.resize(fNumPrimaryTypes);
  Clear();
}

CMSData::~CMSData()
{
  fDataPerPrimaryType.clear();
}

void CMSData::Clear()
{
  for (int ipt = 0; ipt < fNumPrimaryTypes; ++ipt) {
    fDataPerPrimaryType[ipt].Clear();
  }
}

void CMSData::AddDataPerPrimaryType(CMSDataPerPrimary &data, int ptypeindx)
{
  fDataPerPrimaryType[ptypeindx].AddDataPerPrimary(data);
}

//--------------------------------------------------------------------------------------------------------------------//
// CMSDataPerEvent
CMSDataPerEvent::CMSDataPerEvent(int nprimperevent) : fNumPrimaryPerEvent(nprimperevent)
{
  fPerPrimaryData.resize(fNumPrimaryPerEvent);
  Clear();
}

void CMSDataPerEvent::Clear()
{
  for (int i = 0; i < fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i].Clear();
  }
}

CMSDataPerEvent &CMSDataPerEvent::operator+=(const CMSDataPerEvent &other)
{
  for (int i = 0; i < fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i] += other.fPerPrimaryData[i];
  }
  return *this;
}

//--------------------------------------------------------------------------------------------------------------------//
// CMSThreadDataEvents
CMSThreadDataEvents::CMSThreadDataEvents(int nevtbuffered, int nprimperevent) : fNumBufferedEvents(nevtbuffered)
{
  fPerEventData.reserve(fNumBufferedEvents);
  for (int i = 0; i < fNumBufferedEvents; ++i) {
    fPerEventData.push_back(CMSDataPerEvent(nprimperevent));
  }
}

bool CMSThreadDataEvents::Merge(int evtslotindx, const CMSThreadDataEvents &other)
{
  fPerEventData[evtslotindx] += other.GetDataPerEvent(evtslotindx);
  return true;
}

} // namespace cmsapp
