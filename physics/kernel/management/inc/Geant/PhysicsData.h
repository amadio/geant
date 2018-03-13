
#ifndef PHYSICSDATA_H
#define PHYSICSDATA_H

#include <vector>
#include "LightTrack.h"

namespace geantphysics {

class LightTrack;


struct KleinNishinaData{
  static constexpr int dataSize = 2048;
  KleinNishinaData() : fEps(dataSize), fR0(dataSize), fR1(dataSize), fR2(dataSize), fR3(dataSize){}
  std::vector<double> fEps;
  std::vector<double> fR0;
  std::vector<double> fR1;
  std::vector<double> fR2;
  std::vector<double> fR3;
};

class PhysicsData {
public:
  PhysicsData();
  ~PhysicsData() {}

  /** @brief: Insert secondary and return handle to it **/
  LightTrack &InsertSecondary()
  {
    ResizeIfSmall();
    return fListOfSecondaries[fNumUsedSecondaries++];
  }

  /** @brief: Clear list of used secondaries **/
  void ClearSecondaries() { fNumUsedSecondaries = 0; }

  /** @brief: Array of secondaries, Use GetNumOfSecondaries() to get its size **/
  LightTrack *GetListOfSecondaries() { return fListOfSecondaries.data(); }

  /** @brief: Number of used elements in array returned by GetListOfSecondaries **/
  int GetNumOfSecondaries() const { return fNumUsedSecondaries; }

  std::vector<LightTrack> &GetPrimaryTracks() { return fPrimaryTracks; }

  std::vector<int> &GetSecondaryFillVector() { return fSecondaryFillVector; }

  LightTrack_v &GetPrimarySOA() { return fPrimaryLTs; }
  LightTrack_v &GetSecondarySOA() { return fSecondaryLTs; }

  static void ClearAll();
  static std::vector<PhysicsData *> gThePhysicsDataTable;

  KleinNishinaData fKleinNishinaData;
private:
  void ResizeIfSmall()
  {
    if ((int)fListOfSecondaries.size() == fNumUsedSecondaries) fListOfSecondaries.emplace_back();
  }
  int fNumUsedSecondaries; // number of secondary tracks currently used from fListOfSecondaries
  std::vector<LightTrack> fListOfSecondaries;
  std::vector<LightTrack> fPrimaryTracks;
  std::vector<int> fSecondaryFillVector;
  LightTrack_v fPrimaryLTs;
  LightTrack_v fSecondaryLTs;

};

} // namespace geantphysics

#endif
