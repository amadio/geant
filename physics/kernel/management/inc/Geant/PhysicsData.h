
#ifndef PHYSICSDATA_H
#define PHYSICSDATA_H

#include <vector>
#include "LightTrack.h"

namespace geantphysics {

class LightTrack;

struct KleinNishinaData {
  static constexpr int dataSize = 2048;
  KleinNishinaData()
      : fEps((double *)vecCore::AlignedAlloc(kPhysDVAlign, sizeof(double) * dataSize)),
        fR0((double *)vecCore::AlignedAlloc(kPhysDVAlign, sizeof(double) * dataSize)),
        fR1((double *)vecCore::AlignedAlloc(kPhysDVAlign, sizeof(double) * dataSize)),
        fR2((double *)vecCore::AlignedAlloc(kPhysDVAlign, sizeof(double) * dataSize)),
        fR3((double *)vecCore::AlignedAlloc(kPhysDVAlign, sizeof(double) * dataSize)),
        fOneMCos((double *)vecCore::AlignedAlloc(kPhysDVAlign, sizeof(double) * dataSize)),
        fSin2t((double *)vecCore::AlignedAlloc(kPhysDVAlign, sizeof(double) * dataSize))
  {
  }
  ~KleinNishinaData()
  {
    vecCore::AlignedFree(fEps);
    vecCore::AlignedFree(fR0);
    vecCore::AlignedFree(fR1);
    vecCore::AlignedFree(fR2);
    vecCore::AlignedFree(fR3);
    vecCore::AlignedFree(fOneMCos);
    vecCore::AlignedFree(fSin2t);
  }
  double *fEps;
  double *fR0;
  double *fR1;
  double *fR2;
  double *fR3;

  double *fOneMCos;
  double *fSin2t;
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
