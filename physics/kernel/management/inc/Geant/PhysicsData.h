
#ifndef PHYSICSDATA_H
#define PHYSICSDATA_H

#include <vector>
#include "Geant/GSMSCTableSimplified.h"
#include "LightTrack.h"

namespace geantphysics {

class LightTrack;

struct PhysicsModelScratchpad {
  static constexpr int dataSize = 2048;
  PhysicsModelScratchpad()
      : fEps((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fR0((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fR1((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fR2((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fR3((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr2((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr3((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr4((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr5((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr6((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr7((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr8((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr9((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr10((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr11((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr12((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr13((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fDoubleArr14((double *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(double) * dataSize)),
        fIzet((int *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(int) * dataSize)),
        fMatIdx((int *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(int) * dataSize)),
        fNshells((int *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(int) * dataSize)),
        fSampledShells((int *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(int) * dataSize))
        fBoolArr((bool *)vecCore::AlignedAlloc(geant::kVecAlignD, sizeof(bool) * dataSize))
  {
  }
  ~PhysicsModelScratchpad()
  {
    vecCore::AlignedFree(fEps);
    vecCore::AlignedFree(fR0);
    vecCore::AlignedFree(fR1);
    vecCore::AlignedFree(fR2);
    vecCore::AlignedFree(fR3);
    vecCore::AlignedFree(fDoubleArr);
    vecCore::AlignedFree(fDoubleArr2);
    vecCore::AlignedFree(fNshells);
    vecCore::AlignedFree(fSampledShells);
    vecCore::AlignedFree(fDoubleArr3);
    vecCore::AlignedFree(fDoubleArr4);
    vecCore::AlignedFree(fDoubleArr5);
    vecCore::AlignedFree(fDoubleArr6);
    vecCore::AlignedFree(fDoubleArr7);
    vecCore::AlignedFree(fDoubleArr8);
    vecCore::AlignedFree(fDoubleArr9);
    vecCore::AlignedFree(fDoubleArr10);
    vecCore::AlignedFree(fDoubleArr11);
    vecCore::AlignedFree(fDoubleArr12);
    vecCore::AlignedFree(fDoubleArr13);
    vecCore::AlignedFree(fDoubleArr14);
    vecCore::AlignedFree(fIzet);
    vecCore::AlignedFree(fMatIdx);
  }
  double *fEps;
  double *fR0;
  double *fR1;
  double *fR2;
  double *fR3;

  double *fDoubleArr;
  double *fDoubleArr2;
  double *fDoubleArr3;
  double *fDoubleArr4;
  double *fDoubleArr5;
  double *fDoubleArr6;
  double *fDoubleArr7;
  double *fDoubleArr8;
  double *fDoubleArr9;
  double *fDoubleArr10;
  double *fDoubleArr11;
  double *fDoubleArr12;
  double *fDoubleArr13;
  double *fDoubleArr14;

  int *fIzet;
  int *fMatIdx;
  int *fNshells;
  int *fSampledShells;
  bool *fBoolArr;

  std::vector<bool> hasNewDir;
  std::vector<double> truePathLengths;
  std::vector<double> rand0Th1;
  ;
  std::vector<double> expn;
  ;
  std::vector<double> loglabmda;
  ;
  std::vector<double> rand0Th2;
  ;
  std::vector<bool> masked;
  ;

  std::vector<GSMSCTableSimplified::GSMSCAngularDtr *> angDtrCache;
  std::vector<double> transfParCache;

  // Working stack
  std::vector<GSMSCTableSimplified::GSMSCAngularDtr *> angDtr;
  std::vector<double> transfPar;
  std::vector<bool> firstAngle;
  std::vector<int> idx;
  std::vector<double> tempCosStorage;
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

  PhysicsModelScratchpad fPhysicsScratchpad;

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
