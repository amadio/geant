#ifndef GEANTV_VECTORUTILS_H
#define GEANTV_VECTORUTILS_H

namespace geantphysics {

struct SecondariesFillInfo {
  int fTrackId;
  int fNumSecondaries;

  SecondariesFillInfo(int trackIndex, int numSecondaries) : fTrackId(trackIndex), fNumSecondaries(numSecondaries) {}
};
}

#endif // GEANTV_VECTORUTILS_H
