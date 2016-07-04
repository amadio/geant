#ifndef PhotonProcess_H
#define PhotonProcess_H 1

#include "Random.h"
#include "base/PhysicalConstants.h"
#include "base/VecPhys.h"

#include "GUTrack.h"

#include "ComptonKleinNishina.h"
#include "ConversionBetheHeitler.h"
#include "PhotoElectronSauterGavrila.h"

#include "EmProcess.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class PhotonProcess : public EmProcess<PhotonProcess> {
public:
  VECCORE_CUDA_HOST
  PhotonProcess(Random_t *states = 0, int threadId = -1);

  VECCORE_CUDA_HOST_DEVICE
  PhotonProcess(Random_t *states, int threadId, CrossSectionData *data);

  VECCORE_CUDA_HOST_DEVICE
  ~PhotonProcess();

  VECCORE_CUDA_HOST
  void Initialization();

  VECCORE_CUDA_HOST
  void BuildCrossSectionTable();

  template <class Backend>
  inline VECCORE_CUDA_HOST_DEVICE typename Backend::Double_v GetLambda(
      Index_v<typename Backend::Double_v> materialIndex, Index_v<typename Backend::Double_v> ebin,
      typename Backend::Double_v fraction) const;

  template <class Backend>
  inline VECCORE_CUDA_HOST_DEVICE void GetWeightAndAlias(Index_v<typename Backend::Double_v> materialIndex,
                                                         Index_v<typename Backend::Double_v> ebin,
                                                         Index_v<typename Backend::Double_v> iprocess,
                                                         typename Backend::Double_v &weight,
                                                         Index_v<typename Backend::Double_v> &alias) const;

  template <typename Backend>
  inline VECCORE_CUDA_HOST_DEVICE Index_v<typename Backend::Double_v> G3NextProcess(
      Index_v<typename Backend::Double_v> materialIndex, Index_v<typename Backend::Double_v> ebin);

  // the mother is friend in order to access private methods of this
  friend class EmProcess<PhotonProcess>;

private:
  ComptonKleinNishina *fCompton;
  ConversionBetheHeitler *fConversion;
  PhotoElectronSauterGavrila *fPhotoElectron;
};

template <class Backend>
inline VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v PhotonProcess::GetLambda(Index_v<typename Backend::Double_v> matId,
                                                    Index_v<typename Backend::Double_v> ebin,
                                                    typename Backend::Double_v fraction) const
{
  typename Backend::Double_v xlo, xhi;

  // cannot use gather because of usage of struct member fSigma
  for (size_t i = 0; i < VectorSize(ebin); ++i) {
    int ie = LaneAt(ebin, i), im = LaneAt(matId, i);
    AssignLane(xlo, i, fCrossSectionData[im * fNumberOfEnergyBin + ie].fSigma);
    AssignLane(xhi, i, fCrossSectionData[im * fNumberOfEnergyBin + ie + 1].fSigma);
  }

  return xlo + (xhi - xlo) * fraction;
}

template <class Backend>
inline VECCORE_CUDA_HOST_DEVICE void PhotonProcess::GetWeightAndAlias(Index_v<typename Backend::Double_v> matId,
                                                                      Index_v<typename Backend::Double_v> ebin,
                                                                      Index_v<typename Backend::Double_v> iprocess,
                                                                      typename Backend::Double_v &weight,
                                                                      Index_v<typename Backend::Double_v> &alias) const
{
  for (size_t i = 0; i < VectorSize(ebin); ++i) {
    int ie = LaneAt(ebin, i);
    int im = LaneAt(matId, i);
    int ip = LaneAt(iprocess, i);
    double scalar_weight;

    if (ip + 1 == fNumberOfProcess) {
      scalar_weight = 1.0;
      for (size_t j = 0; j < fNumberOfProcess - 1; ++j)
        scalar_weight -= fCrossSectionData[im * fNumberOfEnergyBin + ie].fWeight[j];
    } else {
        scalar_weight = fCrossSectionData[im * fNumberOfEnergyBin + ie].fWeight[ip];
    }

    AssignLane(weight, i, scalar_weight);
    AssignLane(alias, i, fCrossSectionData[im * fNumberOfEnergyBin + ie].fAlias[ip]);
  }
}

template <typename Backend>
inline VECCORE_CUDA_HOST_DEVICE Index_v<typename Backend::Double_v> PhotonProcess::G3NextProcess(
    Index_v<typename Backend::Double_v> matId, Index_v<typename Backend::Double_v> ebin)
{
  // select a physics process randomly based on the weight
  using Double_v = typename Backend::Double_v;

  Double_v rng = UniformRandom<Double_v>(fRandomState, fThreadId);
  Index_v<Double_v> ip(fNumberOfProcess - 1);

  for (size_t i = 0; i < VectorSize(matId); ++i) {
    int ie = LaneAt(ebin, i);
    int im = LaneAt(matId, i);
    double weight = 0.0;

    for (size_t j = 0; j < fNumberOfProcess - 1; ++j) {
      weight += fCrossSectionData[im * fNumberOfEnergyBin + ie].fWeight[j];
      if (weight > LaneAt(rng, i)) {
        AssignLane(ip, i, j);
        break;
      }
    }
  }

  return ip;
}

} // end namespace impl
} // end namespace vecphys

#endif
