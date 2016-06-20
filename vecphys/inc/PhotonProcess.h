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
inline typename Backend::Double_v PhotonProcess::GetLambda(Index_v<typename Backend::Double_v> matId,
                                                           Index_v<typename Backend::Double_v> ebin,
                                                           typename Backend::Double_v fraction) const
{
  int im = (int)(matId);
  int ie = (int)(ebin);

  // linear approximation
  double xlow = fCrossSectionData[im * fNumberOfEnergyBin + ie].fSigma;
  double xhigh = fCrossSectionData[im * fNumberOfEnergyBin + ie + 1].fSigma;

  return xlow + (xhigh - xlow) * fraction;
}

#if !defined(VECCORE_NVCC) && defined(VECCORE_ENABLE_VC)
template <>
inline typename backend::VcVector::Double_v PhotonProcess::GetLambda<backend::VcVector>(
    Index_v<typename backend::VcVector::Double_v> matId, Index_v<typename backend::VcVector::Double_v> ebin,
    typename backend::VcVector::Double_v fraction) const
{
  typedef typename Backend::Double_v Double_t;
  Double_t lambda(0.0);

  for (size_t i = 0; i < VectorSize(ebin); ++i) {
    int im = (int)(matId[i]);
    int ie = (int)(ebin[i]);

    double xlow = fCrossSectionData[im * fNumberOfEnergyBin + ie].fSigma;
    double xhigh = fCrossSectionData[im * fNumberOfEnergyBin + ie + 1].fSigma;
    lambda[i] = xlow + (xhigh - xlow) * fraction[i];
  }
  return lambda;
}
#endif

template <class Backend>
inline VECCORE_CUDA_HOST_DEVICE void PhotonProcess::GetWeightAndAlias(Index_v<typename Backend::Double_v> matId,
                                                                      Index_v<typename Backend::Double_v> ebin,
                                                                      Index_v<typename Backend::Double_v> iprocess,
                                                                      typename Backend::Double_v &weight,
                                                                      Index_v<typename Backend::Double_v> &alias) const
{
  int im = (int)(matId);
  int ie = (int)(ebin);
  int ip = (int)(iprocess);

  if (ip == fNumberOfProcess - 1) {
    weight = 1.0;
    for (int j = 0; j < fNumberOfProcess - 1; ++j)
      weight -= fCrossSectionData[im * fNumberOfEnergyBin + ie].fWeight[j];
  } else {
    weight = fCrossSectionData[im * fNumberOfEnergyBin + ie].fWeight[ip];
  }
  alias = fCrossSectionData[im * fNumberOfEnergyBin + ie].fAlias[ip];
}

#if !defined(VECCORE_NVCC) && defined(VECCORE_ENABLE_VC)
template <>
inline void PhotonProcess::GetWeightAndAlias<backend::VcVector>(
    Index_v<typename backend::VcVector::Double_v> matId, Index_v<typename backend::VcVector::Double_v> ebin,
    Index_v<typename backend::VcVector::Double_v> iprocess, typename backend::VcVector::Double_v &weight,
    Index_v<typename backend::VcVector::Double_v> &alias) const
{
  for (size_t i = 0; i < VectorSize(matId); ++i) {
    int im = (int)(matId[i]);
    int ie = (int)(ebin[i]);
    int ip = (int)(iprocess[i]);

    if (ip == fNumberOfProcess - 1) {
      weight[i] = 1.0;
      for (int j = 0; j < fNumberOfProcess - 1; ++j)
        weight[i] -= fCrossSectionData[im * fNumberOfEnergyBin + ie].fWeight[j];
    } else {
      weight[i] = fCrossSectionData[im * fNumberOfEnergyBin + ie].fWeight[ip];
    }
    alias[i] = fCrossSectionData[im * fNumberOfEnergyBin + ie].fAlias[ip];
  }
}
#endif

template <typename Backend>
inline VECCORE_CUDA_HOST_DEVICE Index_v<typename Backend::Double_v> PhotonProcess::G3NextProcess(
    Index_v<typename Backend::Double_v> matId, Index_v<typename Backend::Double_v> ebin)
{
  // select a physics process randomly based on the weight
  using Double_v = typename Backend::Double_v;

  int im = (int)(matId);
  int ie = (int)(ebin);

  int ip = fNumberOfProcess - 1;

  double weight = 0.0;
  double rp = UniformRandom<Double_v>(fRandomState, fThreadId);

  for (int i = 0; i < fNumberOfProcess - 1; ++i) {
    weight += fCrossSectionData[im * fNumberOfEnergyBin + ie].fWeight[i];
    if (weight > rp) {
      ip = i;
      break;
    }
  }
  return ip;
}

#if !defined(VECCORE_NVCC) && defined(VECCORE_ENABLE_VC)
template <>
inline Index_v<typename backend::VcVector::Double_v> PhotonProcess::G3NextProcess<backend::VcVector>(
    Index_v<typename backend::VcVector::Double_v> matId, Index_v<typename backend::VcVector::Double_v> ebin)
{
  // select a physics process randomly based on the weight
  Index_v<typename backend::VcVector::Double_v> ip = (Index_v<Double_v>)(fNumberOfProcess - 1);

  for (size_t i = 0; i < VectorSize(matId); ++i) {

    int im = (int)(matId[i]);
    int ie = (int)(ebin[i]);

    double weight = 0.0;
    double rp = UniformRandom<double>(fRandomState, fThreadId);

    for (int j = 0; j < fNumberOfProcess - 1; ++j) {
      weight += fCrossSectionData[im * fNumberOfEnergyBin + ie].fWeight[j];
      if (weight > rp) {
        ip[i] = j;
        break;
      }
    }
  }

  return ip;
}
#endif

} // end namespace impl
} // end namespace vecphys

#endif
