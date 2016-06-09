#ifndef ElectronProcess_H
#define ElectronProcess_H 1

#include "Random.h"
#include "base/PhysicalConstants.h"
#include "base/VecPhys.h"

#include "GUTrack.h"

#include "BremSeltzerBerger.h"
#include "IonisationMoller.h"

#include "EmProcess.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

struct ElectronCrossSectionData {
  double fSigma;
  double fWeight[2];
  int fAlias[3];
};

class ElectronProcess : public EmProcess<ElectronProcess> {
public:
  VECCORE_CUDA_HOST
  ElectronProcess(Random_t *states = 0, int threadId = -1);

  VECCORE_CUDA_HOST_DEVICE
  ~ElectronProcess();

  VECCORE_CUDA_HOST
  void Initialization();

  VECCORE_CUDA_HOST
  void BuildCrossSectionTable();

  VECCORE_CUDA_HOST
  void PrintCrossSectionTable();

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
  friend class EmProcess<ElectronProcess>;

private:
  IonisationMoller *fIonisation;
  BremSeltzerBerger *fBremsstrahlung;

  ElectronCrossSectionData **fElectronCrossSectionData;
};

template <class Backend>
inline typename Backend::Double_v ElectronProcess::GetLambda(Index_v<typename Backend::Double_v> matId,
                                                             Index_v<typename Backend::Double_v> ebin,
                                                             typename Backend::Double_v fraction) const
{
  int im = (int)(matId);
  int ie = (int)(ebin);
  // linear approximation
  double xlow = fElectronCrossSectionData[im][ie].fSigma;
  double xhigh = fElectronCrossSectionData[im][ie + 1].fSigma;

  return xlow + (xhigh - xlow) * fraction;
}

#if !defined(VECCORE_NVCC) && defined(VECCORE_ENABLE_VC)
template <>
inline typename backend::VcVector::Double_v ElectronProcess::GetLambda<backend::VcVector>(
    Index_v<typename backend::VcVector::Double_v> matId, Index_v<typename backend::VcVector::Double_v> ebin,
    typename backend::VcVector::Double_v fraction) const
{
  typedef typename Backend::Double_v Double_t;
  Double_t lambda(0.0);

  for (size_t i = 0; i < VectorSize(ebin); ++i) {
    int im = (int)(matId[i]);
    int ie = (int)(ebin[i]);
    // test to call the scalar method: lambda[i] = GetLambda(im,ie,fraction[i]);
    double xlow = fElectronCrossSectionData[im][ie].fSigma;
    double xhigh = fElectronCrossSectionData[im][ie + 1].fSigma;
    lambda[i] = xlow + (xhigh - xlow) * fraction[i];
  }
  return lambda;
}
#endif

template <class Backend>
inline VECCORE_CUDA_HOST_DEVICE void ElectronProcess::GetWeightAndAlias(
    Index_v<typename Backend::Double_v> matId, Index_v<typename Backend::Double_v> ebin,
    Index_v<typename Backend::Double_v> iprocess, typename Backend::Double_v &weight,
    Index_v<typename Backend::Double_v> &alias) const
{
  int im = (int)(matId);
  int ie = (int)(ebin);
  int ip = (int)iprocess;

  if (ip == fNumberOfProcess - 1) {
    for (int j = 0; j < fNumberOfProcess - 1; ++j)
      weight -= fElectronCrossSectionData[im][ie].fWeight[j];
    weight += 1.0;
  }
  else {
    weight = fElectronCrossSectionData[im][ie].fWeight[ip];
  }
  alias = fElectronCrossSectionData[im][ie].fAlias[ip];
}

#if !defined(VECCORE_NVCC) && defined(VECCORE_ENABLE_VC)
template <>
inline void ElectronProcess::GetWeightAndAlias<backend::VcVector>(
    Index_v<typename backend::VcVector::Double_v> matId, Index_v<typename backend::VcVector::Double_v> ebin,
    Index_v<typename backend::VcVector::Double_v> iprocess, typename backend::VcVector::Double_v &weight,
    Index_v<typename backend::VcVector::Double_v> &alias) const
{
  for (size_t i = 0; i < VectorSize(matId); ++i) {
    int im = (int)(matId[i]);
    int ie = (int)(ebin[i]);
    int ip = (int)(iprocess[i]);

    if (ip == fNumberOfProcess - 1) {
      for (int j = 0; j < fNumberOfProcess - 1; ++j)
        weight[i] -= fElectronCrossSectionData[im][ie].fWeight[j];
      weight[i] += 1.0;
    }
    else {
      weight[i] = fElectronCrossSectionData[im][ie].fWeight[ip];
    }
    alias[i] = fElectronCrossSectionData[im][ie].fAlias[ip];
  }
}
#endif

template <typename Backend>
inline VECCORE_CUDA_HOST_DEVICE Index_v<typename Backend::Double_v> ElectronProcess::G3NextProcess(
    Index_v<typename Backend::Double_v> matId, Index_v<typename Backend::Double_v> ebin)
{
  // select a physics process randomly based on the weight
  int im = (int)(matId);
  int ie = (int)(ebin);

  int ip = fNumberOfProcess - 1;

  double weight = 0.0;
  double rp = UniformRandom<double>(fRandomState, fThreadId);

  for (int i = 0; i < fNumberOfProcess - 1; ++i) {
    weight += fElectronCrossSectionData[im][ie].fWeight[i];
    if (weight > rp) {
      ip = i;
      break;
    }
  }
  return ip;
}

#if !defined(VECCORE_NVCC) && defined(VECCORE_ENABLE_VC)
template <>
inline Index_v<typename backend::VcVector::Double_v> ElectronProcess::G3NextProcess<backend::VcVector>(
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
      weight += fElectronCrossSectionData[im][ie].fWeight[j];
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
