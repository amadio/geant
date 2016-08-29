#ifndef UrbanWentzelVI_H
#define UrbanWentzelVI_H 1

#include "base/PhysicalConstants.h"
#include "base/VecPhys.h"

#include "GUAliasSampler.h"
#include "GUTrack.h"

#include "EmModelBase.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class UrbanWentzelVI : public EmModelBase<UrbanWentzelVI> {
public:
  VECCORE_ATT_HOST
  UrbanWentzelVI(Random_t *states = 0, int threadId = -1);

  VECCORE_ATT_HOST_DEVICE
  UrbanWentzelVI(Random_t *states, int threadId, GUAliasSampler *sampler);

  VECCORE_ATT_HOST_DEVICE
  ~UrbanWentzelVI(); //{}

  VECCORE_ATT_HOST
  void Initialization();

  // interfaces for tables
  VECCORE_ATT_HOST
  void BuildCrossSectionTablePerAtom(int Z);

  VECCORE_ATT_HOST
  void BuildPdfTable(int Z, double *p);

public:
  // Auxiliary methods
  VECCORE_ATT_HOST_DEVICE
  GUAliasSampler *GetSampler() { return fAliasSampler; }

  VECCORE_ATT_HOST_DEVICE
  void SetSampler(GUAliasSampler *sampler) { fAliasSampler = sampler; }

private:
  // Implementation methods
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v CrossSectionKernel(typename Backend::Double_v energyIn,
                                                                        Index_v<typename Backend::Double_v> zElement);

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE void InteractKernel(typename Backend::Double_v energyIn,
                                              Index_v<typename Backend::Double_v> zElement,
                                              typename Backend::Double_v &energyOut,
                                              typename Backend::Double_v &sinTheta);

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE void InteractKernelCR(typename Backend::Double_v energyIn,
                                                Index_v<typename Backend::Double_v> zElement,
                                                typename Backend::Double_v &energyOut,
                                                typename Backend::Double_v &sinTheta);

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE void InteractKernelUnpack(typename Backend::Double_v energyIn,
                                                    Index_v<typename Backend::Double_v> zElement,
                                                    typename Backend::Double_v &energyOut,
                                                    typename Backend::Double_v &sinTheta,
                                                    Mask_v<typename Backend::Double_v> &status);

  VECCORE_ATT_HOST_DEVICE
  void SampleByCompositionRejection(int Z, double energyIn, double &energyOut, double &sinTheta);

  VECCORE_ATT_HOST double CrossSectionPerAtom(int Z, double energyIn);
  VECCORE_ATT_HOST double GetG4CrossSection(int Z, double energyIn);

  VECCORE_ATT_HOST_DEVICE
  double CalculateDiffCrossSection(int Zelement, double Ein, double outEphoton) const;

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<UrbanWentzelVI>;
};

template <class Backend>
VECCORE_ATT_HOST_DEVICE typename Backend::Double_v UrbanWentzelVI::CrossSectionKernel(
    typename Backend::Double_v /*energy*/, Index_v<typename Backend::Double_v> /*Z*/)
{
  // vector version of CrossSection
  return 1.0;
}

template <class Backend>
VECCORE_ATT_HOST_DEVICE void UrbanWentzelVI::InteractKernel(typename Backend::Double_v energyIn,
                                                            Index_v<typename Backend::Double_v> /*Z*/,
                                                            typename Backend::Double_v &energyOut,
                                                            typename Backend::Double_v &sinTheta)
{
  energyOut = energyIn;
  sinTheta  = 0.0;
}

template <class Backend>
VECCORE_ATT_HOST_DEVICE void UrbanWentzelVI::InteractKernelCR(typename Backend::Double_v energyIn,
                                                              Index_v<typename Backend::Double_v> /*Z*/,
                                                              typename Backend::Double_v &energyOut,
                                                              typename Backend::Double_v &sinTheta)
{
  energyOut = energyIn;
  sinTheta  = 0.0;
}

template <class Backend>
VECCORE_ATT_HOST_DEVICE void UrbanWentzelVI::InteractKernelUnpack(typename Backend::Double_v energyIn,
                                                                  Index_v<typename Backend::Double_v> /*Z*/,
                                                                  typename Backend::Double_v &energyOut,
                                                                  typename Backend::Double_v &sinTheta,
                                                                  Mask_v<typename Backend::Double_v> & /*status*/)
{
  energyOut = energyIn;
  sinTheta  = 0.0;
}

} // end namespace impl
} // end namespace vecphys

#endif
