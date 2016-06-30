#ifndef GUAliasSampler_H
#define GUAliasSampler_H 1
//
//  Alias Sampler template function implementation
//   For use in implementing methods for scalar, vector, GPU
//   Depends on Backend 'technology' of VecGeom
//
//  First version using 'Backend' - 22 Oct 2014
//   Authors:  Sandro Wenzel, Soon Y. Jun, John Apostolakis
//
//  Desing choice: One Sampler per element
//                    (potentially later per composit material?)
//
//  First version assumes linear 'x' integrand
//   TODO: a measure will be required for the integration for log, theta etc.
//
//  Based on first alias sampler by Soon Y. Jun - July 2014

#include "GUAliasTableManager.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class GUAliasSampler {
public:
  VECCORE_CUDA_HOST
  GUAliasSampler(Random_t *states, int threadId, double incomingMin, double incomingMax,
                 int numEntriesIncoming, // 'energy' (or log) of projectile
                 int numEntriesSampled);

  VECCORE_CUDA_HOST_DEVICE
  GUAliasSampler(Random_t *states, int threadId, double incomingMin, double incomingMax,
                 int numEntriesIncoming, // 'energy' (or log) of projectile
                 int numEntriesSampled, GUAliasTableManager *table);

  VECCORE_CUDA_HOST_DEVICE
  ~GUAliasSampler();

  VECCORE_CUDA_HOST_DEVICE
  void PrintTable();

  VECCORE_CUDA_HOST
  void BuildAliasTable(int z, const double *pdf);

  VECCORE_CUDA_HOST_DEVICE
  GUAliasTableManager *GetAliasTableManager() { return fAliasTableManager; }

  // Backend Implementation:
  template <class Backend>
  VECCORE_CUDA_HOST_DEVICE void SampleBin(typename Backend::Double_v kineticEnergy,
                                          Index_v<typename Backend::Double_v> &index,
                                          Index_v<typename Backend::Double_v> &icol,
                                          typename Backend::Double_v &fraction);

  template <class Backend>
  VECCORE_CUDA_HOST_DEVICE void SampleLogBin(typename Backend::Double_v kineticEnergy,
                                             Index_v<typename Backend::Double_v> &irow,
                                             Index_v<typename Backend::Double_v> &icol,
                                             typename Backend::Double_v &fraction);

  template <class Backend>
  VECCORE_CUDA_HOST_DEVICE typename Backend::Double_v SampleX(typename Backend::Double_v rangeSampled,
                                                              typename Backend::Double_v probNA,
                                                              Index_v<typename Backend::Double_v> aliasInd,
                                                              Index_v<typename Backend::Double_v> icol,
                                                              typename Backend::Double_v fraction);

  template <class Backend>
  VECCORE_CUDA_HOST_DEVICE typename Backend::Double_v SampleXL(Index_v<typename Backend::Double_v> zElement,
                                                               typename Backend::Double_v rangeSampled,
                                                               typename Backend::Double_v probNA,
                                                               Index_v<typename Backend::Double_v> aliasInd,
                                                               Index_v<typename Backend::Double_v> irow,
                                                               Index_v<typename Backend::Double_v> icol);

  template <class Backend>
  inline VECCORE_CUDA_HOST_DEVICE void GatherAlias(Index_v<typename Backend::Double_v> index,
                                                   Index_v<typename Backend::Double_v> zElement,
                                                   typename Backend::Double_v &probNA,
                                                   Index_v<typename Backend::Double_v> &aliasInd) const;

  template <class Backend>
  inline VECCORE_CUDA_HOST_DEVICE typename Backend::Double_v GetPDF(Index_v<typename Backend::Double_v> zElement,
                                                                    Index_v<typename Backend::Double_v> irow,
                                                                    Index_v<typename Backend::Double_v> icol) const;

  // For atomic independent models
  template <class Backend>
  inline VECCORE_CUDA_HOST_DEVICE void GatherAlias(Index_v<typename Backend::Double_v> index,
                                                   typename Backend::Double_v &probNA,
                                                   Index_v<typename Backend::Double_v> &aliasInd) const;

  template <class Backend>
  inline VECCORE_CUDA_HOST_DEVICE typename Backend::Double_v GetPDF(Index_v<typename Backend::Double_v> irow,
                                                                    Index_v<typename Backend::Double_v> icol) const;

  // accessors
  VECCORE_CUDA_HOST_DEVICE
  double GetIncomingMin() const { return fIncomingMin; }

  VECCORE_CUDA_HOST_DEVICE
  double GetIncomingMax() const { return fIncomingMax; }

  VECCORE_CUDA_HOST_DEVICE
  int GetNumEntries() const { return fInNumEntries; }

  VECCORE_CUDA_HOST_DEVICE
  int GetSamplesPerEntry() const { return fSampledNumEntries; }

private:
  Random_t *fRandomState;
  int fThreadId;

  double fIncomingMin; // Min of Incoming - e.g. e_Kinetic or math::Log(E_kinetic)
  double fIncomingMax; // Max
  int fInNumEntries;
  double fLogIncomingMin;
  double fInverseBinIncoming;
  double fInverseLogBinIncoming;

  //  For the sampled variable
  const int fSampledNumEntries; //  Old name fNcol  (number of Columns)
  double fInverseBinSampled;
  GUAliasTableManager *fAliasTableManager;
};

// Backend Implementation

template <class Backend>
VECCORE_CUDA_HOST_DEVICE void GUAliasSampler::SampleBin(
    typename Backend::Double_v kineticEnergy, Index_v<typename Backend::Double_v> & /*index*/, // ~ sampled value
    Index_v<typename Backend::Double_v> &icol,                                                 // ~ input Energy
    typename Backend::Double_v &fraction                                                       //  in sampled variable
    )
{
  using Double_v = typename Backend::Double_v;

  // select the alias table for incoming energy
  Double_v eloc = (kineticEnergy - fIncomingMin) * fInverseBinIncoming;
  Index_v<Double_v> irow = math::Floor(eloc);
  Double_v efrac = eloc - 1.0 * irow;
  // to use fPower2Divisor
  //  Double_v eloc  = (kineticEnergy - fIncomingMin);
  //  Index_v<Double_v> irow = fPower2Divisor->GetBin<Backend>(eloc);
  //  Double_v efrac = fPower2Divisor->FractionWithinBin<Backend>(eloc,irow);

  Double_v u1 = UniformRandom<Double_v>(fRandomState, fThreadId);

  Mask_v<Double_v> useHigh = (u1 <= efrac);

  // irow = useHigh ? irow+1 : irow;
  MaskedAssign(irow, useHigh, irow + 1); // at the upper edge

  Double_v r1 = fSampledNumEntries * UniformRandom<Double_v>(fRandomState, fThreadId);

  // Prepare output values
  icol = math::Floor(r1);
  fraction = r1 - 1.0 * icol;

  // index = irow*fSampledNumEntries  + icol;  // No need to compute - no longer an output
}

template <class Backend>
VECCORE_CUDA_HOST_DEVICE void GUAliasSampler::SampleLogBin(
    typename Backend::Double_v kineticEnergy,
    Index_v<typename Backend::Double_v> &irow, // input energy
    Index_v<typename Backend::Double_v> &icol, // sampled value
    typename Backend::Double_v &fraction       // within the sampled bin
    )
{
  using Double_v = typename Backend::Double_v;

  // select the alias table for incoming energy
  Double_v eloc = (math::Log(kineticEnergy) - fLogIncomingMin) * fInverseLogBinIncoming;

  Double_v irowf = math::Floor(eloc);

  Double_v efrac = eloc - irowf;

  Double_v u1 = UniformRandom<Double_v>(fRandomState, fThreadId);

  Mask_v<Double_v> useHigh = (u1 <= efrac);

  MaskedAssign(irowf, useHigh, irowf + 1.0); // at the upper edge

  // select the sampling bin
  Double_v r1 = fSampledNumEntries * UniformRandom<Double_v>(fRandomState, fThreadId);

  irow = (Index_v<Double_v>)math::Floor(irowf);
  icol = (Index_v<Double_v>)math::Floor(r1);

  fraction = r1 - math::Floor(r1);
}

//    Sample distribution of secondary's 'X' - typically Energy
//      for given zElement ...
//    Feature of this method:  flat distribution within bin

template <class Backend>
VECCORE_CUDA_HOST_DEVICE typename Backend::Double_v GUAliasSampler::SampleX(
    typename Backend::Double_v rangeSampled, typename Backend::Double_v probNA,
    Index_v<typename Backend::Double_v> aliasInd, Index_v<typename Backend::Double_v> icol,
    typename Backend::Double_v fraction)
{
  using Double_v = typename Backend::Double_v;

  Double_v r1 = UniformRandom<Double_v>(fRandomState, fThreadId);

  Double_v xd, xu;
  Double_v binSampled = rangeSampled * fInverseBinSampled;

  Mask_v<Index_v<Double_v>> mask = r1 <= probNA;
  Index_v<Double_v> icolDist = Blend(mask, icol, aliasInd);

  xd = icolDist * binSampled;
  xu = xd + binSampled;

  Double_v x = (1.0 - fraction) * xd + fraction * xu;

  return x;
}

//
//  Lacks a description of what the method does.
//  Since it is a complex method, it will benefit significantly from it.

//  Draft description:
//    Sample distribution of secondary's 'X' - typically Energy
//      for given zElement ...
//    Feature of this method:  linear interpolation using 'PDF'
template <class Backend>
VECCORE_CUDA_HOST_DEVICE typename Backend::Double_v GUAliasSampler::SampleXL(
    Index_v<typename Backend::Double_v> /*zElement*/, typename Backend::Double_v rangeSampled,
    typename Backend::Double_v probNA, Index_v<typename Backend::Double_v> aliasInd,
    Index_v<typename Backend::Double_v> irow, Index_v<typename Backend::Double_v> icol)
{
  using Double_v = typename Backend::Double_v;

  Double_v r1 = UniformRandom<Double_v>(fRandomState, fThreadId);

  Double_v xd, xu;
  Double_v binSampled = rangeSampled * fInverseBinSampled;

  Mask_v<Index_v<Double_v>> mask = r1 <= probNA;
  Index_v<Double_v> icolDist = Blend(mask, icol, aliasInd);

  xd = icolDist * binSampled;
  xu = xd + binSampled;

  // Flat distribution was
  //  Double_v x = (1 - fraction) * xd + fraction* xu;

  // Using pdf of linear interpolation within the sampling bin based on the pdf
  // linear interpolation within the sampling bin based on the pdf
  Double_v pd(0.);
  Double_v pu(0.);

  pd = GetPDF<Backend>(irow, icolDist);
  pu = GetPDF<Backend>(irow, icolDist + 1);

  //* Obtain 'x' in interval [xd, xu] using pdf from linear interpolation
  //    (x,y) from (xd, pd) and (xu, pu)
  //  - Uses two random numbers in order to avoid square root
  Double_v r2 = UniformRandom<Double_v>(fRandomState, fThreadId);
  Double_v r3 = UniformRandom<Double_v>(fRandomState, fThreadId);

  Mask_v<Double_v> below = r2 * (pd + pu) < (1. - r3) * pd + r3 * pu;
  ;

  return Blend(below, (1. - r3) * xd + r3 * xu, r3 * xd + (1. - r3) * xu);
}

template <class Backend>
inline VECCORE_CUDA_HOST_DEVICE
void GUAliasSampler::GatherAlias(Index_v<typename Backend::Double_v> index,
                                 Index_v<typename Backend::Double_v> zElement,
                                 typename Backend::Double_v &probNA,
                                 Index_v<typename Backend::Double_v> &aliasInd) const
{
  for (size_t i = 0; i < VectorSize(index); ++i) {
    int idx  = LaneAt(index, i);
    int Zidx = LaneAt(zElement, i);
    AssignLane(probNA,   i, fAliasTableManager->GetAliasTable(Zidx)->fProbQ[idx]);
    AssignLane(aliasInd, i, fAliasTableManager->GetAliasTable(Zidx)->fAlias[idx]);
  }
}

template <class Backend>
inline VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v GUAliasSampler::GetPDF(Index_v<typename Backend::Double_v> zElement,
                                                  Index_v<typename Backend::Double_v> irow,
                                                  Index_v<typename Backend::Double_v> icol) const
{
  typename Backend::Double_v pdf;
  for (size_t i = 0; i < VectorSize<Double_v>(); ++i) {
    int Zidx = LaneAt(zElement, i);
    int  idx = fSampledNumEntries * LaneAt(irow, i) + LaneAt(icol, i);
    AssignLane(pdf, i, fAliasTableManager->GetAliasTable(Zidx)->fpdf[idx]);
  }
  return pdf;
}

// For atomic independent models

template <class Backend>
inline VECCORE_CUDA_HOST_DEVICE
void GUAliasSampler::GatherAlias(Index_v<typename Backend::Double_v> index, typename Backend::Double_v &probNA,
                                 Index_v<typename Backend::Double_v> &aliasInd) const
{
  for (size_t i = 0; i < VectorSize(index); ++i) {
    int idx = LaneAt(index, i);
    AssignLane(probNA,   i, fAliasTableManager->GetAliasTable(0)->fProbQ[idx]);
    AssignLane(aliasInd, i, fAliasTableManager->GetAliasTable(0)->fAlias[idx]);
  }
}

template <class Backend>
inline VECCORE_CUDA_HOST_DEVICE
typename Backend::Double_v GUAliasSampler::GetPDF(Index_v<typename Backend::Double_v> irow,
                                                  Index_v<typename Backend::Double_v> icol) const
{
  typename Backend::Double_v pdf;
  for (size_t i = 0; i < VectorSize<Double_v>(); ++i) {
    int idx = fSampledNumEntries * LaneAt(irow, i) + LaneAt(icol, i);
    AssignLane(pdf, i, fAliasTableManager->GetAliasTable(0)->fpdf[idx]);
  }
  return pdf;
}

} // end namespace impl
} // end namespace vecphys

#endif
