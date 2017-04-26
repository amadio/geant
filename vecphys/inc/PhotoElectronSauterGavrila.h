#ifndef PhotoElectronSauterGavrila_H
#define PhotoElectronSauterGavrila_H 1

#include "base/PhysicalConstants.h"
#include "base/VecPhys.h"

#include "GUConstants.h"
#include "GUTrack.h"
#include "StaticSandiaData.h"

#include "EmModelBase.h"

#include <stdlib.h>

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class PhotoElectronSauterGavrila : public EmModelBase<PhotoElectronSauterGavrila> {
public:
  VECCORE_ATT_HOST
  PhotoElectronSauterGavrila(Random_t *states, int threadId = -1);

  VECCORE_ATT_HOST_DEVICE
  PhotoElectronSauterGavrila(Random_t *states, int threadId, GUAliasSampler *sampler);

  VECCORE_ATT_HOST_DEVICE
  ~PhotoElectronSauterGavrila() {}

  VECCORE_ATT_HOST
  void Initialization();

  // interfaces for tables
  VECCORE_ATT_HOST
  void BuildCrossSectionTablePerAtom(int Z);

  VECCORE_ATT_HOST
  void BuildPdfTable(int Z, double *p);

  // Alternative Interact method to test energy dependent subtasks for a
  // specific model. Eeventually this method should replace the Interact
  // method of EmBaseModel

  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE void ModelInteract(GUTrack &projectile, const int targetElement, GUTrack &secondary);

// vector
  template <typename Backend>
  void ModelInteract(GUTrack_v &inProjectile, const int *targetElements, GUTrack_v &outSecondaryV);

  // Implementation methods
  template <class Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v CrossSectionKernel(typename Backend::Double_v energyIn,
                                                                        Index_v<typename Backend::Double_v> zElement);

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE void InteractKernel(typename Backend::Double_v energyIn,
                                              Index_v<typename Backend::Double_v> zElement,
                                              typename Backend::Double_v &energyOut,
                                              typename Backend::Double_v &sinTheta);
private:
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
  double GetPhotoElectronEnergyScalar(double E, size_t Z)
  {
    assert(Z > 0 && Z <= 100);

    int i = 0, nShells = fNumberOfShells[Z];

    while (i < nShells && E >= fBindingEnergies[fIndexOfShells[Z] + i] * eV)
      i++;

    return i <= nShells ? E - fBindingEnergies[fIndexOfShells[Z] + i] * eV : 0.0;
  }

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v GetPhotoElectronEnergy(typename Backend::Double_v E,
                                                                            Index_v<typename Backend::Double_v> Z)
  {
    using Double_v = typename Backend::Double_v;
    using DIndex_v = Index_v<Double_v>;
    using Scalar_t = typename ScalarType<DIndex_v>::Type;

    Double_v Eout;
    Double_t *E_ptr    = (Double_t *)&E;
    Double_t *Eout_ptr = (Double_t *)&Eout;
    Scalar_t *Z_ptr    = (Scalar_t *)&Z;

    for (Size_t i = 0; i < VectorSize<Double_v>(); i++) {
      Eout_ptr[i] = GetPhotoElectronEnergyScalar(E_ptr[i], Z_ptr[i]);
    }

    return Eout;
  }

  template <class Backend>
  inline VECCORE_ATT_HOST_DEVICE typename Backend::Double_v SampleSequential(typename Backend::Double_v A,
                                                                             typename Backend::Double_v Ap2,
                                                                             typename Backend::Double_v B,
                                                                             typename Backend::Double_v grej);

  VECCORE_ATT_HOST_DEVICE
  void SampleByCompositionRejection(int Z, double energyIn, double &energyOut, double &sinTheta);

  VECCORE_ATT_HOST double GetG4CrossSection(int Z, double energyIn);

  VECCORE_ATT_HOST_DEVICE
  double CalculateDiffCrossSectionK(int Zelement, double Ein, double outEphoton) const;

  VECCORE_ATT_HOST_DEVICE
  double CalculateDiffCrossSection(int Zelement, double Ein, double outEphoton) const;

  // the mother is friend in order to access private methods of this
  friend class EmModelBase<PhotoElectronSauterGavrila>;

  // private:
};

// Implementation

template <class Backend>
VECCORE_ATT_HOST_DEVICE typename Backend::Double_v PhotoElectronSauterGavrila::CrossSectionKernel(
    typename Backend::Double_v energy, Index_v<typename Backend::Double_v> Z)
{
  using Double_v = typename Backend::Double_v;

  Double_v sigma = 0.;

  // Sandia parameterization for Z < 100
  //  int Z = zElement;

  int fCumulInterval[101] = {0};
  double fSandiaCof[4]    = {0.0};

  fCumulInterval[0] = 1;

  // scan - move it to constructor or use pre-built table
  for (int iz = 1; iz < 101; ++iz) {
    fCumulInterval[iz] = fCumulInterval[iz - 1] + fNbOfIntervals[iz];
  }

  double Emin = fSandiaTable[fCumulInterval[Z - 1]][0] * keV;

  int interval = fNbOfIntervals[Z] - 1;
  int row      = fCumulInterval[Z - 1] + interval;

  while ((interval > 0) && (energy < fSandiaTable[row][0] * keV)) {
    --interval;
    row = fCumulInterval[Z - 1] + interval;
  }

  if (energy >= Emin) {
    double AoverAvo = Z * amu / fZtoAratio[Z];
    fSandiaCof[0]   = AoverAvo * funitc[1] * fSandiaTable[row][1];
    fSandiaCof[1]   = AoverAvo * funitc[2] * fSandiaTable[row][2];
    fSandiaCof[2]   = AoverAvo * funitc[3] * fSandiaTable[row][3];
    fSandiaCof[3]   = AoverAvo * funitc[4] * fSandiaTable[row][4];
  } else {
    fSandiaCof[0] = fSandiaCof[1] = fSandiaCof[2] = fSandiaCof[3] = 0.;
  }

  Double_v energy2 = energy * energy;
  Double_v energy3 = energy * energy2;
  Double_v energy4 = energy2 * energy2;

  Double_v sgima = fSandiaCof[0] / energy + fSandiaCof[1] / energy2 + fSandiaCof[2] / energy3 + fSandiaCof[3] / energy4;

  return sigma;
}

template <class Backend>
VECCORE_ATT_HOST_DEVICE void PhotoElectronSauterGavrila::InteractKernel(typename Backend::Double_v energyIn,
                                                                        Index_v<typename Backend::Double_v> zElement,
                                                                        typename Backend::Double_v &energyOut,
                                                                        typename Backend::Double_v &sinTheta)
{
  using Double_v = typename Backend::Double_v;

  // energy of photo-electron: Sandia parameterization
  energyOut = GetPhotoElectronEnergy<Backend>(energyIn, zElement);

  // sample angular distribution of photo-electron

  Index_v<Double_v> irow;
  Index_v<Double_v> icol;
  Double_v fraction;

  fAliasSampler->SampleLogBin<Backend>(energyIn, irow, icol, fraction);

  Double_v probNA;
  Index_v<Double_v> aliasInd;

  Index_v<Double_v> ncol(fAliasSampler->GetSamplesPerEntry());
  Index_v<Double_v> index = ncol * irow + icol;
  fAliasSampler->GatherAlias<Backend>(index, probNA, aliasInd);

  Double_v mininum = -1.0;
  Double_v deltaE  = 2.0;

  Double_v cosTheta = mininum + fAliasSampler->SampleX<Backend>(deltaE, probNA, aliasInd, icol, fraction);

  sinTheta = math::Sqrt((1 + cosTheta) * (1 - cosTheta));
}

template <class Backend>
VECCORE_ATT_HOST_DEVICE void PhotoElectronSauterGavrila::InteractKernelCR(typename Backend::Double_v energyIn,
                                                                          Index_v<typename Backend::Double_v> zElement,
                                                                          typename Backend::Double_v &energyOut,
                                                                          typename Backend::Double_v &sinTheta)
{
  using Double_v = typename Backend::Double_v;

  // energy of photo-electron: Sandia parameterization
  energyOut = GetPhotoElectronEnergy<Backend>(energyIn, zElement);

  // sample angular direction according to SauterGavrilaAngularDistribution
  Double_v tau = energyIn / electron_mass_c2;

  /*
  const double taulimit = 50.0;
  Mask_v<Double_v> highE = tau > taulimit;
  cosTheta = 1.0;
  if(EarlyReturnAllowed() && MaskFull(highE)) return;
  */

  Double_v gamma = tau + 1.0;
  Double_v beta  = math::Sqrt(tau * (tau + 2.0)) / gamma;

  Double_v A    = (1 - beta) / beta;
  Double_v Ap2  = A + 2;
  Double_v B    = 0.5 * beta * gamma * (gamma - 1.) * (gamma - 2.);
  Double_v grej = 2 * (1 + A * B) / A;

  Double_v z = SampleSequential<Backend>(A, Ap2, B, grej);

  //  MaskedAssign(!highE, 1.0 - z , &cosTheta);
  sinTheta = math::Sqrt(z * (2 - z)); // cosTheta = 1 -z
}
template <class Backend>
inline VECCORE_ATT_HOST_DEVICE typename Backend::Double_v PhotoElectronSauterGavrila::SampleSequential(
    typename Backend::Double_v A, typename Backend::Double_v Ap2, typename Backend::Double_v B,
    typename Backend::Double_v grej)
{
  using Double_v = typename Backend::Double_v;

  Double_v z;
  Double_v g;
  Mask_v<Double_v> done(false);

  do {
    Double_v q = UniformRandom<Double_v>(fRandomState, fThreadId);
    MaskedAssign(z, !done, 2 * A * (2 * q + Ap2 * math::Sqrt(q)) / (Ap2 * Ap2 - 4 * q));
    MaskedAssign(g, !done, (2 - z) * (1.0 / (A + z) + B));
    done |= g < UniformRandom<Double_v>(fRandomState, fThreadId) * grej;
  } while (!MaskFull(done));

  return z;
}

template <class Backend>
VECCORE_ATT_HOST_DEVICE void PhotoElectronSauterGavrila::InteractKernelUnpack(
    typename Backend::Double_v energyIn, Index_v<typename Backend::Double_v> /*zElement*/,
    typename Backend::Double_v &energyOut, typename Backend::Double_v &sinTheta,
    Mask_v<typename Backend::Double_v> & /*status*/)
{
  // dummy for now
  energyOut = energyIn;
  sinTheta  = 0;
}

// Alternative Interact method

template <typename Backend>
VECCORE_ATT_HOST_DEVICE void PhotoElectronSauterGavrila::ModelInteract(GUTrack &inProjectile, const int targetElement,
                                                                       GUTrack &outSecondary)
{
  using Double_v = typename Backend::Double_v;

  Double_v energyIn = inProjectile.E;

  // check for the validity of energy
  if (energyIn < fLowEnergyLimit || energyIn > fHighEnergyLimit) return;

  Double_v energyOut = 0;
  Double_v sinTheta  = 0;

  // a good upper bound of photon energy to apply the alias method for
  // the SauterGavrila angular distribution (above this energy region,
  // dsigma/dcos(theta) is squeezed toward 1
  const double aliaslimit = 1.0 * MeV;

  // lower bound for the approximation dsigma/dcos(theta) =1 driven by Geant4
  //(note that the (geant4) composition and rejection method is very inefficient
  // for the model above this region)
  const double taulimit = 50. * electron_mass_c2;

  if (energyIn < aliaslimit) {
    InteractKernel<Backend>(energyIn, targetElement, energyOut, sinTheta);
  } else if (energyIn < taulimit) {
    InteractKernelCR<Backend>(energyIn, targetElement, energyOut, sinTheta);
  } else {
    energyOut = GetPhotoElectronEnergy<Backend>(energyIn, targetElement);
    sinTheta  = 0; // cosTheta = 1.0;
  }

  // update final states of the primary and store the secondary
  ConvertXtoFinalState<Backend>(energyIn, energyOut, sinTheta, inProjectile, outSecondary);
}

template <typename Backend>
void PhotoElectronSauterGavrila::ModelInteract(GUTrack_v &inProjectile, const int *targetElements,
                                               GUTrack_v &outSecondary)
{
  using Double_v = typename Backend::Double_v;

  // filtering energy regions for sampling methods - setable if necessary
  const double aliaslimit = 1.0 * MeV;
  const double taulimit   = 50. * electron_mass_c2;

  int nTracks          = inProjectile.numTracks;
  double *start        = inProjectile.E;
  auto indexAliasLimit = std::lower_bound(start, start + nTracks, aliaslimit) - start;
  auto indexTauLimit   = std::lower_bound(start, start + nTracks, taulimit) - start;

  for (int j = 0; j < nTracks; ++j) {
    assert((targetElements[j] > 0) && (targetElements[j] <= maximumZ));
  }

  int ibase     = 0;
  int numChunks = (nTracks / VectorSize<Double_v>());

  for (int i = 0; i < numChunks; ++i) {
    Double_v energyIn = FromPtr<Double_v>(&inProjectile.E[ibase]);
    Double_v sinTheta(0.);
    Double_v energyOut;

    Index_v<Double_v> zElement(targetElements[ibase]);

    if (ibase < indexAliasLimit) {
      InteractKernel<Backend>(energyIn, zElement, energyOut, sinTheta);
    } else if (ibase < indexTauLimit) {
      InteractKernelCR<Backend>(energyIn, zElement, energyOut, sinTheta);
    } else {
      energyOut = GetPhotoElectronEnergy<Backend>(energyIn, zElement);
      sinTheta  = 0.; // cosTheta = 1.0;
    }

    ConvertXtoFinalState<Backend>(energyIn, energyOut, sinTheta, ibase, inProjectile, outSecondary);

    ibase += VectorSize<Double_v>();
  }

  // leftover - do scalar (temporary)
  for (int i = numChunks * VectorSize<Double_v>(); i < nTracks; ++i) {

    double senergyIn = inProjectile.E[i];
    double senergyOut, ssinTheta;
    // use InteractKernel for any leftover to be consistent with EmBaseModel
    InteractKernel<ScalarBackend>(senergyIn, targetElements[i], senergyOut, ssinTheta);
    ConvertXtoFinalState_Scalar<ScalarBackend>(senergyIn, senergyOut, ssinTheta, i, inProjectile, outSecondary);
  }
}

} // end namespace impl
} // end namespace vecphys

#endif
