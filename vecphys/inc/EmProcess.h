#ifndef EmProcess_H
#define EmProcess_H 1

#include "base/PhysicalConstants.h"
#include "base/VecPhys.h"

#include "GUTrack.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

struct CrossSectionData {
  double fSigma;
  double fWeight[2];
  int fAlias[3];
};

template <class Process>
class EmProcess {

public:
  VECCORE_ATT_HOST
  EmProcess(Random_t *states = 0, int threadId = -1);

  VECCORE_ATT_HOST_DEVICE
  EmProcess(Random_t *states, int threadId, CrossSectionData *data);

  VECCORE_ATT_HOST
  ~EmProcess() = default;

  VECCORE_ATT_HOST
  void Initialization();

  VECCORE_ATT_HOST
  void BuildCrossSectionTable();

  VECCORE_ATT_HOST
  void BuildAlias();

  VECCORE_ATT_HOST
  CrossSectionData *GetCrossSectionData() { return fCrossSectionData; }

  VECCORE_ATT_HOST_DEVICE
  void SetCrossSectionData(CrossSectionData *table) { fCrossSectionData = table; }

  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE void GetStepLengthAndProcess(GUTrack &track, const int materialIndex);

  template <typename Backend>
  void GetStepLengthAndProcess(GUTrack_v &tracks, const int *materialIndex);

  template <typename Backend>
  void GVStepLengthAndProcess(GUTrack_v &tracks, const int *materialIndex);

  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE void G3StepLengthAndProcess(GUTrack &track, const int materialIndex);

  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE void GetNextStep(Index_v<typename Backend::Double_v> materialIndex,
                                            Index_v<typename Backend::Double_v> ebin, typename Backend::Double_v efrac,
                                            typename Backend::Double_v &nint, typename Backend::Double_v &lambda,
                                            typename Backend::Double_v &step);

  template <typename Backend>
  VECCORE_ATT_HOST_DEVICE Index_v<typename Backend::Double_v> GetNextProcess(
      Index_v<typename Backend::Double_v> materialIndex, Index_v<typename Backend::Double_v> ebin);

  VECCORE_ATT_HOST_DEVICE int GetNumberOfProcess() { return fNumberOfProcess; }

  VECCORE_ATT_HOST_DEVICE int GetNumberOfEnergyBin() { return fNumberOfEnergyBin; }

  VECCORE_ATT_HOST_DEVICE int GetNumberOfMaterialBin() { return fNumberOfMaterialBin; }

protected:
  Random_t *fRandomState;
  int fThreadId;

  int fNumberOfProcess;
  int fNumberOfEnergyBin;
  int fNumberOfMaterialBin;

  double fEnergyLowerBound;
  double fEnergyUpperBound;

  double fLogEnergyLowerBound;
  double fInverseLogEnergyBin;

  CrossSectionData *fCrossSectionData;
};

// Implementation

template <class Process>
VECCORE_ATT_HOST EmProcess<Process>::EmProcess(Random_t *states, int tid)
    : fRandomState(states), fThreadId(tid), fNumberOfProcess(0), fNumberOfEnergyBin(100), fNumberOfMaterialBin(0),
      fEnergyLowerBound(1. * 10e-3), fEnergyUpperBound(1. * 10e+6)
{
}

template <class Process>
VECCORE_ATT_HOST_DEVICE EmProcess<Process>::EmProcess(Random_t *states, int tid, CrossSectionData *data)
    : fRandomState(states), fThreadId(tid), fNumberOfProcess(0), fNumberOfEnergyBin(100), fNumberOfMaterialBin(0),
      fEnergyLowerBound(1. * 10e-3), fEnergyUpperBound(1. * 10e+6)
{
  fCrossSectionData = data;
}

template <class Process>
VECCORE_ATT_HOST void EmProcess<Process>::BuildCrossSectionTable()
{
  static_cast<Process *>(this)->BuildCrossSectionTable();
  this->BuildAias();
}

template <class Process>
template <typename Backend>
VECCORE_ATT_HOST_DEVICE void EmProcess<Process>::GetStepLengthAndProcess(GUTrack &track, const int materialIndex)
{
  using Double_v = typename Backend::Double_v;
  typedef Index_v<typename Backend::Double_v> Index_t;

  // get interaction length

  // assuming that validity of energy is already checked by the scheduler
  Double_v energy = track.E;
  Double_v nint = track.nint;
  Double_v lambda = track.lambda;
  Double_v step = track.s;

  Double_v eloc = (math::Log(energy) - fLogEnergyLowerBound) * fInverseLogEnergyBin;
  Double_v elow = math::Floor(eloc);
  Index_t ebin = Convert<Double_v, Index_v<Double_v>>(elow);
  Double_v efrac = (eloc - elow) * fInverseLogEnergyBin;

  GetNextStep<Backend>(materialIndex, ebin, efrac, nint, lambda, step);

  track.nint = nint;
  track.lambda = lambda;
  track.s = step;

  // II. select a physics process
  track.proc = GetNextProcess<Backend>(materialIndex, ebin);
}

template <class Process>
template <typename Backend>
void EmProcess<Process>::GetStepLengthAndProcess(GUTrack_v &tracks, const int *materialIndex)
{
  using Double_v = typename Backend::Double_v;
  typedef Index_v<typename Backend::Double_v> Index_t;

  int nTracks = tracks.numTracks;
  int ibase = 0;
  int numChunks = (nTracks / VectorSize<Double_v>());

  for (int i = 0; i < numChunks; ++i) {

    // I. get the step length

    Index_t matId(materialIndex[ibase]);
    Double_v energy, nint, lambda, step;

    Load(energy, &tracks.E[ibase]);
    Load(nint, &tracks.nint[ibase]);
    Load(lambda, &tracks.lambda[ibase]);
    Load(step, &tracks.s[ibase]);

    Double_v eloc = (math::Log(energy) - fLogEnergyLowerBound) * fInverseLogEnergyBin;
    Double_v elow = math::Floor(eloc);
    Index_t ebin = Convert<Double_v, Index_v<Double_v>>(elow);
    Double_v efrac = (eloc - elow) * fInverseLogEnergyBin;

    GetNextStep<Backend>(matId, ebin, efrac, nint, lambda, step);

    Store(nint, &tracks.nint[ibase]);
    Store(lambda, &tracks.lambda[ibase]);
    Store(step, &tracks.s[ibase]);

    // II. select a physics process

    Index_t iproc = GetNextProcess<Backend>(matId, ebin);
    LoadStoreImplementation<Index_t>::template Store<int>(iproc, &tracks.proc[ibase]);

    // next operand
    ibase += VectorSize<Double_v>();
  }

  // leftover - do scalar
  for (int i = numChunks * VectorSize<Double_v>(); i < tracks.numTracks; ++i) {

    double senergy = tracks.E[i];
    double snint = tracks.nint[i];
    double slambda = tracks.lambda[i];
    double sstep = tracks.s[i];

    double seloc = (math::Log(senergy) - fLogEnergyLowerBound) * fInverseLogEnergyBin;
    double selow = math::Floor(seloc);
    int sebin = (int)selow;
    double sefrac = (seloc - selow) * fInverseLogEnergyBin;

    GetNextStep<ScalarBackend>(materialIndex[i], sebin, sefrac, snint, slambda, sstep);

    tracks.nint[i] = snint;
    tracks.lambda[i] = slambda;
    tracks.s[i] = sstep;

    // II. select a physics process
    tracks.proc[i] = (int)(GetNextProcess<ScalarBackend>(materialIndex[i], sebin));
  }
}

template <class Process>
template <typename Backend>
void EmProcess<Process>::GVStepLengthAndProcess(GUTrack_v &tracks, const int *materialIndex)
{
  using Double_v = typename Backend::Double_v;
  typedef Index_v<typename Backend::Double_v> Index_t;

  int nTracks = tracks.numTracks;
  int ibase = 0;
  int numChunks = (nTracks / VectorSize<Double_v>());

  for (int i = 0; i < numChunks; ++i) {

    // I. get the step length

    Index_t matId(materialIndex[ibase]);
    Double_v energy, nint, lambda, step;

    Load(energy, &tracks.E[ibase]);
    Load(nint, &tracks.nint[ibase]);
    Load(lambda, &tracks.lambda[ibase]);
    Load(step, &tracks.s[ibase]);

    Double_v eloc = (math::Log(energy) - fLogEnergyLowerBound) * fInverseLogEnergyBin;
    Double_v elow = math::Floor(eloc);
    Index_t ebin = Convert<Double_v, Index_v<Double_v>>(elow);
    Double_v efrac = (eloc - elow) * fInverseLogEnergyBin;

    GetNextStep<Backend>(matId, ebin, efrac, nint, lambda, step);

    Store(nint, &tracks.nint[ibase]);
    Store(lambda, &tracks.lambda[ibase]);
    Store(step, &tracks.s[ibase]);

    // II. select a physics process

    Index_t iproc = static_cast<Process *>(this)->template G3NextProcess<Backend>(matId, ebin);
    LoadStoreImplementation<Index_t>::template Store<int>(iproc, &tracks.proc[ibase]);

    // next operand
    ibase += VectorSize<Double_v>();
  }

  // leftover - do scalar
  for (int i = numChunks * VectorSize<Double_v>(); i < tracks.numTracks; ++i) {

    double senergy = tracks.E[i];
    double snint = tracks.nint[i];
    double slambda = tracks.lambda[i];
    double sstep = tracks.s[i];

    double seloc = (math::Log(senergy) - fLogEnergyLowerBound) * fInverseLogEnergyBin;
    double selow = math::Floor(seloc);
    int sebin = (int)selow;
    double sefrac = (seloc - selow) * fInverseLogEnergyBin;

    GetNextStep<ScalarBackend>(materialIndex[i], sebin, sefrac, snint, slambda, sstep);

    tracks.nint[i] = snint;
    tracks.lambda[i] = slambda;
    tracks.s[i] = sstep;

    // II. select a physics process
    tracks.proc[i] = (int)(GetNextProcess<ScalarBackend>(materialIndex[i], sebin));
  }
}

template <class Process>
template <typename Backend>
VECCORE_ATT_HOST_DEVICE void EmProcess<Process>::GetNextStep(Index_v<typename Backend::Double_v> matId,
                                                              Index_v<typename Backend::Double_v> ebin,
                                                              typename Backend::Double_v efrac,
                                                              typename Backend::Double_v &nint,
                                                              typename Backend::Double_v &lambda,
                                                              typename Backend::Double_v &step)
{
  using Double_v = typename Backend::Double_v;
  typedef Mask_v<typename Backend::Double_v> Bool_t;

  Bool_t reset = nint < 0.0;
  nint = Blend(reset, -math::Log(UniformRandom<Double_v>(fRandomState, fThreadId)), nint - step / lambda);

  lambda = static_cast<Process *>(this)->template GetLambda<Backend>(matId, ebin, efrac);
  step = lambda * nint;
}

template <class Process>
template <typename Backend>
VECCORE_ATT_HOST_DEVICE Index_v<typename Backend::Double_v> EmProcess<Process>::GetNextProcess(
    Index_v<typename Backend::Double_v> matId, Index_v<typename Backend::Double_v> ebin)
{
  // select a physics process

  using Double_v = typename Backend::Double_v;
  //  typedef Mask_v<typename Backend::Double_v> Bool_t;

  Double_v u1 = fNumberOfProcess * UniformRandom<Double_v>(fRandomState, fThreadId);
  Index_v<typename Backend::Double_v> ip = Convert<Double_v, Index_v<Double_v>>(math::Floor(u1));

  //  Bool_t last = (ip == Convert<Double_v, Index_v<Double_v>>(fNumberOfProcess - 1));

  Double_v weight(0.0);
  Index_v<Double_v> alias;
  static_cast<Process *>(this)->template GetWeightAndAlias<Backend>(matId, ebin, ip, weight, alias);

  // calculate non-alias probablity based on weight (i.e., pdf[ip]) and the normalization factor (1/fNumberOfProcess)
  // non-alias probablity = weight/(normalization factor)
  Double_v probNA = weight * fNumberOfProcess;

  // get the relative weight of the cross section for the ip-th model
  // if (random < weigth*fNumberOfProcess ) use ip else use alias[ip]

  Double_v u2 = fNumberOfProcess * UniformRandom<Double_v>(fRandomState, fThreadId);
  Mask_v<Index_v<Double_v>> mask = u2 <= probNA;
  ip = Blend(mask, ip, alias);

  return ip;
}

template <class Process>
template <typename Backend>
VECCORE_ATT_HOST_DEVICE void EmProcess<Process>::G3StepLengthAndProcess(GUTrack &track, const int materialIndex)
{
  using Double_v = typename Backend::Double_v;
  typedef Index_v<typename Backend::Double_v> Index_t;

  // get the next interaction length

  Double_v energy = track.E;
  Double_v nint = track.nint;
  Double_v lambda = track.lambda;
  Double_v step = track.s;

  Double_v eloc = (math::Log(energy) - fLogEnergyLowerBound) * fInverseLogEnergyBin;
  Double_v elow = math::Floor(eloc);
  Index_t ebin = (Index_v<Double_v>)elow;
  Double_v efrac = (eloc - elow) * fInverseLogEnergyBin;

  GetNextStep<Backend>(materialIndex, ebin, efrac, nint, lambda, step);

  track.nint = nint;
  track.lambda = lambda;
  track.s = step;

  // II. select a physics process following the Geant3-approach
  track.proc = static_cast<Process *>(this)->template G3NextProcess<Backend>(materialIndex, ebin);
}

} // end namespace impl
} // end namespace vecphys

#endif
