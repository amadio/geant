
#include "Geant/GSMSCTableSimplified.h"

#include "Geant/SystemOfUnits.h"
#include "Geant/PhysicalConstants.h"

// from material
#include "Geant/Types.h"

#include "Geant/Material.h"
#include "Geant/MaterialProperties.h"
#include "Geant/Element.h"
#include "Geant/MaterialCuts.h"

#include "Geant/GSMottCorrection.h"

// from geantV
#include "Geant/TaskData.h"
#include "Geant/math_wrappers.h"

#include <fstream>
//#include <cstdlib>
#include <cmath>
#include <iostream>
#include <Geant/PhysicsData.h>

namespace geantphysics {

using geant::Double_v;
using geant::IndexD_v;
using geant::kVecLenD;
using geant::MaskD_v;
using vecCore::Get;
using vecCore::Set;
using vecCore::AssignMaskLane;
using vecCore::MaskFull;
using vecCore::MaskEmpty;

// perecomputed GS angular distributions, based on the Screened-Rutherford DCS
// are the same for e- and e+ so make sure we load them only onece
bool GSMSCTableSimplified::gIsInitialised = false;
//
std::vector<GSMSCTableSimplified::GSMSCAngularDtr> GSMSCTableSimplified::gGSMSCAngularDistributions1;
std::vector<GSMSCTableSimplified::GSMSCAngularDtr> GSMSCTableSimplified::gGSMSCAngularDistributions2;
//
std::vector<GSMSCTableSimplified::MoliereData> GSMSCTableSimplified::gMoliere;

GSMSCTableSimplified::GSMSCTableSimplified(bool iselectron) : fIsElectron(iselectron)
{
  // set initial values: final values will be set in the Initialize method
  fLogLambda0        = 0.;
  fLogDeltaLambda    = 0.;
  fInvLogDeltaLambda = 0.;
  fInvDeltaQ1        = 0.;
  fDeltaQ2           = 0.;
  fInvDeltaQ2        = 0.;
}

GSMSCTableSimplified::~GSMSCTableSimplified()
{
  gIsInitialised = false;
}

void GSMSCTableSimplified::Initialize(double /*lownergylimit*/, double /* highenergylimit*/,
                                      const std::vector<bool> &/*activeregionv*/)
{
  double lLambdaMin  = Math::Log(gLAMBMIN);
  double lLambdaMax  = Math::Log(gLAMBMAX);
  fLogLambda0        = lLambdaMin;
  fLogDeltaLambda    = (lLambdaMax - lLambdaMin) / (gLAMBNUM - 1.);
  fInvLogDeltaLambda = 1. / fLogDeltaLambda;
  fInvDeltaQ1        = 1. / ((gQMAX1 - gQMIN1) / (gQNUM1 - 1.));
  fDeltaQ2           = (gQMAX2 - gQMIN2) / (gQNUM2 - 1.);
  fInvDeltaQ2        = 1. / fDeltaQ2;
  // load precomputed angular distributions and set up several values used during the sampling
  // these are particle independet => they go to static container: load them only onece
  if (!gIsInitialised) {
    // load pre-computed GS angular distributions (computed based on Screened-Rutherford DCS)
    LoadMSCData();
    gIsInitialised = true;
  }
  InitMoliereMSCParams();
  // Mott-correction: particle(e- or e+) dependet so init them
}

// samplig multiple scattering angles cos(theta) and sin(thata)
//  - including no-scattering, single, "few" scattering cases as well
//  - Mott-correction will be included if it was requested by the user (i.e. if fIsMottCorrection=true)
// lambdaval : s/lambda_el
// qval      : s/lambda_el G1
// scra      : screening parameter
// cost      : will be the smapled cos(theta)
// sint      : will be the smapled sin(theta)
// lekin     : logarithm of the current kinetic energy
// beta2     : the corresponding beta square
// matindx   : index of the current material
// returns true if it was msc
bool GSMSCTableSimplified::Sampling(double lambdaval, double qval, double scra, double &cost, double &sint,
                                    GSMSCAngularDtr **gsDtr, double &transfPar, geant::TaskData *td, bool isfirst)
{
  double rand0 = td->fRndm->uniform();
  double expn  = Math::Exp(-lambdaval);
  //
  // no scattering case
  if (rand0 < expn) {
    cost = 1.0;
    sint = 0.0;
    return false;
  }
  //
  // single scattering case : sample from the single scattering PDF
  // - Mott-correction will be included if it was requested by the user (i.e. if fIsMottCorrection=true)
  if (rand0 < (1. + lambdaval) * expn) {
    // cost is sampled in SingleScattering()
    double rand1 = td->fRndm->uniform();
    // sample cost from the Screened-Rutherford DCS
    cost = 1. - 2.0 * scra * rand1 / (1.0 - rand1 + scra);
    // add protections
    cost = std::max(cost, -1.0);
    cost = std::min(cost, 1.0);
    //    if (cost<-1.0) cost = -1.0;
    //    if (cost>1.0)  cost =  1.0;
    // compute sin(theta) from the sampled cos(theta)
    double dum0 = 1. - cost;
    sint        = Math::Sqrt(dum0 * (2.0 - dum0));
    return false;
  }
  //
  // handle this case:
  //      -lambdaval < 1 i.e. mean #elastic events along the step is < 1 but
  //       the currently sampled case is not 0 or 1 scattering. [Our minimal
  //       lambdaval (that we have precomputed, transformed angular distributions
  //       stored in a form of equally probabe intervalls together with rational
  //       interp. parameters) is 1.]
  //      -probability of having n elastic events follows Poisson stat. with
  //       lambdaval parameter.
  //      -the max. probability (when lambdaval=1) of having more than one
  //       elastic events is 0.2642411 and the prob of having 2,3,..,n elastic
  //       events decays rapidly with n. So set a max n to 10.
  //      -sampling of this cases is done in a one-by-one single elastic event way
  //       where the current #elastic event is sampled from the Poisson distr.
  if (lambdaval < 1.0) {
    double prob, cumprob;
    prob = cumprob = expn;
    double curcost, cursint;
    // init cos(theta) and sin(theta) to the zero scattering values
    cost = 1.0;
    sint = 0.0;
    for (int iel = 1; iel < 10; ++iel) {
      // prob of having iel scattering from Poisson
      prob *= lambdaval / (double)iel;
      cumprob += prob;
      //
      // sample cos(theta) from the singe scattering pdf:
      // - Mott-correction will be included if it was requested by the user (i.e. if fIsMottCorrection=true)
      double rand1 = td->fRndm->uniform();
      // sample cost from the Screened-Rutherford DCS
      curcost     = 1. - 2.0 * scra * rand1 / (1.0 - rand1 + scra);
      double dum0 = 1. - curcost;
      cursint     = dum0 * (2.0 - dum0); // sin^2(theta)
      //
      // if we got current deflection that is not too small
      // then update cos(theta) sin(theta)
      if (cursint > 1.0e-20) {
        cursint       = Math::Sqrt(cursint);
        double curphi = geant::units::kTwoPi * td->fRndm->uniform();
        cost          = cost * curcost - sint * cursint * std::cos(curphi);
        sint          = Math::Sqrt(std::max(0.0, (1.0 - cost) * (1.0 + cost)));
      }
      //
      // check if we have done enough scattering i.e. sampling from the Poisson
      if (rand0 < cumprob) {
        return false;
      }
    }
    // if reached the max iter i.e. 10
    return false;
  }
  //
  // multiple scattering case with lambdavalue >= 1:
  //   - use the precomputed and transformed Goudsmit-Saunderson angular
  //     distributions to sample cos(theta)
  //   - Mott-correction will be included if it was requested by the user (i.e. if fIsMottCorrection=true)
  cost = SampleCosTheta(lambdaval, qval, scra, gsDtr, transfPar, td, isfirst);
  // add protections
  cost = std::max(cost, -1.0);
  cost = std::min(cost, 1.0);
  //  if (cost<-1.0)  cost = -1.0;
  //  if (cost> 1.0)  cost =  1.0;
  // compute cos(theta) and sin(theta) from the sampled 1-cos(theta)
  double dum0 = 1.0 - cost;
  sint        = Math::Sqrt(dum0 * (2.0 - dum0));
  // return true if it was msc
  return true;
}

double GSMSCTableSimplified::SampleMSCWithSingle(double expn, double lambdaval, double scra, double rndTheta,
                                                 geant::TaskData *td)
{
  double prob, cumprob;
  prob = cumprob = expn;
  double curcost, cursint;
  // init cos(theta) and sin(theta) to the zero scattering values
  double cost = 1.0;
  double sint = 0.0;
  for (int iel = 1; iel < 10; ++iel) {
    // prob of having iel scattering from Poisson
    prob *= 0.5 * lambdaval / (double)iel;
    cumprob += prob;
    //
    // sample cos(theta) from the singe scattering pdf:
    // - Mott-correction will be included if it was requested by the user (i.e. if fIsMottCorrection=true)
    double rand1 = td->fRndm->uniform();
    // sample cost from the Screened-Rutherford DCS
    curcost     = 1. - 2.0 * scra * rand1 / (1.0 - rand1 + scra);
    double dum0 = 1. - curcost;
    cursint     = dum0 * (2.0 - dum0); // sin^2(theta)
    //
    // if we got current deflection that is not too small
    // then update cos(theta) sin(theta)
    if (cursint > 1.0e-20) {
      cursint       = Math::Sqrt(cursint);
      double curphi = geant::units::kTwoPi * td->fRndm->uniform();
      cost          = cost * curcost - sint * cursint * std::cos(curphi);
      sint          = Math::Sqrt(std::max(0.0, (1.0 - cost) * (1.0 + cost)));
    }
    //
    // check if we have done enough scattering i.e. sampling from the Poisson
    if (rndTheta < cumprob) {
      break;
    }
  }
  return cost;
}

void GSMSCTableSimplified::SampleTheta12(const double *lambdaval, const double *qval, const double *scra, double *cost1,
                                         double *cost2, int N, geant::TaskData *td)
{

  auto &rand0Th1 = td->fPhysicsData->fPhysicsScratchpad.rand0Th1;
  rand0Th1.resize(N);
  auto &expn = td->fPhysicsData->fPhysicsScratchpad.expn;
  expn.resize(N);
  auto &loglamda = td->fPhysicsData->fPhysicsScratchpad.loglabmda;
  loglamda.resize(N);
  auto &rand0Th2 = td->fPhysicsData->fPhysicsScratchpad.rand0Th2;
  rand0Th2.resize(N);
  auto &masked = td->fPhysicsData->fPhysicsScratchpad.masked;
  masked.resize(N, false);

  auto &angDtrCache = td->fPhysicsData->fPhysicsScratchpad.angDtrCache;
  angDtrCache.resize(N, nullptr);
  auto &transfParCache = td->fPhysicsData->fPhysicsScratchpad.transfParCache;
  transfParCache.resize(N);

  // Working stack
  auto &angDtr = td->fPhysicsData->fPhysicsScratchpad.angDtr;
  angDtr.clear();
  auto &transfPar = td->fPhysicsData->fPhysicsScratchpad.transfPar;
  transfPar.clear();
  auto &firstAngle = td->fPhysicsData->fPhysicsScratchpad.firstAngle;
  firstAngle.clear();
  auto &idx = td->fPhysicsData->fPhysicsScratchpad.idx;
  idx.clear();

  auto &tempCosStorage = td->fPhysicsData->fPhysicsScratchpad.tempCosStorage;
  tempCosStorage.resize(N * 2);

  for (int i = 0; i < N; ++i) {
    if (0.5 * qval[i] > 7.0) {
      cost1[i]  = 1. - 2. * td->fRndm->uniform();
      cost2[i]  = 1. - 2. * td->fRndm->uniform();
      masked[i] = true;
    } else {
      masked[i] = false;
    }
  }

  for (int i = 0; i < N; i += kVecLenD) {
    Double_v lambdavalVec;
    vecCore::Load(lambdavalVec, lambdaval + i);
    lambdavalVec *= 0.5;
    Double_v rand0Vec = td->fRndm->uniformV();
    Double_v expnVec  = Math::Exp(-lambdavalVec);
    vecCore::Store(rand0Vec, rand0Th1.data() + i);
    vecCore::Store(expnVec, expn.data() + i);
    rand0Vec = td->fRndm->uniformV();
    vecCore::Store(rand0Vec, rand0Th2.data() + i);
    Double_v logLambdaVec = Math::Log(lambdavalVec);
    vecCore::Store(logLambdaVec, loglamda.data() + i);
  }

  for (int i = 0; i < N; ++i) {
    if (masked[i]) continue;
    bool noScat      = rand0Th1[i] < expn[i];
    bool oneScat     = rand0Th1[i] < (1. + 0.5 * lambdaval[i]) * expn[i];
    bool onePlusScat = 0.5 * lambdaval[i] < 1.0;

    if (!(noScat || oneScat || onePlusScat)) { // MSC

      double transfPar_tmp;
      auto gsDtr        = GetGSAngularDtr(scra[i], 0.5 * lambdaval[i], loglamda[i], 0.5 * qval[i], transfPar_tmp, td);
      angDtrCache[i]    = gsDtr;
      transfParCache[i] = transfPar_tmp;
      if (gsDtr) {
        angDtr.push_back(gsDtr);
        transfPar.push_back(transfPar_tmp);
        firstAngle.push_back(true);
        idx.push_back(i);
      } else {
        cost1[i] = 1.0 - 2.0 * td->fRndm->uniform();
      }

    } else {
      if (noScat) {
        cost1[i] = 1.0;
      } else if (oneScat) {
        cost1[i] = SampleWithSingle(scra[i], td);
      } else if (onePlusScat) {
        cost1[i] = SampleMSCWithSingle(expn[i], lambdaval[i], scra[i], rand0Th1[i], td);
      }
    }
  }

  for (int i = 0; i < N; ++i) {
    if (masked[i]) continue;
    bool noScat      = rand0Th2[i] < expn[i];
    bool oneScat     = rand0Th2[i] < (1. + 0.5 * lambdaval[i]) * expn[i];
    bool onePlusScat = 0.5 * lambdaval[i] < 1.0;

    if (!(noScat || oneScat || onePlusScat)) { // MSC
      auto gsDtr           = angDtrCache[i];
      double transfPar_tmp = transfParCache[i];
      if (!gsDtr) {
        gsDtr = GetGSAngularDtr(scra[i], 0.5 * lambdaval[i], loglamda[i], 0.5 * qval[i], transfPar_tmp, td);
      }
      if (gsDtr) {
        angDtr.push_back(gsDtr);
        transfPar.push_back(transfPar_tmp);
        firstAngle.push_back(false);
        idx.push_back(i);
      } else {
        cost2[i] = 1.0 - 2.0 * td->fRndm->uniform();
      }

    } else {
      if (noScat) {
        cost2[i] = 1.0;
      } else if (oneScat) {
        cost2[i] = SampleWithSingle(scra[i], td);
      } else if (onePlusScat) {
        cost2[i] = SampleMSCWithSingle(expn[i], lambdaval[i], scra[i], rand0Th2[i], td);
      }
    }
  }

  size_t tail     = angDtr.size() - (angDtr.size() / kVecLenD) * kVecLenD;
  size_t defficit = tail == 0 ? 0 : kVecLenD - tail;
  if (defficit != 0 && angDtr.size() > 0) {
    angDtr.insert(angDtr.end(), defficit, *(angDtr.end() - 1));
    transfPar.insert(transfPar.end(), defficit, *(transfPar.end() - 1));
    firstAngle.insert(firstAngle.end(), defficit, true);
    idx.insert(idx.end(), defficit, *(idx.end() - 1));
  }

  SampleGSSRCosThetaVector(angDtr.data(), transfPar.data(), tempCosStorage.data(), idx.size(), td);

  for (size_t i = 0; i < idx.size() - defficit; ++i) {
    if (firstAngle[i]) {
      cost1[idx[i]] = tempCosStorage[i];
    } else {
      cost2[idx[i]] = tempCosStorage[i];
    }
  }
}

double GSMSCTableSimplified::SampleCosTheta(double lambdaval, double qval, double scra, GSMSCAngularDtr **gsDtr,
                                            double &transfPar, geant::TaskData *td, bool isfirst)
{
  double cost = 1.;
  // determine the base GS angular distribution if it is the first call (when sub-step sampling is used)
  if (isfirst) {
    *gsDtr = GetGSAngularDtr(scra, lambdaval, qval, transfPar, td);
  }
  // sample cost from the GS angular distribution (computed based on Screened-Rutherford DCS)
  cost = SampleGSSRCosTheta(*gsDtr, transfPar, td);

  return cost;
}

// returns with cost sampled from the GS angular distribution computed based on Screened-Rutherford DCS
double GSMSCTableSimplified::SampleGSSRCosTheta(const GSMSCAngularDtr *gsDtr, double transfpar, geant::TaskData *td)
{
  // check if isotropic theta (i.e. cost is uniform on [-1:1])
  if (!gsDtr) {
    return 1. - 2.0 * td->fRndm->uniform();
  }
  //
  // sampling form the selected distribution
  double ndatm1 = gsDtr->fNumData - 1.;
  double delta  = 1.0 / ndatm1;
  // determine lower cumulative bin inidex
  double rndm = td->fRndm->uniform();
  int indxl   = (int)(rndm * ndatm1);
  double aval = rndm - indxl * delta;
  double dum0 = delta * aval;
  //
  double Ai   = gsDtr->fData[indxl].fParamA;
  double Bi   = gsDtr->fData[indxl].fParamB;
  double Ui   = gsDtr->fData[indxl].fUValues;
  double Uip1 = gsDtr->fData[indxl + 1].fUValues;

  double dum1   = (1.0 + Ai + Bi) * dum0;
  double dum2   = delta * delta + Ai * dum0 + Bi * aval * aval;
  double sample = Ui + dum1 / dum2 * (Uip1 - Ui);
  // transform back u to cos(theta) :
  // this is the sampled cos(theta) = (2.0*para*sample)/(1.0-sample+para)
  return 1. - (2.0 * transfpar * sample) / (1.0 - sample + transfpar);
}

void GSMSCTableSimplified::SampleGSSRCosThetaVector(GSMSCTableSimplified::GSMSCAngularDtr **gsDtr,
                                                    const double *transfpar, double *cost, int N, geant::TaskData *td)
{
  for (int i = 0; i < N; i += kVecLenD) {
    Double_v rndm = td->fRndm->uniformV();
    Double_v Ai;
    Double_v Bi;
    Double_v Ui;
    Double_v Uip1;
    Double_v delta;
    Double_v aval;

    Double_v ndatm1;
    for (int l = 0; l < kVecLenD; ++l) {
      Set(ndatm1, l, gsDtr[i + l]->fNumData - 1.);
      Set(delta, l, gsDtr[i + l]->fDelta);
    }

    IndexD_v indxl = (IndexD_v)(rndm * ndatm1);
    aval           = rndm - indxl * delta;

    for (int l = 0; l < kVecLenD; ++l) {
      // sampling form the selected distribution
      // determine lower cumulative bin inidex
      Set(Ai, l, gsDtr[i + l]->fData[Get(indxl, l)].fParamA);
      Set(Bi, l, gsDtr[i + l]->fData[Get(indxl, l)].fParamB);
      Set(Ui, l, gsDtr[i + l]->fData[Get(indxl, l)].fUValues);
      Set(Uip1, l, gsDtr[i + l]->fData[Get(indxl, l) + 1].fUValues);
    }

    Double_v dum0   = delta * aval;
    Double_v dum1   = (1.0 + Ai + Bi) * dum0;
    Double_v dum2   = delta * delta + Ai * dum0 + Bi * aval * aval;
    Double_v sample = Ui + dum1 / dum2 * (Uip1 - Ui);
    // transform back u to cos(theta) :
    // this is the sampled cos(theta) = (2.0*para*sample)/(1.0-sample+para)
    Double_v transfparVec;
    vecCore::Load(transfparVec, transfpar + i);

    vecCore::Store(1. - (2.0 * transfparVec * sample) / (1.0 - sample + transfparVec), cost + i);
  }
}

GSMSCTableSimplified::GSMSCAngularDtr *GSMSCTableSimplified::GetGSAngularDtr(double scra, double lambdaval,
                                                                             double lLambda, double qval,
                                                                             double &transfpar, geant::TaskData *td)
{
  GSMSCAngularDtr *dtr = nullptr;
  bool first           = false;
  // isotropic cost above gQMAX2 (i.e. dtr stays nullptr)
  if (qval > gQMAX2) return dtr;

  int lamIndx = -1; // lambda value index
  int qIndx   = -1; // lambda value index
  // init to second grid Q values
  int numQVal    = gQNUM2;
  double minQVal = gQMIN2;
  double invDelQ = fInvDeltaQ2;
  double pIndxH  = 0.; // probability of taking higher index
  // check if first or second grid needs to be used
  if (qval < gQMIN2) { // first grid
    first = true;
    // protect against qval<gQMIN1
    if (qval < gQMIN1) {
      qval  = gQMIN1;
      qIndx = 0;
      // pIndxH = 0.;
    }
    // set to first grid Q values
    numQVal = gQNUM1;
    minQVal = gQMIN1;
    invDelQ = fInvDeltaQ1;
  }
  // make sure that lambda = s/lambda_el is in [gLAMBMIN,gLAMBMAX)
  // lambda<gLAMBMIN=1 is already handeled before so lambda>= gLAMBMIN for sure
  if (lambdaval >= gLAMBMAX) {
    lambdaval = gLAMBMAX - 1.e-8;
    lamIndx   = gLAMBNUM - 1;
    lLambda   = Math::Log(lambdaval);
  }
  //
  // determine lower lambda (=s/lambda_el) index: linear interp. on log(lambda) scale
  if (lamIndx < 0) {
    pIndxH  = (lLambda - fLogLambda0) * fInvLogDeltaLambda;
    lamIndx = (int)(pIndxH);    // lower index of the lambda bin
    pIndxH  = pIndxH - lamIndx; // probability of taking the higher index distribution
    if (td->fRndm->uniform() < pIndxH) {
      ++lamIndx;
    }
  }
  //
  // determine lower Q (=s/lambda_el G1) index: linear interp. on Q
  if (qIndx < 0) {
    pIndxH = (qval - minQVal) * invDelQ;
    qIndx  = (int)(pIndxH); // lower index of the Q bin
    pIndxH = pIndxH - qIndx;
    if (td->fRndm->uniform() < pIndxH) {
      ++qIndx;
    }
  }
  // set indx
  int indx = lamIndx * numQVal + qIndx;
  if (first) {
    dtr = &gGSMSCAngularDistributions1[indx];
  } else {
    dtr = &gGSMSCAngularDistributions2[indx];
  }
  dtr = dtr->fNumData == 0 ? nullptr : dtr; // Temporary workaround around different data formats

  // dtr might be nullptr that indicates isotropic cot distribution because:
  // - if the selected lamIndx, qIndx correspond to L(=s/lambda_el) and Q(=s/lambda_el G1) such that G1(=Q/L) > 1
  //   G1 should always be < 1 and if G1 is ~1 -> the dtr is isotropic (this can only happen in case of the 2. grid)
  //
  // compute the transformation parameter
  if (lambdaval > 10.0) {
    transfpar = 0.5 * (-2.77164 + lLambda * (2.94874 - lLambda * (0.1535754 - lLambda * 0.00552888)));
  } else {
    transfpar = 0.5 * (1.347 + lLambda * (0.209364 - lLambda * (0.45525 - lLambda * (0.50142 - lLambda * 0.081234))));
  }
  transfpar *= (lambdaval + 4.0) * scra;
  return dtr;
}
// determine the GS angular distribution we need to sample from: will set other things as well ...
GSMSCTableSimplified::GSMSCAngularDtr *GSMSCTableSimplified::GetGSAngularDtr(double scra, double lambdaval, double qval,
                                                                             double &transfpar, geant::TaskData *td)
{
  GSMSCAngularDtr *dtr = nullptr;
  bool first           = false;
  // isotropic cost above gQMAX2 (i.e. dtr stays nullptr)
  if (qval > gQMAX2) return dtr;

  int lamIndx = -1; // lambda value index
  int qIndx   = -1; // lambda value index
  // init to second grid Q values
  int numQVal    = gQNUM2;
  double minQVal = gQMIN2;
  double invDelQ = fInvDeltaQ2;
  double pIndxH  = 0.; // probability of taking higher index
  // check if first or second grid needs to be used
  if (qval < gQMIN2) { // first grid
    first = true;
    // protect against qval<gQMIN1
    if (qval < gQMIN1) {
      qval  = gQMIN1;
      qIndx = 0;
      // pIndxH = 0.;
    }
    // set to first grid Q values
    numQVal = gQNUM1;
    minQVal = gQMIN1;
    invDelQ = fInvDeltaQ1;
  }
  // make sure that lambda = s/lambda_el is in [gLAMBMIN,gLAMBMAX)
  // lambda<gLAMBMIN=1 is already handeled before so lambda>= gLAMBMIN for sure
  if (lambdaval >= gLAMBMAX) {
    lambdaval = gLAMBMAX - 1.e-8;
    lamIndx   = gLAMBNUM - 1;
  }
  double lLambda = Math::Log(lambdaval);
  //
  // determine lower lambda (=s/lambda_el) index: linear interp. on log(lambda) scale
  if (lamIndx < 0) {
    pIndxH  = (lLambda - fLogLambda0) * fInvLogDeltaLambda;
    lamIndx = (int)(pIndxH);    // lower index of the lambda bin
    pIndxH  = pIndxH - lamIndx; // probability of taking the higher index distribution
    if (td->fRndm->uniform() < pIndxH) {
      ++lamIndx;
    }
  }
  //
  // determine lower Q (=s/lambda_el G1) index: linear interp. on Q
  if (qIndx < 0) {
    pIndxH = (qval - minQVal) * invDelQ;
    qIndx  = (int)(pIndxH); // lower index of the Q bin
    pIndxH = pIndxH - qIndx;
    if (td->fRndm->uniform() < pIndxH) {
      ++qIndx;
    }
  }
  // set indx
  int indx = lamIndx * numQVal + qIndx;
  if (first) {
    dtr = &gGSMSCAngularDistributions1[indx];
  } else {
    dtr = &gGSMSCAngularDistributions2[indx];
  }
  dtr = dtr->fNumData == 0 ? nullptr : dtr; // Temporary workaround around different data formats

  // dtr might be nullptr that indicates isotropic cot distribution because:
  // - if the selected lamIndx, qIndx correspond to L(=s/lambda_el) and Q(=s/lambda_el G1) such that G1(=Q/L) > 1
  //   G1 should always be < 1 and if G1 is ~1 -> the dtr is isotropic (this can only happen in case of the 2. grid)
  //
  // compute the transformation parameter
  if (lambdaval > 10.0) {
    transfpar = 0.5 * (-2.77164 + lLambda * (2.94874 - lLambda * (0.1535754 - lLambda * 0.00552888)));
  } else {
    transfpar = 0.5 * (1.347 + lLambda * (0.209364 - lLambda * (0.45525 - lLambda * (0.50142 - lLambda * 0.081234))));
  }
  transfpar *= (lambdaval + 4.0) * scra;
  return dtr;
}

void GSMSCTableSimplified::LoadMSCData()
{
  // get the path to the main physics data directory
  char *path = std::getenv("GEANT_PHYSICS_DATA");
  if (!path) {
    std::cerr << "******   ERROR in GSMSCTableSimplified::LoadMSCData() \n"
              << "         GEANT_PHYSICS_DATA is not defined! Set the GEANT_PHYSICS_DATA\n"
              << "         environmental variable to the location of Geant data directory!\n"
              << std::endl;
    exit(1);
  }
  //
  gGSMSCAngularDistributions1.resize(gLAMBNUM * gQNUM1);
  for (int il = 0; il < gLAMBNUM; ++il) {
    char fname[512];
    sprintf(fname, "%s/msc_GS/GSGrid_1/gsDistr_%d", path, il);
    std::ifstream infile(fname, std::ios::in);
    if (!infile.is_open()) {
      std::string strfname(fname);
      std::cerr << "******   ERROR in GSMSCTableSimplified::LoadMSCData() \n"
                << "         Cannot open file: " << fname << " \n"
                << std::endl;
      exit(1);
    }
    for (int iq = 0; iq < gQNUM1; ++iq) {
      std::unique_ptr<GSMSCAngularDtr> gsd(new GSMSCAngularDtr());
      infile >> gsd->fNumData;
      gsd->fData.resize(gsd->fNumData);
      double ddummy;
      infile >> ddummy;
      infile >> ddummy;
      for (int i = 0; i < gsd->fNumData; ++i) {
        infile >> gsd->fData[i].fUValues;
        infile >> gsd->fData[i].fParamA;
        infile >> gsd->fData[i].fParamB;
      }
      if (gsd->fNumData > 0) gsd->fDelta            = 1.0 / (gsd->fNumData - 1);
      gGSMSCAngularDistributions1[il * gQNUM1 + iq] = *gsd;
    }
    infile.close();
  }
  //
  // second grid
  gGSMSCAngularDistributions2.resize(gLAMBNUM * gQNUM2);
  for (int il = 0; il < gLAMBNUM; ++il) {
    char fname[512];
    sprintf(fname, "%s/msc_GS/GSGrid_2/gsDistr_%d", path, il);
    std::ifstream infile(fname, std::ios::in);
    if (!infile.is_open()) {
      std::string strfname(fname);
      std::cerr << "******   ERROR in GSMSCTableSimplified::LoadMSCData() \n"
                << "         Cannot open file: " << fname << " \n"
                << std::endl;
      exit(1);
    }
    for (int iq = 0; iq < gQNUM2; ++iq) {
      int numData;
      infile >> numData;
      if (numData > 1) {
        std::unique_ptr<GSMSCAngularDtr> gsd(new GSMSCAngularDtr());
        gsd->fNumData = numData;
        gsd->fData.resize(gsd->fNumData);
        double ddummy;
        infile >> ddummy;
        infile >> ddummy;
        for (int i = 0; i < gsd->fNumData; ++i) {
          infile >> gsd->fData[i].fUValues;
          infile >> gsd->fData[i].fParamA;
          infile >> gsd->fData[i].fParamB;
        }
        if (gsd->fNumData > 0) gsd->fDelta            = 1.0 / (gsd->fNumData - 1);
        gGSMSCAngularDistributions2[il * gQNUM2 + iq] = *gsd;
      }
    }
    infile.close();
  }
}

// compute material dependent Moliere MSC parameters at initialisation
void GSMSCTableSimplified::InitMoliereMSCParams()
{
  constexpr double const1   = 7821.6;         // [cm2/g]
  constexpr double const2   = 0.1569;         // [cm2 MeV2 / g]
  constexpr double finstrc2 = 5.325135453E-5; // fine-structure const. square

  const Vector_t<Material *> &theMaterialTable = Material::GetTheMaterialTable();
  // get number of materials in the table
  size_t numMaterials = theMaterialTable.size();
  // make sure that we have long enough vectors
  gMoliere.resize(numMaterials);
  double xi = 1.0;
  int maxZ  = 200;
  // xi   = 1.0;  <= always set to 1 from now on
  maxZ = GSMottCorrection::GetMaxZet();
  //
  for (size_t imat = 0; imat < numMaterials; ++imat) {
    const Material *theMaterial            = theMaterialTable[imat];
    const Vector_t<Element *> &theElemVect = theMaterial->GetElementVector();
    const double *theNbAtomsPerVolVect     = theMaterial->GetMaterialProperties()->GetNumOfAtomsPerVolumeVect();
    const double theTotNbAtomsPerVol       = theMaterial->GetMaterialProperties()->GetTotalNumOfAtomsPerVol();
    int numelems                           = theElemVect.size();
    //
    double zs = 0.0;
    double zx = 0.0;
    double ze = 0.0;
    double sa = 0.0;
    //
    for (int ielem = 0; ielem < numelems; ielem++) {
      double zet = theElemVect[ielem]->GetZ();
      if (zet > maxZ) {
        zet = (double)maxZ;
      }
      double iwa = theElemVect[ielem]->GetA() * geant::units::mole / geant::units::g;
      double ipz = theNbAtomsPerVolVect[ielem] / theTotNbAtomsPerVol;
      double dum = ipz * zet * (zet + xi);
      zs += dum;
      ze += dum * (-2.0 / 3.0) * Math::Log(zet);
      zx += dum * Math::Log(1.0 + 3.34 * finstrc2 * zet * zet);
      sa += ipz * iwa;
    }
    double density = theMaterial->GetDensity() * geant::units::cm3 / geant::units::g; // [g/cm3]
    //
    gMoliere[theMaterial->GetIndex()].Bc =
        const1 * density * zs / sa * Math::Exp(ze / zs) / Math::Exp(zx / zs); //[1/cm]
    gMoliere[theMaterial->GetIndex()].Xc2 = const2 * density * zs / sa;       // [MeV2/cm]
    // change to internal units of 1/length and energ2/length
    gMoliere[theMaterial->GetIndex()].Bc *= 1.0 / geant::units::cm;
    gMoliere[theMaterial->GetIndex()].Xc2 *= geant::units::MeV * geant::units::MeV / geant::units::cm;
  }
}

double GSMSCTableSimplified::SampleWithSingle(double scra, geant::TaskData *td)
{
  // cost is sampled in SingleScattering()
  double rand1 = td->fRndm->uniform();
  // sample cost from the Screened-Rutherford DCS
  double cost = 1. - 2.0 * scra * rand1 / (1.0 - rand1 + scra);
  // add protections
  cost = std::max(cost, -1.0);
  cost = std::min(cost, 1.0);
  return cost;
}

} // namespace geantphysics
