#include "TNudyCore.h"
#include "TNudyAliasCont.h"

#ifdef USE_ROOT
ClassImp(TNudyAlias)
#endif

    //_______________________________________________________________________________
    TNudyAliasCont::TNudyAliasCont()
    : fLen(0), fChooseBin(nullptr), fP(nullptr), fX(nullptr), fInterX(nullptr), fInterP(nullptr), fTx(nullptr),
      fTp(nullptr), fTa(-1), fInterAlpha(10), fRnd(nullptr), fAlpha(10)
#ifdef TNUDYALIAS_MULTITHREAD
      ,
      fMult(nullptr), fMultLen(0)
#endif
{
}

//_______________________________________________________________________________
TNudyAliasCont::TNudyAliasCont(double *data, int len, double alpha, unsigned int seed)
    : fLen(0), fChooseBin(nullptr), fP(nullptr), fX(nullptr), fInterX(nullptr), fInterP(nullptr), fTx(nullptr),
      fTp(nullptr), fTa(-1), fInterAlpha(10), fRnd(nullptr), fAlpha(10)
#ifdef TNUDYALIAS_MULTITHREAD
      ,
      fMult(nullptr), fMultLen(0)
#endif
{
  // double* of x1,p1,x2,p2.....xn,pn
  int i = 0, j = 0;
  if (len % 2 != 0) {
    printf("TNudyAliasCont::TNudyAliasCont: Incorrect array specified\n");
  }
  double *x = new double[len / 2];
  double *p = new double[len / 2];
  for (i = 0, j = 0; i < len; i += 2, j++) {
    x[j] = data[i];
    p[j] = data[i + 1];
  }
  Initialize(p, x, len, alpha, seed);
  delete[] p;
  delete[] x;
}

//_______________________________________________________________________________
TNudyAliasCont::TNudyAliasCont(double *p, double *x, int len, double alpha, unsigned int seed)
    : fLen(0), fChooseBin(nullptr), fP(nullptr), fX(nullptr), fInterX(nullptr), fInterP(nullptr), fTx(nullptr),
      fTp(nullptr), fTa(-1), fInterAlpha(10), fRnd(nullptr), fAlpha(10)
#ifdef TNUDYALIAS_MULTITHREAD
      ,
      fMult(nullptr), fMultLen(0)
#endif
{
  Initialize(p, x, len, alpha, seed);
}

//_______________________________________________________________________________
void TNudyAliasCont::Initialize(double *p, double *x, const int len, double alpha, unsigned int seed)
{
  // Data sorted(asc) in x
  // alpha to be used when using statistical interpolation between two distributions
  int i;
  double *integral = nullptr, *I = nullptr;
  fInterX = fInterP = nullptr;
  //  fTx = new double[len];
  //  fTp = new double[len];
  fTa         = -1;
  fInterAlpha = 0;
#ifdef TNUDYALIAS_MULTITHREAD
  fMult    = nullptr;
  fMultLen = 0;
#endif
  fLen   = len;
  fAlpha = alpha;
#ifdef USE_ROOT
  fRnd = new TRandom(seed);
#else
  fRnd = &RNG::Instance();
#endif
  fP          = new double[fLen];
  fX          = new double[fLen];
  integral    = new double[len - 1];
  I           = new double[len - 1];
  integral[0] = 0;
  I[0]        = 0;
  double h    = 0;
  for (i = 0; i < len; i++) {
    fP[i] = p[i];
    fX[i] = x[i];
    if (i > 0) h += (x[i] - x[i - 1]) * (p[i] + p[i - 1]) * 0.5;
  }
  for (i = 0; i < len; i++) {
    fP[i] /= h;
    if (i > 0) {
      integral[i - 1] = (fX[i] - fX[i - 1]) * (fP[i] + fP[i - 1]) * 0.5;
      I[i - 1]        = i - 1;
    }
  }
  fChooseBin = new TNudyAlias(integral, I, len - 1, seed);
  delete[] integral;
  delete[] I;
}

//_______________________________________________________________________________
TNudyAliasCont::~TNudyAliasCont()
{
  delete fChooseBin;
  delete[] fP;
  delete[] fX;
  delete[] fInterX;
  delete[] fInterP;
  delete[] fTx;
  delete[] fTp;
#ifdef USE_ROOT
  delete fRnd;
#endif
#ifdef TNUDYALIAS_MULTITHREAD
  delete[] fMult;
#endif
}

//_______________________________________________________________________________
void TNudyAliasCont::DumpTable()
{
  int i;
  printf("Continuous probability distribution --> \n");
  for (i = 0; i < fLen; i++) {
    printf("x[%d]=%e, p[%d]=%e\n", i, fX[i], i, fP[i]);
  }
  if (fInterX && fInterP)
    for (i = 0; i < fLen; i++)
      printf("fIx[%i]=%e, fIp[%d]=%e\n", i, fInterX[i], i, fInterP[i]);
  printf("Alias table -->\n");
  fChooseBin->DumpTable();
}

//_______________________________________________________________________________
double TNudyAliasCont::Random(IntScheme_t iScheme, int binNo, AliasDist_t distribution, double r1, double r2)
{
  if (binNo < 0 || binNo >= fLen) // Choose bins if required as per their probabilities
    binNo    = (int)fChooseBin->Random();
  double *xx = 0;
  double *xp = 0;
  switch (distribution) {
  case kOriginal:
    xx = fX;
    xp = fP;
    break;
  case kBuilt:
    xx = (fInterX && fInterP) ? fInterX : fX;
    xp = (fInterX && fInterP) ? fInterP : fP;
    break;
  default:
    printf("TNudyAliasCont::Random: Unknown distribution type\n");
    return 0.;
  }
  double rnd1, rnd2, x1, x2;
  switch (iScheme) {
  case kLinear:
    // Coerce linear probabiliy distributions to equal probability bins
    rnd1 = r1 < 0 ? fRnd->Uniform() : r1;
    rnd2 = r2 < 0 ? fRnd->Uniform() : r2;
    x1   = (1 - rnd1) * xx[binNo] + rnd1 * xx[binNo + 1];
    x2   = rnd1 * xx[binNo] + (1 - rnd1) * xx[binNo + 1];
    if (rnd2 * (xp[binNo] + xp[binNo + 1]) <= (1 - rnd1) * xp[binNo] + rnd1 * xp[binNo + 1])
      return x1;
    else
      return x2;
    break;
  case kDiscrete:
    // Generate discrete random numbers using probability distribution
    return xx[binNo];
    break;
  case kUniform:
    // Generate random numbers between bins using a uniform probability distribution
    x1 = (xx[binNo + 1] - xx[binNo]) * fRnd->Uniform(1);
    return x1 + xx[binNo];
    break;
  default:
    return 0;
    break;
  }
}

//_______________________________________________________________________________
double TNudyAliasCont::StatisticalInterpolation(TNudyAliasCont *dist1, TNudyAliasCont *dist2, double alpha)
{
  // Simple statistical interpolation between distributions ( Not very accurate )
  if (dist1->Uniform(1) <= ((dist2->fAlpha - alpha) / (dist2->fAlpha - dist1->fAlpha)))
    return dist1->Random();
  else
    return dist2->Random();
}

//_______________________________________________________________________________
double TNudyAliasCont::ImprovedInterpolation(double alpha)
{
  // Improved interpolation between distributions - Interpolates between corresonding values in the
  // intergral of probability density
  double x1, x2, X1, X2, P1, P2, rnd1, rnd2, u;
  u = (alpha - fAlpha) / (fInterAlpha - fAlpha);
  if (alpha < fAlpha || alpha > fInterAlpha) {
    printf("TNudyAliasCont::ImprovedInterpolation: %e does not lie between %e and %e\n", alpha, fAlpha, fInterAlpha);
  }
  double h      = 0;
  if (!fTx) fTx = new double[fLen];
  if (!fTp) fTp = new double[fLen];
  if (alpha != fTa) {
    fTa = alpha;
    for (int i = 0; i < fLen; i++) {
      fTx[i] = (1 - u) * fX[i] + u * fInterX[i];
      fTp[i] = (1 - u) * fP[i] + u * fInterP[i];
      if (i > 0) {
        h += 0.5 * (fTx[i] - fTx[i - 1]) * (fTp[i] + fTp[i - 1]);
      }
    }
    if (fabs(1.0 - h) > 1e-5) {
      for (int i = 0; i < fLen; i++) {
        fTp[i] /= h;
      }
    }
  }
  int binNo = (int)fChooseBin->Random();
  X1        = fTx[binNo];
  X2        = fTx[binNo + 1];
  P1        = fTp[binNo];
  P2        = fTp[binNo + 1];
  rnd1      = fRnd->Uniform();
  rnd2      = fRnd->Uniform();
  x1        = (1.0 - rnd1) * X1 + rnd1 * X2;
  x2        = rnd1 * X1 + (1.0 - rnd1) * X2;
  if (rnd2 * (P1 + P2) <= (1.0 - rnd1) * P1 + rnd1 * P2)
    return x1;
  else
    return x2;
}
//_______________________________________________________________________________
double TNudyAliasCont::SelectBin(double /*alpha*/, AliasDist_t distribution)
{
  int binNo = (int)fChooseBin->Random();
  if (distribution == kOriginal)
    return fX[binNo];
  else
    return fInterX[binNo];
}
//_______________________________________________________________________________
void TNudyAliasCont::BuildIntermediate(TNudyAliasCont *dists, const int len)
{
  // Populates data structures to interpolate between distributions
  int i, j, k;
  double *integral   = nullptr;
  double *integralp1 = nullptr;
  double x;
  int lo, hi, mid, min;
  int count;
  double f1, f2, e1, e2, dele, delf, temp, pp;
  double h = 0;
  lo = hi = mid = 0;
  for (i = 0; i < len - 1; i++) {
    dists[i].fInterAlpha = dists[i + 1].fAlpha;
    integral             = new double[dists[i].fLen];     // Integral for dists[i]
    integralp1           = new double[dists[i + 1].fLen]; // Integral for dists[i+1]
    integral[0]          = 0;
    integralp1[0]        = 0;
    // Evaluate integrals
    for (j = 1; j < dists[i].fLen; j++) {
      integral[j] = 0.5 * (dists[i].fX[j] - dists[i].fX[j - 1]) * (dists[i].fP[j] + dists[i].fP[j - 1]);
      integral[j] = integral[j] + integral[j - 1];
    }
    // Store bin areas of
    for (j = 1; j < dists[i + 1].fLen; j++) {
      integralp1[j] =
          0.5 * (dists[i + 1].fX[j] - dists[i + 1].fX[j - 1]) * (dists[i + 1].fP[j] + dists[i + 1].fP[j - 1]);
      integralp1[j] = integralp1[j] + integralp1[j - 1];
    }
    // Add Extra points into distribution
    for (j = 0, count = 0; j < dists[i + 1].fLen; j++) {
      lo = TNudyCore::Instance()->BinarySearch(integral, dists[i].fLen, integralp1[j]);
      hi = lo + 1;
      f1 = dists[i].fP[lo];
      f2 = dists[i].fP[hi];
      e1 = dists[i].fX[lo];
      e2 = dists[i].fX[hi];
      ;
      delf = f2 - f1;
      if (fabs(delf) != 0) {
        dele = e2 - e1;
        temp = 2 * (integralp1[j] - integral[lo]) * delf / dele;
        if (integralp1[j] < integral[hi])
          x = e1 + dele * (sqrt(pow(f1, 2) + temp) - f1) / delf;
        else
          x = e2 + 1e-5;
      } else {
        if (f1 != 0)
          x = e1 + (integralp1[j] - integral[lo]) / f1;
        else
          x = e1 + 1e-5;
      }
      pp = TNudyCore::Instance()->LinearInterpolation(e1, f1, e2, f2, x);
      //      if(fabs(integralp1[j] - integral[lo]) > 1e-5){
      //** Work in progress to improve accuracy
      /*                 count++;
                  dists[i].fP->Set(dists[i].fLen+count);
                  dists[i].fX->Set(dists[i].fLen+count);
                  dists[i].fX->SetAt(x,dists[i].fLen+count-1);
                  dists[i].fP->SetAt(pp,dists[i].fLen+count-1);
 */ //}
    }
    dists[i].fLen += count;
    delete[] integral;
    integral    = new double[dists[i].fLen];
    integral[0] = 0;
    // Sort array
    for (j = 0; j < dists[i].fLen; j++) {
      min = j;
      for (k = j + 1; k < dists[i].fLen; k++) {
        if (dists[i].fX[k] < dists[i].fX[min]) {
          min = k;
        }
      }
      x                = dists[i].fP[min];
      dists[i].fP[min] = dists[i].fP[j];
      dists[i].fP[j]   = x;
      x                = dists[i].fX[min];
      dists[i].fX[min] = dists[i].fX[j];
      dists[i].fX[j]   = x;
    }

    // Recalculate integrals
    for (j = 1; j < dists[i].fLen; j++) {
      integral[j] = 0.5 * (dists[i].fX[j] - dists[i].fX[j - 1]) * (dists[i].fP[j] + dists[i].fP[j - 1]);
      integral[j] = integral[j] + integral[j - 1];
    }
    double *It = new double[dists[i].fLen - 1];
    double *I  = new double[dists[i].fLen - 1];
    for (int k = 0; k < dists[i].fLen; k++) {
      if (k > 0) {
        It[k - 1] = (dists[i].fX[k] - dists[i].fX[k - 1]) * (dists[i].fP[k] + dists[i].fP[k - 1]) * 0.5;
        I[k - 1]  = k - 1;
      }
    }
    delete dists[i].fChooseBin;
    dists[i].fChooseBin = new TNudyAlias(It, I, dists[i].fLen - 1, i);
    delete[] It;
    delete[] I;

    dists[i].fInterP = new double(dists[i].fLen);
    dists[i].fInterX = new double[dists[i].fLen];

    // Fill interpolation vector for ith entry
    dists[i].fInterX[0] = dists[i + 1].fX[0];
    dists[i].fInterP[0] = dists[i + 1].fP[0];

    for (j = 1; j < dists[i].fLen; j++) {
      lo = TNudyCore::Instance()->BinarySearch(integralp1, dists[i + 1].fLen, integral[j]);
      hi = lo + 1;
      f1 = dists[i + 1].fP[lo];
      f2 = dists[i + 1].fP[hi];
      e1 = dists[i + 1].fX[lo];
      e2 = dists[i + 1].fX[hi];
      ;
      delf = f2 - f1;
      if (fabs(delf) != 0) {
        dele = e2 - e1;
        temp = 2 * (integral[j] - integralp1[lo]) * delf / dele;
        if (integral[j] < integralp1[hi])
          x = e1 + dele * (sqrt(pow(f1, 2) + temp) - f1) / delf;
        else
          x = e2 + 1e-5;
      } else {
        if (f1 != 0)
          x = e1 + (integral[j] - integralp1[lo]) / f1;
        else
          x = e1 + 1e-5;
      }

      dists[i].fInterX[j] = x;
      pp                  = TNudyCore::Instance()->LinearInterpolation(e1, f1, e2, f2, x);
      if (pp < 0) pp      = 0;
      dists[i].fInterP[j] = pp;

      if (j > 0) {
        h += 0.5 * (dists[i].fInterX[j] - dists[i].fInterX[j - 1]) * (dists[i].fInterP[j] + dists[i].fInterP[j - 1]);
      }
    }
    for (j = 0; j < dists[i].fLen; j++)
      dists[i].fInterP[j] /= h;
    delete integral;
    delete integralp1;
  }
  dists[len - 1].fInterP     = nullptr;
  dists[len - 1].fInterX     = nullptr;
  dists[len - 1].fInterAlpha = 0;
}

#ifdef TNUDYALIAS_MULTITHREAD

//_______________________________________________________________________________
double *TNudyAliasCont::Randoms(int n, IntScheme_t iScheme)
{
  delete[] fMult;
  double *bins = fChooseBin->Randoms(n);
  int i;
  fMultLen = n;
  fMult    = new double[n];
  for (i     = 0; i < n; i++)
    fMult[i] = Random(iScheme, bins[i]);
  return fMult;
}

#endif
