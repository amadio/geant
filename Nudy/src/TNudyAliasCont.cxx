#include "TNudyAliasCont.h"

ClassImp(TNudyAliasCont)

    //_______________________________________________________________________________
    TNudyAliasCont::TNudyAliasCont()
    : fLen(0), fChooseBin(NULL), fP(NULL), fX(NULL), fInterX(NULL), fInterP(NULL), fTx(NULL), fTp(NULL), fTa(-1),
      fInterAlpha(10), fRan(NULL), fRnd(NULL), fAlpha(10)
#ifdef TNUDYALIAS_MULTITHREAD
      ,
      fMult(NULL), fMultLen(0)
#endif
{
}

//_______________________________________________________________________________
TNudyAliasCont::TNudyAliasCont(TArrayD *data, double alpha, unsigned int seed)
    : fLen(0), fChooseBin(NULL), fP(NULL), fX(NULL), fInterX(NULL), fInterP(NULL), fTx(NULL), fTp(NULL), fTa(-1),
      fInterAlpha(10), fRan(NULL), fRnd(NULL), fAlpha(10)
#ifdef TNUDYALIAS_MULTITHREAD
      ,
      fMult(NULL), fMultLen(0)
#endif
{
  // TArrayD of x1,p1,x2,p2.....xn,pn
  int len = data->GetSize();
  int i = 0, j = 0;
  if (len % 2 != 0) {
    Error("TNudyAliasCont", "Incorrect TArrayD specified");
  }
  double *x = new double[len / 2];
  double *p = new double[len / 2];
  for (i = 0, j = 0; i < len; i += 2, j++) {
    x[j] = data->GetAt(i);
    p[j] = data->GetAt(i + 1);
  }
  Initialize(p, x, len, alpha, seed);
  delete[] p;
  p = NULL;
  delete[] x;
  x = NULL;
}

//_______________________________________________________________________________
TNudyAliasCont::TNudyAliasCont(double *p, double *x, const int len, double alpha, unsigned int seed)
   : fLen(0), fChooseBin(NULL), fP(NULL), fX(NULL), fInterX(NULL), fInterP(NULL), fTx(NULL), fTp(NULL), fTa(-1),
   fInterAlpha(10), fRan(NULL), fRnd(NULL), fAlpha(10)
#ifdef TNUDYALIAS_MULTITHREAD
   ,
   fMult(NULL), fMultLen(0)
#endif
{
  Initialize(p, x, len, alpha, seed);
}

//_______________________________________________________________________________
void TNudyAliasCont::Initialize(double *p, double *x, const int len, double alpha, unsigned int seed) {
  // Data sorted(asc) in x
  // alpha to be used when using statistical interpolation between two distributions
  int i;
  double *integral = NULL, *I = NULL;
  fInterX = fInterP = NULL;
  //  fTx = new TArrayD(len);
  //  fTp = new TArrayD(len);
  fTa = -1;
  fInterAlpha = 0;
#ifdef TNUDYALIAS_MULTITHREAD
  fMult = NULL;
  fMultLen = 0;
#endif
  fLen = len;
  fAlpha = alpha;
  fRnd = new TRandom(seed);
  fRan = new TRandom(seed >> 2);
  fP = new TArrayD(len);
  fX = new TArrayD(len);
  integral = new double[len - 1];
  I = new double[len - 1];
  integral[0] = 0;
  I[0] = 0;
  double h = 0;
  for (i = 0; i < len; i++) {
    fP->SetAt(p[i], i);
    fX->SetAt(x[i], i);
    if (i > 0)
      h += (x[i] - x[i - 1]) * (p[i] + p[i - 1]) * 0.5;
  }
  for (i = 0; i < len; i++) {
    fP->SetAt(fP->At(i) / h, i);
    if (i > 0) {
      integral[i - 1] = (fX->At(i) - fX->At(i - 1)) * (fP->At(i) + fP->At(i - 1)) * 0.5;
      I[i - 1] = i - 1;
    }
  }
  fChooseBin = new TNudyAlias(integral, I, len - 1, seed);
  delete[] integral;
  delete[] I;
}

//_______________________________________________________________________________
TNudyAliasCont::~TNudyAliasCont() {
  delete fChooseBin;
  delete fP;
  delete fX;
  delete fInterX;
  delete fInterP;
  if (fTx)
    delete fTx;
  if (fTp)
    delete fTp;
  delete fRnd;
  delete fRan;
  fP = fX = fInterX = fInterP = NULL;
  fChooseBin = NULL;
  fRnd = fRan = NULL;
  fLen = 0;
#ifdef TNUDYALIAS_MULTITHREAD
  delete[] fMult;
  fMult = NULL;
  fMultLen = 0;
#endif
}

//_______________________________________________________________________________
void TNudyAliasCont::DumpTable() {
  int i;
  printf("Continuous probability distribution --> \n");
  for (i = 0; i < fLen; i++) {
    printf("x[%d]=%e, p[%d]=%e\n", i, fX->At(i), i, fP->At(i));
  }
  if (fInterX && fInterP)
    for (i = 0; i < fLen; i++)
      printf("fIx[%i]=%e, fIp[%d]=%e\n", i, fInterX->At(i), i, fInterP->At(i));
  printf("Alias table -->\n");
  fChooseBin->DumpTable();
}

//_______________________________________________________________________________
double TNudyAliasCont::Random(IntScheme_t iScheme, int binNo, AliasDist_t distribution, double r1, double r2) {
  if (binNo < 0 || binNo >= fLen) // Choose bins if required as per their probabilities
    binNo = (int)fChooseBin->Random();
  double *xx;
  double *xp;
  switch (distribution) {
  case kOriginal:
    xx = fX->GetArray();
    xp = fP->GetArray();
    break;
  case kBuilt:
    xx = (fInterX && fInterP) ? fInterX->GetArray() : fX->GetArray();
    xp = (fInterX && fInterP) ? fInterP->GetArray() : fP->GetArray();
  }
  double rnd1, rnd2, x1, x2;
  switch (iScheme) {
  case kLinear:
    // Coerce linear probabiliy distributions to equal probability bins
    rnd1 = r1 < 0 ? fRnd->Uniform() : r1;
    rnd2 = r2 < 0 ? fRan->Uniform() : r2;
    x1 = (1 - rnd1) * xx[binNo] + rnd1 * xx[binNo + 1];
    x2 = rnd1 * xx[binNo] + (1 - rnd1) * xx[binNo + 1];
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
double TNudyAliasCont::StatisticalInterpolation(TNudyAliasCont *dist1, TNudyAliasCont *dist2, double alpha) {
  // Simple statistical interpolation between distributions ( Not very accurate )
  if (dist1->Uniform(1) <= ((dist2->fAlpha - alpha) / (dist2->fAlpha - dist1->fAlpha)))
    return dist1->Random();
  else
    return dist2->Random();
}

//_______________________________________________________________________________
double TNudyAliasCont::ImprovedInterpolation(double alpha) {
  // Improved interpolation between distributions - Interpolates between corresonding values in the
  // intergral of probability density
  double x1, x2, X1, X2, P1, P2, rnd1, rnd2, u;
  u = (alpha - fAlpha) / (fInterAlpha - fAlpha);
  if (alpha < fAlpha || alpha > fInterAlpha) {
    Error("ImprovedInterpolation", "%e does not lie between %e and %e", alpha, fAlpha, fInterAlpha);
  }
  double h = 0;
  if (!fTx)
    fTx = new TArrayD(fLen);
  if (!fTp)
    fTp = new TArrayD(fLen);
  if (alpha != fTa) {
    fTa = alpha;
    for (int i = 0; i < fLen; i++) {
      fTx->SetAt((1 - u) * fX->At(i) + u * fInterX->At(i), i);
      fTp->SetAt((1 - u) * fP->At(i) + u * fInterP->At(i), i);
      if (i > 0) {
        h += 0.5 * (fTx->At(i) - fTx->At(i - 1)) * (fTp->At(i) + fTp->At(i - 1));
      }
    }
    if (fabs(1.0 - h) > 1e-5) {
      for (int i = 0; i < fLen; i++) {
        fTp->SetAt(fTp->At(i) / h, i);
      }
    }
  }
  int binNo = (int)fChooseBin->Random();
  X1 = fTx->At(binNo);
  X2 = fTx->At(binNo + 1);
  P1 = fTp->At(binNo);
  P2 = fTp->At(binNo + 1);
  rnd1 = fRnd->Uniform();
  rnd2 = fRan->Uniform();
  x1 = (1.0 - rnd1) * X1 + rnd1 * X2;
  x2 = rnd1 * X1 + (1.0 - rnd1) * X2;
  if (rnd2 * (P1 + P2) <= (1.0 - rnd1) * P1 + rnd1 * P2)
    return x1;
  else
    return x2;
}
//_______________________________________________________________________________
double TNudyAliasCont::SelectBin(double alpha, AliasDist_t distribution) {
  int binNo = (int)fChooseBin->Random();
  if (distribution == kOriginal)
    return fX->At(binNo);
  else
    return fInterX->At(binNo);
}
//_______________________________________________________________________________
void TNudyAliasCont::BuildIntermediate(TNudyAliasCont *dists, const int len) {
  // Populates data structures to interpolate between distributions
  int i, j, k;
  TArrayD *integral;
  TArrayD *integralp1;
  double x;
  int lo, hi, mid, min;
  int count;
  double f1, f2, e1, e2, dele, delf, temp, pp;
  double h = 0;
  lo = hi = mid = 0;
  for (i = 0; i < len - 1; i++) {
    dists[i].fInterAlpha = dists[i + 1].fAlpha;
    integral = new TArrayD(dists[i].fLen);       // Integral for dists[i]
    integralp1 = new TArrayD(dists[i + 1].fLen); // Integral for dists[i+1]
    integral->SetAt(0, 0);
    integralp1->SetAt(0, 0);
    // Evaluate integrals
    for (j = 1; j < dists[i].fLen; j++) {
      integral->SetAt(
          0.5 * (dists[i].fX->At(j) - dists[i].fX->At(j - 1)) * (dists[i].fP->At(j) + dists[i].fP->At(j - 1)), j);
      integral->SetAt(integral->At(j) + integral->At(j - 1), j);
    }
    // Store bin areas of
    for (j = 1; j < dists[i + 1].fLen; j++) {
      integralp1->SetAt(0.5 * (dists[i + 1].fX->At(j) - dists[i + 1].fX->At(j - 1)) *
                            (dists[i + 1].fP->At(j) + dists[i + 1].fP->At(j - 1)),
                        j);
      integralp1->SetAt(integralp1->At(j) + integralp1->At(j - 1), j);
    }
    // Add Extra points into distribution
    for (j = 0, count = 0; j < dists[i + 1].fLen; j++) {
      lo = TNudyCore::Instance()->BinarySearch(integral->GetArray(), dists[i].fLen, integralp1->At(j));
      hi = lo + 1;
      f1 = dists[i].fP->At(lo);
      f2 = dists[i].fP->At(hi);
      e1 = dists[i].fX->At(lo);
      e2 = dists[i].fX->At(hi);
      ;
      delf = f2 - f1;
      if (fabs(delf) != 0) {
        dele = e2 - e1;
        temp = 2 * (integralp1->At(j) - integral->At(lo)) * delf / dele;
        if (integralp1->At(j) < integral->At(hi))
          x = e1 + dele * (sqrt(pow(f1, 2) + temp) - f1) / delf;
        else
          x = e2 + 1e-5;
      } else {
        if (f1 != 0)
          x = e1 + (integralp1->At(j) - integral->At(lo)) / f1;
        else
          x = e1 + 1e-5;
      }
      pp = TNudyCore::Instance()->LinearInterpolation(e1, f1, e2, f2, x);
      //      if(fabs(integralp1->At(j) - integral->At(lo)) > 1e-5){
      //** Work in progress to improve accuracy
      /*                 count++;
                  dists[i].fP->Set(dists[i].fLen+count);
                  dists[i].fX->Set(dists[i].fLen+count);
                  dists[i].fX->SetAt(x,dists[i].fLen+count-1);
                  dists[i].fP->SetAt(pp,dists[i].fLen+count-1);
 */ //}
    }
    dists[i].fLen += count;
    integral->Set(dists[i].fLen);
    // Sort array
    for (j = 0; j < dists[i].fLen; j++) {
      min = j;
      for (k = j + 1; k < dists[i].fLen; k++) {
        if (dists[i].fX->At(k) < dists[i].fX->At(min)) {
          min = k;
        }
      }
      x = dists[i].fP->At(min);
      dists[i].fP->SetAt(dists[i].fP->At(j), min);
      dists[i].fP->SetAt(x, j);
      x = dists[i].fX->At(min);
      dists[i].fX->SetAt(dists[i].fX->At(j), min);
      dists[i].fX->SetAt(x, j);
    }

    // Recalculate integrals
    for (j = 1; j < dists[i].fLen; j++) {
      integral->SetAt(
          0.5 * (dists[i].fX->At(j) - dists[i].fX->At(j - 1)) * (dists[i].fP->At(j) + dists[i].fP->At(j - 1)), j);
      integral->SetAt(integral->At(j) + integral->At(j - 1), j);
    }
    double *It = new double[dists[i].fLen - 1];
    double *I = new double[dists[i].fLen - 1];
    for (int k = 0; k < dists[i].fLen; k++) {
      if (k > 0) {
        It[k - 1] = (dists[i].fX->At(k) - dists[i].fX->At(k - 1)) * (dists[i].fP->At(k) + dists[i].fP->At(k - 1)) * 0.5;
        I[k - 1] = k - 1;
      }
    }
    delete dists[i].fChooseBin;
    dists[i].fChooseBin = new TNudyAlias(It, I, dists[i].fLen - 1, i);
    delete[] It;
    delete[] I;

    dists[i].fInterP = new TArrayD(dists[i].fLen);
    dists[i].fInterX = new TArrayD(dists[i].fLen);

    // Fill interpolation vector for ith entry
    dists[i].fInterX->SetAt(dists[i + 1].fX->At(0), 0);
    dists[i].fInterP->SetAt(dists[i + 1].fP->At(0), 0);

    for (j = 1; j < dists[i].fLen; j++) {
      lo = TNudyCore::Instance()->BinarySearch(integralp1->GetArray(), dists[i + 1].fLen, integral->At(j));
      hi = lo + 1;
      f1 = dists[i + 1].fP->At(lo);
      f2 = dists[i + 1].fP->At(hi);
      e1 = dists[i + 1].fX->At(lo);
      e2 = dists[i + 1].fX->At(hi);
      ;
      delf = f2 - f1;
      if (fabs(delf) != 0) {
        dele = e2 - e1;
        temp = 2 * (integral->At(j) - integralp1->At(lo)) * delf / dele;
        if (integral->At(j) < integralp1->At(hi))
          x = e1 + dele * (sqrt(pow(f1, 2) + temp) - f1) / delf;
        else
          x = e2 + 1e-5;
      } else {
        if (f1 != 0)
          x = e1 + (integral->At(j) - integralp1->At(lo)) / f1;
        else
          x = e1 + 1e-5;
      }

      dists[i].fInterX->SetAt(x, j);
      pp = TNudyCore::Instance()->LinearInterpolation(e1, f1, e2, f2, x);
      if (pp < 0)
        pp = 0;
      dists[i].fInterP->SetAt(pp, j);

      if (j > 0) {
        h += 0.5 * (dists[i].fInterX->At(j) - dists[i].fInterX->At(j - 1)) *
             (dists[i].fInterP->At(j) + dists[i].fInterP->At(j - 1));
      }
    }
    for (j = 0; j < dists[i].fLen; j++)
      dists[i].fInterP->SetAt(dists[i].fInterP->At(j) / h, j);
    delete integral;
    delete integralp1;
  }
  dists[len - 1].fInterP = NULL;
  dists[len - 1].fInterX = NULL;
  dists[len - 1].fInterAlpha = 0;
}

#ifdef TNUDYALIAS_MULTITHREAD

//_______________________________________________________________________________
double *TNudyAliasCont::Randoms(int n, IntScheme_t iScheme) {
  if (fMult)
    delete[] fMult;
  double *bins = fChooseBin->Randoms(n);
  int i;
  fMultLen = n;
  fMult = new double[n];
  for (i = 0; i < n; i++)
    fMult[i] = Random(iScheme, bins[i]);
  return fMult;
}

#endif
