/**
 * @file TNudyAlias.cxx
 * @brief Implementation file for Alias Sampling
*/

#include "TNudyAlias.h"
#include "TRandom.h"

#ifdef TNUDYALIAS_MULTITHREAD
#include <TThread.h>
#endif

#ifdef USE_ROOT
ClassImp(TNudyAlias)
#endif

/**
 * Dummy constructor
 */
//_______________________________________________________________________________
TNudyAlias::TNudyAlias()
  : fLen(0), fP(nullptr), fX(nullptr), fA(nullptr), fR(nullptr), fRnd(nullptr), fMult(nullptr), fMultLen(0) {}

/**
 * @brief Full constructor
 * @param[in] p pointer to the discrete probability density function
 * @param[in] x pointer to the values of the random variable
 * @param[in] len number of points of the discrete probability function
 * @param[in] seed seed for the random number generator
 */
//_______________________________________________________________________________
TNudyAlias::TNudyAlias(double *p, double *x, int len, unsigned int seed)
    : fLen(len), fP(new double[fLen]), fX(new double[fLen]), fA(new double[fLen]), fR(new double[fLen]),
      fRnd(new TRandom(seed)), fMult(nullptr), fMultLen(0) {
  // Improve algorithm for building table
  int i, j;
  double sum, c, d, mean;
  int k, l;
  double *b = new double[len];
  mean = 1.0 / len;
  sum = c = d = k = l = 0;
  // Normalize
  for (i = 0; i < len; i++)
    sum += p[i];
  if (fabs(1.0 - sum) > ERROR_MARGIN) {
    printf("TNudyAlias::TNudyAlias: Data not normalized, Integral = %e \n", sum);
    for (i = 0; i < len; i++)
      p[i] /= sum;
  }
  // Build table
  for (i = 0; i < len; i++) {
    fP[i] = p[i];
    fX[i] = fA[i] = x[i];
    fR[i] = 0;
    b[i] = p[i] - mean;
  }
  for (i = 0; i < len; i++) {
    d = c = 0;
    for (j = 0; j < len; j++) {
      if (b[j] <= c) {
        c = b[j];
        k = j;
      } else if (b[j] >= d) {
        d = b[j];
        l = j;
      }
    }
    sum = 0;
    for (j = 0; j < len; j++)
       sum += std::abs(b[j]);
    if (sum < 1e-9) {
      break;
    } else {
      fA[k] = fX[l];
      fR[k] = 1 + c * (double)len;
      b[k] = 0;
      b[l] = c + d;
    }
  }
  delete[] b;
}

/** 
 * @brief Destructor
 */
//_______________________________________________________________________________
TNudyAlias::~TNudyAlias() {
  delete[] fP;
  delete[] fX;
  delete[] fA;
  delete[] fR;
  delete fRnd;
}

/**
 * @brief Dump the alias table
 */
//_______________________________________________________________________________
void TNudyAlias::DumpTable() const {
  int i, j;
  i = j = 0;
  // Reconstruct probability table
  double *prob = new double[fLen];
  for (i = 0; i < fLen; i++) {
    prob[i] = fR[i] / fLen;
    for (j = 0; j < fLen; j++) {
      if (fA[j] == fX[i])
        prob[i] += (1 - fR[j]) / fLen;
    }
  }
  for (i = 0; i < fLen; i++) {
    printf("p=%e, x=%e, r=%e, a=%e, error=%e \n", fP[i], fX[i], fR[i], fA[i], prob[i] - fP[i]);
  }
  delete[] prob;
}

/** 
 * Sample the distribution
 */
//_______________________________________________________________________________
double TNudyAlias::Random() const {
  double ua = fRnd->Uniform(1);
  double ub = fRnd->Uniform(1);
  int x = (int)(ua * fLen);
  double rx = fX[x];
  if (ub > fR[x])
    rx = fA[x];
  return rx;
}

#ifdef TNUDYALIAS_MULTITHREAD

//_______________________________________________________________________________
void *TNudyAlias::ThreadHandle(void *ptr) {
  TNudyComStruct *com = (TNudyComStruct *)ptr;
  if (!com->fAl->fMult) {
    delete com;
    return (void *)1;
  }
  int i = 0;
  int t = com->fAl->fMultLen / com->fAl->GetLen();
  for (i = com->fI * t; i < (com->fI + 1) * t; i++) {
    com->fAl->fMult[i] = com->fAl->Random();
  }
  delete com;
  return (void *)0;
}

//_______________________________________________________________________________
double *TNudyAlias::Randoms(int n) {
  int i;
  if (fMult)
    delete[] fMult;
  fMult = new double[n];
  fMultLen = n;
  void *(*funPtr)(void *) = &TNudyAlias::ThreadHandle;
  TThread **threads = new TThread *[fLen];
  for (i = 0; i < fLen; i++) {
    TNudyComStruct *thData = new TNudyComStruct(this, i);
    threads[i] = new TThread(Form("TNudyAlias%d", i), funPtr, (void *)thData);
    threads[i]->Run();
  }
  for (i = 0; i < fLen; i++) {
    threads[i]->Join();
  }
  for (i = 0; i < fLen; i++) {
    delete threads[i];
  }
  delete[] threads;
  return fMult;
}

#endif
