#include "TNudyAlias.h"

#ifdef TNUDYALIAS_MULTITHREAD
#include <TThread.h>
#endif

#define ABS(X) X >= 0 ? X : -X

ClassImp(TNudyAlias)

    //_______________________________________________________________________________
    TNudyAlias::TNudyAlias()
    : fLen(0), fP(NULL), fX(NULL), fA(NULL), fR(NULL), fRnd(NULL), fMult(NULL), fMultLen(0) {}

//_______________________________________________________________________________
TNudyAlias::TNudyAlias(double *p, double *x, const int len, unsigned int seed)
    : fLen(len), fP(new double[fLen]), fX(new double[fLen]), fA(new double[fLen]), fR(new double[fLen]),
      fRnd(new TRandom(seed)), fMult(NULL), fMultLen(0) {
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
    Error("TNudyAlias", "Data not normalized, Integral = %e \n", sum);
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
      sum += ABS(b[j]);
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

//_______________________________________________________________________________
TNudyAlias::~TNudyAlias() {
  delete[] fP;
  delete[] fX;
  delete[] fA;
  delete[] fR;
  delete fRnd;
  fP = fX = fA = fR = NULL;
  fRnd = NULL;
  fLen = 0;
}

//_______________________________________________________________________________
void TNudyAlias::DumpTable() {
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

//_______________________________________________________________________________
double TNudyAlias::Random() {
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
