#include "TNudyAlias.h"

#ifdef TNUDYALIAS_MULTITHREAD
#include <TThread.h>
#endif

ClassImp(TNudyAlias)

    //_______________________________________________________________________________
    TNudyAlias::TNudyAlias() {
  fLen = fMultLen = 0;
  fMult = fP = fX = fR = fA = NULL;
  fRnd = NULL;
}

//_______________________________________________________________________________
TNudyAlias::TNudyAlias(Double_t *p, Double_t *x, const Int_t len, UInt_t seed) {
  // Improve algorithm for building table
  int i, j;
  Double_t sum, c, d, mean;
  Int_t k, l;
  Double_t *b = new Double_t[len];
  fP = new Double_t[len];
  fA = new Double_t[len];
  fR = new Double_t[len];
  fX = new Double_t[len];
  fMult = NULL;
  fMultLen = 0;
  mean = 1.0 / len;
  sum = c = d = k = l = 0;
  fRnd = new TRandom(seed);
  fLen = len;
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
      sum += fabs(b[j]);
    if (sum < 1e-9) {
      break;
    } else {
      fA[k] = fX[l];
      fR[k] = 1 + c * (Double_t)len;
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
  Double_t *prob = new Double_t[fLen];
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
Double_t TNudyAlias::Random() {
  Double_t ua = fRnd->Uniform(1);
  Double_t ub = fRnd->Uniform(1);
  Int_t x = (int)(ua * fLen);
  Double_t rx = fX[x];
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
  Int_t i = 0;
  Int_t t = com->fAl->fMultLen / com->fAl->GetLen();
  for (i = com->fI * t; i < (com->fI + 1) * t; i++) {
    com->fAl->fMult[i] = com->fAl->Random();
  }
  delete com;
  return (void *)0;
}

//_______________________________________________________________________________
Double_t *TNudyAlias::Randoms(Int_t n) {
  int i;
  if (fMult)
    delete[] fMult;
  fMult = new Double_t[n];
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
  return fMult;
}

#endif
