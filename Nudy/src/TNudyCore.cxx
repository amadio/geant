#include "TNudyCore.h"
#include "TROOT.h"
#include <stdlib.h>
using std::max;

ClassImp(TNudyCore)

    TNudyCore *TNudyCore::fgInstance = 0;

//______________________________________________________________________________
TNudyCore::~TNudyCore() {
  delete fGeom;
  fGeom = NULL;
  fTable = NULL;
  delete fPdgDB;
  fPdgDB = NULL;
  fListOfObjects->Delete();
  delete fListOfObjects;
  fListOfObjects = NULL;
  gROOT->GetListOfSpecials()->Remove(this);
  fgInstance = NULL;
}

//______________________________________________________________________________
TNudyCore::TNudyCore() : TNamed("NudyCore", "Core of Nudy ENDF Framework") {
  fGeom = new TGeoManager("NudyGeoManager", "Geometry Manager for the TNudy Framework");
  fTable = fGeom->GetElementTable();
  fTable->BuildDefaultElements();
  fTable->ImportElementsRN();
  fPdgDB = new TDatabasePDG();
  fPdgDB->ReadPDGTable();
  fListOfObjects = new TList();
  if (fgInstance) {
    Warning("TNudyCore", "object already instantiated");
  } else {
    fgInstance = this;
    gROOT->GetListOfSpecials()->Add(this);
  }
}

//______________________________________________________________________________
TNudyCore *TNudyCore::Instance() {
  if (!fgInstance)
    fgInstance = new TNudyCore();
  return fgInstance;
}

//______________________________________________________________________________
char *TNudyCore::ExpandReaction(Reaction_t reac) {
  switch (reac) {

  default:
    return Form("%d", (int)reac);
  }
}

//______________________________________________________________________________
TParticle *TNudyCore::GetParticle(int pdgCode) {
  TParticle *particle = new TParticle();
  particle->SetPdgCode(pdgCode);
  return particle;
}

//______________________________________________________________________________
TParticle *TNudyCore::GetParticle(const char *name) {
  TParticle *particle = new TParticle();
  particle->SetPdgCode(fPdgDB->GetParticle(name)->PdgCode());
  return particle;
}

//_______________________________________________________________________________
TParticlePDG *TNudyCore::GetParticlePDG(int pdgCode) { return fPdgDB->GetParticle(pdgCode); }

//______________________________________________________________________________
TParticlePDG *TNudyCore::GetParticlePDG(const char *name) { return fPdgDB->GetParticle(name); }

//______________________________________________________________________________
int TNudyCore::IsMaterial(const TGeoElementRN *endf, const char *key) {
  if (endf == NULL)
    return 1;
  char *matkey = Form("%07d", endf->ENDFCode());
  const char *pos = &key[9];
  while (*pos == matkey[pos - key - 9]) {
    pos++;
  }
  if (pos < key + 16)
    return 0;
  else
    return 1;
}

//______________________________________________________________________________
int TNudyCore::IsTemperature(const ULong_t temp, const char *key) {
  if (temp == 0)
    return 1;
  char *matkey = Form("%12ld", temp);
  const char *pos = key;
  while (*pos == matkey[pos - key]) {
    pos++;
  }
  if (pos < key + 6)
    return 0;
  else
    return 1;
}

//______________________________________________________________________________
int TNudyCore::IsReaction(const Reaction_t r, const char *key) {
  if (r == kNoReaction)
    return 1;
  char *matkey = Form("%03d", r);
  const char *pos = &key[6];
  while (*pos == matkey[pos - key - 6]) {
    pos++;
  }
  if (pos < key + 9)
    return 0;
  else
    return 1;
}

//______________________________________________________________________________
void TNudyCore::MemProfile() {
  MemInfo_t mem;
  gSystem->GetMemInfo(&mem);
  printf("RAM\nTotal: %d\nUsed: %d\nFree: %d\nSWAP\nTotal: %d\nUsed: %d\nFree: %d\n", mem.fMemTotal, mem.fMemUsed,
         mem.fMemFree, mem.fSwapTotal, mem.fSwapUsed, mem.fSwapFree);
}

//______________________________________________________________________________
char *TNudyCore::GetKey(const TGeoElementRN *mat, Reaction_t reac, ULong_t temp) {
  return Form("%012ld%03d%07d", temp, reac, mat->ENDFCode());
}

//______________________________________________________________________________
int TNudyCore::BinarySearch(double *array, int len, double val) {
  int min = 0;
  int max = len - 1;
  int mid = 0;
  if (val <= array[min])
    return 0;
  else if (val >= array[max])
    return max - 1;
  else {
    while (max - min > 1) {
      mid = (min + max) / 2;
      if (val < array[mid])
        max = mid;
      else
        min = mid;
    }
  }
  return min;
}

//______________________________________________________________________________
double TNudyCore::InterpolateScale(double x[2], double y[2], int law, double xx) {
  double yy = -1;
  double small = 1e-20;
  //  Info("InterpolateScale","x1,y1 = %e,%e x2,y2 = %e,%e INT = %d xx = %e",x[0],y[0],x[1],y[1],law,xx);
  if (law == 1 || x[1] <= x[0]) {
    yy = y[0];
    if (xx == x[1])
      yy = y[1];
  } else if (law == 2) {
    yy = y[0] + (y[1] - y[0]) * (xx - x[0]) / (x[1] - x[0]);
  } else if (law == 3) {
    x[0] = max<double>(x[0], small);
    x[1] = max<double>(x[1], small);
    yy = y[0] + (y[1] - y[0]) * log(xx / x[0]) / log(x[1] / x[0]);
  } else if (law == 4) {
    y[0] = max<double>(y[0], small);
    y[1] = max<double>(y[1], small);
    yy = exp(log(y[0]) + log(y[1] / y[0]) * (xx - x[0]) / (x[1] - x[0]));
  } else if (law == 5) {
    x[0] = max<double>(x[0], small);
    x[1] = max<double>(x[1], small);
    y[0] = max<double>(y[0], small);
    y[1] = max<double>(y[1], small);
    yy = exp(log(y[0]) + log(y[1] / y[0]) * log(xx / x[0]) / log(x[1] / x[0]));
  }
  return yy;
}

//______________________________________________________________________________
double TNudyCore::Interpolate(int *nbt, int *interp, int nr, double *x, double *y, int np, double xx) {
  double yy = 0;
  //  Info("Interpolation","E = %e:%e, xx = %e , P = %e:%e",x[0],x[np-1],xx,y[0],y[np-1]);
  //  Info("Interpolation limits","xx = %e min = %e max = %e",xx,x[0],x[np-1]);
  if (xx < x[0]) {
    yy = 0;
    return yy;
  } else if (xx > x[np - 1]) {
    yy = 0;
    return yy;
  } else if (xx == x[np - 1]) {
    yy = y[np - 1];
    return yy;
  }
  int index = BinarySearch(x, np, xx);
  if (xx < x[index] || xx > x[index + 1]) {
    Error("Interpolate", "Error in the interpolation xx = %e does not lie between %e and %e Index = %d", xx, x[index],
          x[index + 1], index);
    return 0;
  }
  int intlaw = 0;
  for (int jnt = 1; jnt <= nr; jnt++) {
    if (index < nbt[jnt - 1]) {
      intlaw = interp[jnt - 1];
      yy = InterpolateScale(x + index, y + index, intlaw, xx);
    }
  }
  //  Info("Value","Interpolated Value = %e",yy);
  return yy;
}

//______________________________________________________________________________
double TNudyCore::LinearInterpolation(double x1, double y1, double x2, double y2, double x) {
  if (x2 == x1) {
    //   Error("TNudyCore::LinearInterpolation","Points specified have same coordinate in X");
    return y1;
  }
  return y1 + ((y2 - y1) / (x2 - x1)) * (x - x1);
}

//______________________________________________________________________________
double TNudyCore::BilinearInterploation(double x1, double y1, double x2, double y2, double z11,
                                          double z12, double z21, double z22, double x, double y) {
  if ((x2 == x1) || (y2 == y1)) {
    Error("TNudyCore::BilinearInterpolation", "Points specified have same coordinate in X or Y");
    return 0;
  }
  // Storing partial sums and products so that they dont have to be recalculated
  double d = (x2 - x1) * (y2 - y1);
  double n2x = x2 - x;
  double n2y = y2 - y;
  double n1x = x - x1;
  double n1y = y - y1;
  // printf("x1 = %e, y1 = %e, x2 = %e, y2 = %e\n",x1,y1,x2,y2);
  // printf("%e, %e, %e, %e, %e\n",d,n2x,n2y,n1x,n1y);
  return (z11 * n2x * n2y + z21 * n1x * n2y + z12 * n2x * n1y + z22 * n1x * n1y) / d;
}
//______________________________________________________________________________
void TNudyCore::CumulativeIntegral(double *x, double *y, double *q, int len) {
  for (int i = 0; i < len; i++) {
    if (i > 0) {
      q[i - 1] = 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1]);
      if (i > 1)
        q[i - 1] += q[i - 2];
    }
  }
}
//______________________________________________________________________________
void TNudyCore::TrapezoidalIntegral(double *xpts, double *ypts, const int npts, double *out) {
  // This function evaluates the integral of discrete points using the trapezoidal rule
  // and returns the value of the integral at the same points in x
  if (!xpts || !ypts || (npts == 0)) {
    Error("TrapezoidalIntegral()", "Incorrect input");
    return;
  }
  if (out == NULL) {
    out = new double[npts];
  }
  out[0] = 0;
  for (int i = 1; i < npts; i++) {
    out[i] = out[i - 1] + (xpts[i] - xpts[i - 1]) * (ypts[i] + ypts[i - 1]) * 0.5;
  }
}
