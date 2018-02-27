#ifndef TNudySampling_H
#define TNudySampling_H

#include "Geant/Particle.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TGraph.h"
#define PI acos(-1.0)
#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
class TNudyEndfRecoPoint;
#endif

class TNudySampling {

public:
  TNudySampling();
  TNudySampling(Particle *, TNudyEndfRecoPoint *recoPoint);
  virtual ~TNudySampling();

private:
  void GetSecParameter(Particle *, TNudyEndfRecoPoint *recoPoint);
  void FillHisto(double icosLab, double isecEnergyLab);
  std::vector<double> crs;
  double kineticE;
  double cosCM = 0, cosLab = 0, secEnergyCM = 0, secEnergyLab = 0;
  double x[1000000], y[1000000];
  double ene[1000000], nu1[1000000], nu2[1000000], nu3[1000000];
  double residueA, residueZ;
  int elemId   = 0;
  int isel     = 0;
  int counter  = 0;
  int ecounter = 0;
  int LCT, MF, MT, MF4, MF5, MF6;
  int events;
  TH2D *h;
  TH2D *hist[10];
  TH1D *h1;
  TH1D *h2;
  TH1D *fissZ1[10];
  TH1D *fissA1[10];
  TH1D *fissA;
  TGraph *gr1;
  TGraph *gr[5];
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudySampling, 1) // class for sampling
#endif
};
#endif
