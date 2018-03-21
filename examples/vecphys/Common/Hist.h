#ifndef HIST_H
#define HIST_H

#include <iostream>

// IT IS JUST FOR TESTING

class Hist {
public:
  Hist(double min, double max, int numbin)
  {
    fMin     = min;
    fMax     = max;
    fNumBins = numbin;
    fDelta   = (fMax - fMin) / (numbin);
    fx       = new double[fNumBins];
    fy       = new double[fNumBins];
    for (int i = 0; i < fNumBins; ++i) {
      fx[i] = fMin + i * fDelta;
      fy[i] = 0.0;
    }
    fSum = 0.0;
  }

  Hist(double min, double max, double delta)
  {
    fMin     = min;
    fMax     = max;
    fDelta   = delta;
    fNumBins = (int)((fMax - fMin) / (delta)) + 1.0;
    fx       = new double[fNumBins];
    fy       = new double[fNumBins];
    for (int i = 0; i < fNumBins; ++i) {
      fx[i] = fMin + i * fDelta;
      fy[i] = 0.0;
    }
    fSum = 0.0;
  }

  void Fill(double x)
  {
    int indx = (int)((x - fMin) / fDelta);
    if (indx < 0) {
      std::cerr << "\n ***** ERROR in Hist::FILL  =>  x = " << x << " < fMin = " << fMin << std::endl;
      exit(1);
    }

    fy[indx] += 1.0;
  }

  void Fill(double x, double w)
  {
    int indx = (int)((x - fMin) / fDelta);
    if (indx < 0) {
      std::cerr << "\n ***** ERROR in Hist::FILL  =>  x = " << x << " < fMin = " << fMin << std::endl;
      exit(1);
    }

    fy[indx] += 1.0 * w;
  }

  int GetNumBins() const { return fNumBins; }
  double GetDelta() const { return fDelta; }
  double *GetX() const { return fx; }
  double *GetY() const { return fy; }

  Hist operator/(Hist &b)
  {
    Hist &a = *this;
    assert(a.fMax == b.fMax);
    assert(a.fMin == b.fMin);
    assert(a.fNumBins = b.fNumBins);
    Hist res(a.fMin, a.fMax, a.fNumBins);
    for (int xi = 0; xi < fNumBins; ++xi) {
      res.fx[xi] = a.fx[xi];
      res.fy[xi] = a.fy[xi] / b.fy[xi];
      if (a.fy[xi] == 0.0 && b.fy[xi] == 0.0) {
        res.fy[xi] = 1.0;
      }
      res.fy[xi] = (res.fy[xi] - 1.0) * 100.0;
    }
    return res;
  }

  void Print()
  {
    for (int xi = 0; xi < fNumBins; ++xi) {
      Printf("%10f  --  %10f", fx[xi], fy[xi]);
    }
  }

private:
  double *fx;
  double *fy;
  double fMin;
  double fMax;
  double fDelta;
  double fSum;
  int fNumBins;
};

#endif
