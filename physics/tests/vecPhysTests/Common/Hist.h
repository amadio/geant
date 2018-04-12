#ifndef HIST_H
#define HIST_H

#include <iostream>
#include <cmath>

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

  /*
   * Compare this histogram with other one, assuming N[bin] ~ Poisson(N_observed[bin])
   * Then Mean[ N[bin] ] = N_observed[bin]
   * Variance[ N[bin] ] = Sqrt[ N_observed[bin] ]
   *
   * Sigma of difference N_hist1[bin] - N_hist2[bin]  approx. = Sqrt[N_hist1 + N_hist2]
   *
   * Reports #bins with difference less than 4 sigma(defined above), and >= 4 sigma
   */
  void Compare(Hist &hist)
  {
    int fourSigmaDev = 0, bigDev = 0;
    Hist &a = *this;
    Hist &b = hist;
    assert(a.fMax == b.fMax);
    assert(a.fMin == b.fMin);
    assert(a.fNumBins = b.fNumBins);
    Hist diff(a.fMin, a.fMax, a.fNumBins);
    Hist diffSigma(a.fMin, a.fMax, a.fNumBins);
    for (int xi = 0; xi < fNumBins; ++xi) {
      diff.fx[xi]      = a.fx[xi];
      diffSigma.fx[xi] = a.fx[xi];
      if (a.fy[xi] == 0.0 && b.fy[xi] == 0.0) {
        diff.fy[xi]      = 0.0;
        diffSigma.fy[xi] = 0.0;
      } else {
        diff.fy[xi]      = ComputeDiff(a.fy[xi], b.fy[xi]);
        diffSigma.fy[xi] = ComputeDiffSigma(a.fy[xi], b.fy[xi]);
      }

      if (std::abs(diff.fy[xi]) <= 4.0 * diffSigma.fy[xi]) {
        fourSigmaDev++;
      } else {
        Printf("Big deviation in bin x: %f H1: %f H2: %f dev: %f dev_sigma %f", diff.fx[xi], a.fy[xi], b.fy[xi],
               diff.fy[xi], diffSigma.fy[xi]);
        bigDev++;
      }
    }
    Printf("4 sigma dev: %d | big dev: %d", fourSigmaDev, bigDev);
  }

  void Print()
  {
    for (int xi = 0; xi < fNumBins; ++xi) {
      Printf("%10f  --  %10f", fx[xi], fy[xi]);
    }
  }

private:
  double ComputeDiff(double a, double b)
  {
    assert(a >= 0.0);
    assert(b >= 0.0);
    return a - b;
  }
  double ComputeDiffSigma(double a, double b)
  {
    assert(a >= 0.0);
    assert(b >= 0.0);
    return std::sqrt(a + b);
  }

  double *fx;
  double *fy;
  double fMin;
  double fMax;
  double fDelta;
  double fSum;
  int fNumBins;
};

#endif
