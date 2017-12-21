
#ifndef HIST_H
#define HIST_H

#include <iostream>

// IT IS JUST FOR TESTING

namespace userapplication {

class Hist{
public:

  Hist(double min, double max, int numbin) {
    fMin     = min;
    fMax     = max;
    fNumBins = numbin;
    fDelta = (fMax-fMin)/(numbin);
    fx = new double[fNumBins];
    fy = new double[fNumBins];
    for (int i=0; i<fNumBins; ++i) {
      fx[i] = fMin+i*fDelta;
      fy[i] = 0.0;
    }
    fSum = 0.0;
  }

  Hist(double min, double max, double delta) {
    fMin     = min;
    fMax     = max;
    fDelta   = delta;
    fNumBins = (int)((fMax-fMin)/(delta))+1.0;
    fx = new double[fNumBins];
    fy = new double[fNumBins];
    for (int i=0; i<fNumBins; ++i) {
      fx[i] = fMin+i*fDelta;
      fy[i] = 0.0;
    }
    fSum = 0.0;
  }

  void   Fill(double x) {
    int indx = (int)((x-fMin)/fDelta);
    if (indx<0) {
      std::cerr<<"\n ***** ERROR in Hist::FILL  =>  x = "<<x << " < fMin = "<<fMin<< std::endl;
      exit(1);
    }
  
    fy[indx] += 1.0;
   
  }


  void   Fill(double x, double w) {
    int indx = (int)((x-fMin)/fDelta);
    if (indx<0 || indx>fNumBins-1) {
      std::cerr<<"\n ***** ERROR in Hist::FILL  =>  x = "<<x << " < fMin = "<<fMin<<" or > fMaxn = "<<fMax<< std::endl;
      exit(1);
    }
    fy[indx] += 1.0*w;
  }

  int    GetNumBins() const { return fNumBins;}
  double GetDelta()   const { return fDelta;}
  double*  GetX() const {return fx;}
  double*  GetY() const {return fy;}


 private:
   double *fx;
   double *fy;
   double fMin;
   double fMax;
   double fDelta;
   double fSum;
   int    fNumBins;
};

}

#endif
