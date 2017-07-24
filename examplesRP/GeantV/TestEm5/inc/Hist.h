
#ifndef HIST_H
#define HIST_H

#include <iostream>

namespace userapplication {

/**
 * @brief Simple histogram object.
 *
 * @class   Hist
 * @author  M Novak
 * @date    june 2016
 */

class Hist {
public:

  Hist(double min, double max, int numbin);
  Hist(double min, double max, double delta);

  void Fill(double x);
  void Fill(double x, double w);

  int      GetNumBins() const { return fNumBins;}
  double   GetDelta()   const { return fDelta;}
  double*  GetX() const {return fx;}
  double*  GetY() const {return fy;}

  Hist& operator+=(const Hist &other);

 private:
   double *fx;
   double *fy;
   double fMin;
   double fMax;
   double fDelta;
   int    fNumBins;
};

}      // namespace userapplication

#endif // HIST_H
