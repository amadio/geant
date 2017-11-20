#ifndef XSECTIONSVECTOR_H
#define XSECTIONSVECTOR_H

#include "EMModel.h"

namespace geantphysics {

class Spline;
/**
 * @brief   Class to handle tabulated cross-sections data -> one object per element for total cross-sections and one
 * object per subshell for partial cross-sections data.
 * @class   XSectionsVector
 * @author  M Bandieramonte
 * @date    September 2017
 *
 *
 * \cite
 */

class XSectionsVector {

public:
  XSectionsVector();
  ~XSectionsVector();

  size_t FindCSBinLocation(double energy, size_t idx) const;
  double GetValueAt(double energy) const;
  double GetValue(double energy, size_t &shellIdx) const;
  inline double LinearInterpolation(double energy, size_t idx) const;

  std::vector<double> fBinVector;  // Cross sections bin vector (i.e. x coordinate)
  std::vector<double> fDataVector; // Cross sections data vector (i.e. y coordinate)

  size_t fNumberOfNodes; // Number of elements
  double fEdgeMin;       // Energy of first point
  double fEdgeMax;       // Energy of last point
  Spline *fSplineInt;    // Spline interpolator
  int fCompID;           // id of the object
};
} // namespace geantphysics

#endif // XSECTIONSVECTOR_H
