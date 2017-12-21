//
//
//   Author: J. Apostolakis, 28 Jan 2016

#ifndef STEPCHECKER_H
#define STEPCHECKER_H

#include "base/Vector3D.h"

#include "Geant/Config.h"

// namespace Geant {
// inline namespace MagField {

class StepChecker
{
  public:

   StepChecker( double eps, double maxLengthDiff= 0.0, bool verbose= false );
   ~StepChecker() {}

   /** @brief Compare two different 'solutions' - each provides a endpoint position/direction */
   bool
   CompareStep( vecgeom::Vector3D<double> const & position,
                vecgeom::Vector3D<double> const & direction,
                                   double         charge,
                                   double         momentum,
                                   double         step,
                vecgeom::Vector3D<double> const & endPosition1,
                vecgeom::Vector3D<double> const & endDirection1,
                vecgeom::Vector3D<double> const & endPositionRef,
                vecgeom::Vector3D<double> const & endDirectionRef ) const;

   /** @brief Check against helical solution */
   bool
   CheckStep( vecgeom::Vector3D<double> const & position,  
              vecgeom::Vector3D<double> const & direction,
                                 double         charge,
                                 double         momentum,
                                 double         step,
              vecgeom::Vector3D<double> const & endPosition,
              vecgeom::Vector3D<double> const & endDirection,
              vecgeom::Vector3D<double> const & BfieldVec ) const;

   // Potential to generalise methods above using template:
   // template<typename VectorType, typename BaseFltType, typename BaseIntType>
   //  in place of i) vecgeom::Vector3D<double> , ii) double,  iii) int
  private:
   double fMaxRelDiff;  /** max relative difference accepted - displacement or direction */
   double fMaxLengthDf; /** max displacement accepted without warning / error */
   bool   fVerbose;     /** Print comparison or not */
};

inline
StepChecker::StepChecker( double eps, double maxLengthDiff, bool verbose )
{
   fMaxRelDiff=   eps;
   fMaxLengthDf=  maxLengthDiff;
   fVerbose=      verbose;
}

//  } // End inline namespace MagField
// } // namespace Geant

#endif
