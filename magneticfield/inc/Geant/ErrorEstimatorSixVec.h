//
//  ErrorEstimatorSixVec: Simple class to Estimate the overall error for each lane
//    
//
//  Created by japost on 12.12.2018
//

#ifndef ErrorEstimatorSixVec_h
#define ErrorEstimatorSixVec_h

// class ScalarFieldTrack;
// struct FieldTrack;
// #include "FieldTrack.h"

//  Attempt to provide a method that works for a single track using Flex/Simple Integration Driver
#define EXTEND_SINGLE 1

class ErrorEstimatorSixVec {
  public:
    static constexpr int    fNoComponents= 6;

    ErrorEstimatorSixVec( double eps_rel_max, double minimumStep ) :
        fEpsRelMax (eps_rel_max ),
        fInvEpsilonRelSq( 1.0 / (eps_rel_max * eps_rel_max) ),
        fMinimumStep( minimumStep )
    {
    } 

    template <class Real_v>
      Real_v EstimateError( const Real_v yError[fNoComponents],
                            const Real_v hStep,
                            // const Real_v yValue[fNoComponents],
                            const Real_v magMomentumSq //   Initial momentum square (used for rel. error)
         ) const;
    //  Returns the Maximum Square Error: the maximum between 
    //    - position relative error square (magnitude^2): (momentum error vec)^2 / (initial momentum)^2
    //    - momentum relative error square (magnitude^2): (position error vec)^2 / (step length)^2
    //
    //  Last argument enables the use of initial momentum square in calculating the relative error

    double GetMaxRelativeError() const { return fEpsRelMax; }

  public:
    static constexpr double tinyValue = 1.0e-80; // Just to ensure there is no division by zero
    
  private:
    const double fEpsRelMax;
    const double fInvEpsilonRelSq; // = 1.0 / (eps_rel_max * eps_rel_max);
    const double fMinimumStep;
};

template <class Real_v> 
   Real_v ErrorEstimatorSixVec::EstimateError(
             const Real_v yEstError[fNoComponents],
             const Real_v hStep,
             // const Real_v yValues[fNoComponents],
             const Real_v magInitMomentumSq  //   Initial momentum square (used for rel. error)
      ) const
{
    // const  double fInvEpsilonRelSq = 1.0 / (eps_rel_max * eps_rel_max);
   
    Real_v epsPosition = fEpsRelMax * vecCore::math::Max(hStep, Real_v(fMinimumStep));
    // Note: it uses the remaining step 'h'
    //       Could change it to use full step size ==> move it outside loop !! 2017.11.10 JA

    Real_v invEpsPositionSq = 1.0 / (epsPosition * epsPosition);

    // Evaluate accuracy
    Real_v errpos_sq;
    errpos_sq = yEstError[0] * yEstError[0] + yEstError[1] * yEstError[1] + yEstError[2] * yEstError[2];
    errpos_sq *= invEpsPositionSq; // Scale relative to required tolerance
    
    // Accuracy for momentum
    Real_v invMagMomentumSq = 1.0 / (magInitMomentumSq + tinyValue);
    // Old choice: use  ( |momentum_final|^2 + tiny ) as divisor
    //   Real_v magmom_sq = yValues[3] * yValues[3] + yValues[4] * yValues[4] + yValues[5] * yValues[5];
    //   Real_v errmom_sq = sumerr_sq / (magmom_sq + tinyValue);
    // Oldest code: 
    //   vecCore::CondAssign(magmom_sq > 0.0, sumerr_sq/magmom_sq, sumerr_sq, &errmom_sq);    
    Real_v sumerr_sq = yEstError[3] * yEstError[3] + yEstError[4] * yEstError[4] + yEstError[5] * yEstError[5];
    Real_v errmom_sq = invMagMomentumSq * sumerr_sq ;

    errmom_sq *= fInvEpsilonRelSq;
    
    return vecCore::math::Max(errpos_sq, errmom_sq); // Maximum Square Error
}

#endif /* ErrorEstimatorSixVec_h */
