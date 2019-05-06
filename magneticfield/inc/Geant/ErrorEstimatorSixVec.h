//
//  ErrorEstimatorSixVec: Simple class to Estimate the overall error for each lane
//    
//
//  Created by japost on 12.12.2018
//

#define RETURN_ERRPOSMEM 1

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

    template <class Real_v>
      Real_v EstimateError( const Real_v   yError[fNoComponents],
                            const Real_v   hStep,
                            //nst Real_v   yValue[fNoComponents],
                            const Real_v   magMomentumSq, // Initial momentum square (used for rel. error)
                            Real_v       & epsPosition, 
                            Real_v       & errpos_sq,
                            Real_v       & errmom_sq                            
         ) const;
    //  Same as above, but returns intermediate values
   
    template <class Real_v>
      Real_v EstimateErrorWithFinalP( const Real_v yError[fNoComponents],
                                      const Real_v hStep,
                                      const Real_v yValue[fNoComponents] // Use it for momentum square (used for rel. error)
         ) const;
    // Same as above, but use final momentum value as divisor for momentum relative error
     
    double GetMaxRelativeError() const { return fEpsRelMax; }

  public:
    // static constexpr
    const double tinyValue = 1.0e-80; // Just to ensure there is no division by zero
    
  private:
    const double fEpsRelMax;
    const double fInvEpsilonRelSq; // = 1.0 / (eps_rel_max * eps_rel_max);
    const double fMinimumStep;
};

template <class Real_v> 
   Real_v ErrorEstimatorSixVec::EstimateError(
             const Real_v   yEstError[fNoComponents],
             const Real_v   hStep,
             const Real_v   magInitMomentumSq, // (Initial) momentum square (used for rel. error)
             Real_v       & epsPosition, 
             Real_v       & errpos_sq,
             Real_v       & errmom_sq
      ) const
{
   Real_v invMagMomentumSq = 1.0 / (magInitMomentumSq + tinyValue);
   
   epsPosition = fEpsRelMax * vecCore::math::Max(hStep, Real_v(fMinimumStep));
   // Note: it uses the remaining step 'h'
   //       Could change it to use full step size ==> move it outside loop !! 2017.11.10 JA

   Real_v invEpsPositionSq = 1.0 / (epsPosition * epsPosition);

   // Evaluate accuracy
   errpos_sq = yEstError[0] * yEstError[0] + yEstError[1] * yEstError[1] + yEstError[2] * yEstError[2];
   errpos_sq *= invEpsPositionSq; // Scale relative to required tolerance
   // Accuracy for momentum
   
   // Old choice: use  ( |momentum_final|^2 + tiny ) as divisor
   Real_v sumerr_sq = yEstError[3] * yEstError[3] + yEstError[4] * yEstError[4] + yEstError[5] * yEstError[5];
   errmom_sq = fInvEpsilonRelSq * invMagMomentumSq * sumerr_sq ;
   // Oldest code: 
   //   vecCore::CondAssign(magmom_sq > 0.0, sumerr_sq/magmom_sq, sumerr_sq, &errmom_sq);    
   
   // ReportOneLanePart2 ( epsPosition, errpos_sq, errmom_sq ); 
    
   return vecCore::math::Max(errpos_sq, errmom_sq); // Maximum Square Error
}

template <class Real_v> 
   Real_v ErrorEstimatorSixVec::EstimateError(
             const Real_v yEstError[fNoComponents],
             const Real_v hStep,
             const Real_v magInitMomentumSq  //   (Initial) momentum square (used for rel. error)
      ) const
{
   Real_v epsPosition=0.0, errpos_sq=0.0, errmom_sq= 0.0;
   return EstimateError( yEstError, hStep, magInitMomentumSq, epsPosition, errpos_sq, errmom_sq );
}

template <class Real_v> 
   Real_v ErrorEstimatorSixVec::EstimateErrorWithFinalP(
             const Real_v yEstError[fNoComponents],
             const Real_v hStep,
             const Real_v yValues[fNoComponents]
      ) const
{
    Real_v magmom_sq = yValues[3] * yValues[3] + yValues[4] * yValues[4] + yValues[5] * yValues[5];

    return EstimateErrorWithFinalP( yEstError, hStep, magmom_sq );
}
 
#endif /* ErrorEstimatorSixVec_h */
