#include "StepChecker.h"

// For Geant::Print,  Printf
#include "Geant/Error.h"
#include "ConstFieldHelixStepper.h"
#include "ConstBzFieldHelixStepper.h"

VECCORE_ATT_HOST_DEVICE          
bool
StepChecker::
 CompareStep( vecgeom::Vector3D<double> const & Position,  
              vecgeom::Vector3D<double> const & Direction,
                                 double         charge,              
                                 double         momentum,
                                 double         step,                            
              vecgeom::Vector3D<double> const & endPosition,
              vecgeom::Vector3D<double> const & endDirection,
              vecgeom::Vector3D<double> const & endPositionRef, //  Reference solution / position
              vecgeom::Vector3D<double> const & endDirectionRef //  Reference solution / direction
      ) const
{
  using ThreeVector = vecgeom::Vector3D<double>;

  if( fVerbose )
  {
     printf("Particle with charge %g  momentum = %8.4g Position= %10.6f %10.6f %10.6f  Momentum= %10.6f %10.6f %10.6f \n",
            charge, momentum, Position[0], Position[1], Position[2], 
            Direction[0],  Direction[1], Direction[2] );
     /* Printf("                  End    Position= %10.6f  %10.6f  %10.6f  Momentum= %10.6f %10.6f %10.6f \n",
         PositionNew[0], PositionNew[1], PositionNew[2],
         DirectionNew[0],  DirectionNew[1], DirectionNew[2] ); */
     ThreeVector PositionShift=  endPosition  - Position;
     ThreeVector DirectionShift= endDirection - Direction;  
     printf(" Step=%8f  Move/RK4:  Pos/shift= %10.6f  %10.6f  %10.6f  Mom/shft= %10.3g %10.3g %10.3g \n",
            step,
            PositionShift[0],  PositionShift[1],  PositionShift[2],
            DirectionShift[0],  DirectionShift[1],  DirectionShift[2]);         
  }

  const double epsTrigger = fMaxRelDiff;
  const double triggerLen = epsTrigger * step; 
  const double triggerDir = epsTrigger;

  //************
  ThreeVector PositionDiff= endPositionRef   - endPosition;
  ThreeVector DirectionDiff= endDirectionRef - endDirection;
  bool differ = PositionDiff.Mag() > triggerLen ||
     DirectionDiff.Mag() > triggerDir;
  if( differ )
  {  
     printf(" Step=%8f    DIFFERS   Pos/diff= %10.6f  %10.6f  %10.6f  Mom/diff= %10.3g %10.3g %10.3g \n",
            step,
            PositionDiff[0],  PositionDiff[1],  PositionDiff[2],
            DirectionDiff[0],  DirectionDiff[1],  DirectionDiff[2]);
     
     ThreeVector PositionShiftH=  endPositionRef   - Position;
     ThreeVector DirectionShiftH= endDirectionRef  - Direction;  
     printf(" Step=%8f  Helix/Mv:  Pos/shift= %10.6f  %10.6f  %10.6f  Mom/shft= %10.3g %10.3g %10.3g \n",
            step,
            PositionShiftH[0],  PositionShiftH[1],  PositionShiftH[2],
            DirectionShiftH[0],  DirectionShiftH[1],  DirectionShiftH[2]);
  }
  return differ;
}


   /** @brief Check against helical solution */
bool 
StepChecker::CheckStep( vecgeom::Vector3D<double> const & Position,
                        vecgeom::Vector3D<double> const & Direction,
                                           double         charge,
                                           double         momentum,
                                           double         step,
                        vecgeom::Vector3D<double> const & endPosition,
                        vecgeom::Vector3D<double> const & endDirection,
                        vecgeom::Vector3D<double> const & BfieldVec ) const
{
  using ThreeVector = vecgeom::Vector3D<double>;
  
  ThreeVector PositionNewHelix(0.,0.,0.);
  ThreeVector DirectionNewHelix(0.,0.,0.);
#if 1
  // Simpler version
  Geant::ConstFieldHelixStepper stepper(BfieldVec[0], BfieldVec[1], BfieldVec[2]);
  stepper.DoStep<double>(Position, Direction, charge,
                         momentum, step,
                         PositionNewHelix, DirectionNewHelix);
#else
  // More complicated version ...
  if ( std::fabs( BfieldVec[2] ) > 1e6 * std::max( std::fabs(BfieldVec[0]), std::fabs(BfieldVec[1]) ) ) {
     // Old - constant field in Z-direction
     Geant::ConstBzFieldHelixStepper stepper( BfieldVec[2] ); // z-component
     stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  charge, momentum, step,
                                            PositionNewHelix, DirectionNewHelix);
  } else {
     Geant::ConstFieldHelixStepper stepper( Bfield ); // double Bfield[3] );
     stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  charge, momentum, step,
                                            PositionNewHelix, DirectionNewHelix);
  }
#endif
  
  return CompareStep( Position,         Direction,  charge,  momentum,  step,
                      endPosition,   endDirection,
                      PositionNewHelix, DirectionNewHelix);               
}
