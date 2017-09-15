//
// class GULineSection
//
// Class description:
//
// A utility class that calculates the distance of a point from a 
// line section.

// History:
// - Created. J. Apostolakis.
// --------------------------------------------------------------------

#ifndef GULineSection_hh
#define GULineSection_hh

// #include "G4Types.hh"
// #include "ThreeVector.h"
#include <base/Vector3D.h> 

class GULineSection
{
   using ThreeVector = vecgeom::Vector3D<double>;
   
  public:  // with description

     inline GULineSection( const ThreeVector& PntA, const ThreeVector& PntB );

     double Dist( ThreeVector OtherPnt ) const;

     inline double GetABdistanceSq() const;

     inline static double   Distline( const ThreeVector& OtherPnt, 
                                      const ThreeVector& LinePntA, 
                                      const ThreeVector& LinePntB );
  private:

     ThreeVector   EndpointA;
     ThreeVector   VecAtoB;
     double fABdistanceSq ;
};

// Inline methods implementations

inline
GULineSection::GULineSection( const ThreeVector& PntA, 
                              const ThreeVector& PntB )
  : EndpointA(PntA), VecAtoB(PntB-PntA)
{ 
  fABdistanceSq = VecAtoB.Mag2();  
}

inline
double GULineSection::GetABdistanceSq() const
{
  return fABdistanceSq;
}

inline
double GULineSection::Distline( const ThreeVector& OtherPnt, 
                                const ThreeVector& LinePntA, 
                                const ThreeVector& LinePntB )
{
  GULineSection LineAB( LinePntA, LinePntB );  // Line from A to B
  return LineAB.Dist( OtherPnt );
}


double GULineSection::Dist( ThreeVector OtherPnt ) const
{
  double       dist_sq;  
  ThreeVector  VecAZ;
  double sq_VecAZ, inner_prod, unit_projection ; 

  VecAZ= OtherPnt - EndpointA;
  sq_VecAZ = VecAZ.Mag2();

  inner_prod= VecAtoB.Dot( VecAZ );
   
  //  Determine  Projection(AZ on AB) / Length(AB) 
  //
  if( fABdistanceSq != 0.0 )
  {
    //  unit_projection= inner_prod * InvsqDistAB();
    unit_projection = inner_prod/fABdistanceSq;

    if( (0. <= unit_projection ) && (unit_projection <= 1.0 ) )
    {
      dist_sq= sq_VecAZ -  unit_projection * inner_prod ;
    }
    else
    {
     //  The perpendicular from the point to the line AB meets the line
     //   in a point outside the line segment!
     
      if( unit_projection < 0. ) // A is the closest point
      {
        dist_sq= sq_VecAZ;  
      }
      else                       // B is the closest point
      {
        ThreeVector   EndpointB = EndpointA + VecAtoB;
        ThreeVector   VecBZ =     OtherPnt - EndpointB;
        dist_sq =  VecBZ.Mag2();
      }
    }
  }
  else
  {
     dist_sq = (OtherPnt - EndpointA).Mag2() ;   
  }  
  if( dist_sq < 0.0 ) dist_sq = 0.0 ;

  return std::sqrt(dist_sq) ;  
}
#endif
