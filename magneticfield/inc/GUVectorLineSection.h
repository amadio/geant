//
// class GUVectorLineSection
//
// Class description:
//
// A utility class that calculates the distance of a point from a 
// line section.

// History:
// - Created. J. Apostolakis.
// --------------------------------------------------------------------

#ifndef GUVectorLineSection_hh
#define GUVectorLineSection_hh


#include <base/Vector3D.h> 
#include <Geant/VectorTypes.h>

class GUVectorLineSection
{
  using Double_v= Geant::Double_v;
   
  public:  // with description

    template<typename T>
       using Vector3D = vecgeom::Vector3D<T>;
    typedef Vector3D<Double_v>  ThreeVectorSimd; 

    inline GUVectorLineSection( const ThreeVectorSimd& PntA, 
                                const ThreeVectorSimd& PntB );

    Double_v Dist( ThreeVectorSimd OtherPnt ) const;

    inline Double_v GetABdistanceSq() const;

    inline static Double_v Distline( const ThreeVectorSimd& OtherPnt, 
                                     const ThreeVectorSimd& LinePntA, 
                                     const ThreeVectorSimd& LinePntB );
  private:

     ThreeVectorSimd   EndpointA;
     ThreeVectorSimd   VecAtoB;
     Double_v          fABdistanceSq;
};

// Inline methods implementations

inline
GUVectorLineSection::GUVectorLineSection( const ThreeVectorSimd& PntA, 
                                          const ThreeVectorSimd& PntB )
  : EndpointA(PntA), VecAtoB(PntB-PntA)
{ 
  fABdistanceSq = VecAtoB.Mag2();  
}

inline
Geant::Double_v GUVectorLineSection::GetABdistanceSq() const
{
  return fABdistanceSq;
}

inline
Geant::Double_v GUVectorLineSection::Distline
              ( const ThreeVectorSimd& OtherPnt, 
                const ThreeVectorSimd& LinePntA, 
                const ThreeVectorSimd& LinePntB )
{
  GUVectorLineSection LineAB( LinePntA, LinePntB );  // Line from A to B
  return LineAB.Dist( OtherPnt );
}


Geant::Double_v GUVectorLineSection::Dist
             ( ThreeVectorSimd OtherPnt ) const
{
  Double_v  dist_sq;  
  ThreeVectorSimd  VecAZ;
  Double_v sq_VecAZ, inner_prod, unit_projection(10.0) ; 

  VecAZ= OtherPnt - EndpointA;
  sq_VecAZ = VecAZ.Mag2();

  inner_prod= VecAtoB.Dot( VecAZ );
   
  //  Determine  Projection(AZ on AB) / Length(AB) 

  // veccore::MaskedAssign( fABdistanceSq != 0.0, inner_prod/fABdistanceSq, &unit_projection );
  vecCore::MaskedAssign( unit_projection, fABdistanceSq != 0.0, inner_prod/fABdistanceSq );
  
  // vecCore::MaskedAssign( (0. <= unit_projection ) && (unit_projection <= 1.0 ), sq_VecAZ - unit_projection*inner_prod, &dist_sq );
  // Mask<Geant::Double_v> goodProjection = (0. <= unit_projection ) && (unit_projection <= 1.0 );
  vecCore::MaskedAssign( dist_sq,
                         // goodProjection,
                         (0. <= unit_projection ) && (unit_projection <= 1.0 ),                       
                         sq_VecAZ - unit_projection*inner_prod );
  
  // -- vecCore::MaskedAssign( unit_projection < 0.0, sq_VecAZ, &dist_sq);
  // Mask<Geant::Double_v> negativeProj = unit_projection < 0.0 ;  
  // vecCore::MaskedAssign( dist_sq, negativeProj, sq_VecAZ );
  vecCore::MaskedAssign( dist_sq, unit_projection < 0.0, sq_VecAZ );
  
  
  // vecCore::MaskedAssign( (fABdistanceSq != 0.0) && (unit_projection > 1.0), (OtherPnt -(EndpointA + VecAtoB)).Mag2(), &dist_sq);
  // Mask<Geant::Double_v> condDistProj=  (fABdistanceSq != 0.0) && (unit_projection > 1.0),;
  vecCore::MaskedAssign( dist_sq, 
                         // condDistProj,
                         (fABdistanceSq != 0.0) && (unit_projection > 1.0),                         
                         (OtherPnt -(EndpointA + VecAtoB)).Mag2());
  
  // if( fABdistanceSq != 0.0 )
  // {
  //   unit_projection = inner_prod/fABdistanceSq;

  //   if( (0. <= unit_projection ) && (unit_projection <= 1.0 ) )
  //   {
  //     dist_sq= sq_VecAZ -  unit_projection * inner_prod ;
  //   }
  //   else
  //   {
     
  //     if( unit_projection < 0. )
  //     {
  //       dist_sq= sq_VecAZ;  
  //     }
  //     else                       
  //     {
  //       ThreeVectorSimd   EndpointB = EndpointA + VecAtoB;
  //       ThreeVectorSimd   VecBZ =     OtherPnt - EndpointB;
  //       dist_sq =  VecBZ.Mag2();
  //     }
  //   }
  // }

  vecCore::MaskedAssign( dist_sq,  !(fABdistanceSq != 0.0), (OtherPnt - EndpointA).Mag2() );
  // else
  // {
  //    dist_sq = (OtherPnt - EndpointA).Mag2() ;   
  // }  

  // vecgeom::MaskedAssign( dist_sq < 0.0, 0.0, &dist_sq );  
  vecCore::MaskedAssign( dist_sq, dist_sq < Double_v(0.0), Double_v(0.0) );
  // dist_sq = vecgeom::Max ( dist_sq, 0.0 );
  
  return vecCore::math::Sqrt(dist_sq) ;
}
#endif
