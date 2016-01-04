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


class GUVectorLineSection
{
  public:  // with description

    typedef vecgeom::Vector3D<typename vecgeom::kVc::precision_v>  ThreeVectorSimd; 
    typedef typename vecgeom::kVc::precision_v Double_v;

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
typename vecgeom::kVc::precision_v GUVectorLineSection::GetABdistanceSq() const
{
  return fABdistanceSq;
}

inline
typename vecgeom::kVc::precision_v GUVectorLineSection::Distline
              ( const vecgeom::Vector3D<typename vecgeom::kVc::precision_v> & OtherPnt, 
                const vecgeom::Vector3D<typename vecgeom::kVc::precision_v> & LinePntA, 
                const vecgeom::Vector3D<typename vecgeom::kVc::precision_v> & LinePntB )
{
  GUVectorLineSection LineAB( LinePntA, LinePntB );  // Line from A to B
  return LineAB.Dist( OtherPnt );
}


typename vecgeom::kVc::precision_v GUVectorLineSection::Dist
             ( vecgeom::Vector3D<typename vecgeom::kVc::precision_v> OtherPnt ) const
{
  typename vecgeom::kVc::precision_v  dist_sq;  
  vecgeom::Vector3D<typename vecgeom::kVc::precision_v>  VecAZ;
  typename vecgeom::kVc::precision_v sq_VecAZ, inner_prod, unit_projection(10.0) ; 

  VecAZ= OtherPnt - EndpointA;
  sq_VecAZ = VecAZ.Mag2();

  inner_prod= VecAtoB.Dot( VecAZ );
   
  //  Determine  Projection(AZ on AB) / Length(AB) 


  vecgeom::MaskedAssign( fABdistanceSq != 0.0, inner_prod/fABdistanceSq, &unit_projection );
  vecgeom::MaskedAssign( (0. <= unit_projection ) && (unit_projection <= 1.0 ), sq_VecAZ - unit_projection*inner_prod, &dist_sq );
  vecgeom::MaskedAssign( unit_projection < 0.0, sq_VecAZ, &dist_sq);
  vecgeom::MaskedAssign( (fABdistanceSq != 0.0) && (unit_projection > 1.0), (OtherPnt -(EndpointA + VecAtoB)).Mag2(), &dist_sq);
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

  vecgeom::MaskedAssign( !(fABdistanceSq != 0.0), (OtherPnt - EndpointA).Mag2(), &dist_sq);
  // else
  // {
  //    dist_sq = (OtherPnt - EndpointA).Mag2() ;   
  // }  

  vecgeom::MaskedAssign( dist_sq < 0.0, 0.0, &dist_sq );

  return Vc::sqrt(dist_sq) ;  
}
#endif
