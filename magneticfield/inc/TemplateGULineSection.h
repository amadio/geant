//
// class TemplateGULineSection
//
// Class description:
//
// A utility class that calculates the distance of a point from a 
// line section.

// History:
// - Created. J. Apostolakis.
// --------------------------------------------------------------------

#ifndef TemplateGULineSection_hh
#define TemplateGULineSection_hh


#include <base/Vector3D.h> 


template <class Backend>
class TemplateGULineSection
{
  public:  // with description

    typedef vecgeom::Vector3D<typename Backend::precision_v>  ThreeVectorSimd; 
    typedef typename Backend::precision_v Double_v;

    inline TemplateGULineSection( const ThreeVectorSimd& PntA, 
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

template <class Backend>
inline
TemplateGULineSection<Backend>
  ::TemplateGULineSection( const ThreeVectorSimd& PntA, 
                           const ThreeVectorSimd& PntB )
  : EndpointA(PntA), VecAtoB(PntB-PntA)
{ 
  fABdistanceSq = VecAtoB.Mag2();  
}

template <class Backend>
inline
typename Backend::precision_v 
TemplateGULineSection<Backend>
  ::GetABdistanceSq() const
{
  return fABdistanceSq;
}

template <class Backend>
inline
typename Backend::precision_v 
TemplateGULineSection<Backend>
  ::Distline( const vecgeom::Vector3D<typename Backend::precision_v> & OtherPnt, 
              const vecgeom::Vector3D<typename Backend::precision_v> & LinePntA, 
              const vecgeom::Vector3D<typename Backend::precision_v> & LinePntB )
{
  TemplateGULineSection<Backend> LineAB( LinePntA, LinePntB );  // Line from A to B
  return LineAB.Dist( OtherPnt );
}

template <class Backend>
typename Backend::precision_v 
TemplateGULineSection<Backend>
  ::Dist( vecgeom::Vector3D<typename Backend::precision_v> OtherPnt ) const
{
  typename Backend::precision_v  dist_sq(0.);  
  vecgeom::Vector3D<typename Backend::precision_v>  VecAZ;
  typename Backend::precision_v sq_VecAZ, inner_prod, unit_projection(10.0) ; 

  VecAZ= OtherPnt - EndpointA;
  sq_VecAZ = VecAZ.Mag2();

  inner_prod= VecAtoB.Dot( VecAZ );
   
  //  Determine  Projection(AZ on AB) / Length(AB) 

  unit_projection = inner_prod/fABdistanceSq; //becomes nan when A and B are same point
  // vecgeom::MaskedAssign( fABdistanceSq != 0.0, inner_prod/fABdistanceSq, &unit_projection );
  vecgeom::MaskedAssign( (0. <= unit_projection ) && (unit_projection <= 1.0 ), sq_VecAZ - unit_projection*inner_prod, &dist_sq );
  vecgeom::MaskedAssign( unit_projection < 0.0, sq_VecAZ, &dist_sq);

  //below works because if fAbdist=0, then vecAtoB is zero.
  vecgeom::MaskedAssign( (fABdistanceSq == 0.0) || (unit_projection > 1.0), (OtherPnt -(EndpointA + VecAtoB)).Mag2(), &dist_sq); 

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

  // vecgeom::MaskedAssign( fABdistanceSq == 0.0, (OtherPnt - EndpointA).Mag2(), &dist_sq);
  // else
  // {
  //    dist_sq = (OtherPnt - EndpointA).Mag2() ;   
  // }  

  //Ananya: Can't see where dist_sq might be negative. Confirm. Can remove the maskedassign below in that case.
  vecgeom::MaskedAssign( dist_sq < 0.0, 0.0, &dist_sq );

  return vecgeom::VECGEOM_IMPL_NAMESPACE::Sqrt(dist_sq) ;  
}
#endif
