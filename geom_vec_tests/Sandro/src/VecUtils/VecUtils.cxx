#include "VecUtils.h"


#ifdef HAVE_SSE2
#  include <emmintrin.h>
/** 1.0 double for sse2 */
  __m128d sse2_double_one=_mm_set_pd(1.,1.);
#elif HAVE_AVX
#include <immintrin.h>
//#include <emmintrin.h>
  __m256d avx_double_one=_mm256_set_pd(1.,1.,1.,1.);
#endif

#ifdef HAVE_SSE2
union w128_t
{
  __m128d w;
  double d[2];
};
#elif HAVE_AVX
union w256_t
{
  __m256d w;
  double d[4];
};
#endif


void VecUtils::min_v( const  double * x, const  double * y, double * result, const int n )
{
  // most naive implementation using intrinsics, no check for alignement etc...
#ifdef HAVE_AVX
  // truly vectorized version (intrinsics)?
  w256_t xl,yl,r;

  // load the registers and perform min operation
  for(unsigned int i=0;i<n;i+=4)
    {
      xl.w = _mm256_set_pd(x[i],x[i+1],x[i+2],x[i+3]);
      yl.w = _mm256_set_pd(y[i],y[i+1],y[i+2],y[i+3]);
      r.w=_mm256_min_pd (xl.w, yl.w);
      result[i]=r.d[0];
      result[i+1]=r.d[1];
      result[i+2]=r.d[2];
      result[i+3]=r.d[3];
    }
#endif
#ifdef HAVE_SSE2
  // truly vectorized version (intrinsics)?
  w128_t xl,yl,r;

  // load the registers and perform min operation
  for(unsigned int i=0;i<n;i+=2)
    {
      xl.w = _mm_set_pd(x[i],x[i+1]);
      yl.w = _mm_set_pd(y[i],y[i+1]);
      r.w=_mm_min_pd (xl.w, yl.w);
      result[i]=r.d[0];
      result[i+1]=r.d[1];
    }
#endif
}


double VecUtils::min(double x, double y)
{
  double d;
  d=x;
  if (y<x) d=y;
  return d;
  //  return (x>y)? x : y; 
}


// loop version of min function
void VecUtils::min_l( const  double * x, const  double * y, double * result, const int n )
{
  for(unsigned int i=0;i<n;i++)
    {
      result[i]=(x[i] > y[i])? y[i] : x[i];
    }
}


void VecUtils::min_Vc( const  double *, const  double *, double * result, const int n )
{
// version that dispatches to Vc library 

}



void VecUtils::abs_v( const double * x, double * y, int n)
{
  const unsigned long long signmask = 1UL << 63;
  for(int unsigned i=0; i<n; i++)
    {
      y[i] =  ( 2.* ( (( (unsigned long long) x[i] ) & signmask ) - 1./2 ) ) * x[i];
    }
}



double VecUtils::abs1(double x)
{
  return (x>0)? x: -x;
}



double VecUtils::abs2(double x)
{
  unsigned char s=127;
  union
  {
    unsigned char c[7];
    double x;
  } data;
  data.x=x;
  data.c[7] = data.c[7] & s;
  return data.x;
}


