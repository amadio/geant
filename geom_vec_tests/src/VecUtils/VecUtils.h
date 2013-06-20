//
// utility class / functions to test and provide vectorization of often occuring math operations
// started: 2.5.2013; Sandro Wenzel
//

#ifndef VecUtils_H
#define VecUtils_H

class VecUtils
{
 public:

  static double min(double x, double y);
  // normal version

  static void min_v( const  double *, const  double *, double * result, const int n );
  // truly vectorized version (intrinsics)?

  static void min_l( const  double *, const  double *, double * result, const int n );
  // loop version of min function

  static void min_Vc( const double *, const  double *, double * result, const int n );
  // version that dispatches to Vc library 

  // ------------------------------------------------------------------------------------
  static void abs_v( const double *, double * result, const int n);

  static double abs2( double );
  static double abs1( double );

};


#endif // VecUtils_H
