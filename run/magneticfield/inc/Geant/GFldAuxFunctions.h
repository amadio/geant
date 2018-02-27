//
//  Separate declarations of necessary or useful methods
//
//   Must be used only in implementation files (ie .cxx) - not in headers
// 
//  Author: J. Apostolakis
//

#ifndef GFldAuxFunctions_H
#define GFldAuxFunctions_H

// add the sincos function on MAC because sincos is not part of math.h
#ifdef __APPLE__ // possibly other conditions
inline void sincos(double x, double *s, double *c){
  __sincos(x,s,c);
}
#endif

#endif
