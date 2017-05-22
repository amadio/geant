#ifndef R123NameSpace_H
#define R123NameSpace_H 1

/**
  R123: NameSpace for SIMD/SIMT Random124 classes based on Random123 
 */

//namespace for Random-123 parameters

#include <array>
#include <stdint.h>
#include "VecRNG.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

namespace R123 {

  //typedef

  template<typename ReturnTypeBackendT> using array4_t = typename ReturnTypeBackendT::UInt32_v [4];
  template<typename ReturnTypeBackendT, unsigned int N> using array_t = typename ReturnTypeBackendT::UInt32_v [N];

  //Threefry parameters 

  VECPHYS_GLOBAL unsigned int R_64x4[8][2] = {
    {14,16},
    {52,57},
    {23,40},
    { 5,37},
    {25,33},
    {46,12},
    {58,22},
    {32,32}
  };  

  VECPHYS_GLOBAL unsigned int R_32x4[8][2] = {
    {10,26},
    {11,21},
    {13,27},
    {23, 5},
    { 6,20},
    {17,11},
    {25,10},
    {18,20}
  };  

  //Philox parameters
  /*
  constexpr uint32_t PHILOX_M4x32_0 = ((uint32_t)0xD2511F53);
  constexpr uint32_t PHILOX_M4x32_1 = ((uint32_t)0xCD9E8D57);
  constexpr uint32_t PHILOX_W32_0   = ((uint32_t)0x9E3779B9);
  constexpr uint32_t PHILOX_W32_1   = ((uint32_t)0xBB67AE85);
  */

  constexpr long long THREEFRY_SKIP_AHEAD = 2^64;
  constexpr long long PHILOX_SKIP_AHEAD = 2^64;

} // end of R123 namespace

} // end namespace impl
} // end namespace vecphys

#endif
