// Auxiliary methods, unrelated to printing

#ifndef AuxVecMethods_Def
#define AuxVecMethods_Def

template <class Real_v>
Real_v PowerSameIf(const Real_v inpValue, double exponent, vecCore::Mask_v<Real_v> condition = true); // const
  // Return the power in each 'lane':
  //  if( condition[i] ) { inpValue[i]^exponent[i] } else { 1.0 }

template <class Real_v>
Real_v PowerDiffIf(const Real_v inpValue, const Real_v exponent, vecCore::Mask_v<Real_v> condition = true);
// Same - but with varying exponent 
  // Return the power in each 'lane':
  //  if( condition[i] ) { inpValue[i]^exponent[i] } else { 1.0 }

// Definitions
// ===========---------------------------------------------------
template <class Real_v>
inline Real_v PowerSameIf(const Real_v inpValue,
                          double       exponent,
                          vecCore::Mask_v<Real_v> condition) // const
{
  using vecCore::Get;
  using vecCore::Set;
  using vecCore::math::Exp;
  using vecCore::math::Log;

  Real_v result(1.0);

  bool allNeeded = vecCore::MaskFull(condition);
  if (allNeeded) // All of the steps failed
  {
     // std::cout << " PowerSameIf: all Needed " << std::endl;
     result = Math::Pow( inpValue, Real_v(exponent) ); //
             // Exp(exponent * Log(inpValue));
  } else {
     // std::cout << " PowerSameIf: selection branch  " << std::endl;     
     // Do expensive 'pow' only for continuing ('condition') lanes
     for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
        if (vecCore::Get(condition, i)) {
           double resValue = Math::Pow(Get(inpValue, i), exponent);
           vecCore::Set(result, i, resValue);
        }
     }
  }

  return result;
}
// ===========---------------------------------------------------

template <class Real_v>
inline Real_v PowerDiffIf(const Real_v inpValue,
                          const Real_v exponent,
                          vecCore::Mask_v<Real_v> condition) // const
{
  using vecCore::Get;
  using vecCore::Set;
  using vecCore::math::Exp;
  using vecCore::math::Log;

  Real_v result(1.0);

  bool allNeeded = vecCore::MaskFull(condition);
  if (allNeeded) // All of the steps failed
  {
     // std::cout << " PowerDiffIf: all Needed " << std::endl;     
     result = Math::Pow( inpValue, exponent );
           // Exp(exponent * Log(inpValue));
  } else {
     // std::cout << " PowerDiffIf: selection branch  " << std::endl;          
    // Do expensive 'pow' only for needed ('condition') lanes
    for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
      if (vecCore::Get(condition, i)) {
        double resValue = std::pow( vecCore::Get(inpValue, i), vecCore::Get(exponent,i) );
           // std::exp( std::log( vecCore::Get(inpValue, i)) * Get(exponent, i));
            // Math::Pow(Get(inpValue, i), Get(exponent,i) );
        vecCore::Set(result, i, resValue);
      }
    }
  }

  return result;
}
// ===========---------------------------------------------------


template <class Real_v>
inline int countMaskTrue( vecCore::Mask_v<Real_v> flagLane )
{
   // using Bool_v = vecCore::Mask_v<Real_v>;

   constexpr unsigned int VecSize = vecCore::VectorSize<Real_v>();

   int     count=VecSize;
#if   0   
   if( ! vecCore::MaskFull(flagLane) ){
      count = 0;
      if( ! vecCore::MaskFull(flagLane) ){
         for( unsigned int i=0; i<VecSize; i++ ){
            if( vecCore::Get( flagLane, i ) )
               count++;
         }
      }
   }
#else   
   count = 0;
   for( unsigned int i=0; i<VecSize; i++ ){
      if( vecCore::Get( flagLane, i ) )
         count++;
   }   
#endif   
   return count;
}
#endif
