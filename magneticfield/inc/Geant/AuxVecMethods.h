// Auxiliary methods, unrelated to printing

template <class Real_v>
Real_v PowerIf(const Real_v value, double exponent, vecCore::Mask_v<Real_v> condition = true); // const
  // Return the power in each 'lane':
  //  if( condition[i] ) { value[i]^exponent[i] } else { 1.0 }

// Definitions
// ===========---------------------------------------------------
template <class Real_v>
inline Real_v PowerIf(const Real_v value,
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
    result = Exp(exponent * Log(value));
  } else {
    // Do expensive 'pow' only for continuing ('condition') lanes
    for (size_t i = 0; i < vecCore::VectorSize<Real_v>(); ++i) {
      if (vecCore::Get(condition, i)) {
        double redFactor = Math::Pow(Get(value, i), exponent);
        vecCore::Set(result, i, redFactor);
      }
    }
  }

  return result;
}
