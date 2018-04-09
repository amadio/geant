#ifndef GEANTV_FASTEXP_H
#define GEANTV_FASTEXP_H

/*
 * exp.h
 * The basic idea is to exploit Pade polynomials.
 * A lot of ideas were inspired by the cephes math library (by Stephen L. Moshier
 * moshier@na-net.ornl.gov) as well as actual code.
 * The Cephes library can be found here:  http://www.netlib.org/cephes/
 *
 *  Created on: Jun 23, 2012
 *      Author: Danilo Piparo, Thomas Hauth, Vincenzo Innocente
 */

/*
 * VDT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.
 *
 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <limits>
#include "FastMathCommon.h"

namespace geant {

namespace details {

const double EXP_LIMIT = 708;

const double PX1exp = 1.26177193074810590878E-4;
const double PX2exp = 3.02994407707441961300E-2;
const double PX3exp = 9.99999999999999999910E-1;
const double QX1exp = 3.00198505138664455042E-6;
const double QX2exp = 2.52448340349684104192E-3;
const double QX3exp = 2.27265548208155028766E-1;
const double QX4exp = 2.00000000000000000009E0;

const double LOG2E = 1.4426950408889634073599; // 1/log(2)
}

// Exp double precision --------------------------------------------------------

/// Exponential Function double precision
inline double Exp(double initial_x)
{

  double x  = initial_x;
  double px = details::fpfloor(details::LOG2E * x + 0.5);

  const int32_t n = int32_t(px);

  x -= px * 6.93145751953125E-1;
  x -= px * 1.42860682030941723212E-6;

  const double xx = x * x;

  // px = x * P(x**2).
  px = details::PX1exp;
  px *= xx;
  px += details::PX2exp;
  px *= xx;
  px += details::PX3exp;
  px *= x;

  // Evaluate Q(x**2).
  double qx = details::QX1exp;
  qx *= xx;
  qx += details::QX2exp;
  qx *= xx;
  qx += details::QX3exp;
  qx *= xx;
  qx += details::QX4exp;

  // e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
  x = px / (qx - px);
  x = 1.0 + 2.0 * x;

  // Build 2^n in double.
  x *= details::uint642dp((((uint64_t)n) + 1023) << 52);

  if (initial_x > details::EXP_LIMIT) x  = std::numeric_limits<double>::infinity();
  if (initial_x < -details::EXP_LIMIT) x = 0.;

  return x;
}

} // end namespace geant

#endif // GEANTV_FASTEXP_H
