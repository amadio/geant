/*
 * log.h
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

#ifndef GEANTV_FASTLOG_H
#define GEANTV_FASTLOG_H

#include <limits>
#include "FastMathCommon.h"

namespace geant {

// local namespace for the constants/functions which are necessary only here
namespace details {

const double LOG_UPPER_LIMIT = 1e307;
const double LOG_LOWER_LIMIT = 0;

const double SQRTH = 0.70710678118654752440;

inline double get_log_px(const double x)
{
  const double PX1log = 1.01875663804580931796E-4;
  const double PX2log = 4.97494994976747001425E-1;
  const double PX3log = 4.70579119878881725854E0;
  const double PX4log = 1.44989225341610930846E1;
  const double PX5log = 1.79368678507819816313E1;
  const double PX6log = 7.70838733755885391666E0;

  double px = PX1log;
  px *= x;
  px += PX2log;
  px *= x;
  px += PX3log;
  px *= x;
  px += PX4log;
  px *= x;
  px += PX5log;
  px *= x;
  px += PX6log;
  return px;
}

inline double get_log_qx(const double x)
{
  const double QX1log = 1.12873587189167450590E1;
  const double QX2log = 4.52279145837532221105E1;
  const double QX3log = 8.29875266912776603211E1;
  const double QX4log = 7.11544750618563894466E1;
  const double QX5log = 2.31251620126765340583E1;

  double qx = x;
  qx += QX1log;
  qx *= x;
  qx += QX2log;
  qx *= x;
  qx += QX3log;
  qx *= x;
  qx += QX4log;
  qx *= x;
  qx += QX5log;
  return qx;
}
}

// Log double precision --------------------------------------------------------
inline double Log(double x)
{

  const double original_x = x;

  /* separate mantissa from exponent */
  double fe;
  x = details::getMantExponent(x, fe);

  // blending
  x > details::SQRTH ? fe += 1. : x += x;
  x -= 1.0;

  /* rational form */
  double px = details::get_log_px(x);

  // for the final formula
  const double x2 = x * x;
  px *= x;
  px *= x2;

  const double qx = details::get_log_qx(x);

  double res = px / qx;

  res -= fe * 2.121944400546905827679e-4;
  res -= 0.5 * x2;

  res = x + res;
  res += fe * 0.693359375;

  if (original_x > details::LOG_UPPER_LIMIT) res = std::numeric_limits<double>::infinity();
  if (original_x < details::LOG_LOWER_LIMIT) // THIS IS NAN!
    res = -std::numeric_limits<double>::quiet_NaN();

  return res;
}

inline double Log10(double x)
{
  const double invlog10 = 0.434294481903251827651128;
  return Log(x) * invlog10;
}
}

#endif // GEANTV_FASTLOG_H
