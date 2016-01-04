#ifndef FIELD_CONSTANTS_H
#define FIELD_CONSTANTS_H 1

#include <cfloat>
#include "Units.h"

namespace Constants {

static constexpr double     pi = 3.1415926535897932384626;
static constexpr double  twopi = 2.0 * pi;

static constexpr double c_light   = 2.99792458e+8 * fieldUnits::meter/fieldUnits::second;
static constexpr double c_squared = c_light * c_light;

};
#endif
