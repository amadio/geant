// ApproxEqual Functions for geometry test programs
//
// History:
// 20.07.95 P.Kent Translated from old code

#ifndef APPROXEQUAL_HH
#define APPROXEQUAL_HH

#include <iomanip>

const float kApproxEqualTolerance = 1E-6;

// Return true if the double x is approximately equal to y
//
// Process:
//
// Return true is x if less than kApproxEqualTolerance from y

bool ApproxEqual(const float x, const float y, const float r, const float z, const int i)
{
  using std::cout;
  using std::setw;
  using std::endl;
  
  if (x == y) {
    // std::cout<<"case1"<<std::endl;
    return true;
  } else if ( false ) { // x * y == 0.0) {
    // std::cout<<"case2"<<std::endl;
    float diff = std::fabs(x - y);
    cout << "Diff : " << diff << std::endl;
    return diff < kApproxEqualTolerance;
    return true;
  } else {
    // cout<<"case3"<<std::endl;
    float diff  = std::fabs(x - y);
    float abs_x = std::fabs(x), abs_y = std::fabs(y);
    if ( diff > kApproxEqualTolerance * (abs_x + abs_y) ) {
       cout << "Difference at  r: " << setw(6) << r << " and z: " << setw(6) << z << " : ";
       if (i == 1) cout << "Midpoint btwn r values";
       if (i == 2) cout << "Midpoint btwn z values";
       if (i == 3) cout << "Middle of cell        ";
       cout << "  val1 = " << setw( 14 ) << x << " val2= " << setw(14) << y
            << " Relative diff = " << diff / (abs_x + abs_y)
            << std::endl;
      // cout <<"\n"<<std::endl;
    }
    // return true;
    return diff < kApproxEqualTolerance * (abs_x + abs_y) ;
  }
}

// Return true if the 3vector check is approximately equal to target
template <class Vec_t>
bool ApproxEqual(const Vec_t &check, const Vec_t &target, const float r, const float z, const int i)
{
  return (ApproxEqual(check.x(), target.x(), r, z, i) && ApproxEqual(check.y(), target.y(), r, z, i) &&
          ApproxEqual(check.z(), target.z(), r, z, i))
             ? true
             : false;
}

bool ApproxEqual(const double x, const double y)
{
  if (x == y) {
    return true;
  } else if (x * y == 0.0) {
    double diff = std::fabs(x - y);
    return diff < kApproxEqualTolerance;
  } else {
    double diff  = std::fabs(x - y);
    double abs_x = std::fabs(x), abs_y = std::fabs(y);
    return diff / (abs_x + abs_y) < kApproxEqualTolerance;
  }
}

bool ApproxEqual(const float x, const float y)
{
  if (x == y) {
    return true;
  } else if (x * y == 0.0) {
    float diff = std::fabs(x - y);
    return diff < kApproxEqualTolerance;
  } else {
    float diff  = std::fabs(x - y);
    float abs_x = std::fabs(x), abs_y = std::fabs(y);
    return diff / (abs_x + abs_y) < kApproxEqualTolerance;
  }
}

// Return true if the 3vector check is approximately equal to target
template <class Vec_t>
bool ApproxEqual(const Vec_t &check, const Vec_t &target)
{
  return (ApproxEqual(check.x(), target.x()) && ApproxEqual(check.y(), target.y()) &&
          ApproxEqual(check.z(), target.z()))
             ? true
             : false;
}

#endif
