// ApproxEqual Functions for geometry test programs
//
// History:
// 20.07.95 P.Kent Translated from old code

#ifndef APPROXEQUAL_HH
#define APPROXEQUAL_HH

const float kApproxEqualTolerance = 1E-6;

// Return true if the double x is approximately equal to y
//
// Process:
//
// Return true is x if less than kApproxEqualTolerance from y

bool ApproxEqual(const float x, const float y, const float r, const float z, const int i)
{
    if (x == y) {
        //std::cout<<"case1"<<std::endl;
        return true;
    } else if (x * y == 0.0) {
        //std::cout<<"case2"<<std::endl;
        float diff = std::fabs(x - y);
        std::cout<<"Diff : "<<diff<<std::endl;
        return diff < kApproxEqualTolerance;
        return true;
    } else {
        //std::cout<<"case3"<<std::endl;
        float diff = std::fabs(x - y);
        float abs_x = std::fabs(x), abs_y = std::fabs(y);
        if((diff / (abs_x + abs_y)) > kApproxEqualTolerance){
            std::cout<<"\nFor r: "<<r<<" and z: "<<z<<"\nRelative error is:  "<<diff / (abs_x + abs_y)*1e+9<<std::endl;
            if(i==1) std::cout<<"On edge, between 2 r values"<<std::endl;
            if(i==2) std::cout<<"On edge, between 2 z values"<<std::endl;
            if(i==3) std::cout<<"Middle of cell"<<std::endl;
            //std::cout<<"\n"<<std::endl;
        }
        return true;
        return diff / (abs_x + abs_y) < kApproxEqualTolerance;
    }
}

// Return true if the 3vector check is approximately equal to target
template <class Vec_t>
bool ApproxEqual(const Vec_t& check, const Vec_t& target, const float r, const float z, const int i)
{
  return (ApproxEqual(check.x(),target.x(), r, z, i)&&
      ApproxEqual(check.y(),target.y(), r, z, i)&&
      ApproxEqual(check.z(),target.z(), r, z, i)) ? true : false;
}


bool ApproxEqual(const double x, const double y)
{
    if (x == y) {
        return true;
    } else if (x * y == 0.0) {
        double diff = std::fabs(x - y);
        return diff < kApproxEqualTolerance;
    } else {
        double diff = std::fabs(x - y);
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
        float diff = std::fabs(x - y);
        float abs_x = std::fabs(x), abs_y = std::fabs(y);
        return diff / (abs_x + abs_y) < kApproxEqualTolerance;
    }
}

// Return true if the 3vector check is approximately equal to target
template <class Vec_t>
bool ApproxEqual(const Vec_t& check, const Vec_t& target)
{
  return (ApproxEqual(check.x(),target.x())&&
      ApproxEqual(check.y(),target.y())&&
      ApproxEqual(check.z(),target.z())) ? true : false;
}

#endif









