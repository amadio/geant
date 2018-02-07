#include "base/Vector3D.h"
#include "base/Global.h"
#include "base/SOA3D.h"
#include <iostream>
#include <Vc/Vc>
#include "backend/vc/Backend.h"
#include "backend/vcfloat/Backend.h"
#include "backend/scalarfloat/Backend.h"
#include "backend/Backend.h"


using namespace vecgeom;
int main(){
Vector3D<double> vec(1,2,3);
std::cout<<"Vec3d vector  : "<<vec<<std::endl;

SOA3D<double> soa(2);
int b = soa.size();
std::cout<<"SOA size is: "<< b <<std::endl;


typedef vecgeom::Vector3D<float> Vec_t;
typedef vecgeom::Vector3D<vecgeom::kVcFloat::precision_v> ThreeVectorVcFloat;
typedef vecgeom::Vector3D<vecgeom::kVc::precision_v> ThreeVectorVc;

ThreeVectorVc v1(1.,1.,1.), v3;

float f1, f2, f3 =4;
ThreeVectorVcFloat v2(f1,f2,f3);
std::cout<<"ThreeVectorVc uninitialized is: "<<v1<<std::endl;
std::cout<<"ThreeVectorVcFloat uninitialized is: "<<v2<<std::endl;


return 0;
}
