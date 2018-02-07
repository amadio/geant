
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/AOS3D.h"
#include "base/Global.h"
// #include "backend/Backend.h"
#include <cassert>
#include <ctime>
#include "base/Vector.h"
#ifdef USE_ROOT
#include "TRandom1.h"
#else
#include "base/RNG.h"
using VECGEOM_NAMESPACE::RNG;
#endif

using namespace std;

int main(){
  
  vecgeom:: AOS3D<double> testAOS;
  // testAOS.push_back(2.,3.,10.);
  // testAOS.push_back(1.,2.,3.);
  //double c = testAOS[1][1];
  // vecgeom::Vector3D<double> vec = testAOS[1];
  // cout<<vec<<endl;
  //cout<<c<<endl;
  //cout<<testAOS.x(0)<<endl;

 // vecgeom::Vector3D<int> v1(1,2,3);
  vecgeom::Vector3D<double> v2(2.,3.,4.);
  vecgeom::Vector3D<float> v3(2.,3.,4.);
  vecgeom::Vector3D<double> diff = v2 - v3;
  cout<<"diff of v2: "<<v2<<" and v3: "<<v3<<" is: "<<diff<<endl;

/*  typedef vecgeom::Vector3D<double> ThreeVectorD;
  typedef vecgeom::Vector3D<float> ThreeVectorF;
  ThreeVectorD v5 = (ThreeVectorD) v3;
  ThreeVectorD v4 = v5 - v2 ;
  cout<<v4<<endl;*/

  // v1.push_back(10);
  //cout<<v1[0]<<endl;

  //Vc::SimdArray<double,4> a;

/*  vecgeom::kVc::precision_v   Br[2];
  Br[0] = 1;
  Br[1] = 2;
  //Br[2] = 3;
  cout<<Br[0]<<" "<<Br[1]<<" "<<Br[2]<<endl;

  //cout<<Br[10]<<endl;
  cout<<Br[1][3]<<endl;*/

#ifdef USE_ROOT
  // TRandom1::TRandom1  (   UInt_t seed,Int_t   lux = 3 )
  TRandom1 r1(0,4);
  TRandom1 *r2= new TRandom1(0,4);
  double x = r2->Rndm(1);
  double y = r2->Rndm(2);
  float z = r2->Rndm();
  cout<<r1.Rndm()<<endl;
  cout<<r1.Rndm()<<endl;
  cout<<x<<" "<<y<<" "<<z<<endl;
#else
  // TRandom1::TRandom1  (   UInt_t seed,Int_t   lux = 3 )
  RNG r1; // (0,4);
  RNG *r2 = new RNG; // (0,4);
  double x = r2->uniform(1);
  double y = r2->uniform(2);
  float z = r2->uniform();
  cout<<r1.uniform()<<endl;
  cout<<r1.uniform()<<endl;
  cout<<x<<" "<<y<<" "<<z<<endl;
#endif

  return 0;
}
