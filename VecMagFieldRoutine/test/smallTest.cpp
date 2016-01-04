//===--- BenchmarkTiming.cpp - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file  BenchmarkTiming.cpp
 * @brief Benchmark of implementations of bi-linear interpolation of CMS field
 * @author Ananya
 */
//===----------------------------------------------------------------------===//

#include <iostream>

#include <string>
#include <vector>
#include <ctime>
#include <cmath> //for sqrt
// #include <stdlib.h>
#include <cstdlib>

#include <numeric>
#include <string>
#include <functional>

#include <Vc/Vc>

#include "backend/vc/Backend.h"
// #include "backend/vcfloat/Backend.h"
#include "VcFloatBackend.h"
#include "base/Vector.h"

#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"

// #include "MagField.h"

#include "TemplateCMSmagField.h"

using namespace std;

typedef vecgeom::Vector3D<double> ThreeVector; //normal Vector3D
typedef vecgeom::Vector3D<vecgeom::kVc::precision_v> ThreeVecSimd_t;
typedef vecgeom::Vector<double> VcVectorFloat;




int main(){

  TemplateCMSmagField<vecgeom::kVc> m1;
  TemplateCMSmagField<vecgeom::kScalar> m2;

  m1.ReadVectorData("../VecMagFieldRoutine/cmsmagfield2015.txt");
  m2.ReadVectorData("../VecMagFieldRoutine/cmsmagfield2015.txt");
  ThreeVector position, xyzField;
  ThreeVecSimd_t position_v, xyzField_v;

  std::cout<<"Give x,y and z"<<std::endl;
  std::cin>>position.x()>>position.y()>>position.z();

  for (int i = 0; i < 3; ++i)
  {
  	position_v[i] = position[i];
  	xyzField_v[i] = xyzField[i];
  }

  m2.GetFieldValue(position, xyzField);
  m1.GetFieldValue(position_v, xyzField_v);


  std::cout<<"Magnetic Field at "<<position<<" is: "<<xyzField<<" tesla."<<std::endl;
  cout<<"Vector mag field is: " << xyzField_v << endl;

  for (int i = 0; i < 4; ++i)
  {
  	position[0] = i;
  	position[1] = i*10;
  	position[2] = i*100;
  	position_v[0][i] = i;
  	position_v[1][i] = i*10;
  	position_v[2][i] = i*100;

  	m2.GetFieldValue(position, xyzField);
	  cout << "Magnetic Field at " << position <<" is: "<< xyzField <<" tesla."<< endl;
  }
  m1.GetFieldValue(position_v, xyzField_v);
  cout<<"Vector mag field at position: "<< position_v  << " is: " << xyzField_v << endl;


  return 0;
}


