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

#include "base/Vector3D.h"
#include <Geant/VectorTypes.h>

// #include "VC_NO_MEMBER_GATHER"

#include "CMSmagField.h"
// using MagField=MagField3<float>

using namespace std;

using Double_v = Geant::Double_v;
using Float_v = Geant::Float_v;

typedef vecgeom::Vector3D<double> ThreeVector; // normal Vector3D
typedef vecgeom::Vector3D<Double_v> ThreeVecSimd_t;
typedef vecgeom::Vector3D<Float_v> ThreeVecSimdF_t;

// typedef MagVector3<float>         MagField;
// using MagField=MagVector3<float>

using MagField = CMSmagField;

// constexpr float tesla = geant::tesla;
// constexpr float kilogauss = geant::kilogauss;
constexpr float millimeter = geant::millimeter;

const float kRMax = 9000 * millimeter;
const float kZMax = 16000 * millimeter;

float RandR()
{
  float rnd = (float)rand() / (RAND_MAX);
  return rnd * kRMax;
}

float RandZ()
{
  float rnd = (float)rand() / (RAND_MAX);
  return -kZMax + 2 * rnd * kZMax;
}

void GenVecCart(ThreeVector &pos)
{
  float rnd = (float)rand() / (RAND_MAX);
  float phi = 2. * M_PI * rnd;
  float r   = RandR();
  float x   = r * vecCore::math::Cos(phi);
  float y   = r * vecCore::math::Sin(phi);
  float z   = RandZ();
  pos.Set(x, y, z);
}

void GenVecCart(vector<ThreeVector> &posVec, const size_t &n)
{
  for (size_t i = 0; i < n; ++i) {
    ThreeVector pos;
    GenVecCart(pos);
    posVec.push_back(pos);
  }
}

float TimeScalar(MagField &m1, const vector<ThreeVector> &posVec, vector<ThreeVector> &outputVec, const size_t &n)
{
  ThreeVector xyzField, sumField(0.);
  vector<float> scaTimePerRepitition;

  size_t noRunsAvg = 16;

  cout << "=== Scalar field start: " << endl;
  float tmin  = FLT_MAX;
  float tmax  = -FLT_MAX;
  size_t imin = 0, imax = 0;

  for (size_t k = 0; k < noRunsAvg; k++) {
    clock_t clock1 = clock();
    for (size_t i = 0; i < n; ++i) {
      m1.GetFieldValue(posVec[i], xyzField);
      // std::cout << i << ": " << posVec[i] << " => " << xyzField << std::endl;
      sumField += xyzField;
      outputVec.push_back(xyzField);
    }
    clock1              = clock() - clock1;
    float clock1InFloat = ((float)clock1) / CLOCKS_PER_SEC;
    float tpercall      = clock1InFloat / n;
    if (tpercall > tmax) {
      tmax = tpercall;
      imax = k;
    }
    if (tpercall < tmin) {
      tmin = tpercall;
      imin = k;
    }
    scaTimePerRepitition.push_back(tpercall);
  }
  // Remove imin and imax measurements from the sample
  scaTimePerRepitition.erase(scaTimePerRepitition.begin() + std::max(imin, imax));
  scaTimePerRepitition.erase(scaTimePerRepitition.begin() + std::min(imin, imax));

  float timeSum  = std::accumulate(scaTimePerRepitition.begin(), scaTimePerRepitition.end(), 0.0);
  float timeMean = timeSum / scaTimePerRepitition.size();
  float timeSqSum =
      std::inner_product(scaTimePerRepitition.begin(), scaTimePerRepitition.end(), scaTimePerRepitition.begin(), 0.0);
  float timeStDev = std::sqrt(timeSqSum / scaTimePerRepitition.size() - timeMean * timeMean);

  cout << "   Scalar sumField is: " << sumField << endl;
  cout << "   totScaTime is: " << timeSum << endl;
  cout << "   Mean time is: " << timeMean * 1e+9 << "ns" << endl;
  cout << "   Standard devi. is: " << timeStDev * 1e+9 << "ns" << endl;
  // return clock1InFloat;
  return timeSum;
}

float TimeVector(MagField &m1, const vector<ThreeVector> &posVec, vector<ThreeVector> &outputVec, const size_t &n)
{
  cout << "\n=== Vector field start: " << endl;
  float tmin  = FLT_MAX;
  float tmax  = -FLT_MAX;
  size_t imin = 0, imax = 0;
  ThreeVector xyzFieldS, sumField(0.);
  // decides no. of doubles that one Vc vector can contain.
  // depends on architecture. 4 for avx. Later can be modified
  // to take the value itself from architecture
  vector<float> vecTimePerRepitition;
  size_t noRunsAvg = 16;

  size_t inputVcLen = ceil(((float)n) / Geant::kVecLenF);
  ThreeVecSimdF_t inputForVec;
  // We read the field in float
  ThreeVecSimdF_t xyzField;
  // Then convert to double for later use
  ThreeVecSimd_t xyzField1, xyzField2;
  
  for (size_t k = 0; k < noRunsAvg; ++k) {
    clock_t clock1 = clock();
    for (size_t i = 0; i < inputVcLen; ++i) {
      // We benchmark also the AOS->SOA
      for (size_t lane = 0; lane < Geant::kVecLenF; ++lane) {
        vecCore::Set(inputForVec.x(), lane, posVec[i * Geant::kVecLenF + lane].x());
        vecCore::Set(inputForVec.y(), lane, posVec[i * Geant::kVecLenF + lane].y());
        vecCore::Set(inputForVec.z(), lane, posVec[i * Geant::kVecLenF + lane].z());
      }

      // We need the field in Double_v for further computations. We do a trick:
      // extract field values with a Float_v which we then copy to two Double_v
      m1.GetFieldValue(inputForVec, xyzField);
      // std::cout << i << ": " << inputForVec << " => " << xyzField << std::endl;
     // We benchmark also writing to the Double_v
      Geant::CopyFltToDbl(xyzField, xyzField1, xyzField2);
      for (size_t lane = 0; lane < Geant::kVecLenD; ++lane) {
        xyzFieldS.Set(vecCore::Get(xyzField1.x(), lane), vecCore::Get(xyzField1.y(), lane),
                      vecCore::Get(xyzField1.z(), lane));
        sumField += xyzFieldS;
        outputVec.push_back(xyzFieldS);
        xyzFieldS.Set(vecCore::Get(xyzField2.x(), lane), vecCore::Get(xyzField2.y(), lane),
                      vecCore::Get(xyzField2.z(), lane));
        sumField += xyzFieldS;
        outputVec.push_back(xyzFieldS);
      }
    }
    clock1              = clock() - clock1;
    float clock1InFloat = ((float)clock1) / CLOCKS_PER_SEC;
    float tpercall      = clock1InFloat / n;
    if (tpercall > tmax) {
      tmax = tpercall;
      imax = k;
    }
    if (tpercall < tmin) {
      tmin = tpercall;
      imin = k;
    }
    vecTimePerRepitition.push_back(tpercall);
  }

  // Remove imin and imax measurements from the sample
  vecTimePerRepitition.erase(vecTimePerRepitition.begin() + std::max(imin, imax));
  vecTimePerRepitition.erase(vecTimePerRepitition.begin() + std::min(imin, imax));
  float timeSum  = std::accumulate(vecTimePerRepitition.begin(), vecTimePerRepitition.end(), 0.0);
  float timeMean = timeSum / vecTimePerRepitition.size();
  float timeSqSum =
      std::inner_product(vecTimePerRepitition.begin(), vecTimePerRepitition.end(), vecTimePerRepitition.begin(), 0.0);
  float timeStDev = std::sqrt(timeSqSum / vecTimePerRepitition.size() - timeMean * timeMean);

  cout << "   Vector sumField is: " << sumField << endl;
  cout << "   totVecTime is: " << timeSum << endl;
  cout << "   Mean time is: " << timeMean * 1e+9 << "ns" << endl;
  cout << "   Standard devi. is: " << timeStDev * 1e+9 << "ns" << endl;
  return timeSum;
}

int main() // int argc, char**argv)
{
  CMSmagField m1("../VecMagFieldRoutine/cmsmagfield2015.txt");
  // m1.ReadVectorData("/home/ananya/Work/MagFieldRoutine/cms2015.txt");
  // No absolute path required now.
  // input file copied to build/VecMagFieldRoutine
  /// m1.ReadVectorData("../VecMagFieldRoutine/cms2015.txt");
  // vector<ThreeVector> posVec;
  vector<ThreeVector> posVec;
  vector<ThreeVector> outputVec;

  size_t n = 1E5;
  // cout << "Give input vector size: ";
  // cin >> n;

  // srand(time(NULL));
  srand(2);
  GenVecCart(posVec, n);
  cout << "Size of posVec is: " << posVec.size() << endl;

  float Ts = TimeScalar(m1, posVec, outputVec, n);
  float Tv = TimeVector(m1, posVec, outputVec, n);

  cout << "Vector speedup: " << Ts / Tv << endl;
}
