
//#include <string>
#include <iostream>
#include <vector>
#include <cassert>
#include <ctime>
#include <cmath>

#include <ApproxEqual.h>
#include <Geant/VectorTypes.h>
#include <base/Vector3D.h>
#include <base/Global.h>
#include "CMSmagField.h"

#undef NDEBUG
//#define VERBOSE 1

using namespace std;

using Double_v      = Geant::Double_v;
using ThreeVector   = vecgeom::Vector3D<double>;
using ThreeVector_v = vecgeom::Vector3D<Double_v>;

template <typename T>
using vector_t        = std::vector<T>;
constexpr float tesla = geant::tesla;
// constexpr float kilogauss = geant::kilogauss;
constexpr float millimeter = geant::millimeter;

const double kRMax = 9000 * millimeter;
const double kZMax = 16000 * millimeter;

double RandR()
{
  double rnd = (double)rand() / (RAND_MAX);
  return rnd * kRMax;
}

double RandZ()
{
  double rnd = (double)rand() / (RAND_MAX);
  return -kZMax + 2 * rnd * kZMax;
}

void GenVecCart(ThreeVector &pos)
{
  double rnd = (double)rand() / (RAND_MAX);
  double phi = 2. * M_PI * rnd;
  double r   = RandR();
  double x   = r * vecCore::math::Cos(phi);
  double y   = r * vecCore::math::Sin(phi);
  double z   = RandZ();
  pos.Set(x, y, z);
}

void GenVecCart(vector_t<ThreeVector> &posVec, const int &n)
{
  for (int i = 0; i < n; ++i) {
    ThreeVector pos;
    GenVecCart(pos);
    posVec.push_back(pos);
  }
}

int main()
{
  CMSmagField m1;
  m1.ReadVectorData("../VecMagFieldRoutine/cmsmagfield2015.txt");
  vector_t<ThreeVector> posVec;

  size_t n = 4;

  srand(time(NULL));
  // srand(2);
  GenVecCart(posVec, n);

  ThreeVector sumXYZField(0.), sumXYZFieldVec(0.), xyzField;
  vector_t<ThreeVector> outputScalar;

#ifdef VERBOSE
  cout << "Size of posVec is: " << posVec.size() << endl;
  cout << "Scalar fields start: " << endl;
#endif

  for (size_t i = 0; i < n; ++i) {
    m1.GetFieldValue(posVec[i], xyzField);
#ifdef VERBOSE
    cout << "   point: " << posVec[i] << "   field: " << xyzField / tesla << endl;
#endif
    sumXYZField += xyzField;
    outputScalar.push_back(xyzField);
  }
  cout << "   Scalar field sum: " << sumXYZField / tesla << " [Tesla]" << endl;

#ifdef VERBOSE
  cout << "\nVector fields start: " << endl;
#endif
  size_t inputVcLen          = ceil(((double)n) / Geant::kVecLenD);
  ThreeVector_v *inputForVec = new ThreeVector_v[inputVcLen];
  size_t init                = 0;
  for (size_t i = 0; i < n; i += Geant::kVecLenD) {
    for (size_t j = 0; j < Geant::kVecLenD; ++j) {
      vecCore::Set(inputForVec[init].x(), j, posVec[i + j].x());
      vecCore::Set(inputForVec[init].y(), j, posVec[i + j].y());
      vecCore::Set(inputForVec[init].z(), j, posVec[i + j].z());
    }
    init++;
  }

  // vector_t<ThreeVector_v> outputVec;
  ThreeVector_v *outputVec = new ThreeVector_v[inputVcLen];
  ThreeVector_v sumXYZField_v(Double_v(0));
  for (size_t i = 0; i < inputVcLen; ++i) {
    m1.GetFieldValue(inputForVec[i], outputVec[i]);
#ifdef VERBOSE
    cout << "   point: " << inputForVec[i] << "   field: " << outputVec[i] / tesla << endl;
#endif
    sumXYZField_v += outputVec[i];
  }

  sumXYZFieldVec.Set(vecCore::ReduceAdd(sumXYZField_v.x()), vecCore::ReduceAdd(sumXYZField_v.y()),
                     vecCore::ReduceAdd(sumXYZField_v.z()));
  cout << "   Vector field sum (after ReduceAdd): " << sumXYZFieldVec / tesla << " [Tesla]" << endl;
  assert(ApproxEqual(sumXYZField / tesla, sumXYZFieldVec / tesla));

  // Now compare the results scalar/vector
  for (size_t i = 0; i < inputVcLen; ++i) {
    for (size_t lane = 0; lane < Geant::kVecLenD; ++lane) {
      // ThreeVector testVec2(xyzField_v[0][j], xyzField_v[1][j], xyzField_v[2][j]);
      size_t k = i * Geant::kVecLenD + lane;
      ThreeVector testVec(vecCore::Get(outputVec[i].x(), lane), vecCore::Get(outputVec[i].y(), lane),
                          vecCore::Get(outputVec[i].z(), lane));
#ifdef VERBOSE
      cout << k << ": " << testVec / tesla << " being tested against " << outputScalar[k] / tesla << endl;
#endif
      assert(ApproxEqual(testVec / tesla, outputScalar[k] / tesla));
      k++;
    }
  }

  std::cout << "=== VectorAgainstScalar: success\n";
  return 0;
}
