#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath> //for sqrt
#include <iostream>

#include "base/Vector3D.h"
#include "base/Global.h"
//#include "test/unit_tests/ApproxEqual.h"
#include "Geant/ApproxEqual.h"

#include "CMSmagField.h"
#include <Geant/VectorTypes.h>

// ensure asserts are compiled in
#undef NDEBUG
#include <cassert>

using namespace std;

class ReadVectorData {
public:
  vector<float> fRadius, fPhi, fZ, fBr, fBz, fBphi;
  ReadVectorData(string inputMap)
  {
    dataFile = inputMap;
    PleaseReadData();
  };

  ~ReadVectorData() {}

private:
  string dataFile;
  void PleaseReadData()
  {
     // int lineNum=0;
    string line;
    string s1, s2, s3, s4, s5, s0;
    float d1, d2, d3, d4, d5, d0;
    ifstream pFile(dataFile);
    if (pFile.is_open()) {
      while (getline(pFile, line)) {
        stringstream ss(line);
        ss >> s0 >> d1 >> s1 >> d0 >> s2 >> d2 >> s3 >> d3 >> s4 >> d4 >> s5 >> d5;
        fRadius.push_back(d1);
        fPhi.push_back(d0);
        fZ.push_back(d2);
        fBz.push_back(d3);
        fBr.push_back(d4);
        fBphi.push_back(d5);

        // if( lineNum++ < 10 ) std::cout << " Line# " << lineNum << " : " << s0 << " " << d1 << " " << s1 << " " << s2 << " " << d2 << " " << s3 << " " << d3 << " " << s4 << " " << d4 << " " << s5 << " " << d5 << endl;
      }
      pFile.close();
    } else {
      cout << "Unable to open file" << endl;
    }
  }
};

bool ReportDifference( double value, double checkVal, vecgeom::Vector3D<float> pos,                           
                       int kComp,    const char* descr );

int main()
{

  CMSmagField m1;
  // input file is copied to build/examples/magneticfield/simplifiedCMS using CMakeLists
  std::string inputMap("cmsmagfield2015.txt");
  m1.ReadVectorData(inputMap); 
  ReadVectorData dataMap(inputMap);

  static constexpr float millimeter = geant::units::millimeter;  
  const float kRDiff    = 50. * millimeter;
  const float kZDiff    = 200. * millimeter;
  const float kRDiffInv = 1.0 / kRDiff;
  const float kZDiffInv = 1.0 / kZDiff;
  const float kRMax     =  9000 * millimeter;
  const float kZMax     = 16000 * millimeter;
  const int noZValues   = 161;
  const int noRValues   = 181;
  const int halfZValues = ( noZValues - 1 ) / 2;
  // const float invSqrt2 = 1/sqrt(2);
  const float halfWeight = 0.5;
  const float zero       = 0.;
  const float threshold  = 2.0e-7;  // Above 10^-8 that corresponds to 23 bit accuracy
  using vecgeom::kTiny;
  
  assert( ApproxEqual( (double) (noZValues-1) * kZDiff, (double) 2. * kZMax ) ) ; 
  assert( ApproxEqual( (double) (noRValues-1) * kRDiff, (double)   kRMax ) ) ;
  
  const float toTesla    = 1.0 / geant::units::tesla;
  std::cout << " toTesla = " << toTesla << std::endl;

  // using TypeFlt = float;
  using TypeFlt = double;  
  
  //(r,0,z) corresponds exactly to (r,z) in terms that xyzField obtained is same as rzField since
  // theta=0 in this case. Hence can check GetFieldValue<vecgeom::kScalar> in place of GetFieldValueTest
  // Limitation however is that can't check for points with non zero y.
  // for (float r = 0; r <= kRMax; r = r + kRDiff) {

  bool allGoodPoints = true;
  int  badPoints= 0;

// #define  NEW_LOOP   1
  
  for (int jr = 0; jr <= noRValues; jr++ ) {
    float r = jr * kRDiff; // r <= kRMax; r = r + kRDiff) {

    for (int lz = 0; lz <= noZValues; lz++ ) {
      float z = -kZMax + lz * kZDiff ; // z <= kZMax; z = z + kZDiff) {

      // Checks for (r,0,z) and (r,0,z) against (r,z)
      // cout << " r= " << r << " z =  " << z << endl;
      vecgeom::Vector3D<TypeFlt> pos1(r, zero, z);
      vecgeom::Vector3D<TypeFlt> xyzField1;
      m1.EstimateFieldValues<TypeFlt>(pos1, xyzField1);

      // std::cout << "Field r= " << r << " z= " << z << " vec= " << xyzField1[1] << " " << xyzField1[2] << " " << xyzField1[0] << " internal GeantV unites" << std::endl;;
      xyzField1 *= toTesla;
      // cout << "Field r= " << r << " z= " << z << " vec= " << xyzField1[1] << " " << xyzField1[2] << " " << xyzField1[0] << " Tesla " << endl;

      int i = jr * noZValues + lz;
      
      // cout << "Expected index is: " << i << endl;
      vecgeom::Vector3D<float> rzCheckField1(dataMap.fBr[i], dataMap.fBphi[i], dataMap.fBz[i]);
      // cout << "xyzField1: " << xyzField1 << " vs rzCheckField1: " << rzCheckField1 << endl;
      bool goodPoint= 
         ApproxEqual(xyzField1, rzCheckField1, r, z, 0); // Working for floats
      // assert( goodPoint && "Comparison on location of r, z node - i.e. provided values." ); 
      allGoodPoints = allGoodPoints && goodPoint;
      if( !goodPoint ) badPoints++;
      // cout << endl;
    }
  }

  // Check for points on mid of cell lines i.e. (r/2,0,z) , (r,0,z/2)

// #define OLD_CODE 1
  
#ifdef OLD_CODE  
  for (float r = 0; r < kRMax; r = r + kRDiff) {
    for (float z = -kZMax; z < kZMax; z = z + kZDiff) {
#else  
  for (int jr = 0; jr < noRValues; jr++ ) {
    float r = jr * kRDiff; // r <= kRMax; r = r + kRDiff) {

    for (int lz = 0; lz < noZValues; lz++ ) {
      float z = -kZMax + lz * kZDiff ; // z <= kZMax; z = z + kZDiff) {
#endif  
      cout<<"r: "<<r<<" and z: "<<z<<endl;

      vecgeom::Vector3D<float> pos2(r + kRDiff * halfWeight, zero, z), pos3(r, zero, z + kZDiff * halfWeight);
      vecgeom::Vector3D<float> xyzField2, xyzField3;
      m1.EstimateFieldValues<float>(pos2, xyzField2);
      xyzField2 *= toTesla;      

      // Say i1, i2, i3, i4
      vecgeom::Vector3D<float> rzCheckField2, rzCheckField3;
      // Now need i1, i2, i3, i4
      // for pos2 and pos5, take i1 and i2. i4 = i3 + 161. Same z, different r. so skip through as many values of z as
      // for one r
#ifdef OLD_CODE      
      int i1 = std::round( r * kRDiffInv * noZValues + halfZValues + z * kZDiffInv );
#else      
      int i1 = std::max( jr  , noRValues - 1)  * noZValues + lz;
#endif
      // if( jr < noRValues - 1 ) {  // Trial control for last line
      int i2 = i1 + noZValues;
      // int i2 = std::max( jr+1, noRValues - 1)  * noZValues + lz;
      
      // for pos3 and pos4, take i3 and i4. Then i4 = i3+1 because same r
      // int i3 = r*kRDiffInv*noZValues + halfZValues + z*kZDiffInv;
      int i3 = i1;
      int i4 = i3 + 1;
      // int i4 = std::max( jr  , noRValues - 1)  * noZValues + std::max(lz+1, noZValues );
      
      rzCheckField2.x() = 0.5 * (dataMap.fBr[i1]   + dataMap.fBr[i2]) ;
      rzCheckField2.y() = 0.5 * (dataMap.fBphi[i1] + dataMap.fBphi[i2]) ;
      rzCheckField2.z() = 0.5 * (dataMap.fBz[i1]   + dataMap.fBz[i2]) ;
      
      // auto diffVec = rzCheckField2 - xyzField2;
      
      if( ( rzCheckField2 - xyzField2 ).Mag() > threshold * (rzCheckField2.Mag() + xyzField2.Mag()) ) {
         cout << "Differences seen in check of values at half-point in r.  r = " << r + 0.5 * kRDiff
              << " z= " << z << endl;         
         cout << "    xyzField:      " << xyzField2     << " vs " << endl
              << "    rzCheckField:  " << rzCheckField2 << endl
              << "    diff(val-chk): " << xyzField2-rzCheckField2 << endl;
         cout << "    relative diff: ( ";
         for( int j=0; j< 3; j++ ) {
            cout << (xyzField2[j]-rzCheckField2[j]) / ( fabs(xyzField2[j]) + fabs(rzCheckField2[j]) + kTiny );
            if( j < 2 ) cout << ", ";
         }
         cout << ")" << endl;
         
         cout << "Checked against: " << endl;
         cout << "B for i1 is: " << dataMap.fBr[i1] << " " << dataMap.fBphi[i1] << " " << dataMap.fBz[i1] << endl;
         cout << "B for i3 is: " << dataMap.fBr[i2] << " " << dataMap.fBphi[i2] << " " << dataMap.fBz[i2] << endl;
         cout << "Direct indices are: " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;
      }
      // assert(ApproxEqual(xyzField2, rzCheckField2, r, z, 1));
      // ApproxEqual(xyzField2, rzCheckField2, r, z, 1);
      for( int k= 0; k<3 ; k++ )
         ApproxEqual(xyzField2[k], rzCheckField2[k], r, z, 1) ||
            ReportDifference( xyzField2[k], rzCheckField2[k], pos2, k, "Half-point in r");
      // }  // Trial control for last line
      
      m1.EstimateFieldValues<float>(pos3, xyzField3);
      xyzField3 *= toTesla;      

      rzCheckField3.x() = (dataMap.fBr[i3] + dataMap.fBr[i4]) * halfWeight;
      rzCheckField3.y() = (dataMap.fBphi[i3] + dataMap.fBphi[i4]) * halfWeight;
      rzCheckField3.z() = (dataMap.fBz[i3] + dataMap.fBz[i4]) * halfWeight;
      // cout << "xyzField3: " << xyzField3 << " vs rzCheckField3: " << rzCheckField3 << endl;
      // assert(
      ApproxEqual(xyzField3, rzCheckField3, r, z, 2); // );
      // cout << "\n" << endl;
    }
  }

  // For point in middle of cell
  for (float r = 0; r < kRMax; r = r + kRDiff) {
    for (float z = -kZMax; z < kZMax; z = z + kZDiff) {
      // cout << "r: " << r << " and z: " << z << endl;
      vecgeom::Vector3D<float> pos4(r + kRDiff * halfWeight, zero, z + kZDiff * halfWeight);
      vecgeom::Vector3D<float> xyzField4, rzCheckField4;
      m1.EstimateFieldValues<float>(pos4, xyzField4);
      xyzField4 *= toTesla;
      
      // need to get rzcheckfield4
      // going to be average of 4 points
      int i1 = r * kRDiffInv * noZValues + halfZValues + z * kZDiffInv;
      int i2 = i1 + noZValues;
      int i3 = i1 + 1;
      int i4 = i2 + 1;

      rzCheckField4.x() = (dataMap.fBr[i1] + dataMap.fBr[i2] + dataMap.fBr[i3] + dataMap.fBr[i4]) * 0.25;
      rzCheckField4.y() = (dataMap.fBphi[i1] + dataMap.fBphi[i2] + dataMap.fBphi[i3] + dataMap.fBphi[i4]) * 0.25;
      rzCheckField4.z() = (dataMap.fBz[i1] + dataMap.fBz[i2] + dataMap.fBz[i3] + dataMap.fBz[i4]) * 0.25;

      // cout << "Direct indices are: " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;

      // cout << "xyzField4: " << xyzField4 << " vs rzCheckField4: " << rzCheckField4 << endl;
      // assert(
      // ApproxEqual(xyzField4, rzCheckField4, r, z, 3); // );
      for( int k= 0; k<3 ; k++ ) {
         ApproxEqual(xyzField4[k], rzCheckField4[k], r, z, 1) ||
            ReportDifference( xyzField4[k], rzCheckField4[k], pos4, k, "Mid-point in r, z cell");
      }
      // cout << "\n" << endl;
    }
  }

  return 0;
}


// void CheckAndReportDiffs( const Vec_t &target, 
//                           const Vec_t &baseline, const float r, const float z, const int i);

bool ReportDifference( double value,
                       double checkVal,
                       vecgeom::Vector3D<float> pos,                           
                       int kComp,
                       const char* description )
{
   vecgeom::Vector3D<float> xy ( pos.x(), pos.y(), 0.0 );
   double  z= pos.z();
   double  r= xy.Mag();
   cout << "Differences seen in check of values " << description
        << "  near point   r = " << r  << " z= " << z << endl;         
   cout << "  component " << kComp << "  field-val:  " << setw(12) << value << " vs " 
        << "    check-val:  " << setw(12) << checkVal  
        << "    diff(val-chk): " << setw(12) << value - checkVal ;

   double relativeDiff = (value-checkVal) / ( fabs(value) + fabs(checkVal) + vecgeom::kTiny );   
   cout << "    relative diff:  " << relativeDiff << endl;
        
      


   return relativeDiff < 1.0e-5;   
}
