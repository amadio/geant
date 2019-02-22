#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "RZMagFieldFromMap.h"
#include "base/Global.h"
#include "base/Vector3D.h"
#include "base/SOA3D.h"
using namespace std;

// typedef vecgeom::Vector3D<double> Vector3D;
using  Vector3D_d = vecgeom::Vector3D<double>;

RZMagFieldFromMap::RZMagFieldFromMap()
   : VVectorField( 3, true)   //  meaning: number of components,   changes Energy
{
}

RZMagFieldFromMap::RZMagFieldFromMap(std::string inputMap)
   : RZMagFieldFromMap()
{
  ReadVectorData(inputMap);   
}

RZMagFieldFromMap::~RZMagFieldFromMap()
{
}

void RZMagFieldFromMap::ReadVectorData(string inputMap)
{
  string line;
  string s1, s2, s3, s4, s5, s0;
  double d1, d2, d3, d4, d5, d0;
  ifstream pFile(inputMap);
  if (pFile.is_open()) {
    // getline() returns the stream. testing the stream with while returns error such as EOF
    while (getline(pFile, line)) {
      // so here we know that the read was a success and that line has valid data
      stringstream ss(line);
      // parsing all the parts. s0's store the string names which are of no use to us.
      ss >> s0 >> d1 >> s1 >> d0 >> s2 >> d2 >> s3 >> d3 >> s4 >> d4 >> s5 >> d5;
      fRadius.push_back(d1);
      fPhi.push_back(d0);
      fZ.push_back(d2);
      fBz.push_back(d3);
      fBr.push_back(d4);
      fBphi.push_back(d5);
    }
    pFile.close();
  } else {
    cout << "Unable to open file";
  }
}

void RZMagFieldFromMap::
  GetFieldValueRZ(const double r, const double Z, Vector3D_d &rzField)
                  // vecgeom::Vector3D<double> &rzField)
{

  // Take care that radius and z for out of limit values take values at end points
  double radius = min(r, kRMax);
  double z      = max(min(Z, kZMax), -kZMax); // max(min(Z,Zmax), Zmin )

  // to make sense of the indices, consider any particular instance e.g. (25,-200)
  int rFloor   = floor(radius * kRDiffInv);
  int rIndLow  = rFloor * kNoZValues;
  int rIndHigh = rIndLow + kNoZValues;

  // if we use z-z0 in place of two loops for Z<0 and Z>0
  // z-z0 = [0,32000]
  // so indices 0 to 160 : total 161 indices for (z-z0)/200
  // i.e. we are saying:
  int zInd = floor((z - kZ0) * kZDiffInv);
  // need i1,i2,i3,i4 for 4 required indices
  int i1            = rIndLow + zInd;
  int i2            = i1 + 1;
  int i3            = rIndHigh + zInd;
  int i4            = i3 + 1;
  double zLow       = (zInd - kHalfZValues) * kZDiff; // 80 because it's the middle index in 0 to 160
  double zHigh      = zLow + kZDiff;
  double radiusLow  = rFloor * kRDiff;
  double radiusHigh = radiusLow + kRDiff;
  // cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<endl;

  // now write function
  double a1 = (radiusHigh - radius) * (zHigh - z); // area to be multiplied with i1
  double a2 = (radiusHigh - radius) * (z - zLow);
  double a3 = (radius - radiusLow) * (zHigh - z);
  double a4 = (radius - radiusLow) * (z - zLow);

  double BR   = (fBr[i1] * a1 + fBr[i2] * a2 + fBr[i3] * a3 + fBr[i4] * a4) * kAInverse;
  double BZ   = (fBz[i1] * a1 + fBz[i2] * a2 + fBz[i3] * a3 + fBz[i4] * a4) * kAInverse;
  double BPhi = (fBphi[i1] * a1 + fBphi[i2] * a2 + fBphi[i3] * a3 + fBphi[i4] * a4) * kAInverse;

  // To make it thread safe. Because the previous predicted_B* vectors weren't threadsafe
  rzField.x() = BR;
  rzField.y() = BPhi;
  rzField.z() = BZ;
}

void RZMagFieldFromMap::GetFieldValueRZ(std::vector<double> radius, std::vector<double> z)
{
  int len= std::min( radius.size(), z.size() );
  for (int i = 0; i < len; ++i) {
    Vector3D_d rzField;
    GetFieldValueRZ(radius[i], z[i], rzField);
  }
}

// Sidenote: For theta =0; xyzField = rzField.
// theta =0 corresponds to y=0
void RZMagFieldFromMap::GetFieldValueXYZ(const Vector3D_d &pos, Vector3D_d &xyzField)
{

  double cyl[2];
  CartesianToCylindrical(pos, cyl);
  Vector3D_d rzField;
  GetFieldValueRZ(cyl[0], cyl[1], rzField); // cyl[2] =[r,z]

  double sinTheta = 0.0, cosTheta = 1.0; // initialize as theta=0
  // To take care of r =0 case
  if (cyl[0] != 0.0) {
    double rInv = 1 / cyl[0];
    sinTheta    = pos.y() * rInv;
    cosTheta    = pos.x() * rInv;
  }

  CylindricalToCartesian(rzField, sinTheta, cosTheta, xyzField);
}

void RZMagFieldFromMap::GetFieldValueTest(const Vector3D_d &pos, Vector3D_d &rzField)
{
  double cyl[2];
  CartesianToCylindrical(pos, cyl);
  GetFieldValueRZ(cyl[0], cyl[1], rzField); // cyl[] =[r,z]
}

void RZMagFieldFromMap::GetFieldValues(const vecgeom::SOA3D<double> &posVec,
                                                   vecgeom::SOA3D<double> &fieldVec)
{
  int len= posVec.size();
  for (int i = 0; i < len; ++i) {
    // fill a vector3D with ith triplet for input to getFieldValue
    Vector3D_d pos(posVec.x(i), posVec.y(i), posVec.z(i));
    Vector3D_d xyzField;
    GetFieldValueXYZ(pos, xyzField); // runs for 1 triplet
    // Fill SOA3D field with single field values
    fieldVec.x(i) = xyzField.x();
    fieldVec.y(i) = xyzField.y();
    fieldVec.z(i) = xyzField.z();
  }
}

/** @brief Vector interface for field retrieval */
void RZMagFieldFromMap::
  ObtainFieldValueSIMD(const Vector3D<Double_v> &positionVec, Vector3D<Double_v> &fieldValueVec)
{
   Vector3D<double> position, fieldValue;
   int vecsz = vecCore::VectorSize<Double_v>();  // Deprecated !?
   // Double_v dv;
   // int vecsz = VectorSize<decltype(dv)>();
   Double_v fieldValueArr[3];
   
   for( int i= 0; i < vecsz ; ++i)           
   {
      position= Vector3D_d( vecCore::Get(positionVec[0], i),
                            vecCore::Get(positionVec[1], i),
                            vecCore::Get(positionVec[2], i) ) ;
      
      ObtainFieldValue( position, fieldValue);
      for( int j= 0; j < 3; j++)
        vecCore::Set( fieldValueArr[j], i, fieldValue[j] );

      vecCore::Set( fieldValueVec[0], i, fieldValue.x() );
      vecCore::Set( fieldValueVec[0], i, fieldValue.x() );
      vecCore::Set( fieldValueVec[0], i, fieldValue.x() );            
   }
}

RZMagFieldFromMap::
  RZMagFieldFromMap(const RZMagFieldFromMap &right)
   : RZMagFieldFromMap()
{
   fRadius = right.fRadius;
   fPhi    = right.fPhi;
   fZ      = right.fZ;
   fBr     = right.fBr;
   fBz     = right.fBz;
   fBphi   = right.fBphi;
}
