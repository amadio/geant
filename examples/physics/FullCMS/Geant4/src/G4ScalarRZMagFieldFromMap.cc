//  Created by J. Apostolakis,  18 May 2018
//
//    based on ScalarRZMagFieldFromMap
//          which was based on the work of Ananya Fall 2015
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>

#include <cassert>

#include "G4Types.hh"
#include "G4ScalarRZMagFieldFromMap.hh"

using namespace std;

G4ScalarRZMagFieldFromMap::G4ScalarRZMagFieldFromMap()
   : G4MagneticField()  
{
}

G4ScalarRZMagFieldFromMap::G4ScalarRZMagFieldFromMap(std::string inputMap)
   : G4ScalarRZMagFieldFromMap()
{
  ReadVectorData(inputMap);   
}

G4ScalarRZMagFieldFromMap::~G4ScalarRZMagFieldFromMap()
{
}

void G4ScalarRZMagFieldFromMap::ReadVectorData(string inputMap)
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
    G4cout << "ReadVectorData> Data read - closing file " << inputMap << G4endl;
  } else {
    G4cerr << "G4ScalarRZMagFieldFromMap::ReadVectorData> Unable to open file"
           << inputMap << G4endl;
    exit(1);
  }
}

void G4ScalarRZMagFieldFromMap::
  GetFieldValueRZ(double r, double Zin, G4ThreeVector &rzField) const
{

  // Take care that radius and z for out of limit values take values at end points
  double radius = min(r, kRMax);
  double z      = max(min(Zin, kZMax), -kZMax); // max(min(Z,Zmax), Zmin )

  // to make sense of the indices, consider any particular instance e.g. (25,-200)
  int rFloor   = std::min( (int) floor(radius * kRDiffInv), kNoRValues - 2 );
  int rIndLow  = rFloor * kNoZValues;
  int rIndHigh = rIndLow + kNoZValues;  
  // int rIndHigh = (radius < kRMax) ? rIndLow + kNoZValues ? rIndLow;

  // if we use z-z0 in place of two loops for Z<0 and Z>0
  // z-z0 = [0,32000]
  // so indices 0 to 160 : total 161 indices for (z-z0)/200

  int zInd = std::min( (int) std::floor(z * kZDiffInv) + kHalfZValues, kNoZValues - 2);

  // cout << " r:  floor = " << rFloor << " Index = " << rIndLow << " hi= " << rIndHigh << endl;
  // cout << " z-Index = " << zInd    << endl;

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

  unsigned long minSzUL= std::min( std::min( fBr.size(), fBphi.size()) , fBz.size() );

  assert( 0 <= i1 );
  assert( i4 <= (int) minSzUL ) ;
  assert( 0. <= a1  && a1 * kAInverse <= 1.0 );
  assert( 0. <= a2  && a2 * kAInverse <= 1.0 );
  assert( 0. <= a3  && a3 * kAInverse <= 1.0 );
  assert( 0. <= a4  && a4 * kAInverse <= 1.0 );  
  
  double BR   = (fBr[i1] * a1 + fBr[i2] * a2 + fBr[i3] * a3 + fBr[i4] * a4) * kAInverse;
  double BZ   = (fBz[i1] * a1 + fBz[i2] * a2 + fBz[i3] * a3 + fBz[i4] * a4) * kAInverse;
  double BPhi = (fBphi[i1] * a1 + fBphi[i2] * a2 + fBphi[i3] * a3 + fBphi[i4] * a4) * kAInverse;

  rzField = G4ThreeVector( BR, BPhi, BZ );
}

// Sidenote: For theta =0; xyzField = rzField.
// theta =0 corresponds to y=0
void G4ScalarRZMagFieldFromMap::GetFieldValueXYZ(const G4ThreeVector &pos, G4ThreeVector &xyzField) const
{
  double cyl[2];
  CartesianToCylindrical(pos, cyl);
  G4ThreeVector rzField;
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

void G4ScalarRZMagFieldFromMap::GetFieldValueTest(const G4ThreeVector &pos, G4ThreeVector &rzField)
{
  double cyl[2];
  CartesianToCylindrical(pos, cyl);
  GetFieldValueRZ(cyl[0], cyl[1], rzField); // cyl[] =[r,z]
}

/*******
void G4ScalarRZMagFieldFromMap::GetFieldValues(const vecgeom::SOA3D<double> &posVec,
                                                   vecgeom::SOA3D<double> &fieldVec)
{
  int len= posVec.size();
  for (int i = 0; i < len; ++i) {
    // fill a vector3D with ith triplet for input to getFieldValue
    G4ThreeVector pos(posVec.x(i), posVec.y(i), posVec.z(i));
    G4ThreeVector xyzField;
    GetFieldValueXYZ(pos, xyzField); // runs for 1 triplet
    // Fill SOA3D field with single field values
    fieldVec.x(i) = xyzField.x();
    fieldVec.y(i) = xyzField.y();
    fieldVec.z(i) = xyzField.z();
  }
}
  ************/

/** @brief Vector interface for field retrieval */
/************

void G4ScalarRZMagFieldFromMap::
  ObtainFieldValueSIMD(const Vector3D<Double_v> &positionVec, Vector3D<Double_v> &fieldValueVec)
{
   Vector3D<double> position, fieldValue;
   int vecsz = vecCore::VectorSize<Double_v>();  // Deprecated !?
   // Double_v dv;
   // int vecsz = VectorSize<decltype(dv)>();
   Double_v fieldValueArr[3];
   
   for( int i= 0; i < vecsz ; ++i)           
   {
      position= G4ThreeVector( vecCore::Get(positionVec[0], i),
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
  ************/
  
G4ScalarRZMagFieldFromMap::
  G4ScalarRZMagFieldFromMap(const G4ScalarRZMagFieldFromMap &right)
   : G4ScalarRZMagFieldFromMap()
{
   fRadius = right.fRadius;
   fPhi    = right.fPhi;
   fZ      = right.fZ;
   fBr     = right.fBr;
   fBz     = right.fBz;
   fBphi   = right.fBphi;
}
