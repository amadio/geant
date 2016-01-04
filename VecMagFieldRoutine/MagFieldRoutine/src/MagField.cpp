#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "MagField.h"
#include "base/Global.h"
#include "base/Vector3D.h"
#include "base/SOA3D.h"
using namespace std;


typedef vecgeom::Vector3D<double> Vector3D;

MagField::MagField(){

}

MagField::~MagField(){

}

void MagField::ReadVectorData(string inputMap){
    string line;
    string s1,s2,s3,s4,s5,s0;
    double d1,d2,d3,d4,d5,d0;
    ifstream pFile(inputMap);
    if (pFile.is_open()){
        // getline() returns the stream. testing the stream with while returns error such as EOF
        while(getline(pFile,line)){
             // so here we know that the read was a success and that line has valid data
            stringstream ss(line);
            //parsing all the parts. s0's store the string names which are of no use to us. 
            ss>> s0>> d1>> s1>> d0>> s2>> d2>> s3>> d3>> s4>> d4>> s5>> d5;
            fRadius.push_back(d1);
            fPhi.push_back(d0);
            fZ.push_back(d2);
            fBz.push_back(d3);
            fBr.push_back(d4);
            fBphi.push_back(d5);
        }
        pFile.close();
    }
    else{
        cout<<"Unable to open file";
    }
}


void MagField::GetFieldValueRZ(const double r, const double Z, Vector3D &rzField){

    //Take care that radius and z for out of limit values take values at end points 
    double radius = min(r, kRMax);
    double z = max(min(Z, kZMax), -kZMax); //max(min(Z,Zmax), Zmin ) 

    //to make sense of the indices, consider any particular instance e.g. (25,-200)
    int rFloor = floor(radius*kRDiffInv);
    int rIndLow = rFloor*kNoZValues;
    int rIndHigh = rIndLow + kNoZValues;

    
    //if we use z-z0 in place of two loops for Z<0 and Z>0
    //z-z0 = [0,32000]
    //so indices 0 to 160 : total 161 indices for (z-z0)/200
    //i.e. we are saying:
    int zInd = floor((z-kZ0)*kZDiffInv);
    //need i1,i2,i3,i4 for 4 required indices
    int i1 = rIndLow + zInd;
    int i2 = i1 + 1;
    int i3 = rIndHigh + zInd;
    int i4 = i3 + 1;
    double zLow = (zInd- kHalfZValues)*kZDiff; //80 because it's the middle index in 0 to 160
    double zHigh = zLow + kZDiff;
    double radiusLow = rFloor*kRDiff;
    double radiusHigh = radiusLow + kRDiff;
    //cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<endl;

    //now write function
    double a1 = (radiusHigh - radius)*(zHigh - z); //area to be multiplied with i1
    double a2 = (radiusHigh - radius)*(z - zLow);
    double a3 = (radius - radiusLow)*(zHigh - z);
    double a4 = (radius - radiusLow)*(z- zLow);

    double BR   = (fBr[i1]  *a1 + fBr[i2]  *a2 + fBr[i3]  *a3 + fBr[i4]  *a4)*kAInverse;
    double BZ   = (fBz[i1]  *a1 + fBz[i2]  *a2 + fBz[i3]  *a3 + fBz[i4]  *a4)*kAInverse;
    double BPhi = (fBphi[i1]*a1 + fBphi[i2]*a2 + fBphi[i3]*a3 + fBphi[i4]*a4)*kAInverse;

    //To make it thread safe. Because the previous predicted_B* vectors weren't threadsafe
    rzField.x()= BR;
    rzField.y()= BPhi;
    rzField.z()= BZ;
}

void MagField::GetFieldValueRZ(std::vector<double> radius, std::vector<double> z){
    for (int i = 0; i < radius.size(); ++i)
    {   
        Vector3D rzField;
        GetFieldValueRZ(radius[i], z[i], rzField);
    }
}

//Sidenote: For theta =0; xyzField = rzField. 
//theta =0 corresponds to y=0
void MagField::GetFieldValueXYZ(const Vector3D &pos, Vector3D &xyzField){

    double cyl[2];
    CartesianToCylindrical(pos, cyl); 
    Vector3D rzField;
    GetFieldValueRZ(cyl[0], cyl[1], rzField); //cyl[2] =[r,z]
    
    double sinTheta=0.0, cosTheta=1.0; //initialize as theta=0
    //To take care of r =0 case 
    if (cyl[0] != 0.0){
        double rInv = 1/cyl[0]; 
        sinTheta = pos.y()*rInv;
        cosTheta = pos.x()*rInv;
    }

    CylindricalToCartesian(rzField, sinTheta, cosTheta, xyzField);
}


void MagField::GetFieldValueTest(const Vector3D &pos, Vector3D &rzField){
    double cyl[2];
    CartesianToCylindrical(pos, cyl); 
    GetFieldValueRZ(cyl[0], cyl[1], rzField); //cyl[] =[r,z]
}



void MagField::GetFieldValues(const vecgeom::SOA3D<double> &posVec, vecgeom::SOA3D<double> &fieldVec){
    for (int i = 0; i < posVec.size(); ++i)
    {
        //fill a vector3D with ith triplet for input to getFieldValue 
        Vector3D pos(posVec.x(i), posVec.y(i) ,posVec.z(i));
        Vector3D xyzField;
        GetFieldValueXYZ(pos, xyzField); // runs for 1 triplet
        //Fill SOA3D field with single field values 
        fieldVec.x(i) = xyzField.x();
        fieldVec.y(i) = xyzField.y();
        fieldVec.z(i) = xyzField.z();

    }
}
