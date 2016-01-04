#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath> //for sqrt
#include <iostream>

#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"
//#include "test/unit_tests/ApproxEqual.h"
#include "ApproxEqual.h"

#include <Vc/Vc>
#include "backend/vc/Backend.h"
#include "backend/vcfloat/Backend.h"
#include "backend/scalarfloat/Backend.h"

// #include "MagField.h"
#include "CMSmagField.h"

// ensure asserts are compiled in                                                                                           
#undef NDEBUG
#include <cassert>

typedef vecgeom::Vector3D<double> ThreeVectorD;
typedef vecgeom::Vector3D<float > ThreeVectorF;

using namespace std;

class ReadVectorData
{
public: 
    vector<float> fRadius, fPhi, fZ, fBr, fBz, fBphi;
    ReadVectorData(string inputMap){
        dataFile = inputMap;
        PleaseReadData();
    };
    
    ~ReadVectorData(){}

private:
    string dataFile;
    void PleaseReadData(){
        string line;
        string s1,s2,s3,s4,s5,s0;
        float d1,d2,d3,d4,d5,d0;
        ifstream pFile(dataFile);
        if (pFile.is_open()){
            while(getline(pFile,line)){
                stringstream ss(line);
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
            cout<<"Unable to open file"<<endl;
        }
    }

};


ThreeVectorD diff(ThreeVectorD doubleVec, ThreeVectorF floatVec){
    ThreeVectorD floatToDouble((double) floatVec.x(),(double) floatVec.y(), (double) floatVec.z() );
    ThreeVectorD diff = doubleVec - floatToDouble;

    //Calculate relative error:
    double abs_x = abs(diff.x());
    double abs_y = abs(diff.y());
    double abs_z = abs(diff.z());
    double norm_x = abs(floatVec.x()) + abs(floatToDouble.x());
    double norm_y = abs(floatVec.y()) + abs(floatToDouble.y());
    double norm_z = abs(floatVec.z()) + abs(floatToDouble.z());
    ThreeVectorD rel_diff(abs_x/norm_x, abs_y/norm_y, abs_z/norm_z);
    return rel_diff*1e+6;
};

int main(){

    MagField<float> m1;
    MagField<double> m2;
    //input file is copied to build/VecMagFieldRoutine using CMakeLists
    m1.ReadVectorData("../VecMagFieldRoutine/cms2015.txt");
    m2.ReadVectorData("../VecMagFieldRoutine/cms2015.txt");

    ReadVectorData dataMap("../VecMagFieldRoutine/cms2015.txt");

    const float kRDiff = 50;
    const float kZDiff = 200;
    const float kRDiffInv = 1.0/kRDiff;
    const float kZDiffInv = 1.0/kZDiff;
    const float kRMax= 9000;
    const float kZMax= 16000;
    const int noZValues = 161;
    const int halfZValues = 80;
    //const float invSqrt2 = 1/sqrt(2);
    const float halfWeight = 0.5;
    const float zero = 0.;

    vector<ThreeVectorD> vRelDiffR;
    vector<ThreeVectorD> vRelDiffZ;
    vector<ThreeVectorD> vRelDiffRZ;
    ofstream myfile1, myfile2, myfile3;
    myfile1.open("vRelDiffR.txt");
    myfile2.open("vRelDiffZ.txt");
    myfile3.open("vRelDiffRZ.txt");

    //(r,0,z) corresponds exactly to (r,z) in terms that xyzField obtained is same as rzField since
    //theta=0 in this case. Hence can check GetFieldValue<vecgeom::kScalar> in place of GetFieldValueTest
    //Limitation however is that can't check for points with non zero y.
    for (float r = 0; r <= kRMax; r =r+kRDiff)
    {
        for (float z = -kZMax; z <= kZMax; z = z+kZDiff)
        {
            //Checks for (r,0,z) and (r,0,z) against (r,z)
            vecgeom::Vector3D<float> pos1(r,zero, z);
            vecgeom::Vector3D<float> xyzField1;
            m1.GetFieldValue<vecgeom::kScalarFloat>(pos1,xyzField1);

            int i = r*kRDiffInv*noZValues + halfZValues + z*kZDiffInv;

            vecgeom::Vector3D<float> rzCheckField1(dataMap.fBr[i], dataMap.fBphi[i], dataMap.fBz[i]);
            assert(ApproxEqual(xyzField1, rzCheckField1, r ,z, 0 )); //Working for floats 
        }
    }
    

    //Check for points on mid of cell lines i.e. (r/2,0,z) , (r,0,z/2) 

    for (float r = 0; r < kRMax; r =r+kRDiff)
    {
        for (float z = -kZMax; z < kZMax; z = z+kZDiff)
        {
            //cout<<"r: "<<r<<" and z: "<<z<<endl;

            vecgeom::Vector3D<float> pos2(r+kRDiff*halfWeight,zero,z),  pos3(r,zero,z+kZDiff*halfWeight);
            vecgeom::Vector3D<float> xyzField2,  xyzField3;
            m1.GetFieldValue<vecgeom::kScalarFloat>(pos2, xyzField2);

            vecgeom::Vector3D<double> xyzField2d,  xyzField3d;
            vecgeom::Vector3D<double> pos2d(r+kRDiff*halfWeight,zero,z), pos3d(r,zero,z+kZDiff*halfWeight);
            m2.GetFieldValue<vecgeom::kScalar>(pos2d, xyzField2d);
            
            ThreeVectorD relDiffR;
            relDiffR = diff(xyzField2d, xyzField2);
            vRelDiffR.push_back(relDiffR);
            myfile1 << "For r= "<<r<<" and z= "<<z<<" , rel. diff. is: "<<relDiffR<<"\n";

            vecgeom::Vector3D<float> rzCheckField2, rzCheckField3;
          
            int i1 = r*kRDiffInv*noZValues + halfZValues + z*kZDiffInv;
            int i2 = i1 + noZValues;

            int i3 = i1;
            int i4 = i3 + 1;

            rzCheckField2.x() = (dataMap.fBr[i1]   + dataMap.fBr[i2]  )*halfWeight;
            rzCheckField2.y() = (dataMap.fBphi[i1] + dataMap.fBphi[i2])*halfWeight;
            rzCheckField2.z() = (dataMap.fBz[i1]   + dataMap.fBz[i2]  )*halfWeight;

            assert(ApproxEqual(xyzField2, rzCheckField2, r, z, 1)); 

            m1.GetFieldValue<vecgeom::kScalarFloat>(pos3, xyzField3);
            m2.GetFieldValue<vecgeom::kScalar>(pos3d, xyzField3d);
           
            ThreeVectorD relDiffZ;
            relDiffZ = diff(xyzField3d, xyzField3);
            vRelDiffZ.push_back(relDiffZ);
            myfile2 << "For r= "<<r<<" and z= "<<z<<" , rel. diff. is: "<<relDiffZ<<"\n";

            rzCheckField3.x() = (dataMap.fBr[i3]   + dataMap.fBr[i4]  )*halfWeight;
            rzCheckField3.y() = (dataMap.fBphi[i3] + dataMap.fBphi[i4])*halfWeight;
            rzCheckField3.z() = (dataMap.fBz[i3]   + dataMap.fBz[i4]  )*halfWeight;
            assert(ApproxEqual(xyzField3, rzCheckField3, r, z, 2));
        }
    }

    
    //For point in middle of cell
    for (float r = 0; r < kRMax; r =r+kRDiff)
    {
        for (float z = -kZMax; z < kZMax; z = z+kZDiff)
        {
            vecgeom::Vector3D<float> pos4(r + kRDiff*halfWeight, zero, z + kZDiff*halfWeight);
            vecgeom::Vector3D<float> xyzField4, rzCheckField4;
            m1.GetFieldValue<vecgeom::kScalarFloat>(pos4, xyzField4);

            vecgeom::Vector3D<double> pos4d(r + kRDiff*halfWeight, zero, z + kZDiff*halfWeight);
            vecgeom::Vector3D<double> xyzField4d;
            m2.GetFieldValue<vecgeom::kScalar>(pos4d, xyzField4d);

            ThreeVectorD relDiffRZ;
            relDiffRZ = diff(xyzField4d, xyzField4);
            vRelDiffRZ.push_back(relDiffRZ);
            myfile3 << "For r= "<<r<<" and z= "<<z<<" , rel. diff. is: "<<relDiffRZ<<"\n";

            //need to get rzcheckfield4
            //going to be average of 4 points
            int i1 = r*kRDiffInv*noZValues + halfZValues + z*kZDiffInv;
            int i2 = i1 + noZValues;
            int i3 = i1 + 1;
            int i4 = i2 + 1;

            rzCheckField4.x() = (dataMap.fBr[i1]   + dataMap.fBr[i2]   + dataMap.fBr[i3]   + dataMap.fBr[i4]  )*0.25;
            rzCheckField4.y() = (dataMap.fBphi[i1] + dataMap.fBphi[i2] + dataMap.fBphi[i3] + dataMap.fBphi[i4])*0.25;
            rzCheckField4.z() = (dataMap.fBz[i1]   + dataMap.fBz[i2]   + dataMap.fBz[i3]   + dataMap.fBz[i4]  )*0.25;

            assert(ApproxEqual(xyzField4, rzCheckField4, r, z, 3));  
        }
    }
        
    myfile1.close();
    myfile2.close();
    myfile3.close();
    return 0;
}
