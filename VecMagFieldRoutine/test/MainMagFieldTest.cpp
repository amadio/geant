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
#include "ApproxEqual.h"

// #include "MagField.h"
#include "CMSmagField.h"
#include <Geant/VectorTypes.h>

// ensure asserts are compiled in                                                                                           
#undef NDEBUG
#include <cassert>

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

int main(){

    CMSmagField m1;
    //input file is copied to build/VecMagFieldRoutine using CMakeLists
    std::string inputMap("../VecMagFieldRoutine/cms2015.txt");
    
    m1.ReadVectorData(inputMap); // "../VecMagFieldRoutine/cms2015.txt");
    ReadVectorData dataMap(inputMap); // "../VecMagFieldRoutine/cms2015.txt");

    const float kRDiff = 50.;
    const float kZDiff = 200.;
    const float kRDiffInv = 1.0/kRDiff;
    const float kZDiffInv = 1.0/kZDiff;
    const float kRMax= 9000;
    const float kZMax= 16000;
    const int noZValues = 161;
    const int halfZValues = 80;
    //const float invSqrt2 = 1/sqrt(2);
    const float halfWeight = 0.5;
    const float zero = 0.;


    //(r,0,z) corresponds exactly to (r,z) in terms that xyzField obtained is same as rzField since
    //theta=0 in this case. Hence can check GetFieldValue<vecgeom::kScalar> in place of GetFieldValueTest
    //Limitation however is that can't check for points with non zero y.
    for (float r = 0; r <= kRMax; r =r+kRDiff)
    {
        for (float z = -kZMax; z <= kZMax; z = z+kZDiff)
        {
            //Checks for (r,0,z) and (r,0,z) against (r,z)
            // cout<<r<<endl;
            // cout<<z<<endl;
            vecgeom::Vector3D<float> pos1(r,zero, z);
            vecgeom::Vector3D<float> xyzField1;
            m1.GetFieldValue<float>(pos1,xyzField1);

            int i = r*kRDiffInv*noZValues + halfZValues + z*kZDiffInv;

            //cout<<"Correct index is: "<<i<<endl;
            vecgeom::Vector3D<float> rzCheckField1(dataMap.fBr[i], dataMap.fBphi[i], dataMap.fBz[i]);
            //cout<<"xyzField1: "<<xyzField1<<" vs rzCheckField1: "<<rzCheckField1<<endl;
            assert(ApproxEqual(xyzField1, rzCheckField1, r ,z, 0)); //Working for floats 
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
            m1.GetFieldValue<float>(pos2, xyzField2);

            //Say i1, i2, i3, i4
            vecgeom::Vector3D<float> rzCheckField2, rzCheckField3;
            //Now need i1, i2, i3, i4
            //for pos2 and pos5, take i1 and i2. i4 = i3 + 161. Same z, different r. so skip through as many values of z as for one r
            int i1 = r*kRDiffInv*noZValues + halfZValues + z*kZDiffInv;
            int i2 = i1 + noZValues;

            //for pos3 and pos7, take i3 and i4. Then i4 = i3+1 because same r
            //int i3 = r*kRDiffInv*noZValues + halfZValues + z*kZDiffInv;
            int i3 = i1;
            int i4 = i3 + 1;

            rzCheckField2.x() = (dataMap.fBr[i1]   + dataMap.fBr[i2]  )*halfWeight;
            rzCheckField2.y() = (dataMap.fBphi[i1] + dataMap.fBphi[i2])*halfWeight;
            rzCheckField2.z() = (dataMap.fBz[i1]   + dataMap.fBz[i2]  )*halfWeight;

            // cout<<"Checked against: "<<endl;
            // cout<<"B for i1 is: "<<dataMap.fBr[i1]<<" "<<dataMap.fBphi[i1]<<" "<<dataMap.fBz[i1]<<endl;
            // cout<<"B for i3 is: "<<dataMap.fBr[i2]<<" "<<dataMap.fBphi[i2]<<" "<<dataMap.fBz[i2]<<endl;
            // cout<<"Direct indices are: "<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<endl;
            // cout<<"xyzField2: "<<xyzField2<<" vs rzCheckField2: "<<rzCheckField2<<endl;
            assert(ApproxEqual(xyzField2, rzCheckField2, r, z, 1)); 

            m1.GetFieldValue<float>(pos3, xyzField3);

            rzCheckField3.x() = (dataMap.fBr[i3]   + dataMap.fBr[i4]  )*halfWeight;
            rzCheckField3.y() = (dataMap.fBphi[i3] + dataMap.fBphi[i4])*halfWeight;
            rzCheckField3.z() = (dataMap.fBz[i3]   + dataMap.fBz[i4]  )*halfWeight;
            // cout<<"xyzField3: "<<xyzField3<<" vs rzCheckField3: "<<rzCheckField3<<endl;
            assert(ApproxEqual(xyzField3, rzCheckField3, r, z, 2));
            //cout<<"\n"<<endl;
        }
    }

    
    //For point in middle of cell
    for (float r = 0; r < kRMax; r =r+kRDiff)
    {
        for (float z = -kZMax; z < kZMax; z = z+kZDiff)
        {
            //cout<<"r: "<<r<<" and z: "<<z<<endl;
            vecgeom::Vector3D<float> pos4(r + kRDiff*halfWeight, zero, z + kZDiff*halfWeight);
            vecgeom::Vector3D<float> xyzField4, rzCheckField4;
            m1.GetFieldValue<float>(pos4, xyzField4);

            //need to get rzcheckfield4
            //going to be average of 4 points
            int i1 = r*kRDiffInv*noZValues + halfZValues + z*kZDiffInv;
            int i2 = i1 + noZValues;
            int i3 = i1 + 1;
            int i4 = i2 + 1;

            rzCheckField4.x() = (dataMap.fBr[i1]   + dataMap.fBr[i2]   + dataMap.fBr[i3]   + dataMap.fBr[i4]  )*0.25;
            rzCheckField4.y() = (dataMap.fBphi[i1] + dataMap.fBphi[i2] + dataMap.fBphi[i3] + dataMap.fBphi[i4])*0.25;
            rzCheckField4.z() = (dataMap.fBz[i1]   + dataMap.fBz[i2]   + dataMap.fBz[i3]   + dataMap.fBz[i4]  )*0.25;

            // cout<<"Direct indices are: "<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<endl;
            
            //cout<<"xyzField4: "<<xyzField4<<" vs rzCheckField4: "<<rzCheckField4<<endl;
            assert(ApproxEqual(xyzField4, rzCheckField4, r, z, 3));  
            //cout<<"\n"<<endl;

        }
    }
    
    
    return 0;
}
