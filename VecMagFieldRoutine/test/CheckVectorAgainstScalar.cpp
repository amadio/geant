
#include <string>
#include <iostream>
#include <vector>
#include <cassert>
#include <ctime>
#include <cmath> //for sqrt
#include <stdlib.h>

#include "base/Vector.h"

//#include "test/unit_tests/ApproxEqual.h"
#include "ApproxEqual.h"

#include "base/Vector3D.h"
#include "base/Global.h"

#include "CMSmagField.h"
#include <Geant/VectorTypes.h>

#undef NDEBUG

typedef float dataType;
//typedef double dataType;

using namespace std;
typedef vecgeom::Vector3D<dataType> ThreeVector; //normal Vector3D
using Double_v = Geant::Double_v;
using Float_v = Geant::Float_v;

typedef vecgeom::Vector3D<float>    ThreeVector_f;
typedef vecgeom::Vector3D<double>   ThreeVector_d;

typedef vecgeom::Vector3D<Float_v> ThreeVecSimd_v;
typedef vecgeom::Vector<dataType> VcVectorFloat;
typedef vecgeom::Vector<ThreeVecSimd_v> VecGeomVector;



const dataType kRMax=9000;
const dataType kZMax= 16000;

dataType RandR(){
    dataType r = (dataType) rand()/(RAND_MAX) ;
    r = r*kRMax; //because r is in range (0,9000) mm                                                                          
    return r;
}

dataType RandZ(){
    dataType z = (dataType) rand()/(RAND_MAX) ;
    z = z*kZMax; //range of z is between -16k and 16k                                                                         
    int sign = rand()%2; //to define the sign, since it can be both positive and negative                                     
    if (sign==0){
            z= -z;
    }
    return z;
}

void GenVecCartSubR(dataType &x, dataType &y){
    x = RandR();
    y = RandR();
    if((x*x + y*y)> kRMax*kRMax){
        GenVecCartSubR(x,y);
    }
}

void GenVecCart(ThreeVector &pos){
    dataType x=0,y=0;
    dataType z = RandZ();
    GenVecCartSubR(x, y);
    pos.x()=x;
    pos.y()=y;
    pos.z()=z;
}

void GenVecCart(vecgeom::Vector<ThreeVector> &posVec, const int &n){
    for (int i = 0; i < n; ++i)
    {       
        ThreeVector pos;
        GenVecCart(pos);
        posVec.push_back(pos);


    }
}

int main()
{
    CMSmagField m1;
    m1.ReadVectorData("../VecMagFieldRoutine/cms2015.txt");
    vecgeom::Vector<ThreeVector> posVec;
    
    int n = 1e+4;

    srand(time(NULL));
    //srand(2);
    GenVecCart(posVec, n);
    cout<<"Size of posVec is: "<<posVec.size()<<endl;

    ThreeVector sumXYZField(0., 0., 0.), xyzField;
    vector<ThreeVector> outputScalar;

    vecgeom::Vector<ThreeVector> outputScalar2;
    cout<<"Scalar fields start: "<<endl;

    for (int i = 0; i < n; ++i)
    {
        m1.GetFieldValue<float>(posVec[i], xyzField);
        // m1.GetFieldValue(posVec[i], xyzField);        
        sumXYZField += xyzField;
        outputScalar.push_back(xyzField);
        outputScalar2.push_back(xyzField);
    }
    cout<<sumXYZField<<endl;
    for (int i = 0; i < 8; ++i)
    {
        cout<<outputScalar2[i]<<endl;
    }
 

    cout<<"\nVector fields start: "<<endl;
    Float_v vX;
    Float_v vY;
    Float_v vZ;

    int inputVcLen = ceil(((dataType)n)/Geant::kVecLenF);
    ThreeVecSimd_v *inputForVec = new ThreeVecSimd_v[inputVcLen];
    //std::vector<ThreeVecSimd_v> outputVec;
    int init = 0;
    
    for (int i = 0; i < n; i=i+Geant::kVecLenF){
       for (size_t j = 0; j < Geant::kVecLenF; ++j){
            vX[j]= posVec[i+j].x();
            vY[j]= posVec[i+j].y();
            vZ[j]= posVec[i+j].z();
        }
        ThreeVecSimd_v Pos;
        Pos[0] = vX;
        Pos[1] = vY;
        Pos[2] = vZ;

        inputForVec[init] = Pos;
        init++;
    }

    //==================================================
    //=================Test Block=======================
    // ThreeVecSimd_v v1, v2, v3;
    // // outputVec.push_back(v1);
    // // outputVec.push_back(v2);
    // // outputVec.push_back(v3);
    // VecGeomVec VecGeomVec1;
    // VecGeomVec1.push_back(v1);
    // VecGeomVec1.push_back(v2);
    // VecGeomVec1.push_back(v3);

    //==================================================

    VecGeomVector outputVec;
    ThreeVecSimd_v sumXYZField_v, xyzField_v;
    for (int i = 0; i < inputVcLen; ++i){
        m1.GetFieldValue<Float_v>(inputForVec[i], xyzField_v);
        outputVec.push_back(xyzField_v);
        sumXYZField_v += xyzField_v;
    }
    cout<<sumXYZField<<endl;

    cout<<outputVec[0]<<endl;

    for (int i = 0, k=0; i < 256 ; ++i)
    {
        for (size_t j = 0; j < Geant::kVecLenF; ++j)
        {
            //ThreeVector testVec2(xyzField_v[0][j], xyzField_v[1][j], xyzField_v[2][j]);
            cout<<k<<endl;
            ThreeVector testVec(outputVec[i][0][j],outputVec[i][1][j] ,outputVec[i][2][j] );
            cout<<testVec<<" being tested against "<<outputScalar[k]<<endl;
            assert(ApproxEqual(testVec, outputScalar[k] ));
            k++;
        }
       
    }
    


    return 0;

}


