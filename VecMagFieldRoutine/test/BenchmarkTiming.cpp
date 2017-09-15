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

#include "base/Vector.h"

#include "base/Vector3D.h"
#include "base/Global.h"
#include <Geant/VectorTypes.h>


// #include "VC_NO_MEMBER_GATHER"

#include "CMSmagField.h"
// using MagField=MagField3<float>

using namespace std;

using Double_v = Geant::Double_v;
using Float_v = Geant::Float_v;

typedef vecgeom::Vector3D<float> ThreeVector; //normal Vector3D
typedef vecgeom::Vector3D<Float_v> ThreeVecSimd_t;
typedef vecgeom::Vector<float>   VcVectorFloat;

// typedef MagVector3<float>         MagField;
// using MagField=MagVector3<float>

using MagField= CMSmagField;

const float kRMax=9000;
const float kZMax= 16000;

float RandR(){
    float r = (float) rand()/(RAND_MAX) ;
    r = r*kRMax; //because r is in range (0,9000) mm                                                                          
    return r;
}

float RandZ(){
    float z = (float) rand()/(RAND_MAX) ;
    z = z*kZMax; //range of z is between -16k and 16k                                                                         
    int sign = rand()%2; //to define the sign, since it can be both positive and negative                                     
    if (sign==0){
        z= -z;
    }
    return z;
}

void GenVecCartSubR(float &x, float &y){
    x = RandR();
    y = RandR();
    if((x*x + y*y)> kRMax*kRMax){
        GenVecCartSubR(x,y);
    }
}

void GenVecCart(ThreeVector &pos){
    float x=0,y=0;
    float z = RandZ();
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

float TimeScalar(CMSmagField &m1, const vecgeom::Vector<ThreeVector> &posVec, const int &n, const int &nRepetitions){
    ThreeVector sumXYZField(0., 0., 0.), xyzField;
    float totScaTime = 0.;
    vector<float> scaTimePerRepitition; 

    int noRunsAvg = 16;

    cout<<"Scalar fields start: "<<endl;

    for(int k=0; k< noRunsAvg; k++){
        clock_t clock1= clock();
        for(int j=0;j<nRepetitions;j++){
            for (int i = 0; i < n; ++i)
            {
                m1.GetFieldValue<float>(posVec[i], xyzField);
                sumXYZField += xyzField;
            }
        }
        clock1 = clock() - clock1;
        float clock1InFloat = ((float)clock1)/CLOCKS_PER_SEC;
        scaTimePerRepitition.push_back(clock1InFloat/n/nRepetitions);
        totScaTime += clock1InFloat;
    }
    cout<<sumXYZField<<endl;
    
    float timeSum   = std::accumulate(scaTimePerRepitition.begin(), scaTimePerRepitition.end(), 0.0);
    float timeMean  = timeSum/scaTimePerRepitition.size();
    float timeSqSum = std::inner_product(scaTimePerRepitition.begin(), scaTimePerRepitition.end(), scaTimePerRepitition.begin(), 0.0);
    float timeStDev = std::sqrt(timeSqSum/scaTimePerRepitition.size() - timeMean*timeMean);
 
    cout<<"\nScalar: "<<endl;
    // cout<<"Total time is: "<<clock1InFloat <<endl;
    // cout<<"Time per field value is : "<<clock1InFloat/(n*nRepetitions)*1e+9 << " ns "<<endl;
    cout<<"totScaTime is: "<<totScaTime<<endl;
    //cout<<"Time per call inside loop: "<<totScaTime/(n*nRepetitions*noRunsAvg)*1e+9 << " ns "<<endl;
    cout<<"Mean time is: "<<timeMean*1e+9<<"ns"<<endl;
    cout<<"Standard devi. is: "<<timeStDev*1e+9<<"ns"<<endl;
    //return clock1InFloat;
    return totScaTime/noRunsAvg;
}

    

float TimeVector(MagField &m1, const vecgeom::Vector<ThreeVector> &posVec, const int &n, const int &nRepetitions){
    cout<<"\nVector fields start: "<<endl;
    Float_v vX;
    Float_v vY;
    Float_v vZ;

    //decides no. of doubles that one Vc vector can contain.
    //depends on architecture. 4 for avx. Later can be modified
    //to take the value itself from architecture
    float totVecTime= 0.;
    vector<float> vecTimePerRepitition; 
    int noRunsAvg = 16;

    int inputVcLen = ceil(((float)n)/Geant::kVecLenF);
    ThreeVecSimd_t *inputForVec = new ThreeVecSimd_t[inputVcLen];
    int init = 0;
    
    for (int i = 0; i < n; i=i+Geant::kVecLenF){
       for (size_t j = 0; j < Geant::kVecLenF; ++j){
            vX[j]= posVec[i+j].x();
            vY[j]= posVec[i+j].y();
            vZ[j]= posVec[i+j].z();
        }
        ThreeVecSimd_t Pos;
        Pos[0] = vX;
        Pos[1] = vY;
        Pos[2] = vZ;

        inputForVec[init] = Pos;
        init++;
    }
    
    ThreeVecSimd_t sumXYZField, xyzField;

    for (int k = 0; k < noRunsAvg; ++k)
    {
        clock_t clock1= clock();
        for (int j = 0; j < nRepetitions; ++j){
            for (int i = 0; i < inputVcLen; ++i){
                m1.GetFieldValue<Float_v>(inputForVec[i], xyzField);
                sumXYZField += xyzField;
            }
        }
        clock1 = clock() - clock1;
        float clock1InFloat = ((float)clock1)/CLOCKS_PER_SEC;
        vecTimePerRepitition.push_back(clock1InFloat/n/nRepetitions);
        totVecTime += clock1InFloat;

    }


    float timeSum   = std::accumulate(vecTimePerRepitition.begin(), vecTimePerRepitition.end(), 0.0);
    float timeMean  = timeSum/vecTimePerRepitition.size();
    float timeSqSum = std::inner_product(vecTimePerRepitition.begin(), vecTimePerRepitition.end(), vecTimePerRepitition.begin(), 0.0);
    float timeStDev = std::sqrt(timeSqSum/vecTimePerRepitition.size() - timeMean*timeMean);


    cout<<sumXYZField<<endl;  
    //float clock1InFloat = ((float)clock1)/CLOCKS_PER_SEC;  
        
    cout<<"\nVector: "<<endl;
    // cout<<"Total time is: "<<clock1InFloat<<endl;
    // cout<<"Time per field value is : "<<clock1InFloat/(n*nRepetitions)*1e+9 << " ns "<<endl;
    cout<<"totVecTime is: "<<totVecTime<<endl;
    //cout<<"Time per call inside loop: "<<totVecTime/(n*nRepetitions*noRunsAvg)*1e+9 << " ns "<<endl;
    cout<<"Mean time is: "<<timeMean*1e+9<<"ns"<<endl;
    cout<<"Standard devi. is: "<<timeStDev*1e+9<<"ns"<<endl;
    return totVecTime/noRunsAvg; 
}

int main() // int argc, char**argv)
{
    CMSmagField m1("../VecMagFieldRoutine/cmsmagfield2015.txt");
    //m1.ReadVectorData("/home/ananya/Work/MagFieldRoutine/cms2015.txt");
    //No absolute path required now. 
    //input file copied to build/VecMagFieldRoutine
    /// m1.ReadVectorData("../VecMagFieldRoutine/cms2015.txt");
    //vector<ThreeVector> posVec;
    vecgeom::Vector<ThreeVector> posVec;
    
    // int n = 1e+5;
    // int nRepetitions =100;

    int n;
    cout<<"Give input vector size: ";
    cin>>n;
    int nRepetitions;
    cout<<"Give nRepetitions: ";
    cin>>nRepetitions;

    //srand(time(NULL));
    srand(2);
    GenVecCart(posVec, n);
    cout<<"Size of posVec is: "<<posVec.size()<<endl;

    float Ts= TimeScalar(m1,posVec,n,nRepetitions);
    float Tv= TimeVector(m1,posVec,n,nRepetitions);

    cout<<"Vector speedup: " << Ts/ Tv <<endl;  

}


