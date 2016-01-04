
/*

Generates some random points following the previous points and some new points

*/

#include "iostream"

#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"
#include <string>
#include <vector>
#include <ctime>
#include <cmath> //for sqrt
#include <stdlib.h>
#include <Vc/Vc>
#include "backend/vc/Backend.h"
#include "backend/vcfloat/Backend.h"
#include "base/Vector.h"

#include "CMSmagField.h"

using namespace std;
typedef vecgeom::Vector3D<float> ThreeVector; //normal Vector3D
typedef vecgeom::Vector3D<vecgeom::kVcFloat::precision_v> ThreeVecSimd_t;
typedef vecgeom::Vector<float> VcVectorFloat;


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


float ExpoRandom(float mean){
    //Generate a random number between 0 and 1.
    float rLog;
    rLog = (float)rand()/(float)(RAND_MAX);

    //Put it as a parameter to the inverse distribution function.
    return  -mean*log(rLog);
}

void GenUnitVec(ThreeVector &unitVec){
    float x = (float)rand()/(float)(RAND_MAX);
    float y = (float)rand()/(float)(RAND_MAX);
    float z = (float)rand()/(float)(RAND_MAX);
    unitVec.x() = x;
    unitVec.y() = y;
    unitVec.z() = z;
    if((x*x + y*y + z*z) > 1 ){
        GenVecCart(unitVec);
    }
}

ThreeVector DelPos(int t){
    float R = 150.f ; //150 mm , 15 cm and 2.5 cm option
    float v = 50.f; //5cm 

    //exponential random variable 
    float randomExp = ExpoRandom(R/t);
    ThreeVector unitVec1, unitVec2;
    GenUnitVec(unitVec1);
    GenUnitVec(unitVec2);
    ThreeVector delta = v*t*unitVec1 + randomExp*unitVec2;
/*    cout<<"randomExp: "<<randomExp<<endl;
    cout<<"unitVec1: "<<unitVec1<<endl;
    cout<<"unitVec2: "<<unitVec2<<endl;
    cout<<"Delta is: "<<delta<<endl;*/
    return delta;
}

void GenVecCartContinuous(vecgeom::Vector<ThreeVector> &posVec, int &n){
    n = floor(n/8.);

    for (int i = 0; i < 8; ++i)
    {       
        ThreeVector pos;
        //8 new particles
        GenVecCart(pos);
        // cout<<i<<" : "<<pos<<endl;
        posVec.push_back(pos);
    }

    int t = 1;
    for (int i = 8; i < n; ++i)
    {       
        ThreeVector pos;

        //4 old particles
        for (int j = 0; j < 4; ++j)
        {
            pos = posVec[i-8] + DelPos(t);
            //cout<<i+j<<" : "<<DelPos(t)<<endl;
            // cout<<i+j<<" : "<<pos<<endl;
            //GenVecCart(pos);
            posVec.push_back(pos);
        }

        //4 new particles 
        for (int j = 0; j < 4; ++j)
        {
            GenVecCart(pos);
            // cout<<i+j+4<<" : "<<pos<<endl;
            posVec.push_back(pos);
        }
        t++;
    }

    n = n*8;
}


float TimeScalar(MagField &m1, const vecgeom::Vector<ThreeVector> &posVec, const int &n, const int &nRepetitions, ThreeVector &sumXYZField){
    ThreeVector  xyzField; //sumXYZField(0., 0., 0.),

    //cout<<"\nScalar fields start: "<<endl;
    clock_t clock1= clock();
    for(int j=0;j<nRepetitions;j++){
        for (int i = 0; i < n; ++i)
        {
            m1.GetFieldValue<vecgeom::kScalarFloat>(posVec[i], xyzField);
            sumXYZField += xyzField;
        }
    }

    clock1 = clock() - clock1;
    
    //cout<<sumXYZField<<endl;
    float clock1InFloat = ((float)clock1)/CLOCKS_PER_SEC;

/*    cout<<"\nScalar: "<<endl;
    cout<<"Total time is: "<<clock1InFloat <<endl;
    cout<<"Time per field value is : "<<clock1InFloat/(n*nRepetitions)*1e+9 << " ns "<<endl;*/

    return clock1InFloat/(n*nRepetitions)*1e+9;
}

float TimeVector(MagField &m1, const vecgeom::Vector<ThreeVector> &posVec, const int &n, const int &nRepetitions, ThreeVecSimd_t &sumXYZField){
    //cout<<"\nVector fields start: "<<endl;
    vecgeom::kVcFloat::precision_v vX;
    vecgeom::kVcFloat::precision_v vY;
    vecgeom::kVcFloat::precision_v vZ;

    //decides no. of doubles that one Vc vector can contain.
    //depends on architecture. 4 for avx. Later can be modified
    //to take the value itself from architecture
    int noOfDoubles = 8;

    int inputVcLen = ceil(((float)n)/noOfDoubles);
    ThreeVecSimd_t *inputForVec = new ThreeVecSimd_t[inputVcLen];
    int init = 0;
    
    for (int i = 0; i < n; i=i+noOfDoubles){
       for (int j = 0; j < noOfDoubles; ++j){
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
    
    clock_t clock1= clock();
    ThreeVecSimd_t xyzField; //, sumXYZField;
    for (int k = 0; k < nRepetitions; ++k){
        for (int i = 0; i < inputVcLen; ++i){
            m1.GetFieldValue<vecgeom::kVcFloat>(inputForVec[i], xyzField);
            sumXYZField += xyzField;
        }
    }
    clock1 = clock() - clock1;
    
    //cout<<sumXYZField<<endl;  
    float clock1InFloat = ((float)clock1)/CLOCKS_PER_SEC;  
        
    // cout<<"\nVector: "<<endl;
    // cout<<"Total time is: "<<clock1InFloat<<endl;
    // cout<<"Time per field value is : "<<clock1InFloat/(n*nRepetitions)*1e+9 << " ns "<<endl;

    return clock1InFloat/(n*nRepetitions)*1e+9; 
}


void CalculateMeanStDev(const vector<float> timeVec, float &mean, float &stDev){
    float sum   = std::accumulate(timeVec.begin(), timeVec.end(), 0.0);
    mean  = sum/timeVec.size();
    float sqSum = std::inner_product(timeVec.begin(), timeVec.end(), timeVec.begin(), 0.0);
    stDev = std::sqrt(sqSum/timeVec.size() - mean*mean);
    cout<<"Mean Time: "<<mean<<"ns"<<endl;
    cout<<"Standard dev: "<<mean<<"ns"<<endl;
}



int main(){

    MagField m1;
    //m1.ReadVectorData("/home/ananya/Work/MagFieldRoutine/cms2015.txt");
    //No absolute path required now. 
    //input file copied to build/VecMagFieldRoutine
    m1.ReadVectorData("../VecMagFieldRoutine/cms2015.txt");
    //vector<ThreeVector> posVec;
    
    vecgeom::Vector<ThreeVector> posVec2;
    
    int n = 1e+5;
    int nRepetitions =100;

    // int n;
    // cout<<"Give input vector size: ";
    // cin>>n;
    // int nRepetitions;
    // cout<<"Give nRepetitions: ";
    // cin>>nRepetitions;

    srand(time(NULL));
    //srand(2);
/*   vecgeom::Vector<ThreeVector> posVec1; 
    GenVecCart(posVec1, n);
    float Ts= TimeScalar(m1,posVec1,n,nRepetitions);
    float Tv= TimeVector(m1,posVec1,n,nRepetitions);
    cout<<"Vector speedup: " << Ts/ Tv <<endl;*/ 
    
    // Mean with different random numbers

    int noOfRuns = 16;
    float totVecTime= 0.;
    float totScaTime= 0.;
    vector<float> vecTimePerRepitition, scaTimePerRepitition; 
    ThreeVector sumFieldScalar; 
    ThreeVecSimd_t sumFieldVector;

    for (int i = 0; i < noOfRuns; ++i)
    {
        vecgeom::Vector<ThreeVector> posVec1;
        //srand(time(NULL));
        GenVecCart(posVec1, n);
        float Ts= TimeScalar(m1,posVec1,n,nRepetitions,sumFieldScalar);
        float Tv= TimeVector(m1,posVec1,n,nRepetitions,sumFieldVector);
        cout<<"Vector speedup: " << Ts/ Tv <<endl;

        vecTimePerRepitition.push_back(Tv);
        totVecTime += Tv;

        scaTimePerRepitition.push_back(Ts);
        totScaTime += Ts;

    }

    cout<<sumFieldScalar<<endl;
    cout<<sumFieldVector<<endl;

    float meanScalar, meanVector, stDevScalar, stDevVector;

    cout<<"\nScalar"<<endl;
    CalculateMeanStDev(scaTimePerRepitition, meanScalar, stDevScalar);
    cout<<"\nVector"<<endl;
    CalculateMeanStDev(vecTimePerRepitition, meanVector, stDevVector);
    cout<<"Mean speedup: "<<meanScalar/meanVector<<endl;


    



/*    GenVecCartContinuous(posVec2,n);

    cout<<"Size of posVec is: "<<posVec1.size()<<endl;
    // for (int i = 0; i < 100; ++i)
    // {
    //     cout<<i<<" : "<<posVec[i]<<endl;
    // }

    float Ts2= TimeScalar(m1,posVec2,n,nRepetitions);
    float Tv2= TimeVector(m1,posVec2,n,nRepetitions);

    cout<<"Vector speedup: " << Ts2/ Tv2 <<endl;  
*/
}


