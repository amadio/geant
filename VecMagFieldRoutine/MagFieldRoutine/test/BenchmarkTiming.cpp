#include "iostream"
#include "MagField.h"
#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"
#include <string>
#include <vector>
#include <ctime>
#include <cmath> //for sqrt
#include <stdlib.h>

using namespace std;
typedef vecgeom::Vector3D<double> Vector3D;

const double kRMax=9000;
const double kZMax= 16000;

double RandR(){
        double r = (double) rand()/(RAND_MAX) ;
        r = r*kRMax; //because r is in range (0,9000) mm                                                                          
        return r;
}

double RandZ(){
        double z = (double) rand()/(RAND_MAX) ;
        z = z*kZMax; //range of z is between -16k and 16k                                                                         
        int sign = rand()%2; //to define the sign, since it can be both positive and negative                                     
        if (sign==0){
                z= -z;
        }
        return z;
}

void GenVecCartSubR(double &x, double &y){
        x = RandR();
        y = RandR();
        if((x*x + y*y)> kRMax*kRMax){
            GenVecCartSubR(x,y);
        }
}

void GenVecCart(Vector3D &pos){
        //cout<<"in genVecCart"<<endl;                                                                                            
        double x=0,y=0;
        double z = RandZ();
        GenVecCartSubR(x, y);
        pos.x()=x;
        pos.y()=y;
        pos.z()=z;
}

void GenVecCart(vector<Vector3D> &posVec, const int &n){
    for (int i = 0; i < n; ++i)
    {       
        Vector3D pos;
        GenVecCart(pos);
        posVec.push_back(pos);

    }
}



int main(){

    MagField m1;
    m1.ReadVectorData("VecMagFieldRoutine/cms2015.txt");
    vector<Vector3D> posVec, fieldVec;
    vector<double> testVec;
    int n = 1e+6;
    GenVecCart(posVec, n);

    clock_t c1= clock();
    for (int i = 0; i < n; ++i)
    {
        Vector3D xyzField;
        m1.GetFieldValueXYZ(posVec[i], xyzField);
        fieldVec.push_back(xyzField);
    }
    c1 = clock() - c1;

    float c = ((float)c1)/CLOCKS_PER_SEC;
    cout<<"Time for cartesian coordinates is : "<<c<<endl;
}


