#include "iostream"
#include "MagField.h"
#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"
#include <ctime>
using namespace std;

typedef vecgeom::Vector3D<double> Vector3D;


int main(){

    //cout<<"Print statements working"<<endl;
    
    MagField m1;
    //cout<<"here1"<<endl;
    m1.ReadVectorData("VecMagFieldRoutine/cms2015.txt");
    
    /*clock_t t, t2,t_; 
    int n = 1e+7; //size of randomR and Z vectors 
    t_= clock();
    m1.genVecR(n);
    m1.genVecZ(n);
    t = clock();
    m1.biLinearInterpolationRZ(m1.randomR, m1.randomZ);
    t= clock() - t;
    t2= clock() - t_ ;
    
    
    int n2= 1e+7; //size of vectors for cartesian coordinates
    m1.genVecCart(n2);
    clock_t c1= clock();
    m1.biLinearInterpolation(m1.random_x, m1.random_y, m1.random_z);
    c1 = clock() - c1;
    
    float t3 = ((float) t)/CLOCKS_PER_SEC;
    float t4 = ((float) t2)/CLOCKS_PER_SEC;
    float c = ((float)c1)/CLOCKS_PER_SEC;
    cout<<"Time taken(including time for generation of vectors) is "<<t4<<endl;
    cout<<"Time taken is "<<t3<<endl;
    cout<<"Time for cartesian coordinates is : "<<c<<endl;
    */
    /*  
    for (int i = 0; i < 50 ; ++i)  //or i < m1.predicted_Br.size()
    {
        //cout<<"Entering loop"<<endl;
        //cout<<m1.allR[i]<<" "<<m1.allZ[i]<<" "<<m1.predicted_BZ[i]<< " "<<m1.predicted_Br[i]<<" "<<m1.predicted_Bphi[i]<<endl;
        cout<<m1.all_x[i]<<" "<<m1.all_y[i]<<" "<<m1.all_z[i]<<" "<<m1.predicted_Bx[i]<<" "<<m1.predicted_By[i]<<" "<<m1.predicted_Bz[i]<<endl;
    }
    */
    

    Vector3D Pos ={1849,0,-8199}; 
    //or Vector3D Pos{1849,0,-8199}; 
    Vector3D rzField, xyzField;
    m1.GetFieldValueXYZ(Pos, xyzField);
    //cout<<Pos<<" "<<Field<<" "<<m1.predicted_BZ.back()<<" "<<m1.predicted_Br.back()<<" "<<m1.predicted_Bphi.back()<<endl;
    //m1.TestMagField(200, 400);
    return 0;
}
