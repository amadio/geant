#include "iostream"
#include <vector>
#include "stdlib.h"
#include <cmath>
#include <ctime>
//#include "/home/ananya/Softwares/iaca-lin64/include/iacaMarks.h" 

using namespace std;

int main(){


    //int x =rand()%10;
    double r = (double) rand()/(RAND_MAX) ;
    //cout<<r<<endl;
    //cout<< -220%200 <<endl;
    cout<< -221.1/200 << endl;
    cout<< -220/200 <<endl;
    cout<<25/50 <<endl;
    //vector<int> v1; v1={2, 34, 1};
    //cout<<v1[2]<<endl;
    double a = fmod(20.,3); //requires cmath library
    vector<int> v1;
    //IACA_START
    v1.push_back(10);
    v1.push_back(5);
    int b = v1.back();
    cout<<b<<endl;
    
    cout<<v1[0]<<" "<<v1[1]<<endl;
    //cout<< sqrt(4)<<endl;
    if(2<=3){cout<<"True"<<endl;}
    //IACA_END


    clock_t c1;
    cout<<c1<<endl;

    float fa = 2;
    double db =3 ;
    double dc = db-fa;

    // std::array<int,2> a={1,2};
    // std::array<int,2> b = a;

    int arr1[2] ={1,2};
    int arr2[2];
    copy(arr1,arr1+2,arr2);
}