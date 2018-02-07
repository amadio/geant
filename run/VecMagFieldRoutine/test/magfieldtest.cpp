#include <iostream>
#include <ctime>

#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"

#include <Vc/Vc>
#include "backend/vc/Backend.h"
#include "backend/vcfloat/Backend.h"
#include "backend/scalarfloat/Backend.h"
#include "backend/Backend.h"

// #include "MagField.h"
#include "CMSMagField.h"

using namespace std;

typedef vecgeom::Vector3D<float> Vec_t;
typedef vecgeom::Vector3D<vecgeom::kVcFloat::precision_v> ThreeVector;

typedef vecgeom::SOA3D<double> SOA3D;

int main(){
/*    
    MagField m1;
    m1.ReadVectorData("/home/ananya/Work/MagFieldRoutine/cms2015.txt");*/
/*
    Vector3D v1(2.,2.,2.);
    cout<<"vc Vector3D "<<v1<<endl;*/

/*    vecgeom::kVcFloat::precision_v v1;
    v1[0] = 2.;
    v1[1] = 3.;
    v1[2] = 4.,
    v1[3] = 5.;
    cout<<"Vc vector is: "<<v1<<endl;

    float f1=1, f2=1;
    Vec_t Pos2(f1,f1,f1);
    cout<<"Input vector is: "<<Pos2<<endl;
    Vec_t xyzField2;
    m1.GetFieldValue<vecgeom::kScalarFloat>(Pos2, xyzField2);

    cout<<"Vector First "<<endl;
    ThreeVector Pos(f2,f2,f2);
    ThreeVector Pos1;
    Pos1[0]= v1;
    cout<<"Pos1[0] is: "<<Pos1[0]<<endl;
    cout<<"Pos1 is: "<<Pos1<<endl;
    cout<<"first element : "<<Pos[0]<<endl;
    cout<<"Input vector is: "<<Pos<<endl;
    ThreeVector xyzField;
    m1.GetFieldValue<vecgeom::kVcFloat>(Pos, xyzField);
    cout<<"xyzField is: "<<xyzField<<endl;
*/


/*    cout<<"\n"<<endl;
    SOA3D a(10);
    cout<<"SOA3D "<<a[1]<<" " <<a[0]<<endl;
    cout<<"kSize is : "<<vecgeom::kVectorSize<<endl;

    vecgeom::kVc::precision_v A;
    vecgeom::kVc::precision_v B;
    for (int i = 0; i < 4; ++i)
    {
        A[i]=1;
        B[i]=2;
    }
    vecgeom::kVc::precision_v C;
    C = A+B;
    cout<<"Addition in C is: "<<C<<endl;
*/


    //Testing double template Magfield
    
    MagField<float> m1;
    m1.ReadVectorData("/home/ananya/Work/MagFieldRoutine/cms2015.txt");

    vecgeom::kVcFloat::precision_v v1;
    v1[0] = 2.;
    v1[1] = 3.;
    v1[2] = 4.,
    v1[3] = 5.;

    float f2 = 1.f;
    ThreeVector Pos(f2,f2,f2);
    ThreeVector Pos1;
    Pos1[0]= v1;
    ThreeVector xyzField;
    m1.GetFieldValue<vecgeom::kVcFloat>(Pos, xyzField);



    return 0;
}
