/*

Auto-vec. (Built on Cell-version)

*/


#ifndef _MAGFIELDAUTOVEC_H_
#define _MAGFIELDAUTOVEC_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"
#include <Vc/Vc>
#include "backend/vc/Backend.h"
#include "backend/vcfloat/Backend.h"
#include "backend/scalarfloat/Backend.h"
#include "backend/Backend.h"
#include <cassert>
#include <ctime>
using namespace std;


//typedef double dataType;
typedef float dataType;

#define FAST 

#ifdef FAST 
#define INLINE_CHOICE inline __attribute__ ((always_inline))
#endif
#ifndef FAST
#define INLINE_CHOICE __attribute__ ((noinline))
#endif

#ifndef _MAGFIELD_H_
#include "Structs.h"
#endif


class MagFieldAutoVec
{
public:
    MagFieldAutoVec();
    
    //New stuff 
    //Takes as input x,y,z; Gives output Bx,By,Bz
    template <class Backend>
    void GetFieldValue(const vecgeom::Vector3D<typename Backend::precision_v> &pos, 
                             vecgeom::Vector3D<typename Backend::precision_v> &xyzField);

    //Reads data from given 2D magnetic field map. Can be easily modified to read a given 2D map, in case the file changes
    void ReadVectorData(std::string inputMap);

    void StoreDataInMag4();

    void StoreDataInMag3SOA();
    
    ~MagFieldAutoVec();

private:
    //  Invariants -- parameters of the field 
    const dataType millimeter = 1.0;             // Currently -- to be native GeantV

    const dataType kRMax=  9000. * millimeter;   //  Maximum value of R =  9.00 meters
    const dataType kZMax= 16000. * millimeter;   //  Max value of Z = 16.00 meters
    const int kNoZValues = 161;
    const int kNoRValues = 181;
    const int kHalfZValues = 80;

    // Derived values
    //kRDiff and kZDiff take care of mm because they come from kRMax and kZMax which have mm in them
    const dataType kRDiff = kRMax/(kNoRValues-1);   //  Radius increment between lattice points
    const dataType kZDiff = 2*kZMax/(kNoZValues-1); //  Z increment

    const dataType kZ0 = -kZMax;
    const dataType kRDiffInv = 1.f/kRDiff;
    const dataType kZDiffInv = 1.f/kZDiff;
    const dataType kAInverse = 1/(kRDiff*kZDiff);

    //For (R,Z) pairs : gives field in cylindrical coordinates in rzfield
    template <class Backend>
    void GetFieldValueRZ(const typename Backend::precision_v &radius,
                         const typename Backend::precision_v      &z, 
                         vecgeom::Vector3D<typename Backend::precision_v> &rzField); 

    //Used to convert cartesian coordinates to cylindrical coordinates R-Z-phi
    //Does not calculate phi
    template <class Backend>
    void CartesianToCylindrical(const vecgeom::Vector3D<typename Backend::precision_v> &cart, 
                                                        typename Backend::precision_v cyl[2]);

    //Converts cylindrical magnetic field to field in cartesian coordinates 
    template <class Backend>
    void CylindricalToCartesian(const vecgeom::Vector3D<typename Backend::precision_v> &rzField, 
                                                  const typename Backend::precision_v  sinTheta, 
                                                  const typename Backend::precision_v  cosTheta, 
                                     vecgeom::Vector3D<typename Backend::precision_v> &xyzField);


    //Takes care of indexing into multiple places in AOS. Gather because using gather
    //defined in Vc class. Not self-defined gather like before 
    template <class Backend>
    void SandroGather(const typename Backend::precision_v   index, 
                            typename Backend::precision_v   Br[4],
                            typename Backend::precision_v Bphi[4],
                            typename Backend::precision_v   Bz[4]);


private: 
    std::vector<dataType> fRadius, fPhi, fZ, fBr, fBz, fBphi;
    MagVector      *fMagVector     = new MagVector     [30000];
    MagCellStructs *fMagCellStructs= new MagCellStructs[30000];
    MagCellArrays  *fMagCellArrays = new MagCellArrays [30000];

};



MagFieldAutoVec::MagFieldAutoVec(){
    std::cout<<"Auto-vec"<<std::endl;
}

MagFieldAutoVec::~MagFieldAutoVec(){
    
}


INLINE_CHOICE
void MagFieldAutoVec::ReadVectorData(std::string inputMap){
    std::string line;
    std::string s1,s2,s3,s4,s5,s0;
    dataType d1,d2,d3,d4,d5,d0;
    int ind =0;
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

            fMagVector[ind].SetBr(d4*kAInverse);
            fMagVector[ind].SetBphi(d5*kAInverse);
            fMagVector[ind].SetBz(d3*kAInverse);
            ind++;
        }
        pFile.close();

        StoreDataInMag4();
        StoreDataInMag3SOA();
    }
    else{
        cout<<"Unable to open file";
    }

}

INLINE_CHOICE
void MagFieldAutoVec::StoreDataInMag4(){
    for (int i = 0; i < (kNoRValues-1)*kNoZValues; ++i)
    {
        fMagCellStructs[i].m0 = fMagVector[i];
        fMagCellStructs[i].m1 = fMagVector[i+1];
        fMagCellStructs[i].m2 = fMagVector[i+kNoZValues];
        fMagCellStructs[i].m3 = fMagVector[i+kNoZValues+1];
    }
}

INLINE_CHOICE
void MagFieldAutoVec::StoreDataInMag3SOA(){
    for (int i = 0; i < (kNoRValues-1)*kNoZValues; ++i)
    {
        fMagCellArrays[i].sBr[0]   = fMagVector[i].Br;
        fMagCellArrays[i].sBphi[0] = fMagVector[i].Bphi;
        fMagCellArrays[i].sBz[0]   = fMagVector[i].Bz;

        fMagCellArrays[i].sBr[1]   = fMagVector[i+1].Br;
        fMagCellArrays[i].sBphi[1] = fMagVector[i+1].Bphi;
        fMagCellArrays[i].sBz[1]   = fMagVector[i+1].Bz;

        fMagCellArrays[i].sBr[2]   = fMagVector[i+kNoZValues].Br;
        fMagCellArrays[i].sBphi[2] = fMagVector[i+kNoZValues].Bphi;
        fMagCellArrays[i].sBz[2]   = fMagVector[i+kNoZValues].Bz;

        fMagCellArrays[i].sBr[3]   = fMagVector[i+kNoZValues+1].Br;
        fMagCellArrays[i].sBphi[3] = fMagVector[i+kNoZValues+1].Bphi;
        fMagCellArrays[i].sBz[3]   = fMagVector[i+kNoZValues+1].Bz;

    }
}


template <class Backend>
INLINE_CHOICE
void MagFieldAutoVec::CartesianToCylindrical(const vecgeom::Vector3D<typename Backend::precision_v> &cart, 
                                                              typename Backend::precision_v  cyl[2])
{

    //cyl[] =[r,z]
    cyl[0] = sqrt(cart.x()*cart.x() + cart.y()*cart.y()); // r = sqrt(x^2 + y^2)
    cyl[1] = cart.z(); //z = z 
}

template <class Backend>
INLINE_CHOICE
void MagFieldAutoVec::CylindricalToCartesian(const vecgeom::Vector3D<typename Backend::precision_v> &rzField, 
                                                        const typename Backend::precision_v  sinTheta, 
                                                        const typename Backend::precision_v  cosTheta, 
                                           vecgeom::Vector3D<typename Backend::precision_v> &xyzField)
{
    //rzField[] has r, phi and z

    xyzField.x() = rzField.x()*cosTheta - rzField.y()*sinTheta; // Bx= Br cos(theta) - Bphi sin(theta)
    xyzField.y() = rzField.x()*sinTheta + rzField.y()*cosTheta; // By = Br sin(theta) + Bphi cos(theta)
    xyzField.z() = rzField.z();   //Bz = Bz 
}



// Scalar Backend method 
template <class Backend>
INLINE_CHOICE 
void MagFieldAutoVec::SandroGather(const typename Backend::precision_v   index, 
                                  typename Backend::precision_v   Br[4],
                                  typename Backend::precision_v Bphi[4],
                                  typename Backend::precision_v   Bz[4])
{
    // //snippet 1
    //Better speed than snippet 2 (from vtune)
    int intIndex= (int) index;

    Br  [0] = fMagCellArrays[intIndex].sBr[0];
    Br  [1] = fMagCellArrays[intIndex].sBr[1];
    Br  [2] = fMagCellArrays[intIndex].sBr[2];
    Br  [3] = fMagCellArrays[intIndex].sBr[3];

    Bphi[0] = fMagCellArrays[intIndex].sBphi[0];
    Bphi[1] = fMagCellArrays[intIndex].sBphi[1];
    Bphi[2] = fMagCellArrays[intIndex].sBphi[2];
    Bphi[3] = fMagCellArrays[intIndex].sBphi[3];

    Bz  [0] = fMagCellArrays[intIndex].sBz[0];
    Bz  [1] = fMagCellArrays[intIndex].sBz[1];
    Bz  [2] = fMagCellArrays[intIndex].sBz[2];
    Bz  [3] = fMagCellArrays[intIndex].sBz[3];

    //snippet 2
 /*   Br[0] = fMagCellStructs[intIndex].m0.Br;
    Bphi[0] = fMagCellStructs[intIndex].m0.Bphi;
    Bz[0] = fMagCellStructs[intIndex].m0.Bz;

    Br[1] = fMagCellStructs[intIndex].m1.Br;
    Bphi[1] = fMagCellStructs[intIndex].m1.Bphi;
    Bz[1] = fMagCellStructs[intIndex].m1.Bz;

    Br[2] = fMagCellStructs[intIndex].m2.Br;
    Bphi[2] = fMagCellStructs[intIndex].m2.Bphi;
    Bz[2] = fMagCellStructs[intIndex].m2.Bz;

    Br[3] = fMagCellStructs[intIndex].m3.Br;
    Bphi[3] = fMagCellStructs[intIndex].m3.Bphi;
    Bz[3] = fMagCellStructs[intIndex].m3.Bz;
*/


}


// VcFloat Backend method 
template<>
INLINE_CHOICE
void MagFieldAutoVec::SandroGather<vecgeom::kVcFloat>(const typename vecgeom::kVcFloat::precision_v   index, 
                                                     typename vecgeom::kVcFloat::precision_v   Br[4],
                                                     typename vecgeom::kVcFloat::precision_v Bphi[4],
                                                     typename vecgeom::kVcFloat::precision_v   Bz[4])
{
    typedef typename vecgeom::kVcFloat::Int_t Int_v;
    Int_v ind = (Int_v) index; 

    Br  [0].gather(fMagCellStructs, &MagCellStructs::m0, &MagVector::Br  ,ind);
    Bphi[0].gather(fMagCellStructs, &MagCellStructs::m0, &MagVector::Bphi,ind);
    Bz  [0].gather(fMagCellStructs, &MagCellStructs::m0, &MagVector::Bz  ,ind);


    Br  [1].gather(fMagCellStructs, &MagCellStructs::m1, &MagVector::Br  ,ind);
    Bphi[1].gather(fMagCellStructs, &MagCellStructs::m1, &MagVector::Bphi,ind);
    Bz  [1].gather(fMagCellStructs, &MagCellStructs::m1, &MagVector::Bz  ,ind);


    Br  [2].gather(fMagCellStructs, &MagCellStructs::m2, &MagVector::Br  ,ind);
    Bphi[2].gather(fMagCellStructs, &MagCellStructs::m2, &MagVector::Bphi,ind);
    Bz  [2].gather(fMagCellStructs, &MagCellStructs::m2, &MagVector::Bz  ,ind);

    Br  [3].gather(fMagCellStructs, &MagCellStructs::m3, &MagVector::Br  ,ind);
    Bphi[3].gather(fMagCellStructs, &MagCellStructs::m3, &MagVector::Bphi,ind);
    Bz  [3].gather(fMagCellStructs, &MagCellStructs::m3, &MagVector::Bz  ,ind);

}



template <class Backend>
INLINE_CHOICE
//__attribute__((noinline))
void MagFieldAutoVec::GetFieldValueRZ(const typename Backend::precision_v &r, 
                               const typename Backend::precision_v &Z, 
                               vecgeom::Vector3D<typename Backend::precision_v> &rzField)
{

    typedef typename Backend::precision_v Float_v;

    //Take care that radius and z for out of limit values take values at end points 
    Float_v radius = min(r, kRMax);
    Float_v z      = max(min(Z, kZMax), -kZMax); 

    //to make sense of the indices, consider any particular instance e.g. (25,-200)
    Float_v rFloor   = floor(radius*kRDiffInv);
    Float_v rIndLow  = rFloor*kNoZValues;

    //if we use z-z0 in place of two loops for Z<0 and Z>0
    //z-z0 = [0,32000]
    //so indices 0 to 160 : total 161 indices for (z-z0)/200
    //i.e. we are saying:
    Float_v zInd = floor((z-kZ0)*kZDiffInv);
    //need i1,i2,i3,i4 for 4 required indices
    Float_v i1 = rIndLow + zInd;

    Float_v zLow       = (zInd- kHalfZValues)*kZDiff; //80 because it's the middle index in 0 to 160
    Float_v zHigh      = zLow + kZDiff;
    Float_v radiusLow  = rFloor*kRDiff;
    Float_v radiusHigh = radiusLow + kRDiff;
 
    //To enable auto vectorization
    Float_v area[4];
    area[0] = (radiusHigh - radius)   *(zHigh - z);
    area[1] = (radiusHigh - radius)   *(z - zLow );
    area[2] = (radius     - radiusLow)*(zHigh - z);
    area[3] = (radius     - radiusLow)*(z- zLow  );

    //We need B1[0], B2[0] ... 0s to be in one Br so that we have things aligned in memory.
    Float_v Br[4], Bphi[4], Bz[4];
    SandroGather<Backend>(i1, Br, Bphi, Bz);

    dataType zero = 0.f;

    Float_v BR(zero), BPhi(zero), BZ(zero);
    //cout<<BR<<endl;

    #pragma GCC ivdep
    //#pragma vector aligned 
    for (int i = 0; i < 4; ++i)
    {

       Br  [i] = Br[i]*area[i];
       Bphi[i] = Bphi[i]*area[i];
       Bz  [i] = Bz[i]*area[i];
    }

/*    #pragma GCC ivdep
    //#pragma vector aligned 
    for (int i = 0; i < 4; ++i)
    {
       BR   += Br  [i];
       BPhi += Bphi[i];
       BZ   += Bz  [i];
    }*/

    Br[0]   += Br[2];
    Bz[0]   += Bz[2];       
    Bphi[0] += Bphi[2];

    Br[1]   += Br[3];
    Bz[1]   += Bz[3];
    Bphi[1] += Bphi[3];

    BR  = Br[0]   + Br[1];
    BZ  = Bz[0]   + Bz[1];
    BPhi= Bphi[0] + Bphi[1]; 

    // Float_v BR   = (Br[0]  *area[0] + Br[1]  *area[1] + Br[2]  *area[2] + Br[3]  *area[3]);
    // Float_v BPhi = (Bphi[0]*area[0] + Bphi[1]*area[1] + Bphi[2]*area[2] + Bphi[3]*area[3]);
    // Float_v BZ   = (Bz[0]  *area[0] + Bz[1]  *area[1] + Bz[2]  *area[2] + Bz[3]  *area[3]);

    rzField.x()= BR;
    rzField.y()= BPhi;
    rzField.z()= BZ;
}


template <class Backend>
INLINE_CHOICE
//__attribute__ ((noinline))
//Sidenote: For theta =0; xyzField = rzField. 
//theta =0 corresponds to y=0

void MagFieldAutoVec::GetFieldValue(const vecgeom::Vector3D<typename Backend::precision_v>      &pos, 
                                   vecgeom::Vector3D<typename Backend::precision_v> &xyzField)
{

    typedef typename Backend::precision_v Float_v;
    typedef typename Backend::bool_v      Bool_v;

    Float_v cyl[2];
    CartesianToCylindrical<Backend>(pos, cyl); 
    vecgeom::Vector3D<Float_v> rzField;
    GetFieldValueRZ<Backend>(cyl[0], cyl[1], rzField); //cyl[2] =[r,z]
    
    //Did Float_v instead of dataType in order to generalize later 
    Float_v zero(0.);
    Float_v one(1.);
    //dataType zero = 0.;
    //dataType one = 1.;
    Float_v sinTheta(zero), cosTheta(one); //initialize as theta=0
    //To take care of r =0 case 

    //MaskedAssign(cond, value , var );
    //where cond is Bool_v, value is value calculated, var is the variable taking value 
    Bool_v  nonZero = (cyl[0] != zero); 
    Float_v rInv    = zero;
    vecgeom::MaskedAssign<dataType>(nonZero, 1/cyl[0]    , &rInv    );
    sinTheta = pos.y()*rInv;
    vecgeom::MaskedAssign<dataType>(nonZero, pos.x()*rInv, &cosTheta);

    CylindricalToCartesian<Backend>(rzField, sinTheta, cosTheta, xyzField);

}

#endif
