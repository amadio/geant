//===--- (CMS)CMSmagField.h - Geant-V ------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file   (CMS)CMSmagField.h
 * @brief  Bi-linear interpolation of CMS field
 * @author Ananya
 */
//===----------------------------------------------------------------------===//

/*
 *  Details of current version / choices:
 *   - Reordered the way in which Gather was being used:
 *     Gathered elements for i1 and then i3 (i.e. the next r value)
 *     Then i2, then i4.
 *        The idea is to ensure that the different cache lines are accessed early, 
 *        so that they are available when the remaining values are needed, without
 *        further waiting.
 *   - Floats used 
 * 
 *  Note about ordering of memory:
 *         row  (stagerred in memory) 
 *          || 
 *          \/    Column ( consecutive in memory )
 *                  ===>
 *          i        i1       i2 (= i1 +1 )
 * 
 *         i+1       i3       i4
 */

#ifndef _CMSMAGFIELD_H__
#define _CMSMAGFIELD_H__

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cassert>
#include <ctime>

#include <base/Vector3D.h>
#include <base/Global.h>
#include <SystemOfUnits.h>
#include <Geant/VectorTypes.h>

//#include "Units.h"

#define FORCE_INLINE   1

// #include "VVectorField.h"

template<typename T>
struct MagVector3 {
  T fBr   =0.;
  T fBphi =0.;
  T fBz   =0.;
};

class CMSmagField // : public VVectorField
{
  using Double_v = Geant::Double_v;
  using Float_v = Geant::Float_v;
  
  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

public:
  CMSmagField();   
  CMSmagField(std::string inputMap);
  CMSmagField(const CMSmagField &right);

  /** @brief Scalar interface for field retrieval */
  // virtual
  void  GetFieldValue( const Vector3D<double> &position, 
                             Vector3D<double>  &fieldValue ) // override
  {
    GetFieldValue<double>(position, fieldValue);
  }

  /** @brief Vector interface for field retrieval */
  // virtual
  void GetFieldValueSIMD( const Vector3D<Double_v> &position, 
                                Vector3D<Double_v> &fieldValue ) // override
  {
    GetFieldValue<Double_v>(position, fieldValue);
  }

  /** @brief Templated field interface */
  template <typename Real_v>
  void GetFieldValue( const Vector3D<Real_v> &position,
                            Vector3D<Real_v> &fieldValue );

  // Reads data from specific 2D magnetic field map.
  //    ( Must be modified if used to read a different 2D map. )
  bool ReadVectorData(std::string inputMap);
   // Return value: success of finding and reading file.

  ~CMSmagField();

  void ReportVersion();

public: 
    //  Invariants -- parameters of the field 
    // static constexpr float millimeter = 0.1;             // Equal to Native GeantV unit
    // static constexpr float tesla = 10.0;                 // Navite unit = KiloGauss
    static constexpr float tesla     = geant::tesla; 
    static constexpr float kilogauss = geant::kilogauss;
    static constexpr float millimeter = geant::millimeter;
   
    const float kRMax   = 9000.  * millimeter;   //  Maximum value of R =  9.00 meters
    const float kZMax   = 16000. * millimeter;   //  Max value of Z = 16.00 meters
    static constexpr int kNoZValues   = 161;
    static constexpr int kNoRValues   = 181;
    static constexpr int kHalfZValues = 80;
    static constexpr int   gNumFieldComponents= 3;
    static constexpr bool  gFieldChangesEnergy= false;

    // Derived values
    //kRDiff and kZDiff take care of mm because they come from kRMax and kZMax which have mm in them
    const float kRDiff    = kRMax/(kNoRValues-1);   //  Radius increment between lattice points
    const float kZDiff    = 2*kZMax/(kNoZValues-1); //  Z increment

    const float kZ0       = -kZMax;
    const float kRDiffInv = 1.0/kRDiff;
    const float kZDiffInv = 1.0/kZDiff;
    const float kAInverse = tesla/(kRDiff*kZDiff);  // Values in files are Tesla

    //For (R,Z) pairs : gives field in cylindrical coordinates in rzfield
    template <typename Real_v>
    void GetFieldValueRZ(const    Real_v  &radius,
                         const    Real_v  &z, 
                         Vector3D<Real_v> &rzField); 

protected:
    //Used to convert cartesian coordinates to cylindrical coordinates R-Z-phi
    //Does not calculate phi
    template <typename Real_v>
    GEANT_FORCE_INLINE
    void CartesianToCylindrical(const Vector3D<Real_v> &cart, 
                                               Real_v   cyl[2]);

    //Converts cylindrical magnetic field to field in cartesian coordinates 
    template <typename Real_v>
    GEANT_FORCE_INLINE
    void CylindricalToCartesian(const Vector3D<Real_v> &rzField, 
                                const          Real_v  &sinTheta, 
                                const          Real_v  &cosTheta, 
                                      Vector3D<Real_v> &xyzField);

    //Gets the field array pointer in the appropriate form
    template <typename Real_v>
    GEANT_FORCE_INLINE
    const typename vecCore::Scalar<Real_v> *GetFieldArray() const;

    //Takes care of indexing into multiple places in AOS. 
    template <typename Real_v>
    GEANT_FORCE_INLINE
    void Gather2(const vecCore::Index<Real_v> index,
                       Real_v B1[3],
                       Real_v B2[3]);

public:
  // Methods for Multi-treading
  CMSmagField* CloneOrSafeSelf( bool* pSafe );
   // VVectorField*    Clone() const override;

  enum  kIndexRPhiZ { kNumR = 0, kNumPhi= 1, kNumZ = 2 } ;
private: 
   // MagVector3<float> *fMagvArray; //  = new MagVector3<float>[30000];
  float   fMagLinArray[3*kNoZValues*kNoRValues];
  double  fMagLinArrayD[3*kNoZValues*kNoRValues];
  bool    fReadData;
  bool    fVerbose;
  bool    fPrimary;  /** Read in and own the data arrays */
};

CMSmagField::CMSmagField()
  : // VVectorField(gNumFieldComponents, gFieldChangesEnergy),
    fReadData(false), fVerbose(true), fPrimary(false)
{
   // fMagvArray = new MagVector3<float>[kNoZValues*kNoRValues];
   if( fVerbose ) {   ReportVersion();  }
}

CMSmagField::CMSmagField(std::string inputMap) : CMSmagField() 
{
   if( fVerbose ) {   
     // ReportVersion();  
     std::cout<<"- CMSmagField c-tor #2" << std::endl;
   }
   // std::cout<<" Version: Reorder2 (floats) (with VC_NO_MEMBER_GATHER enabled if required)"<<std::endl;
   fReadData= CMSmagField::ReadVectorData(inputMap);
   if( fVerbose ) {
     std::cout<<"- CMSmagField c-tor #2: data has been read." << std::endl;
   }
   fPrimary= true;   // Own the data!
}

void CMSmagField::ReportVersion()
{
      printf( "\n%s", "CMSmagField class: Version: Reorder2 (floats)");
#ifdef VC_NO_MEMBER_GATHER
      printf( "%s", ", with VC_NO_MEMBER_GATHER enabled." );
#endif
}

CMSmagField::CMSmagField(const CMSmagField &right)
  : // VVectorField(gNumFieldComponents, gFieldChangesEnergy),
    fReadData(right.fReadData),
    fVerbose(right.fVerbose),
    fPrimary(false)
{
    // fMagvArray= right.fMagvArray;
}

CMSmagField::~CMSmagField()
{
   // if( fPrimary )
   //    delete[] fMagvArray;
}

bool CMSmagField::ReadVectorData(std::string inputMap)
{
   std::cout << "CMSmagField::ReadVectorData called with filename= " << inputMap << std::endl;
   std::string line;
   std::string s1,s2,s3,s4,s5,s0;
   float d1,d2,d3,d4,d5,d0;
   int ind =0;
   std::ifstream pFile(inputMap);
   if (pFile.is_open())
   {
      // getline() returns the stream. testing the stream with while returns error such as EOF
      while(getline(pFile,line)){
         // so here we know that the read was a success and that line has valid data
         std::stringstream ss(line);
         //parsing all the parts. s0's store the string names which are of no use to us. 
         ss>> s0>> d1>> s1>> d0>> s2>> d2>> s3>> d3>> s4>> d4>> s5>> d5;
      
         fMagLinArray[ind+kNumR]   = d4 * kAInverse;
         fMagLinArray[ind+kNumPhi] = d5 * kAInverse;
         fMagLinArray[ind+kNumZ]   = d3 * kAInverse;
#if    VERBOSE
         if( ind % 10 == 0 ) std::cout << "Read in line " << ind
                                       << " Values= " << d3 << " " << d4 << " "
                                       << d5 << std::endl;
#endif         
         ind += 3;
      }
      pFile.close();
      for (size_t i = 0; i < 3*kNoZValues*kNoRValues; ++i)
        fMagLinArrayD[i] = fMagLinArray[i];
   }
   else
   {
      std::cerr << "Unable to open file (for CMS mag field). Name = '" << inputMap
                << "'" << std::endl;
      exit(1);
   }
   return true;
}

template <typename Real_v>
GEANT_FORCE_INLINE
void CMSmagField::CartesianToCylindrical(const Vector3D<Real_v> &cart, 
                                               Real_v            cyl[2])
{

    //cyl[] =[r,z]
    cyl[0] = cart.Perp(); // calling sqrt at every call...
    cyl[1] = cart.z();
}

template <typename Real_v>
GEANT_FORCE_INLINE
void CMSmagField::CylindricalToCartesian(const Vector3D<Real_v>  &rzField, 
                                         const Real_v            &sinTheta, 
                                         const Real_v            &cosTheta, 
                                         Vector3D<Real_v>        &xyzField)
{
    //rzField[] has r, phi and z

    xyzField.x() = rzField.x()*cosTheta - rzField.y()*sinTheta; // Bx= Br cos(theta) - Bphi sin(theta)
    xyzField.y() = rzField.x()*sinTheta + rzField.y()*cosTheta; //By = Br sin(theta) + Bphi cos(theta)
    xyzField.z() = rzField.z();   //Bz = Bz 
}

template <typename Real_v>
GEANT_FORCE_INLINE
const typename vecCore::Scalar<Real_v> *CMSmagField::GetFieldArray() const
{
  return nullptr;
}

template <>
GEANT_FORCE_INLINE
const float *CMSmagField::GetFieldArray<Geant::Float_v>() const
{
  return fMagLinArray;
}

template <>
GEANT_FORCE_INLINE
const double *CMSmagField::GetFieldArray<Geant::Double_v>() const
{
  return fMagLinArrayD;
}

template <typename Real_v>
GEANT_FORCE_INLINE
void CMSmagField::Gather2(const vecCore::Index<Real_v> index, 
                                Real_v B1[3],
                                Real_v B2[3])
{
  using Index_v = vecCore::Index<Real_v>;
  using Real_s = vecCore::Scalar<Real_v>;
  
  Index_v ind1 = 3 * index;           // 3 components per 'location'
  Index_v ind2 = ind1   + 3 * kNoZValues;
  Real_s const *addr = GetFieldArray<Real_v>();
  //float const *addr = (Real_s const *)fMagLinArray;
  
  //Fetch one component of each point first, then the rest.
  B1[0] = vecCore::Gather<Real_v>( addr, ind1 );
  B2[0] = vecCore::Gather<Real_v>( addr, ind2 );

  Index_v ind1phi= ind1 + kNumPhi;
  Index_v ind2phi= ind2 + kNumPhi;
  B1[1] = vecCore::Gather<Real_v>( addr, ind1phi );
  B2[1] = vecCore::Gather<Real_v>( addr, ind2phi );
  
  Index_v ind1z = ind1 + kNumZ;
  Index_v ind2z = ind2 + kNumZ;
  B1[2] = vecCore::Gather<Real_v>( addr, ind1z );
  B2[2] = vecCore::Gather<Real_v>( addr, ind2z );
}

// Scalar specialization
template <>
GEANT_FORCE_INLINE
void CMSmagField::Gather2<double>(const vecCore::Index<double> index,
                                        double B1[3],
                                        double B2[3])
{
    int ind1= 3 * int(index);
    int ind2 = ind1 + 3 * kNoZValues;

    //Fetch one component of each point first, then the rest. 
    B1[0] = fMagLinArrayD[ind1+kNumR];
    B2[0] = fMagLinArrayD[ind2+kNumR];

    B1[1] = fMagLinArrayD[ind1+kNumPhi];
    B2[1] = fMagLinArrayD[ind2+kNumPhi];
    
    B1[2] = fMagLinArrayD[ind1+kNumZ];
    B2[2] = fMagLinArrayD[ind2+kNumZ];
}

template <typename Real_v>
GEANT_FORCE_INLINE
void CMSmagField::GetFieldValueRZ(const Real_v &r, 
                                  const Real_v &Z, 
                                  Vector3D<Real_v> &rzField)
{

    using namespace vecCore::math;
    using Index_v = vecCore::Index<Real_v>;

    //Take care that radius and z for out of limit values take values at end points 
    Real_v radius = Min(r, Real_v(kRMax));
    Real_v z = Max(Min(Z, Real_v(kZMax)), Real_v(-kZMax)); 

    //to make sense of the indices, consider any particular instance e.g. (25,-200)
    Real_v rFloor = Floor(radius * kRDiffInv);
    Real_v rIndLow = rFloor * kNoZValues;
    // Real_v rIndHigh = rIndLow + kNoZValues;

    //if we use z-z0 in place of two loops for Z<0 and Z>0
    //z-z0 = [0,32000]
    //so indices 0 to 160 : total 161 indices for (z-z0)/200
    //i.e. we are saying:
    Real_v zInd = Floor((z - Real_v(kZ0)) * Real_v(kZDiffInv));
    //need i1,i2,i3,i4 for 4 required indices
    Index_v i1 = Index_v(rIndLow + zInd);
    Index_v i2 = i1 + 1;

    Real_v zLow       = (zInd - Real_v(kHalfZValues)) * Real_v(kZDiff); //80 because it's the middle index in 0 to 160
    Real_v zHigh      = zLow + Real_v(kZDiff);
    Real_v radiusLow  = rFloor * Real_v(kRDiff);
    Real_v radiusHigh = radiusLow + Real_v(kRDiff);

    Real_v a1 = (radiusHigh - radius) * (zHigh - z); //area to be multiplied with i1
    Real_v a2 = (radiusHigh - radius) * (z - zLow);
    Real_v a3 = (radius - radiusLow) * (zHigh - z);
    Real_v a4 = (radius - radiusLow) * (z- zLow);

    
    Real_v B1[3], B2[3], B3[3], B4[3];
    Gather2<Real_v>(i1, B1, B3);
    Gather2<Real_v>(i2, B2, B4);

    Real_v BR   = B1[0]  *a1 + B2[0]  *a2 + B3[0]  *a3 + B4[0]  *a4; 
    Real_v BPhi = B1[1]  *a1 + B2[1]  *a2 + B3[1]  *a3 + B4[1]  *a4; 
    Real_v BZ   = B1[2]  *a1 + B2[2]  *a2 + B3[2]  *a3 + B4[2]  *a4; 

    rzField.x()= BR;
    rzField.y()= BPhi;
    rzField.z()= BZ;
}


template <typename Real_v>
GEANT_FORCE_INLINE
void CMSmagField::GetFieldValue(const Vector3D<Real_v>      &pos,
                                      Vector3D<Real_v> &xyzField)
{

//Sidenote: For theta =0; xyzField = rzField. 
//theta =0 corresponds to y=0
    Real_v cyl[2];
    CartesianToCylindrical<Real_v>(pos, cyl); 
    vecgeom::Vector3D<Real_v> rzField;
    GetFieldValueRZ<Real_v>(cyl[0], cyl[1], rzField); //cyl[2] =[r,z]

    using vecCore::Mask_v;

    Mask_v<Real_v> nonZero = (cyl[0] != Real_v(0.) );
    Real_v rInv      = vecCore::Blend(nonZero, Real_v(1.) / cyl[0],  Real_v(0.) );
    Real_v sinTheta  = pos.y() * rInv;
    Real_v cosTheta = vecCore::Blend(nonZero, pos.x() * rInv, Real_v(1.) );

    CylindricalToCartesian<Real_v>(rzField, sinTheta, cosTheta, xyzField);
}

// This class is thread safe.  So other threads can use the same instance
//
/*
CMSmagField* CMSmagField::CloneOrSafeSelf( bool* pSafe )
{
   if( pSafe ) *pSafe= true;
   return this;
}

VVectorField* CMSmagField::Clone() const
{
   return new CMSmagField( *this );
}
*/
#endif
