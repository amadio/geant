#ifndef GUVVectorMagneticField_H
#define GUVVectorMagneticField_H

#include "GUVVectorField.h"
#include "GUVField.h"
#include "GUVMagneticField.h"

class GUVVectorMagneticField :  public GUVVectorField 
{
  typedef typename vecgeom::kVc::precision_v      Double_v;
//  typedef typename vecgeom::kVc::precision_v      Float_v;     // Was vecgeom::kVcFloat::precision_v 
  typedef typename vecgeom::kVcFloat::precision_v Float_v;
  
  public:
    static constexpr int   fNumFieldComponents= 3;
    static constexpr bool  fFieldChangesEnergy= false;
  
    GUVVectorMagneticField():  
     GUVVectorField( fNumFieldComponents, fFieldChangesEnergy) 
    {std::cout<<"--- GUVVectorMagneticField entered here ---"<<std::endl;}

    virtual ~GUVVectorMagneticField(){}; 

    void  GetFieldValue( const Double_v  Point[4],     // The old interface
                               Double_v* Field );
    
    virtual void GetFieldValue( const vecgeom::Vector3D<Double_v> &Position, 
                                      vecgeom::Vector3D<Float_v> &FieldValue ) = 0;

    GUVVectorMagneticField& operator = (const GUVVectorMagneticField &p);
    //  Copy 'standard' components ...
};

void
GUVVectorMagneticField::GetFieldValue( const typename vecgeom::kVc::precision_v  Point[4],     // The old interface
                                             typename vecgeom::kVc::precision_v* FieldArr )
{
   typedef typename vecgeom::kVc::precision_v Double_v;
   // typedef typename vecgeom::kVc::precision_v Float_v;     // Was  vecgeom::kVcFloat::precision_v;
   typedef typename vecgeom::kVcFloat::precision_v Float_v;
   
   vecgeom::Vector3D<Double_v> PositionV3D( Point[0], Point[1], Point[2]);
   vecgeom::Vector3D<Float_v>  Field_v3f;
   this->GetFieldValue( PositionV3D, Field_v3f );
   FieldArr[0]= (Double_v) Field_v3f.x();
   FieldArr[1]= (Double_v) Field_v3f.y();
   FieldArr[2]= (Double_v) Field_v3f.z();
}
#endif
