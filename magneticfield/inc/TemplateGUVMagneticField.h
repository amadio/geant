//===----------------------------------------------------------------------===//
/**
 * @file TemplateGUVMagneticField.h
 * @brief  Abstract field class for Geant-V prototype
 * @author Ananya
 *         Based on GUVMagneticField from John Apostolakis
 */
//===----------------------------------------------------------------------===//

#ifndef TemplateGUVMagneticField_H
#define TemplateGUVMagneticField_H

#include "TemplateGUVField.h"

template <class Backend>
class TemplateGUVMagneticField :  public TemplateGUVField<Backend> 
{

  typedef typename Backend::precision_v      Double_v;

  public:
    static constexpr int   fNumFieldComponents= 3;
    static constexpr bool  fFieldChangesEnergy= false;
  
    TemplateGUVMagneticField():  
     TemplateGUVField<Backend>( fNumFieldComponents, fFieldChangesEnergy) 
    {
      // std::cout<<"--- TemplateGUVMagneticField entered here ---"<<std::endl;
    }

    virtual ~TemplateGUVMagneticField(){}; 

    /*** void  GetFieldValue( const Double_v  Point[4],     // The old interface
     ***                            Double_v* Field );
     ***/
    
    virtual void GetFieldValue( const vecgeom::Vector3D<Double_v> &Position, 
                                      vecgeom::Vector3D<Double_v> &FieldValue ) = 0;

    TemplateGUVMagneticField& operator = (const TemplateGUVMagneticField &p);
    //  Copy 'standard' components ...
};

/***
template <class Backend>
void
TemplateGUVMagneticField<Backend>
  ::GetFieldValue( const typename Backend::precision_v  Point[4],     // The old interface
                         typename Backend::precision_v* FieldArr )
{
   typedef typename vecgeom::kVc::precision_v Double_v;

   vecgeom::Vector3D<Double_v> PositionV3D( Point[0], Point[1], Point[2]);
   vecgeom::Vector3D<Double_v>  Field_v3f;
   this->GetFieldValue( PositionV3D, Field_v3f );
   FieldArr[0]= Field_v3f.x();
   FieldArr[1]= Field_v3f.y();
   FieldArr[2]= Field_v3f.z();
}
 ***/
#endif
