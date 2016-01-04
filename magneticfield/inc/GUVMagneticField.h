#ifndef GUVMagneticField_H
#define GUVMagneticField_H

#include "GUVField.h"

class GUVMagneticField : public GUVField

{
  public:
    static constexpr int   fNumFieldComponents= 3;
    static constexpr bool  fFieldChangesEnergy= false;
    GUVMagneticField():  GUVField( fNumFieldComponents, fFieldChangesEnergy) {}

    virtual ~GUVMagneticField(); // {}

    void  GetFieldValue( const double  Point[4],     // The old interface
                               double* Field );

    virtual void  GetFieldValue( const vecgeom::Vector3D<double> &Position, 
                                       vecgeom::Vector3D<float>  &FieldValue ) = 0;

    /*
     * The expected vector interface is: 
     *
     * virtual void GetFieldValue( const vecgeom::Vector3D<typename Backend::precision_v> &Position, 
     *                               vecgeom::Vector3D<typename Backend::precision_v> &FieldValue ) = 0;
     */

    // virtual GUVField* Clone() const;
    //   Concrete subclasses can (should?) implement it!

    GUVMagneticField& operator = (const GUVMagneticField &p);
    //  Copy 'standard' components ...
};

void
GUVMagneticField::GetFieldValue( const double  Point[4],     // The old interface
                                       double* FieldArr )
{
   vecgeom::Vector3D<double> PositionV3D( Point[0], Point[1], Point[2]);
   vecgeom::Vector3D<float>  Field_v3f;
   this->GetFieldValue( PositionV3D, Field_v3f );
   FieldArr[0]= Field_v3f.x();
   FieldArr[1]= Field_v3f.y();
   FieldArr[2]= Field_v3f.z();
}
#endif
