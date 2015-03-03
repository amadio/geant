#ifndef GUVMagneticField_H
#define GUVMagneticField_H

#include "GUVField.h"

class GUVMagneticField : public GUVField

{
  public: 
    GUVMagneticField():  GUVField( 3, false) {}   // 3  components, Not change energy
       // GUVField:fChangesEnergy(true), GUVField::fNumberOfComponents(3)  {}
   virtual ~GUVMagneticField() {}

   virtual void  GetFieldValue( const  double  Point[4],
                                       double* Field ) const = 0;
   virtual GUVField* Clone() const;

   GUVMagneticField& operator = (const GUVMagneticField &p);
};

#endif
