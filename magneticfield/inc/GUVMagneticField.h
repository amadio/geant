#ifndef GUVMagneticField_H
#define GUVMagneticField_H

#include "GUVField.h"

class GUVMagneticField : public GUVField

{
  public: 
    GUVMagneticField():  GUVField( 3, false) {}   // 3  components, Not change energy
      //   GUVField::fNumberOfComponents(3)     GUVField:fChangesEnergy(true),{}
    virtual ~GUVMagneticField(); // {}

    virtual void  GetFieldValue( const  double  Point[4],
                                       double* Field ) const = 0;
    // virtual GUVField* Clone() const;
    //   Concrete subclasses can (should?) implement it!

    GUVMagneticField& operator = (const GUVMagneticField &p);
    //  Copy 'standard' components ...
};

#endif
