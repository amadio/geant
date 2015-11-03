//  
//  First version:      (Josh) - GSoC 2014 project
//  Current version:  J. Apostolakis

#ifndef TUniformMagField_H
#define TUniformMagField_H

#include "GUVMagneticField.h"
#include <iostream>

#include "ThreeVector.h"  // Or whatever defines such a class
#include "Constants.h"  //   For pi & twopi - Temporary solution ..

// using fieldConstants::pi;
// using fieldConstants::twopi;

class TUniformMagField : public GUVMagneticField

{
    public:  // with description

        TUniformMagField(const ThreeVector& FieldVector )
           : GUVMagneticField() //NumberOfComponents(3)
            // A field with value equal to FieldVector.
        {
            fFieldComponents[0] = FieldVector.x();
            fFieldComponents[1] = FieldVector.y();
            fFieldComponents[2] = FieldVector.z();
        }

        TUniformMagField(double vField,
                         double vTheta,
                         double vPhi     );

        // virtual 
        ~TUniformMagField() {}

        TUniformMagField(const TUniformMagField &p)   // : G4MagneticField(p)
           
        {
            for (int i=0; i<3; i++)
                fFieldComponents[i] = p.fFieldComponents[i];
        }

        TUniformMagField& operator = (const TUniformMagField &p)
            // Copy constructor and assignment operator.
        {
            if (&p == this) return *this;
            for (int i=0; i<3; i++)
                fFieldComponents[i] = p.fFieldComponents[i];
            return *this;
        }
        
        // virtual
        void GetFieldValue(const double [4], // yTrack[4],
                                 double *B) const 
        {
            B[0]= fFieldComponents[0] ;
            B[1]= fFieldComponents[1] ;
            B[2]= fFieldComponents[2] ;
        }

        void SetFieldValue(const ThreeVector& fieldValue)
        {
            fFieldComponents[0] = fieldValue.x();
            fFieldComponents[1] = fieldValue.y();
            fFieldComponents[2] = fieldValue.z();
        }

        ThreeVector GetConstantFieldValue() const
        {
            ThreeVector B(fFieldComponents[0],
                    fFieldComponents[1],
                    fFieldComponents[2]);
            return B;
        }
        // Return the field value

        // virtual
        TUniformMagField* Clone() const
        { 
            return new TUniformMagField( ThreeVector(this->fFieldComponents[0],
                        this->fFieldComponents[1],
                        this->fFieldComponents[2]) );
        }

    private:
        double fFieldComponents[3];
};

TUniformMagField::TUniformMagField(double vField,
                                   double vTheta,
                                   double vPhi     ) 
{
   if ( (vField<0) || (vTheta<0) || (vTheta>Constants::pi) || (vPhi<0) || (vPhi>Constants::twopi) )
   {
      // Exception("TUniformMagField::TUniformMagField()",
      //     "GeomField0002", FatalException, "Invalid parameters.") ;
      std::cerr << "ERROR in TUniformMagField::TUniformMagField()"
                << "Invalid parameter(s): expect " << std::endl;
      std::cerr << " - Theta angle: Value = " << vTheta
                << "  Expected between 0 <= theta <= pi = " << Constants::pi << std::endl;
      std::cerr << " - Phi   angle: Value = " << vPhi
                << "  Expected between 0 <=  phi  <= 2*pi = " << Constants::twopi << std::endl;
      std::cerr << " - Magnitude vField: Value = " << vField
                << "  Expected vField > 0 " << Constants::twopi << std::endl;
   }
   fFieldComponents[0] = vField*std::sin(vTheta)*std::cos(vPhi) ;
   fFieldComponents[1] = vField*std::sin(vTheta)*std::sin(vPhi) ;
   fFieldComponents[2] = vField*std::cos(vTheta) ;
}
#endif
