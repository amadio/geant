//
//  First version:      (Josh) - GSoC 2014 project
//  Current version:  J. Apostolakis

#ifndef TUniformMagField_H
#define TUniformMagField_H

#include "GUVMagneticField.h"
#include <iostream>

#include "base/Vector3D.h"

#include "Constants.h"  //   For pi & twopi - Temporary solution ..

// using fieldConstants::pi;
// using fieldConstants::twopi;

class TUniformMagField : public GUVMagneticField
{
    public:  // with description

        TUniformMagField(const vecgeom::Vector3D<float>& FieldVector )
           : GUVMagneticField() //NumberOfComponents(3)
            // A field with value equal to FieldVector.
        {
           fFieldComponents = FieldVector;
        }

        TUniformMagField(double vField,
                         double vTheta,
                         double vPhi  );

        // virtual
        ~TUniformMagField() {}

        TUniformMagField(const TUniformMagField &p)   // : G4MagneticField(p)
        {
           fFieldComponents = p.fFieldComponents;
        }

        TUniformMagField& operator = (const TUniformMagField &p);

        // virtual
        void GetFieldValue( const vecgeom::Vector3D<double> &, // Position,
                                  vecgeom::Vector3D<float> &FieldValue )
        {
           FieldValue= fFieldComponents;
        }

        void SetFieldValue(const vecgeom::Vector3D<float>& fieldValue)
        {
           fFieldComponents= fieldValue;
        }

        vecgeom::Vector3D<float> GetConstantFieldValue() const
        {
           return fFieldComponents;
        }
        // Return the field value

        // virtual
        // TUniformMagField*
        GUVField* Clone() const
        {  return new TUniformMagField( *this ); }

        TUniformMagField* CloneOrSafeSelf( bool* pSafe )
        {
           if( pSafe ) *pSafe= true;
           return this; // ->CloneOrSafeSelf(*pSafe);
        }
        //  Class is thread-safe, can use 'self' instead of clone

    private:
        vecgeom::Vector3D<float> fFieldComponents;
};

TUniformMagField& 
TUniformMagField:: operator = (const TUniformMagField &p)
   // Copy constructor and assignment operator.
{
   if (&p == this) return *this;
   // for (int i=0; i<3; i++) fFieldComponents[i] = p.fFieldComponents[i];
   fFieldComponents = p.fFieldComponents;
   return *this;
}

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
   fFieldComponents.Set( vField*std::sin(vTheta)*std::cos(vPhi),
                         vField*std::sin(vTheta)*std::sin(vPhi),
                         vField*std::cos(vTheta)                );
}
#endif
