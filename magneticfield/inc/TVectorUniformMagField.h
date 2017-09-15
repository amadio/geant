//
//  First version:      (Josh) - GSoC 2014 project
//  Current version:  J. Apostolakis

#ifndef TVectorUniformMagField_H
#define TVectorUniformMagField_H

#include <iostream>
#include <base/Vector3D.h>
#include <Geant/VectorTypes.h>

#include "GUVVectorMagneticField.h"

#include "Constants.h"  //   For pi & twopi - Temporary solution ..

// using fieldConstants::pi;
// using fieldConstants::twopi;

class TVectorUniformMagField : public GUVVectorMagneticField
{
   using Double_v = Geant::Double_v;
   using Float_v = Geant::Float_v;
   
   template <typename T>
   using Vector3D = vecgeom::Vector3D<T>;

    public:  // with description
     
        TVectorUniformMagField(const Vector3D<float>& FieldVector )
           : GUVVectorMagneticField() //NumberOfComponents(3)
            // A field with value equal to FieldVector.
        {
           fFieldComponents = FieldVector;
        }

        TVectorUniformMagField(double vField,
                         double vTheta,
                         double vPhi  );

        // virtual
        ~TVectorUniformMagField() {}

        TVectorUniformMagField(const TVectorUniformMagField &p)   // : G4MagneticField(p)
        {
           fFieldComponents = p.fFieldComponents;
        }

        TVectorUniformMagField& operator = (const TVectorUniformMagField &p)
            // Copy constructor and assignment operator.
        {
           if (&p == this) return *this;
           // for (int i=0; i<3; i++) fFieldComponents[i] = p.fFieldComponents[i];
           fFieldComponents = p.fFieldComponents;
           return *this;
        }

        void GetFieldValue( Vector3D<double> const &, // Position,
                            Vector3D<float> &FieldValue )  override final
        {
           FieldValue= fFieldComponents;
        }

        void GetFieldValueSIMD( Vector3D<Double_v> const &, // Position,
                                Vector3D<Double_v> &FieldValue ) override final
        {
           // for (int i=0; i<3; i++) FieldValue[i] = Float_v( fFieldComponents[i] );
           FieldValue[0] = Double_v( fFieldComponents[0] );
           FieldValue[1] = Double_v( fFieldComponents[1] );
           FieldValue[2] = Double_v( fFieldComponents[2] );           
           // FieldValue= Vector3D<Float_v> (fFieldComponents);
        }

        void SetFieldValue(const Vector3D<float>& fieldValue)
        {
           fFieldComponents= fieldValue;
        }

        Vector3D<float> GetConstantFieldValue() const
        {
           return fFieldComponents;
        }
        // Return the field value

        // virtual
        TVectorUniformMagField* Clone() const override final
        {
           return new TVectorUniformMagField( *this );
        }

        TVectorUniformMagField* CloneOrSafeSelf( bool /*Safe = 0*/ )
        // {  Safe= true; return this; }  //  Class is thread-safe, can use 'self' instead of clone
        // { Safe= false; return new TVectorUniformMagField( this ); }  // Check ...
        { /*Safe= false;*/ return Clone(); }  // Check ...

        TVectorUniformMagField* CloneOrSafeSelf( bool* pSafe )
        {
           if( pSafe ) *pSafe= true;
           return this; // ->CloneOrSafeSelf(*pSafe);
        }
        //  Class is thread-safe, can use 'self' instead of clone

    private:
        Vector3D<float> fFieldComponents;
};

TVectorUniformMagField::TVectorUniformMagField(double vField,
                                   double vTheta,
                                   double vPhi     )
{
   if ( (vField<0) || (vTheta<0) || (vTheta>Constants::pi) || (vPhi<0) || (vPhi>Constants::twopi) )
   {
      // Exception("TVectorUniformMagField::TVectorUniformMagField()",
      //     "GeomField0002", FatalException, "Invalid parameters.") ;
      std::cerr << "ERROR in TVectorUniformMagField::TVectorUniformMagField()"
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
