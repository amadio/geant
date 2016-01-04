//
//  First version:      (Josh) - GSoC 2014 project
//  Current version:  J. Apostolakis

#ifndef TemplateTUniformMagField_H
#define TemplateTUniformMagField_H

#include "TemplateGUVMagneticField.h"
#include <iostream>

#include "base/Vector3D.h"

#include "Constants.h"  


template <class Backend>
class TemplateTUniformMagField : public TemplateGUVMagneticField<Backend>
{
  public:  // with description

    typedef typename Backend::precision_v Double_v;

    // Create a field, with value equal to FieldVector.
    TemplateTUniformMagField(const vecgeom::Vector3D<double>& FieldVector )
       : TemplateGUVMagneticField<Backend>()
    {
      fFieldComponents = FieldVector;
    }

    TemplateTUniformMagField(double vField, double vTheta, double vPhi, char );

    ~TemplateTUniformMagField() {}

    TemplateTUniformMagField(const TemplateTUniformMagField &p)   // : G4MagneticField(p)
    {
      fFieldComponents = p.fFieldComponents;
    }

    TemplateTUniformMagField& operator = (const TemplateTUniformMagField &p)
        // Copy constructor and assignment operator.
    {
      if (&p == this) return *this;
      // for (int i=0; i<3; i++) fFieldComponents[i] = p.fFieldComponents[i];
      fFieldComponents = p.fFieldComponents;
      return *this;
    }

/*        // virtual
    //  need to comment GetFieldValue, otherwise function overloading/repetition with Backend = kScalar
    void GetFieldValue( const vecgeom::Vector3D<double> &, // Position,
                              vecgeom::Vector3D<double> &FieldValue )
    {
       FieldValue= fFieldComponents;
    }*/

    void GetFieldValue( const vecgeom::Vector3D<Double_v> &, // Position,
                              vecgeom::Vector3D<Double_v> &FieldValue )
    {
      // for (int i=0; i<3; i++) FieldValue[i] = fFieldComponents[i];
      FieldValue[0] = fFieldComponents[0];
      FieldValue[1] = fFieldComponents[1];
      FieldValue[2] = fFieldComponents[2];
    }

    void SetFieldValue(const vecgeom::Vector3D<double>& fieldValue)
    {
      fFieldComponents= fieldValue;
    }

    vecgeom::Vector3D<double> GetConstantFieldValue() const
    {
      return fFieldComponents;
    }
    // Return the field value

    // virtual
    TemplateTUniformMagField* Clone() const
    {
      return new TemplateTUniformMagField( *this );
    }

    TemplateTUniformMagField* CloneOrSafeSelf( bool /*Safe = 0*/ )
    // {  Safe= true; return this; }  //  Class is thread-safe, can use 'self' instead of clone
    // { Safe= false; return new TemplateTUniformMagField( this ); }  // Check ...
    { /*Safe= false;*/ return Clone(); }  // Check ...

    TemplateTUniformMagField* CloneOrSafeSelf( bool* pSafe )
    {
      if( pSafe ) *pSafe= true;
      return this; // ->CloneOrSafeSelf(*pSafe);
    }
    //  Class is thread-safe, can use 'self' instead of clone

  private:
      vecgeom::Vector3D<double> fFieldComponents;
};


template <class Backend>
TemplateTUniformMagField<Backend>::TemplateTUniformMagField(double vField,
                                                            double vTheta,
                                                            double vPhi,
                                                            char   )
{
 if ( (vField<0) || (vTheta<0) || (vTheta>Constants::pi) || (vPhi<0) || (vPhi>Constants::twopi) )
 {
    // Exception("TemplateTUniformMagField::TemplateTUniformMagField()",
    //     "GeomField0002", FatalException, "Invalid parameters.") ;
    std::cerr << "ERROR in TemplateTUniformMagField::TemplateTUniformMagField()"
              << "Invalid parameter(s): expect " << std::endl;
    std::cerr << " - Theta angle: Value = " << vTheta
              << "  Expected between 0 <= theta <= pi = " << Constants::pi << std::endl;
    std::cerr << " - Phi   angle: Value = " << vPhi
              << "  Expected between 0 <=  phi  <= 2*pi = " << Constants::twopi << std::endl;
    std::cerr << " - Magnitude vField: Value = " << vField
              << "  Expected vField > 0 " << Constants::twopi << std::endl;
 }

 //std::sin and cos should work since vTheta etc are always scalar, but what the heck
 fFieldComponents.Set( vField*Backend::Sin(vTheta)*Backend::Cos(vPhi),
                       vField*Backend::Sin(vTheta)*Backend::Sin(vPhi),
                       vField*Backend::Cos(vTheta)                  );
}
#endif
