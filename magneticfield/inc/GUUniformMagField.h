// 
// class GUUniformMagField
//
// Class description:
//
// Class for creation of Uniform Magnetic Field.

#ifndef GUUNIFORMMAGFIELD_HH
#define GUUNIFORMMAGFIELD_HH

#include "ThreeVector.h"
#include "GUField.h"

class GUUniformMagField : public GUField
{
  public:  // with description
  
    GUUniformMagField(const ThreeVector& FieldVector );
      // A field with value equal to FieldVector.

    GUUniformMagField(double vField,
                      double vTheta,
                      double vPhi     ) ;

    virtual ~GUUniformMagField() ;

    GUUniformMagField(const GUUniformMagField &p);
    GUUniformMagField& operator = (const GUUniformMagField &p);
      // Copy constructor and assignment operator.

    virtual void GetFieldValue(const double yTrack[4],
                                     double *MagField) const ;

    void SetFieldValue(const ThreeVector& newFieldValue);

    ThreeVector GetConstantFieldValue() const;
      // Return the field value
    
    virtual GUUniformMagField* Clone() const;

  private:

    double fFieldComponents[3] ;
};

#endif
///      Imlementation 

// #include "GUUniformMagField.h"
// #include "GUPhysicalConstants.h"

inline GUUniformMagField::GUUniformMagField(const ThreeVector& FieldVector )
{
      fFieldComponents[0] = FieldVector.x();
      fFieldComponents[1] = FieldVector.y();
      fFieldComponents[2] = FieldVector.z();
}

inline GUUniformMagField* GUUniformMagField::Clone() const
{
    return new GUUniformMagField( ThreeVector(this->fFieldComponents[0],
                                                this->fFieldComponents[1],
                                                this->fFieldComponents[2]) );
}

inline void
GUUniformMagField::SetFieldValue(const ThreeVector& newFieldVector )
{
      fFieldComponents[0] = newFieldVector.x();
      fFieldComponents[1] = newFieldVector.y();
      fFieldComponents[2] = newFieldVector.z();
}
   
inline GUUniformMagField::GUUniformMagField(double vField,
                                     double vTheta,
                                     double vPhi    )
{
   assert( vField >= 0); 
   assert( vTheta >= 0); 
   assert( vPhi >= 0); 
   if ( (vField<0) || (vTheta<0) || (vTheta>pi) || (vPhi<0) || (vPhi>twopi) )
   {
      // GUException("GUUniformMagField::GUUniformMagField()",
      //          "GeomField0002", FatalException, "Invalid parameters.") ;
      std::cerr << "GUUniformMagField::GUUniformMagField(): " 
                << "Fatal ERROR Invalid parameters."  ;
   }
   fFieldComponents[0] = vField*std::sin(vTheta)*std::cos(vPhi) ;
   fFieldComponents[1] = vField*std::sin(vTheta)*std::sin(vPhi) ;
   fFieldComponents[2] = vField*std::cos(vTheta) ;
}

GUUniformMagField::~GUUniformMagField()
{
}

GUUniformMagField::GUUniformMagField (const GUUniformMagField &p)
   : GUMagneticField(p)
{
   for (int i=0; i<3; i++)
      fFieldComponents[i] = p.fFieldComponents[i];
}

GUUniformMagField& GUUniformMagField::operator = (const GUUniformMagField &p)
{
   if (&p == this) return *this;
   for (int i=0; i<3; i++)
      fFieldComponents[i] = p.fFieldComponents[i];
   return *this;
}

// ------------------------------------------------------------------------

inline
void GUUniformMagField::GetFieldValue (const double [4],
                                             double *B  ) const 
{
   B[0]= fFieldComponents[0] ;
   B[1]= fFieldComponents[1] ;
   B[2]= fFieldComponents[2] ;
}

inline
ThreeVector GUUniformMagField::GetConstantFieldValue() const
{
   ThreeVector B(fFieldComponents[0],
                   fFieldComponents[1],
                   fFieldComponents[2]);
  return B;
}
