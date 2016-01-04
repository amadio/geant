//===----------------------------------------------------------------------===//
/**
 * @file GUVField.h
 * @brief  Abstract field class for Geant-V prototype
 * @author John Apostolakis
 */
//===----------------------------------------------------------------------===//

//
//
// class GUVField
//
// Class description:
//
// Abstract class for any kind of Field.
// It allows any kind of field (vector, scalar, tensor and any set of them)
// to be defined by implementing the inquiry function interface.
//
// The key method is  GetFieldValue( const  double Point[4],
//                    *************         double *fieldArr ) 
// Given an input position/time vector 'Point', 
// this method must return the value of the field in "fieldArr".
//
// A field must also specify whether it changes a track's energy:
//                    DoesFieldChangeEnergy() 
//                    *********************
// A field must co-work with a corresponding Equation of Motion, to
// enable the integration of a particle's position, momentum and, optionally, 
// spin.  For this a field and its equation of motion must follow the
// same convention for the order of field components in the array "fieldArr"
// -------------------------------------------------------------------

#ifndef GUVFIELD_HH
#define GUVFIELD_HH

#include <vector>
#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"
#include "backend/Backend.h"

// #include "GUVTypes.hh"
// #include "globals.hh"

/**
 * @brief Class GUVField
 */

class GUVField
{
  public:  // with description

      /**
       * @brief GeantTrack parametrized constructor
       *
       * @param Position - position (0,1,2=x,y,z)   [Input]   - Note: time is suppressed => B(t)=B(0)
       * @param fieldArr - output values of field. Usual convention:
       *                   0,1,2 = B_x, B_y, B_z
       *                   3,4,5 = E_x, E_y, E_z  (foreseen extension)
       *        Units are expected to be native GeantV units.
       */
      // virtual void  GetFieldValue( const double Position[4],
      //                                    double *fieldArr ) const = 0;
      virtual void GetFieldValue( const vecgeom::Vector3D<double> &Position,
                                        vecgeom::Vector3D<float>  &FieldValue ) = 0;

      /*
       * The expected vector interface is: 
       *
       * virtual void GetFieldValue( const vecgeom::Vector3D<typename Backend::precision_v> &Position, 
       *                                   vecgeom::Vector3D<typename Backend::precision_v> &FieldValue ) = 0;
       */

      inline
      GUVField( int NumberOfComponents, bool changesEnergy );
      inline
      GUVField( const GUVField & );
      virtual ~GUVField();
      // inline GUVField& operator = (const GUVField &p); 

     // A field signature function that can be used to insure
     // that the Equation of motion object and the GUVField object
     // have the same "field signature"?

      bool   DoesFieldChangeEnergy() const { return fChangesEnergy; } 
      int    GetNumberOfComponents() const { return fNumberOfComponents; } 

      GUVField& operator = (const GUVField &p); // Useful ?
      
      virtual GUVField* Clone() const = 0;
        // Implements cloning, likely needed for MT 

      // Expect a method of the following signature
      //  [Derived-Field-type] * CloneOrSafeSelf( bool* pSafe ) const
      // to be implemented for each derived class.
      // If the class is thread-safet, it can be implemented as:
      //  { if( pSafe ) { *pSafe= false; } ; return Clone(); } 
      
private:
      const int  fNumberOfComponents; // E.g.  B -> N=3 , ie x,y,z 
                                     //       E+B -> N=6 
      bool fChangesEnergy; 
       //  Each type/class of field set this accordingly:
       //    - an electric field     - "true"
       //    - a pure magnetic field - "false"
};

inline GUVField::GUVField( int NumberOfComponents, bool changesEnergy )
   : fNumberOfComponents(NumberOfComponents),
     fChangesEnergy(changesEnergy)
{
}

inline GUVField::GUVField( const GUVField &field )
   : fNumberOfComponents(field.fNumberOfComponents)
{
   // *this = field; 
   fChangesEnergy= field.fChangesEnergy;
}
#endif /* GUVFIELD_HH */
