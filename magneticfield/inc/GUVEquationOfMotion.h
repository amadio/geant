//
// class GUVEquationOfMotion
//
// Class description:
//
// Abstract Base Class for the right hand size of the equation of
// motion of a particle in a field.

// History:
// - Created. J.Apostolakis     Dec 2014/Jan 2015
// -------------------------------------------------------------------

#ifndef GUV_EquationOfMotion_H
#define GUV_EquationOfMotion_H

#include <cassert>
#include <iostream>

// #include "GUVTypes.hh"      // "globals.hh"
#include "GUVField.h"   // required in inline method implementations

class GUVEquationOfMotion 
{
  public:  // with description

     GUVEquationOfMotion( GUVField *Field, unsigned short verbose=0 );
     virtual ~GUVEquationOfMotion();
       // Constructor and virtual destructor. No operations, just checks

     virtual void EvaluateRhsGivenB( const  double y[],
                                     const  double B[3],
                                     /*  double charge, */
                                            double dydx[] ) const = 0;
       // Given the value of the  field "B", this function 
       // calculates the value of the derivative dydx.
       // --------------------------------------------------------
       // This is the _only_ function a subclass must define.
       // The other two functions use Rhs_givenB.

      virtual void InitializeCharge(double particleCharge)=0;
       // Must be called to correctly initialise and provide charge
      virtual void InvalidateParameters()=0;
      
      inline void InformReady(); // All parameters have been set (charge+)
      inline void InformDone();  // Invalidate charge, other parameters
      inline void CheckInitialization() const; // Ensure initialization
      inline void CheckDone() const;

     // virtual void SetChargeMomentumMass(double particleCharge,
     //                                    double MomentumXc,
     //                                    double MassXc2) = 0;
     //   // Set the charge, momentum and mass of the current particle
     //   // --> used to set the equation's coefficients ...

      inline void RightHandSide( const  double y[],
                                     // double charge ,
                                        double dydx[] ) const;
       // This calculates the value of the derivative dydx at y.
       // It is the usual enquiry function.
       // ---------------------------
       // It uses the virtual function EvaluateRhsGivenB

     void EvaluateRhsReturnB( const  double y[],
                              double dydx[],
                           // double charge,
                              double Field[] ) const;
       // Same as RHS above, but also returns the value of B.
       // Should be made the new default ? after putting dydx & B in a class.

     void GetFieldValue( const double Point[4],
                               double Field[] )  const;
       // Obtain only the field - the stepper assumes it is pure Magnetic.
       // Not protected, because GUVRKG3_Stepper uses it directly.

     const GUVField* GetFieldObj() const {return fField;}
           GUVField* GetFieldObj()       {return fField;}
     void            SetFieldObj(GUVField* pField){fField=pField;}

     bool         Initialised() const { return fInitialised; } 
     unsigned int GetId() const       { return fEquationId; }
     static unsigned int GetNumCreated() { return fNumObjectsCreated; }
     static unsigned int GetNumLive() { return fNumObjectsCreated - fNumObjectsDeleted; }
       // For debugging, checking

     friend std::ostream&
             operator<<( std::ostream& os, const GUVEquationOfMotion& eq);

  public:
     static const unsigned int idxTime=3;  // Convention for location of time 't' in vector

  private:
     static unsigned int fNumObjectsCreated;
     static unsigned int fNumObjectsDeleted;
     // const int GUVmaximum_number_of_field_components = 24;
     enum { GUVmaximum_number_of_field_components = 24 } ;

     GUVField *     fField;
     unsigned int   fEquationId;  //
     unsigned short fVerbose;
     bool           fInitialised;
};

// #include "GUVEquationOfMotion.icc"

//  Inline implementation
//
// -------------------------------------------------------------------

inline
GUVEquationOfMotion::GUVEquationOfMotion(GUVField* pField, unsigned short verbose)
   : fField(pField), fEquationId(fNumObjectsCreated++),
     fVerbose(verbose), fInitialised(false)
{
   if( fVerbose )
      std::cout << " Created Equation " << this << " info= " << *this << std::endl;
}

inline
void GUVEquationOfMotion::InformReady() // was Initialize()
{
   // std::cout << " Called GUVEquationOfMotion::InformReady() " << std::endl;

   // assert( ! fInitialised ); // Sanity checking - assumes Clear() is always called!
                      // BUT: Will signal problem if two steppers share an equation
   fInitialised= true;
}

inline
void GUVEquationOfMotion::InformDone()  // was Clear() and before Finished();
{
   if(fVerbose)
   {
     std::cout << " Called GUVEquationOfMotion::InformDone() " << std::endl;
     std::cout << *this << std::endl;
   }
   assert( fInitialised );
   fInitialised= false;
}

inline
void GUVEquationOfMotion::GetFieldValue( const  double Point[4],
                             double Field[] ) const
{
    fField-> GetFieldValue( Point, Field );
}

inline
void
GUVEquationOfMotion::RightHandSide( const  double y[],
                                       //  double charge,
                                           double dydx[]  ) const
{
   CheckInitialization();

   double Field[GUVmaximum_number_of_field_components];
   double PositionAndTime[4];

     // Position
   PositionAndTime[0] = y[0];
   PositionAndTime[1] = y[1];
   PositionAndTime[2] = y[2];
   // Global Time
   PositionAndTime[3] = y[idxTime];  // See GUVFieldTrack::LoadFromArray

   GetFieldValue(PositionAndTime, Field) ;
   EvaluateRhsGivenB( y, Field, /*charge,*/ dydx );
}

#include <iostream>

void GUVEquationOfMotion::CheckInitialization() const
{
#ifdef GUVERBOSE
   if( fVerbose && !fInitialised ){
      std::cerr << "GUVEquationOfMotion is not Initialised" << std::endl;
   }
#endif
   assert( fInitialised );
}

void GUVEquationOfMotion::CheckDone() const
{
#ifdef GUVERBOSE
   if( fVerbose && fInitialised ){
      std::cerr << "GUVEquationOfMotion was NOT told it is Done!" << std::endl;
   }
#endif
   assert( !fInitialised );
}

#endif /* GUV_EquationOfMotion_DEF */
