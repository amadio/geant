//
//  Simple interface class to ScalarIntegrationDriver (with does Runge Kutta integration)
//   that follows the interface of ConstFieldHelixStepper.h
//
#ifndef GUFIELDPROPAGATOR_H
#define GUFIELDPROPAGATOR_H 1

#include "Geant/Config.h"
// #include "ThreeVector.h"
#include "base/Vector3D.h"
// typedef vecgeom::Vector3D<double>  ThreeVector; 
// #include "ScalarIntegrationDriver.h"

// namespace Geant {
// inline namespace GEANT_IMPL_NAMESPACE {

// template <class Backend> class TemplateGUIntegrationDriver;

class VVectorField;
class ScalarIntegrationDriver;
class FlexIntegrationDriver;

class GUFieldPropagator
{
  public:
    GUFieldPropagator(ScalarIntegrationDriver* scalarDriver,
                      double                   epsilon, 
                      FlexIntegrationDriver*   flexDriver = nullptr);

    // template <typename Backend>
    // GUFieldPropagator(TemplateGUIntegrationDriver<Backend>* driver, double epsilon);

    template<typename FieldType>  // , typename StepperType>
      GUFieldPropagator(FieldType* magField, double epsilon, double hminimum= 1.0e-4);

    ~GUFieldPropagator() {}  // Not virtual anymore.

    /**
      * Propagate track along in a field for length 'step'
      *    input: current position, current direction, particle properties
      *   output: success(returned), new position, new direction of particle
      */
  // VECCORE_ATT_HOST_DEVICE                 
    bool DoStep( vecgeom::Vector3D<double> const & position,  
                 vecgeom::Vector3D<double> const & direction,
                          int const & charge,         
                       double const & momentum,
                       double const & step,
                  vecgeom::Vector3D<double>      & endPosition,
                  vecgeom::Vector3D<double>      & endDiretion
        ) ;   //  Goal => make it 'const';  -- including all classes it uses

    ScalarIntegrationDriver* GetScalarIntegrationDriver(){ return fScalarDriver; }
    const ScalarIntegrationDriver* GetScalarIntegrationDriver() const { return fScalarDriver; }

    // static FlexIntegrationDriver* GetFlexibleIntegrationDriver(){ return fVectorDriver; }
    static const FlexIntegrationDriver* GetFlexibleIntegrationDriver() /*const*/ { return fVectorDriver; }
    static void SetFlexIntegrationDriver( FlexIntegrationDriver * flexDrv);
    
    double GetEpsilon() { return fEpsilon; }

    VVectorField* GetField();
    GUFieldPropagator* Clone() const; // { return this; }
         // Choice 1:  No longer allowing cloning !!  -- later solution
         // Choice 2:  Clone only the scalar 'old' stepper.  Share the flexible stepper 

  /******
    template<typename Vector3D, typename DblType, typename IntType>
    inline
    __attribute__((always_inline))
    VECCORE_ATT_HOST_DEVICE          
    template<typename Vector3D, typename DblType, typename IntType>
       void DoStep( Vector3D  const & pos,    Vector3D const & dir,
                    IntType   const & charge, DblType  const & momentum,
                    DblType   const & step,
                    Vector3D        & newpos,
                    Vector3D        & newdir
          ) const;

    
    //  Single value 
    // 
    template<typename DblType, typename IntType>
    inline
    __attribute__((always_inline))    
    VECCORE_ATT_HOST_DEVICE
       void DoStep( DblType const & posx, DblType const & posy, DblType const & posz,
                    DblType const & dirx, DblType const & diry, DblType const & dirz,
                    IntType const & charge, DblType const & momentum, DblType const & step,
                    DblType & newsposx, DblType  & newposy, DblType  & newposz,
                    DblType & newdirx, DblType  & newdiry, DblType  & newdirz
                  ) const ;
   *****/

private:
    ScalarIntegrationDriver*       fScalarDriver= nullptr;
    static FlexIntegrationDriver*  fVectorDriver;
    double                         fEpsilon;
};

// } // GEANT_IMPL_NAMESPACE
// } // Geant

#endif
