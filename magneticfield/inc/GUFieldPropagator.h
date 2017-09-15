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
// #include "TemplateGUIntegrationDriver.h"

// namespace Geant {
// inline namespace GEANT_IMPL_NAMESPACE {

template <class Backend>
class TemplateGUIntegrationDriver;

class ScalarIntegrationDriver;
class GUVField;

class GUFieldPropagator
{
  public:
    // GUFieldPropagator(GUVField *); // First idea -- sidelined, at least for now
    GUFieldPropagator(ScalarIntegrationDriver* driver, double epsilon); // (GUVField* field)

    template <typename Backend>
    GUFieldPropagator(TemplateGUIntegrationDriver<Backend>* driver, double epsilon);

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

    ScalarIntegrationDriver* GetIntegrationDriver(){ return fDriver; }
    const ScalarIntegrationDriver* GetIntegrationDriver() const { return fDriver; }
    double GetEpsilon() { return fEpsilon; }

    GUVField* GetField();
    GUFieldPropagator* Clone() const;

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
    ScalarIntegrationDriver* fDriver;
    double                  fEpsilon;
};

// } // GEANT_IMPL_NAMESPACE
// } // Geant

#endif
