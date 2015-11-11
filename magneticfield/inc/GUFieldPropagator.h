//
//  Simple interface class to GUIntegrationDriver (with does Runge Kutta integration)
//   that follows the interface of ConstFieldHelixStepper.h
//
#include "ThreeVector.h"

class GUIntegrationDriver;
class GUVField;

class GUFieldPropagator
{
  public:
    // GUFieldPropagator(GUVField *); // First idea -- sidelined, at least for now
    GUFieldPropagator(GUIntegrationDriver* driver, double epsilon); // (GUVField* field)

    template<typename FieldType>  // , typename StepperType>
       GUFieldPropagator(FieldType* magField, double epsilon, double hminimum= 1.0e-4);

    virtual ~GUFieldPropagator() {}   //  Likely needed - to enable use of templated classes

     /**
       * Propagate track along in a field for length 'step'
       *    input: current position, current direction, particle properties
       *   output: success(returned), new position, new direction of particle
       */
  // GEANT_CUDA_BOTH_CODE          
      bool DoStep( ThreeVector const & position,  ThreeVector const & direction,
                           int const & charge,         double const & momentum,
                        double const & step,
                   ThreeVector       & endPosition,
                   ThreeVector       & endDiretion
         ) ;   //  Goal => make it 'const';  -- including all classes it uses

  /******
    template<typename Vector3D, typename DblType, typename IntType>
    inline
    __attribute__((always_inline))
    GEANT_CUDA_BOTH_CODE          
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
    GEANT_CUDA_BOTH_CODE
       void DoStep( DblType const & posx, DblType const & posy, DblType const & posz,
                    DblType const & dirx, DblType const & diry, DblType const & dirz,
                    IntType const & charge, DblType const & momentum, DblType const & step,
                    DblType & newsposx, DblType  & newposy, DblType  & newposz,
                    DblType & newdirx, DblType  & newdiry, DblType  & newdirz
                  ) const ;
   *****/

private:
    GUIntegrationDriver* fDriver;
    double               fEpsilon;
};

//  For Multi-threaded version only -- ie not for CUDA

#include <vector>

class GUFieldPropagatorPool
{
  public:
    // Access methods
    // static GUFieldPropagator* CreateOrFind(int numThreads);
      // It can be called from many threads -- same value must be returned
      //  numThreads must be constant between calls

    GUFieldPropagator* GetPropagator( int num );

  private:
    static GUFieldPropagatorPool* Instance();

    GUFieldPropagatorPool() {} // , void** banks=0 );  // Ensure one per thread
    ~GUFieldPropagatorPool() {} 

    // void Extend(int); 
  private:
    // static std::vector<GUFieldPropagator*> fFieldPropagatorVec;

};
