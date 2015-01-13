//
//  Simple interface class to GUIntegrationDriver (with does Runge Kutta integration)
//   that follows the interface of TGeoHelix  (for now)
//
#include "ThreeVector.h"

class GUIntegrationDriver;
class GUVField;

class GUFieldPropagator
{
  public: 
    GUFieldPropagator(GUVField *);
    virtual ~GUFieldPropagator() {}   //  Likely needed - to enable use of templated classes
    void SetCharge(double charge)  { fCharge= charge;} 
    void InitPoint(double x, double y, double z) { fInitialPosition= ThreeVector(x,y,z);}
    void InitDirection(double dx, double dy, double dz) { fInitialDirection= ThreeVector(dx,dy,dz);}

    void SetMomentum(double momentum) { fMomentumMag= momentum; } 
    
    // Auxiliary methods -- 
    void SetXYcurvature(double curvature) { fInitialCurvature= curvature; }

    //  methods -- names and actions to be reviewed
    void Step(double length);

    // Output methods
    const double *GetCurrentPoint()     { return fCurrentPoint; } 
    const double *GetCurrentDirection() { return fCurrentDirection; }

    // Null methods - needed to have same interface as Helix (for now)
    void UpdateHelix() {}
    void SetHelixStep(double length) { fStepLength= length; }

private:
    double      fCharge;
    ThreeVector fInitialPosition;
    ThreeVector fInitialDirection;
    double      fMomentumMag;          // Assume constant value (for now)
    double      fInitialCurvature;
    double      fStepLength;
    double      fCurrentPoint[3];
    double      fCurrentDirection[3];

    GUIntegrationDriver* fDriver;
};

//  For Multi-threaded version only -- ie not for CUDA

#include <vector>

class GUFieldPropagatorPool
{
  public:
    // Access methods
    static GUFieldPropagatorPool* CreateOrFind(int numThreads);
      // It can be called from many threads -- same value must be returned
      //  numThreads must be constant between calls

    GUFieldPropagator* GetPropagator( int num );

  private:
    static GUFieldPropagatorPool* Instance();

    GUFieldPropagatorPool() {} // , void** banks=0 );  // Ensure one per thread
    ~GUFieldPropagatorPool() {} 

  private:
    static std::vector<GUFieldPropagator*> fFieldPropagatorVec;

};
