//
// Embedded explicit Runge-Kutta Stepper using Dormand Prince's RK tableau
//

// First templated version:   J. Apostolakis, May 2018
//
//  DormandPrince7 - 5(4) embedded RK method
//
//  Class desription: 
//    An implementation of the 5th order embedded RK method from the paper
//    J. R. Dormand and P. J. Prince, “A family of embedded Runge-Kutta formulae,”
//	    Journal of computational and applied …, vol. 6, no. 1, pp. 19–26, 1980.
//
//  Design & Implementation by Somnath Banerjee
//  Supervision & code review: John Apostolakis
//
// Templated implementation adapted from 'CashKarp.h' by Ananya & J. Apostolakis

#ifndef __DefinedDormandPrince5RK_H
#define __DefinedDormandPrince5RK_H

#include "Geant/GUVectorLineSection.h"
// #include "VVectorIntegrationStepper.h"

// #include "AlignedBase.h"  // ==> Ensures alignment of storage for Vector objects

// #define Outside_DormandPrince5RK     1

template <class T_Equation, unsigned int Nvar>
class DormandPrince5RK
// : public VVectorIntegrationStepper
{
public:
  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

  using Double_v        = geant::Double_v;
  using ThreeVectorSimd = Vector3D<Double_v>;

  static constexpr unsigned int kOrderMethod = 4;
  static constexpr unsigned int sNstore // How many variables the full state entails
      = Nvar > 6 ? Nvar : 6;            // = std::max( 6, Nvar );
                                        // (GUIntegrationNms::NumVarBase > Nvar) ? GUIntegrationNms::NumVarBase : Nvar;
  // std::max( GUIntegrationNms::NumVarBase,  Nvar);
  // static const double IntegratorCorrection = 1./((1<<4)-1);
  inline int    GetIntegratorOrder() { return kOrderMethod; }
  inline double    IntegratorCorrection() { return 1. / ((1 << kOrderMethod) - 1); }

public:
  inline DormandPrince5RK(T_Equation *EqRhs, unsigned int numStateVariables = 0);

  DormandPrince5RK(const DormandPrince5RK &);

  virtual ~DormandPrince5RK();

  // GUVVectorIntegrationStepper* Clone() const;

  template <typename Real_v>
  struct ScratchSpaceDormandPrince5RK; // defined below

#ifdef Outside_DormandPrince5RK
  template <typename Real_v>
  // GEANT_FORCE_INLINE -- large method => do not force inline
  void StepWithErrorEstimate(const Real_v yInput[], // Consider __restrict__
                             const Real_v dydx[], const Real_v &charge, const Real_v &hStep, Real_v yOut[],
                             Real_v yErr[]
                             //, ScratchSpaceDormandPrince5RK<Real_v>* sp
                             ) const;
#endif

  //  ------Start of mandatory methods ( for transitional period. ) ------------
  //  To continue to inherit (for now) need to define:
  GEANT_FORCE_INLINE
  void StepWithErrorEstimate(const Double_v yInput[], // Consider __restrict__
                             const Double_v dydx[], const Double_v &charge, const Double_v &hStep, Double_v yOut[],
                             Double_v yErr[] ) const
  {
    StepWithErrorEstimate<Double_v>(yInput, dydx, charge, hStep, yOut, yErr);
  }

  Double_v DistChord() const { return Double_v(0.0); };
//  -------- End of mandatory methods ( for transitional period. ) ------------

#if ENABLE_CHORD_DIST
  template <typename Real_v>
  Real_v DistChord() const;
#endif

  template <typename Real_v>
  GEANT_FORCE_INLINE void RightHandSideInl(Real_v y[], const Real_v &charge, Real_v dydx[]) const
  {
    assert(fEquation_Rhs != nullptr);
    fEquation_Rhs->T_Equation::template RightHandSide<Real_v>(y, charge, dydx);
  }

  void SetEquationOfMotion(T_Equation *equation);

private:
  DormandPrince5RK &operator=(const DormandPrince5RK &) = delete;
  // private assignment operator.

public:
  template <typename Real_v>
  struct ScratchSpaceDormandPrince5RK {
    // State -- intermediate values used during RK step
    // -----
    Real_v ak2[sNstore];
    Real_v ak3[sNstore];
    Real_v ak4[sNstore];
    Real_v ak5[sNstore];
    Real_v ak6[sNstore];
    Real_v ak7[sNstore];
    Real_v ak8[sNstore];     
    Real_v yTemp2[sNstore]; // Separate temporaries per step - to aid compiler
    Real_v yTemp3[sNstore]; //   tradeoff benefit to be evaluated
    Real_v yTemp4[sNstore];
    Real_v yTemp5[sNstore];
    Real_v yTemp6[sNstore];
    // Real_v yTemp7[sNstore]; // Also the output as yOut = yTemp7
     
    Real_v yIn[sNstore];
// scratch space

#if ENABLE_CHORD_DIST
    // State -- values used ONLY for subsequent call to DistChord
    // -----
    Real_v fLastStepLength;
    Real_v fLastInitialVector[sNstore];
    Real_v fLastFinalVector[sNstore];
    Real_v fLastDyDx[sNstore];
    Real_v fMidVector[sNstore];
    Real_v fMidError[sNstore];
// for DistChord calculations
#endif
  public:
    ScratchSpaceDormandPrince5RK() {}
    ~ScratchSpaceDormandPrince5RK() {}
  };

  template <typename Real_v>
  ScratchSpaceDormandPrince5RK<Real_v> *ObtainScratchSpace()
  // Obtain object which can hold the scratch space for integration
  //   ( Should be re-used between calls - preferably long time
  {
    return new ScratchSpaceDormandPrince5RK<Real_v>();
  }

  // How to use it:
  //   auto = stepper->CreatedScratchSpace<Double_v>();

private:
  // 'Invariant' during integration - the pointers must not change
  // -----------
  T_Equation *fEquation_Rhs;
  bool fOwnTheEquation; //  --> indicates ownership of Equation object

  bool fDebug = false;

#ifdef Outside_DormandPrince5RK
};
#endif

// -------------------------------------------------------------------------------

#ifdef Outside_DormandPrince5RK
// template <class Real_v, class T_Equation, unsigned int Nvar>
template <class Real_v>
template <class T_Equation, unsigned int Nvar>
void DormandPrince5RK<T_Equation, Nvar>::
     StepWithErrorEstimate (   
        //  template StepWithErrorEstimate<Real_v> (
#else
public:
  template <typename Real_v>
  void           StepWithErrorEstimate (
#endif
        const Real_v yInput[],
        const Real_v dydx[], const Real_v &charge, const Real_v &Step, Real_v yOut[], Real_v yErr[]
        //, DormandPrince5RK<T_Equation,Nvar>::template ScratchSpaceDormandPrince5RK<Real_v>& sp
        ) const
{
  // const double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
  typename DormandPrince5RK<T_Equation, Nvar>::template ScratchSpaceDormandPrince5RK<Real_v> sp;

  // const double
  static constexpr double
    b21 = 0.2 ,
    
    b31 = 3.0/40.0, b32 = 9.0/40.0 ,
    
    b41 = 44.0/45.0, b42 = -56.0/15.0, b43 = 32.0/9.0,
    
    b51 = 19372.0/6561.0, b52 = -25360.0/2187.0, b53 = 64448.0/6561.0,
    b54 = -212.0/729.0 ,
    
    b61 = 9017.0/3168.0 , b62 =   -355.0/33.0,
    b63 =  46732.0/5247.0    , b64 = 49.0/176.0 ,
    b65 = -5103.0/18656.0 ,
    
    b71 = 35.0/384.0, b72 = 0.,
    b73 = 500.0/1113.0, b74 = 125.0/192.0,
    b75 = -2187.0/6784.0, b76 = 11.0/84.0;

  const double
    dc1 = -( b71 - 5179.0/57600.0),
    dc2 = -( b72 - .0),
    dc3 = -( b73 - 7571.0/16695.0),
    dc4 = -( b74 - 393.0/640.0),
    dc5 = -( b75 + 92097.0/339200.0),
    dc6 = -( b76 - 187.0/2100.0),
    dc7 = -( - 1.0/40.0 ); 
  
  // Initialise time to t0, needed when it is not updated by the integration.
  //       [ Note: Only for time dependent fields (usually electric)
  //                 is it neccessary to integrate the time.]
  // yOut[7] = yTemp[7]   = yIn[7];

  //  Saving yInput because yInput and yOut can be aliases for same array
  for (unsigned int i = 0; i < Nvar; i++) {
    sp.yIn[i] = yInput[i];
  }
  // RightHandSideInl(yIn, charge,  dydx) ;              // 1st Step

  for (unsigned int i = 0; i < Nvar; i++) {
    sp.yTemp2[i] = sp.yIn[i] + b21 * Step * dydx[i];
  }
  this->RightHandSideInl(sp.yTemp2, charge, sp.ak2); // 2nd Step

  for (unsigned int i = 0; i < Nvar; i++) {
    sp.yTemp3[i] = sp.yIn[i] + Step * (b31 * dydx[i] + b32 * sp.ak2[i]);
  }
  this->RightHandSideInl(sp.yTemp3, charge, sp.ak3); // 3rd Step

  for (unsigned int i = 0; i < Nvar; i++) {
    sp.yTemp4[i] = sp.yIn[i] + Step * (b41 * dydx[i] + b42 * sp.ak2[i] + b43 * sp.ak3[i]);
  }
  this->RightHandSideInl(sp.yTemp4, charge, sp.ak4); // 4th Step

  for (unsigned int i = 0; i < Nvar; i++) {
    sp.yTemp5[i] = sp.yIn[i] + Step * (b51 * dydx[i] + b52 * sp.ak2[i] + b53 * sp.ak3[i] + b54 * sp.ak4[i]);
  }
  this->RightHandSideInl(sp.yTemp5, charge, sp.ak5); // 5th Step

  for (unsigned int i = 0; i < Nvar; i++) {
    sp.yTemp6[i] =
        sp.yIn[i] + Step * (b61 * dydx[i] + b62 * sp.ak2[i] + b63 * sp.ak3[i] + b64 * sp.ak4[i] + b65 * sp.ak5[i]);
  }
  this->RightHandSideInl(sp.yTemp6, charge, sp.ak6); // 6th Step

  for(unsigned int i=0;i< Nvar;i++)
  {
     yOut[i] = yInput[i] + Step*(b71*dydx[i]   + b72*sp.ak2[i] + b73*sp.ak3[i] +
                                 b74*sp.ak4[i] + b75*sp.ak5[i] + b76*sp.ak6[i] );
  }
  this->RightHandSideInl(yOut, charge, sp.ak7); //7th and Final stage
  
  for (unsigned int i = 0; i < Nvar; i++) {
    // Estimate error as difference between 4th and 5th order methods
    //
    yErr[i] = Step * ( dc1 * dydx[i]   + dc2 * sp.ak2[i] + dc3 * sp.ak3[i] + dc4 * sp.ak4[i]
                     + dc5 * sp.ak5[i] + dc6 * sp.ak6[i] + dc7 * sp.ak7[i] );
    // std::cout<< "----In Stepper, yerrr is: "<<yErr[i]<<std::endl;
  }
#if ENABLE_CHORD_DIST
  for (unsigned int i = 0; i < Nvar; i++) {
    // Store Input and Final values, for possible use in calculating chord
    fLastInitialVector[i] = sp.yIn[i];
    fLastFinalVector[i]   = yOut[i];
    fLastDyDx[i]          = dydx[i];
  }
  fLastStepLength = Step;
#endif

  return;
}

#ifndef Outside_DormandPrince5RK
}
; // End of class declaration

//  The remaining functions / methods are defined below
#endif

// -------------------------------------------------------------------------------

template <class T_Equation, unsigned int Nvar>
inline DormandPrince5RK<T_Equation, Nvar>::DormandPrince5RK(T_Equation *EqRhs, unsigned int numStateVariables)
    : fEquation_Rhs(EqRhs),
      // fLastStepLength(0.),
      fOwnTheEquation(false)
{
  if (fDebug) {
    std::cout << "DormandPrince5RK constructor: Nvar is: " << Nvar << std::endl;
  }
// assert( dynamic_cast<TemplateVScalarEquationOfMotion<Backend>*>(EqRhs) != 0 );
#if ENABLE_CHORD_DIST
  fLastStepLength = Double_v(0.);
#endif
#ifndef GEANT_DEBUG
  (void)numStateVariables;
#endif
  assert((numStateVariables == 0) || (numStateVariables >= Nvar));
  assert(fEquation_Rhs != nullptr);
  // std::cout << "----end of constructor of DormandPrince5RK" << std::endl;
}

// -------------------------------------------------------------------------------

template <class T_Equation, unsigned int Nvar>
void DormandPrince5RK<T_Equation, Nvar>::SetEquationOfMotion(T_Equation *equation)
{
  fEquation_Rhs = equation;
  assert(fEquation_Rhs != nullptr);
}

// -------------------------------------------------------------------------------

//  Copy - Constructor
//
template <class T_Equation, unsigned int Nvar>
inline DormandPrince5RK<T_Equation, Nvar>::DormandPrince5RK(const DormandPrince5RK &right)
    : // fEquation_Rhs( (T_Equation*) nullptr ),
      fOwnTheEquation(false)
{
  if (fDebug) {
    std::cout << "----Entered *copy* constructor of DormandPrince5RK " << std::endl;
  }
  SetEquationOfMotion(new T_Equation(*(right.fEquation_Rhs)));
  assert(fEquation_Rhs != nullptr);
  // fEquation_Rhs= right.GetEquationOfMotion()->Clone());

  // assert( dynamic_cast<GUVVectorEquationOfMotion*>(fEquation_Rhs) != 0 );   // No longer Deriving
  assert(this->GetNumberOfStateVariables() >= Nvar);

#if ENABLE_CHORD_DIST
  fLastStepLength = Double_v(0.);
#endif

  if (fDebug)
    std::cout << " DormandPrince5RK - copy constructor: " << std::endl
              << " Nvar = " << Nvar << " Nstore= " << sNstore << " Own-the-Equation = " << fOwnTheEquation << std::endl;
}

// -------------------------------------------------------------------------------

template <class T_Equation, unsigned int Nvar>
GEANT_FORCE_INLINE DormandPrince5RK<T_Equation, Nvar>::~DormandPrince5RK()
{
  std::cout << "----- Vector DormandPrince5RK destructor" << std::endl;
  if (fOwnTheEquation)
    delete fEquation_Rhs; // Expect to own the equation, except if auxiliary (then sharing the equation)
  fEquation_Rhs = nullptr;
  std::cout << "----- VectorDormandPrince5RK destructor (ended)" << std::endl;
}

// -------------------------------------------------------------------------------

#ifdef Inheriting_DormandPrince5RK
template <class T_Equation, unsigned int Nvar>
GUVVectorIntegrationStepper *DormandPrince5RK<T_Equation, Nvar>::Clone() const
{
  // return new DormandPrince5RK( *this );
  return new DormandPrince5RK<T_Equation, Nvar>(*this);
}
#endif

// -------------------------------------------------------------------------------

#if ENABLE_CHORD_DIST
template <class Real_v, class T_Equation, unsigned int Nvar>
inline geant::Real_v DormandPrince5RK<T_Equation, Nvar>::DistChord() const
{
  Real_v distLine, distChord;
  ThreeVectorSimd initialPoint, finalPoint, midPoint;

  // Store last initial and final points (they will be overwritten in self-Stepper call!)
  initialPoint = ThreeVectorSimd(fLastInitialVector[0], fLastInitialVector[1], fLastInitialVector[2]);
  finalPoint   = ThreeVectorSimd(fLastFinalVector[0], fLastFinalVector[1], fLastFinalVector[2]);

  // Do half a step using StepNoErr
  fAuxStepper->StepWithErrorEstimate(fLastInitialVector, fLastDyDx, 0.5 * fLastStepLength, fMidVector, fMidError);

  midPoint = ThreeVectorSimd(fMidVector[0], fMidVector[1], fMidVector[2]);

  // Use stored values of Initial and Endpoint + new Midpoint to evaluate
  //  distance of Chord

  distChord = GUVectorLineSection::Distline(midPoint, initialPoint, finalPoint);
  return distChord;
}
#endif

// -------------------------------------------------------------------------------

#ifdef Outside_DormandPrince5RK
#undef Outside_DormandPrince5RK
#endif

#endif /*GUV Vector DormandPrince5RK */
