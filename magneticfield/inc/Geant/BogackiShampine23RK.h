//
// Embedded explicit Runge-Kutta Stepper using Cash Karp's RK tableau
//
// Adaptations of template interface: J. Apostolakis, Oct/Nov 2017
// First templated version:  Ananya, Feb/March 2016
//     ( commit 95e1316bcc156a04c876d6ea0fc9e60a15eeac4f )
//
// Adapted from 'GUTBogackiShampine23RK.hRKF45' by John Apostolakis, Nov 2015
//
// Adapted from 'GUTBogackiShampine23RK.hRKF45' by Qieshen Xie, GSoC 2014
//         (derived from G4BogackiShampine23RK.hRKF45)
//
//
#ifndef __BogackiShampine23RK_STEPPER_h
#define __BogackiShampine23RK_STEPPER_h

#include "Geant/GUVectorLineSection.h"
// #include "VVectorIntegrationStepper.h"

// #include "AlignedBase.h"  // ==> Ensures alignment of storage for Vector objects

// #define Outside_BogackiShampine23RK     1

template <class T_Equation, unsigned int Nvar>
class BogackiShampine23RK
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
  inline int GetIntegratorOrder() { return kOrderMethod; }
  inline double IntegratorCorrection() { return 1. / ((1 << kOrderMethod) - 1); }

public:
  inline BogackiShampine23RK(T_Equation *EqRhs, unsigned int numStateVariables = 0);

  BogackiShampine23RK(const BogackiShampine23RK &);

  virtual ~BogackiShampine23RK();

  // GUVVectorIntegrationStepper* Clone() const;

  template <typename Real_v>
  struct ScratchSpaceBogackiShampine23RK; // defined below

#ifdef Outside_BogackiShampine23RK
  template <typename Real_v>
  // GEANT_FORCE_INLINE -- large method => do not force inline
  void StepWithErrorEstimate(const Real_v yInput[], // Consider __restrict__
                             const Real_v dydx[], const Real_v &charge, const Real_v &hStep, Real_v yOut[],
                             Real_v yErr[]
                             //, ScratchSpaceBogackiShampine23RK<Real_v>* sp
                             ) const;
#endif

  //  ------Start of mandatory methods ( for transitional period. ) ------------
  //  To continue to inherit (for now) need to define:
  GEANT_FORCE_INLINE
  void StepWithErrorEstimate(const Double_v yInput[], // Consider __restrict__
                             const Double_v dydx[], const Double_v &charge, const Double_v &hStep, Double_v yOut[],
                             Double_v yErr[]) const
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
  BogackiShampine23RK &operator=(const BogackiShampine23RK &) = delete;
  // private assignment operator.

public:
  template <typename Real_v>
  struct ScratchSpaceBogackiShampine23RK {
    // State -- intermediate values used during RK step
    // -----
    Real_v ak2[sNstore];
    Real_v ak3[sNstore];
    Real_v ak4[sNstore];
    Real_v ak5[sNstore];
    Real_v ak6[sNstore];
    Real_v ak7[sNstore];
    Real_v yTemp2[sNstore]; // Separate temporaries per step - to aid compiler
    Real_v yTemp3[sNstore]; //   tradeoff benefit to be evaluated
    Real_v yTemp4[sNstore];
    Real_v yTemp5[sNstore];
    Real_v yTemp6[sNstore];

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
    ScratchSpaceBogackiShampine23RK() {}
    ~ScratchSpaceBogackiShampine23RK() {}
  };

  template <typename Real_v>
  ScratchSpaceBogackiShampine23RK<Real_v> *ObtainScratchSpace()
  // Obtain object which can hold the scratch space for integration
  //   ( Should be re-used between calls - preferably long time
  {
    return new ScratchSpaceBogackiShampine23RK<Real_v>();
  }

  // How to use it:
  //   auto = stepper->CreatedScratchSpace<Double_v>();

private:
  // 'Invariant' during integration - the pointers must not change
  // -----------
  T_Equation *fEquation_Rhs;
  bool fOwnTheEquation; //  --> indicates ownership of Equation object

  bool fDebug = false;

#ifdef Outside_BogackiShampine23RK
};
#endif

// -------------------------------------------------------------------------------

#ifdef Outside_BogackiShampine23RK
// template <class Real_v, class T_Equation, unsigned int Nvar>
template <class Real_v>
template <class T_Equation, unsigned int Nvar>
void BogackiShampine23RK<T_Equation, Nvar>::
    /*template*/ StepWithErrorEstimate /*<Real_v>*/ (
        const Real_v yInput[],
#else
public:
  template <typename Real_v>
  void StepWithErrorEstimate(const Real_v yInput[],
#endif

        const Real_v dydx[], const Real_v &charge, const Real_v &Step, Real_v yOut[], Real_v yErr[]
        //, BogackiShampine23RK<T_Equation,Nvar>::template ScratchSpaceBogackiShampine23RK<Real_v>& sp
        ) const
{
  // const double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
  typename BogackiShampine23RK<T_Equation, Nvar>::template ScratchSpaceBogackiShampine23RK<Real_v> sp;

  unsigned int i;

  const double  b21 = 0.5 ,
                b31 = 0. , b32 = 3.0/4.0 ,
                b41 = 2.0/9.0, b42 = 1.0/3.0 , b43 = 4.0/9.0;

  const double  dc1 = b41 - 7.0/24.0 ,  dc2 = b42 - 1.0/4.0 ,
                dc3 = b43 - 1.0/3.0 , dc4 = - 0.125 ;

  //  Saving yInput because yInput and yOut can be aliases for same array
  for (i = 0; i < Nvar; i++) {
    sp.yIn[i] = yInput[i];
  }
  // RightHandSideInl(yIn, charge,  dydx) ;          // 1st Stage

  for (i = 0; i < Nvar; i++) {
    sp.yTemp2[i] = sp.yIn[i] + b21 * Step * dydx[i];
  }
  this->RightHandSideInl(sp.yTemp2, charge, sp.ak2); // 2nd Stage

  for (i = 0; i < Nvar; i++) {
    sp.yTemp3[i] = sp.yIn[i] + Step * (b31 * dydx[i] + b32 * sp.ak2[i]);
  }
  this->RightHandSideInl(sp.yTemp3, charge, sp.ak3); // 3rd Stage

  for (i = 0; i < Nvar; i++) {
    // Accumulate increments with correct weights
    yOut[i] = sp.yIn[i] + Step * (b41 * dydx[i] + b42 * sp.ak2[i] + b43 * sp.ak3[i] );
  }
  this->RightHandSideInl(yOut, charge, sp.ak4);      // 4th Stage
  // Derivative and end-point already calculated in 'ak4' ! => Can be used in FSAL version
  
  for (i = 0; i < Nvar; i++) {
    // Estimate error as difference between 3rd and 2nd order methods
    //
    yErr[i] = Step * (dc1 * dydx[i] + dc2 * sp.ak2[i] +
                      dc3 * sp.ak3[i] + dc4 * sp.ak4[i]);
    // std::cout<< "----In Stepper, yerrr is: "<<yErr[i]<<std::endl;
  }
#if ENABLE_CHORD_DIST
  for (i = 0; i < Nvar; i++) {
    // Store Input and Final values, for possible use in calculating chord
    fLastInitialVector[i] = sp.yIn[i];
    fLastFinalVector[i]   = yOut[i];
    fLastDyDx[i]          = dydx[i];
  }
  fLastStepLength = Step;
#endif

  return;
}

#ifndef Outside_BogackiShampine23RK
}
; // End of class declaration

//  The remaining functions / methods are defined below
#endif

// -------------------------------------------------------------------------------

template <class T_Equation, unsigned int Nvar>
inline BogackiShampine23RK<T_Equation, Nvar>::BogackiShampine23RK(T_Equation *EqRhs, unsigned int numStateVariables)
    : fEquation_Rhs(EqRhs),
      // fLastStepLength(0.),
      fOwnTheEquation(false)
{
  if (fDebug) {
    std::cout << "\n----Entered constructor of BogackiShampine23RK " << std::endl;
    std::cout << "----In BogackiShampine23RK constructor, Nvar is: " << Nvar << std::endl;
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
  std::cout << "----end of constructor of BogackiShampine23RK" << std::endl;
}

// -------------------------------------------------------------------------------

template <class T_Equation, unsigned int Nvar>
void BogackiShampine23RK<T_Equation, Nvar>::SetEquationOfMotion(T_Equation *equation)
{
  std::cout << "Constructed BogackiShampine23RK stepper" << std::endl;
  fEquation_Rhs = equation;
  assert(fEquation_Rhs != nullptr);
}

// -------------------------------------------------------------------------------

//  Copy - Constructor
//
template <class T_Equation, unsigned int Nvar>
inline BogackiShampine23RK<T_Equation, Nvar>::BogackiShampine23RK(const BogackiShampine23RK &right)
    : // fEquation_Rhs( (T_Equation*) nullptr ),
      fOwnTheEquation(false)
{
  if (fDebug) {
    std::cout << "----Entered *copy* constructor of BogackiShampine23RK " << std::endl;
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
    std::cout << " BogackiShampine23RK - copy constructor: " << std::endl
              << " Nvar = " << Nvar << " Nstore= " << sNstore << " Own-the-Equation = " << fOwnTheEquation << std::endl;
}

// -------------------------------------------------------------------------------

template <class T_Equation, unsigned int Nvar>
GEANT_FORCE_INLINE BogackiShampine23RK<T_Equation, Nvar>::~BogackiShampine23RK()
{
  std::cout << "----- (Flexible) BogackiShampine23RK destructor" << std::endl;
  if (fOwnTheEquation)
    delete fEquation_Rhs; // Expect to own the equation, except if auxiliary (then sharing the equation)
  fEquation_Rhs = nullptr;
  std::cout << "----- (Flexible) BogackiShampine23RK destructor (ended)" << std::endl;
}

// -------------------------------------------------------------------------------

#ifdef Inheriting_BogackiShampine23RK
template <class T_Equation, unsigned int Nvar>
GUVVectorIntegrationStepper *BogackiShampine23RK<T_Equation, Nvar>::Clone() const
{
  // return new BogackiShampine23RK( *this );
  return new BogackiShampine23RK<T_Equation, Nvar>(*this);
}
#endif

// -------------------------------------------------------------------------------

#if ENABLE_CHORD_DIST
template <class Real_v, class T_Equation, unsigned int Nvar>
inline geant::Real_v BogackiShampine23RK<T_Equation, Nvar>::DistChord() const
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

#ifdef Outside_BogackiShampine23RK
#undef Outside_BogackiShampine23RK
#endif

#endif
