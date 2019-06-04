//
// Runge-Kutta Stepper using Cash Karp's RK tableau
//
// Adapted from 'GUTCashKarpRKF45' by Qieshen Xie, GSoC 2014
//         (derived from G4CashKarpRKF45)
//
// First version:  John Apostolakis,  4 Nov 2015
//
#ifndef TDORMANDPRINCERK45_H
#define TDORMANDPRINCERK45_H

#include <iomanip> // For  C++ style output (debug)
#include <iostream>

#include "Geant/GULineSection.h"
#include "Geant/VScalarIntegrationStepper.h"
#include "base/Vector3D.h"

#define INLINERHS 1

#ifdef INLINERHS
#define REALLY_INLINE inline __attribute__((always_inline))
#else
#define REALLY_INLINE inline
#endif

template <class T_Equation, unsigned int Nvar>
class GUTDormandPrinceRK45 : public VScalarIntegrationStepper {
  using ThreeVector = vecgeom::Vector3D<double>;

public:
  static constexpr unsigned int kOrderMethod = 4;
  static constexpr unsigned int sNumVarBase  = 6; // Expected min number of Vars
  static constexpr unsigned int sNstore      = (sNumVarBase > Nvar) ? sNumVarBase : Nvar;
  // std::max( GUIntegrationNms::NumVarBase,  Nvar);
  // static const IntegratorCorrection = 1./((1<<4)-1);
  inline double IntegratorCorrection() { return 1. / ((1 << 4) - 1); }
  inline void SetVerbose(bool v) { fVerbose = v; }

public:
  inline GUTDormandPrinceRK45(T_Equation *EqRhs, unsigned int numStateVariables = 0, bool verbose = false);

  GUTDormandPrinceRK45(const GUTDormandPrinceRK45 &);

  virtual ~GUTDormandPrinceRK45();

  VScalarIntegrationStepper *Clone() const override;

  REALLY_INLINE
  void StepWithErrorEstimate(const double *yInput, // Consider __restrict__
                             const double *dydx, double charge, double Step, double *yOut, double *yErr) override final;

  double DistChord() const override;

  REALLY_INLINE
  void RightHandSideInl(const double y[], double charge, double dydx[]) const
  {
#ifndef VERBOSE_RHS
    fEquation_Rhs->T_Equation::RightHandSide(y, charge, dydx);
#else     
    using geant::units::tesla;
    std::cout << "Scalar GUT-DormandPrinceRK45::RHS called with q= " << charge
              << " at Position = " << y[0] << " y= " << y[1] << " z= " << y[2]
              << " with Momentum = " << y[3] << " y= " << y[4] << " z= " << y[5] << " ";
    vecgeom::Vector3D<double> Bfield;
    
    // fEquation_Rhs->T_Equation::RightHandSide(y, charge, dydx);
    fEquation_Rhs->T_Equation::EvaluateRhsReturnB(y, charge, dydx, Bfield);

    std::cout << " from B-field,  Bx= " << Bfield.x() / tesla << " By= " << Bfield.y() / tesla << " Bz= " << Bfield.z() / tesla << " ";
    std::cout << " gives Derivs dydx= :  x = " << dydx[0] << " y = " << dydx[1] << " z = " << dydx[2]
              << " px= " << dydx[3] << " py= " << dydx[4] << " pz= " << dydx[5]
              << std::endl;
#endif
  }

  REALLY_INLINE
  void RightHandSideInl(const double y[], double charge, double dydx[], vecgeom::Vector3D<double> &Bfield

                        )
  {
    // vecgeom::Vector3D<double> Position = { y[0], y[1], y[2] } ;
    // PositionTmp.Set( y[0], y[1], y[2] );
    // fEquation_Rhs->RightHandSide(y, /*PositionTmp,*/ dydx, charge, Bfield);

    //-- fEquation_Rhs->GetField()->T_Field::GetFieldValue(Point, Bfield);

    // fEquation_Rhs->T_Equation::RightHandSide(y, dydx, Bfield);
    // fEquation_Rhs->TEvaluateRhsGivenB( y, dydx, charge, Bfield);

    fEquation_Rhs->TEvaluateRhsReturnB(y, charge, dydx, Bfield);
  }

  void SetEquationOfMotion(T_Equation *equation);

  void PrintField(const char *label, const double y[6], const vecgeom::Vector3D<double> &Bfield) const;
  void PrintDyDx(const char *label, const double dydx[Nvar], const double y[Nvar]) const;
  void PrintDyDxLong(const char *label, const double dydx[Nvar], const double y[Nvar]) const;

private:
  GUTDormandPrinceRK45 &operator=(const GUTDormandPrinceRK45 &) = delete;
  // private assignment operator.

private:
  // 'Invariant' during integration - the pointers must not change
  // -----------
  T_Equation *fEquation_Rhs;
  bool fOwnTheEquation;
  GUTDormandPrinceRK45 *fAuxStepper;

  // State -- intermediate values used during RK step
  // -----
  double ak2[sNstore];
  double ak3[sNstore];
  double ak4[sNstore];
  double ak5[sNstore];
  double ak6[sNstore];
  double ak7[sNstore];
  double yTemp2[sNstore];
  double yTemp3[sNstore];
  double yTemp4[sNstore];
  double yTemp5[sNstore];
  double yTemp6[sNstore];
  double yIn[sNstore];
  vecgeom::Vector3D<double> Bfield2, Bfield3, Bfield4, Bfield5, Bfield6;
  vecgeom::Vector3D<double> PositionTmp;
  // scratch space

  // State -- values used for subsequent call to DistChord
  // -----
  double fLastStepLength;
  double *fLastInitialVector;
  double *fLastFinalVector;
  double *fInitialDyDx;
  /*volatile*/ double *fMidVector;
  /*volatile*/ double *fMidError;
  // for DistChord calculations

  // Parameters - for debugging etc
  bool fVerbose;
};

template <class T_Equation, unsigned int Nvar>
inline GUTDormandPrinceRK45<T_Equation, Nvar>::GUTDormandPrinceRK45(T_Equation *EqRhs,
                                                            // unsigned int noIntegrationVariables,
                                                            unsigned int numStateVariables,
                                                            bool verbose)
    : VScalarIntegrationStepper(EqRhs, // dynamic_cast<VScalarEquationOfMotion*>(EqRhs),
                                kOrderMethod, Nvar, ((numStateVariables > 0) ? numStateVariables : sNstore)),
      fEquation_Rhs(EqRhs), fOwnTheEquation(true), fAuxStepper(0), fLastStepLength(0.), fVerbose(verbose)
{
  assert(dynamic_cast<VScalarEquationOfMotion *>(EqRhs) != 0);
  assert((numStateVariables == 0) || (numStateVariables >= Nvar));
  assert(IntegratorOrder() == kOrderMethod);
  assert(GetNumberOfVariables() == Nvar);

  fLastInitialVector = new double[sNstore];
  fLastFinalVector   = new double[sNstore];
  fInitialDyDx          = new double[sNstore];

  fMidVector = new double[sNstore];
  fMidError  = new double[sNstore];
  if (verbose)
    std::cout << " GUTDormandPrinceRK45 - constructed class. " << std::endl
              << " Nvar = " << Nvar << " Nstore= " << sNstore << std::endl;
}

template <class T_Equation, unsigned int Nvar>
void GUTDormandPrinceRK45<T_Equation, Nvar>::SetEquationOfMotion(T_Equation *equation)
{
  fEquation_Rhs = equation;
  this->VScalarIntegrationStepper::SetEquationOfMotion(fEquation_Rhs);
}

//  Copy - Constructor
//
template <class T_Equation, unsigned int Nvar>
inline GUTDormandPrinceRK45<T_Equation, Nvar>::GUTDormandPrinceRK45(const GUTDormandPrinceRK45 &right)
    : VScalarIntegrationStepper((VScalarEquationOfMotion *)0, kOrderMethod, Nvar, right.GetNumberOfStateVariables()),
      fEquation_Rhs((T_Equation *)0), fOwnTheEquation(true), fAuxStepper(0), //  May overwrite below
      fLastStepLength(0.), fVerbose(right.fVerbose)
{
  SetEquationOfMotion(new T_Equation(*(right.fEquation_Rhs)));
  fOwnTheEquation = true;
  // fEquation_Rhs= right.GetEquationOfMotion()->Clone());

  assert(dynamic_cast<VScalarEquationOfMotion *>(fEquation_Rhs) != 0);
  assert(GetNumberOfStateVariables() >= Nvar);

  fLastInitialVector = new double[sNstore];
  fLastFinalVector   = new double[sNstore];
  fInitialDyDx          = new double[sNstore];

  fMidVector = new double[sNstore];
  fMidError  = new double[sNstore];

  if (fVerbose)
    std::cout << " GUTDormandPrinceRK45 - copy constructor: " << std::endl
              << " Nvar = " << Nvar << " Nstore= " << sNstore << " Own-the-Equation = " << fOwnTheEquation << std::endl;

  if (right.fAuxStepper) {
    // Reuse the Equation of motion in the Auxiliary Stepper
    fAuxStepper = new GUTDormandPrinceRK45(fEquation_Rhs, GetNumberOfStateVariables(), false);
  }
}

template <class T_Equation, unsigned int Nvar>
// inline
REALLY_INLINE GUTDormandPrinceRK45<T_Equation, Nvar>::~GUTDormandPrinceRK45()
{
  delete[] fLastInitialVector;
  delete[] fLastFinalVector;
  delete[] fInitialDyDx;
  delete[] fMidVector;
  delete[] fMidError;

  delete fAuxStepper;
  if (fOwnTheEquation)
    delete fEquation_Rhs; // Expect to own the equation, except if auxiliary (then sharing the equation)
}

template <class T_Equation, unsigned int Nvar>
VScalarIntegrationStepper *GUTDormandPrinceRK45<T_Equation, Nvar>::Clone() const
{
  // return new GUTDormandPrinceRK45( *this );
  return new GUTDormandPrinceRK45<T_Equation, Nvar>(*this);
}

template <class T_Equation, unsigned int Nvar>
inline void
   GUTDormandPrinceRK45<T_Equation, Nvar>::StepWithErrorEstimate(const double * yInput, 
                                                                 const double * dydx,   
                                                                 double   charge,
                                                                 double   Step,
                                                                 double * yOut, 
                                                                 double * yErr) 
{
  // const int nvar = 6 ;
  // const double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
  // std::cout << " Entered StepWithErrorEstimate of scalar " << std::endl;

  unsigned int i;

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

  static constexpr double
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
  for (i = 0; i < Nvar; i++) {
    yIn[i] = yInput[i];
  }
  // RightHandSideInl(yIn, charge,  dydx) ;              // 1st Step

  for (unsigned int i = 0; i < Nvar; i++) {
    yTemp2[i] = yIn[i] + b21 * Step * dydx[i];
  }
  this->RightHandSideInl(yTemp2, charge, ak2); // 2nd Step

  for (unsigned int i = 0; i < Nvar; i++) {
    yTemp3[i] = yIn[i] + Step * (b31 * dydx[i] + b32 * ak2[i]);
  }
  this->RightHandSideInl(yTemp3, charge, ak3); // 3rd Step

  for (unsigned int i = 0; i < Nvar; i++) {
    yTemp4[i] = yIn[i] + Step * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
  }
  this->RightHandSideInl(yTemp4, charge, ak4); // 4th Step

  for (unsigned int i = 0; i < Nvar; i++) {
    yTemp5[i] = yIn[i] + Step * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
  }
  this->RightHandSideInl(yTemp5, charge, ak5); // 5th Step

  for (unsigned int i = 0; i < Nvar; i++) {
    yTemp6[i] =
        yIn[i] + Step * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
  }
  this->RightHandSideInl(yTemp6, charge, ak6); // 6th Step

  for(unsigned int i=0;i< Nvar;i++)
  {
     yOut[i] = yInput[i] + Step*(b71*dydx[i]   + b72*ak2[i] + b73*ak3[i] +
                                 b74*ak4[i] + b75*ak5[i] + b76*ak6[i] );
  }
  this->RightHandSideInl(yOut, charge, ak7); //7th and Final stage
  
  for (unsigned int i = 0; i < Nvar; i++) {
    // Estimate error as difference between 4th and 5th order methods
    //
    yErr[i] = Step * ( dc1 * dydx[i]   + dc2 * ak2[i] + dc3 * ak3[i] + dc4 * ak4[i]
                     + dc5 * ak5[i] + dc6 * ak6[i] + dc7 * ak7[i] );
    // std::cout<< "----In Stepper, yerrr is: "<<yErr[i]<<std::endl;
  }
// #if ENABLE_CHORD_DIST
  for (unsigned int i = 0; i < Nvar; i++) {
    // Store Input and Final values, for possible use in calculating chord
    fLastInitialVector[i] = yIn[i];
    fLastFinalVector[i]   = yOut[i];
    fInitialDyDx[i]       = dydx[i];   // At initial point 
  }
// #endif
  fLastStepLength = Step;
  // std::cout << " Exiting StepWithErrorEstimate of scalar " << std::endl;

  return;
}

// #if ENABLE_CHORD_DIST
template <class T_Equation, unsigned int Nvar>
inline double GUTDormandPrinceRK45<T_Equation, Nvar>::DistChord() const
{
  // Coefficients were taken from Some Practical Runge-Kutta Formulas by Lawrence F. Shampine, page 149, c*
  static constexpr double hf1 = 6025192743.0 / 30085553152.0,
      hf2 = 0.0,
      hf3 = 51252292925.0 / 65400821598.0,
      hf4 = - 2691868925.0 / 45128329728.0,
      hf5 = 187940372067.0 / 1594534317056.0,
      hf6 = - 1776094331.0 / 19743644256.0,
      hf7 = 11237099.0 / 235043384.0;

  double midVector[3];

  for(int i = 0; i < 3; ++i) {
     midVector[i] = fLastInitialVector[i] + 0.5 * fLastStepLength * 
          (hf1 * fInitialDyDx[i] + hf2 * ak2[i] + hf3 * ak3[i] + 
           hf4 * ak4[i] + hf5 * ak5[i] + hf6 * ak6[i] + hf7 * ak7[i]);
  }
  double  distChord;
  ThreeVector initialPoint, finalPoint, midPoint;

  initialPoint = ThreeVector(fLastInitialVector[0], fLastInitialVector[1], fLastInitialVector[2]);
  finalPoint   = ThreeVector(fLastFinalVector[0], fLastFinalVector[1], fLastFinalVector[2]);
  midPoint     = ThreeVector(midVector[0], midVector[1], midVector[2]);

  // Use stored values of Initial and Endpoint + new Midpoint to evaluate
  //  distance of Chord
  distChord  = GULineSection::Distline(midPoint, initialPoint, finalPoint);
  
  return distChord;
}
// #endif

template <class T_Equation, unsigned int Nvar>
inline void GUTDormandPrinceRK45<T_Equation, Nvar>::PrintField(const char *label, const double y[Nvar],
                                                           const vecgeom::Vector3D<double> &Bfield) const
{
  std::cout << " PrintField/Stepper>  Field " << label << " "
            << "at x,y,z= ( " << y[0] << " , " << y[1] << " , " << y[2] << " ) "
            << " is ( " << Bfield.x() << " , " << Bfield.y() << " , " << Bfield.z()
            << " ) kGauss - mag = " << Bfield.Mag() << std::endl;
}

template <class T_Equation, unsigned int Nvar>
inline void GUTDormandPrinceRK45<T_Equation, Nvar>::PrintDyDx(const char *label, const double dydx[Nvar],
                                                          const double y[Nvar]) const
{
  using std::cout;

  if (fVerbose > 0) {
    vecgeom::Vector3D<double> dir(dydx[0], dydx[1], dydx[2]);
    vecgeom::Vector3D<double> dpds(dydx[3], dydx[4], dydx[5]);
    vecgeom::Vector3D<double> p(y[3], y[4], y[5]);
    int oldPrec = cout.precision(3);
    cout << " DyDx " << std::setw(4) << label << "> "
         << " xyz: " << std::setw(12) << dydx[0] << " , " << std::setw(12) << dydx[1] << " , " << std::setw(12)
         << dydx[2] << " ) "
         << " - mag = " << std::setw(12) << dir.Mag() << " dp/ds: " << std::setw(12) << dydx[3] << " , "
         << std::setw(12) << dydx[4] << " , " << std::setw(12) << dydx[5] << " "
         << " - mag = " << std::setw(5) << dpds.Mag() << " p-mag= " << p.Mag() << std::endl;
    cout.precision(oldPrec);
  }
}

template <class T_Equation, unsigned int Nvar>
inline void GUTDormandPrinceRK45<T_Equation, Nvar>::PrintDyDxLong(const char *label, const double dydx[Nvar],
                                                              const double y[Nvar]) const
{
  vecgeom::Vector3D<double> dir(dydx[0], dydx[1], dydx[2]);
  vecgeom::Vector3D<double> dpds(dydx[3], dydx[4], dydx[5]);
  std::cout << " PrintDyDx/Stepper>  dy/dx '" << std::setw(4) << label << "' "
            << " for x,y,z= ( " << dydx[0] << " , " << dydx[1] << " , " << dydx[2] << " ) "
            << " - mag = " << dir.Mag() << std::endl
            << "                              "
            << " dp/ds(x,y,z) = ( " << dydx[3] << " , " << dydx[4] << " , " << dydx[5] << " ) "
            << " ) - mag = " << dpds.Mag() << std::endl;
  vecgeom::Vector3D<double> p(y[3], y[4], y[5]);
  double pMag           = p.Mag();
  if (pMag == 0.0) pMag = 1.0;
  std::cout << "                                 "
            << " 1/p dp/ds = " << dydx[3] / pMag << " , " << dydx[4] / pMag << " , " << dydx[5] / pMag << " ) "
            << std::endl;
  std::cout << "                                 "
            << "         p = " << y[3] << " , " << y[4] << " , " << y[5] << " ) " << std::endl;
}
#endif /*TCashKARP_RKF45 */
