//
//
#include <iomanip>

#include "base/Vector3D.h"
#include <Geant/VectorTypes.h>

#include "Geant/PhysicalConstants.h"

// #include "GUVVectorEquationOfMotion.h"
// #include "TVectorMagFieldEquation.h"
#include "Geant/MagFieldEquation.h"
#include "Geant/ScalarUniformMagField.h"

#include "Geant/VVectorField.h"

#include "Geant/MagFieldEquation.h"

// #include "TMagFieldEquation.h"
// #include "Geant/ScalarFieldEquation.h"
#include "Geant/FieldEquationFactory.h"

// #include "Geant/ScalarUniformMagField.h"

#define CMS_FIELD 1

#ifdef CMS_FIELD
#include "Geant/Utils.h"
#include "CMSmagField.h"

using VectorFieldType = CMSmagField;
// using ScalarFieldType = CMSmagField;

#else
#include "Geant/UniformMagField.h"
#include "Geant/ScalarUniformMagField.h"

using VectorFieldType = UniformMagField;
// using ScalarFieldType = ScalarUniformMagField;
#endif

using Double_v = geant::Double_v;

using ThreeVector_f = vecgeom::Vector3D<float>;
using ThreeVector_d = vecgeom::Vector3D<double>;

using ThreeVector_dVec = vecgeom::Vector3D<Double_v>;

// constexpr unsigned int SixComp = 6; // Number of components in state of equation: position / momentum 
constexpr unsigned int gNposmom = 6; // Position 3-vec + Momentum 3-vec

// using ScalarEquationType = /* geant:: */ ScalarMagFieldEquation<ScalarFieldType, gNposmom>;
using VectorEquationType = /* geant:: */ MagFieldEquation<VectorFieldType>;

// using Float_v  = geant::Float_v;

using ThreeVector_DblV = vecgeom::Vector3D<Double_v>;
// using ThreeVector_FltV = vecgeom::Vector3D<Float_v>;

using std::cout;
using std::cerr;
using std::endl;
using std::setw;

VectorEquationType *CreateUniformFieldAndEquation(ThreeVector_f constFieldValue);
VectorEquationType *CreateFieldAndEquation(const char *filename);

bool TestEquation(VectorEquationType *);

const char *defaultFieldFileName = "cmsmagfield2015.txt";

int gVerbose = 0;  // Control global verbosity

int main(int, char **)
{
  ThreeVector_f fieldValue(0.0, 0.0, 1.0);
  bool bad = false;
  bool testCMSfield= false;

#ifdef CMS_FIELD
  testCMSfield= true;
  // VectorEquationType  or ScalarEquationType ?
  auto *eq2 = CreateFieldAndEquation(defaultFieldFileName); // ("cmsMagneticField2015.txt");
  bool errCMSfield               = TestEquation(eq2);
  bad = errCMSfield;
#else
  // VectorEquationType  or ScalarEquationType ?  
  auto *eq = CreateUniformFieldAndEquation(fieldValue);
  
  bool errUniform = TestEquation(eq);
  bad             = errUniform;
#endif

  if (bad) {
    cout << "*** TEST FAILED ***\n";
    cout << "    version tested:  " << ( testCMSfield ? "simplified CMS-field" : "Uniform" ) << endl;
  } else {
    // if( gVerbose ) 
    cout << "*** TEST PASSED ***\n";
  }

  return bad;
}

#if !defined(CMS_FIELD)

VectorEquationType *CreateUniformFieldAndEquation(ThreeVector_f fieldValue)
{
  // using Field_t = ScalarUniformMagField;

  auto *pConstBfield = // new ScalarUniformMagField(fieldValue);
      new UniformMagField(fieldValue);

  // 1. Original way of creating an equation
  // using VectorEquationType = ScalarMagFieldEquation<ScalarUniformMagField, gNposmom>;
                  // Before: MagFieldEquation<ScalarUniformMagField, gNposmom>;
  auto magEquation   = new VectorEquationType(pConstBfield);
  return magEquation;

  //  2. Different method of creating equation:  Factory
  // auto vecEquation = FieldEquationFactory::CreateMagEquation<ScalarUniformMagField>(pConstBfield);
  // return vecEquation;
}

#else

VectorEquationType *CreateFieldAndEquation(const char *filename)
{
  // const char *defaultFieldFileName= "cmsMagneticField2015.txt";

  //  3. Equation for CMS field
  auto cmsField    = new CMSmagField(filename ? filename : defaultFieldFileName);
  // auto equationCMS = FieldEquationFactory::CreateMagEquation<CMSmagField>(cmsField);
  auto    equationCMS = new VectorEquationType(cmsField);
  
  return equationCMS;
}

#endif

// Auxiliary methods (define below) to check results & report problems.
//
bool SanityCheckMagFieldEquation(double charge, ThreeVector_d Momentum, ThreeVector_d Field, ThreeVector_d ForceVec,
                                 int lane);

bool CheckDerivativeInLanesAndReport(const Double_v         & chargeVec,
                                     const Double_v           PositionMomentum[gNposmom],
                                     const ThreeVector_DblV /*vecgeom::Vector3D<Double_v>*/ & FieldVec,
                                     const Double_v           dydxVec[gNposmom],
                                     bool                     printContents = false);

bool TestEquation(VectorEquationType *equation)
{

  bool hasError = false; // Return value

  ThreeVector_d initialPosition3d(1., 2., 3.); // initial
  ThreeVector_d Momentum(0., 0.1, 1.);
  ThreeVector_d BFieldValue(0., 0., 1.); // Magnetic field value (constant)

  // double PositionTime[4] = { Position.x(), Position.y(), Position.z(), 0.0};
  Double_v PositionMomentum[gNposmom];
  Double_v dydxVec[gNposmom];
  int chargeVecScalar[8] = {-1, 1, 2, -2, -3, -3, -3, -3};

  // Revise the values, so that they are no longer equal
  Double_v chargeVec = Double_v(0.0); // { -1.0, 1.0, 2.0, -2.0 } ;
  for (size_t lane = 0; lane < geant::kVecLenD; ++lane)
    vecCore::Set(chargeVec, lane, (double)chargeVecScalar[lane % 8]);

  vecgeom::Vector3D<Double_v> FieldVec = {Double_v(BFieldValue[0]), Double_v(BFieldValue[1]), Double_v(BFieldValue[2])};
  for (int i = 0; i < 3; i++) {
    PositionMomentum[i]     = Double_v(initialPosition3d[i]);
    PositionMomentum[3 + i] = Double_v(Momentum[i]);
  }

  // Input check
  //
  const bool printInput = false;

  if (printInput) {
    for (int i = 0; i < 3; i++)
      cout << " pos[" << i << "] = " << setw(6) << initialPosition3d[i] << " PositionMomentum[] = " << PositionMomentum[i]
           << endl;
    for (int i = 0; i < 3; i++)
      cout << " mom[" << i << "] = " << setw(6) << Momentum[i] << " PositionMomentum[] = " << PositionMomentum[3 + i]
           << endl;

    cout << "Charge Vec = " << chargeVec << "  expected  -1, 1, -2, 2. " << endl;
    cout << "Field Vec = " << FieldVec << endl;
  }

  bool printContents = false; // Verbose output of  dy/dx, y, etc.

  // 1.) Simple case: Use the equation with a given field value
  // ------------------------------------------------------------
  if( gVerbose )
     cout << "Checking with fixed values of Field Vec = { 0, 1, 2 }. " << endl;
  
  equation->EvaluateRhsGivenB(PositionMomentum, chargeVec, FieldVec, dydxVec);

  //   Simple printing for visual cross check of output
  if (printContents) {
    cout << endl;
    cout << " ============================================ " << endl;
    cout << " Output - dy/dx Vec " << endl;
    cout << " dy/dx Vec [3] : " << dydxVec[3] << endl;
    cout << " dy/dx Vec [4] : " << dydxVec[4] << endl;
    cout << " dy/dx Vec [5] : " << dydxVec[5] << endl;
    cout << " ============================================ " << endl;
  }

  // a.) Test each output individually first
  // ======================================
  bool laneError1 = CheckDerivativeInLanesAndReport(chargeVec, PositionMomentum, FieldVec, dydxVec, printContents);

  hasError = laneError1;

  // b.) Check lanes against each other ( know relations from charge, momentum etc. )
  //  TO DO

  // 2.) Regular case: Use the full equation -- allow it to obtain the field too !
  // -------------------------------------------------------------------------
  cout << "Checking using RightHandSide -- Field Vec from " << endl;

  // equation->GetFieldValue( PositionMomentum, FieldVec );
  VVectorField *magField = equation->GetField();
  vecgeom::Vector3D<Double_v> Position( PositionMomentum[0], PositionMomentum[1], PositionMomentum[2] );
  // ObtainFieldValue(const Vector3D<double> &position, Vector3D<double> &fieldValue);
  // ObtainFieldValueSIMD( const Vector3D<Double_v> &position, Vector3D<Double_v> &fieldValue);  
  magField-> ObtainFieldValueSIMD( Position, FieldVec );
  
  Double_v dydxVecRegular[gNposmom];
  equation->RightHandSide(PositionMomentum, chargeVec, dydxVecRegular);

  //  Note: Check below assumes that 'FieldVec' is the result of obtaining the field
  bool laneError2 =
      CheckDerivativeInLanesAndReport(chargeVec, PositionMomentum, FieldVec, dydxVecRegular, printContents);

  hasError = hasError || laneError2;

  // 3.) Full case: Use the full equation, and even get back the value of the field
  // ------------------------------------------------------------------------------
  Double_v dydxVecFull[gNposmom];
  vecgeom::Vector3D<Double_v> FieldVecEval;
  equation->EvaluateRhsReturnB(PositionMomentum, chargeVec, dydxVecFull, FieldVecEval);

  bool laneError3 =
      CheckDerivativeInLanesAndReport(chargeVec, PositionMomentum, FieldVecEval, dydxVecFull, printContents);

  hasError = hasError || laneError3;

  return hasError;
}

bool CheckDerivativeInLanesAndReport(const Double_v              & chargeVec,
                                     const Double_v                PositionMomentum[gNposmom],
                                     const ThreeVector_DblV /*vecgeom::Vector3D<Double_v>*/ & FieldVec,
                                     const Double_v                dydxVec[gNposmom],
                                     bool                          printContents)
{
  bool hasError = false;

  for (unsigned int lane = 0; lane < geant::kVecLenD; lane++)
  // int lane= 0;
  {
    double dydxArr[6]     = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // To ensure zeroes at each iteration
    double yArr[6]        = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double MomentumArr[3] = {0.0, 0.0, 0.0};
    double fieldVal[3]    = {0.0, 0.0, 0.0};

    for (int i = 0; i < 6; i++) {
      dydxArr[i] = vecCore::Get(dydxVec[i], lane);
      yArr[i]    = vecCore::Get(PositionMomentum[i], lane);
    }
    ThreeVector_d ForceVec(dydxArr[3], dydxArr[4], dydxArr[5]);

    for (int i = 0; i < 3; i++) {
      fieldVal[i]    = vecCore::Get(FieldVec[i], lane);
      MomentumArr[i] = vecCore::Get(PositionMomentum[3 + i], lane);
    }
    ThreeVector_d BFieldVal(fieldVal[0], fieldVal[1], fieldVal[2]);
    ThreeVector_d Momentum3v(MomentumArr[0], MomentumArr[1], MomentumArr[2]);

    // PositionMomentum[3+i] = Double_v( Momentum[i] );

    if (printContents) {
      cout << " Vectors:      Y   dy/dx   ( lane = " << lane << " ) " << endl;
      for (int i = 0; i < 6; i++) {
        cout << "  [" << i << "] "
             << " " << setw(10) << yArr[i] // PositionMomentum[i]
             << " " << setw(10) << dydxArr[i] << endl;
      }
    }
    double charge = vecCore::Get(chargeVec, lane);

    bool badLane = SanityCheckMagFieldEquation(charge, Momentum3v, BFieldVal, ForceVec, lane);
    if( badLane ) {
       cerr << "Problem seen by Sanity Check in lane " << lane << std::endl;
    }
    hasError     = hasError || badLane;
  }
  return hasError;
}

bool SanityCheckMagFieldEquation(double charge, ThreeVector_d Momentum, ThreeVector_d Field, ThreeVector_d ForceVec,
                                 int lane //  To print (if error)
                                 )
{
  using geant::units::kCLight;  
  constexpr double perMillion = 1e-6;
  bool hasError= false;

  // Check result
  double MdotF = Momentum.Dot(ForceVec);
  double BdotF = Field.Dot(ForceVec);

  double momentumMag = Momentum.Mag();
  double fieldMag    = Field.Mag();
  double sineAngle   = Field.Cross(Momentum).Mag() / (momentumMag * fieldMag);

  double ForceMag = ForceVec.Mag();
  double expectedFmag = kCLight * std::fabs(charge) * fieldMag * sineAngle;
   
  // Tolerance of difference in values (used below)
  double tolerance = perMillion;

  if (gVerbose) {
    cout << "Test output: ( lane = " << lane << " ). " << endl;
  }
  
  if ( std::fabs(ForceMag - expectedFmag ) /*( kCLight * std::fabs(charge) * fieldMag * sineAngle ) */
              > 0.5 * tolerance * (ForceMag+expectedFmag) )
  {
    cerr << "ERROR: Force magnitude is NOT equal to   c * |charge| * |field| * sin( p, B )." << endl;
    cerr << "    Force magnitude = " << ForceMag << endl;
    cerr << "    Expected        = " << expectedFmag << endl; // kCLight * std::fabs(charge) * fieldMag * sineAngle;
    cerr << "    Components: ";
    cerr << " charge = " << charge << " field-Mag= " << fieldMag << " sin(p, B) = " << sineAngle
         << " c_light = " << kCLight
         << endl;
    cerr << "       Force = " << ForceVec[0] << " " << ForceVec[1] << " " << ForceVec[2] << " " << endl;

    // exit(1);
  }

  assert(ForceMag != momentumMag * fieldMag); // Must add coefficient !!

  if (std::fabs(MdotF) > tolerance * Momentum.Mag() * ForceVec.Mag()) {
    cerr << "ERROR: Force due to magnetic field is not perpendicular to momentum!!" << endl;
    cerr << "       Lane = " << lane << endl;
    hasError = true;
  } else if (gVerbose) {
    cout << " Success:  Good (near zero) dot product momentum . force " << endl;
  }
  if (std::fabs(BdotF) > tolerance * Field.Mag() * ForceVec.Mag()) {
    cerr << "ERROR: Force due to magnetic field is not perpendicular to B field!" << std::endl;
    cerr << " Vectors:  BField   Momentum   Force " << std::endl;
    for (int i = 0; i < 3; i++)
      cerr << "  [" << i << "] "
           << " " << setw(10) << Field[i] << " " << setw(10) << Momentum[i] << " " << setw(10) << ForceVec[i]
           << std::endl;

    hasError = true;
  } else if (gVerbose) {
    cout << " Success:  Good (near zero) dot product magnetic-field . force " << std::endl;
  }

  return hasError;
}
