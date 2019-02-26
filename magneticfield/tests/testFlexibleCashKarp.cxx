//
//  Compare the output of a selected steppers (scalar and vector.)
//
//   - adapted for the 'flexible' version of Cash Karp stepper (and similar)
//     which can work as either vector or scalar (using VecCore).
//
//  Configurable output - currently compilation time configuration to choose
//    * variables to output ( position x, y, z ; momentum x, y, z; dX/ds x/y/z dP/ds (x/y/z)
//    * 
//
//  Initially based on testStepperFixed.cc created by Somnath Banerjee in GSoC 2015
//
#include <iomanip>

#include "base/Vector3D.h"
#include <Geant/VectorTypes.h>

using namespace vecCore::math;

#include "Geant/UniformMagField.h" // New type (universal) class
#include "Geant/MagFieldEquation.h"

// #include "VIntegrationStepper.h"   // Supressed ?
#include "Geant/CashKarp.h"

#define BASELINESTEPPER  1

// ---- Future 'general' steppers (i.e. scalar + vector )
// #include "SimpleRunge.h"
// #include "ExactHelixStepper.h"

#ifdef BASELINESTEPPER
#include "Geant/ScalarUniformMagField.h"   // Old type (scalar only) class

#include "Geant/VScalarEquationOfMotion.h"
#include "Geant/ScalarMagFieldEquation.h"

#include "Geant/VScalarIntegrationStepper.h"
// #include "Geant/TClassicalRK4.h"
#include "Geant/GUTCashKarpRKF45.h"
#endif
  
#include "Geant/SystemOfUnits.h"

// using geant::meter;
using geant::units::millimeter;
using geant::units::second;
using geant::units::eplus;
using geant::units::tesla;
using geant::units::degree;

// using namespace std;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;

template <typename Real_v, typename Stepper_t> // , typename Equation_t>
bool TestFlexibleStepper(Stepper_t *stepper);  // , Equation_t *equation);

/* -----------------------------SETTINGS-------------------------------- */

// Parameters of test - values can be modified here (or in arguments to main)
int no_of_steps    = 100;  // No. of Steps for the stepper
int stepper_no     = 5;    // Choose stepper no., for refernce see above
double step_len_mm = 200.; // meant as millimeter;  //Step length
double z_field_in  = DBL_MAX;

std::string datafile("cmsmagfield2015.txt");

bool debug = false;       // true;

int    magExp10 = 6;      // Differences are magnified by 10^this (magExp10)

double facStepLen= 1.05;  // Factor for progressive lengthening of step size
double magFactor = 1.0 * pow(10,magExp10);  // Magnification factor for difference(s)

/* ------------------ Parameters for PRINTING --------------------------- */
// Choice of output coordinates
int columns[] = {
#ifdef CHECK_POSITION
    1, 1, 1, // position  x, y, z
    0, 0, 0, // momentum  x, y, z
    1, 1, 1, // dydx pos  x, y, z
    0, 0, 0  // dydx mom  x, y, z
#else   
    0, 0, 0, // position  x, y, z
    1, 1, 1, // momentum  x, y, z
    0, 0, 0, // dydx pos  x, y, z
    1, 1, 1  // dydx mom  x, y, z
#endif    
};           // Variables in yOut[] & dydx[] we want to display - 0 for No, 1 for yes

bool printDiff   = 1, // Print the diffrence
     printRef    = 0, // Print the reference Solution
     printDrvRef = 1  // Print the derivative of Reference ( ORed with above to give Derivative )
       ;
bool printSep = 1; // separator  '|'
bool printInp = 0; // print the input values
bool printInpX= 0;  // print the input values for Ref

// #if PRINT_ERR_VALS
bool printErr    = true,    // Print the predicted Error ? 
     printErrRef = false,
     printErrDiff= true;

const unsigned int nwdf = 12; // Width for float/double
// -------   End of Parameters for PRINTING ------------------------------ */

void processArguments(int argc, char *args[])
{
  // Checking for command line values :
  if (argc > 1) stepper_no  = atoi(args[1]);
  if (argc > 2) step_len_mm = (float)(atof(args[2])); // *mm);   //  stof( std::string, size_t ... )
  if (argc > 3) no_of_steps = atoi(args[3]);
  if (argc > 4) z_field_in  = (float)(atof(args[4])); // tesla
}

void printBanner();

constexpr unsigned int Nposmom = 6; // Position 3-vec + Momentum 3-vec

// #define UNIFORM_FIELD   1

#ifdef UNIFORM_FIELD
#include "Geant/UniformMagField.h" // New type (universal) class
#include "Geant/ScalarUniformMagField.h" // Old type - for single track, i.e. one position, one B-value

using FieldType      = UniformMagField;
using Field_Type_Scalar= ScalarUniformMagField;
#else
#include "CMSmagField.h"          
#include "RZMagFieldFromMap.h"   // Scalar field in sub-directory (for testing)
// #include "ScalarCMSMagField.h"

using FieldType      =   CMSmagField;
using Field_Type_Scalar= RZMagFieldFromMap; // ScalarCMSMagField;
#endif

FieldType      * gvField = nullptr;
using GvEquationType = MagFieldEquation<FieldType>;
using StepperType = CashKarp<GvEquationType, Nposmom>;

vecgeom::Vector3D<float> fieldValInput(0.0, 0.0, 0.0);   // labelled ThreeVector_f below (within main)


// ---------------  Main ------------------------------------------------ */

int main(int argc, char *args[])
{

  // using Backend = vecgeom::kVc ;
  // typedef typename Backend::precision_v Double_v;

  using Double_v = geant::Double_v;
  // template <typename Tpod>
  // using Vector3D  = vecgeom::Vector3D<Tpod>;
  using ThreeVector_f = vecgeom::Vector3D<float>;
  using ThreeVector_d = vecgeom::Vector3D<double>;
  // using ThreeVectorSimd = vecgeom::Vector3D<Double_v>;

  if (debug) cout << "Running with debug (progress) messages." << endl;

  processArguments(argc, args);

  double x_field = 0., y_field = 0., z_field = 0.;      // Uniform Magnetic Field (x,y,z)
  z_field = (z_field_in < DBL_MAX) ? z_field_in : -4.0; //  Tesla // *tesla ;

  if (debug) cout << "---- Creating UniformMagField object" << endl;
  // Vector3D<float>
  // ThreeVector_f
  fieldValInput= ThreeVector_f(x_field, y_field, z_field); // Default = ( 0.0, 0.0, -1.0);
  fieldValInput *= geant::units::tesla;

  // Field
#ifdef UNIFORM_FIELD  
  gvField = new UniformMagField(fieldValInput); // New class
  //  UniformMagField( geant::units::tesla * ThreeVector_f(x_field, y_field, z_field) );
#else
  cerr << "Using CMS field from file " << "cmsmagfield2015.txt" << endl;
  gvField = new FieldType(datafile.c_str());
#endif  

  std::cout << " Simulation parameters:  stepper_no= " << stepper_no 
            << " step length (mm) = " << step_len_mm 
            << " number of steps  = " << no_of_steps 
            << " B_z value (Tesla) = " << z_field  << std::endl;
  
  if (debug) {
    cout << "#  Initial  Field strength (GeantV) = " << x_field << " , " << y_field << " , " << z_field << " T (tesla) "
         << endl;

    ThreeVector_d origin(0.0, 0.0, 0.0), fieldValue;
    gvField->GetFieldValue(origin, fieldValue);
    cout << "#    Values in object created       = " << (1.0 / geant::units::tesla) * fieldValue.x() << ",  "
         << (1.0 / geant::units::tesla) * fieldValue.y() << ",  " << (1.0 / geant::units::tesla) * fieldValue.z()
         << endl;
    cout << endl;
  }

  // Create Equation :
  if (debug) cout << "Create Equation" << endl;

  // 1. Original way of creating an equation
  auto magEquation = new GvEquationType(gvField);
  // new MagFieldEquation<UniformMagField>(gvField);

  if (debug) cout << "----Equation instantiated. " << endl;

  // Create a stepper :
  if (debug) cout << "---- Preparing to create (Vector) CashKarpRKF45 Stepper " << endl;

  // 2. Create stepper

  // auto myStepper= new StepperType(magEquation);
  StepperType myStepper2(magEquation);
  auto myStepper = &myStepper2;
  // myStepper = new VectorCashKarp<GvEquationType,Nposmom>(gvEquation);
  
  if (debug) {
    cout << "---- constructed (flexible) CashKarp" << endl;
    cout << " Stepper  information: ptr= " << myStepper << endl;
    // cout << "       full " << *myStepper << endl;
  }

  // Phase 1 - get it to work without cloning

// #ifdef CREATE_SCALAR_STEPPER
//   cout << endl;
//   // cout << " =============================================================================" << endl;  
//   cout << " 1. Testing scalar field, equation and stepper     " << endl;
//   cout << " =============================================================================" << endl;
  
//   using ScalarFieldType=    
//   using ScalarEquationType= ScalarMagFieldEquation<ScalarFieldType, Nposmom>;
//   using ScalarStepperType=   CashKarp<ScalarEquationType, Nposmom>;
//   bool okScalar = TestFlexibleStepper<double, ScalarStepperType>(myStepper);
// #endif
  
  cout << endl;
  cout << " =============================================================================" << endl;
  cout << " 2. Testing scalar version of 'flexible' field, equation and stepper." << endl;
  cout << " =============================================================================" << endl;
  
  bool okScalar = TestFlexibleStepper<double, StepperType>(myStepper);
  
  // cout << " Testing Vector Float_v.   " << endl;
  // cout << " 3. a) Testing Vector Float_v version of 'flexible' field, equation and stepper." << endl;  
  // bool okVecFloat  = TestFlexibleStepper<Float_v, StepperType >(myStepper); // , magEquation);
  
  // cout << " Testing Vec Double . " << endl;
  cout << endl;
  cout << " =============================================================================" << endl;
  cout << " 3. Testing Vector Double_v version of 'flexible' field, equation and stepper." << endl;
  cout << " =============================================================================" << endl;  
  bool okVecDouble = TestFlexibleStepper<Double_v, StepperType>(myStepper);

  bool good = okScalar && // okVecFloat &&
              okVecDouble;

  // delete myStepper;  // Only if object was new'ed !!

  if (debug) cout << "----deletion of stepper done " << endl;

  // delete gvEquation;  // The stepper now takes ownership of the equation
  // delete gvEquation2;
  delete gvField;

  if (debug) cout << "\n\n#-------------End of all tests -----------------\n";

  return good;
}

const double mmGVf = geant::units::millimeter;
const double ppGVf = geant::units::GeV; //   it is really  momentum * c_light
                                        //   Else it must be divided by geant::c_light;

template <typename Real_v, typename scalar_t, int Ncomp>
void SelectOutput(const Real_v arrVec[Ncomp], scalar_t arrLane[Ncomp], int lane)
{
  for (unsigned int i = 0; i < Ncomp; i++) {
    arrLane[i] = vecCore::Get(arrVec[i], lane);
  }
}

// template < typename scalar_t, typename scalar_t, int Ncomp >
// void SelectOutput( const double arrVec[Ncomp] , scalar_t arrLane[Ncomp], int lane )
// template < typename double, typename double, int Ncomp >
template <int Ncomp>
void SelectOutput(const double arrVec[Ncomp], double arrLane[Ncomp], int /*lane*/)
{
  for (unsigned int i = 0; i < Ncomp; i++) {
    arrLane[i] = arrVec[i];
  }
}

template <typename Real_v, typename Stepper_t> // , typename Equation_t>
bool TestFlexibleStepper(Stepper_t *stepper)   // , Equation_t *equation)
{
  // Initialising 'parameters'
  Real_v chargeVec(-1.);
  double stepLengthValue = step_len_mm * geant::units::millimeter;

  // auto scratch = stepper->ObtainScratchSpace<Double_v>();
  // Scratch space needed for this thread / task ...

  // Initial coordinates , momentum
  double x_pos = 0., y_pos = 0., z_pos = 0.; // pos - position  : input unit = mm
  double x_mom = 0., y_mom = 1., z_mom = 1.; // mom - momentum  : input unit = GeV / c

  if (debug) {
    cout << "TestFlexibleStepper called again ----------  " << endl;
    cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;
  }
  // Double_v yIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
  //                  x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};
  Real_v yInVec[] = {Real_v(x_pos * mmGVf), Real_v(y_pos * mmGVf), Real_v(z_pos * mmGVf),
                     Real_v(x_mom * ppGVf), Real_v(y_mom * ppGVf), Real_v(z_mom * ppGVf)};
  double yInVecX[] = {x_pos * mmGVf, y_pos * mmGVf, z_pos * mmGVf,
                     x_mom * ppGVf, y_mom * ppGVf, z_mom * ppGVf };  
  if (debug) cout << "Initialized yIn: values [0]= " << yInVec[0] << endl;

  // double yInX[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
  //                 x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};

  // Should be able to share the Equation -- eventually
  // For now, it checks that it was Done() -- and fails an assert

  // Empty buckets for results
  Real_v dydxVec[8] = {Real_v(0.), 0., 0., 0., 0., 0., 0., 0.}, // 2 extra safety buffer
      yOutVec[8] = {0., 0., 0., 0., 0., 0., 0., 0.}, yErrVec[8] = {0., 0., 0., 0., 0., 0., 0., 0.};

#ifdef BASELINESTEPPER
  const double mmRef = mmGVf; // Unit for reference of lenght   - milli-meter
  const double ppRef = ppGVf; // Unit for reference of momentum - GeV / c^2

  // auto gvEquation2 = new GvEquationType(gvField);
  // new TMagFieldEquation<ScalarUniformMagField, Nposmom>(gvField);  

#ifdef UNIFORM_FIELD
  // using Field_Type_Scalar= ScalarUniformMagField;
  auto gvScalarField = new ScalarUniformMagField(fieldValInput);
#else
  auto gvScalarField = new RZMagFieldFromMap(datafile);
#endif

  using GvScalarEquationType = ScalarMagFieldEquation<Field_Type_Scalar, Nposmom>;
  auto  gvScalarEquation2 = new GvScalarEquationType(gvScalarField);
  
  // Creating the baseline stepper
  //  auto exactStepperGV = new TClassicalRK4<GvEquationType, Nposmom>(gvEquation2);
  // cout << "#  Reference stepper is: TClassicalRK4<GvEquationType,Nposmom>(gvScalarEquation2);" << endl;  
  auto exactStepperGV = new GUTCashKarpRKF45<GvScalarEquationType, Nposmom>(gvScalarEquation2);  
  cout << "#  Reference stepper is Cash Karp, i.e. GUTCashKarpRKF45 <GvEquationType,Nposmom> (gvScalarEquation2)" << endl;

  // new TSimpleRunge<GvEquationType,Nposmom>(gvEquation2);
  // new GUExactHelixStepper(gvEquation2);

  auto exactStepper = exactStepperGV;

  double dydxRef[8]  = {0., 0., 0., 0., 0., 0., 0., 0.},
         yOutVecX[8] = {0., 0., 0., 0., 0., 0., 0., 0.},
         yErrVecX[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
  cout << " mmRef= " << mmRef << "   ppRef= " << ppRef << endl;

  // Real_v yInX[] = {x_pos * mmRef, y_pos * mmRef, z_pos * mmRef, x_mom * ppRef, y_mom * ppRef, z_mom * ppRef};

  // Simple arrays for outputing and/or checking values
  double yinX[Nposmom], youtX[Nposmom], yerrX[Nposmom];
#endif
  cout << "# step_len_mm = " << step_len_mm;

  // double stepLengthRef = step_len_mm * mmRef;
  /*-----------------------END PREPARING STEPPER---------------------------*/

  int lane = 0; // Lane for printing
  cout << "# Results of lane = " << lane << endl;
  cout << "#-------------------------------------------------------------------#" << endl;  
  printBanner();

  // Units
  const char *nameUnitLength   = "mm";
  const char *nameUnitMomentum = "GeV/c";
  cout << setw(6) << "#Numbr" << " ";
  cout << setw(6) << "(mm) ";
  for (int i = 0; i < 6; i++)
  {
    if (columns[i])
    {
      if (printSep) cout << " | "; // Separator
      const char *nameUnit = (i < 3) ? nameUnitLength : nameUnitMomentum;
      cout << setw(nwdf) << nameUnit;
      if (printRef) cout << setw(nwdf) << nameUnit; // Reference output
      if (printErr) cout << setw(nwdf) << nameUnit; // Estim. error
      if (printDiff) cout << setw(15) << nameUnit;  // Diff  new - ref
    }
  }    
  cout << "\n";

  //-> Then print the data

  // cout.setf (ios_base::scientific);
  // cout.setf (ios_base::scientific);
  cout.precision(3);

  /*----------------NOW STEPPING-----------------*/
  // no_of_steps = 25;
    
  for (int j = 0; j < no_of_steps; j++) {
    cout << setw(6) << j << " " ; // Printing Step number
    cout << setw(8) << stepLengthValue << " ";
    
    Real_v step_len(stepLengthValue);

    stepper->RightHandSideInl(   yInVec, chargeVec, dydxVec); // compute dydx - to supply the stepper
#ifdef BASELINESTEPPER
    double charge= vecCore::Get( chargeVec, 0 );
    exactStepper->RightHandSideInl( yInVecX, charge,   dydxRef);
    // exactStepper->RightHandSideVIS( yInVecX, dydxRef); // compute the value of dydx for the exact stepper
#endif

    using scalar_t        = double;
    scalar_t yin[Nposmom] = {-999., -99.9, -9.99, 10.0, 1.0, 0.1};
    scalar_t yout[Nposmom], dydx[Nposmom], yerr[Nposmom];
    if (j > 0) // Do nothing for j=0, so we can print the initial points!
    {
      stepper->StepWithErrorEstimate(yInVec, dydxVec, chargeVec, step_len, yOutVec, yErrVec /*, scratch */);
//             *********************
#ifdef BASELINESTEPPER
      exactStepperGV->StepWithErrorEstimate(yInVecX, dydxRef, charge, stepLengthValue, yOutVecX,
                                            yErrVecX /* ,scratchRef*/);
//                    *********************
#endif
    }

    // Try1 :
    // SelectoOutput( yInVec,  yin,  lane);
    // Try2 :
    // void CopyToArr(Real_v realArr[], scalar_t scalArr[], int lane ) {
    // SelectOutput<Real_v, scalar_t, N>(realArr, scalArr, lane); }
    // CopyToArr( yInVec,  yin,  lane);

    // Select one lane for printing.

    // General code - for vecCore types or built-in types

    // yin[i] = vecCore::Get(yInVec[i], lane);    
    SelectOutput<Real_v, scalar_t, Nposmom>(yInVec, yin, lane);

    if (j > 0)  //  yOut is not yet set for j=0    
      SelectOutput<Real_v, scalar_t, Nposmom>(yOutVec, yout, lane); // yout[i] = vecCore::Get(yOutVec[i], lane) : 
    else
      SelectOutput<Real_v, scalar_t, Nposmom>(yInVec, yout, lane);  // yout[i] = vecCore::Get(yInVec[i], lane); 
    SelectOutput<Real_v, scalar_t, Nposmom>(yErrVec, yerr, lane);   // yerr[i] = vecCore::Get(yErrVec[i], lane);
    SelectOutput<Real_v, scalar_t, Nposmom>(dydxVec, dydx, lane);   // dydx[i] = vecCore::Get(dydxVec[i], lane);

#ifdef BASELINESTEPPER
    // using Real_vRef= double;
    // int   laneRef= 0;  //  This is scalar - there is no other lane !!
    // SelectOutput<Real_vRef, scalar_t, Nposmom>(yInVecX, yinX, laneRef);
    /* if (j > 0)
      SelectOutput<Real_vRef, scalar_t, Nposmom>(yOutVecX, youtX, laneRef);
    else
      SelectOutput<Real_vRef, scalar_t, Nposmom>(yInVecX, youtX, laneRef);  **/
    // SelectOutput<Real_vRef, scalar_t, Nposmom>(yErrVecX, yerrX, laneRef);
    // SelectOutput<Real_vRef, scalar_t, Nposmom>(dydxRef, dydx, laneRef);

    for (int i = 0; i < 6; i++)
    {
       yinX[i]= yInVecX[i];
       if (j > 0)       
          youtX[i] = yOutVecX[i]; // ( j > 0 ) ? yOutVecX[i] : yInVecX[i] ;
       else
          youtX[i] = yInVecX[i]; // ( j > 0 ) ? yOutVecX[i] : yInVecX[i] ;       
       yerrX[i] = yErrVecX[i];
    }
#endif
    // -->> End of Selecting 'lane' for printing

    //-> Then print the data
    cout.setf(std::ios_base::fixed);
    cout.precision(2);
    for (int i = 0; i < 3; i++)
    {
      // Check the difference first
      double magDiff = yout[i] / mmGVf - youtX[i] / mmRef;
      if( magFactor * fabs( magDiff ) * mmRef > 0.01 * fabs( youtX[i] ) ) {
        cerr << " Non-zero difference in step # " << j << " len= " << stepLengthValue
             << "column[ " << i << " ] "
             << " Diff = " << magDiff << " vs Reference " << youtX[i] / mmRef << endl;
      }
       
      if (columns[i])
      {
        cout << std::defaultfloat;         
        if (printSep)  cout << " | "; // Separator
        if (printInp)  cout << setw(nwdf - 2) << yin[i] / mmGVf << " ";
#ifdef BASELINESTEPPER
        if (printInpX)   cout << setw(nwdf - 2) << yinX[i] / mmRef << " ";
#endif        
        cout << setw(nwdf) << yout[i] / mmGVf;
#ifdef BASELINESTEPPER        
        if (printRef)    cout << setw(nwdf) << youtX[i] / mmRef << " "; // Reference Solution
        if (printDiff)   cout << /* "  Df= " << */ setw(nwdf) << magFactor * (yout[i] / mmGVf - youtX[i] / mmRef) << " ";
        // if (printDrvRef) cout << setw(nwdf - 2) << dydxRef[i] / mmRef;
#endif
        if (printErr)    cout << std::scientific << setw(nwdf) << yerr[i] / ppGVf;
           
// #if PRINT_ERR_VALS
        // cout.unsetf(std::ios_base::fixed);
        // cout.setf(std::ios_base::scientific);
        if (printErrRef)  cout << setw(nwdf) << yerrX[i] / mmGVf << " ";
        cout << std::defaultfloat;
        double errDiff= ( yerr[i] / mmGVf - yerrX[i] / mmRef );
        if( errDiff != 0.0 && ( errDiff < 1.e-4 || errDiff > 1.e+5 ) )
           cout << std::scientific;        
        if (printErrDiff) cout // << "mag*ErrDif= " 
                               << setw(nwdf-4) << magFactor * errDiff << " ";
// #endif
        cout << std::defaultfloat;
      }
    }      
    cout.precision(3);

    cout.unsetf(std::ios_base::fixed);
    cout.setf(std::ios_base::scientific);
    for (int i = 3; i < 6; i++)
      if (columns[i]) {
        if (printSep) cout << " | "; // Separator
        if (printInp) cout << setw(nwdf - 1) << yin[i] / ppGVf << " ";

        double outVal=  yout[i] / ppGVf;
        if ( outVal == 0.0 || ( fabs(outVal) > 1.0e-4 && fabs(outVal) < 1.0e+4) )
        {
           cout.unsetf(std::ios_base::scientific);           
           cout.setf(std::ios_base::fixed);
           cout << setw(nwdf) << outVal; // yout[i] / ppGVf;
           cout.unsetf(std::ios_base::fixed);           
           cout.setf(std::ios_base::scientific);                   
        } else {
           // cout.unsetf(std::ios_base::fixed);           
           // cout.setf(std::ios_base::scientific);
           cout << setw(nwdf) << outVal; // yout[i] / ppGVf;
        }
        // cout << setw(nwdf) << outVal; // yout[i] / ppGVf;

#ifdef BASELINESTEPPER
        if (printInpX) cout << setw(nwdf - 1) << yinX[i] / ppRef << " ";
        if (printRef)  cout << setw(nwdf) << youtX[i] / ppRef; // Reference Solution
        double rawDiffGeVc = ((yout[i] / ppGVf) - (youtX[i] / ppRef));        
        if (printDiff)
        {
           // cout << "  Df= " << setw(nwdf + 2) << magFactor * rawDiffGeVc;
           double scaledDf = magFactor * std::fabs(rawDiffGeVc);
           if ( scaledDf == 0.0 ) { // || ( scaledDf > 1.0e-4 && scaledDf < 1.0e+4) ) {
              // cout.setf(std::ios_base::fixed);
              // cout << setw(nwdf + 2) << magFactor * rawDiffGeVc;
              cout << setw(nwdf + 2) << "0.00";
              // cout.setf(std::ios_base::scientific);
           } else {
              cout << setw(nwdf + 2) << magFactor * rawDiffGeVc;
           }
        }
#endif
        if (printErr) {
           cout << std::scientific;           
           cout << setw(nwdf) << yerr[i] / ppGVf;
        }
// #if PRINT_ERR_VALS
        // cout.unsetf(std::ios_base::fixed);
        // cout.setf(std::ios_base::scientific);
        if (printErrRef)  cout << setw(nwdf) << yerrX[i] / ppGVf << " ";

        double errDiff= ( yerr[i] / ppGVf - yerrX[i] / ppRef );

        cout << std::defaultfloat;            
        if( errDiff != 0.0 && ( errDiff < 1.e-4 || errDiff > 1.e+5 ) )
           cout << std::scientific;
        if (printErrDiff) cout // << "mag*ErrDif= " 
                               << setw(nwdf-4) << magFactor * errDiff << " ";
// #endif
        cout << std::defaultfloat;    
        
      }
    cout.unsetf(std::ios_base::scientific);

    for (int i = 0; i < 6; i++) // Print auxiliary components
    {
      double unitGVf = 1;
      // if( i < 3 )             // length / length
      if (i >= 3) {
        unitGVf = ppGVf / mmGVf; //  dp / ds
      }
#ifdef BASELINESTEPPER
      double unitRef=1;
      if (i >= 3) {        
        unitRef = ppRef / mmRef; //  dp / ds
#endif
      }

      if (i == 0) { // length / length
        // cout.unsetf (ios_base::scientific);
        cout.setf(std::ios_base::fixed);
      } else if (i == 3) {
        cout.unsetf(std::ios_base::fixed);
        cout.setf(std::ios_base::scientific);
      }

      if (columns[i + 6]) {
        // cout << " dy/dx["<<i<<"] = ";
        if (printSep) cout << " | "; // Separator
        // cout << std::defaultfloat;
        cout << std::scientific;
        cout << setw(nwdf) << dydx[i] / unitGVf;
        cout << std::defaultfloat;        
#ifdef BASELINESTEPPER
        if (printRef | printDrvRef ) cout << setw(nwdf) << dydxRef[i] / unitRef << " "; // Reference Solution
        if (printDiff)                                            // Deriv )
           cout /* << "  Df= " */
                << setw(nwdf-3) <<  magFactor * ( (dydx[i] / unitGVf) - (dydxRef[i] / unitRef) );
        // cout << setw(nwdf) << (dydx[i] / unitGVf) - (dydxRef[i] / unitRef);        
#endif
        // bool printDiffDeriv = true;
      }
      // if( i == 2 )     { cout.unsetf(std::ios_base::fixed);      }
      // else if ( i==5 ) { cout.unsetf(std::ios_base::scientific); }
    }
    cout.unsetf(std::ios_base::scientific);
    if (j > 0) // Step 0 did not move -- printed the starting values
    {
      // cout.unsetf(std::ios_base::scientific);
      cout.setf(std::ios_base::fixed);
      cout.precision(2);
      cout << setw(nwdf) << atan2(yout[1], yout[0]) / degree;

#ifdef BASELINESTEPPER // atan(yout[1]/yout[0])/degree;
      if (printRef) cout << setw(nwdf) << atan2(youtX[1], youtX[0]) / degree;
#endif

      // Prepare the state of yIn(Vec)
      // Copy yout into yIn
      for (unsigned int i = 0; i < Nposmom; i++) {
        yInVec[i] = yOutVec[i];
#ifdef BASELINESTEPPER
        yInVecX[i] = yOutVecX[i];
#endif
      }
    }

    stepLengthValue *= facStepLen;
    
    cout << "\n";
  }

  /*-----------------END-STEPPING------------------*/

  if (debug) cout << "----Stepping done " << endl;

/*------ Clean up ------*/
#ifdef BASELINESTEPPER
  // exactStepper->InformDone(); // Old protocol ... 
  delete exactStepper;
#endif

  return true; // Should perform some checks -- for not just PASS
}

void printBanner()
{
  /*---------------------------------------------------*/
  //        -> First Print the (commented) title header
  cout << "\n#";
  cout << setw(5) << "Step";
  cout << setw(9) << "Size";
  for (int i = 0; i < 6; i++)
    if (columns[i]) {
      if (printSep)  cout << " | "; // Separator
      if (printInp)  cout << setw(nwdf - 2) << "yIn[" << i << "]";
#ifdef BASELINESTEPPER            
      if (printInpX) cout << setw(nwdf - 2) << "yInX[" << i << "]";
#endif      
      cout << setw(nwdf - 2) << "yOut[" << i << "]";
#ifdef BASELINESTEPPER      
      if (printRef) cout << setw(nwdf - 2) << "yOutX[" << i << "]";
      if (printDiff) cout << " 10^" << setw(1) << magExp10
                          << setw(nwdf-9) <<  "*(yO-Ref)"; // [" << i << "]";
#endif
      // if (printDeriv)  cout << setw(nwdf - 2) << "dydx[" << i << "]";
      if (printErr) cout << setw(nwdf - 1) << "yErr[" << i << "]";
      if (printErrRef)  cout << setw(nwdf) << "ErrRef" << " ";
      if (printErrDiff) cout << setw(nwdf) << "m*(yE-Ref)" << " ";
                               
    }
  for (int i = 0; i < 6; i++)
    if (columns[i + 6]) {
      if (printSep) cout << " | "; // Separator
      if (printInp) cout << setw(nwdf - 2) << "yIn[" << i << "]";
      cout << setw(nwdf - 2) << "dydx[" << i << "]";
#ifdef BASELINESTEPPER      
      if (printRef | printDrvRef ) 
         cout << setw(nwdf - 2) << "dydxRef[" << i << "]";               
#endif      
      // if( printErr ) cout << setw(nwdf-2) << "yErr[" << i << "]";
      if (printDiff) //  printDiffDeriv )
        cout << setw(nwdf - 2) << "dydx-Ref[" << i << "]";
    }
  cout << setw(nwdf) << "tan-1(y/x)";
  cout << "\n";
}
