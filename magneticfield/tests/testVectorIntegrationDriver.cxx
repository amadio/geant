//
//  Compare the output of different steppers
//
//  Based on testStepperFixed.cc
//    was the work of Somnath Banerjee in GSoC 2015
//
#include <iomanip>
#include <ctime>

#include "Geant/SystemOfUnits.h"

// using fieldUnits::meter;
using geant::units::millimeter;
using geant::units::second;
using geant::units::eplus;
using geant::units::tesla;
using geant::units::degree;

#include <base/Vector3D.h>
#include <Geant/VectorTypes.h>

#include "Geant/UniformMagField.h"
#include "Geant/MagFieldEquation.h"

// #include "IntegrationStepper.h"

#include "Geant/CashKarp.h"
#include "Geant/SimpleIntegrationDriver.h"
#include "Geant/RollingIntegrationDriver.h"

#include "Geant/FieldTrack.h"

// #define  NEW_SCALAR_FIELD 1

#define USECMSFIELD 1

#ifdef USECMSFIELD
#include "CMSmagField.h"
// #include "RZMagFieldFromMap.h"
// using ScalarCMSmagField = RZMagFieldFromMap;
using ScalarCMSmagField = CMSmagField;   // Trial 2019.03.20 16:30
#include "Geant/Utils.h"
// #else
#endif

#ifndef NEW_SCALAR_FIELD
//  Transition measure --- compare to old Scalar field types 2017.11.16
#include "Geant/ScalarUniformMagField.h"
#include "Geant/ScalarMagFieldEquation.h"
#include "Geant/StepperFactory.h"
#include "Geant/ScalarFieldTrack.h"
#include "Geant/ScalarIntegrationDriver.h"
#endif

#include <stdlib.h>

// using namespace std;
using std::cout;
using std::cerr;
using std::endl;

using Double_v = geant::Double_v;
using Bool_v   = vecCore::Mask_v<Double_v>;

constexpr unsigned int Nposmom = 6; // Position 3-vec + Momentum 3-vec


//  A few parameters ...
//
int   gVerbose = 0;  // Global verbosity level ... 



int CompareSixVectors( const FieldTrack  & outVec,
                       const ScalarFieldTrack & yTrackOut,
                       const std::string & varName,
                       int                 varNum,
                       double              hStep,
                       FieldTrack        & yVecInput,
                       ScalarFieldTrack  & yScalInput,
                       double            & maxRelDiff                      
   )
   
{
   const int  hiPrec= 8;  // High precision output
   const int lowPrec= 5;  
   const int    Npm = Nposmom;

   constexpr double threshold = 1.0e-6;    
   constexpr double kTiny=  1.0e-30;    

   int oldPrecErr= -1;
   double outSer[8] = { 0., 0.1, 0.2, 0.3, 0.4, 0.5 };

   yTrackOut.DumpToArray(outSer);
   
   // Return number of differences found between lanes
   // 
   int     diffFound = 0;
   maxRelDiff = 0.0;

   for (unsigned int j = 0; j < Npm; ++j)
   {
      double  dif   = outVec[j] - outSer[j];
      double midAbs = 0.5 * ( std::fabs(outVec[j]) + std::fabs(outSer[j]) ) + kTiny ;
      if( std::fabs(dif) > threshold * midAbs )
      {
         if( ! diffFound ) {
            oldPrecErr= cerr.precision(hiPrec);
            cerr << "##-------------------------------------------------------------------------------" << endl;
            cerr << " Compare-Six> Difference found in " << varName 
                 << " : vector [ i= " << varNum << " ] and for hStep = " << hStep << ": " << endl;
            diffFound++;
         }
         cerr.precision(hiPrec);
         cerr << "   [" << j << " ]:"
              << "  vec = " << std::setw(6+hiPrec) << outVec[j]
              << "  ser = " << std::setw(6+hiPrec) << outSer[j]
              << " diff = " << std::setw(6+hiPrec) << dif
              << std::setprecision(lowPrec)
              << " rel dif = " << std::setw(6+hiPrec) << dif / midAbs
              << endl;
      }
   }
   
   if( diffFound ) {
      cerr.precision(oldPrecErr);
      cerr << "##-------------------------------------------------------------------------------" << endl;
   }
   
   if( gVerbose || diffFound ) {
      int oldPrecOut = cout.precision(hiPrec);
      cout << " Vector  len = " <<  outVec.GetCurveLength()    << " yOut [" << varNum << "]  is : " << outVec
           << " for yInput: " << yVecInput << endl;
      cout << " Serial  len = " <<  yTrackOut.GetCurveLength() << " yTrackOut  is : " << yTrackOut << endl
           << " for yTrackIn: " << yScalInput
           << " for hstep= " << hStep << endl << endl;
      cout.precision(oldPrecOut);
      cout << "##-------------------------------------------------------------------------------" << endl;
   }
   
   return diffFound;
}

int main(int argc, char *argv[])
{
  // template <typename T> using Vector3D = vecgeom::Vector3D<T>;
  using ThreeVector_d = vecgeom::Vector3D<double>;

  bool verbose= true;
  
  
#ifdef USECMSFIELD
  using Field_Type        = CMSmagField;
  using Field_Type_Scalar = ScalarCMSmagField;
  // using Field_Type_Scalar = RZMagFieldFromMap; // Use same class for 'scalar' - i.e. one track at a time
  // using Field_Type_Scalar = TemplateCMSmagField<vecgeom::kScalar>;
#else
  using Field_Type = UniformMagField; // TemplateScalarUniformMagField<Backend>;
#ifdef NEW_SCALAR_FIELD
  // New types ... under development  2017.11.16
  using Field_Type_Scalar = UniformMagField;
#else
  using Field_Type_Scalar = ScalarUniformMagField;
#endif
#endif

  using GvEquationType = MagFieldEquation<Field_Type>; // , Nposmom>;

  /* -----------------------------SETTINGS-------------------------------- */

  /* Parameters of test - Modify values  */

  // int no_of_steps = 20;         // No. of Steps for the stepper
  int stepper_no     = 5;    // Choose stepper no., for refernce see above
  double step_len_mm = 200.; // meant as millimeter;  //Step length
  double z_field_in  = DBL_MAX;
  int    trackToPrint = 11;   // If # < vecLength, print that track in vector mode
  
  double x_field = 0., // Uniform Magnetic Field (x,y,z)
      y_field    = 0., //  Tesla // *tesla ;
      z_field;

  if (z_field_in < DBL_MAX)
    z_field = z_field_in;
  else
    z_field = -1.0; //  * Tesla

  if (argc > 2) {
     int candTrack2pr = atoi(argv[2]);
     if( candTrack2pr >= 0 ) {
       std::cout << "- Chosen track= " << candTrack2pr << " for printing. " << std::endl;
       trackToPrint = candTrack2pr;
     } else {
        std::cout << "- Details of no track will be printed. ";
        if( trackToPrint != -1 ) 
           std::cout << "Candidate # = " << candTrack2pr << " (expected = -1) .";
        std::cout << std::endl;
       trackToPrint = -1; // None printed
     }
  }
  
// Field
#ifdef USECMSFIELD
  std::string datafile(geant::GetDataFileLocation(argc, argv, "cmsmagfield2015.txt"));

  cerr << "Using CMS field from file " << "cmsmagfield2015.txt" << endl;
  auto gvField = new Field_Type(datafile.c_str());

  (void) x_field;
  (void) y_field;
  (void) z_field;
#else
  (void)argc;
  (void)argv;
  cerr << "Using Uniform Magnetic field with strength  ( " << x_field << " , " << y_field << " , "
       << z_field << " )  Tesla. " << endl;
  auto gvField = new UniformMagField(geant::units::tesla * ThreeVector_d(x_field, y_field, z_field));

  cout << "#  Initial  Field strength (GeantV) = " << x_field << " , " << y_field << " , " << z_field
       // << (1.0/geant::units::tesla) * gvField->GetValue()->X() << ",  "
       << " Tesla " << endl;
#endif  
  // cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;

  // Create an Equation :
  auto gvEquation = new MagFieldEquation<Field_Type>(gvField);
  // TemplateFieldEquationFactory<Double_v>::CreateMagEquation<Field_Type >(gvField);

  if( verbose ) { cout << " Equation object created - address = " << gvEquation << endl; } 
  cout << "# step lenght (mm) = " << step_len_mm << endl;

  //=======================Test part for Integration driver====================
  double hminimum = 0.2;
  double epsTol = 1.0e-5;

  //==========  Vector Driver: start preparation ===============
  using StepperType = CashKarp<GvEquationType, Nposmom>;
  auto myStepper    = new StepperType(gvEquation);

  // FlexIntegrationDriver *vectorDriver;
// #define SIMPLE_DRIVER   1 
#ifdef SIMPLE_DRIVER  
  using DriverType = SimpleIntegrationDriver<StepperType, Nposmom>;
#else  
  using DriverType = RollingIntegrationDriver<StepperType, Nposmom>;
#endif  
  
  auto vectorDriver =
      new DriverType(hminimum, myStepper, epsTol,  Nposmom );
  if( verbose ) { cout << " Vector Driver created." << endl; } 
  // ========== Vector Driver prepared ========================

  //========= Preparing scalar Integration Driver ============

#ifdef USECMSFIELD
  // auto gvFieldScalar         = new Field_Type_Scalar(datafile.c_str());
  auto gvFieldScalar         = new ScalarCMSmagField(datafile.c_str());  
  // using GvEquationTypeScalar = MagFieldEquation<Field_Type_Scalar>; // New equation ...
  // auto gvEquationScalar      = new GvEquationTypeScalar(gvFieldScalar);
  cout << "Using Scalar CMS Field type' . " << endl;  
#else
#ifdef NEW_SCALAR_FIELD
// #define gvFieldScalar gvField;
  auto gvFieldScalar         = new UniformMagField(fieldValueVec);
  cout << "Using new Scalar Field types and 'MagFieldEquation' . " << endl;  
#else
  // If we plan to test against 'plain' scalar objects: field, equation, stepper, ...
  auto fieldValueVec = geant::units::tesla * ThreeVector_d(x_field, y_field, z_field);
  auto gvFieldScalar = new ScalarUniformMagField(fieldValueVec);
  //                      new Field_Type_Scalar( fieldValueVec );
  cout << "Using old Scalar Field types and ScalarMagFieldEquation. " << endl;
#endif
#endif  

#ifdef NEW_SCALAR_FIELD
  using GvScalarEquationType = MagFieldEquation<Field_Type_Scalar>; // New scalar
#else
  using GvScalarEquationType = ScalarMagFieldEquation<Field_Type_Scalar, Nposmom>;
#endif
  
  auto gvEquationScalar = new GvScalarEquationType(gvFieldScalar); // Same as vector, yes
  //                           ^was MagFieldEquation

#ifdef NEW_SCALAR_FIELD  
  auto myStepperScalar = CashKarp<GvEquationTypeScalar, Nposmom>(gvEquationScalar);
  cout << "Using new Scalar Stepper types - i.e. flexible CashKarp " << endl;
#else
  VScalarIntegrationStepper *myStepperScalar;
  myStepperScalar = StepperFactory::CreateStepper<GvScalarEquationType>(gvEquationScalar, stepper_no);
  cout << "Scalar> Using old Scalar Stepper types - i.e. Scalar stepper from factory.  Stepper type #= " << stepper_no << endl;
#endif
  
  int statisticsVerbosity = 1;

  auto refScalarDriver = new ScalarIntegrationDriver(hminimum, myStepperScalar, epsTol, 
                                                     Nposmom, statisticsVerbosity);
  cout << "Scalar> Driver created:  Type= ScalarIntegrationDriver.  hminimum= " << hminimum << " , Nvar = " << Nposmom << " , statsVerb= " << statisticsVerbosity << endl;
  //==========  Scalar Driver prepared =========================

  // -- Call internal Test code of Driver
  // int numTracks;
  #if 0
    bool ok= vectorDriver->TestInitializeLanes<Double_v>();
    if( !ok ) {
       std::cerr << "WARNING> vectorDriver<Double_v>->TestInitializeLanes() returned " << ok << std::endl;
    }
  #endif

  // double total_step = 0.;

  Bool_v goodAdvance(true);

  // Double_v  hStep1( 10.0 );
  // goodAdvance = testDriver->AccurateAdvance( yTrackIn, hStep1, epsTol, yTrackOut );

  constexpr int nTracks = 16;
  FieldTrack yInput[nTracks], yOutput[nTracks];

  // std::fill_n(hstep, nTracks, 20);
  // double hstep[nTracks] = {  0.001,  0.01,    0.1,    1.00,   2.0,    4.0,     8.0,     16.0,
  //                             32.0,   64.0,   200.0,  600.,  2000.0, 10000.,  30000.0,  100000.0};
  // double hstep[nTracks] = {11.0,       20.0,  33.0,  44.44, 55.555, 66.6666, 77.77777, 88.888888,  
  double hstep[nTracks] = {  0.001,  0.01,    0.1,    1.00,   2.0,    4.0,     8.0,     16.0,  
                         99.9999999, 101.0, 111.0, 122.0, 133.0,  144.0,   155.0,    166.0};  
  // double hstep[] = { 11.0, 200.0, 330.0, 444.40, 555.55, 666.666, 777.7777, 888.88888, 999.999999, 100.0, 11.0, 12.0,
  // 13.0, 14.0, 15.3, 16.9 };

  double charge[nTracks] = {1.0, -2.0,  3.0,  -4.0,  5.0,  -6.0,  7.0, -8.0,
                            9.0, -10.0, 11.0, -12.0, 13.0, -14.0, 15., -16.};
  // = {0, 0, 0, 1, -.3, .4, 20, 178., 920.};
  bool succeeded[nTracks];

  bool firstTrial = true;
#ifdef MAINTESTING
  firstTrial = false; // Just go to the loop of calls.
#endif

  if (firstTrial) {
    double x_pos = 20., // pos - position  : input unit = mm
        y_pos = 20., z_pos = 20.;
    double x_mom       = 0.; // mom - momentum  : input unit = GeV / c
    double y_mom       = 1.;
    double z_mom       = 1.;
    const double mmGVf = geant::units::millimeter;
    const double ppGVf = geant::units::GeV;

    double posMom[] = {x_pos * mmGVf, y_pos * mmGVf, z_pos * mmGVf, x_mom * ppGVf, y_mom * ppGVf, z_mom * ppGVf};
    // double posMom[] = {0.0513401, 0.095223, 0.0916195, 0.635712, 0.717297, 0.141603 };

    for (int j = 0; j < nTracks; ++j) {
      yInput[j].LoadFromArray(posMom);
      yInput[j].ResetCurveLength();
      yOutput[j].LoadFromArray(posMom); // Not strictly needed
    }

#ifdef CHECK_ONE_LANE
// If compiled with     
    int trackToCheck = ( trackToPrint < nTracks ) ? trackToPrint : -1; 
    // vectorDriver->CheckTrack( trackToCheck );
    IntegrationDriverConstants::GetInstance()->SetTrackToCheck(trackToCheck);
#endif
  
    vectorDriver->AccurateAdvance<Double_v>(yInput, hstep, charge, /* epsTol, */ yOutput, succeeded, nTracks);
    // ==========================

    // if( verbose ) cout << "-- Vector Driver done (advanced)." << endl;

    const ThreeVector_d startPosition(posMom[0], posMom[1], posMom[2]);
    const ThreeVector_d startMomentum(posMom[3], posMom[4], posMom[5]);

    ScalarFieldTrack yScalTrackIn(startPosition, startMomentum, charge[0]);  // constant !! yStart
    ScalarFieldTrack yScalTrackOut(startPosition, startMomentum, charge[0]); // yStart
    
#if 0    
    refScalarDriver->SetPrintDerived(false);
    bool scalarResult = refScalarDriver->AccurateAdvance(yScalTrackIn, hstep[0], /* epsTol, */ yScalTrackOut);

    if( verbose )
    {
       cout << "-- Scalar Driver done (advanced)." << endl;
       cout << " yOutput is   : " << yOutput[0] << " for yInput: " << yInput[0] << endl;
       cout << " yTrackOut is : " << yScalTrackOut << " for yTrackIn: " << yScalTrackIn << endl;
       cout << " Success of Vector: " << succeeded[0] << endl;
       cout << " Success of scalar: " << scalarResult << endl;
    }
#endif

    for (int i = 0; i < nTracks; ++i)
    {
      yScalTrackIn= ScalarFieldTrack (startPosition, startMomentum, charge[i]);
      yScalTrackOut= ScalarFieldTrack(startPosition, startMomentum, charge[i]); // yStart

#ifdef CHECK_ONE_LANE      
      refScalarDriver->SetPrintDerived(i==trackToCheck);
      refScalarDriver->SetTrackNumber(i);  // For info only
#endif      
      refScalarDriver->AccurateAdvance(yScalTrackIn, hstep[i], /* epsTol, */  yScalTrackOut);
      // ************  ***************

      // bool diffFound;
      double maxRelDiff= 0.0; // Max (Absolute) relative diff
      // int dfLane =
         CompareSixVectors( yOutput[i], yScalTrackOut, "yOut/x,p - (Vector vs. Scalar)", i, hstep[i],
                                      yInput[i],  yScalTrackIn, maxRelDiff );
         std::cout << "testVecIntDrv -  hDid :  vector = " << yOutput[i].GetCurveLength()
                   << " scalar  = " << yScalTrackOut.GetCurveLength() << std::endl;
      // diffFound= dfLane > 0;
    }
  }
// #define CALCULATETIME

#ifdef MAINTESTING
#ifdef CALCULATETIME
  std::vector<double> ratioVector;
#endif

  srand(time(NULL));

  for (int step = 0; step < no_of_steps; ++step)
  // for (int step = 0; step < 20; ++step)
  {
    total_step += step_len;
    std::fill_n(hstep, nTracks, total_step);

    // srand(9);

    x_pos           = (float)rand() / (RAND_MAX);
    y_pos           = (float)rand() / (RAND_MAX);
    z_pos           = (float)rand() / (RAND_MAX);
    x_mom           = (float)rand() / (RAND_MAX);
    y_mom           = (float)rand() / (RAND_MAX);
    z_mom           = (float)rand() / (RAND_MAX);
    double posMom[] = {x_pos * mmGVf, y_pos * mmGVf, z_pos * mmGVf, x_mom * ppGVf, y_mom * ppGVf, z_mom * ppGVf};

    const ThreeVector_d startPosition(posMom[0], posMom[1], posMom[2]);
    const ThreeVector_d startMomentum(posMom[3], posMom[4], posMom[5]);
    ScalarFieldTrack yTrackIn(startPosition, startMomentum);  // yStart
    ScalarFieldTrack yTrackOut(startPosition, startMomentum); // yStart

    for (int j = 0; j < nTracks; ++j) {
      yInput[j].LoadFromArray(posMom);
      yOutput[j].LoadFromArray(posMom);
    }

    /*      double X_Pos[nTracks], Y_Pos[nTracks], Z_Pos[nTracks];
          double X_Mom[nTracks], Y_Mom[nTracks], Z_Mom[nTracks];
          double posMomMatrix[nTracks][6];
          for (int i = 0; i < nTracks; ++i)
          {
            X_Pos[i] = (float) rand()/(RAND_MAX) ;
            Y_Pos[i] = (float) rand()/(RAND_MAX) ;
            Z_Pos[i] = (float) rand()/(RAND_MAX) ;
            X_Mom[i] = (float) rand()/(RAND_MAX) ;
            Y_Mom[i] = (float) rand()/(RAND_MAX) ;
            Z_Mom[i] = (float) rand()/(RAND_MAX) ;

            // posMomMatrix[i] = {X_Pos[i] * mmGVf, Y_Pos[i] * mmGVf ,Z_Pos[i] * mmGVf,
            //                    X_Mom[i] * ppGVf ,Y_Mom[i] * ppGVf ,Z_Mom[i] * ppGVf};
            posMomMatrix[i][0] = X_Pos[i] * mmGVf;
            posMomMatrix[i][1] = Y_Pos[i] * mmGVf;
            posMomMatrix[i][2] = Z_Pos[i] * mmGVf;
            posMomMatrix[i][3] = X_Mom[i] * ppGVf;
            posMomMatrix[i][4] = Y_Mom[i] * ppGVf;
            posMomMatrix[i][5] = Z_Mom[i] * ppGVf;
            yInput [i].LoadFromArray(posMomMatrix[i]);
            yOutput[i].LoadFromArray(posMomMatrix[i]);
          }
    */
    for (int i = 1; i < nTracks; ++i) {
      hstep[i] = hstep[i] + i * hstep[i];
    }
// #define DebuggingSection
#ifdef CALCULATETIME
    clock_t biasClock  = clock();
    biasClock          = clock() - biasClock;
    clock_t biasClock2 = clock();
    biasClock2         = clock() - biasClock2;
    cout << "biasClock2 is: " << biasClock2 << " and 1 is : " << biasClock << endl;
    cout << " and diff. in float is : " << ((float)(biasClock - biasClock2)) / CLOCKS_PER_SEC;

    clock_t clock1 = clock();
#endif

#ifndef DebuggingSection
    vectorDriver->AccurateAdvance<Double_v>(yInput, hstep, /* epsTol, */ yOutput, nTracks, succeeded);
#endif

#ifdef CALCULATETIME
    clock1              = clock() - clock1 - biasClock;
    float clock1InFloat = ((float)clock1) / CLOCKS_PER_SEC;
    cout << "Vector time is: " << clock1InFloat << endl;
#endif
//      refScalarDriver->AccurateAdvance( yTrackIn, hstep[0], /*epsTol,*/ yTrackOut );

      cout<<" yOutput[0] is: "<< yOutput[0]<<" for yInput: "  <<yInput[0]<< endl;
      cout<<" yTrackOut is: " << yTrackOut <<" for yTrackIn: "<<yTrackIn << endl;*/

/*      for (int i = 0; i < nTracks; ++i)
      {
        const ThreeVector_d  startPosition( posMomMatrix[i][0], posMomMatrix[i][1], posMomMatrix[i][2]);
        const ThreeVector_d  startMomentum( posMomMatrix[i][3], posMomMatrix[i][4], posMomMatrix[i][5]);
        ScalarFieldTrack yTrackIn ( startPosition, startMomentum );
        ScalarFieldTrack yTrackOut( startPosition, startMomentum );

        refScalarDriver->AccurateAdvance( yTrackIn, hstep[i], // epsTol, 
            yTrackOut );

        cout<<" yOutput["<<i<<"] is: "<< yOutput[i]<<" for yInput: "  <<yInput[i]<< endl;
        cout<<" yTrackOut is: " << yTrackOut <<" for yTrackIn: "<<yTrackIn << endl;
      }
*/
#ifdef CALCULATETIME
    clock_t clock2 = clock();
#endif
    for (int i = 0; i < nTracks; ++i) {
#ifndef DebuggingSection
      refScalarDriver->AccurateAdvance(yTrackIn, hstep[i], /*epsTol,*/ yTrackOut);

      cout << " yOutput[" << i << "] is: " << yOutput[i] << " for yInput: " << yInput[i] << endl;
      cout << " yTrackOut is : " << yTrackOut << " for yTrackIn: " << yTrackIn << " for hstep: " << hstep[i] << endl;
#endif
    }
#ifdef CALCULATETIME
    clock2              = clock() - clock2 - biasClock;
    float clock2InFloat = ((float)clock2) / CLOCKS_PER_SEC;
    cout << "scalar time is: " << clock2InFloat << endl;
    // cout<<"ratio is: "<<clock2InFloat/clock1InFloat<<endl;
    ratioVector.push_back(clock2InFloat / clock1InFloat);
#endif
  }

#ifdef CALCULATETIME
  int sizeOfRatioVector = ratioVector.size(); // no_steps;
  cout << "Size of ratioVector is: " << ratioVector.size() << endl;
  cout << "Time ratios are: " << endl;
  for (int i = 0; i < sizeOfRatioVector; ++i) {
    cout << i << " " << ratioVector.back() << " " << endl;
    ratioVector.pop_back();
  }
  cout << endl;
#endif

#endif

#ifdef DebuggingSection
  if (true) {
    // double posMomt[] = { 0.0394383, 0.0783099, 0.079844, 0.911647, 0.197551, 0.335223};
    // double posMomt[] = { 0.018234, 0.0576585, 0.0694335, 0.857714, 0.205414, 0.186222 };
    // double posMomt[] = {0.0627625, 0.0184593, 0.0449569, 0.68876, 0.163077, 0.841872};
    // double posMomt[] = {0.0614384, 0.0736116, 0.0124955, 0.236072, 0.737118, 0.0545562};
    double posMomt[] = {0.00668596, 0.0106866, 0.000127922, 0.126255, 0.101243, 0.384278};
    // double posMomt[] = {0.0513401, 0.095223, 0.0916195, 0.635712, 0.717297, 0.141603 };

    total_step = 120;
    std::fill_n(hstep, nTracks, total_step);

    for (int i = 0; i < nTracks; ++i) {
      hstep[i] = hstep[i] + i * hstep[i];
      // hstep[i] = (float) rand()/(RAND_MAX) *200.;
    }

    for (int j = 0; j < nTracks; ++j) {
      yInput[j].LoadFromArray(posMomt);
      yOutput[j].LoadFromArray(posMomt);
    }

    const ThreeVector_d startPosition(posMomt[0], posMomt[1], posMomt[2]);
    const ThreeVector_d startMomentum(posMomt[3], posMomt[4], posMomt[5]);
    ScalarFieldTrack yTrackIn(startPosition, startMomentum, charge);
    ScalarFieldTrack yTrackOut(startPosition, startMomentum, charge);

    vectorDriver->AccurateAdvance(yInput, charge, hstep, /*epsTol,*/ yOutput, nTracks, succeeded);
    // refScalarDriver->AccurateAdvance( yTrackIn, hstep[11], /* epsTol,*/ yTrackOut );

    for (int i = 0; i < nTracks; ++i) {
      refScalarDriver->AccurateAdvance(yTrackIn[i], hstep[i], /* epsTol, */ yTrackOutScalar[i]);

      // cout<<" yOutput["<<i<<"] is: "<< yOutput[i]<<" for hstep: "<<hstep[i]<< endl ;
      cout << " yOutput[" << i << "] is: " << yOutput[i] << " for hstep: " << hstep[i] << endl;
      cout << " yTrackOut is : " << yTrackOut << endl;
    }

    // ScalarFieldTrack randNew = yTrackIn + yTrackOut;

    // cout<<yOutput[0]<<endl;
    // cout<<yOutput[1]<<endl;
    // cout<<yOutput[2]<<hstep[2]<<endl;
  }
#endif

  cout << " Scalar Stepper function calls are: " << refScalarDriver->fStepperCalls << " and OneGoodStep calls are "
       << refScalarDriver->fNoTotalSteps << endl;
  cout << " Vector Stepper function calls are: " << vectorDriver->GetNumberOfStepperCalls() << " and Driver (OneGoodStep/KeepStepping) steps are "
       << vectorDriver->GetNumberOfTotalSteps() << endl;

  //========================End testing IntegrationDriver=======================

  /*------ Clean up ------*/
  delete myStepper;
  delete gvField;

  // deleting IntegrationDriver
  delete vectorDriver;
  delete refScalarDriver;

  cout << "\n\n#-------------End of output-----------------\n";
}
