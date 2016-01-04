//
//  Compare the output of different steppers
// 
//  Based on testStepperFixed.cc
//    was the work of Somnath Banerjee in GSoC 2015
//
#include <iomanip>
#include <ctime>

#include "Units.h"

// using fieldUnits::meter;
using fieldUnits::millimeter;   
using fieldUnits::second;  
using fieldUnits::eplus;  
using fieldUnits::tesla;
using fieldUnits::degree;

#include <Vc/Vc>
#include "base/Vector3D.h"

#include "TUniformMagField.h"
#include "TMagFieldEquation.h"
#include "GUVIntegrationStepper.h"
#include "StepperFactory.h"
#include "GUFieldTrack.h"
#include "GUIntegrationDriver.h"

#include "TemplateTUniformMagField.h"
#include "TemplateTMagFieldEquation.h"
#include "TemplateFieldEquationFactory.h"
#include "TemplateGUVIntegrationStepper.h"
#include "TemplateGUTCashKarpRKF45.h"
#include "TemplateGUIntegrationDriver.h"
#include "FieldTrack.h"

#include "TemplateCMSmagField.h"
#include "ScalarCMSmagField.h"

#include <stdlib.h>

using namespace std;

int main(/*int argc, char *args[]*/)
{
    constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

    using Backend = vecgeom::kVc ;
    // typedef typename Backend::precision_v Double;
    // typedef vecgeom::Vector3D<Double> ThreeVectorSimd;
    typedef vecgeom::Vector3D<double> ThreeVector_d;
  // #define USECMSFIELD
  #ifdef USECMSFIELD
    using Field_Type        = TemplateCMSmagField<Backend>;
    using Field_Type_Scalar = ScalarCMSmagField;
    // using Field_Type_Scalar = TemplateCMSmagField<vecgeom::kScalar>;
  #else
    using Field_Type_Scalar = TUniformMagField;
    using Field_Type        = TemplateTUniformMagField<Backend>;
  #endif 

    using GvEquationType    = TemplateTMagFieldEquation<Backend, Field_Type, Nposmom>;
   
    /* -----------------------------SETTINGS-------------------------------- */
    
    /* Parameters of test
     - Modify values  */
    
    // int no_of_steps = 20;         // No. of Steps for the stepper
    int stepper_no  =  5;         // Choose stepper no., for refernce see above
    double step_len_mm = 200.;    // meant as millimeter;  //Step length 
    double z_field_in  = DBL_MAX;
    
    //Checking for command line values :
/*    if(argc>1)
        stepper_no = atoi(args[1]);
    if(argc > 2)
       step_len_mm = (float)(stof(args[2]));   // *mm);
    if(argc > 3)
        no_of_steps = atoi(args[3]);
    if(argc > 4)
       z_field_in = (float) (stof(args[4]));     // tesla*/
    // double step_len = step_len_mm * fieldUnits::millimeter;
    
    //Set Charge etc.
    double particleCharge = +1.0;      // in e+ units
    
    //Choice of output coordinates
/*    int
    columns[] =
    {
       1 , 1 , 1 ,  // position  x, y, z 
       1 , 1 , 1 ,  // momentum  x, y, z
       0 , 0 , 0 ,  // dydx pos  x, y, z
       1 , 1 , 0    // dydx mom  x, y, z
    };*/ //Variables in yOut[] & dydx[] we want to display - 0 for No, 1 for yes

    // const unsigned int nwdf= 12;  // Width for float/double
    
    //Set coordinates here
/*    double
       x_pos = 0.,                 //pos - position  : input unit = mm
       y_pos = 0.,
       z_pos = 0.;
    double   
       x_mom = 0.,                 //mom - momentum  : input unit = GeV / c
       y_mom = 1.,
       z_mom = 1.;*/
    double          
       x_field = 0.,               //Uniform Magnetic Field (x,y,z)
       y_field = 0.;               //  Tesla // *tesla ;
    double z_field;
    // vecgeom::MaskedAssign(z_field_in < DBL_MAX, z_field_in, &z_field);

    if( z_field_in < DBL_MAX )
       z_field = z_field_in;
    else
       z_field = -1.0;  //  Tesla // *tesla ;

    // Field
  #ifdef USECMSFIELD
    auto gvField = new Field_Type ("../VecMagFieldRoutine/cmsmagfield2015.txt");
  #else 
    auto gvField = new Field_Type( fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field) );
  #endif


    cout << "#  Initial  Field strength (GeantV) = "
         << x_field << " , " << y_field << " , " << z_field 
       // << (1.0/fieldUnits::tesla) * gvField->GetValue()->X() << ",  "
         << " Tesla " << endl;
    // cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;
    //Create an Equation :
 
    auto gvEquation =
       TemplateFieldEquationFactory<Backend>::CreateMagEquation<Field_Type >(gvField);

    // gvEquation->InitializeCharge( particleCharge );  // Send it via Stepper instead    
    

/*    double yIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                    x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};*/

    std::cout << "# step_len_mm = " << step_len_mm;
    
    
    //Empty buckets for results
/*    Double dydx[8] = {0.,0.,0.,0.,0.,0.,0.,0.},  // 2 extra safety buffer
           yout[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
           yerr[8] = {0.,0.,0.,0.,0.,0.,0.,0.};*/
  

    //=======================Test part for Integration driver====================
    double hminimum = 0.2;
    // auto testVectorDriver = new TemplateGUIntegrationDriver<Backend>(hminimum, myStepper);

    //========= Preparing scalar Integration Driver ============
    using  GvEquationTypeScalar=  TMagFieldEquation<Field_Type_Scalar, Nposmom>;
    
  #ifdef USECMSFIELD
    auto gvFieldScalar    = new Field_Type_Scalar("../VecMagFieldRoutine/cmsmagfield2015.txt");
  #else
    auto gvFieldScalar    = new Field_Type_Scalar( fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field) );
  #endif
    auto gvEquationScalar = new GvEquationTypeScalar(gvFieldScalar);

    GUVIntegrationStepper *myStepperScalar; 
    myStepperScalar= StepperFactory::CreateStepper<GvEquationTypeScalar>(gvEquationScalar, stepper_no);

    int statisticsVerbosity = 1;
    auto testScalarDriver= new GUIntegrationDriver( hminimum,
                                                    myStepperScalar,
                                                    Nposmom,
                                                    statisticsVerbosity); 
    testScalarDriver->InitializeCharge( particleCharge );
    //==========  Scalar Driver prepared =========================

    auto myStepper = new TemplateGUTCashKarpRKF45<Backend,GvEquationType,Nposmom>(gvEquation);

    // myStepper->InitializeCharge( particleCharge );

    // auto testVectorDriver = new TemplateGUIntegrationDriver<Backend>(hminimum, myStepper, myStepperScalar);
    auto testVectorDriver = new TemplateGUIntegrationDriver<Backend>(hminimum, myStepper, myStepperScalar, testScalarDriver);

    bool debugValue ; 
    cout<< "Debug? " << endl;
    cin >> debugValue;
    testVectorDriver->SetPartDebug(debugValue);
    // ========== Vector Driver prepared ========================


    double total_step = 0.;

    typedef typename Backend::bool_v Bool;
    Bool goodAdvance(true);
    double epsTol = 1.0e-5;

    // goodAdvance = testDriver->AccurateAdvance( yTrackIn, total_step, epsTol, yTrackOut );

    constexpr int nTracks = 16;
    FieldTrack yInput[nTracks], yOutput[nTracks];
    // double posMom[] ={0., 0., 0., 0., 1., 1.};

    double hstep[nTracks] = { 11.0, 2.0, 33.0, 4.0, 55.0, 6.0, 77.0, 8.0, 9.0, 100.0, 11.0, 12.0, 13.0, 14.0, 15.3, 16.9 };
    double charge[nTracks] = { 1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0, -10.0, 11.0, -12.0, 13.0, -14.0, 15., -16. };    
        // = {0, 0, 0, 1, -.3, .4, 20, 178., 920.}; 
    bool   succeeded[nTracks];

#define DebuggingSection
// #define CALCULATETIME 

#ifndef DebuggingSection
#ifndef MAINTESTING    
    double
      x_pos = 20.,                 //pos - position  : input unit = mm
      y_pos = 20.,
      z_pos = 20.;  */
    //  double x_mom = 0. ;               //mom - momentum  : input unit = GeV / c
    //  double y_mom = 1. ;
    //  double z_mom = 1. ;
    // const double mmGVf = fieldUnits::millimeter;
    // const double ppGVf = fieldUnits::GeV ; 

    double posMom[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                       x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};
    // double posMom[] = {0.0513401, 0.095223, 0.0916195, 0.635712, 0.717297, 0.141603 };

    std::fill_n(hstep, nTracks, 20);
    const ThreeVector_d  startPosition( posMom[0], posMom[1], posMom[2]);
    const ThreeVector_d  startMomentum( posMom[3], posMom[4], posMom[5]);
    GUFieldTrack yTrackIn ( startPosition, startMomentum );  // yStart
    GUFieldTrack yTrackOut( startPosition, startMomentum );  // yStart

    for (int j = 0; j < nTracks; ++j)
    {
      yInput [j].LoadFromArray(posMom);
      yOutput[j].LoadFromArray(posMom);
    }

    testVectorDriver->AccurateAdvance( yInput,   hstep,    epsTol, yOutput, nTracks, succeeded );
    bool scalarResult = testScalarDriver->AccurateAdvance( yTrackIn, hstep[0], epsTol, yTrackOut );
    cout<<" yOutput is   : "<< yOutput[0]<<" for yInput: "  <<yInput[0]<< endl;
    cout<<" yTrackOut is : "<< yTrackOut <<" for yTrackIn: "<<yTrackIn << endl;
    cout<<" Success of Vector: "<< succeeded[0] << endl;
    cout<<" Success of scalar: "<< scalarResult << endl;
#endif 
#endif 

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

      x_pos = (float) rand()/(RAND_MAX) ;
      y_pos = (float) rand()/(RAND_MAX) ;
      z_pos = (float) rand()/(RAND_MAX) ;
      x_mom = (float) rand()/(RAND_MAX) ;
      y_mom = (float) rand()/(RAND_MAX) ;
      z_mom = (float) rand()/(RAND_MAX) ;
      double posMom[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                         x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};

      const ThreeVector_d  startPosition( posMom[0], posMom[1], posMom[2]);
      const ThreeVector_d  startMomentum( posMom[3], posMom[4], posMom[5]);
      GUFieldTrack yTrackIn ( startPosition, startMomentum );  // yStart
      GUFieldTrack yTrackOut( startPosition, startMomentum );  // yStart

      for (int j = 0; j < nTracks; ++j)
      {
        yInput [j].LoadFromArray(posMom);
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
      for (int i = 1; i < nTracks; ++i)
      {
         hstep[i] = hstep[i] + i*hstep[i];
      }
// #define DebuggingSection
#ifdef CALCULATETIME
      clock_t biasClock= clock();
      biasClock = clock() - biasClock;
      clock_t biasClock2 = clock();
      biasClock2 = clock() - biasClock2;
      cout<<"biasClock2 is: "<<biasClock2<<" and 1 is : "<<biasClock<<endl;
      cout<<" and diff. in float is : "<< ((float)(biasClock-biasClock2))/CLOCKS_PER_SEC;

      clock_t clock1= clock();
#endif 

#ifndef DebuggingSection
      testVectorDriver->AccurateAdvance( yInput, hstep, epsTol, yOutput, nTracks, succeeded );
#endif 

#ifdef CALCULATETIME
      clock1 = clock() - clock1 - biasClock;
      float clock1InFloat = ((float)clock1)/CLOCKS_PER_SEC;
      cout<<"Vector time is: "<<clock1InFloat<<endl;
#endif 
/*      testScalarDriver->AccurateAdvance( yTrackIn, hstep[0], epsTol, yTrackOut );

      cout<<" yOutput[0] is: "<< yOutput[0]<<" for yInput: "  <<yInput[0]<< endl;
      cout<<" yTrackOut is: " << yTrackOut <<" for yTrackIn: "<<yTrackIn << endl;*/

/*      for (int i = 0; i < nTracks; ++i)
      {
        const ThreeVector_d  startPosition( posMomMatrix[i][0], posMomMatrix[i][1], posMomMatrix[i][2]);
        const ThreeVector_d  startMomentum( posMomMatrix[i][3], posMomMatrix[i][4], posMomMatrix[i][5]);
        GUFieldTrack yTrackIn ( startPosition, startMomentum ); 
        GUFieldTrack yTrackOut( startPosition, startMomentum ); 

        testScalarDriver->AccurateAdvance( yTrackIn, hstep[i], epsTol, yTrackOut );

        cout<<" yOutput["<<i<<"] is: "<< yOutput[i]<<" for yInput: "  <<yInput[i]<< endl;
        cout<<" yTrackOut is: " << yTrackOut <<" for yTrackIn: "<<yTrackIn << endl;
      }
*/
#ifdef CALCULATETIME
      clock_t clock2= clock();
#endif 
      for (int i = 0; i < nTracks; ++i)
      {
#ifndef DebuggingSection
        testScalarDriver->AccurateAdvance( yTrackIn, hstep[i], epsTol, yTrackOut );

        cout<<" yOutput["<<i<<"] is: "<< yOutput[i]<<" for yInput: "  <<yInput[i]<< endl;
        cout<<" yTrackOut is : "      << yTrackOut <<" for yTrackIn: "<<yTrackIn <<" for hstep: "<<hstep[i]<< endl;
#endif 
      }
#ifdef CALCULATETIME
      clock2 = clock() - clock2 - biasClock;
      float clock2InFloat = ((float)clock2)/CLOCKS_PER_SEC;
      cout<<"scalar time is: "<<clock2InFloat<<endl;
      // cout<<"ratio is: "<<clock2InFloat/clock1InFloat<<endl;
      ratioVector.push_back(clock2InFloat/clock1InFloat);
#endif 
    }

#ifdef CALCULATETIME
    int sizeOfRatioVector = ratioVector.size(); //no_steps;
    cout<<"Size of ratioVector is: "<<ratioVector.size()<<endl;
    cout<<"Time ratios are: "<<endl;
    for (int i = 0; i < sizeOfRatioVector; ++i)
    {
      cout<<i<<" "<<ratioVector.back()<< " " <<endl;
      ratioVector.pop_back();
    }
    cout<<endl;
#endif 

#endif 

#ifdef DebuggingSection
    if(true){
      // double posMomt[] = { 0.0394383, 0.0783099, 0.079844, 0.911647, 0.197551, 0.335223};  
      // double posMomt[] = { 0.018234, 0.0576585, 0.0694335, 0.857714, 0.205414, 0.186222 };
      // double posMomt[] = {0.0627625, 0.0184593, 0.0449569, 0.68876, 0.163077, 0.841872};
      // double posMomt[] = {0.0614384, 0.0736116, 0.0124955, 0.236072, 0.737118, 0.0545562};
      double posMomt[] = {0.00668596, 0.0106866, 0.000127922, 0.126255, 0.101243, 0.384278};
      // double posMomt[] = {0.0513401, 0.095223, 0.0916195, 0.635712, 0.717297, 0.141603 };

      
      total_step = 120;
      std::fill_n(hstep, nTracks, total_step);

      for (int i = 0; i < nTracks; ++i)
      {
         hstep[i] = hstep[i] + i*hstep[i];
         // hstep[i] = (float) rand()/(RAND_MAX) *200.; 
      }

      for (int j = 0; j < nTracks; ++j)
      {
        yInput [j].LoadFromArray(posMomt);
        yOutput[j].LoadFromArray(posMomt);
      }

      const ThreeVector_d  startPosition( posMomt[0], posMomt[1], posMomt[2]);
      const ThreeVector_d  startMomentum( posMomt[3], posMomt[4], posMomt[5]);
      GUFieldTrack yTrackIn ( startPosition, startMomentum );  // yStart
      GUFieldTrack yTrackOut( startPosition, startMomentum );  // yStart

      testVectorDriver->AccurateAdvance( yInput,
                                         charge,
                                         hstep,
                                         epsTol,
                                         yOutput,
                                         nTracks,
                                         succeeded );
      // testScalarDriver->AccurateAdvance( yTrackIn, hstep[11], epsTol, yTrackOut );

      for (int i = 0; i < nTracks; ++i)
      {
        // testScalarDriver->AccurateAdvance( yTrackIn, hstep[i], epsTol, yTrackOut );

        // cout<<" yOutput["<<i<<"] is: "<< yOutput[i]<<" for hstep: "<<hstep[i]<< endl ;
        cout<<" yOutput["<<i<<"] is: "<< yOutput[i]<<" for hstep: "<<hstep[i]<< endl ;
        cout<<" yTrackOut is : "      << yTrackOut << endl;
      }

      // GUFieldTrack randNew = yTrackIn + yTrackOut;

      // cout<<yOutput[0]<<endl;
      // cout<<yOutput[1]<<endl;
      // cout<<yOutput[2]<<hstep[2]<<endl;

    }
#endif 

    cout<<" Scalar Stepper function calls are: "<< testScalarDriver->fStepperCalls <<" and OneGoodStep calls are "<<testScalarDriver->fNoTotalSteps << endl;
    cout<<" Vector Stepper function calls are: "<< testVectorDriver->fStepperCalls <<" and OneStep calls are "<<testVectorDriver->fNoTotalSteps << endl;


    //========================End testing IntegrationDriver=======================

    /*------ Clean up ------*/
    delete myStepper; 
    delete gvField;

    // deleting IntegrationDriver
    delete testVectorDriver;
    delete testScalarDriver;      
    
    cout<<"\n\n#-------------End of output-----------------\n";
}
