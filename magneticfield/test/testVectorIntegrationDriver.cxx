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

#include <base/Vector3D.h>
#include <Geant/VectorTypes.h>

#include "UniformMagField.h"
#include "MagFieldEquation.h"

// #include "IntegrationStepper.h"

// #include "TemplateTMagFieldEquation.h"
// #include "TemplateFieldEquationFactory.h"
// #include "TemplateGUVIntegrationStepper.h"

#include "CashKarp.h"
#include "SimpleIntegrationDriver.h"

#include "FieldTrack.h"
// #include "TemplateFieldTrack.h"

// #define  NEW_SCALAR_FIELD 1

// #define USECMSFIELD
#ifdef   USECMSFIELD
#include "TemplateCMSmagField.h"
#include "ScalarCMSmagField.h"
#else
  #ifndef NEW_SCALAR_FIELD
  //  Transition measure --- compare to old Scalar field types 2017.11.16
    #include "TUniformMagField.h"
    #include "TMagFieldEquation.h"
    #include "StepperFactory.h"
    #include "ScalarFieldTrack.h"
    #include "ScalarIntegrationDriver.h"
  #endif
#endif

#include <stdlib.h>

// using namespace std;
using std::cout;
using std::cerr;
using std::endl;

using Double_v = Geant::Double_v;
using Bool_v   = vecCore::Mask_v<Double_v>;

int main(/*int argc, char *args[]*/)
{
    constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec
    // template <typename T> using Vector3D = vecgeom::Vector3D<T>;
    using ThreeVector_d = vecgeom::Vector3D<double>;
    
  #ifdef USECMSFIELD
    using Field_Type        = FlexibleCMSmagField; // TO-DO 
    using Field_Type_Scalar = ScalarCMSmagField;
    // using Field_Type_Scalar = TemplateCMSmagField<vecgeom::kScalar>;
  #else
    using Field_Type        = UniformMagField;  // TemplateTUniformMagField<Backend>;
  #ifdef NEW_SCALAR_FIELD
    // New types ... under development  2017.11.16
    using Field_Type_Scalar = UniformMagField;
  #else
    using Field_Type_Scalar = TUniformMagField;
  #endif
  #endif 

    using GvEquationType    = MagFieldEquation<Field_Type>; // , Nposmom>;
   
    /* -----------------------------SETTINGS-------------------------------- */
    
    /* Parameters of test - Modify values  */
    
    // int no_of_steps = 20;         // No. of Steps for the stepper
    int    stepper_no  =  5;         // Choose stepper no., for refernce see above
    double step_len_mm = 200.;    // meant as millimeter;  //Step length 
    double z_field_in  = DBL_MAX;
    
    double x_field = 0.,               //Uniform Magnetic Field (x,y,z)
           y_field = 0.,               //  Tesla // *tesla ;
           z_field;

    if( z_field_in < DBL_MAX ) z_field = z_field_in;
    else                       z_field = -1.0;  //  * Tesla

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
    auto gvEquation = new MagFieldEquation<Field_Type>(gvField);
    // TemplateFieldEquationFactory<Double_v>::CreateMagEquation<Field_Type >(gvField);

    cout << " Equation object created - address = " << gvEquation << endl;
    cout << "# step_len_mm = " << step_len_mm << endl;
    
    //=======================Test part for Integration driver====================
    double hminimum = 0.2;

    //==========  Vector Driver: start preparation ===============
    using StepperType = CashKarp<GvEquationType,Nposmom>;
    auto myStepper = new StepperType(gvEquation);
       // new CashKarp<GvEquationType,Nposmom>(gvEquation);
       // new TemplateGUTCashKarpRKF45<Double_v,GvEquationType,Nposmom>(gvEquation);       

    // myStepper->InitializeCharge( particleCharge );
    int statsVerbose=1;

    using  DriverType = SimpleIntegrationDriver<StepperType,Nposmom>;
    auto vectorDriver =
       // new SimpleIntegrationDriver<StepperType,Nposmom>
       new DriverType
                                                        (hminimum,
                                                         myStepper,
                                                         Nposmom,
                                                         statsVerbose);
    cout << " Vector Driver created." << endl;
    // ========== Vector Driver prepared ========================
    
    //========= Preparing scalar Integration Driver ============
#ifdef USECMSFIELD
    auto gvFieldScalar    = new Field_Type_Scalar("../VecMagFieldRoutine/cmsmagfield2015.txt");
    using  GvEquationTypeScalar=  MagFieldEquation<Field_Type_Scalar>;   // New equation ...
    auto gvEquationScalar = new GvEquationTypeScalar(gvFieldScalar);    
#else
#ifdef NEW_SCALAR_FIELD
    using  GvEquationTypeScalar=  MagFieldEquation<Field_Type_Scalar>;   // New scalar
    #define gvFieldScalar gvField;
    auto    gvFieldScalar =  new UniformMagField( fieldValueVec );
    auto gvEquationScalar = new  GvEquationTypeScalar(gvFieldScalar);    // Same as vector, yes
    //                           ^was MagFieldEquation
    auto  myStepperScalar= CashKarp<GvEquationTypeScalar,Nposmom>(gvEquationScalar);
#else
    using GvEquationTypeScalar=  TMagFieldEquation<Field_Type_Scalar, Nposmom>;
    // If we plan to test against 'plain' scalar objects: field, equation, stepper, ... 
    auto fieldValueVec = fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field);
    auto gvFieldScalar = new TUniformMagField( fieldValueVec );
//                      new Field_Type_Scalar( fieldValueVec );
    auto gvEquationScalar = new GvEquationTypeScalar(gvFieldScalar);

    GUVIntegrationStepper *myStepperScalar; 
    myStepperScalar= StepperFactory::CreateStepper<GvEquationTypeScalar>(gvEquationScalar, stepper_no);
#endif
#endif

    int statisticsVerbosity = 1;
    auto refScalarDriver= new ScalarIntegrationDriver( hminimum,
                                                    myStepperScalar,
                                                    Nposmom,
                                                    statisticsVerbosity); 
    // refScalarDriver->InitializeCharge( particleCharge );
    //==========  Scalar Driver prepared =========================



    // -- Call internal Test code of Driver
    // int numTracks;
#if 0        
    bool ok= vectorDriver->TestInitializeLanes<Double_v>();
    if( !ok ) {
       std::cerr << "WARNING> vectorDriver<Double_v>->TestInitializeLanes() returned " << ok << std::endl;
    }
#endif    
    // 
    // double total_step = 0.;

    Bool_v goodAdvance(true);
    double epsTol = 1.0e-5;

    // Double_v  hStep1( 10.0 ); 
    // goodAdvance = testDriver->AccurateAdvance( yTrackIn, hStep1, epsTol, yTrackOut );

    constexpr int nTracks = 16;
    FieldTrack yInput[nTracks], yOutput[nTracks];
    // double posMom[] ={0., 0., 0., 0., 1., 1.};

    // std::fill_n(hstep, nTracks, 20);    
    double hstep[nTracks] = { 11.0, 20.0, 33.0, 44.44, 55.555, 66.6666, 77.77777, 88.888888, 99.9999999, 101.0, 111.0, 122.0, 133.0, 144.0, 155.0, 166.0 };
    // double hstep[] = { 11.0, 200.0, 330.0, 444.40, 555.55, 666.666, 777.7777, 888.88888, 999.999999, 100.0, 11.0, 12.0, 13.0, 14.0, 15.3, 16.9 };

    double charge[nTracks] = { 1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0, -10.0, 11.0, -12.0, 13.0, -14.0, 15., -16. };    
        // = {0, 0, 0, 1, -.3, .4, 20, 178., 920.}; 
    bool   succeeded[nTracks];

    bool firstTrial= true;
#ifdef MAINTESTING
    firstTrial= false;  // Just go to the loop of calls.
#endif    
    
    if( firstTrial )
    {
       double
          x_pos = 20.,                 //pos - position  : input unit = mm
          y_pos = 20.,
          z_pos = 20.;
       double x_mom = 0. ;               //mom - momentum  : input unit = GeV / c
       double y_mom = 1. ;
       double z_mom = 1. ;
       const double mmGVf = fieldUnits::millimeter;
       const double ppGVf = fieldUnits::GeV ; 

       double posMom[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                          x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};
    // double posMom[] = {0.0513401, 0.095223, 0.0916195, 0.635712, 0.717297, 0.141603 };

       for (int j = 0; j < nTracks; ++j)
       {
          yInput [j].LoadFromArray(posMom);
          yOutput[j].LoadFromArray(posMom); // Not strictly needed
       }
          
       vectorDriver
          ->AccurateAdvance<Double_v>( yInput, hstep, charge, epsTol, yOutput, nTracks, succeeded );
       // ==========================
       cout<<"-- Vector Driver done (advanced)." << endl;

       
       const ThreeVector_d  startPosition( posMom[0], posMom[1], posMom[2]);
       const ThreeVector_d  startMomentum( posMom[3], posMom[4], posMom[5]);
       ScalarFieldTrack yTrackIn ( startPosition, startMomentum, charge[0] );  // yStart
       ScalarFieldTrack yTrackOut( startPosition, startMomentum, charge[0] );  // yStart

       bool scalarResult =
             refScalarDriver->AccurateAdvance( yTrackIn, hstep[0], epsTol, yTrackOut );
       
       cout<<"-- Scalar Driver done (advanced)." << endl;
       cout<<" yOutput is   : "<< yOutput[0]<<" for yInput: "  <<yInput[0]<< endl;
       cout<<" yTrackOut is : "<< yTrackOut <<" for yTrackIn: "<<yTrackIn << endl;
       cout<<" Success of Vector: "<< succeeded[0] << endl;
       cout<<" Success of scalar: "<< scalarResult << endl;
       
       for (int i = 0; i < nTracks; ++i)
       {
          refScalarDriver->AccurateAdvance( yTrackIn, hstep[i], epsTol, yTrackOut );
          
          cout<<" yOutput["<<i<<"] is: "<< yOutput[i]<<" for yInput: "  <<yInput[i]<< endl;
          cout<<" yTrackOut is : "      << yTrackOut <<" for yTrackIn: "<<yTrackIn <<" for hstep: "<<hstep[i]<< endl;
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
      ScalarFieldTrack yTrackIn ( startPosition, startMomentum );  // yStart
      ScalarFieldTrack yTrackOut( startPosition, startMomentum );  // yStart

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
      vectorDriver->AccurateAdvance<Double_v>( yInput, hstep, epsTol, yOutput, nTracks, succeeded );
#endif 

#ifdef CALCULATETIME
      clock1 = clock() - clock1 - biasClock;
      float clock1InFloat = ((float)clock1)/CLOCKS_PER_SEC;
      cout<<"Vector time is: "<<clock1InFloat<<endl;
#endif 
/*      refScalarDriver->AccurateAdvance( yTrackIn, hstep[0], epsTol, yTrackOut );

      cout<<" yOutput[0] is: "<< yOutput[0]<<" for yInput: "  <<yInput[0]<< endl;
      cout<<" yTrackOut is: " << yTrackOut <<" for yTrackIn: "<<yTrackIn << endl;*/

/*      for (int i = 0; i < nTracks; ++i)
      {
        const ThreeVector_d  startPosition( posMomMatrix[i][0], posMomMatrix[i][1], posMomMatrix[i][2]);
        const ThreeVector_d  startMomentum( posMomMatrix[i][3], posMomMatrix[i][4], posMomMatrix[i][5]);
        ScalarFieldTrack yTrackIn ( startPosition, startMomentum ); 
        ScalarFieldTrack yTrackOut( startPosition, startMomentum ); 

        refScalarDriver->AccurateAdvance( yTrackIn, hstep[i], epsTol, yTrackOut );

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
        refScalarDriver->AccurateAdvance( yTrackIn, hstep[i], epsTol, yTrackOut );

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
      ScalarFieldTrack yTrackIn ( startPosition, startMomentum, charge );
      ScalarFieldTrack yTrackOut( startPosition, startMomentum, charge );

      vectorDriver->AccurateAdvance( yInput,
                                     charge,
                                     hstep,
                                     epsTol,
                                     yOutput,
                                     nTracks,
                                     succeeded );
      // refScalarDriver->AccurateAdvance( yTrackIn, hstep[11], epsTol, yTrackOut );

      for (int i = 0; i < nTracks; ++i)
      {
        refScalarDriver->AccurateAdvance( yTrackIn[i], hstep[i], epsTol, yTrackOutScalar[i] );

        // cout<<" yOutput["<<i<<"] is: "<< yOutput[i]<<" for hstep: "<<hstep[i]<< endl ;
        cout<<" yOutput["<<i<<"] is: "<< yOutput[i]<<" for hstep: "<<hstep[i]<< endl ;
        cout<<" yTrackOut is : "      << yTrackOut << endl;
      }

      // ScalarFieldTrack randNew = yTrackIn + yTrackOut;

      // cout<<yOutput[0]<<endl;
      // cout<<yOutput[1]<<endl;
      // cout<<yOutput[2]<<hstep[2]<<endl;

    }
#endif 

    cout<<" Scalar Stepper function calls are: "<< refScalarDriver->fStepperCalls <<" and OneGoodStep calls are "<<refScalarDriver->fNoTotalSteps << endl;
    cout<<" Vector Stepper function calls are: "<< vectorDriver->GetNumberOfStepperCalls() <<" and OneStep calls are "<<vectorDriver->GetNumberOfTotalSteps() << endl;


    //========================End testing IntegrationDriver=======================

    /*------ Clean up ------*/
    delete myStepper; 
    delete gvField;

    // deleting IntegrationDriver
    delete vectorDriver;
    delete refScalarDriver;      
    
    cout<<"\n\n#-------------End of output-----------------\n";
}

