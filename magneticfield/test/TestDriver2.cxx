//
//  Compare the output of different steppers
//
//  Created by Ananya <ananya.ananya@cern.ch>
//*****************************************************************
//  Based on testStepperFixed.cc
//    was the work of Somnath Banerjee in GSoC 2015
//
#include <iomanip>
#include <ctime>

#include <numeric>

#include "Units.h"

// using fieldUnits::meter;
using fieldUnits::millimeter;   
using fieldUnits::second;  
using fieldUnits::eplus;  
using fieldUnits::tesla;
using fieldUnits::degree;

#include <Vc/Vc>
#include "base/Vector3D.h"
#include "base/Vector.h"

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

#include <stdlib.h>

#include "ScalarCMSmagField.h"
#include "TemplateCMSmagField.h"

using namespace std;

typedef vecgeom::Vector3D<double> ThreeVector_d;

const double kRMax = 9000  * fieldUnits::millimeter; 
const double kZMax = 16000 * fieldUnits::millimeter;

double RandR(){
    double r = (float) rand()/(RAND_MAX) ;
    r = r*kRMax; //because r is in range (0,9000) mm                               
    return r;
}

double RandZ(){
    double z = (float) rand()/(RAND_MAX) ;
    z = z*kZMax; //range of z is between -16k and 16k                                                                         
    int sign = rand()%2; //to define the sign, since it can be both positive and negative                                     
    if (sign==0){
        z= -z;
    }
    return z;
}

void GenVecCartSubR(double &x, double &y){
    x = RandR();
    y = RandR();
    if((x*x + y*y)> kRMax*kRMax){
        GenVecCartSubR(x,y);
    }
}

void GenVecCart(ThreeVector_d &pos){
    double x=0,y=0;
    double z = RandZ();
    GenVecCartSubR(x, y);
    pos.x()=x;
    pos.y()=y;
    pos.z()=z;
}

void GenVecCart(vecgeom::Vector<ThreeVector_d> &posVec, const int &n){
    for (int i = 0; i < n; ++i)
    {       
        ThreeVector pos;
        GenVecCart(pos);
        posVec.push_back(pos);
    }
}

void MeanAndStDev(std::vector<double> timeVec, double &mean, double &stdev)
{
  double sum   = std::accumulate(timeVec.begin(), timeVec.end(), 0.0);
         mean  = sum/timeVec.size();
  double sqSum = std::inner_product(timeVec.begin(), timeVec.end(), timeVec.begin(), 0.0);
         stdev = std::sqrt(sqSum/timeVec.size() - mean*mean);
}


int main(int argc, char *args[])
{
  constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

  using Backend1 = vecgeom::kVc ;
  using Backend2 = vecgeom::kScalar;
  typedef vecgeom::Vector3D<double> ThreeVector_d;
  
// #define USECMSFIELD
#ifdef USECMSFIELD
  using Field_Type        = TemplateCMSmagField<Backend1>;
  using Field_Type_Scalar = TemplateCMSmagField<Backend2>;
  // using Field_Type_Scalar = ScalarCMSmagField;
#else
  using Field_Type        = TemplateTUniformMagField<Backend1>;
  using Field_Type_Scalar = TemplateTUniformMagField<Backend2>;
#endif 


  /* -----------------------------SETTINGS-------------------------------- */
  
  /* Parameters of test
   - Modify values  */
  
  int no_of_steps = 20;         // No. of Steps for the stepper
  // int stepper_no  =  5;         // Choose stepper no., for refernce see above
  double step_len_mm = 200.;    // meant as millimeter;  //Step length 
  double z_field_in = DBL_MAX;
  
  // Checking for command line values :
  // if(argc>1) stepper_no = atoi(args[1]);
  if(argc > 2)
     step_len_mm = (float)(stof(args[2]));   // *mm);
  if(argc > 3)
      no_of_steps = atoi(args[3]);
  if(argc > 4)
     z_field_in = (float) (stof(args[4]));     // tesla

  // double step_len = step_len_mm * fieldUnits::millimeter;
  
  // Set coordinates here
  /***  double
     x_pos = 0.,                 //pos - position  : input unit = mm
     y_pos = 0.,
     z_pos = 0.;
   ***/
  double   
     x_mom = 0.,                 //mom - momentum  : input unit = GeV / c
     y_mom = 1.,
     z_mom = 1.;
  double          
     x_field = 0.,               //Uniform Magnetic Field (x,y,z)
     y_field = 0.;               //  Tesla // *tesla ;
  double z_field;

  if( z_field_in < DBL_MAX )
     z_field = z_field_in;
  else
     z_field = -1.0;  //  Tesla // *tesla ;


  // Field
#ifdef USECMSFIELD
  auto gvField = new Field_Type ("../VecMagFieldRoutine/cmsmagfield2015.txt");
#else 
  auto gvField= new Field_Type( fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field) );
#endif

  cout << "#  Initial  Field strength (GeantV) = "
       << x_field << " , " << y_field << " , " << z_field 
     // << (1.0/fieldUnits::tesla) * gvField->GetValue()->X() << ",  "
       << " Tesla " << endl;
  cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;

  // Create an Equation :
  using GvEquationType    = TemplateTMagFieldEquation<Backend1, Field_Type, Nposmom>;
  auto gvEquation =
     TemplateFieldEquationFactory<Backend1>::CreateMagEquation< Field_Type >(gvField);


  /*-------------------------PREPARING STEPPER-----------------------------*/
  
  // Create a stepper :

  TemplateGUVIntegrationStepper<Backend1> *myStepper = new TemplateGUTCashKarpRKF45<Backend1,GvEquationType,Nposmom>(gvEquation);

//  const double mmGVf = fieldUnits::millimeter;
  const double ppGVf = fieldUnits::GeV ;  //   it is really  momentum * c_light
                                       //   Else it must be divided by fieldUnits::c_light;

  std::cout << "# step_len_mm = " << step_len_mm;
  
  
  /*-----------------------END PREPARING STEPPER---------------------------*/


  //=======================Test part for Integration driver====================
  double hminimum = 0.2;


  // Preparing scalar Integration Driver
  using  GvEquationTypeScalar=  TemplateTMagFieldEquation<Backend2, Field_Type_Scalar, Nposmom>;

#ifdef USECMSFIELD
  auto gvFieldScalar    = new Field_Type_Scalar("../VecMagFieldRoutine/cmsmagfield2015.txt");
#else
  auto gvFieldScalar    = new Field_Type_Scalar( fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field) );
#endif 

  auto gvEquationScalar = TemplateFieldEquationFactory<Backend2>::CreateMagEquation< Field_Type_Scalar >(gvFieldScalar);

  // GUVIntegrationStepper *myStepperScalar; 
  auto myStepperScalar= new TemplateGUTCashKarpRKF45<Backend2, GvEquationTypeScalar, Nposmom>(gvEquationScalar);

  // int statisticsVerbosity = 1;
  auto testScalarDriver= new TemplateGUIntegrationDriver<Backend2>
                                                ( hminimum,
                                                  myStepperScalar);

  auto testVectorDriver = new TemplateGUIntegrationDriver<Backend1>(hminimum, myStepper);

  bool chooseSteppingMethod;
  cout<<"Give 1 for OneStep and 0 for KeepStepping" << endl;
  cin >> chooseSteppingMethod;
  testVectorDriver->SetSteppingMethod(chooseSteppingMethod); 

  //  double total_step = 0.;

  typedef typename Backend1::bool_v Bool;
  Bool goodAdvance(true);
  double epsTol = 1.0e-5;

  // double charge1 = -1.0;
  // goodAdvance = testDriver->AccurateAdvance( yTrackIn, charge1, total_step, epsTol, yTrackOut );

  constexpr int nTracks = 16;

  FieldTrack yInput[nTracks];
  // FieldTrack yOutput[nTracks];

  // double posMom[] ={0., 0., 0., 0., 1., 1.};
  // double hstep[nTracks] = {0}; // = {0, 0, 0, 1, -.3, .4, 20, 178., 920.}; 

  bool   succeeded[nTracks];

#define TIMINGTESTING 
#define CALCULATETIME

 //=======================Proper timing calculation===========================
 // Storing random input data in array of arrays and 
 // vector of GUFieldTrack for vector and sequential version respectively. 

#ifdef TIMINGTESTING 
  int nRepititions = 1;

  constexpr int noOfVectorCalls = 32; // scalarcalls = nTracks*noOfVectorCalls

  no_of_steps = 1;

  // bool debugValue ; 
  // cout<< "Debug? " << endl;
  // cin >> debugValue;
  // testVectorDriver->SetPartDebug(debugValue);
  cout << "Give no_of_steps: "     << endl;
  cin >> no_of_steps;
  cout << "Give nRepititions: "    << endl;
  cin >> nRepititions;
  // cout << "Give noOfVectorCalls: " << endl;
  // cin >> noOfVectorCalls;
  // cout << "Using noOfVectorCalls: " << noOfVectorCalls << endl;
  cout << "Compiled parameter: noOfVectorCalls = " << noOfVectorCalls << endl;

  
  std::vector<double> speedUp, scalarTime, vectorTime;
  // std::vector<GUFieldTrack> vectorGUFieldTrack;
  long double outputVarForScalar = 0, outputVarForVector = 0;
  int indOutputVar = 2;
  
  int randVal = time(NULL);
  // srand(time(NULL));
  cout<<"Give seed for rng" << endl;
  cin >> randVal;
  srand(randVal);
  // srand(1458229725);
  cout<< "Random value initializer is : "<< randVal << endl;
  // srand(19);
  
  vecgeom::Vector<ThreeVector_d> posVec; 
  // GenVecCart( posVec, noOfVectorCalls * nTracks*no_of_steps);

  // int indPosVec = 0;

  for (int step = 0; step < no_of_steps; ++step)
  {
     // double X_Pos[nTracks], Y_Pos[nTracks], Z_Pos[nTracks];
    double X_Mom[nTracks], Y_Mom[nTracks], Z_Mom[nTracks];
    double posMomMatrix[nTracks][6];
    FieldTrack yInputMatrix[noOfVectorCalls][nTracks]; // [6];
    FieldTrack yOutput[nTracks]; // [6];    
    // std::vector<GUFieldTrack> vectorGUFieldTrack;
    std::vector<TemplateGUFieldTrack<Backend2> > vectorGUFieldTrack;

    int indPosVec = 0;
    GenVecCart( posVec, noOfVectorCalls * nTracks);

    for (int j = 0; j < noOfVectorCalls; ++j)
    {
      for (int i = 0; i < nTracks; ++i)
      {
/*        X_Pos[i] = (float) rand()/(RAND_MAX) ;
        Y_Pos[i] = (float) rand()/(RAND_MAX) ;
        Z_Pos[i] = (float) rand()/(RAND_MAX) ;*/
        X_Mom[i] = (float) rand()/(RAND_MAX) ;
        Y_Mom[i] = (float) rand()/(RAND_MAX) ;
        Z_Mom[i] = (float) rand()/(RAND_MAX) ;

/*        posMomMatrix[i][0] = X_Pos[i] * mmGVf;
        posMomMatrix[i][1] = Y_Pos[i] * mmGVf;
        posMomMatrix[i][2] = Z_Pos[i] * mmGVf;*/
        posMomMatrix[i][0] = posVec[indPosVec][0];
        posMomMatrix[i][1] = posVec[indPosVec][1];
        posMomMatrix[i][2] = posVec[indPosVec][2];       
        posMomMatrix[i][3] = X_Mom[i] * ppGVf;
        posMomMatrix[i][4] = Y_Mom[i] * ppGVf;
        posMomMatrix[i][5] = Z_Mom[i] * ppGVf;
        yInput [i].LoadFromArray(posMomMatrix[i]);
        // yOutput[i].LoadFromArray(posMomMatrix[i]);

        yInputMatrix[j][i] = yInput[i];

        const ThreeVector_d  startPosition( posMomMatrix[i][0], posMomMatrix[i][1], posMomMatrix[i][2]);
        const ThreeVector_d  startMomentum( posMomMatrix[i][3], posMomMatrix[i][4], posMomMatrix[i][5]);
        TemplateGUFieldTrack<Backend2> yTrackIn ( startPosition, startMomentum ); 
        vectorGUFieldTrack.push_back(yTrackIn); 

        // cout << "yInput["<<i<<"] is: " <<yInput[i] << " vs yTrackIn: "<< yTrackIn << endl;

        indPosVec++;
      }
    }

    // Random hstep between 0 and 200 cm (2m)
    // x, y, z values are multiplied with mmRef before being passed to function
    // the value of which is 0.1, so passing 200 directly would be in cm
    double hstepMatrix[noOfVectorCalls][nTracks];
    double charge[noOfVectorCalls][nTracks];    
    for (int j = 0; j < noOfVectorCalls; ++j)
    {
      for (int i = 0; i < nTracks; ++i)
      {
        hstepMatrix[j][i] = (float) rand()/(RAND_MAX) *200.;
        charge[j][i] = ( i < j ) ? -1.0 : +1.0 ;
      }
    }

    // for (int i = 0; i < nTracks; ++i)
    // {
    //   hstep[i] = (float) rand()/(RAND_MAX) *200.; 
    // }

    clock_t clock1 = clock();
    for (int repeat = 0; repeat < nRepititions; ++repeat)
    {
      for (int j = 0; j < noOfVectorCalls; ++j)
      {
        testVectorDriver->AccurateAdvance( yInputMatrix[j],
                                           charge[j],
                                           hstepMatrix[j],
                                           epsTol,
                                           yOutput,
                                           nTracks,
                                           succeeded );
        // testVectorDriver->AccurateAdvance( yInputMatrix[j], charge[j], hstep, epsTol, yOutput, nTracks, succeeded );
        for (int i = 0; i < nTracks; ++i)
        {
          // cout<<" yOutput["<<i<<"] is: "<< yOutput[i]<<" for yInput: "  <<yInput[i]<< " for hstep: " << hstepMatrix[j][i] << endl;
          outputVarForVector += yOutput[i].GetComponent(indOutputVar);// .PosMomVector[indOutputVar];          
        }      
      }
    }
    clock1 = clock() - clock1 ;
    float clock1InFloat = ((float)clock1)/CLOCKS_PER_SEC;
    vectorTime.push_back(clock1InFloat);
    cout<<"Vector time is: "<<clock1InFloat/(nRepititions*noOfVectorCalls*nTracks)*1e+6<<" ms"<<endl;


    const ThreeVector_d  startPos( 0, 0, 0);
    const ThreeVector_d  startMom( 0, 0, 0);
    TemplateGUFieldTrack<Backend2> yTrackOut( startPos, startMom ); 

    clock_t clock2 = clock();
    for (int repeat = 0; repeat < nRepititions; ++repeat)
    {
      int indScalar = 0;
      for (int j = 0; j < noOfVectorCalls; ++j)
      {
        for (int i = 0; i < nTracks; ++i)
        {
          // testScalarDriver->AccurateAdvance( vectorGUFieldTrack[i], hstep[i%nTracks], epsTol, yTrackOut );
          // testScalarDriver->AccurateAdvance( vectorGUFieldTrack[indScalar], hstep[i], epsTol, yTrackOut );
           testScalarDriver->AccurateAdvance( vectorGUFieldTrack[indScalar], charge[j][i], hstepMatrix[j][i], epsTol, yTrackOut );
          // cout<<" yTrackOut is : " << yTrackOut <<" for yTrackIn: "<<vectorGUFieldTrack[indScalar] <<" for hstep: "<< hstepMatrix[j][i] << endl;

          outputVarForScalar += yTrackOut.fPositionMomentum[indOutputVar];
          indScalar++;
        }      
      }
    }
    clock2 = clock() - clock2 ;
    float clock2InFloat = ((float)clock2)/CLOCKS_PER_SEC;
    scalarTime.push_back(clock2InFloat);
    cout<<"scalar time is: "<<clock2InFloat/(nRepititions*noOfVectorCalls*nTracks)*1e+6<<" ms"<< endl;
    speedUp.push_back(clock2InFloat/clock1InFloat);

  }

  // cout<<"indPosVec at end is: " << indPosVec <<" should be equal to: " << nTracks*noOfVectorCalls*no_of_steps << endl;
  cout<<"outputVarForScalar: "<< outputVarForScalar<< endl;
  cout<<"outputVarForVector: "<< outputVarForVector<< endl;
  int sizeOfRatioVector = speedUp.size(); // no_steps;
  cout<<"Size of speedUp is: "<<speedUp.size()<<endl;
  cout<<"Time ratios are: "<<endl;
  for (int i = 0; i < sizeOfRatioVector; ++i)
  {
    cout<< sizeOfRatioVector - i<<": "<<speedUp.back()<< " " <<endl;
    speedUp.pop_back();
  }
  cout<<endl;

  double vecMean, scaMean, vecStDev, scaStDev;

  MeanAndStDev(vectorTime, vecMean, vecStDev);
  MeanAndStDev(scalarTime, scaMean, scaStDev);

  cout<<" Ratio of mean timings is: "<< scaMean/vecMean << endl;

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
