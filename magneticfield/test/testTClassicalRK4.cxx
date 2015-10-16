//
//
#include "GUVEquationOfMotion.h"
#include "GUVIntegrationStepper.h"

#include "TMagFieldEquation.h"
#include "TUniformMagField.h"
#include "TClassicalRK4.h"
// #include "TCashKarpRKF45.h"

bool  TestStepper(GUVIntegrationStepper*);

const unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

ThreeVector    FieldValue(0.0, 0.0, 1.0);

// GUVField*       pField;

GUVIntegrationStepper* CreateStepper(ThreeVector          FieldValue, 
                                     int                  stepperId, 
                                     GUVEquationOfMotion* equation)
{
  TUniformMagField*     ConstBfield= new TUniformMagField( FieldValue );

  // pField = ConstBfield; 
  // GUVEquationOfMotion*

  // typedef typename TMagFieldEquation<TUniformMagField, Nposmom>  EquationType;  // Fails 
  using EquationType = TMagFieldEquation<TUniformMagField, Nposmom>;
     
  // TMagFieldEquation<TUniformMagField, Nposmom> *
  // auto
  EquationType *magEquation= new EquationType(ConstBfield);
  GUVIntegrationStepper*  stepper;
  switch ( stepperId )
  {
    case 4 : 
       stepper = new TClassicalRK4<EquationType, Nposmom>(magEquation); // , Nposmom);
       break;
    // case 3 : 
    //    stepper = new TSimpleRunge<EquationType, Nposmom>;
    //    break;
    // case 5 : 
    //   stepper = new TCashKarpRK45<EquationType, Nposmom>;    
    //   break;       
    default: 
       stepper = new TClassicalRK4<EquationType, Nposmom>(magEquation); // , Nposmom);
       break;              
  }
  equation= magEquation;

  return stepper;
}

int gVerbose= 1;

bool TestStepper(GUVIntegrationStepper* stepper)
{
  const double perMillion = 1e-6;
  bool   hasError= false;  // Return value
  
  ThreeVector PositionVec( 1., 2.,  3.);  // initial
  ThreeVector MomentumVec( 0., 0.1, 1.);
  ThreeVector FieldVec( 0., 0., 1.);  // Magnetic field value (constant)

  double PositionTime[4] = { PositionVec.x(), PositionVec.y(), PositionVec.z(), 0.0};
  // double PositionTime[4] = { 1., 2., 3., 4.};
  PositionTime[0]= PositionVec.x();
  PositionTime[1]= PositionVec.y();
  PositionTime[2]= PositionVec.z();

  // double magField[3];

  double dydx[Nposmom];
  double PositionMomentum[Nposmom];
  double yout[Nposmom];
  double yerr[Nposmom];  

  double charge= -1;

  for( int i = 0; i< 3; ++i)
  { 
    PositionMomentum[i]=   PositionVec[i];
    PositionMomentum[3+i]= MomentumVec[i];
  }
  /***
  PositionMomentum[0]= PositionVec[0];
  PositionMomentum[1]= PositionVec[1];  
  PositionMomentum[2]= PositionVec[2];
  PositionMomentum[3]= MomentumVec[0];
  PositionMomentum[4]= MomentumVec[1];  
  PositionMomentum[5]= MomentumVec[2];
   ***/

  // double FieldArr[3]= { FieldVec.x(), FieldVec.y(), FieldVec.z() };
  // equation->EvaluateRhsGivenB( PositionTime, FieldArr, charge, dydx );

  double step_len = 100.0; //  * millimeter;

  stepper->GetEquationOfMotion()->InitializeCharge( charge );

  stepper->RightHandSide(PositionMomentum, /*charge,*/ dydx);
  stepper->StepWithErrorEstimate(PositionMomentum, dydx, step_len, yout, yerr); 

  //  Check the output HERE
  //



  
  // 
  stepper->GetEquationOfMotion()->InformDone(); 
     //InformDone();  // Finished for now - value of charge is 'revoked'  

  
  ThreeVector  ForceVec( dydx[3], dydx[4], dydx[5]);

  // Check result
  double MdotF= MomentumVec.Dot(ForceVec);
  double BdotF= FieldVec.Dot(ForceVec); 
  
  if( std::fabs(MdotF) > perMillion * MomentumVec.Mag() * ForceVec.Mag() )
  { 
     std::cerr << "ERROR: Force due to magnetic field is not perpendicular to momentum!!"  << std::endl;
     hasError= true;
  }
  else if ( gVerbose )
  {
     std::cout << " Success:  Good (near zero) dot product momentum . force " << std::endl;
  }
  if( std::fabs(BdotF) > perMillion * FieldVec.Mag() * ForceVec.Mag() )
  { 
    std::cerr << "ERROR: Force due to magnetic field is not perpendicular to B field!"
              << std::endl; 
    std::cerr << " Vectors:  BField   Force " << std::endl;
    for ( int i=0; i<3; i++ )
      std::cerr << "   [" << i << "] " << FieldVec[i] << " " << ForceVec[i] << std::endl;

    hasError= true;
  }
  else if ( gVerbose )
  {
    std::cout << " Success:  Good (near zero) dot product magnetic-field . force " << std::endl;
  }

  return hasError;
}

int
main( int argc, char** argv)
{
  GUVEquationOfMotion*  equation= 0;
  int    stepperType = 4;   // Default: Classical RK4  - Later Cash-Karp
  int    maxStepper =  8;
  double step_len = 100.0;

  if( argc > 0){
      int type= atoi(argv[0]);
      if( type > 0 && type <= maxStepper ) { 
         stepperType = type; 
      } else {
         std::cerr << "ERROR: Invalid value for Stepper type = " << type
                   << std::endl
                   << "       It must be between 0 and " << maxStepper 
                   << std::endl;
         exit(1);
      }
  }   
  if (argc > 2)
     step_len = (float) (atof(argv[2]) ); // * millimeter);

  GUVIntegrationStepper* stepper= CreateStepper( FieldValue, 
                                                 stepperType, // 
                                                 equation );  // Returns equation 

  TestStepper(stepper);

  return 0;
}
