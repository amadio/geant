//
//  Compare the output of different steppers
// 
//  Based on testStepperFixed.cc
//    was the work of Somnath Banerjee in GSoC 2015
//
#include <cfloat>
#include <iomanip>

// #include "G4UniformMagField.hh"
// #include "G4SystemOfUnits.hh"
#include "Units.h"

// using fieldUnits::meter;
using fieldUnits::millimeter;   
using fieldUnits::second;  
using fieldUnits::eplus;  
using fieldUnits::tesla;
using fieldUnits::degree;

// #include "G4ios.hh"

// #define  COMPARE_TO_G4  1

#include "TUniformMagField.h"

#include "TMagFieldEquation.h"
#include "FieldEquationFactory.h"

#include "GUVIntegrationStepper.h"
#include "StepperFactory.h"

#include "TClassicalRK4.h"
#include "GUTCashKarpRKF45.h"
#include "TSimpleRunge.h"
#include "GUExactHelixStepper.h"

// #include "BogackiShampine23.hh"
// #include "DormandPrince745.hh"
// #include "BogackiShampine45.h"

// #include "G4SimpleHeum.hh"

#ifdef COMPARE_TO_G4
//#include <system.h>
//#include "G4Types.h"
#include "G4UniformMagField.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ClassicalRK4.hh"
#endif

using namespace std;
// using namespace CLHEP;

//Version 2.0 - Includes direct comparison with ExactHelix

/* Stepper No.
 0. GUExactHelix
 4. TClassicalRK4
 5. TCashKarpRKF45

Potential expansion:
 2. G4SimpleHeum
 3. BogackiShampine23
 6. BogackiShampine45
 7. DormandPrince745
 */


int main(int argc, char *args[])
{
    constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

    using  GvEquationType=  TMagFieldEquation<TUniformMagField, Nposmom>;
    void Usage();

    using ThreeVectorF = vecgeom::Vector3D<float>;
    using ThreeVectorD = vecgeom::Vector3D<double>;
    
    /* -----------------------------SETTINGS-------------------------------- */
    
    /* Parameters of test
     - Modify values  */
    
    int no_of_steps = 250;         // No. of Steps for the stepper
    int stepper_no =  5;         // Choose stepper no., for refernce see above
    double step_len_mm = 200.;    // meant as millimeter;  //Step length 
    double z_field_in = DBL_MAX;
    
    //Checking for command line values :
    if(argc==1)
       Usage();

    if(argc>1) {
       if(   (strcmp(args[1],"-h")==0)     || (strcmp(args[1],"-u")==0)
          || (strcmp(args[1],"--help")==0) || (strcmp(args[1],"--usage")==0) ) {
         Usage();
         exit(0);
       } else { 
         stepper_no = atoi(args[1]);
       }
    }    
    if(argc > 2)
       step_len_mm = (float)(stof(args[2]));   // *mm);
    if(argc > 3)
        no_of_steps = atoi(args[3]);
    if(argc > 4)
       z_field_in = (float) (stof(args[4]));     // tesla
    double step_len = step_len_mm * fieldUnits::millimeter;
    
    //Set Charge etc.
    double particleCharge = +1.0;      // in e+ units
    
    //Choice of output coordinates
    int
    columns[] =
    {
       1 , 1 , 1 ,  // position  x, y, z 
       1 , 1 , 1 ,  // momentum  x, y, z
       0 , 0 , 0 ,  // dydx pos  x, y, z
       1 , 1 , 0    // dydx mom  x, y, z
    }; //Variables in yOut[] & dydx[] we want to display - 0 for No, 1 for yes

    bool printErr= 0,   // Print the predicted Error
         printRef= 0,   // Print the reference Solution 
         printDiff= 0;  // Print the diffrence 
    bool printSep = 0;  // separator  '|'
    bool printInp = 0;  // print the input values
    bool printInpX= 0;  // print the input values for Ref 
    
    const unsigned int nwdf= 12;  // Width for float/double
    
    //Set coordinates here
    double
       x_pos = 0.,                 //pos - position  : input unit = mm
       y_pos = 0.,
       z_pos = 0.;
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
       z_field = -1.0;  //  Tesla

    cout << "#  Input  Field strength (GeantV) = " << x_field << " , " << y_field << " , " << z_field 
       // << (1.0/fieldUnits::tesla) * gvField->GetValue()->X() << ",  "
         << " Tesla " << endl;

    // Field
    auto gvField= new TUniformMagField( fieldUnits::tesla * ThreeVectorF(x_field, y_field, z_field) );

    // double position[3] = { 0.0, 0.0, 0.0 }; 
    // double fieldArr[3] = { 0.0, 0.0, 0.0 };
    ThreeVectorD positionVec( 0.0, 0.0, 0.0 );
    ThreeVectorD fieldVec( 0.0, 0.0, 0.0 );
    
    gvField->GetFieldValue( positionVec, fieldVec );
    cout << "#DEBUG> Field value from TUniformMagField = " << fieldVec[0] / fieldUnits::kilogauss
         << " ,  " << fieldVec[1] / fieldUnits::kilogauss
         << " , "  << fieldVec[2] / fieldUnits::kilogauss << " KGauss " << endl; 

    cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;
    //Create an Equation :
    auto gvEquation =
       FieldEquationFactory::CreateMagEquation<TUniformMagField>(gvField);

       // new GvEquationType(gvField);
       // new TMagFieldEquation<TUniformMagField, Nposmom>(gvField);

    // gvEquation->InitializeCharge( particleCharge );  // Send it via Stepper instead    

    /*-------------------------PREPARING STEPPER-----------------------------*/
    
    //Create a stepper :
    GUVIntegrationStepper *myStepper; // , *exactStepper;
    // G4MagIntegrationStepper *g4refStepper;    
    const int cloneBump= 10;
    bool useClonedStepper= (stepper_no > cloneBump);
    if(  useClonedStepper )
       stepper_no -= cloneBump;

    myStepper= StepperFactory::CreateStepper<GvEquationType>(gvEquation, stepper_no);
    // myStepper= StepperFactory::CreateStepper<decltype(gvEquation)*>(gvEquation, stepper_no);

    if( useClonedStepper ){
       auto baseStepper = myStepper;
       auto cloneStepper = myStepper->Clone();
       delete baseStepper;
       myStepper = cloneStepper;
    }

    //myStepper->InitializeCharge( particleCharge );
    
    //Initialising coordinates
    const double mmGVf = fieldUnits::millimeter;
    const double ppGVf = fieldUnits::GeV ;  //   it is really  momentum * c_light
                                         //   Else it must be divided by fieldUnits::c_light;
    // const double ppGVf = fieldUnits::GeV / Constants::c_light;     // OLD

    // double yIn[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    double yIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                    x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};
    
#if COMPARE_TO_G4
    double mass = 1;
    
    const double mmG4 = CLHEP::millimeter;
    const double ppG4 = CLHEP::GeV ;  //  In G4 too 'p' means p*c -- so no division  / CLHEP::c_light;

    const double mmRef = mmG4;  // Unit for reference of lenght   - milli-meter
    const double ppRef = ppG4;  // Unit for reference of momentum - GeV / c^2
    
    // G4double yInX[] = {x_pos * mmG4, y_pos * mmG4 ,z_pos * mmG4,
    //                   x_mom * ppG4 ,y_mom * ppG4 ,z_mom * ppG4};

    G4UniformMagField myField( CLHEP::tesla * G4ThreeVector(x_field, y_field, z_field));
    cout << " Field strength = " << (1.0/CLHEP::tesla) * myField << " Tesla " 
         << endl;
    
    G4Mag_UsualEqRhs *g4Equation = new G4Mag_UsualEqRhs(&myField);
    
    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0);

    g4Equation->SetChargeMomentumMass( chargeState,
                                       G4ThreeVector(x_mom, y_mom, z_mom).mag(), //momentum magnitude
                                       mass);  // unused
//  auto g4exactStepper = new G4ExactHelixStepper(g4Equation);
    auto g4exactStepper = new G4ClassicalRK4(g4Equation);
    
    auto exactStepper = g4ExactStepperGV;
#else
    // double yInX[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
    //                 x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};    
    const double mmRef = mmGVf; // Unit for reference of lenght   - milli-meter
    const double ppRef = ppGVf; // Unit for reference of momentum - GeV / c^2
    
    auto gvEquation2 = new GvEquationType(gvField);
                   // new TMagFieldEquation<TUniformMagField, Nposmom>(gvField);
    // gvEquation2->InitializeCharge( particleCharge ); // Let's make sure
    
    // Should be able to share the Equation -- eventually
    // For now, it checks that it was Done() -- and fails an assert

    //Creating the baseline stepper
    auto exactStepperGV =
        new TClassicalRK4<GvEquationType,Nposmom>(gvEquation2);
    cout << "#  Reference stepper is: TClassicalRK4<GvEquationType,Nposmom>(gvEquation2);" << endl;
       // new TSimpleRunge<GvEquationType,Nposmom>(gvEquation2);    
       // new GUExactHelixStepper(gvEquation2);

    // Configure Stepper for current particle
    //exactStepperGV->InitializeCharge( particleCharge ); // Passes to Equation, is cached by stepper
    // gvEquation2->InitializeCharge( particleCharge ); //  Different way - in case this works
    
    auto exactStepper = exactStepperGV;
#endif
    std::cout << "# step_len_mm = " << step_len_mm;
    std::cout << " mmRef= " << mmRef << "   ppRef= " << ppRef << std::endl;
    
    double yInX[] = {x_pos * mmRef, y_pos * mmRef ,z_pos * mmRef,
                     x_mom * ppRef ,y_mom * ppRef ,z_mom * ppRef};

    double stepLengthRef = step_len_mm * mmRef;
    
    //Empty buckets for results
    double dydx[8] = {0.,0.,0.,0.,0.,0.,0.,0.},  // 2 extra safety buffer
        dydxRef[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
           yout[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          youtX[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
           yerr[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          yerrX[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    
    /*-----------------------END PREPARING STEPPER---------------------------*/
    


    /*---------------------------------------------------*/
    //        -> First Print the (commented) title header
    cout<<"\n#";
    cout<<setw(5)<<"Step";
    for (int i=0; i<6;i++)
        if (columns[i])
        {
           if( printSep ) cout << " | " ;  // Separator
           if( printInp )
              cout << setw(nwdf-2)<< "yIn[" << i << "]";
           if( printInp )
              cout << setw(nwdf-2)<< "yInX[" << i << "]";           
           cout << setw(nwdf-2)<< "yOut[" << i << "]";
           if( printRef )
              cout << setw(nwdf-2) << "yOutX[" << i << "]";
           if( printErr )
              cout << setw(nwdf-1) << "yErr[" << i << "]";
           if( printDiff )
              cout << setw(nwdf) << " yOut-yOutX[" << i << "]";
        }
    for (int i=0; i<6;i++)
        if (columns[i+6])
        {
           if( printSep ) cout << " | " ;  // Separator
           if( printInp )                               
              cout << setw(nwdf-2)<< "yIn[" << i << "]";           
           cout << setw(nwdf-2)<< "dydx[" << i << "]";
           if( printRef )                    
              cout << setw(nwdf-2) << "dydxRef[" << i << "]";
           // if( printErr ) cout << setw(nwdf-2) << "yErr[" << i << "]";
           if( printDiff ) //  printDiffDeriv )
              cout << setw(nwdf-2) << "dydx-Ref[" << i << "]";
        }    
    cout<<setw(nwdf)<<"tan-1(y/x)";
    cout<<"\n";
    
    // Units
    const char *nameUnitLength= "mm";
    const char *nameUnitMomentum= "GeV/c";
    cout<<setw(6)<<"#Numbr";
    for (int i=0; i<6;i++)
        if (columns[i])
        {
           if( printSep ) cout << " | " ;  // Separator                      
           const char* nameUnit = ( i<3 ) ? nameUnitLength : nameUnitMomentum ; 
           cout << setw(nwdf)<< nameUnit;
           if( printRef )  cout << setw(nwdf) << nameUnit;  // Reference output
           if( printErr )  cout << setw(nwdf) << nameUnit;  // Estim. error
           if( printDiff ) cout << setw(15) << nameUnit;    // Diff  new - ref
        }    
    cout<<"\n";
    
    //-> Then print the data

    
    // cout.setf (ios_base::scientific);
    // cout.setf (ios_base::scientific);    
    cout.precision(3);
    
    /*----------------NOW STEPPING-----------------*/
    
    for(int j=0; j<no_of_steps; j++)
    {
        cout<<setw(6)<<j ;           //Printing Step number

        myStepper->RightHandSideVIS(yIn, particleCharge, dydx);               //compute dydx - to supply the stepper
        exactStepper->RightHandSideVIS(yInX, particleCharge, dydxRef);   //compute the value of dydx for the exact stepper

        if( j > 0 )  // Let's print the initial points!
        {
           myStepper->StepWithErrorEstimate(yIn,particleCharge,dydx,step_len,yout,yerr);   //Call the 'trial' stepper
        
#ifdef COMPARE_TO_G4        
           g4ExactStepper->Stepper(yInX,dydxRef,stepLengthRef,youtX,yerrX); //call the reference stepper
#else
           exactStepperGV->StepWithErrorEstimate(yInX,particleCharge,dydxRef,stepLengthRef,youtX,yerrX); //call the reference stepper
#endif
        }
        //-> Then print the data
        cout.setf (ios_base::fixed);
        cout.precision(4);
        for(int i=0; i<3;i++)
            if(columns[i]){
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(nwdf-2)<< yIn[i] / mmGVf;
               if( printInpX ) cout << setw(nwdf-2)<< yInX[i] / mmRef;
               cout<<setw(nwdf)<< yout[i] / mmGVf ;
               if( printRef ) cout<<setw(nwdf)<< youtX[i] / mmRef; // Reference Solution
               if( printErr ) cout<<setw(nwdf)<< yerr[i] / mmGVf ;
               cout.precision(3);
               if( printDiff )                
                  cout<<setw(nwdf)<< yout[i] /  mmGVf - youtX[i] / mmRef ;
            }

        cout.unsetf (ios_base::fixed);        
        cout.setf (ios_base::scientific);
        for(int i=3; i<6;i++)
            if(columns[i]){
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(nwdf-1)<< yIn[i] / ppGVf << " ";
               if( printInpX ) cout << setw(nwdf-1)<< yInX[i] / ppRef << " ";
               cout<<setw(nwdf) << yout[i] / ppGVf ;
               if( printRef ) cout<<setw(nwdf)<< youtX[i] / ppRef; // Reference Solution
               if( printErr ) cout<<setw(nwdf)<< yerr[i] / ppGVf ;
               if( printDiff ) 
                  cout<<setw(nwdf+2)<< ( yout[i] /  ppGVf )
                                     - ( youtX[i] / ppRef );
            }
        cout.unsetf (ios_base::scientific);
        
        for(int i=0; i<6;i++)   // Print auxiliary components
        {
           double unitGVf=1;  
           double unitRef=1;
           // if( i < 3 )             // length / length

           if( i >= 3 ){
              unitGVf = ppGVf / mmGVf; //  dp / ds
              unitRef = ppRef / mmRef; //  dp / ds
           }
           
           if( i == 0 ){            // length / length                      
              // cout.unsetf (ios_base::scientific);
              cout.setf (ios_base::fixed);
           }else if( i == 3 ){
              cout.unsetf (ios_base::fixed);              
              cout.setf (ios_base::scientific);
           }
           
           if(columns[i+6])
           {
               // cout << " dy/dx["<<i<<"] = ";
               if( printSep ) cout << " | " ;  // Separator
               cout<<setw(nwdf)<< dydx[i] / unitGVf ;
               if( printRef )                
                 cout<<setw(nwdf)<< dydxRef[i] / unitRef; // Reference Solution
               // bool printDiffDeriv = true;
               if( printDiff ) // Deriv ) 
                  cout<<setw(nwdf)<< ( dydx[i] / unitGVf )
                                   - ( dydxRef[i] / unitRef );
           }
           // if( i == 2 )     { cout.unsetf(ios_base::fixed);      }
           // else if ( i==5 ) { cout.unsetf(ios_base::scientific); }
        }
        cout.unsetf(ios_base::scientific);
        if( j > 0 )  // Step 0 did not move -- printed the starting values
        {
           // cout.unsetf(ios_base::scientific);
           cout.setf (ios_base::fixed);                         
           cout.precision(2);
           cout<<setw(nwdf) << atan2(yout[1],yout[0])/degree;
                             // atan(yout[1]/yout[0])/degree;
           if( printRef ) 
             cout<<setw(nwdf) << atan2(youtX[1],youtX[0])/degree;
           
           //Copy yout into yIn
           for(int i=0; i<6;i++){
              yIn[i] = yout[i];
              yInX[i] = youtX[i];
           }
        }
        
        cout<<"\n";
    }
    
    /*-----------------END-STEPPING------------------*/

    /*------ Clean up ------*/
    myStepper->InformDone(); 
    
#ifndef COMPARE_TO_G4
    exactStepper->InformDone();
#endif
    delete myStepper;
    delete exactStepper;
    // delete gvEquation;  // The stepper now takes ownership of the equation
    // delete gvEquation2;    
    delete gvField;
    
    cout<<"\n\n#-------------End of output-----------------\n";
    
}

void Usage()
{
  cout << endl;
  cout <<   " This test cross-checks the output of a Stepper class against a referrence stepper (ClassicalRK4) " 
       <<   "   using a series of steps." << endl << endl;
  cout <<   " Usage of this test: " << endl;
  cout <<   "   arg[1]:  stepper type - values: 4 - ClassicalRK4,  5 - Cash-Karp  1 - Simple Runge (order 2)" << endl
       <<   "   arg[2]:  step_len_mm  - step length size (in mm) " << endl
       <<   "   arg[3]:  no_of_steps  - number of steps. " << endl
       <<   "   arg[4]:  z_field      - Magnetic field value (Tesla)" << endl;
}
