//
//  Compare the output of different steppers
//     TClassicalRK4
//     GUTCashKarpRKF45
//  against each other.
//  NOTE: *Scalar* implementations only
//
//  Based on testStepperFixed.cc
//    was the work of Somnath Banerjee in GSoC 2015
//
#include <iomanip>


#include "Units.h"

// using fieldUnits::meter;
using fieldUnits::millimeter;   
using fieldUnits::second;  
using fieldUnits::eplus;  
using fieldUnits::tesla;
using fieldUnits::degree;


#include "TUniformMagField.h"

#include "TMagFieldEquation.h"
#include "FieldEquationFactory.h"

#include "GUVIntegrationStepper.h"
#include "StepperFactory.h"

#include "TClassicalRK4.h"
#include "GUTCashKarpRKF45.h"
#include "TSimpleRunge.h"
#include "GUExactHelixStepper.h"

using namespace std;


int main(int argc, char *args[])
{
    constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

    using  GvEquationType=  TMagFieldEquation<TUniformMagField, Nposmom>;
   
    /* -----------------------------SETTINGS-------------------------------- */
    
    /* Parameters of test
     - Modify values  */
    
    int no_of_steps = 250;         // No. of Steps for the stepper
    int stepper_no =  5;         // Choose stepper no., for refernce see above
    double step_len_mm = 200.;    // meant as millimeter;  //Step length 
    double z_field_in = DBL_MAX;
    
    //Checking for command line values :
    if(argc>1)
        stepper_no = atoi(args[1]);
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
    #ifdef BASELINESTEPPER
      bool printInpX= 0;  // print the input values for Ref 
    #endif 
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
       z_field = -1.0;  //  Tesla // *tesla ;

    // Field
    auto gvField= new TUniformMagField( fieldUnits::tesla * ThreeVector(x_field, y_field, z_field) );

    cout << "#  Initial  Field strength (GeantV) = "
         << x_field << " , " << y_field << " , " << z_field 
       // << (1.0/fieldUnits::tesla) * gvField->GetValue()->X() << ",  "
         << " Tesla " << endl;
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

    myStepper = new GUTCashKarpRKF45<GvEquationType,Nposmom>(gvEquation);
    // myStepper= StepperFactory::CreateStepper<GvEquationType>(gvEquation, stepper_no);
    // myStepper= StepperFactory::CreateStepper<decltype(gvEquation)*>(gvEquation, stepper_no);

    if( useClonedStepper ){
       auto baseStepper = myStepper;
       auto cloneStepper = myStepper->Clone();
       delete baseStepper;
       myStepper = cloneStepper;
    }

    myStepper->InitializeCharge( particleCharge );
    
    //Initialising coordinates
    const double mmGVf = fieldUnits::millimeter;
    const double ppGVf = fieldUnits::GeV ;  //   it is really  momentum * c_light
                                         //   Else it must be divided by fieldUnits::c_light;
    // const double ppGVf = fieldUnits::GeV / Constants::c_light;     // OLD

    // double yIn[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    double yIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                    x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};
    

    // double yInX[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
    //                 x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};    
    const double mmRef = mmGVf; // Unit for reference of lenght   - milli-meter
    const double ppRef = ppGVf; // Unit for reference of momentum - GeV / c^2
    
    // auto gvEquation2 = new GvEquationType(gvField);
                   // new TMagFieldEquation<TUniformMagField, Nposmom>(gvField);
    // gvEquation2->InitializeCharge( particleCharge ); // Let's make sure
    
    // Should be able to share the Equation -- eventually
    // For now, it checks that it was Done() -- and fails an assert

    //Creating the baseline stepper
  #ifdef BASELINESTEPPER
    auto exactStepperGV =
        new TClassicalRK4<GvEquationType,Nposmom>(gvEquation2);
    cout << "#  Reference stepper is: TClassicalRK4<GvEquationType,Nposmom>(gvEquation2);" << endl;

       // new TSimpleRunge<GvEquationType,Nposmom>(gvEquation2);    
       // new GUExactHelixStepper(gvEquation2);

    // Configure Stepper for current particle
    exactStepperGV->InitializeCharge( particleCharge ); // Passes to Equation, is cached by stepper
    // gvEquation2->InitializeCharge( particleCharge ); //  Different way - in case this works
    
    auto exactStepper = exactStepperGV;
  #endif 
    std::cout << "# step_len_mm = " << step_len_mm;
    std::cout << " mmRef= " << mmRef << "   ppRef= " << ppRef << std::endl;
    
    // double yInX[] = {x_pos * mmRef, y_pos * mmRef ,z_pos * mmRef,
                     // x_mom * ppRef ,y_mom * ppRef ,z_mom * ppRef};

    // double stepLengthRef = step_len_mm * mmRef;
    
    //Empty buckets for results
    double dydx   [8] = {0.,0.,0.,0.,0.,0.,0.,0.};  // 2 extra safety buffer
    // double dydxRef[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    double yout   [8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    // double youtX  [8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    double yerr   [8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    // double yerrX  [8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    
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
    no_of_steps = 5;
    for(int j=0; j<no_of_steps; j++)
    {
        cout<<setw(6)<<j ;           //Printing Step number

        myStepper->RightHandSideVIS(yIn, dydx);               //compute dydx - to supply the stepper
        #ifdef baseline
        exactStepper->RightHandSideVIS(yInX, dydxRef);   //compute the value of dydx for the exact stepper
        #endif

        if( j > 0 )  // Let's print the initial points!
        {
           myStepper->StepWithErrorEstimate(yIn,dydx,step_len,yout,yerr);   //Call the 'trial' stepper
/*           cout.setf (ios_base::fixed);
           cout.precision(4);
           cout<<"\n   Ananya yout: "<<endl;
           for (int tempIndex = 0; tempIndex < 6; ++tempIndex)
           {
             cout<<yout[tempIndex]*10<<" ";
           }
           cout<<"\n   Ananya dydx: "<<endl; 
           for (int tempIndex = 0; tempIndex < 6; ++tempIndex)
           {
             cout<<dydx[tempIndex]<<" ";
           }    */       
        
          #ifdef  BASELINESTEPPER
           exactStepperGV->StepWithErrorEstimate(yInX,dydxRef,stepLengthRef,youtX,yerrX); //call the reference stepper
          #endif
        }
        //-> Then print the data
        cout.setf (ios_base::fixed);
        cout.precision(4);
        for(int i=0; i<3;i++)
            if(columns[i]){
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(nwdf-2)<< yIn[i] / mmGVf;
              #ifdef BASELINESTEPPER
               if( printInpX ) cout << setw(nwdf-2)<< yInX[i] / mmRef;
               if( printRef ) cout<<setw(nwdf)<< youtX[i] / mmRef; // Reference Solution
               if( printDiff )                
                  cout<<setw(nwdf)<< yout[i] /  mmGVf - youtX[i] / mmRef ;
              #endif
               cout<<setw(nwdf)<< yout[i] / mmGVf ;
               
               if( printErr ) cout<<setw(nwdf)<< yerr[i] / mmGVf ;
               cout.precision(3);
               
            }

        cout.unsetf (ios_base::fixed);        
        cout.setf (ios_base::scientific);
        for(int i=3; i<6;i++)
            if(columns[i]){
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(nwdf-1)<< yIn[i] / ppGVf << " ";
              #ifdef BASELINESTEPPER
               if( printInpX ) cout << setw(nwdf-1)<< yInX[i] / ppRef << " ";
               if( printRef ) cout<<setw(nwdf)<< youtX[i] / ppRef; // Reference Solution
               if( printDiff ) 
                  cout<<setw(nwdf+2)<< ( yout[i] /  ppGVf )
                                     - ( youtX[i] / ppRef );
              #endif
               cout<<setw(nwdf) << yout[i] / ppGVf ;
               if( printErr ) cout<<setw(nwdf)<< yerr[i] / ppGVf ;

            }
        cout.unsetf (ios_base::scientific);
        
        for(int i=0; i<6;i++)   // Print auxiliary components
        {
           double unitGVf=1;  
           // double unitRef=1;
           // if( i < 3 )             // length / length

           if( i >= 3 ){
              unitGVf = ppGVf / mmGVf; //  dp / ds
            #ifdef BASELINESTEPPER
              unitRef = ppRef / mmRef; //  dp / ds
            #endif
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
               #ifdef BASELINESTEPPER
               if( printRef )
                 cout<<setw(nwdf)<< dydxRef[i] / unitRef; // Reference Solution
                 if( printDiff ) // Deriv ) 
                  cout<<setw(nwdf)<< ( dydx[i] / unitGVf )
                                   - ( dydxRef[i] / unitRef );
                #endif 
               // bool printDiffDeriv = true;

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
           
          #ifdef BASELINESTEPPER                   // atan(yout[1]/yout[0])/degree;
           if( printRef ) 
             cout<<setw(nwdf) << atan2(youtX[1],youtX[0])/degree;
           #endif
           //Copy yout into yIn
           for(int i=0; i<6;i++){
              yIn[i] = yout[i];
              #ifdef BASELINESTEPPER
              yInX[i] = youtX[i];
              #endif
           }
        }
        
        cout<<"\n";
    }
    
    /*-----------------END-STEPPING------------------*/

    /*------ Clean up ------*/
    myStepper->InformDone(); 
    
    #ifdef BASELINESTEPPER
    exactStepper->InformDone();
    #endif 

    delete myStepper;

    #ifdef BASELINESTEPPER
    delete exactStepper;
    #endif 
    // delete gvEquation;  // The stepper now takes ownership of the equation
    // delete gvEquation2;    
    delete gvField;
    
    cout<<"\n\n#-------------End of output-----------------\n";
    
}
