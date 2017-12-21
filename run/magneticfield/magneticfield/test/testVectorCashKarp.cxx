//
//  Compare the output of different steppers
// 
//  Based on testStepperFixed.cc
//    which was part of the work of Somnath Banerjee in GSoC 2015
//
#include <iomanip>

#include "base/Vector3D.h"
#include <Geant/VectorTypes.h>

using namespace vecCore::math;

#include "UniformMagField.h"          // New type (universal) class

// #include "GUVVectorEquationOfMotion.h"

#include "MagFieldEquation.h"

#include "GUVVectorIntegrationStepper.h"
#include "VectorCashKarpRKF45.h"

// #include "TemplateTSimpleRunge.h"
// #include "TemplateGUExactHelixStepper.h"

#include "Units.h"

// using fieldUnits::meter;
using fieldUnits::millimeter;   
using fieldUnits::second;  
using fieldUnits::eplus;  
using fieldUnits::tesla;
using fieldUnits::degree;

// using namespace std;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;

int main(int argc, char *args[])
{
    constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

    // using Backend = vecgeom::kVc ;
    // typedef typename Backend::precision_v Double_v;

    using Double_v  = Geant::Double_v;
    // template <typename Tpod>
    // using Vector3D  = vecgeom::Vector3D<Tpod>;
    using ThreeVector_f   = vecgeom::Vector3D<float>;
    using ThreeVector_d   = vecgeom::Vector3D<double>;
    // using ThreeVectorSimd = vecgeom::Vector3D<Double_v>;

    using  FieldType      =  UniformMagField;
    using  GvEquationType =  MagFieldEquation<FieldType>;
   
    /* -----------------------------SETTINGS-------------------------------- */
    
    // Parameters of test - values can be modified
    int no_of_steps = 250;        // No. of Steps for the stepper
    int stepper_no =  5;          // Choose stepper no., for refernce see above
    double step_len_mm = 200.;    // meant as millimeter;  //Step length 
    double z_field_in = DBL_MAX;

    bool debug= true;
    if( debug ) cout << "Running with debug (progress) messages." << endl;
    
    //Checking for command line values :
    if(argc>1)
        stepper_no = atoi(args[1]);
    if(argc > 2)
       step_len_mm = (float)(atof(args[2]));   // *mm);   //  stof( std::string, size_t ... )
    if(argc > 3)
        no_of_steps = atoi(args[3]);
    if(argc > 4)
       z_field_in = (float) (atof(args[4]));     // tesla
    double stepLengthValue = step_len_mm * fieldUnits::millimeter;
    
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
    // bool printInpX= 0;  // print the input values for Ref 
    
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
    // vecgeom::MaskedAssign(z_field_in < DBL_MAX, z_field_in, &z_field);

    if( z_field_in < DBL_MAX )
       z_field = z_field_in;
    else
       z_field = -1.0;  //  Tesla // *tesla ;

    if( debug ) cout<<"----Just before making UniformMagField object"<<endl;

    // Field
    auto gvUniformField= new
          UniformMagField( fieldUnits::tesla * ThreeVector_f(x_field, y_field, z_field) ); // New classes
     //    VectorUniformMagField<Backend>( fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field) );
     // TemplateScalarUniformMagField<Backend>( fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field) );
    
    if( debug ) cout<<"----UniformMagField Object constructed"<<endl;
    cout << "#  Initial  Field strength (GeantV) = "
         << x_field << " , " << y_field << " , " << z_field
         << " Tesla " << endl;

    ThreeVector_d origin(0.0, 0.0, 0.0), fieldValue;
    gvUniformField->GetFieldValue( origin, fieldValue );
    cout << "#    Values in object created       = "
         << (1.0/fieldUnits::tesla) * fieldValue.x() << ",  "
         << (1.0/fieldUnits::tesla) * fieldValue.y() << ",  "
         << (1.0/fieldUnits::tesla) * fieldValue.z() << endl;
    cout << endl;    
    cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;

    //Create Equation :
    if( debug ) cout << "Create Equation" << endl;

    // 1. Original way of creating an equation
    auto   magEquation = new GvEquationType(gvUniformField);
    if( debug ) cout<<"----Equation instantiated. "<<endl;    

    //Create a stepper :
    if( debug ) cout<<"---- Preparing to create (Vector) CashKarpRKF45 Stepper "<<endl;

    VectorCashKarpRKF45< /*Backend,*/ GvEquationType,Nposmom> myStepper2(magEquation);

    auto myStepper = &myStepper2;
    // myStepper = new VectorCashKarpRKF45<Backend,GvEquationType,Nposmom>(gvEquation);
                                                                                  
    if( debug )
    {
       cout << "---- constructed VectorCashKarpRKF45" << endl;
       cout << " Stepper  information: ptr= " << myStepper << endl;
       // cout << "       full " << *myStepper << endl;
    }
                                                              
    // Phase 1 - get it to work without cloning

    //Initialising coordinates
    const double mmGVf = fieldUnits::millimeter;
    const double ppGVf = fieldUnits::GeV ;  //   it is really  momentum * c_light
                                         //   Else it must be divided by fieldUnits::c_light;

    // Double_v yIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
    //                  x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};
    Double_v yInVec[]= { Double_v(x_pos * mmGVf), Double_v(y_pos * mmGVf) ,
                      Double_v(z_pos * mmGVf),
                      Double_v(x_mom * ppGVf), Double_v(y_mom * ppGVf) ,
                      Double_v(z_mom * ppGVf) };
    if( debug ) cout << "Initialized yIn: values [0]= " << yInVec[0] << endl;

    // double yInX[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
    //                 x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};    
    const double mmRef = mmGVf; // Unit for reference of lenght   - milli-meter
    const double ppRef = ppGVf; // Unit for reference of momentum - GeV / c^2
    
    // auto gvEquation2 = new GvEquationType(gvUniformField);
                   // new TMagFieldEquation<ScalarUniformMagField, Nposmom>(gvUniformField);
    // gvEquation2->InitializeCharge( particleCharge ); // Let's make sure
    
    // Should be able to share the Equation -- eventually
    // For now, it checks that it was Done() -- and fails an assert

    //Empty buckets for results
    Double_v dydxVec[8] = { Double_v(0.),0.,0.,0.,0.,0.,0.,0.},  // 2 extra safety buffer
             yOutVec[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
             yErrVec[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    
    //Creating the baseline stepper
  #ifdef BASELINESTEPPER
    auto exactStepperGV =
        new TClassicalRK4<GvEquationType,Nposmom>(gvEquation2);
    cout << "#  Reference stepper is: TClassicalRK4<GvEquationType,Nposmom>(gvEquation2);" << endl;

       // new TSimpleRunge<GvEquationType,Nposmom>(gvEquation2);    
       // new GUExactHelixStepper(gvEquation2);

    auto exactStepper = exactStepperGV;

    Double_v dydxVecRef[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
             yOutVecX[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
             yErrVecX[8] = {0.,0.,0.,0.,0.,0.,0.,0.};    
  #endif 
    cout << "# step_len_mm = " << step_len_mm;
    cout << " mmRef= " << mmRef << "   ppRef= " << ppRef << endl;
    
    // Double_v yInX[] = {x_pos * mmRef, y_pos * mmRef ,z_pos * mmRef,
                     // x_mom * ppRef ,y_mom * ppRef ,z_mom * ppRef};

    // double stepLengthRef = step_len_mm * mmRef;
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

    const double particleCharge= -1.0; 
    cout << " particle charge = " << particleCharge << endl;        
    
    /*----------------NOW STEPPING-----------------*/
    no_of_steps = 25;
    for(int j=0; j<no_of_steps; j++)
    {
        cout<<setw(6)<<j ;           //Printing Step number
        Double_v chargeVec( particleCharge );
        Double_v step_len( stepLengthValue ); 
        myStepper->RightHandSideInl(yInVec, chargeVec, dydxVec);  //compute dydx - to supply the stepper
#ifdef  BASELINESTEPPER
        double yinX[Nposmom], youtX[Nposmom], yerrX[Nposmom];
        exactStepper->RightHandSideVIS(yInVecX, dydxVecRef);   //compute the value of dydx for the exact stepper
#endif

        int lane=0; // Lane for printing
        double yin[Nposmom] = { -999., -99.9, -9.99, 10.0, 1.0, 0.1 };
        double yout[Nposmom], dydx[Nposmom], yerr[Nposmom];                  
        if( j > 0 )  // Do nothing for j=0, so we can print the initial points!
        {
           myStepper->StepWithErrorEstimate( yInVec, dydxVec, chargeVec, step_len, yOutVec, yErrVec );   //Call the 'trial' stepper
           //         *********************
#ifdef  BASELINESTEPPER
           exactStepperGV->StepWithErrorEstimate(yInVecX,dydxVecRef,chargeVec,stepLengthRef,yOutVecX,yErrVecX); //call the reference stepper
           //              *********************
#endif
        }

        // Select one lane for printing.
        for(unsigned int i=0; i<Nposmom;i++) {
           yin[i]   = vecCore::Get( yInVec[i], lane );           
           yout[i]  = ( j > 0 ) ? vecCore::Get( yOutVec[i], lane ) :
                                  vecCore::Get( yInVec[i], lane ); //  yOut is not yet set for j=0
           yerr[i]  = vecCore::Get( yErrVec[i], lane );
           dydx[i]  = vecCore::Get( dydxVec[i], lane );
#ifdef BASELINESTEPPER
           yinX[i]  = vecCore::Get( yInVecX[i], lane  );
           youtX[i] = ( j> 0 ) ? vecCore::Get( yOutVecX[i], lane ) :
                                 vecCore::Get( yInVecX[i], lane ); //  yOut is not yet set for j=0           
           yerrX[i] = vecCore::Get( yErrVecX[i], lane );
           dydxX[i] = vecCore::Get( dydxVecRef[i], lane );           
#endif           
        }

        //-> Then print the data
        cout.setf (std::ios_base::fixed);
        cout.precision(4);
        for(int i=0; i<3;i++)
            if(columns[i]){
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(nwdf-2)<< yin[i] / mmGVf;
              #ifdef BASELINESTEPPER
               if( printInpX ) cout << setw(nwdf-2)<< yinX[i] / mmRef;
               if( printRef ) cout<<setw(nwdf)<< youtX[i] / mmRef; // Reference Solution
               if( printDiff )                
                  cout<<setw(nwdf)<< yout[i] /  mmGVf - youtX[i] / mmRef ;
              #endif
               cout<<setw(nwdf)<< yout[i] / mmGVf ;
               
               if( printErr ) cout<<setw(nwdf)<< yerr[i] / mmGVf ;
               cout.precision(3);
               
            }

        cout.unsetf (std::ios_base::fixed);        
        cout.setf (std::ios_base::scientific);
        for(int i=3; i<6;i++)
            if(columns[i]){
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(nwdf-1)<< yin[i] / ppGVf << " ";
              #ifdef BASELINESTEPPER
               if( printInpX ) cout << setw(nwdf-1)<< yinX[i] / ppRef << " ";
               if( printRef ) cout<<setw(nwdf)<< youtX[i] / ppRef; // Reference Solution
               if( printDiff ) 
                  cout<<setw(nwdf+2)<< ( yout[i] /  ppGVf )
                                     - ( youtX[i] / ppRef );
              #endif
               cout<<setw(nwdf) << yout[i] / ppGVf ;
               if( printErr ) cout<<setw(nwdf)<< yerr[i] / ppGVf ;

            }
        cout.unsetf (std::ios_base::scientific);
        
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
              cout.setf (std::ios_base::fixed);
           }else if( i == 3 ){
              cout.unsetf (std::ios_base::fixed);              
              cout.setf (std::ios_base::scientific);
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
           // if( i == 2 )     { cout.unsetf(std::ios_base::fixed);      }
           // else if ( i==5 ) { cout.unsetf(std::ios_base::scientific); }
        }
        cout.unsetf(std::ios_base::scientific);
        if( j > 0 )  // Step 0 did not move -- printed the starting values
        {
           // cout.unsetf(std::ios_base::scientific);
           cout.setf (std::ios_base::fixed);                         
           cout.precision(2);
           cout<<setw(nwdf) << atan2(yout[1],yout[0])/degree;
           
          #ifdef BASELINESTEPPER                   // atan(yout[1]/yout[0])/degree;
           if( printRef ) 
             cout<<setw(nwdf) << atan2(youtX[1],youtX[0])/degree;
           #endif

           // Prepare the state of yIn(Vec)
           //Copy yout into yIn
           for(unsigned int i=0; i<Nposmom;i++){
              yInVec[i]  = yOutVec[i];
              #ifdef BASELINESTEPPER
              yInVecX[i] = yOutVecX[i];
              #endif
           }
        }
        
        cout<<"\n";
    }
    
    /*-----------------END-STEPPING------------------*/

    if( debug ) cout<<"----Stepping done "<<endl;

    /*------ Clean up ------*/
    #ifdef BASELINESTEPPER
    exactStepper->InformDone();
    #endif 

    // delete myStepper;  // Only if object was new'ed !!

    if( debug )  cout<<"----deletion of stepper done "<<endl;

    #ifdef BASELINESTEPPER
    delete exactStepper;
    #endif 
    // delete gvEquation;  // The stepper now takes ownership of the equation
    // delete gvEquation2;    
    delete gvUniformField;
    
    cout<<"\n\n#-------------End of output-----------------\n";
}
