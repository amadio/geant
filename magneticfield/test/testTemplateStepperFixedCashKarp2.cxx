//
//  Compare the output of different steppers
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

#include <Vc/Vc>
#include "base/Vector3D.h"

#include "TemplateTUniformMagField.h"
#include "TemplateTMagFieldEquation.h"
#include "TemplateFieldEquationFactory.h"
#include "TemplateGUVIntegrationStepper.h"
#include "TemplateGUTCashKarpRKF45.h"

#include "TemplateGUIntegrationDriver.h"
#include "FieldTrack.h"

using namespace std;

#define DEBUGAnanya

int main(/*int argc, char *args[]*/)
{
    constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

    using Backend = vecgeom::kVc ;
    typedef typename Backend::precision_v Double;
    typedef vecgeom::Vector3D<Double> ThreeVectorSimd;
    typedef vecgeom::Vector3D<double> ThreeVector_d;

    using  GvEquationType=  TemplateTMagFieldEquation<Backend, TemplateTUniformMagField<Backend>, Nposmom>;
   
    /* -----------------------------SETTINGS-------------------------------- */
    
    /* Parameters of test
     - Modify values  */
    
    int no_of_steps = 250;         // No. of Steps for the stepper
    // int stepper_no =  5;         // Choose stepper no., for refernce see above
    double step_len_mm = 200.;    // meant as millimeter;  //Step length 
    double z_field_in = DBL_MAX;
    
    //Checking for command line values :
/*    if(argc>1)
        stepper_no = atoi(args[1]);
    if(argc > 2)
       step_len_mm = (float)(stof(args[2]));   // *mm);
    if(argc > 3)
        no_of_steps = atoi(args[3]);
    if(argc > 4)
       z_field_in = (float) (stof(args[4]));     // tesla*/
    double step_len = step_len_mm * fieldUnits::millimeter;
    
    //Set Charge etc.
    // double particleCharge = +1.0;      // in e+ units
    
    //Choice of output coordinates
    int
    columns[] =
    {
       1 , 1 , 1 ,  // position  x, y, z 
       1 , 1 , 1 ,  // momentum  x, y, z
       0 , 0 , 0 ,  // dydx pos  x, y, z
       1 , 1 , 0    // dydx mom  x, y, z
    }; //Variables in yOut[] & dydx[] we want to display - 0 for No, 1 for yes

    // bool printErr= 0;   // Print the predicted Error
    // bool printSep = 0;  // separator  '|'
    // bool printInp = 0;  // print the input values
    
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

     #ifdef DEBUGAnanya
      cout<<"----Just before making TemplateTUniformMagField"<<endl;
     #endif 

    // Field
    auto gvField= new TemplateTUniformMagField<Backend>( fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field) );
    #ifdef DEBUGAnanya
     cout<<"----TemplateTUniformMagField Object constructed"<<endl;
    #endif
    cout << "#  Initial  Field strength (GeantV) = "
         << x_field << " , " << y_field << " , " << z_field 
       // << (1.0/fieldUnits::tesla) * gvField->GetValue()->X() << ",  "
         << " Tesla " << endl;
    cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;
    //Create an Equation :
    #ifdef DEBUGAnanya
      cout<<"----Just before making EquationFactory"<<endl;
    #endif 
    auto gvEquation =
       TemplateFieldEquationFactory<Backend>::CreateMagEquation<TemplateTUniformMagField<Backend> >(gvField);
    #ifdef DEBUGAnanya
       cout<<"----EquationFactory made "<<endl;
    #endif 

    // gvEquation->InitializeCharge( particleCharge );  // Send it via Stepper instead    

    /*-------------------------PREPARING STEPPER-----------------------------*/
    
    //Create a stepper :

  #ifdef DEBUGAnanya
     cout<<"---- "<<endl;
     cout<<"---- Making TemplateGUTCashKarpRKF45"<<endl;
  #endif   
    // TemplateGUTCashKarpRKF45<Backend,GvEquationType,Nposmom> myStepper2(gvEquation);
  #ifdef DEBUGAnanya
    cout<<"---- constructed TemplateGUTCashKarpRKF45"<<endl;
  #endif
/*    TemplateGUVIntegrationStepper<Backend> *myStepper;
    myStepper = &myStepper2;*/

    TemplateGUVIntegrationStepper<Backend> *myStepper = new TemplateGUTCashKarpRKF45<Backend,GvEquationType,Nposmom>(gvEquation);

    // myStepper->InitializeCharge( particleCharge );

    //Initialising coordinates
    const double mmGVf = fieldUnits::millimeter;
    const double ppGVf = fieldUnits::GeV ;  //   it is really  momentum * c_light
                                         //   Else it must be divided by fieldUnits::c_light;
    // const double ppGVf = fieldUnits::GeV / Constants::c_light;     // OLD

    Double yIn[6]; //  = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                   //  x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};

    yIn[0] = x_pos * mmGVf ;
    yIn[1] = y_pos * mmGVf ;
    yIn[2] = z_pos * mmGVf ;
    yIn[3] = x_mom * ppGVf ;
    yIn[4] = y_mom * ppGVf ;
    yIn[5] = z_mom * ppGVf ;

    Double X_pos, Y_pos, Z_pos, X_mom, Y_mom, Z_mom;

    for (int i = 0; i < 4; ++i)
    {
      X_pos[i] = i-1.;
      Y_pos[i] = i-1.;
      Z_pos[i] = i-1.;
      X_mom[i] = i-1.;
      Y_mom[i] = i+1.-1.;
      Z_mom[i] = i+1.-1.;

    }
    cout<<"New X position is: "<<X_pos<<endl;

    // Double yIn[] = { X_pos * mmGVf, Y_pos * mmGVf ,Z_pos * mmGVf,
    //                  X_mom * ppGVf ,Y_mom * ppGVf ,Z_mom * ppGVf};
    double yOneIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                       x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};

    #ifdef DEBUGAnanya
      cout << yIn[0] << endl;
    #endif 
       

    // auto gvEquation2 = new GvEquationType(gvField);
                   // new TMagFieldEquation<TUniformMagField, Nposmom>(gvField);
    // gvEquation2->InitializeCharge( particleCharge ); // Let's make sure
    
    // Should be able to share the Equation -- eventually
    // For now, it checks that it was Done() -- and fails an assert


    std::cout << "# step_len_mm = " << step_len_mm;
    
    //Empty buckets for results
    Double dydx[8] = {0.,0.,0.,0.,0.,0.,0.,0.},  // 2 extra safety buffer
           yout[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
           yerr[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    
    /*-----------------------END PREPARING STEPPER---------------------------*/


    //=======================Test part for Integration driver====================
    auto testDriver = new TemplateGUIntegrationDriver<Backend>(0.2, myStepper);

    const ThreeVectorSimd  startPosition( yIn[0], yIn[1], yIn[2]);
    const ThreeVectorSimd  startMomentum( yIn[3], yIn[4], yIn[5]);
    TemplateGUFieldTrack<Backend>  yTrackIn(  startPosition, startMomentum );  // yStart
    TemplateGUFieldTrack<Backend>  yTrackOut( startPosition, startMomentum );  // yStart
    // Double total_step = 0.;

    typedef typename Backend::bool_v Bool_v;
    Bool_v succeded(true);
    double epsTol = 1.0e-5;

    // goodAdvance = testDriver->AccurateAdvance( yTrackIn, total_step, epsTol, yTrackOut );

    constexpr int nTracks = 16;
    FieldTrack yInput[nTracks];
    FieldTrack yOutput[nTracks];
    double hstep[nTracks];
    double charge[nTracks];
    bool   succeeded[nTracks];
    for (int i=0; i < nTracks; i++){
       charge[i] =  2.0 * ( i % 2 ) - 1.0;
       hstep[i]  = i;
       yInput[i] = FieldTrack( yOneIn ); // startPosition, startMomentum, 0.0 );
       yOutput[i]= FieldTrack( yOneIn ); // startPosition, startMomentum, 0.0 );
       succeded[i] = false;
    }
    
    testDriver->AccurateAdvance( yInput, charge, hstep, epsTol, yOutput, nTracks, succeeded );

    //========================End testing IntegrationDriver=======================


    /*---------------------------------------------------*/
    //        -> First Print the (commented) title header
    cout<<"\n#";
    cout<<setw(5)<<"Step";
    for (int i=0; i<6;i++)
        if (columns[i])
        {
           cout << setw(nwdf-2)<< "yOut[" << i << "]";
        }
    for (int i=0; i<6;i++)
        if (columns[i+6])
        {
           cout << setw(nwdf-2)<< "dydx[" << i << "]";
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
           const char* nameUnit = ( i<3 ) ? nameUnitLength : nameUnitMomentum ; 
           cout << setw(nwdf)<< nameUnit;
        }    
    cout<<"\n";
    
    //-> Then print the data

    /*--------------------Making a varied first input vector-------*/
/*    cout<<"Starting "<<endl;
    using  GvEquationTypeScalar=  TemplateTMagFieldEquation<vecgeom::kScalar, TemplateTUniformMagField<vecgeom::kScalar>, Nposmom>;
    auto gvFieldScalar= new TemplateTUniformMagField<vecgeom::kScalar>( fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field) );
    auto gvEquationScalar =
       TemplateFieldEquationFactory<vecgeom::kScalar>::CreateMagEquation<TemplateTUniformMagField<vecgeom::kScalar> >(gvFieldScalar);
    TemplateGUVIntegrationStepper<vecgeom::kScalar> *myStepperScalar;
    TemplateGUTCashKarpRKF45<vecgeom::kScalar,GvEquationTypeScalar,Nposmom> myStepper2Scalar(gvEquationScalar);
    myStepperScalar = &myStepper2Scalar;


    double yInScalar[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                    x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};
    double dydxScalar[8] = {0.,0.,0.,0.,0.,0.,0.,0.},  // 2 extra safety buffer
           youtScalar[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
           yerrScalar[8] = {0.,0.,0.,0.,0.,0.,0.,0.};

    double Charge = -1.;
    vecgeom::kVc::precision_v vY1, vY2, vY3, vY4, vY5;
    vecgeom::kVc::precision_v vDydx1, vDydx2, vDydx3, vDydx4, vDydx5;
    vecgeom::kVc::precision_v vY[6];
    vecgeom::kVc::precision_v vDydx[6];

    // vY[a][b] gets ath row, bth column


    for (int i = 0; i < 3; ++i)
    {
      cout<<"Made almost all objects 0"<<endl;
      myStepperScalar->RightHandSideVIS(yInScalar, Charge, dydxScalar );
      cout<<"Made almost all objects"<<endl;
      myStepperScalar->StepWithErrorEstimate(yInScalar, dydxScalar, step_len, youtScalar, yerrScalar);
      for (int k = 0; k < 6; ++k)
      {
        yInScalar[k] = youtScalar[k];
        cout<<"youtScalar is: "<<k<<" "<<youtScalar[k]<<endl;
        vY[k][i+1] = youtScalar[k];
        vDydx[k][i+1] = dydxScalar[k]; 
      }
    }

    for (int i = 0; i < 6; ++i)
    {
      yIn[i]  = vY[i];
      dydx[i] = vDydx[i];
      cout<<i<<" "<<yIn[i]<<endl;
      cout<<dydx[i]<<endl;
    }
*/

/*    vecgeom::kVc::precision_v vY0, vY1, vY2, vY3, vY4, vY5;
    vecgeom::kVc::precision_v vDydx1, vDydx2, vDydx3, vDydx4, vDydx5, vDydx6;

    vY0[0]= 0.;
    vY0[1]= -2.9975/10.;
    vY0[2]= -11.9845/10.;
    vY0[3]= -26.9450/10.;
    vY1[0]= 0./10.;
    vY1[1]= 141.379/10.;
    vY1[2]= 282.504/10.;
    vY1[3]=423.121/10.;
    vY2[0]= 0./10.;
    vY2[1]=141.421/10.;
    vY2[2]=282.843/10.;
    vY2[3]=424.264/10.;
    vY3[0]= 0./10.;
    vY3[1]=  -.04238/10.;
    vY3[2]=  -.08469/10.;
    vY3[3]=  -.01268/10.;
    vY4[0]=   0./10.;
    vY4[1]=  .9991/10.;
    vY4[2]=  .9964/10.;
    vY4[3]=  .9919/10.;
    vY5[0]=   0./10.;
    vY5[1]=  1./10.;
    vY5[2]=  1./10.;
    vY5[3]=  1./10.;

    yIn[0] = vY0;
    yIn[1] = vY1;
    yIn[2] = vY2;
    yIn[3] = vY3;
    yIn[4] = vY4;
    yIn[5] = vY5;*/
    



    /*-------------------Done making it-----------------------*/
    
    // cout.setf (ios_base::scientific);
    // cout.setf (ios_base::scientific);    
    cout.precision(3);
    
    /*----------------NOW STEPPING-----------------*/
    no_of_steps = 1;
    for(int j=0; j<no_of_steps; j++)
    {
        cout<<setw(6)<<j ;           //Printing Step number
        Double charge(-1.);


/*        #ifdef DEBUGAnanya
          cout<<"\n----y before evaluating RightHandSideVIS was: "<<yIn[3]<<endl;
        #endif */

        myStepper->RightHandSideVIS(yIn, charge, dydx);               //compute dydx - to supply the stepper

/*        #ifdef DEBUGAnanya
          cout<<"----dydx from RightHandSideVIS is: "<<dydx[3]<<endl;
          cout<<"----y given to RightHandSideVIS was: "<<yIn[3]<<endl;
        #endif */

        if( j > 0 )  // Let's print the initial points!
        {
           myStepper->StepWithErrorEstimate(yIn, dydx, charge, step_len, yout, yerr);   //Call the 'trial' stepper
        }
        //-> Then print the data
        cout.setf (ios_base::fixed);
        cout.precision(4);
        for(int i=0; i<3;i++)
            if(columns[i]){
               cout<<setw(nwdf)<< yout[i] / mmGVf ;
               cout.precision(3);
            }

        cout.unsetf (ios_base::fixed);        
        cout.setf (ios_base::scientific);
        for(int i=3; i<6;i++)
            if(columns[i]){
               cout<<setw(nwdf) << yout[i] / ppGVf ;
            }
        cout.unsetf (ios_base::scientific);
        
        for(int i=0; i<6;i++)   // Print auxiliary components
        {
           double unitGVf=1;  
           // double unitRef=1;
           // if( i < 3 )             // length / length

           if( i >= 3 ){
              unitGVf = ppGVf / mmGVf; //  dp / ds
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
               cout<<setw(nwdf)<< dydx[i] / unitGVf ;
           }
        }
        cout.unsetf(ios_base::scientific);
        if( j > 0 )  // Step 0 did not move -- printed the starting values
        {
           // cout.unsetf(ios_base::scientific);
           cout.setf (ios_base::fixed);                         
           cout.precision(2);
           cout<<setw(nwdf) << atan2(yout[1],yout[0])/degree;
           
           //Copy yout into yIn
           for(int i=0; i<6;i++){
              yIn[i] = yout[i];

           }
        }
        
        cout<<"\n";
    }
    
    /*-----------------END-STEPPING------------------*/

    #ifdef DEBUGAnanya
      cout<<"----Stepping done "<<endl;
    #endif 


    /*------ Clean up ------*/
    #ifdef DEBUGAnanya
      cout<<"----Informing done "<<endl;
    #endif 

    delete myStepper;

    #ifdef DEBUGAnanya
      cout<<"----deletion of stepper done "<<endl;
    #endif 
 
    delete gvField;
    
    cout<<"\n\n#-------------End of output-----------------\n";
}
