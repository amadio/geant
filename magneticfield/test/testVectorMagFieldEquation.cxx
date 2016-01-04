
#include "base/Vector3D.h"

#include "GUVEquationOfMotion.h"

#include "GUVVectorEquationOfMotion.h"
#include "GUVVectorField.h"
#include "TVectorMagFieldEquation.h"
#include "GUVMagneticField.h"
#include "GUVVectorMagneticField.h"
#include "TVectorUniformMagField.h"
#include "GUVVectorIntegrationStepper.h"
#include "GUVectorLineSection.h"
#include "ConstVectorFieldHelixStepper.h"
#include "GUVVectorHelicalStepper.h"

#include "GUVField.h"
#include "TMagFieldEquation.h"
#include "TUniformMagField.h"

//include template header files for testing
#include "TemplateGUVField.h"
#include "TemplateGUVMagneticField.h"
#include "TemplateTUniformMagField.h"
#include "TemplateGUVEquationOfMotion.h"
#include "TemplateTMagFieldEquation.h"
#include "TemplateGULineSection.h"
#include "TemplateGUVIntegrationStepper.h"
// #include "TemplateConstFieldHelixStepper.h" //nothing templatized really
#include "TemplateGUFieldTrack.h"
#include "TemplateGUVHelicalStepper.h"
#include "TemplateGUTCashKarpRKF45.h"
// #include "TemplateGUIntegrationDriver.h" //nothing done yet. Also, need to append stuff from .cxx file
#include "TemplateTMagErrorStepper.h"
#include "TemplateGUExactHelixStepper.h"
#include "TemplateTClassicalRK4.h"
#include "TemplateTSimpleRunge.h"
#include "TemplateStepperFactory.h"
#include "TemplateFieldEquationFactory.h"

// #include "ToyClass.h"

#include "Units.h"

// using fieldUnits::meter;
using fieldUnits::millimeter;   
using fieldUnits::second;  
using fieldUnits::eplus;  
using fieldUnits::tesla;
using fieldUnits::degree;

using namespace std;

typedef typename vecgeom::kVc::precision_v Double_v;
typedef typename vecgeom::kVcFloat::precision_v Float_v;
// typedef typename vecgeom::kVc::precision_v Float_v;      // Was: vecgeom::kVcFloat::precision_v 

typedef typename vecgeom::kVc::int_v Int_v;
typedef typename vecgeom::kVc::bool_v Bool_v;

using ThreeVector_f = vecgeom::Vector3D<float>;
using ThreeVector_d = vecgeom::Vector3D<double>;

using ThreeVectorSimd_f = vecgeom::Vector3D<Float_v>;
using ThreeVectorSimd_d = vecgeom::Vector3D<Double_v>;


GUVVectorEquationOfMotion*  CreateFieldAndEquation(ThreeVector_f constFieldValue);
bool  TestEquation(GUVVectorEquationOfMotion* );

constexpr unsigned int gNposmom= 6; // Position 3-vec + Momentum 3-vec

ThreeVectorSimd_f  FieldValue (0.0, 0.0, 1.0);
ThreeVector_f      FieldValue1(0.0, 0.0, 1.0);
ThreeVector_f      FieldValue2(1.0, 2.0, 3.0);
ThreeVectorSimd_f  FieldValueV(1.0, 2.0, 3.0);


int main(int /*argc*/, char ** /* *argv[] */ )
{
  GUVVectorEquationOfMotion* eq = CreateFieldAndEquation( FieldValue1 );
  TestEquation(eq);

  return 1;
}


GUVVectorEquationOfMotion* CreateFieldAndEquation(ThreeVector_f FieldValue)
{
  // std::cout<<"Creating field and equation..."<<std::endl;
  ThreeVector_f FieldValue2(1.0, 2.0, 3.0);

  TVectorUniformMagField* ConstBfield = new TVectorUniformMagField( FieldValue );
  // cout<<"--- TVectorUniformMagField object created"<<endl;
  using EquationType = TVectorMagFieldEquation<TVectorUniformMagField, gNposmom>;
  // cout<<"--- EquationType declared"<<endl;
  // GUVVectorEquationOfMotion*  magEquation= new TVectorMagFieldEquation<TUniformMagField, gNposmom>(ConstBfield);
  GUVVectorEquationOfMotion*  magEquation = new EquationType(ConstBfield);
  // cout<<"--- GUVVectorEquationOfMotion object created"<<endl;
  return magEquation;
  // pField = ConstBfield; 
  
}

int gVerbose= 1;

bool TestEquation(GUVVectorEquationOfMotion* equation)
{
  constexpr double perMillion = 1e-6;
  bool   hasError = false;  // Return value
  
  ThreeVectorSimd_d PositionVec( 1., 2.,  3.);  // initial
  ThreeVectorSimd_d MomentumVec( 0., 0.1, 1.);
  ThreeVectorSimd_f FieldVec( 0., 0., 1.);  // Magnetic field value (constant)

  // Double_v PositionTime[4] = { PositionVec.x(), PositionVec.y(), PositionVec.z(), 0.0};

  Double_v dydx[gNposmom];
  Double_v PositionMomentum[gNposmom];

  double charge= -1;  
  Double_v Charge(-1);

  PositionMomentum[0] = PositionVec[0];
  PositionMomentum[1] = PositionVec[1];
  PositionMomentum[2] = PositionVec[2];
  PositionMomentum[3] = MomentumVec[0];
  PositionMomentum[4] = MomentumVec[1];
  PositionMomentum[5] = MomentumVec[2];

  
  equation->InitializeCharge( charge );

  equation->EvaluateRhsGivenB( PositionMomentum, FieldVec,  Charge,   dydx );

  ThreeVectorSimd_d  ForceVec( dydx[3], dydx[4], dydx[5]);

  // Check result
  Double_v MdotF = MomentumVec.Dot(ForceVec);
  ThreeVectorSimd_d doubleFieldVec;
  for(int i=0; i<3;i++) doubleFieldVec[i] = (Double_v) FieldVec[i];

  Double_v BdotF = doubleFieldVec.Dot(ForceVec);

  Double_v momentumMag = MomentumVec.Mag();
  Double_v fieldMag =   doubleFieldVec.Mag();
  Double_v ForceMag =   ForceVec.Mag();

  cout<<"\n";
  cout<<" momentumMag: "<<momentumMag<<endl;
  cout<<" fieldMag:    "<<fieldMag<<endl;
  cout<<" ForceMag:    "<<ForceMag<<endl;
  cout<<" abs(BdotF):  "<<Vc::abs(BdotF)<<endl;
  cout<<" abs(MdotF):  "<<Vc::abs(MdotF)<<endl;

  bool cond1 = !( vecgeom::IsEmpty( ForceMag != momentumMag*fieldMag ) );
  bool cond2 = !( vecgeom::IsEmpty( Vc::abs(MdotF) > perMillion * MomentumVec.Mag()    * ForceVec.Mag() ) );
  bool cond3 = !( vecgeom::IsEmpty( Vc::abs(BdotF) > perMillion * doubleFieldVec.Mag() * ForceVec.Mag() ) );

  Bool_v print1 = (ForceMag != momentumMag*fieldMag) ;
  Bool_v print2 = Vc::abs(MdotF) > perMillion * MomentumVec.Mag()    * ForceVec.Mag();
  Bool_v print3 = Vc::abs(BdotF) > perMillion * doubleFieldVec.Mag() * ForceVec.Mag();

  cout<<"\nPrinting Bool_v "<<endl;
  cout<< print1 <<" and "<<cond1<<endl;
  cout<< print2 <<" and "<<cond2<<endl;
  cout<< print3 <<" and "<<cond3<<endl;
  cout<<"\n"<<endl;


  if( cond1 ) {
     std::cerr << "ERROR: Force magnitude is not equal to momentum * field."  << std::endl;  
     std::cout<< ForceMag<<" vs "<< momentumMag*fieldMag<<"\n"<<endl;   
  }
     
  assert( cond1 );  // Must add coefficient !!
  
  if( cond2 )
  { 
     std::cerr << "ERROR: Force due to magnetic field is not perpendicular to momentum!!"  << std::endl;
     hasError= true;
     cout<<Vc::abs(MdotF)<<" vs "<< perMillion * MomentumVec.Mag()    * ForceVec.Mag() <<"\n"<<endl;
  }
  else if ( gVerbose )
  {
     std::cout << " Success:  Good (near zero) dot product momentum . force " << std::endl;
  }

  if( cond3 )
  { 
    std::cerr << "ERROR: Force due to magnetic field is not perpendicular to B field!"
              << std::endl; 
    std::cerr << " Vectors:  BField   Force " << std::endl;
    for ( int i = 0; i < 3; i ++ )
      std::cerr << "   [" << i << "] " << FieldVec[i] << " " << ForceVec[i] << std::endl;

    hasError= true;

    cout<<"\n"<<Vc::abs(BdotF)<<" vs "<< perMillion * doubleFieldVec.Mag()  * ForceVec.Mag() <<endl;
  }
  else if ( gVerbose )
  {
    std::cout << " Success:  Good (near zero) dot product magnetic-field . force " << std::endl;
  }

  return hasError;
  //return true;
}
