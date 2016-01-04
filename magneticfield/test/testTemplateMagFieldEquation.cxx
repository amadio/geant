//
//
#include "base/Vector3D.h"

#include "TemplateGUVEquationOfMotion.h"
#include "TemplateGUVMagneticField.h"
#include "TemplateGUVField.h"
#include "TemplateTMagFieldEquation.h"
#include "TemplateTUniformMagField.h"

//To be able to call vecgeom::kScalar/kVc/kCuda
#include <Vc/Vc>
// #include "backend/vc/Backend.h"
// #include "backend/vcfloat/Backend.h"
// #include "backend/scalarfloat/Backend.h"


using ThreeVector_d = vecgeom::Vector3D<double>;

typedef vecgeom::Vector3D<typename vecgeom::kVc::precision_v> ThreeVectorSimd_d;

template <class Backend>
TemplateGUVEquationOfMotion<Backend>*  CreateFieldAndEquation(ThreeVector_d constFieldValue);

template <class Backend>
bool TestEquation(TemplateGUVEquationOfMotion<Backend>* );

constexpr unsigned int gNposmom= 6; // Position 3-vec + Momentum 3-vec

ThreeVector_d  FieldValue(0.0, 0.0, 1.0);
// ThreeVector_f  FieldValue(1.0, 1.0, 1.0);

int
main( int, char** )
{
  TemplateGUVEquationOfMotion<vecgeom::kVc>* eq = CreateFieldAndEquation<vecgeom::kVc>( FieldValue );

  TestEquation(eq);

  // TemplateTUniformMagField<vecgeom::kScalar>*   ConstBfield = new TemplateTUniformMagField<vecgeom::kScalar>( FieldValue );  
  // using EquationType = TemplateTMagFieldEquation < vecgeom::kScalar, TemplateTUniformMagField<vecgeom::kScalar>, gNposmom>;
  // TemplateGUVEquationOfMotion<vecgeom::kScalar>* magEquation = new EquationType(ConstBfield);

  // TestEquation(magEquation);

  return 1;
}

template <class Type>
Type CustomAbs(Type T ){
  return vecgeom::VECGEOM_IMPL_NAMESPACE::Abs(T);
}


template <class Backend>
TemplateGUVEquationOfMotion<Backend>* CreateFieldAndEquation(ThreeVector_d FieldValue)
{
  TemplateTUniformMagField<Backend>*   ConstBfield = new TemplateTUniformMagField<Backend>( FieldValue );
  
  using EquationType = TemplateTMagFieldEquation < Backend, TemplateTUniformMagField<Backend>, gNposmom>;

  TemplateGUVEquationOfMotion<Backend>* magEquation = new EquationType(ConstBfield);
  
  return magEquation;
}

int gVerbose= 1;

template <class Backend>
bool TestEquation(TemplateGUVEquationOfMotion<Backend>* equation)
{
  using ThreeVector_v = typename vecgeom::Vector3D<typename Backend::precision_v>; 
  using      Double_v = typename Backend::precision_v;
  using        Bool_v = typename Backend::bool_v;
   
  constexpr double perMillion = 1e-6;
  bool   hasError = false;  // Return value
  
  // vecgeom::Vector3D<typename Backend::precision_v> 
  ThreeVector_v PositionVec( 1., 2.,  3.);  // initial
  ThreeVector_v MomentumVec( 0., 0.1, 1.);
  ThreeVector_v FieldVec   ( 0., 0.,  1.);  // Magnetic field value (constant)

  // Double_v PositionTime[4] = { PositionVec.x(), PositionVec.y(), PositionVec.z(), 0.0};

  Double_v dydx[gNposmom];
  Double_v PositionMomentum[gNposmom];

  // double charge= -1;  
  Double_v ChargeVec( 1. );

  PositionMomentum[0] = PositionVec[0];
  PositionMomentum[1] = PositionVec[1];
  PositionMomentum[2] = PositionVec[2];
  PositionMomentum[3] = MomentumVec[0];
  PositionMomentum[4] = MomentumVec[1];
  PositionMomentum[5] = MomentumVec[2];

  // equation->InitializeCharge( charge );

  equation->EvaluateRhsGivenB( PositionMomentum, FieldVec, ChargeVec, dydx );

  vecgeom::Vector3D<typename Backend::precision_v>  ForceVec( dydx[3], dydx[4], dydx[5]);

  // Check result
  Double_v MdotF = MomentumVec.Dot(ForceVec);
  vecgeom::Vector3D<typename Backend::precision_v>  doubleFieldVec;
  // for(int i=0; i<3;i++) doubleFieldVec[i] = (Double_v) FieldVec[i];

  Double_v BdotF = FieldVec.Dot(ForceVec);

  Double_v momentumMag = MomentumVec.Mag();
  Double_v fieldMag =   FieldVec.Mag();
  Double_v ForceMag =   ForceVec.Mag();

  cout<<"\n";
  cout<<" momentumMag: "<<momentumMag<<endl;
  cout<<" fieldMag:    "<<fieldMag<<endl;
  cout<<" ForceMag:    "<<ForceMag<<endl;
  cout<<" abs(BdotF):  "<<CustomAbs(BdotF)<<endl;
  cout<<" abs(MdotF):  "<<CustomAbs(MdotF)<<endl;

  bool cond1 = !( vecgeom::IsEmpty( ForceMag != momentumMag*fieldMag ) );
  bool cond2 = !( vecgeom::IsEmpty( CustomAbs(MdotF) > perMillion * MomentumVec.Mag() * ForceVec.Mag() ) );
  bool cond3 = !( vecgeom::IsEmpty( CustomAbs(BdotF) > perMillion * FieldVec.Mag()    * ForceVec.Mag() ) );

  Bool_v print1 = (ForceMag != momentumMag*fieldMag) ;
  Bool_v print2 = CustomAbs(MdotF) > perMillion * MomentumVec.Mag() * ForceVec.Mag();
  Bool_v print3 = CustomAbs(BdotF) > perMillion * FieldVec.Mag()    * ForceVec.Mag();

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
     cout<<CustomAbs(MdotF)<<" vs "<< perMillion * MomentumVec.Mag()    * ForceVec.Mag() <<"\n"<<endl;
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

    cout<<"\n"<<CustomAbs(BdotF)<<" vs "<< perMillion * FieldVec.Mag()  * ForceVec.Mag() <<endl;
  }
  else if ( gVerbose )
  {
    std::cout << " Success:  Good (near zero) dot product magnetic-field . force " << std::endl;
  }

  return hasError;
  //return true;
}
