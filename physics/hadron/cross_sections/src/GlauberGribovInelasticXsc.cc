#include "GlauberGribovInelasticXsc.h"
#include "Proton.h"
#include "Neutron.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"
#include "Parameterizations.h"

#include <cmath>
#include <iostream>

namespace geantphysics {

  GlauberGribovInelasticXsc::GlauberGribovInelasticXsc() 
  {
    this->SetName("GlauberGribovInelasticXsc");

    std::vector< int > projVec;
    projVec.push_back(1);
    projVec.push_back(3);
    projVec.push_back(10);
    projVec.push_back(11);
    projVec.push_back(12);
    projVec.push_back(13);
    projVec.push_back(14);
    projVec.push_back(15);
    projVec.push_back(16);
    projVec.push_back(17);
    
    this->SetProjectileCodeVec(projVec);
  }


  GlauberGribovInelasticXsc::~GlauberGribovInelasticXsc()
  {}

  ////////////////////////////////////////////////////////////////////////////////////////
  //
  // Calculates the inelastic Xsc according to Glauber model with Gribov correction calculated in the dipole approximation on
  // light cone. Gaussian density of point-like nucleons helps to calculate rest integrals of the model.
  // [1] B.Z. Kopeliovich, nucl-th/0306044 + simplification above

  double GlauberGribovInelasticXsc::GetIsotopeCrossSection(const int particleCode, const double energyKin, const double mass,
							   const int Z, const int A)
  {
    double NucleusInelasticXsc, sigma, cofInelastic, cofTotal, nucleusSquare, ratio;

    int particlePDG = Particle::GetParticleByInternalCode(particleCode)->GetPDGCode(); 

    int N = A - Z;
    if( A > 1 )
      { 
	double R = 0; 

	if( particlePDG == 321   || 
	    particlePDG == -321  || 
	    particlePDG == 310   || 
	    particlePDG == 130    ) 
	  {
	    sigma = Z*GetKaonNucleonTotalXscGG(particlePDG, mass, energyKin, 2212);
	    sigma += N*GetKaonNucleonTotalXscGG(particlePDG, mass, energyKin, 2112);
    
	    cofInelastic = 2.2;
	    cofTotal     = 2.0;
	    R = 1.3*geant::fermi;
	    R *= std::pow(double(A), 0.3333);
	  }
	else
	  {
	    sigma = Z*GetHadronNucleonTotalXscNS(particlePDG, mass, energyKin, 2212);
	    sigma += N*GetHadronNucleonTotalXscNS(particlePDG, mass, energyKin, 2112);

	    cofInelastic = 2.4;
	    cofTotal     = 2.0;
	    R = GetNucleusRadius(A);   
	  }

	nucleusSquare = cofTotal*geant::kPi*R*R;   // basically 2piRR
	ratio = sigma/nucleusSquare;

	double fAxsc2piR2 = cofInelastic*ratio;

	double fModelInLog = std::log( 1. + fAxsc2piR2 );

	NucleusInelasticXsc = nucleusSquare*fModelInLog/cofInelastic;

	NucleusInelasticXsc *= GetParticleBarCorIn(particlePDG, Z);

      }
    else // H
      {
	if( particlePDG == 321   || 
	    particlePDG == -321  || 
	    particlePDG == 310   || 
	    particlePDG == 130    ) 
	  {
	    NucleusInelasticXsc = GetKaonNucleonInelasticXscGG(particlePDG, mass, energyKin, 2212);
	  }
	else
	  {
	    NucleusInelasticXsc = GetHadronNucleonInelasticXscNS(particlePDG, mass, energyKin, 2212);
	  }
      }
    return NucleusInelasticXsc; 
  }


  ///////////////////////////////////////////////////////////////////////////////
  //
  // Correction arrays for GG <-> Bar changea at ~ 90 GeV


  const double GlauberGribovInelasticXsc::fNeutronBarCorrectionIn[93] = {

    1.0, 1.0,     1.167421e+00, 1.156250e+00, 1.205364e+00, 1.154225e+00, 1.120391e+00, // 6
    1.124632e+00, 1.129460e+00, 1.107863e+00, 1.102152e+00, 1.104593e+00, 1.100285e+00, // 12
    1.098450e+00, 1.092677e+00, 1.101124e+00, 1.106461e+00, 1.115049e+00, 1.123903e+00, // 18
    1.126661e+00, 1.131259e+00, 1.133949e+00, 1.134185e+00, 1.133767e+00, 1.132813e+00, // 24
    1.131515e+00, 1.144338e+00, // 1.130338e+00, 
    1.134171e+00, 1.139206e+00, 1.148474e+00, // 1.141474e+00, 
    1.142189e+00, 
    1.140725e+00, 1.140100e+00, 1.139848e+00, 1.137674e+00, 1.138645e+00, 1.136339e+00, 
    1.136439e+00, 1.135946e+00, 1.136431e+00, 1.135702e+00, 1.135703e+00, 1.134113e+00, 
    1.131935e+00, 1.128381e+00, 1.126373e+00, 1.122453e+00, 1.120908e+00, 1.115953e+00, 
    1.115947e+00, 1.114426e+00, 1.111749e+00, 1.106207e+00, 1.107494e+00, 1.103622e+00, 
    1.102576e+00, 1.098816e+00, 1.097889e+00, 1.097306e+00, 1.097130e+00, 1.094578e+00, 
    1.094552e+00, 1.090222e+00, 1.089358e+00, 1.085409e+00, 1.084560e+00, 1.082182e+00, 
    1.080773e+00, 1.079464e+00, 1.078724e+00, 1.076121e+00, 1.075235e+00, 1.073159e+00, 
    1.071920e+00, 1.070395e+00, 1.069503e+00, 1.067525e+00, 1.066919e+00, 1.065779e+00, 
    1.065319e+00, 1.063730e+00, 1.062092e+00, 1.061085e+00, 1.059908e+00, 1.059815e+00, 
    1.059109e+00, 1.051920e+00, 1.051258e+00, 1.049473e+00, 1.048823e+00, 1.045984e+00, 
    1.046435e+00, 1.042614e+00

  };

  const double GlauberGribovInelasticXsc::fProtonBarCorrectionIn[93] = {

    1.0, 1.0,     
    1.147419e+00, // 1.167419e+00, 
    1.156248e+00, 1.205362e+00, 1.154224e+00, 1.120390e+00, 1.124630e+00, // 7 
    1.129459e+00, 1.107861e+00, 1.102151e+00, 1.104591e+00, 1.100284e+00, 1.098449e+00, // 13
    1.092675e+00, 1.101122e+00, 1.106460e+00, 1.115048e+00, 1.123902e+00, 1.126659e+00, // 19
    1.131258e+00, 1.133948e+00, 1.134183e+00, 1.133766e+00, 1.132812e+00, 1.131514e+00, // 25
    1.140337e+00, // 1.130337e+00, 

    1.134170e+00, 1.139205e+00, 1.151472e+00,  // 1.141472e+00, 
    1.142188e+00, 1.140724e+00, 
    1.140099e+00, 1.139847e+00, 1.137672e+00, 1.138644e+00, 1.136338e+00, 1.136438e+00, 
    1.135945e+00, 1.136429e+00, 1.135701e+00, 1.135702e+00, 1.134112e+00, 1.131934e+00, 
    1.128380e+00, 1.126371e+00, 1.122452e+00, 1.120907e+00, 1.115952e+00, 1.115946e+00, 
    1.114425e+00, 1.111748e+00, 1.106205e+00, 1.107493e+00, 1.103621e+00, 1.102575e+00, 
    1.098815e+00, 1.097888e+00, 1.097305e+00, 1.097129e+00, 1.094577e+00, 1.094551e+00, 
    1.090221e+00, 1.089357e+00, 1.085408e+00, 1.084559e+00, 1.082181e+00, 1.080772e+00, 
    1.079463e+00, 1.078723e+00, 1.076120e+00, 1.075234e+00, 1.073158e+00, 1.071919e+00, 
    1.070394e+00, 1.069502e+00, 1.067524e+00, 1.066918e+00, 1.065778e+00, 1.065318e+00, 
    1.063729e+00, 1.062091e+00, 1.061084e+00, 1.059907e+00, 1.059814e+00, 1.059108e+00, 
    1.051919e+00, 1.051257e+00, 1.049472e+00, 1.048822e+00, 1.045983e+00, 1.046434e+00, 
    1.042613e+00 

  };


  const double GlauberGribovInelasticXsc::fPionPlusBarCorrectionIn[93] = {

    1.0, 1.0,    
    1.140246e+00, 1.097872e+00, 1.104301e+00, 1.068722e+00, 1.056495e+00, 1.062622e+00, // 7
    1.047987e+00, 1.037032e+00, 1.035686e+00, 1.042870e+00, 1.052222e+00, 1.075100e+00, // 13
    1.084480e+00, 1.078286e+00, 1.081488e+00, 1.089713e+00, 1.099105e+00, 1.098003e+00, // 19
    1.102175e+00, 1.117707e+00, 1.121734e+00, 1.125229e+00, 1.126457e+00, 1.128905e+00, // 25
    1.163312e+00, 1.126263e+00, 1.126459e+00, 1.135191e+00, 1.116986e+00, 1.117184e+00, // 31
    1.117037e+00, 1.116777e+00, 1.115858e+00, 1.115745e+00, 1.114489e+00, 1.113993e+00, // 37
    1.113226e+00, 1.112818e+00, 1.111890e+00, 1.111238e+00, 1.111209e+00, 1.111775e+00, // 43
    1.110256e+00, 1.109414e+00, 1.107647e+00, 1.106980e+00, 1.106096e+00, 1.107331e+00, // 49
    1.107849e+00, 1.106407e+00, 1.103426e+00, 1.103896e+00, 1.101756e+00, 1.101031e+00, // 55
    1.098915e+00, 1.098260e+00, 1.097768e+00, 1.097487e+00, 1.095964e+00, 1.095773e+00, // 61
    1.093348e+00, 1.092687e+00, 1.090465e+00, 1.089821e+00, 1.088394e+00, 1.087462e+00, // 67
    1.086571e+00, 1.085997e+00, 1.084451e+00, 1.083798e+00, 1.082513e+00, 1.081670e+00, // 73
    1.080735e+00, 1.075659e+00, 1.074341e+00, 1.073689e+00, 1.072787e+00, 1.072237e+00, // 79
    1.071107e+00, 1.069955e+00, 1.074856e+00, 1.065873e+00, 1.065938e+00, 1.065694e+00, 
    1.062192e+00, 1.061967e+00, 1.061180e+00, 1.060960e+00, 1.059646e+00, 1.059975e+00, 
    1.059658e+00
 
  };


  const double GlauberGribovInelasticXsc::fPionMinusBarCorrectionIn[93] = {

    1.0, 1.0,    
    1.463e+00,    1.100898e+00, 1.106773e+00, 1.070289e+00, 1.040514e+00, 1.062628e+00, // 7
    1.047992e+00, 1.038041e+00, 1.035862e+00, 1.043679e+00, 1.052466e+00, 1.065780e+00, // 13
    1.070551e+00, 1.078869e+00, 1.081541e+00, 1.090455e+00, 1.100847e+00, 1.098511e+00, // 19 
    1.102226e+00, 1.118865e+00, 1.123143e+00, 1.126904e+00, 1.127785e+00, 1.130444e+00, // 25
    1.148502e+00, 1.127678e+00, 1.127244e+00, 1.123634e+00, 1.118347e+00, 1.118988e+00, 
    1.118957e+00, 1.118696e+00, 1.118074e+00, 1.117722e+00, 1.116717e+00, 1.116111e+00, 
    1.115311e+00, 1.114745e+00, 1.113814e+00, 1.113069e+00, 1.113141e+00, 1.113660e+00, 
    1.112249e+00, 1.111343e+00, 1.109718e+00, 1.108942e+00, 1.108310e+00, 1.109549e+00, 
    1.110227e+00, 1.108846e+00, 1.106183e+00, 1.106354e+00, 1.104388e+00, 1.103583e+00, 
    1.101632e+00, 1.100896e+00, 1.100296e+00, 1.099873e+00, 1.098420e+00, 1.098082e+00, 
    1.095892e+00, 1.095162e+00, 1.093144e+00, 1.092438e+00, 1.091083e+00, 1.090142e+00, 
    1.089236e+00, 1.088604e+00, 1.087159e+00, 1.086465e+00, 1.085239e+00, 1.084388e+00, 
    1.083473e+00, 1.078373e+00, 1.077136e+00, 1.076450e+00, 1.075561e+00, 1.074973e+00, 
    1.073898e+00, 1.072806e+00, 1.067706e+00, 1.068684e+00, 1.068618e+00, 1.068294e+00, 
    1.065241e+00, 1.064939e+00, 1.064166e+00, 1.063872e+00, 1.062659e+00, 1.062828e+00, 
    1.062699e+00 

  };


} // namespace geantphysics

