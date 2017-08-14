#include "GlauberGribovTotalXsc.h"
#include "Proton.h"
#include "Neutron.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

#include <cmath>
#include <iostream>

namespace geantphysics{
  
  GlauberGribovTotalXsc::GlauberGribovTotalXsc() 
  {
    this->SetName("GlauberGribovTotalXsc");
    
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


  GlauberGribovTotalXsc::~GlauberGribovTotalXsc()
  {}

  ////////////////////////////////////////////////////////////////////////////////////////
  //
  // Calculates total cross section according to the Glauber model with Gribov correction calculated in the dipole approximation on
  // light cone. Gaussian density of point-like nucleons helps to calculate rest integrals of the model.
  // [1] B.Z. Kopeliovich, nucl-th/0306044 + simplification above

  double GlauberGribovTotalXsc::GetIsotopeCrossSection(const int particleCode, const double energyKin, const double mass,
						       const int Z, const int A)
  {
    double xsection, sigma, cofTotal, nucleusSquare, ratio;
    int N = A - Z;

    int particlePDG = Particle::GetParticleByInternalCode(particleCode)->GetPDGCode(); 


    if( A > 1 )
      { 
	double R = 0; 
  
	if( particlePDG == 321   || 
	    particlePDG == -321  || 
	    particlePDG == 310   || 
	    particlePDG == 130    )
	  {
	    sigma = Z * GetKaonNucleonTotalXscGG(particlePDG, mass, energyKin, 2212);
	    sigma += N * GetKaonNucleonTotalXscGG(particlePDG, mass, energyKin, 2112);
    
	    cofTotal = 2.0;
	    R = 1.3 * geant::fermi;
	    R *= std::pow(double(A), 0.3333);
	  }
	else
	  {	    
	    sigma = Z * GetHadronNucleonTotalXscNS(particlePDG, mass, energyKin, 2212);
	    sigma += N * GetHadronNucleonTotalXscNS(particlePDG, mass, energyKin, 2112);

	    cofTotal = 2.0;
	    R = GetNucleusRadius(A); 
	  }

	nucleusSquare = cofTotal * geant::kPi * R * R;   // basically 2piRR
	ratio = sigma / nucleusSquare;

	xsection =  nucleusSquare * std::log(1. + ratio);
	xsection *= GetParticleBarCorTot(particlePDG, Z);
      }
    else // H
      {
	if( particlePDG == 321   || 
	    particlePDG == -321  || 
	    particlePDG == 310   || 
	    particlePDG == 130    )
	  {
	    xsection = GetKaonNucleonTotalXscGG(particlePDG, mass, energyKin, 2212);
	  }
	else
	  {
	    xsection = GetHadronNucleonTotalXscNS(particlePDG, mass, energyKin, 2212);
	  }
      }
    return xsection; 
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Correction arrays for GG <-> Bar changea at ~ 90 GeV

  const double GlauberGribovTotalXsc::fNeutronBarCorrectionTot[93] = {

    1.0, 1.0,     1.42517e+00,  // 1.118517e+00, 
    1.082002e+00, 1.116171e+00, 1.078747e+00, 1.061315e+00, 
    1.058205e+00, 1.082663e+00, 1.068500e+00, 1.076912e+00, 1.083475e+00, 1.079117e+00, 
    1.071856e+00, 1.071990e+00, 1.073774e+00, 1.079356e+00, 1.081314e+00, 1.082056e+00,
    1.090772e+00, 1.096776e+00, 1.095828e+00, 1.097678e+00, 1.099157e+00, 1.103677e+00, 
    1.105132e+00, 1.109806e+00, 1.110816e+00, 1.117378e+00, 1.115165e+00, 1.115710e+00, 
    1.111855e+00, 1.110482e+00, 1.110112e+00, 1.106676e+00, 1.108706e+00, 1.105549e+00, 
    1.106318e+00, 1.106242e+00, 1.107672e+00, 1.107342e+00, 1.108119e+00, 1.106655e+00, 
    1.102588e+00, 1.096657e+00, 1.092920e+00, 1.086629e+00, 1.083592e+00, 1.076030e+00, 
    1.083777e+00, 1.089460e+00, 1.086545e+00, 1.079924e+00, 1.082218e+00, 1.077798e+00, 
    1.077062e+00, 1.072825e+00, 1.072241e+00, 1.072104e+00, 1.072490e+00, 1.069829e+00, 
    1.070398e+00, 1.065458e+00, 1.064968e+00, 1.060524e+00, 1.060048e+00, 1.057620e+00, 
    1.056428e+00, 1.055366e+00, 1.055017e+00, 1.052304e+00, 1.051767e+00, 1.049728e+00, 
    1.048745e+00, 1.047399e+00, 1.045876e+00, 1.042972e+00, 1.041824e+00, 1.039993e+00, 
    1.039021e+00, 1.036627e+00, 1.034176e+00, 1.032526e+00, 1.033633e+00, 1.036107e+00, 
    1.037803e+00, 1.031266e+00, 1.032991e+00, 1.033284e+00, 1.035015e+00, 1.033945e+00, 
    1.037075e+00, 1.034721e+00

  };


  const double GlauberGribovTotalXsc::fProtonBarCorrectionTot[93] = {

    1.0, 1.0,     
    1.118515e+00, 1.082000e+00, 1.116169e+00, 1.078745e+00, 1.061313e+00, 1.058203e+00, 
    1.082661e+00, 1.068498e+00, 1.076910e+00, 1.083474e+00, 1.079115e+00, 1.071854e+00, 
    1.071988e+00, 1.073772e+00, 1.079355e+00, 1.081312e+00, 1.082054e+00, 1.090770e+00, 
    1.096774e+00, 1.095827e+00, 1.097677e+00, 1.099156e+00, 1.103676e+00, 1.105130e+00, 
    1.109805e+00, 1.110814e+00, 1.117377e+00, 1.115163e+00, 1.115708e+00, 1.111853e+00, 
    1.110480e+00, 1.110111e+00, 1.106674e+00, 1.108705e+00, 1.105548e+00, 1.106317e+00, 
    1.106241e+00, 1.107671e+00, 1.107341e+00, 1.108118e+00, 1.106654e+00, 1.102586e+00, 
    1.096655e+00, 1.092918e+00, 1.086628e+00, 1.083590e+00, 1.076028e+00, 1.083776e+00, 
    1.089458e+00, 1.086543e+00, 1.079923e+00, 1.082216e+00, 1.077797e+00, 1.077061e+00, 
    1.072824e+00, 1.072239e+00, 1.072103e+00, 1.072488e+00, 1.069828e+00, 1.070396e+00, 
    1.065456e+00, 1.064966e+00, 1.060523e+00, 1.060047e+00, 1.057618e+00, 1.056427e+00, 
    1.055365e+00, 1.055016e+00, 1.052303e+00, 1.051766e+00, 1.049727e+00, 1.048743e+00, 
    1.047397e+00, 1.045875e+00, 1.042971e+00, 1.041823e+00, 1.039992e+00, 1.039019e+00, 
    1.036626e+00, 1.034175e+00, 1.032525e+00, 1.033632e+00, 1.036106e+00, 1.037802e+00, 
    1.031265e+00, 1.032990e+00, 1.033283e+00, 1.035014e+00, 1.033944e+00, 1.037074e+00, 
    1.034720e+00 

  };


  const double GlauberGribovTotalXsc::fPionPlusBarCorrectionTot[93] = {

    1.0, 1.0,     
    1.075927e+00, 1.074407e+00, 1.126098e+00, 1.100127e+00, 1.089742e+00, 1.083536e+00, 
    1.089988e+00, 1.103566e+00, 1.096922e+00, 1.126573e+00, 1.132734e+00, 1.136512e+00, 
    1.136629e+00, 1.133086e+00, 1.132428e+00, 1.129299e+00, 1.125622e+00, 1.126992e+00, 
    1.127840e+00, 1.162670e+00, 1.160392e+00, 1.157864e+00, 1.157227e+00, 1.154627e+00, 
    1.192555e+00, 1.197243e+00, 1.197911e+00, 1.200326e+00, 1.220053e+00, 1.215019e+00, 
    1.211703e+00, 1.209080e+00, 1.204248e+00, 1.203328e+00, 1.198671e+00, 1.196840e+00, 
    1.194392e+00, 1.193037e+00, 1.190408e+00, 1.188583e+00, 1.206127e+00, 1.210028e+00, 
    1.206434e+00, 1.204456e+00, 1.200547e+00, 1.199058e+00, 1.200174e+00, 1.200276e+00, 
    1.198912e+00, 1.213048e+00, 1.207160e+00, 1.208020e+00, 1.203814e+00, 1.202380e+00, 
    1.198306e+00, 1.197002e+00, 1.196027e+00, 1.195449e+00, 1.192563e+00, 1.192135e+00, 
    1.187556e+00, 1.186308e+00, 1.182124e+00, 1.180900e+00, 1.178224e+00, 1.176471e+00, 
    1.174811e+00, 1.173702e+00, 1.170827e+00, 1.169581e+00, 1.167205e+00, 1.165626e+00, 
    1.180244e+00, 1.177626e+00, 1.175121e+00, 1.173903e+00, 1.172192e+00, 1.171128e+00, 
    1.168997e+00, 1.166826e+00, 1.164130e+00, 1.165412e+00, 1.165504e+00, 1.165020e+00, 
    1.158462e+00, 1.158014e+00, 1.156519e+00, 1.156081e+00, 1.153602e+00, 1.154190e+00, 
    1.152974e+00
 
  };


  const double GlauberGribovTotalXsc::fPionMinusBarCorrectionTot[93] = {

    1.0, 1.0,     
    1.3956e+00, 1.077959e+00, 1.129145e+00, 1.102088e+00, 1.089765e+00, 1.083542e+00,  // 7
    1.089995e+00, 1.104895e+00, 1.097154e+00, 1.127663e+00, 1.133063e+00, 1.137425e+00, // 13
    1.136724e+00, 1.133859e+00, 1.132498e+00, 1.130276e+00, 1.127896e+00, 1.127656e+00, // 19
    1.127905e+00, 1.164210e+00, 1.162259e+00, 1.160075e+00, 1.158978e+00, 1.156649e+00, // 25 
    1.194157e+00, 1.199177e+00, 1.198983e+00, 1.202325e+00, 1.221967e+00, 1.217548e+00, 
    1.214389e+00, 1.211760e+00, 1.207335e+00, 1.206081e+00, 1.201766e+00, 1.199779e+00, 
    1.197283e+00, 1.195706e+00, 1.193071e+00, 1.191115e+00, 1.208838e+00, 1.212681e+00, 
    1.209235e+00, 1.207163e+00, 1.203451e+00, 1.201807e+00, 1.203283e+00, 1.203388e+00, 
    1.202244e+00, 1.216509e+00, 1.211066e+00, 1.211504e+00, 1.207539e+00, 1.205991e+00, 
    1.202143e+00, 1.200724e+00, 1.199595e+00, 1.198815e+00, 1.196025e+00, 1.195390e+00, 
    1.191137e+00, 1.189791e+00, 1.185888e+00, 1.184575e+00, 1.181996e+00, 1.180229e+00, 
    1.178545e+00, 1.177355e+00, 1.174616e+00, 1.173312e+00, 1.171016e+00, 1.169424e+00, 
    1.184120e+00, 1.181478e+00, 1.179085e+00, 1.177817e+00, 1.176124e+00, 1.175003e+00, 
    1.172947e+00, 1.170858e+00, 1.168170e+00, 1.169397e+00, 1.169304e+00, 1.168706e+00, 
    1.162774e+00, 1.162217e+00, 1.160740e+00, 1.160196e+00, 1.157857e+00, 1.158220e+00, 
    1.157267e+00 
  };


} // namespace geantphysics

