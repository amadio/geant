#include "HadronNucleusXsc.h"
#include "Proton.h"
#include "Neutron.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

#include <cmath>

namespace geantphysics {

HadronNucleusXsc::HadronNucleusXsc() 
{
}


HadronNucleusXsc::~HadronNucleusXsc()
{}

////////////////////////////////////////////////////////////////////////////////////////
//
// Calculates total and inelastic Xsc, derives elastic as total - inelastic according to
// Glauber model with Gribov correction calculated in the dipole approximation on
// light cone. Gaussian density of point-like nucleons helps to calculate rest integrals of the model.
// [1] B.Z. Kopeliovich, nucl-th/0306044 + simplification above

double HadronNucleusXsc::CalculateCrossSection(int particlePDG, double mass, double energyKin, int Z, int A)
{
  double xsection, sigma, cofInelastic, cofTotal, nucleusSquare, ratio;
  double hpInXsc(0.), hnInXsc(0.);
  double R = GetNucleusRadius(A); 
  
  int N = A - Z;              // number of neutrons
  if (N < 0) N = 0;

  if( particlePDG == 2212 || 
      particlePDG == 2112 ||
      particlePDG == 211  || 
      particlePDG == -211 )
  {
    sigma = Z*GetHadronNucleonXscNS(particlePDG, mass, energyKin, 2212);

    hpInXsc = fNucleonInelasticXsc;

    sigma += N*GetHadronNucleonXscNS(particlePDG, mass, energyKin, 2112);

    hnInXsc = fNucleonInelasticXsc;  

    cofInelastic = 2.4;
    cofTotal     = 2.0;
  }
  else if( particlePDG == 321   || 
           particlePDG == -321  || 
           particlePDG == 310   || 
           particlePDG == 130    ) 
  {
    sigma = Z*GetKaonNucleonXscGG(particlePDG, mass, energyKin, 2112);

    hpInXsc = fNucleonInelasticXsc;

    sigma += N*GetKaonNucleonXscGG(particlePDG, mass, energyKin, 2112);

    hnInXsc = fNucleonInelasticXsc;
    
    cofInelastic = 2.2;
    cofTotal     = 2.0;
    R = 1.3*geant::fermi;
    R *= std::pow(double(A), 0.3333);
  }
  else
  {
    sigma = Z*GetHadronNucleonXscNS(particlePDG, mass, energyKin, 2212);

    hpInXsc = fNucleonInelasticXsc;

    sigma += N*GetHadronNucleonXscNS(particlePDG, mass, energyKin, 2112);

    hnInXsc = fNucleonInelasticXsc;  

    cofInelastic = 2.4;
    cofTotal     = 2.0;
    
    //   sigma        = GetHadronNucleonXscNS(aParticle, A, Z);
    //   cofInelastic = 2.2;
    //   cofTotal     = 2.0;
  }
  // cofInelastic = 2.0;

  if( A > 1 )
  { 
    nucleusSquare = cofTotal*geant::kPi*R*R;   // basically 2piRR
    ratio = sigma/nucleusSquare;

    xsection =  nucleusSquare*std::log( 1. + ratio );

    xsection *= GetParticleBarCorTot(particlePDG, Z);

    fNucleusTotalXsc = xsection;

    // inelastic xsc

    double fAxsc2piR2 = cofInelastic*ratio;

    double fModelInLog = std::log( 1. + fAxsc2piR2 );

    fNucleusInelasticXsc = nucleusSquare*fModelInLog/cofInelastic;

    fNucleusInelasticXsc *= GetParticleBarCorIn(particlePDG, Z);

    fNucleusElasticXsc   = fNucleusTotalXsc - fNucleusInelasticXsc;

    if( fNucleusElasticXsc < 0. ) fNucleusElasticXsc = 0.;
    
    double difratio = ratio/(1.+ratio);

    fNucleusDiffractionXsc = 0.5*nucleusSquare*( difratio - std::log( 1. + difratio ) );

    sigma = Z*hpInXsc + N*hnInXsc;

    ratio = sigma/nucleusSquare;

    fNucleusProductionXsc = nucleusSquare*std::log( 1. + cofInelastic*ratio )/cofInelastic;

    fNucleusProductionXsc *= GetParticleBarCorIn(particlePDG, Z);

    if (fNucleusElasticXsc < 0.) fNucleusElasticXsc = 0.;
  }
  else // H
  {
    fNucleusTotalXsc = sigma;
    xsection  = sigma;

    fNucleusInelasticXsc = fNucleonInelasticXsc;

    if ( particlePDG != -2212 ) 
    {
      fNucleusElasticXsc = fNucleonElasticXsc;
      
     // sigma         = GetHNinelasticXsc(aParticle, A, Z);
     // fNucleonInelasticXsc = sigma;
     // fNucleonElasticXsc   = fTotalXsc - fNucleonInelasticXsc;      
    }
    else if( particlePDG == 321   || 
             particlePDG == -321  || 
             particlePDG == 310   || 
             particlePDG == 130    ) 
    { 
      fNucleusInelasticXsc = hpInXsc;
      fNucleusElasticXsc   = fNucleusTotalXsc - fNucleusInelasticXsc;
    }   
    else
    {
      fNucleusInelasticXsc = hpInXsc;
      fNucleusElasticXsc   = fNucleusTotalXsc - fNucleusInelasticXsc;
    }
    if (fNucleusElasticXsc < 0.) fNucleusElasticXsc = 0.;
      
  }
  return xsection; 
}



/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon cross-section based on N. Starkov parametrisation of
// data from mainly http://wwwppds.ihep.su:8001/c5-6A.html database
  
double 
HadronNucleusXsc::GetHadronNucleonXscNS(int particlePDG, double mass, double energyKin, int targetPDG)
{
  double xsection(0); 
  
  double A0, B0;
  double hpXsc(0);
  double hnXsc(0);


  double tM = 0.939 * geant::GeV;  // ~mean neutron and proton ???

  double pM   = mass;
  double pE   = energyKin + mass; // total energy!!!!
  double pLab  = std::sqrt(energyKin*(energyKin+2.*pM));

  double sMand = CalcMandelstamS ( pM , tM , pLab );


  double logP = std::log(pLab);


  // General PDG fit constants

  double s0   = 5.38*5.38; // in Gev^2
  double eta1 = 0.458;
  double eta2 = 0.458;
  double B    = 0.308;
  double minLogP = 3.5;       // min of (lnP-minLogP)^2 
  double cofLogE = .0557;     // elastic (lnP-minLogP)^2 
  double cofLogT = .3;        // total (lnP-minLogP)^2 
  double pMin = .1;        // fast LE calculation 
  double pMax = 1000.;     // fast HE calculation 


  // proton = 2212, neutron = 2112
  bool pORn = (targetPDG == 2212 || targetPDG == 2112);  
  bool proton = (targetPDG == 2212);
  bool neutron = (targetPDG == 2112);

  
  if( particlePDG == 2112 && pORn ) 
  {
    if( pLab >= 373.)
    {
      xsection =  GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;

      fNucleonElasticXsc = 6.5 + 0.308*std::pow(std::log(sMand/400.),1.65) + 9.19*std::pow(sMand,-0.458);
    
      fNucleonTotalXsc = xsection;
    
    }
    else if( pLab >= 100.)
    {
      B0 = 7.5;
      A0 = 100. - B0*std::log(3.0e7);

      xsection = A0 + B0*std::log(pE) - 11
	  // + 103*std::pow(2*0.93827*pE + pM*pM+0.93827*0.93827,-0.165);        //  mb
                  + 103*std::pow(sMand,-0.165);        //  mb

      fNucleonElasticXsc = 5.53 + 0.308*std::pow(std::log(sMand/28.9),1.1) + 9.19*std::pow(sMand,-0.458);
      
      fNucleonTotalXsc = xsection;
    }
    else if( pLab >= 10.)
    {
        B0 = 7.5;
        A0 = 100. - B0*std::log(3.0e7);

        xsection = A0 + B0*std::log(pE) - 11
                  + 103*std::pow(2*0.93827*pE + pM*pM+
                     0.93827*0.93827,-0.165);        //  mb      
      fNucleonTotalXsc = xsection;
      fNucleonElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
    }
    else  // pLab < 10 GeV/c
    {
      if( neutron )      // nn to be pp
      {
        if( pLab < 0.4 )
        {
          hnXsc = 23 + 50*( std::pow( std::log(0.73/pLab), 3.5 ) );
          fNucleonElasticXsc = hnXsc;
        }
        else if( pLab < 0.73 )
        {
          hnXsc = 23 + 50*( std::pow( std::log(0.73/pLab), 3.5 ) );
          fNucleonElasticXsc = hnXsc; 
        }
        else if( pLab < 1.05  )
        {
          hnXsc = 23 + 40*(std::log(pLab/0.73))*
                         (std::log(pLab/0.73));
          fNucleonElasticXsc = 23 + 20*(std::log(pLab/0.73))*
                         (std::log(pLab/0.73));
        }
        else    // 1.05 - 10 GeV/c
        {
          hnXsc = 39.0+75*(pLab - 1.2)/(std::pow(pLab,3.0) + 0.15);

          fNucleonElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
        }
        fNucleonTotalXsc = hnXsc;
      }
      if( proton )   // pn to be np
      {
        if( pLab < 0.02 )
        {
          hpXsc = 4100+30*std::pow(std::log(1.3/pLab),3.6); // was as pLab < 0.8
	  fNucleonElasticXsc = hpXsc;
        }      
        else if( pLab < 0.8 )
        {
          hpXsc = 33+30*std::pow(std::log(pLab/1.3),4.0);
	  fNucleonElasticXsc = hpXsc;
        }      
        else if( pLab < 1.05 )
        {
          hpXsc = 33+30*std::pow(std::log(pLab/0.95),2.0);
          fNucleonElasticXsc =  6 + 52/( std::log(0.511/pLab)*std::log(0.511/pLab) + 1.6 );
        }
        else if( pLab < 1.4 )
        {
          hpXsc = 33+30*std::pow(std::log(pLab/0.95),2.0);
          fNucleonElasticXsc =  6 + 52/( std::log(0.511/pLab)*std::log(0.511/pLab) + 1.6 );
        }
        else    // 1.4 < pLab < 10.  )
        {
          hpXsc = 33.3 + 20.8*(std::pow(pLab,2.0) - 1.35)/(std::pow(pLab,2.50) + 0.95);
          
          fNucleonElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
        }
        fNucleonTotalXsc = hpXsc;
      }
    }
  } 
  else if( particlePDG == 2212 && pORn ) ////// proton //////////////////////////////////////////////
  {
    if( pLab >= 373.) // pdg due to TOTEM data
    {
      xsection =  GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;

      fNucleonElasticXsc = 6.5 + 0.308*std::pow(std::log(sMand/400.),1.65) + 9.19*std::pow(sMand,-0.458);
     
      fNucleonTotalXsc = xsection;
    }
    else if( pLab >= 100.)
    {
      B0 = 7.5;
      A0 = 100. - B0*std::log(3.0e7);

      xsection = A0 + B0*std::log(pE) - 11 + 103*std::pow(sMand,-0.165);        //  mb

      fNucleonElasticXsc = 5.53 + 0.308*std::pow(std::log(sMand/28.9),1.1) + 9.19*std::pow(sMand,-0.458);
      
      fNucleonTotalXsc = xsection;
    }
    else if( pLab >= 10.)
    {
      B0 = 7.5;
      A0 = 100. - B0*std::log(3.0e7);

      xsection = A0 + B0*std::log(pE) - 11 + 103*std::pow(sMand,-0.165);        //  mb

      fNucleonElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
      
      fNucleonTotalXsc = xsection;
    }
    else
    {
      // pp

      if( proton )
      {
        if( pLab < 0.4 )
        {
          hpXsc = 23 + 50*( std::pow( std::log(0.73/pLab), 3.5 ) );
          fNucleonElasticXsc = hpXsc;
        }
        else if( pLab < 0.73 )
        {
          hpXsc = 23 + 50*( std::pow( std::log(0.73/pLab), 3.5 ) );
          fNucleonElasticXsc = hpXsc; 
        }
        else if( pLab < 1.05  )
        {
          hpXsc = 23 + 40*(std::log(pLab/0.73))*
                         (std::log(pLab/0.73));
          fNucleonElasticXsc = 23 + 20*(std::log(pLab/0.73))*
                         (std::log(pLab/0.73));
        }
        else    // 1.05 - 10 GeV/c
        {
          hpXsc = 39.0+75*(pLab - 1.2)/(std::pow(pLab,3.0) + 0.15);

          fNucleonElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
        }
        fNucleonTotalXsc = hpXsc;
      }
      if( neutron )     // pn to be np
      {
        if( pLab < 0.02 )
        {
          hnXsc = 4100+30*std::pow(std::log(1.3/pLab),3.6); // was as pLab < 0.8
	  fNucleonElasticXsc = hnXsc;
        }      
        else if( pLab < 0.8 )
        {
          hnXsc = 33+30*std::pow(std::log(pLab/1.3),4.0);
	  fNucleonElasticXsc = hnXsc;
        }      
        else if( pLab < 1.05 )
        {
          hnXsc = 33+30*std::pow(std::log(pLab/0.95),2.0);
          fNucleonElasticXsc =  6 + 52/( std::log(0.511/pLab)*std::log(0.511/pLab) + 1.6 );
        }
        else if( pLab < 1.4 )
        {
          hnXsc = 33+30*std::pow(std::log(pLab/0.95),2.0);
          fNucleonElasticXsc =  6 + 52/( std::log(0.511/pLab)*std::log(0.511/pLab) + 1.6 );
        }
        else    // 1.4 < pLab < 10.  )
        {
          hnXsc = 33.3 + 20.8*(std::pow(pLab,2.0) - 1.35)/(std::pow(pLab,2.50) + 0.95);
          
          fNucleonElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
        }
        fNucleonTotalXsc = hnXsc;
      }
    }    
  } 
  else if( particlePDG == -2212 && pORn ) /////////////////// p_bar ///////////////////////////
  {
    if( proton )
    {
      xsection  = 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) + 33.34*std::pow(sMand,-eta2);
    }
    if( neutron ) // ???
    {
      xsection = 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) + 30.*std::pow(sMand,-eta2);
    }
    fNucleonTotalXsc = xsection;
  } 
  else if( particlePDG == 211 && pORn ) // pi+ /////////////////////////////////////////////
  {
    if( proton ) // pi+ p
    {
      if( pLab < 0.28 )
      {
        hpXsc       = 10./((logP + 1.273)*(logP + 1.273) + 0.05);
        fNucleonElasticXsc = hpXsc;
      }
      else if( pLab < 0.4 )
      {
        hpXsc       = 14./( (logP + 1.273)*(logP + 1.273) + 0.07);
        fNucleonElasticXsc = hpXsc;
      }
      else if( pLab < 0.68 )
      {
        hpXsc       = 14./( (logP + 1.273)*(logP + 1.273) + 0.07);
        fNucleonElasticXsc = hpXsc;
      }
      else if( pLab < 0.85 )
      {
        double Ex4 = 88*(std::log(pLab/0.77))*(std::log(pLab/0.77));
        hpXsc        = Ex4 + 14.9;
        fNucleonElasticXsc = hpXsc*std::exp(-3.*(pLab - 0.68));  
      }
      else if( pLab < 1.15 )
      {
        double Ex4 = 88*(std::log(pLab/0.77))*(std::log(pLab/0.77));
        hpXsc        = Ex4 + 14.9;

        fNucleonElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
      }
      else if( pLab < 1.4) // ns original
      {
        double Ex1 = 3.2*std::exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        double Ex2 = 12*std::exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hpXsc        = Ex1 + Ex2 + 27.5;
        fNucleonElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
      }
      else if( pLab < 2.0 ) // ns original
      {
        double Ex1 = 3.2*std::exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        double Ex2 = 12*std::exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hpXsc        = Ex1 + Ex2 + 27.5;
        fNucleonElasticXsc = 3.0 + 1.36/( (logP - 0.336)*(logP - 0.336) + 0.08);    
      }
      else if( pLab < 3.5 ) // ns original
      {
        double Ex1 = 3.2*std::exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        double Ex2 = 12*std::exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hpXsc        = Ex1 + Ex2 + 27.5;
        fNucleonElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      else if( pLab < 200. ) // my
      {
        hpXsc = 10.6 + 2.*std::log(pE) + 25*std::pow(pE, -0.43 ); // ns original
        fNucleonElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      else //  pLab > 100 // my
      {
        hpXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
        fNucleonElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      fNucleonTotalXsc = hpXsc;
    }    
    if( neutron )  // pi+ n = pi- p??
    {
      if( pLab < 0.28 ) 
      {
        hnXsc       = 0.288/((pLab - 0.28)*(pLab - 0.28) + 0.004);
        fNucleonElasticXsc = 1.8/((logP + 1.273)*(logP + 1.273) + 0.07);
      }
      else if( pLab < 0.395676 ) // first peak
      {
        hnXsc       = 0.648/((pLab - 0.28)*(pLab - 0.28) + 0.009);
        fNucleonElasticXsc = 0.257/((pLab - 0.28)*(pLab - 0.28) + 0.01);
       }
      else if( pLab < 0.5 )
      {
        hnXsc       = 26 + 110*(std::log(pLab/0.48))*(std::log(pLab/0.48));
        fNucleonElasticXsc = 0.37*hnXsc;
      }
      else if( pLab < 0.65 )
      {
        hnXsc       = 26 + 110*(std::log(pLab/0.48))*(std::log(pLab/0.48));
        fNucleonElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 0.72 )
      {
        hnXsc = 36.1 + 10*std::exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*std::exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fNucleonElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 0.88 )
      {
        hnXsc = 36.1 + 10*std::exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*std::exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fNucleonElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 1.03 )
      {
        hnXsc = 36.1 + 10*std::exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*std::exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fNucleonElasticXsc = 2.0 + 0.4/((pLab - 1.03)*(pLab - 1.03) + 0.016);
      }
      else if( pLab < 1.15 )
      {
        hnXsc = 36.1 + 10*std::exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*std::exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fNucleonElasticXsc = 2.0 + 0.4/((pLab - 1.03)*(pLab - 1.03) + 0.016);
      }
      else if( pLab < 1.3 )
      {
        hnXsc = 36.1 + 10*std::exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*std::exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fNucleonElasticXsc = 3. + 13./pLab;
      }
      else if( pLab < 2.6 ) // < 3.0) // ns original
      {
        hnXsc = 36.1 + 0.079-4.313*std::log(pLab)+
                3*std::exp(-(pLab-2.1)*(pLab-2.1)/0.4/0.4)+
                1.5*std::exp(-(pLab-1.4)*(pLab-1.4)/0.12/0.12);
        fNucleonElasticXsc = 3. + 13./pLab; 
      }
      else if( pLab < 20. ) // < 3.0) // ns original
      {
        hnXsc = 36.1 + 0.079 - 4.313*std::log(pLab)+
                3*std::exp(-(pLab-2.1)*(pLab-2.1)/0.4/0.4)+
                1.5*std::exp(-(pLab-1.4)*(pLab-1.4)/0.12/0.12);
        fNucleonElasticXsc = 3. + 13./pLab; 
      }
      else   // mb 
      {
        hnXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
        fNucleonElasticXsc = 3. + 13./pLab;
      }
      fNucleonTotalXsc = hnXsc;
    }
  } 
  else if( particlePDG == -211 && pORn ) /// pi- ////////////////////////////////////////////
  {
    if( neutron )     // pi- n = pi+ p??
    {
      if( pLab < 0.28 )
      {
        hnXsc       = 10./((logP + 1.273)*(logP + 1.273) + 0.05);
        fNucleonElasticXsc = hnXsc;
      }
      else if( pLab < 0.4 )
      {
        hnXsc       = 14./( (logP + 1.273)*(logP + 1.273) + 0.07);
        fNucleonElasticXsc = hnXsc;
      }
      else if( pLab < 0.68 )
      {
        hnXsc       = 14./( (logP + 1.273)*(logP + 1.273) + 0.07);
        fNucleonElasticXsc = hnXsc;
      }
      else if( pLab < 0.85 )
      {
        double Ex4 = 88*(std::log(pLab/0.77))*(std::log(pLab/0.77));
        hnXsc        = Ex4 + 14.9;
        fNucleonElasticXsc = hnXsc*std::exp(-3.*(pLab - 0.68));  
      }
      else if( pLab < 1.15 )
      {
        double Ex4 = 88*(std::log(pLab/0.77))*(std::log(pLab/0.77));
        hnXsc        = Ex4 + 14.9;

        fNucleonElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
      }
      else if( pLab < 1.4) // ns original
      {
        double Ex1 = 3.2*std::exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        double Ex2 = 12*std::exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hnXsc        = Ex1 + Ex2 + 27.5;
        fNucleonElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
      }
      else if( pLab < 2.0 ) // ns original
      {
        double Ex1 = 3.2*std::exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        double Ex2 = 12*std::exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hnXsc        = Ex1 + Ex2 + 27.5;
        fNucleonElasticXsc = 3.0 + 1.36/( (logP - 0.336)*(logP - 0.336) + 0.08);    
      }
      else if( pLab < 3.5 ) // ns original
      {
        double Ex1 = 3.2*std::exp(-(pLab-2.55)*(pLab-2.55)/0.55/0.55);
        double Ex2 = 12*std::exp(-(pLab-1.47)*(pLab-1.47)/0.225/0.225);
        hnXsc        = Ex1 + Ex2 + 27.5;
        fNucleonElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      else if( pLab < 200. ) // my
      {
        hnXsc = 10.6 + 2.*std::log(pE) + 25*std::pow(pE, -0.43 ); // ns original
        fNucleonElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      else //  pLab > 100 // my
      {
        hnXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
        fNucleonElasticXsc = 3.0 + 6.20/( (logP - 0.336)*(logP - 0.336) + 0.8);    
      }
      fNucleonTotalXsc = hnXsc;
    }
    if( proton )    // pi- p
    {
      if( pLab < 0.28 ) 
      {
        hpXsc       = 0.288/((pLab - 0.28)*(pLab - 0.28) + 0.004);
        fNucleonElasticXsc = 1.8/((logP + 1.273)*(logP + 1.273) + 0.07);
      }
      else if( pLab < 0.395676 ) // first peak
      {
        hpXsc       = 0.648/((pLab - 0.28)*(pLab - 0.28) + 0.009);
        fNucleonElasticXsc = 0.257/((pLab - 0.28)*(pLab - 0.28) + 0.01);
       }
      else if( pLab < 0.5 )
      {
        hpXsc       = 26 + 110*(std::log(pLab/0.48))*(std::log(pLab/0.48));
        fNucleonElasticXsc = 0.37*hpXsc;
      }
      else if( pLab < 0.65 )
      {
        hpXsc       = 26 + 110*(std::log(pLab/0.48))*(std::log(pLab/0.48));
        fNucleonElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 0.72 )
      {
        hpXsc = 36.1+
                10*std::exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*std::exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fNucleonElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 0.88 )
      {
        hpXsc = 36.1+
                10*std::exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*std::exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fNucleonElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
      }
      else if( pLab < 1.03 )
      {
        hpXsc = 36.1+
                10*std::exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*std::exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fNucleonElasticXsc = 2.0 + 0.4/((pLab - 1.03)*(pLab - 1.03) + 0.016);
      }
      else if( pLab < 1.15 )
      {
        hpXsc = 36.1+
                10*std::exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*std::exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fNucleonElasticXsc = 2.0 + 0.4/((pLab - 1.03)*(pLab - 1.03) + 0.016);
      }
      else if( pLab < 1.3 )
      {
        hpXsc = 36.1+
                10*std::exp(-(pLab-0.72)*(pLab-0.72)/0.06/0.06)+
                24*std::exp(-(pLab-1.015)*(pLab-1.015)/0.075/0.075);
        fNucleonElasticXsc = 3. + 13./pLab;
      }
      else if( pLab < 2.6 ) // < 3.0) // ns original
      {
        hpXsc = 36.1+0.079-4.313*std::log(pLab)+
                3*std::exp(-(pLab-2.1)*(pLab-2.1)/0.4/0.4)+
                1.5*std::exp(-(pLab-1.4)*(pLab-1.4)/0.12/0.12);
        fNucleonElasticXsc = 3. +13./pLab; // *std::log(pLab*6.79);
      }
      else   // mb
      {
        hpXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
        fNucleonElasticXsc = 3. + 13./pLab;
      }
      fNucleonTotalXsc = hpXsc;
    }
  } 
  else if( (particlePDG ==  -321 || particlePDG == 310 ) && proton )   // Kmp/K0p //////
  {
    if( pLab < pMin)
    {
      double psp = pLab*std::sqrt(pLab);
      fNucleonElasticXsc  = 5.2/psp;
      fNucleonTotalXsc    = 14./psp;
    }
    else if( pLab > pMax )
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      fNucleonElasticXsc           = cofLogE*ld2 + 2.23;
      fNucleonTotalXsc           = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      double sp  = std::sqrt(pLab);
      double psp = pLab*sp;
      double p2  = pLab*pLab;
      double p4  = p2*p2;
      double lm  = pLab - .39;
      double md  = lm*lm + .000356;

      double lh1  = pLab - 0.78;
      double hd1  = lh1*lh1 + .00166;
 
      double lh  = pLab - 1.01;
      double hd  = lh*lh + .011;

      double lh2  = pLab - 1.63;
      double hd2  = lh2*lh2 + .007;

      fNucleonElasticXsc  = 5.2/psp + (1.1*cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) 
	             + .004/md + 0.005/hd1+ 0.01/hd2 +.15/hd; // small peaks were added

      fNucleonTotalXsc    = 14./psp + (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4) 
	             + .006/md  + 0.01/hd1+ 0.02/hd2 + .20/hd ;
    }
  }
  else if( ( particlePDG ==  -321 || particlePDG == 310 ) && neutron )   // Kmn/K0n ///////
  {
    if( pLab > pMax )
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      fNucleonElasticXsc           = cofLogE*ld2 + 2.23;
      fNucleonTotalXsc           = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    {
 
      double lh  = pLab - 0.98;
      double hd  = lh*lh + .021;

      double LogPlab = std::log( pLab );
      double sqrLogPlab = LogPlab * LogPlab;

      fNucleonElasticXsc  = // 5.2/psp + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) + .004/md 
                     5.0 +  8.1*std::pow(pLab,-1.8 ) + 0.16*sqrLogPlab - 1.3*LogPlab + .15/hd;
      fNucleonTotalXsc    = // 14./psp + 
                     //  (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4) 
	// WP                     25.2 +  0. *std::pow(pLab, 0.  ) + 0.38*sqrLogPlab - 2.9*LogPlab	             
                     25.2 +  0.38*sqrLogPlab - 2.9*LogPlab	             
                     //       + .006/md  + 0.01/hd1+ 0.02/hd2 
                        + 0.30/hd ;
    }
  }
  else if(  (particlePDG ==  321 || particlePDG == 130) && proton  )  // Kpp/aKp ////////////////////////
  {
    if( pLab < pMin )
    {
      double lr = pLab - .38;
      double lm = pLab - 1.;
      double md = lm*lm + .392;   
      fNucleonElasticXsc = .7/(lr*lr + .076) + 2./md;
      fNucleonTotalXsc   = .7/(lr*lr + .076) + 2.6/md;
    }
    else if( pLab > pMax )
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      fNucleonElasticXsc           = cofLogE*ld2 + 2.23;
      fNucleonTotalXsc           = cofLogT*ld2 + 19.2;
    }
    else
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      double lr  = pLab - .38;
      double LE  = .7/(lr*lr + .076);
      double sp  = std::sqrt(pLab);
      double p2  = pLab*pLab;
      double p4  = p2*p2;
      double lm  = pLab - 1.;
      double md  = lm*lm + .392;
      fNucleonElasticXsc  = LE + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
      fNucleonTotalXsc    = LE + (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 2.6/md;
    }
  }
  else if(  (particlePDG ==  321 || particlePDG == 130) && neutron  )  // Kpn/aKn ///////////////////////
  {
    if( pLab < pMin )
    {
      double lm = pLab - 0.94;
      double md = lm*lm + .392;   
      fNucleonElasticXsc = 2./md;
      fNucleonTotalXsc   = 4.6/md;
    }
    else if( pLab > pMax )
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      fNucleonElasticXsc           = cofLogE*ld2 + 2.23;
      fNucleonTotalXsc           = cofLogT*ld2 + 19.2;
    }
    else
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      double sp  = std::sqrt(pLab);
      double p2  = pLab*pLab;
      double p4  = p2*p2;
      double lm  = pLab - 0.94;
      double md  = lm*lm + .392;
      fNucleonElasticXsc  = (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
      fNucleonTotalXsc    = (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 4.6/md;
    }
  }
  else if( particlePDG == 3112 && pORn ) 
  {
    xsection  = 35.20 + B*std::pow(std::log(sMand/s0),2.) 
                          - 199.*std::pow(sMand,-eta1) + 264.*std::pow(sMand,-eta2);
  } 
  else if( particlePDG == 22  && pORn ) // modify later on
  {
    xsection  = 0.0 + B*std::pow(std::log(sMand/s0),2.) 
      + 0.032*std::pow(sMand,-eta1); // WP - 0.0*std::pow(sMand,-eta2);
    fNucleonTotalXsc = xsection;   
  } 
  else  // other then p,n,pi+,pi-,K+,K- as proton ??? 
  {
    if( proton )
    {
      xsection  = 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2);
    }
    if( neutron )
    {
      xsection += 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2);
    }
    fNucleonTotalXsc = xsection;
  } 
  fNucleonTotalXsc   *= geant::millibarn; // parametrised in mb
  fNucleonElasticXsc *= geant::millibarn; // parametrised in mb

  if( proton && geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGCharge() > 0. )
  {
    double proton_mass = geant::kProtonMassC2;
    double cB = GetCoulombBarrier(particlePDG, mass, energyKin, targetPDG, proton_mass);
    fNucleonTotalXsc   *= cB;
    fNucleonElasticXsc *= cB; 
  }
  fNucleonInelasticXsc = fNucleonTotalXsc - fNucleonElasticXsc;
  if( fNucleonInelasticXsc < 0. ) fNucleonInelasticXsc = 0.;

 
  return fNucleonTotalXsc;
}

/////////////////////////////////////////////////////////////////////////////////////
//
// Returns hadron-nucleon Xsc according to PDG parametrisation (2005):
// http://pdg.lbl.gov/2006/reviews/hadronicrpp.pdf
  
double 
HadronNucleusXsc::GetHadronNucleonXscPDG(int particlePDG, double mass, double energyKin, int targetPDG)
{
  double xsection(0);
  int Zt=1, Nt=1, At=1;

  double targ_mass = 0.939 * geant::GeV;  // ~mean neutron and proton ???

  double proj_mass     = mass;
  double proj_momentum = std::sqrt(energyKin*(energyKin + 2*proj_mass));
 
  double sMand = CalcMandelstamS ( proj_mass , targ_mass , proj_momentum );

  sMand         /= geant::GeV * geant::GeV;  // in GeV for parametrisation

  // General PDG fit constants

  double s0   = 5.38*5.38; // in Gev^2
  double eta1 = 0.458;
  double eta2 = 0.458;
  double B    = 0.308;

  // proton = 2212, neutron = 2112
  bool pORn = (targetPDG == 2212 || targetPDG == 2112);  
  bool proton = (targetPDG == 2212);
  bool neutron = (targetPDG == 2112);

  
  if(particlePDG == 2112) // proton-neutron fit 
  {
    if ( proton )
    {
      xsection = Zt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
		 + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));// on p
    }
    if ( neutron )
    {
      xsection  = Nt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
		      + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2)); // on n pp for nn
    }
  } 
  else if(particlePDG == 2212) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
    }
  } 
  else if(particlePDG == -2212) // particlePDG == theAProton) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                          + 42.53*std::pow(sMand,-eta1) + 33.34*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                          + 40.15*std::pow(sMand,-eta1) + 30.*std::pow(sMand,-eta2));
    }
  } 
  else if( particlePDG == 211 && pORn ) 
  {
    xsection  = At*( 20.86 + B*std::pow(std::log(sMand/s0),2.) 
                          + 19.24*std::pow(sMand,-eta1) - 6.03*std::pow(sMand,-eta2));
  } 
  else if(particlePDG == -211 && pORn ) 
  {
    xsection  = At*( 20.86 + B*std::pow(std::log(sMand/s0),2.) 
                          + 19.24*std::pow(sMand,-eta1) + 6.03*std::pow(sMand,-eta2));
  } 
  else if(particlePDG == 321) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 17.91 + B*std::pow(std::log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) - 13.45*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 17.87 + B*std::pow(std::log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) - 7.23*std::pow(sMand,-eta2));
    }
  } 
  else if(particlePDG == -321) 
  {
    if ( proton )
    {      
      xsection  = Zt*( 17.91 + B*std::pow(std::log(sMand/s0),2.) 
                          + 7.14*std::pow(sMand,-eta1) + 13.45*std::pow(sMand,-eta2));
    }
    if ( neutron )
    {
      xsection = Nt*( 17.87 + B*std::pow(std::log(sMand/s0),2.) 
                          + 5.17*std::pow(sMand,-eta1) + 7.23*std::pow(sMand,-eta2) );
    }
  }
  else if(particlePDG == 3112 && pORn ) 
  {
    xsection  = At*( 35.20 + B*std::pow(std::log(sMand/s0),2.) 
                          - 199.*std::pow(sMand,-eta1) + 264.*std::pow(sMand,-eta2) );
  } 
  else if(particlePDG == 22 && pORn ) // modify later on
  {
    xsection  = At*( 0.0 + B*std::pow(std::log(sMand/s0),2.) 
                          + 0.032*std::pow(sMand,-eta1) - 0.0*std::pow(sMand,-eta2) );
   
  } 
  else  // as proton ??? 
  {
    if ( proton )
    {      
      xsection  = Zt*( 35.45 + B*std::pow(std::log(sMand/s0),2.) 
                       + 42.53*std::pow(sMand,-eta1) - 33.34*std::pow(sMand,-eta2) );
    }
    if ( neutron )
    {
      xsection = Nt*( 35.80 + B*std::pow(std::log(sMand/s0),2.) 
                      + 40.15*std::pow(sMand,-eta1) - 30.*std::pow(sMand,-eta2));
    }
  } 
  xsection *= geant::millibarn; // parametrised in mb

  fNucleonTotalXsc     = xsection;
  fNucleonInelasticXsc = 0.75*xsection;
  fNucleonElasticXsc   = fNucleonTotalXsc - fNucleonInelasticXsc;
  if (fNucleonElasticXsc < 0.) fNucleonElasticXsc = 0.;

  return xsection;
}

/////////////////////////////////////////////////////////////////////////////////////
//
// Returns kaon-nucleon cross-section based on smoothed NS for GG model

double HadronNucleusXsc::GetKaonNucleonXscGG(int particlePDG, double mass, double energyKin, int targetPDG)
{

  double pLab = std::sqrt(energyKin*(energyKin + 2*mass));
  
  pLab /= geant::GeV;
  double LogPlab = std::log( pLab );
  double sqrLogPlab = LogPlab * LogPlab;

  double minLogP = 3.5;       // min of (lnP-minLogP)^2 
  double cofLogE = .0557;     // elastic (lnP-minLogP)^2 
  double cofLogT = .3;        // total (lnP-minLogP)^2 
  double pMin = .1;        // fast LE calculation 
  double pMax = 1000.;     // fast HE calculation 


  bool proton = (targetPDG == 2212);
  bool neutron = (targetPDG == 2112);

  if(  (particlePDG == -321 || particlePDG == 310) && proton ) // (K-,K0)on p ////////////////////////////
  {

    if( pLab < pMin)
    {
      double psp = pLab*std::sqrt(pLab);
      fNucleonElasticXsc  = 5.2/psp;
      fNucleonTotalXsc    = 14./psp;
    }
    else if( pLab > pMax )
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      fNucleonElasticXsc           = cofLogE*ld2 + 2.23;
      fNucleonTotalXsc           = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      double sp  = std::sqrt(pLab);
      double psp = pLab*sp;
      double p2  = pLab*pLab;
      double p4  = p2*p2;
 
      double lh  = pLab - 0.98;
      double hd  = lh*lh + .045;


      fNucleonElasticXsc  = 5.2/psp + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) // + .004/md 
               + .15/hd;
      fNucleonTotalXsc    = 14./psp + (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4) 
	              //  + .006/md  + 0.01/hd1 + 0.02/hd2 
                     + .60/hd;
    }
  }
  else if( (particlePDG == -321 || particlePDG == 310) && neutron )   // Kmn/K0n /////////////////////////////
  {
    if( pLab > pMax )
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      fNucleonElasticXsc           = cofLogE*ld2 + 2.23;
      fNucleonTotalXsc           = 1.1*cofLogT*ld2 + 19.7;
    }
    else
    {
 
      double lh  = pLab - 0.98;
      double hd  = lh*lh + .045;

      fNucleonElasticXsc  = // 5.2/psp + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) + .004/md 
                     5.0 +  8.1*std::pow(pLab,-1.8 ) + 0.16*sqrLogPlab - 1.3*LogPlab + .15/hd;
      fNucleonTotalXsc    = // 14./psp + 
                     //  (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4) 
	// WP                     25.2 +  0. *std::pow(pLab, 0.  ) + 0.38*sqrLogPlab - 2.9*LogPlab	             
                     25.2 + 0.38*sqrLogPlab - 2.9*LogPlab	             
                     //       + .006/md  + 0.01/hd1+ 0.02/hd2 
                        + 0.60/hd ;
    }
  }
  else if(  (particlePDG == 321 || particlePDG == 130) && proton )  // Kpp/aKp //////////////////////
  {
    if( pLab < pMin )
    {
      double lr = pLab - .38;
      double lm = pLab - 1.;
      double md = lm*lm + .392;   
      fNucleonElasticXsc = .7/(lr*lr + .076) + 2./md;
      fNucleonTotalXsc   = // .7/(lr*lr + .076) + 
                2.6/md;
    }
    else if( pLab > pMax )
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      fNucleonElasticXsc           = cofLogE*ld2 + 2.23;
      fNucleonTotalXsc           = cofLogT*ld2 + 19.2;
    }
    else
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      double lr  = pLab - .38;
      double LE  = .7/(lr*lr + .076);
      double sp  = std::sqrt(pLab);
      double p2  = pLab*pLab;
      double p4  = p2*p2;
      double lm  = pLab - 0.8;
      double md  = lm*lm + .652;
      fNucleonElasticXsc  = LE + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
      fNucleonTotalXsc    = (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 7.6/md; // + LE;
    }
  }
  else if( (particlePDG == 321 || particlePDG == 130) && neutron )  // Kpn/aKn //////////////////////////////////
  {
    if( pLab < pMin )
    {
      double lm = pLab - 0.94;
      double md = lm*lm + .392;   
      fNucleonElasticXsc = 2./md;
      fNucleonTotalXsc   = 4.6/md;
    }
    else if( pLab > pMax )
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      fNucleonElasticXsc           = cofLogE*ld2 + 2.23;
      fNucleonTotalXsc           = cofLogT*ld2 + 19.2;
    }
    else
    {
      double ld  = std::log(pLab) - minLogP;
      double ld2 = ld*ld;
      double sp  = std::sqrt(pLab);
      double p2  = pLab*pLab;
      double p4  = p2*p2;
      double lm  = pLab - 0.8;
      double md  = lm*lm + .652;
      fNucleonElasticXsc  = (cofLogE*ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
      fNucleonTotalXsc    = (cofLogT*ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 7.6/md;
    }
  }
  fNucleonTotalXsc   *= geant::millibarn; // parametrised in mb
  fNucleonElasticXsc *= geant::millibarn; // parametrised in mb

  if( proton && geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGCharge() > 0. )
  {
    double proton_mass = geant::kProtonMassC2;
    double cB = GetCoulombBarrier(particlePDG, mass, energyKin, targetPDG, proton_mass);
    fNucleonTotalXsc   *= cB;
    fNucleonElasticXsc *= cB; 
  }
  fNucleonInelasticXsc = fNucleonTotalXsc - fNucleonElasticXsc;
  if( fNucleonInelasticXsc < 0. ) fNucleonInelasticXsc = 0.;

  return fNucleonTotalXsc;
}




double HadronNucleusXsc::CalcMandelstamS( const double mp , 
					    const double mt , 
					    const double Plab )
{
  double Elab = std::sqrt ( mp * mp + Plab * Plab );
  double sMand  = mp*mp + mt*mt + 2*Elab*mt ;

  return sMand;
}


double HadronNucleusXsc::GetCoulombBarrier(int particlePDG, double proj_mass, double energyKin,
					   int targetPDG, double target_mass)
{
  double ratio;

  double tR = 0.895*geant::fermi, pR;

  if     ( particlePDG == 2212 ) pR = 0.895*geant::fermi;
  else if( particlePDG == 211 )  pR = 0.663*geant::fermi;
  else if( particlePDG == 321 )  pR = 0.340*geant::fermi;
  else                           pR = 0.500*geant::fermi;

  double pZ = geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGCharge();
  double tZ = geantphysics::Particle::GetParticleByPDGCode(targetPDG)->GetPDGCharge(); 

  double pTkin = energyKin;
  
  double pM    = proj_mass;
  double tM    = target_mass; 

  double pElab = pTkin + pM;

  double totEcm  = std::sqrt(pM*pM + tM*tM + 2.*pElab*tM);

  double totTcm  = totEcm - pM -tM;

  double bC    = geant::kFineStructConst * geant::kHBarPlanckCLight *pZ*tZ;
           bC   /= pR + tR;
           bC   /= 2.;  // 4., 2. parametrisation cof ??? vmg

  if( totTcm <= bC ) ratio = 0.;
  else               ratio = 1. - bC/totTcm;

  if( ratio < 0.) ratio = 0.;

  return ratio;
}



double HadronNucleusXsc::GetNucleusRadius(int At)
{
  double oneThird = 1.0/3.0;
  double cubicrAt = std::pow(double(At), oneThird); 
  const double RadiusConst = 1.08 * geant::fermi;
  
  double R;

  R = RadiusConst*cubicrAt;

  double meanA = 20.;
  double tauA  = 20.; 

  if (At > 20) 
  {
    R *= ( 0.8 + 0.2*std::exp( -(double(At) - meanA)/tauA) ); 
  }
  else
  {
    R *= ( 1.0 + 0.1*( 1. - std::exp( (double(At) - meanA)/tauA) ) ); 
  }

  return R;
}


///////////////////////////////////////////////////////////////////////////////
//
// Correction arrays for GG <-> Bar changea at ~ 90 GeV

const double HadronNucleusXsc::fNeutronBarCorrectionTot[93] = {

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

const double HadronNucleusXsc::fNeutronBarCorrectionIn[93] = {

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

const double HadronNucleusXsc::fProtonBarCorrectionTot[93] = {

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

const double HadronNucleusXsc::fProtonBarCorrectionIn[93] = {

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


const double HadronNucleusXsc::fPionPlusBarCorrectionTot[93] = {

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

const double HadronNucleusXsc::fPionPlusBarCorrectionIn[93] = {

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


const double HadronNucleusXsc::fPionMinusBarCorrectionTot[93] = {

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


const double HadronNucleusXsc::fPionMinusBarCorrectionIn[93] = {

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

