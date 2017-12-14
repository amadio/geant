
#include "Parameterizations.h"

#include "SystemOfUnits.h"
#include "PhysicalConstants.h"
#include "Proton.h"
#include "Neutron.h"

#include <cmath>


namespace geantphysics {


  /// Mandelstam
  double CalcMandelstamS(const double mp, const double mt, const double Plab) {
    double Elab   = std::sqrt ( mp * mp + Plab * Plab );
    double sMand  = mp * mp + mt * mt + 2 * Elab * mt ;
    return sMand;
  }

  /// Nucleus radius
  double GetNucleusRadius(int At) {
    constexpr double oneThird    = 1.0 / 3.0;
    constexpr double RadiusConst = 1.08 * geant::fermi;
    double cubicrAt = std::pow(double(At), oneThird);
    double R        = RadiusConst * cubicrAt;
    double meanA    = 20.;
    double tauA     = 20.;
    if (At > 20) {
      R *= ( 0.8 + 0.2 * std::exp( -(double(At) - meanA) / tauA) );
    } else {
      R *= ( 1.0 + 0.1 * ( 1. - std::exp( (double(At) - meanA) / tauA) ) );
    }
    return R;
  }


  /// Coulomb barrier
  // M. Novak: this method can be improved
  double GetCoulombBarrier(int particlePDG, double proj_mass, double energyKin, int targetPDG, double target_mass) {
    double ratio;
    double tR = 0.895 * geant::fermi, pR;

    if     ( particlePDG == 2212 ) pR = 0.895 * geant::fermi;
    else if( particlePDG == 211 )  pR = 0.663 * geant::fermi;
    else if( particlePDG == 321 )  pR = 0.340 * geant::fermi;
    else                           pR = 0.500 * geant::fermi;

    double pZ = geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGCharge();
    double tZ = geantphysics::Particle::GetParticleByPDGCode(targetPDG)->GetPDGCharge();

    double pTkin = energyKin;

    double pM      = proj_mass;
    double tM      = target_mass;

    double pElab   = pTkin + pM;

    double totEcm  = std::sqrt(pM * pM + tM * tM + 2. * pElab * tM);

    double totTcm  = totEcm - pM -tM;

    double bC      = geant::kFineStructConst * geant::kHBarPlanckCLight * pZ * tZ;
    
    bC   /= pR + tR;
    bC   /= 2.;  // 4., 2. parametrisation cof ??? vmg

    if( totTcm <= bC ) ratio = 0.;
    else               ratio = 1. - bC / totTcm;

    if( ratio < 0.) ratio = 0.;

    return ratio;
  }

  
  // PDG paremeterization
  double GetHadronNucleonXscPDG(int particlePDG, double mass, double energyKin, int targetPDG) {
    
    constexpr double targ_mass = 0.939 * geant::GeV;  // ~mean neutron and proton ???
    // General PDG fit constants
    constexpr double s0        = 5.38 * 5.38; // in Gev^2
    constexpr double eta1      = 0.458;
    constexpr double eta2      = 0.458;
    constexpr double B         = 0.308;
    // M.Novak: target Z,N and A properties will be set later when the method will be used
    int   Zt = 1;
    int   Nt = 1;
    int   At = 1;

    double NucleonTotalXsc = 0.;
    double proj_mass       = mass;
    double proj_momentum   = std::sqrt(energyKin * (energyKin + 2 * proj_mass));

    double sMand  = CalcMandelstamS(proj_mass, targ_mass, proj_momentum);
    sMand        /= geant::GeV * geant::GeV;  // in GeV for parametrisation

    double dumy0  = std::log(sMand/s0);
    double term1  = B * dumy0 * dumy0;
    double term2  = std::pow(sMand,-eta1);
    double term3  = std::pow(sMand,-eta2);

    // proton = 2212, neutron = 2112
    bool pORn    = (targetPDG == 2212 || targetPDG == 2112);
    bool proton  = (targetPDG == 2212);
    bool neutron = (targetPDG == 2112);

    if(particlePDG == 2112) { // proton-neutron fit
      if ( proton )  { NucleonTotalXsc  = Zt*( 35.80 + term1 + 40.15*term2 - 30.  *term3); } // on p
      if ( neutron ) { NucleonTotalXsc  = Nt*( 35.45 + term1 + 42.53*term2 - 33.34*term3); } // on n pp for nn
    } else if(particlePDG == 2212) {
      if ( proton )  { NucleonTotalXsc  = Zt*( 35.45 + term1 + 42.53*term2 - 33.34*term3); }
      if ( neutron ) { NucleonTotalXsc  = Nt*( 35.80 + term1 + 40.15*term2 - 30.  *term3); }
    } else if(particlePDG == -2212) { // particlePDG == theAProton)
      if ( proton )  { NucleonTotalXsc  = Zt*( 35.45 + term1 + 42.53*term2 + 33.34*term3); }
      if ( neutron ) { NucleonTotalXsc  = Nt*( 35.80 + term1 + 40.15*term2 + 30.  *term3); }
    } else  if(particlePDG == 321) {
      if ( proton )  { NucleonTotalXsc  = Zt*( 17.91 + term1 +  7.14*term2 - 13.45*term3); }
      if ( neutron ) { NucleonTotalXsc  = Nt*( 17.87 + term1 +  5.17*term2 -  7.23*term3); }
    } else if(particlePDG == -321) {
      if ( proton )  { NucleonTotalXsc  = Zt*( 17.91 + term1 + 7.14*term2 + 13.45*term3); }
      if ( neutron ) { NucleonTotalXsc  = Nt*( 17.87 + term1 + 5.17*term2 +  7.23*term3); }
    } else if( particlePDG == 211 && pORn ) {
      NucleonTotalXsc  = At*( 20.86 + term1 +  19.24*term2 -   6.03*term3);
    } else if(particlePDG == -211 && pORn ) {
      NucleonTotalXsc  = At*( 20.86 + term1 +  19.24*term2 +   6.03*term3);
    } else if(particlePDG == 3112 && pORn ) {
      NucleonTotalXsc  = At*( 35.20 + term1 - 199.  *term2 + 264.  *term3);
    } else if(particlePDG == 22 && pORn ) { // modify later on
      // M. Novak: a bit too many zeros in this line :)
      NucleonTotalXsc  = At*( 0.0 + term1 + 0.032*term2 - 0.0*term3 );
    } else { // as proton ???
      if ( proton )  { NucleonTotalXsc  = Zt*( 35.45 + term1 + 42.53*term2 - 33.34*term3); }
      if ( neutron ) { NucleonTotalXsc  = Nt*( 35.80 + term1 + 40.15*term2 - 30.  *term3); }
    }
    NucleonTotalXsc *= geant::millibarn; // parametrised in mb

    return NucleonTotalXsc;
  }


  // Starkov parameterization
  double GetHadronNucleonTotalXscNS(int particlePDG, double mass, double energyKin, int targetPDG) {

    double NucleonTotalXsc(0);

    double A0, B0;

    double tM = 0.939 * geant::GeV;  // ~mean neutron and proton ???

    double pM   = mass;
    double pE   = energyKin + mass; // total energy!!!!
    double pLab  = std::sqrt(energyKin * (energyKin + 2. * pM));

    double sMand = CalcMandelstamS ( pM , tM , pLab );
    double logP = std::log(pLab);

    // General PDG fit constants
    double s0   = 5.38 * 5.38; // in Gev^2
    double eta1 = 0.458;
    double eta2 = 0.458;
    double B    = 0.308;

    // proton = 2212, neutron = 2112
    bool pORn = (targetPDG == 2212 || targetPDG == 2112);
    bool proton = (targetPDG == 2212);
    bool neutron = (targetPDG == 2112);

    if( particlePDG == 2112 && pORn )
      {
	if( pLab >= 373.)
	  {
	    NucleonTotalXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
	  }
	else if( pLab >= 100.)
	  {
	    B0 = 7.5;
	    A0 = 100. - B0*std::log(3.0e7);
	    
	    NucleonTotalXsc = A0 + B0*std::log(pE) - 11 + 103*std::pow(sMand,-0.165); // mb
	  }
	else if( pLab >= 10.)
	  {
	    B0 = 7.5;
	    A0 = 100. - B0*std::log(3.0e7);

	    NucleonTotalXsc = A0 + B0*std::log(pE) - 11 + 103*std::pow(2*0.93827*pE + pM*pM +0.93827*0.93827,-0.165); //  mb
	  }
	else  // pLab < 10 GeV/c
	  {
	    if( neutron )      // nn to be pp
	      {
		if( pLab < 0.4 )
		  {
		    NucleonTotalXsc = 23 + 50 * ( std::pow( std::log(0.73/pLab), 3.5 ) );
		  }
		else if( pLab < 0.73 )
		  {
		    NucleonTotalXsc = 23 + 50 * ( std::pow( std::log(0.73/pLab), 3.5 ) );
		  }
		else if( pLab < 1.05  )
		  {
		    double lp73 = std::log(pLab/0.73);
		    NucleonTotalXsc = 23 + 40 * lp73 * lp73;
		  }
		else    // 1.05 - 10 GeV/c
		  {
		    NucleonTotalXsc = 39.0 + 75 * (pLab - 1.2)/(std::pow(pLab,3.0) + 0.15);
		  }
	      }
	    if( proton )   // pn to be np
	      {
		if( pLab < 0.02 )
		  {
		    NucleonTotalXsc = 4100 + 30 * std::pow(std::log(1.3/pLab), 3.6); // was as pLab < 0.8
		  }
		else if( pLab < 0.8 )
		  {
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/1.3), 4.0);
		  }
		else if( pLab < 1.05 )
		  {
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/0.95), 2.0);
		  }
		else if( pLab < 1.4 )
		  {
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/0.95), 2.0);
		  }
		else    // 1.4 < pLab < 10.  )
		  {
		    NucleonTotalXsc = 33.3 + 20.8 * (std::pow(pLab,2.0) - 1.35)/(std::pow(pLab,2.50) + 0.95);
		  }
	      }
	  }
      }
    else if( particlePDG == 2212 && pORn ) ////// proton //////////////////////////////////////////////
      {
	if( pLab >= 373.) // pdg due to TOTEM data
	  {
	    NucleonTotalXsc =  GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
	  }
	else if( pLab >= 100.)
	  {
	    B0 = 7.5;
	    A0 = 100. - B0*std::log(3.0e7);

	    NucleonTotalXsc = A0 + B0 * std::log(pE) - 11 + 103 * std::pow(sMand,-0.165);        //  mb
	  }
	else if( pLab >= 10.)
	  {
	    B0 = 7.5;
	    A0 = 100. - B0 * std::log(3.0e7);

	    NucleonTotalXsc = A0 + B0 * std::log(pE) - 11 + 103 * std::pow(sMand,-0.165);        //  mb
    	  }
	else
	  {
	    // pp

	    if( proton )
	      {
		if( pLab < 0.4 )
		  {
		    NucleonTotalXsc = 23 + 50 * ( std::pow( std::log(0.73/pLab), 3.5 ) );
		  }
		else if( pLab < 0.73 )
		  {
		    NucleonTotalXsc = 23 + 50 * ( std::pow( std::log(0.73/pLab), 3.5 ) );
		  }
		else if( pLab < 1.05  )
		  {
		    double lp73 = std::log(pLab/0.73);
		    NucleonTotalXsc = 23 + 40 * lp73 * lp73;
		  }
		else    // 1.05 - 10 GeV/c
		  {
		    NucleonTotalXsc = 39.0 + 75 * (pLab - 1.2)/(std::pow(pLab,3.0) + 0.15);
		  }
	      }
	    if( neutron )     // pn to be np
	      {
		if( pLab < 0.02 )
		  {
		    NucleonTotalXsc = 4100 + 30 * std::pow(std::log(1.3/pLab),3.6); // was as pLab < 0.8
		  }
		else if( pLab < 0.8 )
		  {
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/1.3),4.0);
		  } 
		else if( pLab < 1.05 )
		  {
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/0.95),2.0);
		  }
		else if( pLab < 1.4 )
		  {
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/0.95),2.0);
		  }
		else    // 1.4 < pLab < 10.  )
		  {
		    NucleonTotalXsc = 33.3 + 20.8 * (std::pow(pLab,2.0) - 1.35)/(std::pow(pLab,2.50) + 0.95);
		  }
	      }
	  }
      }
    else if( particlePDG == -2212 && pORn ) /////////////////// p_bar ///////////////////////////
      {
	if( proton )
	  {
	    NucleonTotalXsc = 35.45 + B * std::pow(std::log(sMand/s0),2.) + 42.53 * std::pow(sMand,-eta1) + 33.34 * std::pow(sMand,-eta2);
	  }
	if( neutron ) // ???
	  {
	    NucleonTotalXsc = 35.80 + B * std::pow(std::log(sMand/s0),2.) + 40.15 * std::pow(sMand,-eta1) + 30. * std::pow(sMand,-eta2);
	  }
      }
    else if( particlePDG == 211 && pORn ) // pi+ /////////////////////////////////////////////
      {
	if( proton ) // pi+ p
	  {
	    if( pLab < 0.28 )
	      {
		double lp1273 = logP + 1.273;
		NucleonTotalXsc = 10./(lp1273 * lp1273 + 0.05);
	      }
	    else if( pLab < 0.4 )
	      {
		double lp1273 = logP + 1.273;
		NucleonTotalXsc = 14./(lp1273 * lp1273 + 0.07);
	      }
	    else if( pLab < 0.68 )
	      {
		double lp1273 = logP + 1.273;
		NucleonTotalXsc = 14./(lp1273 * lp1273 + 0.07);
	      }
	    else if( pLab < 1.15 )
	      {
		double lp77 = std::log(pLab/0.77);
		NucleonTotalXsc = 88 * lp77 * lp77 + 14.9;
	      }
	    else if( pLab < 1.4) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
	      }
	    else if( pLab < 2.0 ) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
	      }
	    else if( pLab < 3.5 ) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
	      }
	    else if( pLab < 200. ) // my
	      {
		NucleonTotalXsc = 10.6 + 2.*std::log(pE) + 25 * std::pow(pE, -0.43 ); // ns original
	      }
	    else //  pLab > 100 // my
	      {
		NucleonTotalXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
	      }
	  }
	if( neutron )  // pi+ n = pi- p??
	  {
	    if( pLab < 0.28 )
	      {
		NucleonTotalXsc       = 0.288/((pLab - 0.28) * (pLab - 0.28) + 0.004);
	      }
	    else if( pLab < 0.395676 ) // first peak
	      {
		NucleonTotalXsc       = 0.648/((pLab - 0.28) * (pLab - 0.28) + 0.009);
	      }
	    else if( pLab < 0.5 )
	      {
		NucleonTotalXsc       = 26 + 110 * (std::log(pLab/0.48)) * (std::log(pLab/0.48));
	      }
	    else if( pLab < 0.65 )
	      {
		NucleonTotalXsc       = 26 + 110 * (std::log(pLab/0.48)) * (std::log(pLab/0.48));
	      }
	    else if( pLab < 0.72 )
	      {
		NucleonTotalXsc = 36.1 + 10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
	      }
	    else if( pLab < 0.88 )
	      {
		NucleonTotalXsc = 36.1 + 10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
	      }
	    else if( pLab < 1.03 )
	      {
		NucleonTotalXsc = 36.1 + 10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
	      }
	    else if( pLab < 1.15 )
	      {
		NucleonTotalXsc = 36.1 + 10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
	      }
	    else if( pLab < 1.3 )
	      {
		NucleonTotalXsc = 36.1 + 10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
	      }
	    else if( pLab < 2.6 ) // < 3.0) // ns original
	      {
		NucleonTotalXsc = 36.1 + 0.079-4.313 * std::log(pLab) +
		  3 * std::exp(-(pLab-2.1) * (pLab-2.1)/0.4/0.4) +
		  1.5 * std::exp(-(pLab-1.4) * (pLab-1.4)/0.12/0.12);
	      }
	    else if( pLab < 20. ) // < 3.0) // ns original
	      {
		NucleonTotalXsc = 36.1 + 0.079 - 4.313 * std::log(pLab) +
		  3 * std::exp(-(pLab-2.1) * (pLab-2.1)/0.4/0.4) +
		  1.5 * std::exp(-(pLab-1.4) * (pLab-1.4)/0.12/0.12);
	      }
	    else   // mb
	      {
		NucleonTotalXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
	      }
	  }
      }
    else if( particlePDG == -211 && pORn ) /// pi- ////////////////////////////////////////////
      {
	if( neutron )     // pi- n = pi+ p??
	  {
	    if( pLab < 0.28 )
	      {
		NucleonTotalXsc = 10./((logP + 1.273) * (logP + 1.273) + 0.05);
	      }
	    else if( pLab < 0.4 )
	      {
		NucleonTotalXsc = 14./( (logP + 1.273) * (logP + 1.273) + 0.07);
	      }
	    else if( pLab < 0.68 )
	      {
		NucleonTotalXsc = 14./( (logP + 1.273) * (logP + 1.273) + 0.07);
	      }
	    else if( pLab < 0.85 )
	      {
		double lp77 = std::log(pLab/0.77);
		NucleonTotalXsc = 88 * lp77 * lp77 + 14.9;
	      }
	    else if( pLab < 1.15 )
	      {
		double lp77 = std::log(pLab/0.77);
		NucleonTotalXsc = 88 * lp77 * lp77 + 14.9;
	      }
	    else if( pLab < 1.4) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
	      }
	    else if( pLab < 2.0 ) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc        = Ex1 + Ex2 + 27.5;
	      }
	    else if( pLab < 3.5 ) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
	      }
	    else if( pLab < 200. ) // my
	      {
		NucleonTotalXsc = 10.6 + 2.*std::log(pE) + 25*std::pow(pE, -0.43 ); // ns original
	      }
	    else //  pLab > 100 // my
	      {
		NucleonTotalXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
	      }
	  }
	if( proton )    // pi- p
	  {
	    if( pLab < 0.28 )
	      {
		NucleonTotalXsc = 0.288/((pLab - 0.28) * (pLab - 0.28) + 0.004);
	      }
	    else if( pLab < 0.395676 ) // first peak
	      {
		NucleonTotalXsc = 0.648/((pLab - 0.28) * (pLab - 0.28) + 0.009);
	      }
	    else if( pLab < 0.5 )
	      {
		double lp48 = std::log(pLab/0.48);
		NucleonTotalXsc = 26 + 110 * lp48 * lp48;
	      }
	    else if( pLab < 0.65 )
	      {
		double lp48 = std::log(pLab/0.48);
		NucleonTotalXsc = 26 + 110 * lp48 * lp48;
	      }
	    else if( pLab < 0.72 )
	      {
		NucleonTotalXsc = 36.1 +
		  10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
	      }
	    else if( pLab < 0.88 )
	      {
		NucleonTotalXsc = 36.1 +
		  10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
	      }
	    else if( pLab < 1.03 )
	      {
		NucleonTotalXsc = 36.1 +
		  10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
	      }
	    else if( pLab < 1.15 )
	      {
		NucleonTotalXsc = 36.1 +
		  10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
	      }
	    else if( pLab < 1.3 )
	      {
		NucleonTotalXsc = 36.1 +
		  10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06) +
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
	      }
	    else if( pLab < 2.6 ) // < 3.0) // ns original
	      {
		NucleonTotalXsc = 36.1 + 0.079 - 4.313 * std::log(pLab) +
		  3 * std::exp(-(pLab-2.1) * (pLab-2.1)/0.4/0.4) +
		  1.5 * std::exp(-(pLab-1.4) * (pLab-1.4)/0.12/0.12);
	      }
	    else   // mb
	      {
		NucleonTotalXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
	      }
	  }
      }
    else if( particlePDG == 3112 && pORn )
      {
	NucleonTotalXsc = 35.20 + B * std::pow(std::log(sMand/s0),2.)
	  - 199. * std::pow(sMand,-eta1) + 264. * std::pow(sMand,-eta2);
      }
    else if( particlePDG == 22  && pORn ) // modify later on
      {
	NucleonTotalXsc = 0.0 + B * std::pow(std::log(sMand/s0),2.)
	  + 0.032 * std::pow(sMand,-eta1); // WP - 0.0*std::pow(sMand,-eta2);
      }
    else  // other then p,n,pi+,pi-,K+,K- as proton ???
      {
	if( proton )
	  {
	    NucleonTotalXsc = 35.45 + B * std::pow(std::log(sMand/s0),2.)
	      + 42.53 * std::pow(sMand,-eta1) - 33.34 * std::pow(sMand,-eta2);
	  }
	if( neutron )
	  {
	    NucleonTotalXsc += 35.80 + B * std::pow(std::log(sMand/s0),2.)
	      + 40.15 * std::pow(sMand,-eta1) - 30. * std::pow(sMand,-eta2);
	  }
      }
    NucleonTotalXsc   *= geant::millibarn; // parametrised in mb

    if( proton && geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGCharge() > 0. )
      {
	double proton_mass = geant::kProtonMassC2;
	double cB = GetCoulombBarrier(particlePDG, mass, energyKin, targetPDG, proton_mass);
	NucleonTotalXsc   *= cB;
      }

    return NucleonTotalXsc;
  }

  // M. Novak: this method requires significant cleanup
  double GetKaonNucleonTotalXscGG(int particlePDG, double mass, double energyKin, int targetPDG) {
    double pLab = std::sqrt(energyKin*(energyKin + 2*mass));

    pLab /= geant::GeV;
    double LogPlab = std::log( pLab );
    double sqrLogPlab = LogPlab * LogPlab;

    double minLogP = 3.5;       // min of (lnP-minLogP)^2
    //double cofLogE = .0557;     // elastic (lnP-minLogP)^2
    double cofLogT = .3;        // total (lnP-minLogP)^2
    double pMin = .1;        // fast LE calculation
    double pMax = 1000.;     // fast HE calculation

    double NucleonTotalXsc(0);

    bool proton = (targetPDG == 2212);
    bool neutron = (targetPDG == 2112);

    if(  (particlePDG == -321 || particlePDG == 310) && proton ) // (K-,K0)on p ////////////////////////////
      {

	if( pLab < pMin)
	  {
	    double psp = pLab * std::sqrt(pLab);
	    NucleonTotalXsc = 14./psp;
	  }
	else if( pLab > pMax )
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonTotalXsc = 1.1 * cofLogT * ld2 + 19.7;
	  }
	else
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    double sp  = std::sqrt(pLab);
	    double psp = pLab * sp;
	    double p2  = pLab * pLab;
	    double p4  = p2 * p2;

	    double lh  = pLab - 0.98;
	    double hd  = lh * lh + .045;

	    NucleonTotalXsc = 14./psp + (1.1 * cofLogT * ld2 + 19.5)/(1. - .21/sp + .52/p4)
	      //  + .006/md  + 0.01/hd1 + 0.02/hd2
	      + .60/hd;
	  }
      }
    else if( (particlePDG == -321 || particlePDG == 310) && neutron )   // Kmn/K0n /////////////////////////////
      {
	if( pLab > pMax )
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonTotalXsc = 1.1 * cofLogT * ld2 + 19.7;
	  }
	else
	  {

	    double lh  = pLab - 0.98;
	    double hd  = lh * lh + .045;

	    NucleonTotalXsc    = // 14./psp +
	      //  (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4)
	      // WP                     25.2 +  0. *std::pow(pLab, 0.  ) + 0.38*sqrLogPlab - 2.9*LogPlab
	      25.2 + 0.38 * sqrLogPlab - 2.9 * LogPlab
	      //       + .006/md  + 0.01/hd1+ 0.02/hd2
	      + 0.60/hd ;
	  }
      }
    else if(  (particlePDG == 321 || particlePDG == 130) && proton )  // Kpp/aKp //////////////////////
      {
	if( pLab < pMin )
	  {
	    //double lr = pLab - .38;
	    double lm = pLab - 1.;
	    double md = lm*lm + .392;
	    NucleonTotalXsc   = // .7/(lr*lr + .076) +
	      2.6/md;
	  }
	else if( pLab > pMax )
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonTotalXsc = cofLogT * ld2 + 19.2;
	  }
	else
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    //double lr  = pLab - .38;
	    //double LE  = .7/(lr*lr + .076);
	    double sp  = std::sqrt(pLab);
	    double p2  = pLab * pLab;
	    double p4  = p2 * p2;
	    double lm  = pLab - 0.8;
	    double md  = lm * lm + .652;
	    NucleonTotalXsc = (cofLogT * ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 7.6/md; // + LE;
	  }
      }
    else if( (particlePDG == 321 || particlePDG == 130) && neutron )  // Kpn/aKn //////////////////////////////////
      {
	if( pLab < pMin )
	  {
	    double lm = pLab - 0.94;
	    double md = lm * lm + .392;
	    NucleonTotalXsc = 4.6/md;
	  }
	else if( pLab > pMax )
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonTotalXsc = cofLogT * ld2 + 19.2;
	  }
	else
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    double sp  = std::sqrt(pLab);
	    double p2  = pLab * pLab;
	    double p4  = p2 * p2;
	    double lm  = pLab - 0.8;
	    double md  = lm * lm + .652;
	    NucleonTotalXsc = (cofLogT * ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 7.6/md;
	  }
      }
    NucleonTotalXsc   *= geant::millibarn; // parametrised in mb

    if( proton && geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGCharge() > 0. )
      {
	double proton_mass = geant::kProtonMassC2;
	double cB = GetCoulombBarrier(particlePDG, mass, energyKin, targetPDG, proton_mass);
	NucleonTotalXsc   *= cB;
      }

    return NucleonTotalXsc;
  }

  // 
  double GetHadronNucleonInelasticXscNS(int particlePDG, double mass, double energyKin, int targetPDG) {
    double NucleonInelasticXsc(0), NucleonElasticXsc(0), NucleonTotalXsc(0);

    double A0, B0;

    double tM = 0.939 * geant::GeV;  // ~mean neutron and proton ???

    double pM = mass;
    double pE = energyKin + mass; // total energy!!!!
    double pLab = std::sqrt(energyKin * (energyKin + 2. * pM));

    double sMand = CalcMandelstamS ( pM , tM , pLab );

    double logP = std::log(pLab);

    // General PDG fit constants

    double s0   = 5.38 * 5.38; // in Gev^2
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
	    NucleonTotalXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;

	    NucleonElasticXsc = 6.5 + 0.308 * std::pow(std::log(sMand/400.),1.65) + 9.19 * std::pow(sMand,-0.458);
	  }
	else if( pLab >= 100.)
	  {
	    B0 = 7.5;
	    A0 = 100. - B0 * std::log(3.0e7);

	    NucleonTotalXsc = A0 + B0 * std::log(pE) - 11
	      // + 103*std::pow(2*0.93827*pE + pM*pM+0.93827*0.93827,-0.165);        //  mb
	      + 103 * std::pow(sMand,-0.165);        //  mb

	    NucleonElasticXsc = 5.53 + 0.308*std::pow(std::log(sMand/28.9),1.1) + 9.19*std::pow(sMand,-0.458);
	  }
	else if( pLab >= 10.)
	  {
	    B0 = 7.5;
	    A0 = 100. - B0 * std::log(3.0e7);

	    NucleonTotalXsc = A0 + B0 * std::log(pE) - 11 + 103 * std::pow(2 * 0.93827 * pE + pM * pM + 0.93827 * 0.93827,-0.165); //  mb
	    NucleonElasticXsc =  6 + 20/( (logP-0.182) * (logP-0.182) + 1.0 );
	  }
	else  // pLab < 10 GeV/c
	  {
	    if( neutron )      // nn to be pp
	      {
		if( pLab < 0.4 )
		  {
		    NucleonTotalXsc = 23 + 50 * ( std::pow( std::log(0.73/pLab), 3.5 ) );
		    NucleonElasticXsc = NucleonTotalXsc;
		  }
		else if( pLab < 0.73 )
		  {
		    NucleonTotalXsc = 23 + 50 * ( std::pow( std::log(0.73/pLab), 3.5 ) );
		    NucleonElasticXsc = NucleonTotalXsc;
		  }
		else if( pLab < 1.05  )
		  {
		    double lp73 = std::log(pLab/0.73);
		    NucleonTotalXsc = 23 + 40 * lp73 * lp73;
		    NucleonElasticXsc = 23 + 20 * lp73 * lp73;
		  }
		else    // 1.05 - 10 GeV/c
		  {
		    NucleonTotalXsc = 39.0 + 75 * (pLab - 1.2)/(std::pow(pLab,3.0) + 0.15);

		    NucleonElasticXsc = 6 + 20/((logP-0.182) * (logP-0.182) + 1.0);
		  }
	      }
	    if( proton )   // pn to be np
	      {
		if( pLab < 0.02 )
		  {
		    NucleonTotalXsc = 4100 + 30 * std::pow(std::log(1.3/pLab),3.6); // was as pLab < 0.8
		    NucleonElasticXsc = NucleonTotalXsc;
		  }
		else if( pLab < 0.8 )
		  {
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/1.3),4.0);
		    NucleonElasticXsc = NucleonTotalXsc;
		  }
		else if( pLab < 1.05 )
		  {
		    double l5p = std::log(0.511/pLab);
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/0.95),2.0);
		    NucleonElasticXsc = 6 + 52/(l5p * l5p + 1.6);
		  }
		else if( pLab < 1.4 )
		  {
		    double l5p = std::log(0.511/pLab);
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/0.95),2.0);
		    NucleonElasticXsc = 6 + 52/(l5p * l5p + 1.6 );
		  }
		else    // 1.4 < pLab < 10.  )
		  {
		    double lp0 = (logP-0.182);
		    NucleonTotalXsc = 33.3 + 20.8 * (std::pow(pLab,2.0) - 1.35)/(std::pow(pLab,2.50) + 0.95);
		    NucleonElasticXsc = 6 + 20/(lp0 * lp0 + 1.0);
		  }
	      }
	  }
      }
    else if( particlePDG == 2212 && pORn ) ////// proton //////////////////////////////////////////////
      {
	if( pLab >= 373.) // pdg due to TOTEM data
	  {
	    NucleonTotalXsc =  GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;

	    NucleonElasticXsc = 6.5 + 0.308 * std::pow(std::log(sMand/400.),1.65) + 9.19 * std::pow(sMand,-0.458);
	  }
	else if( pLab >= 100.)
	  {
	    B0 = 7.5;
	    A0 = 100. - B0 * std::log(3.0e7);

	    NucleonTotalXsc = A0 + B0 * std::log(pE) - 11 + 103 * std::pow(sMand,-0.165);        //  mb

	    NucleonElasticXsc = 5.53 + 0.308 * std::pow(std::log(sMand/28.9),1.1) + 9.19 * std::pow(sMand,-0.458);
	  }
	else if( pLab >= 10.)
	  {
	    B0 = 7.5;
	    A0 = 100. - B0 * std::log(3.0e7);

	    NucleonTotalXsc = A0 + B0 * std::log(pE) - 11 + 103 * std::pow(sMand,-0.165);        //  mb

	    NucleonElasticXsc =  6 + 20/( (logP-0.182) * (logP-0.182) + 1.0 );
	  }
	else
	  {
	    // pp

	    if( proton )
	      {
		if( pLab < 0.4 )
		  {
		    NucleonTotalXsc = 23 + 50 * ( std::pow( std::log(0.73/pLab), 3.5 ) );
		    NucleonElasticXsc = NucleonTotalXsc;
		  }
		else if( pLab < 0.73 )
		  {
		    NucleonTotalXsc = 23 + 50 * ( std::pow( std::log(0.73/pLab), 3.5 ) );
		    NucleonElasticXsc = NucleonTotalXsc;
		  }
		else if( pLab < 1.05  )
		  {
		    double lp73 = std::log(pLab/0.73);
		    NucleonTotalXsc = 23 + 40 * lp73 * lp73;
		    NucleonElasticXsc = 23 + 20 * lp73 * lp73;
		  }
		else    // 1.05 - 10 GeV/c
		  {
		    NucleonTotalXsc = 39.0 + 75 * (pLab - 1.2)/(std::pow(pLab,3.0) + 0.15);
		    NucleonElasticXsc = 6 + 20/( (logP-0.182) * (logP-0.182) + 1.0 );
		  }
	      }
	    if( neutron )     // pn to be np
	      {
		if( pLab < 0.02 )
		  {
		    NucleonTotalXsc = 4100 + 30 * std::pow(std::log(1.3/pLab),3.6); // was as pLab < 0.8
		    NucleonElasticXsc = NucleonTotalXsc;
		  }
		else if( pLab < 0.8 )
		  {
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/1.3),4.0);
		    NucleonElasticXsc = NucleonTotalXsc;
		  }
		else if( pLab < 1.05 )
		  {
		    double l5p = std::log(0.511/pLab);
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/0.95),2.0);
		    NucleonElasticXsc =  6 + 52/( l5p * l5p + 1.6 );
		  }
		else if( pLab < 1.4 )
		  {
		    double l5p = std::log(0.511/pLab);
		    NucleonTotalXsc = 33 + 30 * std::pow(std::log(pLab/0.95),2.0);
		    NucleonElasticXsc =  6 + 52/(l5p * l5p + 1.6 );
		  }
		else    // 1.4 < pLab < 10.  )
		  {
		    NucleonTotalXsc = 33.3 + 20.8*(std::pow(pLab,2.0) - 1.35)/(std::pow(pLab,2.50) + 0.95);

		    NucleonElasticXsc =  6 + 20/( (logP-0.182)*(logP-0.182) + 1.0 );
		  }
	      }
	  }
      }
    else if( particlePDG == -2212 && pORn ) /////////////////// p_bar ///////////////////////////
      {
	if( proton )
	  {
	    NucleonTotalXsc  = 35.45 + B * std::pow(std::log(sMand/s0),2.)
	      + 42.53 * std::pow(sMand,-eta1) + 33.34 * std::pow(sMand,-eta2);
	  }
	if( neutron ) // ???
	  {
	    NucleonTotalXsc = 35.80 + B * std::pow(std::log(sMand/s0),2.)
	      + 40.15 * std::pow(sMand,-eta1) + 30. * std::pow(sMand,-eta2);
	  }
      }
    else if( particlePDG == 211 && pORn ) // pi+ /////////////////////////////////////////////
      {
	if( proton ) // pi+ p
	  {
	    if( pLab < 0.28 )
	      {
		NucleonTotalXsc = 10./((logP + 1.273) * (logP + 1.273) + 0.05);
		NucleonElasticXsc = NucleonTotalXsc;
	      }
	    else if( pLab < 0.4 )
	      {
		NucleonTotalXsc = 14./( (logP + 1.273) * (logP + 1.273) + 0.07);
		NucleonElasticXsc = NucleonTotalXsc;
	      }
	    else if( pLab < 0.68 )
	      {
		NucleonTotalXsc = 14./( (logP + 1.273) * (logP + 1.273) + 0.07);
		NucleonElasticXsc = NucleonTotalXsc;
	      }
	    else if( pLab < 0.85 )
	      {
		double lp77 = std::log(pLab/0.77);
		
		NucleonTotalXsc = 88 * lp77 * lp77 + 14.9;
		NucleonElasticXsc = NucleonTotalXsc * std::exp(-3.*(pLab - 0.68));
	      }
	    else if( pLab < 1.15 )
	      {
		double lp77 = std::log(pLab/0.77);

		NucleonTotalXsc = 88 * lp77 * lp77 + 14.9;
		NucleonElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
	      }
	    else if( pLab < 1.4) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
		NucleonElasticXsc = 6.0 + 1.4/(( pLab - 1.4) * ( pLab - 1.4) + 0.1);
	      }
	    else if( pLab < 2.0 ) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
		NucleonElasticXsc = 3.0 + 1.36/( (logP - 0.336) * (logP - 0.336) + 0.08);
	      }
	    else if( pLab < 3.5 ) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
		NucleonElasticXsc = 3.0 + 6.20/((logP - 0.336) * (logP - 0.336) + 0.8);
	      }
	    else if( pLab < 200. ) // my
	      {
		NucleonTotalXsc = 10.6 + 2. * std::log(pE) + 25 * std::pow(pE, -0.43 ); // ns original
		NucleonElasticXsc = 3.0 + 6.20/((logP - 0.336) * (logP - 0.336) + 0.8);
	      }
	    else //  pLab > 100 // my
	      {
		NucleonTotalXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
		NucleonElasticXsc = 3.0 + 6.20/((logP - 0.336) * (logP - 0.336) + 0.8);
	      }
	  }
	if( neutron )  // pi+ n = pi- p??
	  {
	    if( pLab < 0.28 )
	      {
		NucleonTotalXsc = 0.288/((pLab - 0.28) * (pLab - 0.28) + 0.004);
		NucleonElasticXsc = 1.8/((logP + 1.273) * (logP + 1.273) + 0.07);
	      }
	    else if( pLab < 0.395676 ) // first peak
	      {
		NucleonTotalXsc = 0.648/((pLab - 0.28) * (pLab - 0.28) + 0.009);
		NucleonElasticXsc = 0.257/((pLab - 0.28) * (pLab - 0.28) + 0.01);
	      }
	    else if( pLab < 0.5 )
	      {
		double lp48 = std::log(pLab/0.48);
		NucleonTotalXsc = 26 + 110 * lp48 * lp48;
		NucleonElasticXsc = 0.37 * NucleonTotalXsc;
	      }
	    else if( pLab < 0.65 )
	      {
		double lp48 = std::log(pLab/0.48);
		NucleonTotalXsc = 26 + 110 * lp48 * lp48;
		NucleonElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
	      }
	    else if( pLab < 0.72 )
	      {
		NucleonTotalXsc = 36.1 + 10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06) +
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
		NucleonElasticXsc = 0.95/((pLab - 0.72) * (pLab - 0.72) + 0.049);
	      }
	    else if( pLab < 0.88 )
	      {
		NucleonTotalXsc = 36.1 + 10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
		NucleonElasticXsc = 0.95/((pLab - 0.72) * (pLab - 0.72) + 0.049);
	      }
	    else if( pLab < 1.03 )
	      {
		NucleonTotalXsc = 36.1 + 10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06)+
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
		NucleonElasticXsc = 2.0 + 0.4/((pLab - 1.03) * (pLab - 1.03) + 0.016);
	      }
	    else if( pLab < 1.15 )
	      {
		NucleonTotalXsc = 36.1 + 10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06) +
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
		NucleonElasticXsc = 2.0 + 0.4/((pLab - 1.03) * (pLab - 1.03) + 0.016);
	      }
	    else if( pLab < 1.3 )
	      {
		NucleonTotalXsc = 36.1 + 10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06) +
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
		NucleonElasticXsc = 3. + 13./pLab;
	      }
	    else if( pLab < 2.6 ) // < 3.0) // ns original
	      {
		NucleonTotalXsc = 36.1 + 0.079-4.313 * std::log(pLab) + 3 * std::exp(-(pLab-2.1) * (pLab-2.1)/0.4/0.4) +
		  1.5 * std::exp(-(pLab-1.4) * (pLab-1.4)/0.12/0.12);
		NucleonElasticXsc = 3. + 13./pLab;
	      }
	    else if( pLab < 20. ) // < 3.0) // ns original
	      {
		NucleonTotalXsc = 36.1 + 0.079 - 4.313 * std::log(pLab) + 3 * std::exp(-(pLab-2.1) * (pLab-2.1)/0.4/0.4) +
		  1.5 * std::exp(-(pLab-1.4) * (pLab-1.4)/0.12/0.12);
		NucleonElasticXsc = 3. + 13./pLab;
	      }
	    else   // mb
	      {
		NucleonTotalXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
		NucleonElasticXsc = 3. + 13./pLab;
	      }
	  }
      }
    else if( particlePDG == -211 && pORn ) /// pi- ////////////////////////////////////////////
      {
	if( neutron )     // pi- n = pi+ p??
	  {
	    if( pLab < 0.28 )
	      {
		NucleonTotalXsc = 10./((logP + 1.273) * (logP + 1.273) + 0.05);
		NucleonElasticXsc = NucleonTotalXsc;
	      }
	    else if( pLab < 0.4 )
	      {
		NucleonTotalXsc = 14./( (logP + 1.273) * (logP + 1.273) + 0.07);
		NucleonElasticXsc = NucleonTotalXsc;
	      }
	    else if( pLab < 0.68 )
	      {
		NucleonTotalXsc = 14./( (logP + 1.273) * (logP + 1.273) + 0.07);
		NucleonElasticXsc = NucleonTotalXsc;
	      }
	    else if( pLab < 0.85 )
	      {
		double lp77 = std::log(pLab/0.77);
		NucleonTotalXsc = 88 * lp77 * lp77 + 14.9;
		NucleonElasticXsc = NucleonTotalXsc*std::exp(-3. * (pLab - 0.68));
	      }
	    else if( pLab < 1.15 )
	      {
		double lp77 = std::log(pLab/0.77);
		NucleonTotalXsc = 88 * lp77 * lp77 + 14.9;

		NucleonElasticXsc = 6.0 + 1.4/(( pLab - 1.4)*( pLab - 1.4) + 0.1);
	      }
	    else if( pLab < 1.4) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
		NucleonElasticXsc = 6.0 + 1.4/(( pLab - 1.4) * ( pLab - 1.4) + 0.1);
	      }
	    else if( pLab < 2.0 ) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
		NucleonElasticXsc = 3.0 + 1.36/( (logP - 0.336) * (logP - 0.336) + 0.08);
	      }
	    else if( pLab < 3.5 ) // ns original
	      {
		double Ex1 = 3.2 * std::exp(-(pLab-2.55) * (pLab-2.55)/0.55/0.55);
		double Ex2 = 12 * std::exp(-(pLab-1.47) * (pLab-1.47)/0.225/0.225);
		NucleonTotalXsc = Ex1 + Ex2 + 27.5;
		NucleonElasticXsc = 3.0 + 6.20/((logP - 0.336) * (logP - 0.336) + 0.8);
	      }
	    else if( pLab < 200. ) // my
	      {
		NucleonTotalXsc = 10.6 + 2. * std::log(pE) + 25 * std::pow(pE, -0.43 ); // ns original
		NucleonElasticXsc = 3.0 + 6.20/( (logP - 0.336) * (logP - 0.336) + 0.8);
	      }
	    else //  pLab > 100 // my
	      {
		NucleonTotalXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
		NucleonElasticXsc = 3.0 + 6.20/( (logP - 0.336) * (logP - 0.336) + 0.8);
	      }
	  }
	if( proton )    // pi- p
	  {
	    if( pLab < 0.28 )
	      {
		NucleonTotalXsc = 0.288/((pLab - 0.28) * (pLab - 0.28) + 0.004);
		NucleonElasticXsc = 1.8/((logP + 1.273) * (logP + 1.273) + 0.07);
	      }
	    else if( pLab < 0.395676 ) // first peak
	      {
		NucleonTotalXsc = 0.648/((pLab - 0.28) * (pLab - 0.28) + 0.009);
		NucleonElasticXsc = 0.257/((pLab - 0.28) * (pLab - 0.28) + 0.01);
	      }
	    else if( pLab < 0.5 )
	      {
		double lp48 = std::log(pLab/0.48);
		NucleonTotalXsc = 26 + 110 * lp48 * lp48;
		NucleonElasticXsc = 0.37 * NucleonTotalXsc;
	      }
	    else if( pLab < 0.65 )
	      {
		double lp48 = std::log(pLab/0.48);
		NucleonTotalXsc = 26 + 110 * lp48 * lp48;
		NucleonElasticXsc = 0.95/((pLab - 0.72) * (pLab - 0.72) + 0.049);
	      }
	    else if( pLab < 0.72 )
	      {
		NucleonTotalXsc = 36.1 +
		  10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06) +
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
		NucleonElasticXsc = 0.95/((pLab - 0.72) * (pLab - 0.72) + 0.049);
	      }
	    else if( pLab < 0.88 )
	      {
		NucleonTotalXsc = 36.1 +
		  10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06) +
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
		NucleonElasticXsc = 0.95/((pLab - 0.72)*(pLab - 0.72) + 0.049);
	      }
	    else if( pLab < 1.03 )
	      {
		NucleonTotalXsc = 36.1 +
		  10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06) +
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
		NucleonElasticXsc = 2.0 + 0.4/((pLab - 1.03) * (pLab - 1.03) + 0.016);
	      }
	    else if( pLab < 1.15 )
	      {
		NucleonTotalXsc = 36.1 +
		  10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06) +
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
		NucleonElasticXsc = 2.0 + 0.4/((pLab - 1.03) * (pLab - 1.03) + 0.016);
	      }
	    else if( pLab < 1.3 )
	      {
		NucleonTotalXsc = 36.1 +
		  10 * std::exp(-(pLab-0.72) * (pLab-0.72)/0.06/0.06) +
		  24 * std::exp(-(pLab-1.015) * (pLab-1.015)/0.075/0.075);
		NucleonElasticXsc = 3. + 13./pLab;
	      }
	    else if( pLab < 2.6 ) // < 3.0) // ns original
	      {
		NucleonTotalXsc = 36.1 + 0.079 - 4.313 * std::log(pLab) +
		  3 * std::exp(-(pLab-2.1) * (pLab-2.1)/0.4/0.4) +
		  1.5 * std::exp(-(pLab-1.4) * (pLab-1.4)/0.12/0.12);
		NucleonElasticXsc = 3. +13./pLab; // *std::log(pLab*6.79);
	      }
	    else   // mb
	      {
		NucleonTotalXsc = GetHadronNucleonXscPDG(particlePDG, mass, energyKin, targetPDG)/geant::millibarn;
		NucleonElasticXsc = 3. + 13./pLab;
	      }
	  }
      }
    else if( (particlePDG == -321 || particlePDG == 310 ) && proton )   // Kmp/K0p //////
      {
	if( pLab < pMin)
	  {
	    double psp = pLab*std::sqrt(pLab);
	    NucleonElasticXsc = 5.2/psp;
	    NucleonTotalXsc = 14./psp;
	  }
	else if( pLab > pMax )
	  {
	    double ld = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonElasticXsc = cofLogE * ld2 + 2.23;
	    NucleonTotalXsc = 1.1 * cofLogT * ld2 + 19.7;
	  }
	else
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    double sp  = std::sqrt(pLab);
	    double psp = pLab * sp;
	    double p2  = pLab * pLab;
	    double p4  = p2 * p2;
	    double lm  = pLab - .39;
	    double md  = lm * lm + .000356;

	    double lh1  = pLab - 0.78;
	    double hd1  = lh1 * lh1 + .00166;

	    double lh  = pLab - 1.01;
	    double hd  = lh * lh + .011;

	    double lh2  = pLab - 1.63;
	    double hd2  = lh2 * lh2 + .007;

	    NucleonElasticXsc = 5.2/psp + (1.1 * cofLogE * ld2 + 2.23)/(1. - .7/sp + .075/p4)
	      + .004/md + 0.005/hd1+ 0.01/hd2 +.15/hd; // small peaks were added

	    NucleonTotalXsc = 14./psp + (1.1 * cofLogT * ld2 + 19.5)/(1. - .21/sp + .52/p4)
	      + .006/md  + 0.01/hd1+ 0.02/hd2 + .20/hd ;
	  }
      }
    else if( ( particlePDG ==  -321 || particlePDG == 310 ) && neutron )   // Kmn/K0n ///////
      {
	if( pLab > pMax )
	  {
	    double ld = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonElasticXsc = cofLogE * ld2 + 2.23;
	    NucleonTotalXsc = 1.1 * cofLogT * ld2 + 19.7;
	  }
	else
	  {

	    double lh  = pLab - 0.98;
	    double hd  = lh * lh + .021;

	    double LogPlab = std::log( pLab );
	    double sqrLogPlab = LogPlab * LogPlab;

	    NucleonElasticXsc  = // 5.2/psp + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) + .004/md
	      5.0 +  8.1 * std::pow(pLab,-1.8 ) + 0.16 * sqrLogPlab - 1.3 * LogPlab + .15/hd;
	    NucleonTotalXsc    = // 14./psp +
	      //  (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4)
	      // WP                     25.2 +  0. *std::pow(pLab, 0.  ) + 0.38*sqrLogPlab - 2.9*LogPlab
	      25.2 + 0.38*sqrLogPlab - 2.9*LogPlab
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
	    double md = lm * lm + .392;
	    NucleonElasticXsc = .7/(lr * lr + .076) + 2./md;
	    NucleonTotalXsc   = .7/(lr * lr + .076) + 2.6/md;
	  }
	else if( pLab > pMax )
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonElasticXsc = cofLogE * ld2 + 2.23;
	    NucleonTotalXsc = cofLogT * ld2 + 19.2;
	  }
	else
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld*ld;
	    double lr  = pLab - .38;
	    double LE  = .7/(lr*lr + .076);
	    double sp  = std::sqrt(pLab);
	    double p2  = pLab * pLab;
	    double p4  = p2 * p2;
	    double lm  = pLab - 1.;
	    double md  = lm * lm + .392;
	    NucleonElasticXsc  = LE + (cofLogE * ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
	    NucleonTotalXsc    = LE + (cofLogT * ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 2.6/md;
	  }
      }
    else if(  (particlePDG ==  321 || particlePDG == 130) && neutron  )  // Kpn/aKn ///////////////////////
      {
	if( pLab < pMin )
	  {
	    double lm = pLab - 0.94;
	    double md = lm * lm + .392;
	    NucleonElasticXsc = 2./md;
	    NucleonTotalXsc   = 4.6/md;
	  }
	else if( pLab > pMax )
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonElasticXsc = cofLogE * ld2 + 2.23;
	    NucleonTotalXsc = cofLogT * ld2 + 19.2;
	  }
	else
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    double sp  = std::sqrt(pLab);
	    double p2  = pLab * pLab;
	    double p4  = p2 * p2;
	    double lm  = pLab - 0.94;
	    double md  = lm * lm + .392;
	    NucleonElasticXsc  = (cofLogE * ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
	    NucleonTotalXsc    = (cofLogT * ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 4.6/md;
	  }
      }
    else if( particlePDG == 3112 && pORn )
      {
	NucleonTotalXsc  = 35.20 + B * std::pow(std::log(sMand/s0),2.)
	  - 199. * std::pow(sMand,-eta1) + 264. * std::pow(sMand,-eta2);
      }
    else if( particlePDG == 22  && pORn ) // modify later on
      {
	NucleonTotalXsc  = 0.0 + B * std::pow(std::log(sMand/s0),2.)
	  + 0.032 * std::pow(sMand,-eta1); // WP - 0.0*std::pow(sMand,-eta2);
      }
    else  // other then p,n,pi+,pi-,K+,K- as proton ???
      {
	if( proton )
	  {
	    NucleonTotalXsc = 35.45 + B * std::pow(std::log(sMand/s0),2.)
	      + 42.53 * std::pow(sMand,-eta1) - 33.34 * std::pow(sMand,-eta2);
	  }
	if( neutron )
	  {
	    NucleonTotalXsc += 35.80 + B * std::pow(std::log(sMand/s0),2.)
	      + 40.15 * std::pow(sMand,-eta1) - 30. * std::pow(sMand,-eta2);
	  }
      }
    NucleonTotalXsc   *= geant::millibarn; // parametrised in mb
    NucleonElasticXsc *= geant::millibarn; // parametrised in mb

    if( proton && geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGCharge() > 0. )
      {
	double proton_mass = geant::kProtonMassC2;
	double cB = GetCoulombBarrier(particlePDG, mass, energyKin, targetPDG, proton_mass);
	NucleonTotalXsc   *= cB;
	NucleonElasticXsc *= cB;
      }
    NucleonInelasticXsc = NucleonTotalXsc - NucleonElasticXsc;
    if( NucleonInelasticXsc < 0. ) NucleonInelasticXsc = 0.;

    return NucleonInelasticXsc;
  }


  // M. Novak: this method requires significant cleanup
  double GetKaonNucleonInelasticXscGG(int particlePDG, double mass, double energyKin, int targetPDG) {
    double NucleonElasticXsc(0), NucleonTotalXsc(0);

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
	    double psp = pLab * std::sqrt(pLab);
	    NucleonElasticXsc  = 5.2/psp;
	    NucleonTotalXsc    = 14./psp;
	  }
	else if( pLab > pMax )
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld*ld;
	    NucleonElasticXsc = cofLogE * ld2 + 2.23;
	    NucleonTotalXsc = 1.1 * cofLogT * ld2 + 19.7;
	  }
	else
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    double sp  = std::sqrt(pLab);
	    double psp = pLab * sp;
	    double p2  = pLab * pLab;
	    double p4  = p2 * p2;

	    double lh  = pLab - 0.98;
	    double hd  = lh * lh + .045;


	    NucleonElasticXsc  = 5.2/psp + (cofLogE * ld2 + 2.23)/(1. - .7/sp + .075/p4) // + .004/md
	      + .15/hd;
	    NucleonTotalXsc    = 14./psp + (1.1 * cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4)
	      //  + .006/md  + 0.01/hd1 + 0.02/hd2
	      + .60/hd;
	  }
      }
    else if( (particlePDG == -321 || particlePDG == 310) && neutron )   // Kmn/K0n /////////////////////////////
      {
	if( pLab > pMax )
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonElasticXsc = cofLogE * ld2 + 2.23;
	    NucleonTotalXsc = 1.1 * cofLogT * ld2 + 19.7;
	  }
	else
	  {

	    double lh  = pLab - 0.98;
	    double hd  = lh * lh + .045;

	    NucleonElasticXsc  = // 5.2/psp + (cofLogE*ld2 + 2.23)/(1. - .7/sp + .075/p4) + .004/md
	      5.0 +  8.1 * std::pow(pLab,-1.8 ) + 0.16 * sqrLogPlab - 1.3 * LogPlab + .15/hd;
	    NucleonTotalXsc    = // 14./psp +
	      //  (1.1*cofLogT*ld2 + 19.5)/(1. - .21/sp + .52/p4)
	      // WP                     25.2 +  0. *std::pow(pLab, 0.  ) + 0.38*sqrLogPlab - 2.9*LogPlab
	      25.2 + 0.38 * sqrLogPlab - 2.9 * LogPlab
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
	    double md = lm * lm + .392;
	    NucleonElasticXsc = .7/(lr * lr + .076) + 2./md;
	    NucleonTotalXsc   = // .7/(lr*lr + .076) +
	      2.6/md;
	  }
	else if( pLab > pMax )
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonElasticXsc = cofLogE * ld2 + 2.23;
	    NucleonTotalXsc = cofLogT * ld2 + 19.2;
	  }
	else
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    double lr  = pLab - .38;
	    double LE  = .7/(lr * lr + .076);
	    double sp  = std::sqrt(pLab);
	    double p2  = pLab * pLab;
	    double p4  = p2 * p2;
	    double lm  = pLab - 0.8;
	    double md  = lm * lm + .652;
	    NucleonElasticXsc = LE + (cofLogE * ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
	    NucleonTotalXsc = (cofLogT * ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 7.6/md; // + LE;
	  }
      }
    else if( (particlePDG == 321 || particlePDG == 130) && neutron )  // Kpn/aKn //////////////////////////////////
      {
	if( pLab < pMin )
	  {
	    double lm = pLab - 0.94;
	    double md = lm * lm + .392;
	    NucleonElasticXsc = 2./md;
	    NucleonTotalXsc   = 4.6/md;
	  }
	else if( pLab > pMax )
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    NucleonElasticXsc = cofLogE * ld2 + 2.23;
	    NucleonTotalXsc = cofLogT * ld2 + 19.2;
	  }
	else
	  {
	    double ld  = std::log(pLab) - minLogP;
	    double ld2 = ld * ld;
	    double sp  = std::sqrt(pLab);
	    double p2  = pLab * pLab;
	    double p4  = p2 * p2;
	    double lm  = pLab - 0.8;
	    double md  = lm * lm + .652;
	    NucleonElasticXsc  = (cofLogE * ld2 + 2.23)/(1. - .7/sp + .1/p4) + 2./md;
	    NucleonTotalXsc    = (cofLogT * ld2 + 19.5)/(1. + .46/sp + 1.6/p4) + 7.6/md;
	  }
      }
    NucleonTotalXsc   *= geant::millibarn; // parametrised in mb
    NucleonElasticXsc *= geant::millibarn; // parametrised in mb

    if( proton && geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGCharge() > 0. )
      {
	double proton_mass = geant::kProtonMassC2;
	double cB = GetCoulombBarrier(particlePDG, mass, energyKin, targetPDG, proton_mass);
	NucleonTotalXsc   *= cB;
	NucleonElasticXsc *= cB;
      }
    double NucleonInelasticXsc = NucleonTotalXsc - NucleonElasticXsc;
    if( NucleonInelasticXsc < 0. ) NucleonInelasticXsc = 0.;

    return NucleonInelasticXsc;
  }


}
