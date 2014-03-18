#include "GPUrbanMscModel95.h"
#include "GPPhysicalConstants.h"
#include "GPSystemOfUnits.h"
#include "GPRandom.h"

//#include "Randomize.hh"
//#include "G4Electron.hh"
//#include "G4LossTableManager.hh"
//#include "G4ParticleChangeForMSC.hh"
//#include "G4Poisson.hh"
//#include "globals.hh"
//using namespace std;

FQUALIFIER
GPUrbanMscModel95::GPUrbanMscModel95(curandState* devStates, 
				     int threadId)
{
  //curand
  fThreadId = threadId;
  fDevStates = devStates;
  theLambdaTable = 0;

  masslimite    = 0.6*MeV;
  lambdalimit   = 1.*mm;
  fr            = 0.02;
  taubig        = 8.0;
  tausmall      = 1.e-16;
  taulim        = 1.e-6;
  currentTau    = taulim;
  tlimitminfix  = 1.e-6*mm;            
  stepmin       = tlimitminfix;
  smallstep     = 1.e10;
  currentRange  = 0. ;
  rangeinit     = 0.;
  tlimit        = 1.e10*mm;
  tlimitmin     = 10.*tlimitminfix;            
  tgeom         = 1.e50*mm;
  geombig       = 1.e50*mm;
  geommin       = 1.e-3*mm;
  geomlimit     = geombig;
  presafety     = 0.*mm;
  //facsafety     = 0.50 ;
                          
  Zold          = 0.;
  Zeff          = 1.;
  Z2            = 1.;                
  Z23           = 1.;                    
  lnZ           = 0.;
  coeffth1      = 0.;
  coeffth2      = 0.;
  coeffc1       = 0.;
  coeffc2       = 0.;
  coeffc3       = 0.;
  coeffc4       = 0.;
  scr1ini       = fine_structure_const*fine_structure_const*
                  electron_mass_c2*electron_mass_c2/(0.885*0.885*4.*pi);
  scr2ini       = 3.76*fine_structure_const*fine_structure_const;
  scr1          = 0.;
  scr2          = 0.;

  firstStep     = true; 
  inside        = false;  
  insideskin    = false;

 //from base G4VMasModel
  skin = 1.0;
  dtrl = 0.05;
  facrange = 0.04;
  facsafety = 0.3;
  skindepth = skin*stepmin;
  samplez = false;
  latDisplasment = true;
  dedx       = 2.0*MeV*cm2/g;
  localrange = DBL_MAX;
  localtkin  = 0.0;

  mass = proton_mass_c2;
  charge = 1.0;
  currentKinEnergy = currentRadLength = lambda0 = lambdaeff = tPathLength 
    = zPathLength = par1 = par2 = par3 = 0;

  //  currentMaterialIndex = -1;
  //  fParticleChange = 0;
  //  couple = 0;
  SetSampleZ(true);

  //G4VEmModel
  pParticleChange = 0;
  lowLimit = 0.1*keV;
  highLimit = 100.0*TeV;
  eMinActive = 0.0;
  eMaxActive = DBL_MAX;

}

FQUALIFIER
GPUrbanMscModel95::~GPUrbanMscModel95()
{}

FQUALIFIER
void GPUrbanMscModel95::SetSampleZ(G4bool val) {
  samplez = val;
}


FQUALIFIER
void GPUrbanMscModel95::Initialise()
{
  skindepth = skin*stepmin;

  // set values of some data members
  SetParticle();
  //  fParticleChange = GetParticleChangeForMSC(p);

}

FQUALIFIER
void GPUrbanMscModel95::StartTracking()
{
  //  SetParticle(track->GetDynamicParticle()->GetDefinition());
  SetParticle();
  firstStep = true; 
  inside = false;
  insideskin = false;
  tlimit = geombig;
  stepmin = tlimitminfix ;
  tlimitmin = 10.*stepmin ;
}

FQUALIFIER
G4double GPUrbanMscModel95::ComputeTruePathLengthLimit2(
                             GPMaterial* material,
                             G4double kineticEnergy,
                             G4double xlambda,
			     G4double* currentMinimalStep)
{
  tPathLength = *currentMinimalStep;
  //  const G4DynamicParticle* dp = track.GetDynamicParticle();
  
  //  G4StepPoint* sp = track.GetStep()->GetPreStepPoint();
  //  G4StepStatus stepStatus = sp->GetStepStatus();
  //  couple = track.GetMaterialCutsCouple();
  //  SetCurrentCouple(couple); 
  //  currentMaterialIndex = couple->GetIndex();

  currentKinEnergy = kineticEnergy; // dp->GetKineticEnergy();
  //  currentRange = GetRange(particle,currentKinEnergy,couple);
  currentRange = GetRange(material,currentKinEnergy);
  //  lambda0 = GetTransportMeanFreePath(particle,currentKinEnergy);
  //  lambda0 = GetTransportMeanFreePath(currentKinEnergy);
  lambda0 = xlambda;

  if(tPathLength > currentRange) { tPathLength = currentRange; }

  // stop here if small range particle
  if(inside || tPathLength < tlimitminfix) { 
    //    return ConvertTrueToGeom(material,tPathLength, currentMinimalStep); 
    return ConvertTrueToGeom2(material,&tPathLength,currentMinimalStep);
  }
  
  //@@@  presafety = sp->GetSafety();
  presafety = 1.0*mm;

  // far from geometry boundary
  if(currentRange < presafety)
    {
      inside = true;
      //      return ConvertTrueToGeom(material,tPathLength, currentMinimalStep);  
      return ConvertTrueToGeom2(material,&tPathLength,currentMinimalStep);
    }


    {
      // compute presafety again if presafety <= 0 and no boundary
      // i.e. when it is needed for optimization purposes

      //@@@ implement SafetyHelper !!!
      //@@@      if((stepStatus != fGeomBoundary) && (presafety < tlimitminfix)) 
      //@@@        presafety = ComputeSafety(sp->GetPosition(),tPathLength); 
      //@@@ use temporary
      presafety = 1.0*mm;

      // is far from boundary
      if(currentRange < presafety)
        {
          inside = true;
	  //          return ConvertTrueToGeom(material,tPathLength, currentMinimalStep);    
	  return ConvertTrueToGeom2(material,&tPathLength,currentMinimalStep);
        }

      //@@@      if(firstStep || stepStatus == fGeomBoundary)
      if(firstStep)
      {
        rangeinit = currentRange;
        fr = facrange;
        // 9.1 like stepping for e+/e- only (not for muons,hadrons)
        if(mass < masslimite) 
        {
          if(lambda0 > currentRange)
            rangeinit = lambda0;
          if(lambda0 > lambdalimit)
            fr *= 0.75+0.25*lambda0/lambdalimit;
        }

        //lower limit for tlimit
        G4double rat = currentKinEnergy/MeV ;
        rat = 1.e-3/(rat*(10.+rat)) ;
        tlimitmin = 10.*lambda0*rat;
        if(tlimitmin < tlimitminfix) tlimitmin = tlimitminfix;
      }
      //step limit
      tlimit = fr*rangeinit;               

      if(tlimit < facsafety*presafety)
        tlimit = facsafety*presafety;

      //lower limit for tlimit
      if(tlimit < tlimitmin) tlimit = tlimitmin;
      
      if(tPathLength > tlimit) tPathLength = tlimit;

    }
  

    //    return ConvertTrueToGeom(material,tPathLength, currentMinimalStep);
    return ConvertTrueToGeom2(material,&tPathLength,currentMinimalStep);
}

//G4double G4VMscModel::ConvertTrueToGeom(G4double& tlength, 
FQUALIFIER
G4double GPUrbanMscModel95::ConvertTrueToGeom2(GPMaterial* material,
					       G4double* tlength,
					       G4double* glength)
{
  *glength = ComputeGeomPathLength2(material);
  // should return true length 
  return *tlength;
}

G4double GPUrbanMscModel95::ComputeGeomPathLength2(GPMaterial* material)
{
  firstStep = false; 
  lambdaeff = lambda0;
  par1 = -1. ;  
  par2 = par3 = 0. ;  

  //  do the true -> geom transformation
  zPathLength = tPathLength;

  // z = t for very small tPathLength
  if(tPathLength < tlimitminfix) return zPathLength;

  // this correction needed to run MSC with eIoni and eBrem inactivated
  // and makes no harm for a normal run
  // It is already checked
  // if(tPathLength > currentRange)
  //  tPathLength = currentRange ;

  G4double tau   = tPathLength/lambda0 ;

  if ((tau <= tausmall) || insideskin) {
    zPathLength  = tPathLength;
    if(zPathLength > lambda0) zPathLength = lambda0;
    return zPathLength;
  }

  G4double zmean = tPathLength;
  if (tPathLength < currentRange*dtrl) {
    if(tau < taulim) zmean = tPathLength*(1.-0.5*tau) ;
    else             zmean = lambda0*(1.-exp(-tau));
    zPathLength = zmean ;
    return zPathLength;    

  } else if(currentKinEnergy < mass || tPathLength == currentRange)  {
    par1 = 1./currentRange ;
    par2 = 1./(par1*lambda0) ;
    par3 = 1.+par2 ;
    if(tPathLength < currentRange)
      zmean = (1.-exp(par3*log(1.-tPathLength/currentRange)))/(par1*par3) ;
    else {
      zmean = 1./(par1*par3) ;
    }
    zPathLength = zmean ;
    return zPathLength;    

  } else {
    //    G4double T1 = GetEnergy(particle,currentRange-tPathLength,couple);
    G4double T1 = GetEnergy(material,currentRange-tPathLength);
    //    G4double lambda1 = GetTransportMeanFreePath(particle,T1);
    G4double lambda1 = GetTransportMeanFreePath(T1);

    par1 = (lambda0-lambda1)/(lambda0*tPathLength);
    par2 = 1./(par1*lambda0);
    par3 = 1.+par2 ;
    zmean = (1.-exp(par3*log(lambda1/lambda0)))/(par1*par3);
  }

  zPathLength = zmean;

  //  sample z
  if(samplez)
  {
    const G4double  ztmax = 0.999 ;
    G4double zt = zmean/tPathLength ;

    if (tPathLength > stepmin && zt < ztmax)              
    {
      G4double u,cz1;
      if(zt >= 1./3. ) { //third = 1./3.
        G4double cz = 0.5*(3.*zt-1.)/(1.-zt) ;
        cz1 = 1.+cz ;
        G4double u0 = cz/cz1 ;
        G4double grej ;
        do {
            u = exp(log(GPUniformRand(fDevStates,fThreadId))/cz1) ;
            grej = exp(cz*log(u/u0))*(1.-u)/(1.-u0) ;
           } while (grej < GPUniformRand(fDevStates,fThreadId)) ;
      }
      else
      {
        u = 2.*zt*GPUniformRand(fDevStates,fThreadId);
      }
      zPathLength = tPathLength*u ;
    }
  }

  if(zPathLength > lambda0) { zPathLength = lambda0; }
  return zPathLength;
}

G4double GPUrbanMscModel95::ComputeGeomPathLength(GPMaterial* material)
{
  firstStep = false; 
  lambdaeff = lambda0;
  par1 = -1. ;  
  par2 = par3 = 0. ;  

  //  do the true -> geom transformation
  zPathLength = tPathLength;

  // z = t for very small tPathLength
  if(tPathLength < tlimitminfix) return zPathLength;

  // this correction needed to run MSC with eIoni and eBrem inactivated
  // and makes no harm for a normal run
  // It is already checked
  // if(tPathLength > currentRange)
  //  tPathLength = currentRange ;

  G4double tau   = tPathLength/lambda0 ;

  if ((tau <= tausmall) || insideskin) {
    zPathLength  = tPathLength;
    if(zPathLength > lambda0) zPathLength = lambda0;
    return zPathLength;
  }

  G4double zmean = tPathLength;
  if (tPathLength < currentRange*dtrl) {
    if(tau < taulim) zmean = tPathLength*(1.-0.5*tau) ;
    else             zmean = lambda0*(1.-exp(-tau));
    zPathLength = zmean ;
    return zPathLength;    

  } else if(currentKinEnergy < mass || tPathLength == currentRange)  {
    par1 = 1./currentRange ;
    par2 = 1./(par1*lambda0) ;
    par3 = 1.+par2 ;
    if(tPathLength < currentRange)
      zmean = (1.-exp(par3*log(1.-tPathLength/currentRange)))/(par1*par3) ;
    else {
      zmean = 1./(par1*par3) ;
    }
    zPathLength = zmean ;
    return zPathLength;    

  } else {
    //    G4double T1 = GetEnergy(particle,currentRange-tPathLength,couple);
    G4double T1 = GetEnergy(material,currentRange-tPathLength);
    //    G4double lambda1 = GetTransportMeanFreePath(particle,T1);
    G4double lambda1 = GetTransportMeanFreePath(T1);

    par1 = (lambda0-lambda1)/(lambda0*tPathLength);
    par2 = 1./(par1*lambda0);
    par3 = 1.+par2 ;
    zmean = (1.-exp(par3*log(lambda1/lambda0)))/(par1*par3);
  }

  zPathLength = zmean;

  //  sample z
  if(samplez)
  {
    const G4double  ztmax = 0.999 ;
    G4double zt = zmean/tPathLength ;

    if (tPathLength > stepmin && zt < ztmax)              
    {
      G4double u,cz1;
      if(zt >= 1./3.) { //  third = 1./3.
        G4double cz = 0.5*(3.*zt-1.)/(1.-zt) ;
        cz1 = 1.+cz ;
        G4double u0 = cz/cz1 ;
        G4double grej ;
        do {
            u = exp(log(GPUniformRand(fDevStates,fThreadId))/cz1) ;
            grej = exp(cz*log(u/u0))*(1.-u)/(1.-u0) ;
           } while (grej < GPUniformRand(fDevStates,fThreadId)) ;
      }
      else
      {
        u = 2.*zt*GPUniformRand(fDevStates,fThreadId);
      }
      zPathLength = tPathLength*u ;
    }
  }

  if(zPathLength > lambda0) { zPathLength = lambda0; }
  return zPathLength;
}

FQUALIFIER
G4double GPUrbanMscModel95::ComputeTruePathLengthLimit(
                             GPMaterial* material,
                             G4double kineticEnergy,
                             G4double xlambda,
			     G4double& currentMinimalStep)
{
  tPathLength = currentMinimalStep;
  //  const G4DynamicParticle* dp = track.GetDynamicParticle();
  
  //  G4StepPoint* sp = track.GetStep()->GetPreStepPoint();
  //  G4StepStatus stepStatus = sp->GetStepStatus();
  //  couple = track.GetMaterialCutsCouple();
  //  SetCurrentCouple(couple); 
  //  currentMaterialIndex = couple->GetIndex();

  currentKinEnergy = kineticEnergy; // dp->GetKineticEnergy();
  //  currentRange = GetRange(particle,currentKinEnergy,couple);
  currentRange = GetRange(material,currentKinEnergy);
  //  lambda0 = GetTransportMeanFreePath(particle,currentKinEnergy);
  //  lambda0 = GetTransportMeanFreePath(currentKinEnergy);
  lambda0 = xlambda;

  if(tPathLength > currentRange) { tPathLength = currentRange; }

  // stop here if small range particle
  if(inside || tPathLength < tlimitminfix) { 
    return ConvertTrueToGeom(material,tPathLength,currentMinimalStep); 
  }
  
  //@@@  presafety = sp->GetSafety();
  presafety = 1.0*mm;

  // far from geometry boundary
  if(currentRange < presafety)
    {
      inside = true;
      return ConvertTrueToGeom(material,tPathLength,currentMinimalStep);  
    }

  //@@@ steppingAlgorithm = fUseSafety;
  //@@@ see constructor of G4VMscModel.cc

  // standard  version
  //

  /*
  if (steppingAlgorithm == fUseDistanceToBoundary)
    {
      //compute geomlimit and presafety 
      geomlimit = ComputeGeomLimit(track, presafety, currentRange);

      // is it far from boundary ?
      if(currentRange < presafety)
	{
	  inside = true;
	  return ConvertTrueToGeom(tPathLength, currentMinimalStep);   
	}

      smallstep += 1.;
      insideskin = false;

      if(firstStep || (stepStatus == fGeomBoundary))
        {
          rangeinit = currentRange;
          if(firstStep) smallstep = 1.e10;
          else  smallstep = 1.;

          //define stepmin here (it depends on lambda!)
          //rough estimation of lambda_elastic/lambda_transport
          G4double rat = currentKinEnergy/MeV ;
          rat = 1.e-3/(rat*(10.+rat)) ;
          //stepmin ~ lambda_elastic
          stepmin = rat*lambda0;
          skindepth = skin*stepmin;
          //define tlimitmin
          tlimitmin = 10.*stepmin;
          if(tlimitmin < tlimitminfix) tlimitmin = tlimitminfix;

          // constraint from the geometry
          if((geomlimit < geombig) && (geomlimit > geommin))
            {
              // geomlimit is a geometrical step length
              // transform it to true path length (estimation)
              if((1.-geomlimit/lambda0) > 0.)
                geomlimit = -lambda0*log(1.-geomlimit/lambda0)+tlimitmin ;

              if(stepStatus == fGeomBoundary)
                tgeom = geomlimit/facgeom;
              else
                tgeom = 2.*geomlimit/facgeom;
            }
            else
              tgeom = geombig;
        }


      //step limit 
      tlimit = facrange*rangeinit;              

      //lower limit for tlimit
      if(tlimit < tlimitmin) tlimit = tlimitmin;

      if(tlimit > tgeom) tlimit = tgeom;

      // shortcut
      if((tPathLength < tlimit) && (tPathLength < presafety) &&
         (smallstep > skin) && (tPathLength < geomlimit-0.999*skindepth))
	return ConvertTrueToGeom(tPathLength, currentMinimalStep);   

      // step reduction near to boundary
      if(smallstep <= skin)
	{
	  tlimit = stepmin;
	  insideskin = true;
	}
      else if(geomlimit < geombig)
	{
	  if(geomlimit > skindepth)
	    {
	      if(tlimit > geomlimit-0.999*skindepth)
		tlimit = geomlimit-0.999*skindepth;
	    }
	  else
	    {
	      insideskin = true;
	      if(tlimit > stepmin) tlimit = stepmin;
	    }
	}

      if(tlimit < stepmin) tlimit = stepmin;

      // randomize 1st step or 1st 'normal' step in volume
      if(firstStep || ((smallstep == skin) && !insideskin)) 
        { 
          G4double temptlimit = tlimit;
          if(temptlimit > tlimitmin)
          {
            do {
              temptlimit = G4RandGauss::shoot(tlimit,0.3*tlimit);        
               } while ((temptlimit < tlimitmin) || 
                        (temptlimit > 2.*tlimit-tlimitmin));
          }
          else
            temptlimit = tlimitmin;
          if(tPathLength > temptlimit) tPathLength = temptlimit;
        }
      else
        {  
          if(tPathLength > tlimit) tPathLength = tlimit  ; 
        }

    }
    // for 'normal' simulation with or without magnetic field 
    //  there no small step/single scattering at boundaries
  else if(steppingAlgorithm == fUseSafety)
  */
    {
      // compute presafety again if presafety <= 0 and no boundary
      // i.e. when it is needed for optimization purposes

      //@@@ implement SafetyHelper !!!
      //@@@      if((stepStatus != fGeomBoundary) && (presafety < tlimitminfix)) 
      //@@@        presafety = ComputeSafety(sp->GetPosition(),tPathLength); 
      //@@@ use temporary
      presafety = 1.0*mm;

      // is far from boundary
      if(currentRange < presafety)
        {
          inside = true;
	  return ConvertTrueToGeom(material,tPathLength,currentMinimalStep);    
        }

      //@@@      if(firstStep || stepStatus == fGeomBoundary)
      if(firstStep)
      {
        rangeinit = currentRange;
        fr = facrange;
        // 9.1 like stepping for e+/e- only (not for muons,hadrons)
        if(mass < masslimite) 
        {
          if(lambda0 > currentRange)
            rangeinit = lambda0;
          if(lambda0 > lambdalimit)
            fr *= 0.75+0.25*lambda0/lambdalimit;
        }

        //lower limit for tlimit
        G4double rat = currentKinEnergy/MeV ;
        rat = 1.e-3/(rat*(10.+rat)) ;
        tlimitmin = 10.*lambda0*rat;
        if(tlimitmin < tlimitminfix) tlimitmin = tlimitminfix;
      }
      //step limit
      tlimit = fr*rangeinit;               

      if(tlimit < facsafety*presafety)
        tlimit = facsafety*presafety;

      //lower limit for tlimit
      if(tlimit < tlimitmin) tlimit = tlimitmin;
      
      if(tPathLength > tlimit) tPathLength = tlimit;

    }
  
  // version similar to 7.1 (needed for some experiments)
    /*
  else
    {
      if (stepStatus == fGeomBoundary)
	{
	  if (currentRange > lambda0) tlimit = facrange*currentRange;
	  else                        tlimit = facrange*lambda0;

	  if(tlimit < tlimitmin) tlimit = tlimitmin;
	  if(tPathLength > tlimit) tPathLength = tlimit;
	}
    }
    */

    return ConvertTrueToGeom(material,tPathLength,currentMinimalStep);
}

//G4double G4VMscModel::ConvertTrueToGeom(G4double& tlength, 
FQUALIFIER
G4double GPUrbanMscModel95::ConvertTrueToGeom(GPMaterial* material,
					      G4double& tlength,
					      G4double& glength)
{
  glength = ComputeGeomPathLength(material);
  // should return true length 
  return tlength;
}

FQUALIFIER
G4double GPUrbanMscModel95::ComputeTrueStepLength(G4double geomStepLength)
{
  // step defined other than transportation 
  if(geomStepLength == zPathLength)
    { return tPathLength; }

  zPathLength = geomStepLength;

  // t = z for very small step
  if(geomStepLength < tlimitminfix) { 
    tPathLength = geomStepLength; 
  
  // recalculation
  } else {

    G4double tlength = geomStepLength;
    if((geomStepLength > lambda0*tausmall) && !insideskin) {


      if(par1 <  0.) {
	tlength = -lambda0*log(1.-geomStepLength/lambda0) ;

      } else {
	if(par1*par3*geomStepLength < 1.) {
	  tlength = (1.-exp(log(1.-par1*par3*geomStepLength)/par3))/par1 ;
	} else {
	  tlength = currentRange;
	}
      }
      if(tlength < geomStepLength)   { tlength = geomStepLength; }
      else if(tlength > tPathLength) { tlength = tPathLength; }
    }  
    tPathLength = tlength; 
  }

  return tPathLength;
}

FQUALIFIER
G4double GPUrbanMscModel95::ComputeTheta0(G4double trueStepLength,
					  G4double KineticEnergy)
{
  // for all particles take the width of the central part
  //  from a  parametrization similar to the Highland formula
  // ( Highland formula: Particle Physics Booklet, July 2002, eq. 26.10)
  const G4double c_highland = 13.6*MeV ;
  G4double betacp = sqrt(currentKinEnergy*(currentKinEnergy+2.*mass)*
                         KineticEnergy*(KineticEnergy+2.*mass)/
                      ((currentKinEnergy+mass)*(KineticEnergy+mass)));
  G4double y = trueStepLength/currentRadLength;
  //  G4double theta0 = c_highland*abs(charge)*sqrt(y)/betacp;
  G4double theta0 = c_highland*sqrt(y)/betacp;
  y = log(y);
  // correction factor from e- scattering data
  G4double corr = coeffth1+coeffth2*y;                

  theta0 *= corr ;                                               

  return theta0;
}

FQUALIFIER
GPThreeVector
GPUrbanMscModel95::SampleScattering(GPMaterial* material,
				    GPThreeVector oldDirection,
				    G4double safety)
{

  //  fDisplacement.set(0.0,0.0,0.0);
  GPThreeVector fDisplacement = GPThreeVector_create(0.0,0.0,0.0);

  G4double kineticEnergy = currentKinEnergy;

  if (tPathLength > currentRange*dtrl) {
    //    kineticEnergy = GetEnergy(particle,currentRange-tPathLength,couple);
    kineticEnergy = GetEnergy(material,currentRange-tPathLength);
  } else {
    //    kineticEnergy -= tPathLength*GetDEDX(particle,currentKinEnergy,couple);
    kineticEnergy -= tPathLength*GetDEDX(currentKinEnergy);
  }

  if((kineticEnergy <= eV) || (tPathLength <= tlimitminfix) ||
     (tPathLength/tausmall < lambda0)) { return fDisplacement; }

  G4double cth  = SampleCosineTheta(material,tPathLength,kineticEnergy);

  // protection against 'bad' cth values
  if(fabs(cth) > 1.) { return fDisplacement; }

  // extra protection agaist high energy particles backscattered 
  // do Gaussian central scattering
  if(kineticEnergy > 5*GeV && cth < 0.9) {
    //@@@Print warning
  }

  G4double sth  = sqrt((1.0 - cth)*(1.0 + cth));
  G4double phi  = twopi*GPUniformRand(fDevStates,fThreadId);
  G4double dirx = sth*cos(phi);
  G4double diry = sth*sin(phi);

  //  G4ThreeVector newDirection(dirx,diry,cth);
  GPThreeVector newDirection = GPThreeVector_create(dirx,diry,cth);
  //  newDirection.rotateUz(oldDirection);
  GPThreeVector_rotateUz(&newDirection,oldDirection);

  pParticleChange->SetProposedMomentumDirection(newDirection);

  if (latDisplasment && safety > tlimitminfix) {

    G4double r = SampleDisplacement();

    if(r > 0.)
      {
        G4double latcorr = LatCorrelation();
        if(latcorr > r) latcorr = r;

        // sample direction of lateral displacement
        // compute it from the lateral correlation
        G4double Phi = 0.;
        if(abs(r*sth) < latcorr)
          Phi  = twopi*GPUniformRand(fDevStates,fThreadId);
        else
        {
          G4double psi = acos(latcorr/(r*sth));
          if(GPUniformRand(fDevStates,fThreadId) < 0.5)
            Phi = phi+psi;
          else
            Phi = phi-psi;
        }

        dirx = cos(Phi);
        diry = sin(Phi);

	//        fDisplacement.set(r*dirx,r*diry,0.0);
	//        fDisplacement.rotateUz(oldDirection);
        GPThreeVector_set(&fDisplacement,r*dirx,r*diry,0.0);
        GPThreeVector_rotateUz(&fDisplacement,oldDirection);
      }
  }

  return fDisplacement;
}

FQUALIFIER
G4double GPUrbanMscModel95::SampleCosineTheta(GPMaterial* material,
					      G4double trueStepLength,
					      G4double KineticEnergy)
{
  G4double cth = 1. ;
  G4double tau = trueStepLength/lambda0;
  currentTau   = tau;
  lambdaeff    = lambda0;

  //  Zeff = couple->GetMaterial()->GetTotNbOfElectPerVolume()/
  //         couple->GetMaterial()->GetTotNbOfAtomsPerVolume() ;
  Zeff = GPMaterial_GetTotNbOfElectPerVolume(material)/
         GPMaterial_GetTotNbOfAtomsPerVolume(material) ;

  if(Zold != Zeff)  
    UpdateCache();

  if(insideskin)
  {
    //no scattering, single or plural scattering
    G4double mean = trueStepLength/stepmin ;

    G4int n = GPPoisson(fDevStates,fThreadId,mean);
    if(n > 0)
    {
      //screening (Moliere-Bethe)
      G4double mom2 = KineticEnergy*(2.*mass+KineticEnergy);
      G4double beta2 = mom2/((KineticEnergy+mass)*(KineticEnergy+mass));
      G4double ascr = scr1/mom2;
      ascr *= 1.13+scr2/beta2;
      G4double ascr1 = 1.+2.*ascr;
      G4double bp1=ascr1+1.;
      G4double bm1=ascr1-1.;

      // single scattering from screened Rutherford x-section
      G4double ct,st,phi;
      G4double sx=0.,sy=0.,sz=0.;
      for(G4int i=1; i<=n; i++)
      {
        ct = ascr1-bp1*bm1/(2.*GPUniformRand(fDevStates,fThreadId)+bm1);
        if(ct < -1.) ct = -1.;
        if(ct >  1.) ct =  1.; 
        st = sqrt(1.-ct*ct);
        phi = twopi*GPUniformRand(fDevStates,fThreadId);
        sx += st*cos(phi);
        sy += st*sin(phi);
        sz += ct;
      }
      cth = sz/sqrt(sx*sx+sy*sy+sz*sz);
    }
  }
  else
  {
    //    G4double lambda1 = GetTransportMeanFreePath(particle,KineticEnergy);
    G4double lambda1 = GetTransportMeanFreePath(KineticEnergy);
    if(fabs(lambda1/lambda0 - 1) > 0.01 && lambda1 > 0.)
    {
      // mean tau value
      tau = trueStepLength*log(lambda0/lambda1)/(lambda0-lambda1);
    }

    currentTau = tau ;
    lambdaeff = trueStepLength/currentTau;
    currentRadLength = GPMaterial_GetRadlen(material);

    if (tau >= taubig) cth = -1.+2.*GPUniformRand(fDevStates,fThreadId);
    else if (tau >= tausmall)
    {
      const G4double numlim = 0.01;
      G4double xmeanth, x2meanth;
      if(tau < numlim) {
	xmeanth = 1.0 - tau*(1.0 - 0.5*tau);
        x2meanth= 1.0 - tau*(5.0 - 6.25*tau)/3.;
      } else {
	xmeanth = exp(-tau);
	x2meanth = (1.+2.*exp(-2.5*tau))/3.;
      }
      G4double relloss = 1.-KineticEnergy/currentKinEnergy;

      if(relloss > 0.50) { //rellossmax    = 0.50
        return SimpleScattering(xmeanth,x2meanth);
      }
      G4double theta0 = ComputeTheta0(trueStepLength,KineticEnergy);

      // protection for very small angles
      G4double theta2 = theta0*theta0;

      if(theta2 < tausmall) { return cth; }
    
      if(theta0 > pi/6.) { //theta0max     = pi/6.
        return SimpleScattering(xmeanth,x2meanth);
      }

      G4double x = theta2*(1.0 - theta2/12.);
      if(theta2 > numlim) {
	G4double sth = 2*sin(0.5*theta0);
	x = sth*sth;
      }

      // parameter for tail
      G4double ltau= log(tau);
      G4double u   = exp(ltau/6.);
      G4double xx  = log(lambdaeff/currentRadLength);
      G4double xsi = coeffc1+u*(coeffc2+coeffc3*u)+coeffc4*xx;
      G4double   c = xsi;

      //correction of tail for high energy/small step
      if(ltau < -10.63)
        { c *= (0.016*ltau+1.17008); }

      // tail should not be too big
      if(c < 1.9) { c = 1.9; }

      if(abs(c-3.) < 0.001)  c = 3.001;      
      if(abs(c-2.) < 0.001)  c = 2.001;      
      if(abs(c-1.) < 0.001)  c = 1.001;      

      G4double c1 = c-1.;

      G4double ea = exp(-xsi);
      G4double eaa = 1.-ea ;
      G4double xmean1 = 1.-(1.-(1.+xsi)*ea)*x/eaa;
      G4double x0 = 1. - xsi*x;

      if(xmean1 <= 0.999*xmeanth) {
        return SimpleScattering(xmeanth,x2meanth);
      }
      //from continuity of derivatives
      G4double b = 1.+(c-xsi)*x;

      G4double b1 = b+1.;
      G4double bx = c*x;

      G4double eb1 = pow(b1,c1);
      G4double ebx = pow(bx,c1);
      G4double d = ebx/eb1;

      G4double xmean2 = (x0 + d - (bx - b1*d)/(c-2.))/(1. - d);

      G4double f1x0 = ea/eaa;
      G4double f2x0 = c1/(c*(1. - d));
      G4double prob = f2x0/(f1x0+f2x0);

      G4double qprob = xmeanth/(prob*xmean1+(1.-prob)*xmean2);

      // sampling of costheta
      if(GPUniformRand(fDevStates,fThreadId) < qprob)
      {
        G4double var = 0;
        if(GPUniformRand(fDevStates,fThreadId) < prob) {
          cth = 1.+log(ea+GPUniformRand(fDevStates,fThreadId)*eaa)*x;
        } else {
          var = (1.0 - d)*GPUniformRand(fDevStates,fThreadId);
          if(var < numlim*d) {
            var /= (d*c1); 
            cth = -1.0 + var*(1.0 - 0.5*var*c)*(2. + (c - xsi)*x);
	  } else {
	    cth = 1. + x*(c - xsi - c*pow(var + d, -1.0/c1));
	  }
	}
	if(KineticEnergy > 5*GeV && cth < 0.9) {
	  //@@@Print Warning
	}
      }
      else {
        cth = -1.+2.*GPUniformRand(fDevStates,fThreadId);
	if(KineticEnergy > 5*GeV) {
	  //@@@Print Warning
	}
      }
    }
  }
  return cth ;
}

FQUALIFIER
G4double GPUrbanMscModel95::SimpleScattering(G4double xmeanth,G4double x2meanth)
{
  // 'large angle scattering'
  // 2 model functions with correct xmean and x2mean
  G4double a = (2.*xmeanth+9.*x2meanth-3.)/(2.*xmeanth-3.*x2meanth+1.);
  G4double prob = (a+2.)*xmeanth/a;

  // sampling
  G4double cth = 1.;
  if(GPUniformRand(fDevStates,fThreadId) < prob)
    cth = -1.+2.*exp(log(GPUniformRand(fDevStates,fThreadId))/(a+1.));
  else
    cth = -1.+2.*GPUniformRand(fDevStates,fThreadId);
  return cth;
}

FQUALIFIER
G4double GPUrbanMscModel95::SampleDisplacement()
{
  G4double r = 0.0;
  if ((currentTau >= tausmall) && !insideskin) {
    G4double rmax = sqrt((tPathLength-zPathLength)*(tPathLength+zPathLength));
    r = rmax*exp(log(GPUniformRand(fDevStates,fThreadId))/3.);
  }
  return r;
}

FQUALIFIER
G4double GPUrbanMscModel95::LatCorrelation()
{
  const G4double kappa = 2.5;
  const G4double kappami1 = kappa-1.;

  G4double latcorr = 0.;
  if((currentTau >= tausmall) && !insideskin)
  {
    if(currentTau < taulim)
      latcorr = lambdaeff*kappa*currentTau*currentTau*
                (1.-(kappa+1.)*currentTau/3.)/3.;
    else
    {
      G4double etau = 0.;
      if(currentTau < taubig) etau = exp(-currentTau);
      latcorr = -kappa*currentTau;
      latcorr = exp(latcorr)/kappami1;
      latcorr += 1.-kappa*etau/kappami1 ;
      latcorr *= 2.*lambdaeff/3. ;
    }
  }

  return latcorr;
}

FQUALIFIER
void GPUrbanMscModel95::SetParticle()
{
  //  if (p != particle) {
  //    particle = p;
  //@@@paricle -> electron or positron only
  mass = electron_mass_c2; //p->GetPDGMass();
  //charge = p->GetPDGCharge()/CLHEP::eplus;
  charge = -1.0;
  //  }
}

FQUALIFIER
void GPUrbanMscModel95::UpdateCache()
{
  lnZ = log(Zeff);
  // correction in theta0 formula                                             
  coeffth1 = (1. - 8.7780e-2/Zeff)*(0.87 + 0.03*lnZ);
  coeffth2 = (4.0780e-2 + 1.7315e-4*Zeff)*(0.87 + 0.03*lnZ);

  // tail parameters                                                          
  G4double Z13 = exp(lnZ/3.);
  coeffc1  = 2.3785    - Z13*(4.1981e-1 - Z13*6.3100e-2);
  coeffc2  = 4.7526e-1 + Z13*(1.7694    - Z13*3.3885e-1);
  coeffc3  = 2.3683e-1 - Z13*(1.8111    - Z13*3.2774e-1);
  coeffc4  = 1.7888e-2 + Z13*(1.9659e-2 - Z13*2.6664e-3);

  // for single scattering                                                    
  Z2   = Zeff*Zeff;
  Z23  = Z13*Z13;
  scr1 = scr1ini*Z23;
  //  scr2 = scr2ini*Z2*ChargeSquare;
  scr2 = scr2ini*Z2*charge*charge;

  Zold = Zeff;
}

//G4VMscModel::GetDEDX(const G4ParticleDefinition* part,
//                     G4double kinEnergy, const G4MaterialCutsCouple* couple)
FQUALIFIER 
G4double GPUrbanMscModel95::GetDEDX(G4double kinEnergy)
{
  G4double x;
  //  if(ionisation) { x = ionisation->GetDEDX(kinEnergy, couple); }
  //  else { 
  //  G4double q = part->GetPDGCharge()/CLHEP::eplus;
    G4double q = -1.0;
    x = dedx*q*q;
    //  }
  return x;
}

//G4VMscModel::GetRange(const G4ParticleDefinition* part,
//                      G4double kinEnergy, const G4MaterialCutsCouple* couple)
FQUALIFIER 
G4double GPUrbanMscModel95::GetRange(GPMaterial* material,
				     G4double kinEnergy)
{
  //  if(ionisation) { 
  //    localrange = ionisation->GetRangeForLoss(kinEnergy, couple); 
  //  } else { 
  //  G4double q = part->GetPDGCharge()/CLHEP::eplus;
    G4double q = -1.0;
    localrange = kinEnergy/(dedx*q*q*GPMaterial_GetDensity(material)); 
    localtkin  = kinEnergy;
    //  }
  return localrange;
}

//G4VMscModel::GetEnergy(const G4ParticleDefinition* part,
//                       G4double range, const G4MaterialCutsCouple* couple)
FQUALIFIER 
G4double  GPUrbanMscModel95::GetEnergy(GPMaterial* material,
				       G4double range)
{
  G4double e;
  //  if(ionisation) { e = ionisation->GetKineticEnergy(range, couple); }
  //  else { 
    e = localtkin;
    if(localrange > range) {
      //      G4double q = part->GetPDGCharge()/CLHEP::eplus;
      G4double q = -1.0;
      e -= (localrange - range)*dedx*q*q*GPMaterial_GetDensity(material); 
    } 
    //  }
  return e;
}

//G4VMscModel::GetTransportMeanFreePath(const G4ParticleDefinition* part,
//                                      G4double ekin)
FQUALIFIER
G4double GPUrbanMscModel95::GetTransportMeanFreePath(G4double ekin)
{
  G4double x;
  //@@@ using xSectionTable and (*theDensityFactor)[idx] = 1.0 for PbWO4
  //@@@ xSectionTable for msc is the lambda table:
  //  if(xSectionTable) {
  //  G4int idx = CurrentCouple()->GetIndex();
  //  x = (*xSectionTable)[(*theDensityIdx)[idx]]->Value(ekin)
  //    *(*theDensityFactor)[idx]/(ekin*ekin);
  x = theLambdaTable->physicsVectors[1].Value(ekin)*1.0/(ekin*ekin);
  //  } else { 
  //    x = CrossSectionPerVolume(CurrentCouple()->GetMaterial(), part, ekin, 
  //                              0.0, DBL_MAX); 
  //  }
  if(0.0 >= x) { x = DBL_MAX; }
  else { x = 1.0/x; }
  return x;
}

FQUALIFIER
void GPUrbanMscModel95::SetLambdaTable(GPPhysicsTable* lambdaTable)
{
  theLambdaTable = lambdaTable;
}

//---------------------------------------------------------------------------
//
// G4VEmModel
//
//---------------------------------------------------------------------------

FQUALIFIER 
void GPUrbanMscModel95::SetParticleChange(GPVParticleChange* p)
//                                          GPUniversalFluctuation* f)
{
  if(p && pParticleChange != p) { pParticleChange = p; }
  //  flucModel = f;
}

FQUALIFIER 
G4double GPUrbanMscModel95::HighEnergyLimit()
{
  return highLimit;
}

FQUALIFIER 
G4double GPUrbanMscModel95::LowEnergyLimit()
{
  return lowLimit;
}

FQUALIFIER 
G4double GPUrbanMscModel95::HighEnergyActivationLimit()
{
  return eMaxActive;
}

FQUALIFIER 
G4double GPUrbanMscModel95::LowEnergyActivationLimit()
{
  return eMinActive;
}

FQUALIFIER
void GPUrbanMscModel95::SetHighEnergyLimit(G4double val)
{
  highLimit = val;
}

FQUALIFIER
void GPUrbanMscModel95::SetLowEnergyLimit(G4double val)
{
  lowLimit = val;
}

FQUALIFIER 
void GPUrbanMscModel95::SetActivationHighEnergyLimit(G4double val)
{
  eMaxActive = val;
}

FQUALIFIER 
void GPUrbanMscModel95::SetActivationLowEnergyLimit(G4double val)
{
  eMinActive = val;
}

FQUALIFIER 
G4bool GPUrbanMscModel95::IsActive(G4double kinEnergy)
{
  return (kinEnergy >= eMinActive && kinEnergy <= eMaxActive);
}

