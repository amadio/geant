#ifndef GPCLHEPRandom_H
#define GPCLHEPRandom_H

//CLHEP::RandFlat

namespace GPCLHEPRandFlat {

inline
FQUALIFIER double flat() 
{
#ifdef GPNONRANDOM
  return 0.123456;
#else
  return (G4double)rand()/RAND_MAX;
#endif
} 

};

//CLHEP::RandGauss
namespace GPCLHEPRandGauss {

CONSTTYPE static bool set_st = false;
CONSTTYPE static double nextGauss_st = 0.0;

inline FQUALIFIER static bool getFlag() {return set_st;}
inline FQUALIFIER static void setFlag( bool val ) {set_st = val;}
inline FQUALIFIER static double getVal() {return nextGauss_st;}
inline FQUALIFIER static void setVal( double nextVal ) {nextGauss_st = nextVal;}

inline FQUALIFIER double shoot()
{
  // Gaussian random numbers are generated two at the time, so every other
  // time this is called we just return a number generated the time before.

  if ( getFlag() ) {
    setFlag(false);
    double x = getVal();
    return x; 
    // return getVal();
  } 

  double r;
  double v1,v2,fac,val;
  //  HepRandomEngine* anEngine = HepRandom::getTheEngine();

  do {
    //    v1 = 2.0 * anEngine->flat() - 1.0;
    //    v2 = 2.0 * anEngine->flat() - 1.0;
    v1 = 2.0 * GPCLHEPRandFlat::flat() - 1.0;
    v2 = 2.0 * GPCLHEPRandFlat::flat() - 1.0;
    r = v1*v1 + v2*v2;
  } while ( r > 1.0 );

  fac = sqrt(-2.0*log(r)/r);
  val = v1*fac;
  setVal(val);
  setFlag(true);
  return v2*fac;
}

};

// CLHEP::RandPoisson
namespace GPCLHEPRandPoisson {

// Initialisation of static data

CONSTTYPE static double status_st[3] = {0., 0., 0.};
CONSTTYPE static double oldm_st = -1.0;
CONSTTYPE static const double meanMax_st = 2.0E9;

inline FQUALIFIER static  double  getOldMean() {return oldm_st;}
inline FQUALIFIER static  double  getMaxMean() {return meanMax_st;}
inline FQUALIFIER static  void    setOldMean( double val ){oldm_st = val;}
inline FQUALIFIER static  double* getPStatus() {return status_st;}

inline FQUALIFIER static void setPStatus(double sq, double alxm, double g1) 
{
  status_st[0] = sq; status_st[1] = alxm; status_st[2] = g1;
}

inline
FQUALIFIER double gammln(double xx) 
{

// Returns the value ln(Gamma(xx) for xx > 0.  Full accuracy is obtained for 
// xx > 1. For 0 < xx < 1. the reflection formula (6.1.4) can be used first.
// (Adapted from Numerical Recipes in C)

  static double cof[6] = {76.18009172947146,-86.50532032941677,
			  24.01409824083091, -1.231739572450155,
			  0.1208650973866179e-2, -0.5395239384953e-5};
  int j;
  double x = xx - 1.0;
  double tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  double ser = 1.000000000190015;

  for ( j = 0; j <= 5; j++ ) {
    x += 1.0;
    ser += cof[j]/x;
  }
#ifdef GPNONRANDOM
  return 0.123456;
#else
  return -tmp + log(2.5066282746310005*ser);
#endif
}

inline
FQUALIFIER static double normal ()  
{
  double r;
  double v1,v2,fac;
  do {
    //    v1 = 2.0 * eptr->flat() - 1.0;
    //    v2 = 2.0 * eptr->flat() - 1.0;
    v1 = 2.0 * GPCLHEPRandFlat::flat() - 1.0;
    v2 = 2.0 * GPCLHEPRandFlat::flat() - 1.0;
    r = v1*v1 + v2*v2;
  } while ( r > 1.0 );

  fac = sqrt(-2.0*log(r)/r);
  return v2*fac;
}

inline
FQUALIFIER long shoot(double xm) 
{

// Returns as a floating-point number an integer value that is a random
// deviation drawn from a Poisson distribution of mean xm, using flat()
// as a source of uniform random numbers.
// (Adapted from Numerical Recipes in C)

  double em, t, y;
  double sq, alxm, g1;
  double om = getOldMean();
  //  HepRandomEngine* anEngine = HepRandom::getTheEngine();

  double* status = getPStatus();
  sq = status[0];
  alxm = status[1];
  g1 = status[2];

  if( xm == -1 ) return 0;
  if( xm < 12.0 ) {
    if( xm != om ) {
      setOldMean(xm);
      g1 = exp(-xm);
    }
    em = -1;
    t = 1.0;
    do {
      em += 1.0;
      //      t *= anEngine->flat();
      t *= GPCLHEPRandFlat::flat();
    } while( t > g1 );
  }
  else if ( xm < getMaxMean() ) {
    if ( xm != om ) {
      setOldMean(xm);
      sq = sqrt(2.0*xm);
      alxm = log(xm);
      g1 = xm*alxm - gammln(xm + 1.0);
    }
    do {
      do {
	//        y = tan(CLHEP::pi*anEngine->flat());
        y = tan(3.14159265358979323846*GPCLHEPRandFlat::flat());
        em = sq*y + xm;
      } while( em < 0.0 );
      em = floor(em);
      t = 0.9*(1.0 + y*y)* exp(em*alxm - gammln(em + 1.0) - g1);
    } while( GPCLHEPRandFlat::flat() > t );
    //    } while( anEngine->flat() > t );
  }
  else {
    //    em = xm + sqrt(xm) * normal (anEngine); 
    em = xm + sqrt(xm) * normal();
    if ( static_cast<long>(em) < 0 ) 
      em = static_cast<long>(xm) >= 0 ? xm : getMaxMean();
  }    
  setPStatus(sq,alxm,g1);
  return long(em);
}

};

//CLHEP::RandGamma
namespace GPCLHEPRandGamma {

inline
FQUALIFIER double genGamma(double a, double lambda ) 
{
static double aa = -1.0, aaa = -1.0, b, c, d, e, r, s, si, ss, q0,
       q1 = 0.0416666664, q2 =  0.0208333723, q3 = 0.0079849875,
       q4 = 0.0015746717, q5 = -0.0003349403, q6 = 0.0003340332,
       q7 = 0.0006053049, q8 = -0.0004701849, q9 = 0.0001710320,
       a1 = 0.333333333,  a2 = -0.249999949,  a3 = 0.199999867,
       a4 =-0.166677482,  a5 =  0.142873973,  a6 =-0.124385581,
       a7 = 0.110368310,  a8 = -0.112750886,  a9 = 0.104089866,
       e1 = 1.000000000,  e2 =  0.499999994,  e3 = 0.166666848,
       e4 = 0.041664508,  e5 =  0.008345522,  e6 = 0.001353826,
       e7 = 0.000247453;
double gds,p,q,t,sign_u,u,v,w,x;
double v1,v2,v12;

// Check for invalid input values

 if( a <= 0.0 ) return (-1.0);
 if( lambda <= 0.0 ) return (-1.0);

 if (a < 1.0)
   {          // CASE A: Acceptance rejection algorithm gs
    b = 1.0 + 0.36788794412 * a;       // Step 1
    for(;;)
      {
       p = b * GPCLHEPRandFlat::flat();
       if (p <= 1.0)
          {                            // Step 2. Case gds <= 1
           gds = exp(log(p) / a);
           if (log(GPCLHEPRandFlat::flat()) <= -gds) return(gds/lambda);
          }
       else
          {                            // Step 3. Case gds > 1
           gds = - log ((b - p) / a);
           if (log(GPCLHEPRandFlat::flat()) <= ((a - 1.0)*log(gds))){
	     return(gds/lambda);
	   }
          }
      }
   }
 else
   {          // CASE B: Acceptance complement algorithm gd
    if (a != aa)
       {                               // Step 1. Preparations
        aa = a;
        ss = a - 0.5;
        s = sqrt(ss);
        d = 5.656854249 - 12.0 * s;
       }
                                              // Step 2. Normal deviate
    do {
      v1 = 2.0 * GPCLHEPRandFlat::flat() - 1.0;
      v2 = 2.0 * GPCLHEPRandFlat::flat() - 1.0;
      v12 = v1*v1 + v2*v2;
    } while ( v12 > 1.0 );
    t = v1*sqrt(-2.0*log(v12)/v12);
    x = s + 0.5 * t;
    gds = x * x;
    if (t >= 0.0) return(gds/lambda);         // Immediate acceptance

    u = GPCLHEPRandFlat::flat();            // Step 3. Uniform random number
    if (d * u <= t * t * t) return(gds/lambda); // Squeeze acceptance

    if (a != aaa)
       {                               // Step 4. Set-up for hat case
        aaa = a;
        r = 1.0 / a;
        q0 = ((((((((q9 * r + q8) * r + q7) * r + q6) * r + q5) * r + q4) *
                          r + q3) * r + q2) * r + q1) * r;
        if (a > 3.686)
           {
            if (a > 13.022)
               {
                b = 1.77;
                si = 0.75;
                c = 0.1515 / s;
               }
            else
               {
                b = 1.654 + 0.0076 * ss;
                si = 1.68 / s + 0.275;
                c = 0.062 / s + 0.024;
               }
           }
        else
           {
            b = 0.463 + s - 0.178 * ss;
            si = 1.235;
            c = 0.195 / s - 0.079 + 0.016 * s;
           }
       }
    if (x > 0.0)                       // Step 5. Calculation of q
       {
        v = t / (s + s);               // Step 6.
        if (fabs(v) > 0.25)
           {
            q = q0 - s * t + 0.25 * t * t + (ss + ss) * log(1.0 + v);
           }
        else
           {
            q = q0 + 0.5 * t * t * ((((((((a9 * v + a8) * v + a7) * v + a6) *
                            v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
           }                // Step 7. Quotient acceptance
        if (log(1.0 - u) <= q) return(gds/lambda);
       }

    for(;;)
       {                    // Step 8. Double exponential deviate t
        do
        {
         e = -log(GPCLHEPRandFlat::flat());
         u = GPCLHEPRandFlat::flat();
         u = u + u - 1.0;
         sign_u = (u > 0)? 1.0 : -1.0;
         t = b + (e * si) * sign_u;
        }
        while (t <= -0.71874483771719);   // Step 9. Rejection of t
        v = t / (s + s);                  // Step 10. New q(t)
        if (fabs(v) > 0.25)
           {
            q = q0 - s * t + 0.25 * t * t + (ss + ss) * log(1.0 + v);
           }
        else
           {
            q = q0 + 0.5 * t * t * ((((((((a9 * v + a8) * v + a7) * v + a6) *
                            v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
           }
        if (q <= 0.0) continue;           // Step 11.
        if (q > 0.5)
           {
            w = exp(q) - 1.0;
           }
        else
           {
            w = ((((((e7 * q + e6) * q + e5) * q + e4) * q + e3) * q + e2) *
                                     q + e1) * q;
           }                    // Step 12. Hat acceptance
        if ( c * u * sign_u <= w * exp(e - 0.5 * t * t))
           {
            x = s + 0.5 * t;
            return(x*x/lambda);
           }
       }
   }
}

inline
FQUALIFIER double shoot(double k, double lambda) 
{
  return genGamma(k, lambda);
}

};

#endif

