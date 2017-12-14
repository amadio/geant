

#include "Spline.h"

// cstdlib is for std::abs(int)
#include <cstdlib>

namespace geantphysics {

// spline does not own fXdata and fYdata
Spline::Spline(double *xdata, double *ydata, int numdata, bool isacsecderis, bool isconstrained) {
  fNPoints = numdata;
  fXdata   = xdata;
  fYdata   = ydata;
  fA.resize(fNPoints);
  fB.resize(fNPoints);
  fC.resize(fNPoints);
  fD.resize(fNPoints);
  //  fSecondDerivatives.resize(fNPoints);
  fIsConstrained = isconstrained;
  fIsAcSecDerivs = isacsecderis;
  if (fIsConstrained) {
    fIsAcSecDerivs = false;
  }
  // set some members used for the B-search
  if (*(fXdata)>*(fXdata+(fNPoints-1))) {
    fDirection = 1;          // decreasing x grid
    fLowerm    = fNPoints-1;
    fUpperm    = -1;
    if (fIsAcSecDerivs) {  // can be used only for increasing x data
      fIsAcSecDerivs = false;
    }
  } else {
    fDirection = 0;         // increasing x grid
    fLowerm    = -1;
    fUpperm    = fNPoints-1;
  }
  if (fIsConstrained) {
    SetParametersConstrained(); // constrained natural cubic spline
  } else {
    SetParameters();            // natural cubic spline with possible accurate second derivatives
  }
}

// spline does not own fXdata and fYdata
void Spline::SetUpSpline(double *xdata, double *ydata, int numdata, bool isacsecderis, bool isconstrained) {
  fNPoints = numdata;
  fXdata   = xdata;
  fYdata   = ydata;
  fA.resize(fNPoints);
  fB.resize(fNPoints);
  fC.resize(fNPoints);
  fD.resize(fNPoints);
//  fSecondDerivatives.resize(fNPoints);
  fIsConstrained = isconstrained;
  fIsAcSecDerivs = isacsecderis;
  if (fIsConstrained) {
    fIsAcSecDerivs = false;
  }
  // set some members used for the B-search
  if (*(fXdata)>*(fXdata+(fNPoints-1))) {
    fDirection = 1;          // decreasing x grid
    fLowerm    = fNPoints-1;
    fUpperm    = -1;
    if (fIsAcSecDerivs) {  // can be used only for increasing x data
      fIsAcSecDerivs = false;
    }
  } else {
    fDirection = 0;         // increasing x grid
    fLowerm    = -1;
    fUpperm    = fNPoints-1;
  }
  if (fIsConstrained) {
    SetParametersConstrained(); // constrained natural cubic spline
  } else {
    SetParameters();            // natural cubic spline
  }
}

// get interpollated Y value at the given X=val point
double Spline::GetValueAt(double val) {
  int    m,ml,mu,mav;
  // check if 'val' is above/below the highes/lowest value
  if (val>=*(fXdata+(fUpperm+fDirection))) {
    m = fUpperm+2*fDirection-1;
  } else if (val<=*(fXdata+(fLowerm+1-fDirection))) {
    m = fLowerm-2*fDirection+1;
  } else {  // Perform a binary search to find the interval val is in
    ml = fLowerm;
    mu = fUpperm;
    while (std::abs(mu-ml)>1) {
      mav = (ml+mu)/2.;
      if (val<*(fXdata+mav))
        mu = mav;
      else
        ml = mav;
      }
    m = mu+fDirection-1;
  }
  // check if the 2 grid point is 0,0
  if ((*(fYdata+m))+(*(fYdata+m+1))==0.0) {
    return 0.0;
  }
  if (!fIsConstrained)
    val -= (*(fXdata+m));

 return fA[m]+val*(fB[m]+val*(fC[m]+val*fD[m]));
}


// we need an other GetSplie such that if:
// - we know that x-values are increasing
// - and if we know which is the lower bin index j such that x_j < x <x_j+1
// - NOTE: before calling this, one needs to make sure that x_min < x < x_max

double Spline::GetValueAt(double x, int indx){
  // check if the 2 grid point is 0,0
  if ((*(fYdata+indx))+(*(fYdata+indx+1))==0.0) {
    return 0.0;
  }
  if (!fIsConstrained)
    x -= (*(fXdata+indx));
 return fA[indx] + x*(fB[indx] + x*( fC[indx] + x*fD[indx]));
}


// init the spline parameters
// compute the natural cubic spline paraneters a,b,c,d
//
void Spline::SetParameters(){
  if (fIsAcSecDerivs && fNPoints>4) {
    FillSecondDerivatives(); // and compute parameters a,b,c,d
  } else { // use less accurate second derivatives
    fIsAcSecDerivs = false;
    int    m1,m2,m,mr;
    double s,r;
    m1 = 1;
    m2 = fNPoints-1;
    s  = 0.0;
    for (m=0; m < m2; ++m) {
      fD[m] = (*(fXdata+(m+1))) - (*(fXdata+m));
      r     = ((*(fYdata+(m+1))) - (*(fYdata+m)))/fD[m];
      fC[m] = r - s;
      s     = r;
    }
    s     = 0.0;
    r     = 0.0;
    fC[0] = 0.0;
    fC[fNPoints-1] = 0.0;
    for (m=m1; m<m2; ++m) {
      fC[m] = fC[m] + r*fC[m-1];
      fB[m] = 2.0*(  (*(fXdata+(m-1))) - (*(fXdata+(m+1)))  ) - r*s;
      s     = fD[m];
      r     = s/fB[m];
    }
    mr = m2-1;
    for (m=m1; m<m2; ++m) {
      fC[mr] = (fD[mr]*fC[mr+1] - fC[mr])/fB[mr];
      mr     = mr - 1;
    }
    for (m=0; m<m2; ++m) {
     s     = fD[m];
     r     = fC[m+1] - fC[m];
     fD[m] = r/s;
     fC[m] = 3.0*fC[m];
     fB[m] = (  *(fYdata+(m+1))  -  *(fYdata+m)  )/s - (fC[m]+r)*s;
     fA[m] = *(fYdata+m);
    }
  }
}

void Spline::SetParametersConstrained(){
  double f1pxi    = 0.0;
  double f1pxim1  = 0.0;
  double f1ppxi   = 0.0;
  double f1ppxim1 = 0.0;
  for (int i=1;i<fNPoints;++i) {
    f1pxi = 0.0;
    if (i<fNPoints-1) {
      double dum0  = (*(fYdata +(i+1)))-(*(fYdata+i));
      double dum1  = (*(fYdata+i))-(*(fYdata+i-1));
      if (dum0*dum1>0.0) {
        f1pxi = 2.0/( ( (*(fXdata+i+1))-(*(fXdata+i)) )/dum0 + (( (*(fXdata+i))-(*(fXdata+i-1)))/dum1) );
      }
    } else { // last point
      f1pxi = 1.5*( (*(fYdata+fNPoints-1))-(*(fYdata+fNPoints-2)) )/( (*(fXdata+fNPoints-1))-(*(fXdata+fNPoints-2))) - 0.5*f1pxim1;
    }
    if (i==1) { // first point
      f1pxim1 = 1.5*( (*(fYdata+1))-(*(fYdata)) )/( (*(fXdata+1))-(*(fXdata)) ) - 0.5*f1pxi;
    }
    double xdel1 = (*(fXdata+i))-(*(fXdata+i-1));
    double ydel1 = (*(fYdata+i))-(*(fYdata+i-1));
    f1ppxim1     = -2.0*(f1pxi+2.0*f1pxim1)/xdel1 + 6.0*ydel1/(xdel1*xdel1);
    f1ppxi       =  2.0*(2.0*f1pxi+f1pxim1)/xdel1 - 6.0*ydel1/(xdel1*xdel1);
    double xi2   = (*(fXdata+i)) * (*(fXdata+i));
    double xim12 = (*(fXdata+i-1)) * (*(fXdata+i-1));
    fD[i-1] = (f1ppxi-f1ppxim1)/(6.0*xdel1);
    fC[i-1] = 0.5*( (*(fXdata+i))*f1ppxim1-(*(fXdata+i-1))*f1ppxi )/xdel1;
    fB[i-1] = (ydel1-fC[i-1]*(xi2-xim12) - fD[i-1]*( (*(fXdata+i))*xi2 - (*(fXdata+i-1))*xim12 ))/xdel1;
    fA[i-1] = (*(fYdata+i-1)) - fB[i-1]*(*(fXdata+i-1)) - fC[i-1]*xim12 - fD[i-1]*xim12*(*(fXdata+i-1));
    f1pxim1 = f1pxi;
  }
}


///
// Computation of second derivatives using "Not-a-knot" endpoint conditions
// B.I. Kvasov "Methods of shape-preserving spline approximation"
// World Scientific, 2000
void Spline::FillSecondDerivatives() { // will resize it to zero at the end because so far we use only a,b,c,d
  int     n   = fNPoints-1;
  double *u   = new double[n];
  double  p   = 0.0;
  double  sig = 0.0;

  u[1] = ((fYdata[2]-fYdata[1])/(fXdata[2]-fXdata[1]) - (fYdata[1]-fYdata[0])/(fXdata[1]-fXdata[0]));
  u[1] = 6.0*u[1]*(fXdata[2]-fXdata[1]) / ((fXdata[2]-fXdata[0])*(fXdata[2]-fXdata[0]));

  // Decomposition loop for tridiagonal algorithm. secondDerivatives[i]
  // and u[i] are used for temporary storage of the decomposed factors.
  std::vector<double> secondDerivatives;
  secondDerivatives.resize(fNPoints,0.0);
  secondDerivatives[1] = (2.0*fXdata[1]-fXdata[0]-fXdata[2]) / (2.0*fXdata[2]-fXdata[0]-fXdata[1]);
  for (int i=2; i<n-1; ++i) {
    sig = (fXdata[i]-fXdata[i-1]) / (fXdata[i+1]-fXdata[i-1]);
    p   = sig*secondDerivatives[i-1] + 2.0;
    secondDerivatives[i] = (sig - 1.0)/p;
    u[i] = (fYdata[i+1]-fYdata[i])/(fXdata[i+1]-fXdata[i]) - (fYdata[i]-fYdata[i-1])/(fXdata[i]-fXdata[i-1]);
    u[i] = (6.0*u[i]/(fXdata[i+1]-fXdata[i-1])) - sig*u[i-1]/p;
  }

  sig    = (fXdata[n-1]-fXdata[n-2]) / (fXdata[n]-fXdata[n-2]);
  p      = sig*secondDerivatives[n-3] + 2.0;
  u[n-1] = (fYdata[n]-fYdata[n-1])/(fXdata[n]-fXdata[n-1]) - (fYdata[n-1]-fYdata[n-2])/(fXdata[n-1]-fXdata[n-2]);
  u[n-1] = 6.0*sig*u[n-1]/(fXdata[n]-fXdata[n-2]) - (2.0*sig - 1.0)*u[n-2]/p;

  p      = (1.0+sig) + (2.0*sig-1.0)*secondDerivatives[n-2];
  secondDerivatives[n-1] = u[n-1]/p;

  // The back-substitution loop for the triagonal algorithm of solving
  // a linear system of equations.
  for (int k=n-2; k>1; --k) {
    secondDerivatives[k] *= (secondDerivatives[k+1] - u[k]*(fXdata[k+1]-fXdata[k-1])/(fXdata[k+1]-fXdata[k]));
  }
  secondDerivatives[n] = (secondDerivatives[n-1] - (1.0-sig)*secondDerivatives[n-2])/sig;
  sig = 1.0 - ((fXdata[2]-fXdata[1])/(fXdata[2]-fXdata[0]));
  secondDerivatives[1] *= (secondDerivatives[2] - u[1]/(1.0-sig));
  secondDerivatives[0]  = (secondDerivatives[1] - sig*secondDerivatives[2])/(1.0-sig);

  delete [] u;

  // compute a,b,c,d,
  // it would also be possible to use the secondDerivatives on the fly without these a,b,c,d parameters so we might
  // add this later as an option just to compute the secondDerivatives; see the SplineInterpolation method
  for (int i=0; i<n; ++i) {
    double dx = fXdata[i+1]-fXdata[i];
    fA[i] = fYdata[i];
    fB[i] = (fYdata[i+1]-fYdata[i])/dx - (2.0*secondDerivatives[i]+secondDerivatives[i+1])*dx/6.0;
    fC[i] = 0.5*secondDerivatives[i];
    fD[i] = (secondDerivatives[i+1]-secondDerivatives[i])/(6.0*dx);
  }
  secondDerivatives.clear();
}

/*
double Spline::SplineInterpolation(double val) {
  int    m,ml,mu,mav;
  // check if 'val' is above/below the highes/lowest value
  if (val>=*(fXdata+(fUpperm+fDirection))) {
    m = fUpperm+2*fDirection-1;
  } else if (val<=*(fXdata+(fLowerm+1-fDirection))) {
    m = fLowerm-2*fDirection+1;
  } else {  // Perform a binary search to find the interval val is in
    ml = fLowerm;
    mu = fUpperm;
    while (std::abs(mu-ml)>1) {
      mav = (ml+mu)/2.;
      if (val<*(fXdata+mav))
        mu = mav;
      else
        ml = mav;
      }
    m = mu+fDirection-1;
  }
 return  SplineInterpolation(val,m);
}

double Spline::SplineInterpolation(double e, int idx) {
 // Spline interpolation is used to get the value. Before this method
 // is called it is ensured that the energy is inside the bin
 // 0 < idx < numberOfNodes-1
 static const double onesixth = 1.0/6.0;

 // check bin value
 double x1 = fXdata[idx];
 double x2 = fXdata[idx + 1];
 double delta = x2 - x1;

 double a = (x2 - e)/delta;
 double b = (e - x1)/delta;

 // Final evaluation of cubic spline polynomial for return
 double y1 = fYdata[idx];
 double y2 = fYdata[idx + 1];

 double res = a*y1 + b*y2 +
       ( (a*a*a - a)*fSecondDerivatives[idx] +
         (b*b*b - b)*fSecondDerivatives[idx + 1] )*delta*delta*onesixth;

 return res;
}

*/

} // namespace geantphysics
