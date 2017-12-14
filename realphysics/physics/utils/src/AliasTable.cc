
#include "AliasTable.h"

#include "Spline.h"
#include "GLIntegral.h"

// cstdlib is for std::abs(int)
#include <cstdlib>
#include <cmath>

namespace geantphysics {

/**
  * Alias table is built based on the discretized continuous distribution provided as input at a set
  * \f$\{x_i\}_{i=0}^{N-1}\f$ of discrete values of the random variable \f$x\f$. During the sampling, this
  * alias table is used to determine the bin \f$i\f$ such that \f$\mathcal{P}(\xi_i)\leq\xi<\mathcal{P}(\xi_{i+1})\f$
  * where \f$\mathcal{P}(x)\f$ is the cumulative distribution function and \f$ \xi \in \mathbb{U}[0,1]\f$ random
  * number. When numerical inversion of the c.d.f. is used, one needs to solve \f$\mathcal{P}^{-1}(\xi) = x\f$ i.e.
  * interpolation within
  * \f$\mathcal{P}^{-1}(\xi_i) = x_i  \leq \mathcal{P}^{-1}(\xi) = x < \mathcal{P}^{-1}(\xi_{i+1}) = x_{i+1}\f$. The
  * applied interpolation should satisfy \f$\frac{\mathrm{d}\mathcal{P}^{-1}(\xi)}{\mathrm{d}\xi}
  *   = \left( \frac{\mathrm{d}\mathcal{P}(x)}{\mathrm{d}x} \right)^{-1}
  *   =\frac{1}{p(x)}\f$ and \f$\mathcal{P}^{-1}(\xi_i)=x_i\f$, \f$\mathcal{P}^{-1}(\xi_{i+1})=x_{i+1}\f$.
  *
  * When rational function approximation is used in the form of \f$x=\mathcal{P}^{-1}(\xi)
  * \approx \tilde{\mathcal{P}}^{-1}(\xi) =x_i + \frac{(1+a_i+b_i)\alpha}{1+a_i\alpha+b_i\alpha^{2}}[x_{i+1}-x_i] \f$,
  * where \f$\alpha=\frac{\xi-\xi_i}{\xi_{i+1}-\xi_i}\f$:
  *
  *  - \f$\tilde{\mathcal{P}}^{-1}(\xi_i)=x_i\f$ and \f$\tilde{\mathcal{P}}^{-1}(\xi_{i+1})=x_{i+1}\f$
  *        independently form the values \f$a_i,b_i\f$
  *  - \f$\frac{\mathrm{d}\tilde{\mathcal{P}}^{-1}(\xi)}{\mathrm{d}\xi}
  *        =\frac{(1+a_i+b_i)(1-b_i\alpha^{2})}{[1+a_i\alpha+b_i\alpha^{2}]^{2}}\frac{x_{i+1}-x_i}{\xi_{i+1}-\xi_i}\f$
  *  - and the parameters \f$a_i,b_i\f$ can be determined from the requirements:
  *     1. \f$\frac{\mathrm{d}\tilde{\mathcal{P}}^{-1}(\xi)}{\mathrm{d}\xi}_{|\xi=\xi_i}=\frac{1}{p(x_i)}\f$
  *     2. \f$\frac{\mathrm{d}\tilde{\mathcal{P}}^{-1}(\xi)}{\mathrm{d}\xi}_{|\xi=\xi_{i+1}}=\frac{1}{p(x_{i+1})}\f$
  *  - that yields
  *     1. \f$b_i=1-\left[\frac{\xi_{i+1}-\xi_i}{x_{i+1}-x_i}\right]^{2}\frac{1}{p(x_i)p(x_{i+1})}\f$
  *     2. \f$a_i=\frac{\xi_{i+1}-\xi_i}{x_{i+1}-x_i}\frac{1}{p(x_i)}-1-b_i\f$
  *  - then the sampled value \f$x \approx \tilde{x} = \tilde{\mathcal{P}}^{-1}(\xi)
  *       = x_i + \frac{(1+a_i+b_i)\alpha}{1+a_i\alpha+b_i\alpha^{2}}[x_{i+1}-x_i]\f$,
  *       with \f$\alpha=\frac{\xi-\xi_i}{\xi_{i+1}-\xi_i}\f$
  *
  * The cumulative \f$\mathcal{P}(x_i)=\xi_i\f$, parameters \f$a_i,b_i\f$ are computed in this method at the given
  * discrete values \f$\{x_i\}_{i=0}^{N-1}\f$ of the random variable \f$x\f$. The AliasTable::SampleRatin method can
  * be used then to sample \f$x\f$ based on the procedure described above.
  */
void AliasTable::PreparRatinTable(double *xdata, double *ydata, double *comf, double *paradata, double *parbdata,
                                  double *xx, int *binindx, int numdata, bool isconstraintspline, int glnum) {
  // compute (integral) probability within each bin; compute c.d.f.
  GLIntegral *gl = new GLIntegral(glnum, 0.0, 1.0);
  const std::vector<double> &glx = gl->GetAbscissas();
  const std::vector<double> &glw = gl->GetWeights();
  Spline *sp = new Spline(xdata, ydata, numdata, true, isconstraintspline); // use cubic spline
  double sum = 0.0;
  comf[0]    = 0.0;
  for (int i=1; i<numdata; ++i) {
    double dum0 = 0.0;
    for (int j=0; j<glnum; ++j) {
      double xval = glx[j]*(xdata[i]-xdata[i-1])+xdata[i-1];
      dum0 += glw[j]*(xdata[i]-xdata[i-1])*sp->GetValueAt(xval);
    }
    if (dum0<1.e-40)
      dum0 = 1.e-40;
    xx[i-1]      = -dum0;
    binindx[i-1] = i-1;
    sum       += dum0;
    comf[i]    = sum;
  }
  // prepare alias table
  sum /= (double)(numdata-1);
  for (int i=0; i<numdata-1; ++i) {
    int indxh = 0;
    for (indxh = 0; indxh<numdata-1; ++indxh) {
      if (xx[indxh]<0.0 && std::fabs(xx[indxh])>sum)
        break;
    }
    int indxl = 0;
    for (indxl = 0; indxl<numdata-1; ++indxl) {
      if (xx[indxl]<0.0 && std::fabs(xx[indxl])<sum)
        break;
    }
    double dum0 = sum - std::fabs(xx[indxl]);
    xx[indxh] +=  dum0;
    xx[indxl] = -1.0*xx[indxl]/sum;
    binindx[indxl] = indxh;
    if (i==numdata-2) {
      xx[indxh] = 1.0;
      xx[indxl] = 1.0;
    }
  }
  // ensure normality of the p.d.f. (and c.d.f.)
  double norm = 1.0/comf[numdata-1];
  for (int j=0; j<numdata; ++j) {
    comf[j]   *= norm;
    ydata[j]  *= norm;
    if (ydata[j]<1.e-40)
      ydata[j] = 1.e-40;
  }
  // compute parameters a,b of ratin
  for (int j=0; j<numdata-1; ++j) {
    double dum1 = (comf[j+1]-comf[j])/(xdata[j+1]-xdata[j]);
    double parb = 1.0-dum1*dum1/(ydata[j]*ydata[j+1]);
    double para = dum1/ydata[j]-1.0-parb;
    paradata[j] = para;
    parbdata[j] = parb;
  }
  delete gl;
  delete sp;
}

void AliasTable::PreparLinearTable(double *xdata, double *ydata, double *xx, int *binindx, int numdata) {
  // compute (integral) probability within each bin; compute c.d.f.
  double sum = 0.0;
  for (int i=1; i<numdata; ++i) {
    double dum0 = 0.5*(xdata[i]-xdata[i-1])*(ydata[i]+ydata[i-1]);
    if (dum0<1.e-40)
     dum0 = 1.e-40;
    xx[i-1]      = -dum0;
    binindx[i-1] = i-1;
    sum += dum0;
  }
  // prepare alias table
  sum /= (double)(numdata-1);
  for (int i=0; i<numdata-1; ++i) {
    int indxh = 0;
    for (indxh = 0; indxh<numdata-1; ++indxh) {
      if (xx[indxh]<0.0 && std::fabs(xx[indxh])>sum)
        break;
    }
    int indxl = 0;
    for (indxl = 0; indxl<numdata-1; ++indxl) {
      if (xx[indxl]<0.0 && std::fabs(xx[indxl])<sum)
        break;
    }
    double dum0 = sum - std::fabs(xx[indxl]);
    xx[indxh] +=  dum0;
    xx[indxl] = -1.0*xx[indxl]/sum;
    binindx[indxl] = indxh;
    if (i==numdata-2) {
      xx[indxh] = 1.0;
      xx[indxl] = 1.0;
    }
  }
}

double AliasTable::SampleRatin(double *xdata, double *comf, double *paradata, double *parbdata, double *xx, int *binindx,
                               int numdata, double rndm1, double rndm2, int above) {
  // get the lower index of the bin by using the alias part
  const double rest = rndm1*(numdata-1);
  int indxl         = (int) (rest);
  const double dum0 = rest-indxl;
  if (xx[indxl]<dum0)
    indxl = binindx[indxl];

  double res = 0.0;
  if (indxl>above) {
    // sample value within the selected bin by using ratin based numerical inversion
    const double  delta = comf[indxl+1]-comf[indxl];
    const double  aval  =  rndm2 * delta;
    const double  dum1  = (1.0+paradata[indxl]+parbdata[indxl])*delta*aval;
    const double  dum2  = delta*delta+paradata[indxl]*delta*aval+parbdata[indxl]*aval*aval;
    res = xdata[indxl] +  dum1/dum2 *(xdata[indxl+1]-xdata[indxl]);
  } else { // use linear aprx by assuming that pdf[indxl]=0.0
    res = xdata[indxl] +  (xdata[indxl+1]-xdata[indxl])*std::sqrt(rndm2);
  }
  return res;
}


double AliasTable::SampleRatin(const double *xdata, const double *comf, const double *paradata, const double *parbdata,
                               const double *xx, const int *binindx, const int numdata, const double rndm1,
                               const double rndm2, const int above) {
  // get the lower index of the bin by using the alias part
  const double rest = rndm1*(numdata-1);
  int indxl         = (int) (rest);
  const double dum0 = rest-indxl;
  if (xx[indxl]<dum0)
    indxl = binindx[indxl];

  double res = 0.0;
  if (indxl>above) {
    // sample value within the selected bin by using ratin based numerical inversion
    const double  delta = comf[indxl+1]-comf[indxl];
    const double  aval  =  rndm2 * delta;
    const double  dum1  = (1.0+paradata[indxl]+parbdata[indxl])*delta*aval;
    const double  dum2  = delta*delta+paradata[indxl]*delta*aval+parbdata[indxl]*aval*aval;
    res = xdata[indxl] +  dum1/dum2 *(xdata[indxl+1]-xdata[indxl]);
  } else { // use linear aprx by assuming that pdf[indxl]=0.0
    res = xdata[indxl] +  (xdata[indxl+1]-xdata[indxl])*std::sqrt(rndm2);
  }
  return res;
}


double AliasTable::SampleLinear(double *xdata, double *ydata, double *xx, int *binindx, int numdata, double rndm1,
                                double rndm2) {
  // get the lower index of the bin by using the alias part
  const double rest = rndm1*(numdata-1);
  int indxl         = (int) (rest);
  const double dum0 = rest-indxl;
  if (xx[indxl]<dum0)
    indxl = binindx[indxl];
  // sample value within the selected bin by using linear approximation
  const double xval   = xdata[indxl];
  const double xdelta = xdata[indxl+1]-xval;
  if (ydata[indxl] > 0.0) {
    const double dum = (ydata[indxl+1]-ydata[indxl])/ydata[indxl];
    if (std::abs(dum)>0.1)
      return xval - xdelta/dum * (1.0 - std::sqrt(1.0+rndm2*dum*(dum+2.0)));
    else // use second order Taylor around dum = 0.0
      return xval + rndm2*xdelta*(1.0-0.5*dum*(rndm2-1.0)*(1.0+dum*rndm2));
  }
  return xval + xdelta*std::sqrt(rndm2);
}

double AliasTable::SampleLinear(const double *xdata, const double *ydata, const double *xx, const int *binindx,
                                const int numdata, const double rndm1, const double rndm2) {
  // get the lower index of the bin by using the alias part
  const double rest = rndm1*(numdata-1);
  int indxl         = (int) (rest);
  const double dum0 = rest-indxl;
  if (xx[indxl]<dum0)
    indxl = binindx[indxl];
  // sample value within the selected bin by using linear approximation
  const double xval   = xdata[indxl];
  const double xdelta = xdata[indxl+1]-xval;
  if (ydata[indxl] > 0.0) {
    const double dum = (ydata[indxl+1]-ydata[indxl])/ydata[indxl];
    if (std::abs(dum)>0.1)
      return xval - xdelta/dum * (1.0 - std::sqrt(1.0+rndm2*dum*(dum+2.0)));
    else // use second order Taylor around dum = 0.0
      return xval + rndm2*xdelta*(1.0-0.5*dum*(rndm2-1.0)*(1.0+dum*rndm2));
  }
  return xval + xdelta*std::sqrt(rndm2);
}

/**
  * When the rational interpolation based numerical inversion of the c.d.f. is used as described at
  * AliasTable::PreparRatinTable, the approximated p.d.f. at the random variable value \f$x_i\leq x<x_{i+1}\f$ can be
  * written as \cite salvat2006penelope
  * \f[
  *  p(x)\approx\tilde{p}(x)=\left[\frac{\mathrm{d}\tilde{\mathcal{P}}^{-1}(\xi)}{\mathrm{d}\xi}\right]^{-1}
  *   =\frac{(1+a_i\alpha+b_i\alpha^2)^2}{(1+a_i+b_i)(1-b_i\alpha^2)}\frac{\xi_{i+1}-\xi_i}{x_{i+1}-x_i}
  * \f]
  * where \f$\mathcal{P(x_i)}=\xi_i,\;\mathcal{P(x_{i+1})}=\xi_{i+1}\f$ and the parameters \f$a_i,b_i\f$ are the same
  * as described at AliasTable::PreparRatinTable. In order to compute the approximated p.d.f. one needs to determine
  * the value of \f$alpha\f$. This can be done by using the sampling formula given at AliasTable::PreparRatinTable
  * and the roots that satisfy the conditions \f$ x=x_i \to \alpha=0\f$ and \f$x=x_{i+1} \to \alpha=1\f$
  * \cite salvat2006penelope
  * \f[
  *  \alpha = \frac{1+a_i+b_i-a_i\beta}{2b_i\beta}\left[1- \sqrt{1-\frac{4b_i\beta^2}{(1+a_i+b_i-a_i\beta)^2}}\right]
  * \f]
  * where \f$\beta \equiv (x-x_i)/(x_{i+1}-x_i)\f$.
  *
  * The parameters of the rational interpolation are computed in this method and can be used in
  * AliasTable::GetRatinForPDF or AliasTable::GetRatinForPDF1 to get the approximated p.d.f. value
  * \f$p(x)\approx\tilde{p}(x)\f$ at a given random variable value \f$x\f$. This can be used to compute the
  * approximation error at a given \f$x_j,x_{j+1}\f$ interval by computing \cite salvat2006penelope
  * \f[
  *   \varepsilon_j = \int_{x_j}^{x_{j+1}} \mid p(x)-\tilde{p}(x) \mid \mathrm{d}x
  * \f]
  * Knowing these approximation error values for all intervals, one can adaptively insert new discrete sample points
  * of the random variable into the interval with the highest approximation error gradualay decreasing the overall
  * approximation error.
  */
double AliasTable::PreparRatinForPDF(double *xdata, double *ydata, double *comf, double *paradata, double *parbdata,
                                     int numdata, bool isconstraintspline, int glnum) {
  // compute (integral) probability within each bin; compute c.d.f.
  GLIntegral *gl = new GLIntegral(glnum, 0.0, 1.0);
  const std::vector<double> &glx = gl->GetAbscissas();
  const std::vector<double> &glw = gl->GetWeights();
  Spline     *sp = new Spline(xdata, ydata, numdata, true, isconstraintspline);
  double sum   = 0.0;
  comf[0]      = 0.0;
  for (int i=1; i<numdata; ++i) {
    double dum0 = 0.0;
    for (int j=0; j<glnum; ++j) {
      const double xval = glx[j]*(xdata[i]-xdata[i-1])+xdata[i-1];
      dum0 += glw[j]*(xdata[i]-xdata[i-1])*sp->GetValueAt(xval);
    }
    sum     += dum0;
    comf[i]  = sum;
  }
  // ensure normality of the p.d.f. (and c.d.f.)
  const double norm = 1.0/comf[numdata-1];
  for (int j=0; j<numdata; ++j) {
    comf[j]   *= norm;
    ydata[j]  *= norm;
    if (ydata[j]<1.e-40)
      ydata[j] = 1.e-40;
  }
  // compute parameters a,b of ratin
  for (int j=0; j<numdata-1; ++j) {
    const double dum1 = (comf[j+1]-comf[j])/(xdata[j+1]-xdata[j]);
    const double parb = 1.0-dum1*dum1/(ydata[j]*ydata[j+1]);
    const double para = dum1/ydata[j]-1.0-parb;
    paradata[j] = para;
    parbdata[j] = parb;
  }
  delete gl;
  delete sp;
  return norm;
}

double AliasTable::GetRatinForPDF(double x, double *xdata, double *comf, double *paradata, double *parbdata, int numdata) {
  // find indx of lower bin of xdata where x is located
  const int lindx = BSearch(x, xdata, numdata);
  // compute the approximated p.d.f. value
  const double lx  = xdata[lindx];
  const double la  = paradata[lindx];
  const double lb  = parbdata[lindx];
  const double invDeltX = 1.0/(xdata[lindx+1]-xdata[lindx]);
  const double deltCum  = comf[lindx+1]-comf[lindx];
  const double tau  = (x-lx)*invDeltX;
  const double dum0 = 1.0+la+lb-la*tau;
  double eta  = dum0/(2.0*lb*tau)*(1.0-std::sqrt(1.0-4.0*lb*tau*tau/(dum0*dum0)));
  if (x==lx)
    eta = 0.0;
  const double dum1  = 1.0+la*eta+lb*eta*eta;
  return dum1*dum1*deltCum*invDeltX/((1.0+la+lb)*(1.0-lb*eta*eta));
}

double AliasTable::GetRatinForPDF1(double x, double *xdata, double *comf, double *paradata, double *parbdata, int lindx) {
  const double lx  = xdata[lindx];
  const double la  = paradata[lindx];
  const double lb  = parbdata[lindx];
  const double invDeltX = 1.0/(xdata[lindx+1]-xdata[lindx]);
  const double deltCum  = comf[lindx+1]-comf[lindx];
  const double tau  = (x-lx)*invDeltX;
  const double dum0 = 1.0+la+lb-la*tau;
  double eta  = dum0/(2.0*lb*tau)*(1.0-std::sqrt(1.0-4.0*lb*tau*tau/(dum0*dum0)));
  if (x==lx)
    eta = 0.0;
  const double dum1 = 1.0+la*eta+lb*eta*eta;
  return dum1*dum1*deltCum*invDeltX/((1.0+la+lb)*(1.0-lb*eta*eta));
}

double AliasTable::BSearch(double val, double *xdata, int npoints) {
  int lowerm,upperm,direction,m,ml,mu,mav;
  if (xdata[0]>xdata[npoints-1]) {
    direction = 1;
    lowerm    = npoints-1;
    upperm    = -1;
  } else {
    direction = 0;
    lowerm    = -1;
    upperm    = npoints-1;
  }
  // check if 'val' is above/below the highes/lowest value
  if (val>=xdata[upperm+direction]) {
    m = upperm+2*direction-1;
  } else if (val<=xdata[lowerm+1-direction]) {
    m = lowerm-2*direction+1;
  } else {  // Perform a binary search to find the interval val is in
    ml = lowerm;
    mu = upperm;
    while (std::abs(mu-ml)>1) {
      mav = (ml+mu)/2.;
      if (val<xdata[mav])
        mu = mav;
      else
        ml = mav;
      }
    m = mu+direction-1;
  }
 return m;
}


} // namespace geantphysics
