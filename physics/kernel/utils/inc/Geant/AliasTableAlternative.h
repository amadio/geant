#ifndef GEANTV_ALIASTABLEALTERNATIVE_H
#define GEANTV_ALIASTABLEALTERNATIVE_H

#include <vector>
#include <cstdlib>
#include <cmath>

namespace geantphysics {

/** @brief This is table with data needed for linear alias sampling.
 *  The difference is that this one:
 *  1) Stores data in form of SOA(only one memory access is required to get all the data)
 *  2) Stores some additional values for simplified calculations(there is memory overhead)
 */
struct LinAliasCached {
  LinAliasCached(int num) : fNumData(num) { fLinAliasData.resize(num); }

  struct LinAliasData {
    long long int fAliasIndx;
    double fAliasW;
    double fX;
    double fXdelta;
    double fYdataDelta;
    double fXdivYdelta;
    double fYdata;
  };
  int fNumData;
  std::vector<LinAliasData> fLinAliasData;

  /**
   * @param r1 random number [0,1)
   * @param r2 random number [0,1)
   * @return sampled value
   */
  double Sample(double r1, double r2)
  {
    LinAliasCached::LinAliasData *data    = fLinAliasData.data();
    const double rest                     = r1 * (fNumData - 1);
    int indxl                             = (int)(rest);
    const double dum0                     = rest - indxl;
    if (data[indxl].fAliasW < dum0) indxl = (int)data[indxl].fAliasIndx;

    if (data[indxl].fYdata > 0) {
      const double dum = data[indxl].fYdataDelta;
      if (std::abs(dum) > 0.1)
        return data[indxl].fX - data[indxl].fXdivYdelta * (1.0 - std::sqrt(1.0 + r2 * dum * (dum + 2.0)));
      else // use second order Taylor around dum = 0.0
        return data[indxl].fX + r2 * data[indxl].fXdelta * (1.0 - 0.5 * dum * (r2 - 1.0) * (1.0 + dum * r2));
    }
    return data[indxl].fX + data[indxl].fXdelta * std::sqrt(r2);
  }
};

/** @brief This is table with data needed for ratin alias sampling.
 *  The difference is that this one:
 *  1) Stores data in form of SOA(only one memory access is required to get all the data)
 */
struct RatinAliasDataTrans {
  RatinAliasDataTrans(int n) : fNumdata(n) { fData.resize(n); }
  RatinAliasDataTrans(){};
  struct Data {
    long long int fAliasIndx;
    double fAliasW;
    double fParaA;
    double fParaB;
    double fCumulative;
    double fXdata;
  };
  int fNumdata;
  std::vector<Data> fData;

  /**
   * @param r1 random number [0,1)
   * @param r2 random number [0,1)
   * @param above sample with linear for bins > above
   * @return
   */
  double Sample(double r1, double r2, int above)
  {
    RatinAliasDataTrans::Data *data       = fData.data();
    const double rest                     = r1 * (fNumdata - 1);
    int indxl                             = (int)(rest);
    const double dum0                     = rest - indxl;
    if (data[indxl].fAliasW < dum0) indxl = (int)data[indxl].fAliasIndx;

    double res = 0.0;
    if (indxl > above) {
      // sample value within the selected bin by using ratin based numerical inversion
      const double delta = data[indxl + 1].fCumulative - data[indxl].fCumulative;
      const double aval  = r2 * delta;
      const double dum1  = (1.0 + data[indxl].fParaA + data[indxl].fParaB) * delta * aval;
      const double dum2  = delta * delta + data[indxl].fParaA * delta * aval + data[indxl].fParaB * aval * aval;
      res                = data[indxl].fXdata + dum1 / dum2 * (data[indxl + 1].fXdata - data[indxl].fXdata);
    } else { // use linear aprx by assuming that pdf[indxl]=0.0
      res = data[indxl].fXdata + (data[indxl + 1].fXdata - data[indxl].fXdata) * std::sqrt(r2);
    }
    return res;
  }
};
}
#endif // GEANTV_ALIASTABLEALTERNATIVE_H
