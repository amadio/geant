#ifndef GEANTV_ALIASTABLEALTERNATIVE_H
#define GEANTV_ALIASTABLEALTERNATIVE_H

#include <vector>
#include <cstdlib>
#include <cmath>

namespace geantphysics {

struct LinAliasCached {
  LinAliasCached(int num) { fLinAliasData.resize(num); }

  struct LinAliasData {
    long long int fAliasIndx;
    double fAliasW;
    double fX;
    double fXdelta;
    double fYdataDelta;
    double fXdivYdelta;
    double fYdata;
  };
  std::vector<LinAliasData> fLinAliasData;
};

struct RatinAliasDataTrans {
  RatinAliasDataTrans(size_t n) { fData.resize(n); }
  RatinAliasDataTrans(){};
  struct Data {
    long long int fAliasIndx;
    double fAliasW;
    double fParaA;
    double fParaB;
    double fCumulative;
    double fXdata;
  };
  std::vector<Data> fData;
};

class AliasTableAlternative {
public:
  static double SampleLinear(LinAliasCached &linAliasCached, int numdata, double r1, double r2)
  {
    LinAliasCached::LinAliasData *data    = linAliasCached.fLinAliasData.data();
    const double rest                     = r1 * (numdata - 1);
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

  static double SampleRatin(RatinAliasDataTrans &ratinData, int numdata, double r1, double r2, int above)
  {
    RatinAliasDataTrans::Data *data       = ratinData.fData.data();
    const double rest                     = r1 * (numdata - 1);
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
