#include "TVNudyModel.h"

ClassImp(TVNudyModel)

    //______________________________________________________________________________
    void TVNudyModel::ReadFile3(TNudyEndfFile *file) {
  TIter secIter(file->GetSections());
  TNudyEndfSec *sec;
  while ((sec = (TNudyEndfSec *)secIter.Next())) {
    if (sec->GetMT() == (int)fReaction) {
      // Read the data required for GetXSect(E)
      TNudyEndfTab1 *record = (TNudyEndfTab1 *)(sec->GetRecords()->At(0));
      fEXSect_length = record->GetNP();
      fE_file3 = new double[fEXSect_length]();
      fXSect_file3 = new double[fEXSect_length]();
      for (int i = 0; i < fEXSect_length; i++) {
        fE_file3[i] = record->GetX(i);
        fXSect_file3[i] = record->GetY(i);
      }
    }
  }
}

//______________________________________________________________________________
double TVNudyModel::GetXSect(double e) {
  if (!fE_file3 || !fXSect_file3 || fEXSect_length == 0) {
    Error("GetXSect", "Energy-Cross Section data for model is not set\n");
    return 0;
  }
  if (e <= fE_file3[0])
    return fXSect_file3[0];
  if (e >= fE_file3[fEXSect_length - 1])
    return fXSect_file3[fEXSect_length - 1];
  int low, high, mid;
  low = 0;
  high = fEXSect_length;
  while (high - low > 1) {
    mid = (low + high) / 2;
    if (fE_file3[mid] >= e)
      high = mid;
    else
      low = mid;
  }
  if (high == low)
    return fXSect_file3[low];
  else
    return TNudyCore::Instance()->LinearInterpolation(fE_file3[low], fXSect_file3[low], fE_file3[high],
                                                      fXSect_file3[high], e);
}
