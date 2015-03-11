#ifndef GUHistogram_H
#define GUHistogram_H 1

#include "TFile.h"
#include "TH1F.h"

namespace vecphys {

class GUHistogram 
{
public:
 
  GUHistogram(std::string fileName);
  ~GUHistogram();

  TH1F*  ftime;
  TH1F*  fenergy;
  TH1F*  fangle;

private:
  void BookHistograms();

private:
  TFile* fHistFile;
};

} // end namespace vecphys

#endif
