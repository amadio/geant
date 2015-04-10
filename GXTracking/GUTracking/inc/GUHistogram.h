#ifndef GUHistogram_H
#define GUHistogram_H 1

#ifdef VECPHYS_ROOT
#include "TFile.h"
#include "TH1F.h"
#endif

namespace vecphys {

class GUHistogram 
{
public:
 
  GUHistogram(std::string fileName, double maxEnergy );
  ~GUHistogram();

  void RecordTime( double );
  void RecordHistos( double EnPrimary,
                     double EnFinalGamma,
                     double angleGamma,    
                     double EnElectron,
                     double angleElectron);

#ifdef VECPHYS_ROOT
private:
  void BookHistograms( double maxEnergy );
  TFile* fHistFile;

private:
  TH1F*  ftime;
  TH1F*  fenergyPrimary;
  TH1F*  fenergyGam;
  TH1F*  fangleGam;
  TH1F*  fenergyElec;
  TH1F*  fangleElec;  
  
#endif
};

} // end namespace vecphys

#endif
