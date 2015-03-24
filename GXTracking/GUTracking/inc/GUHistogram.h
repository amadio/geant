#ifndef GUHistogram_H
#define GUHistogram_H 1

#include "TFile.h"
#include "TH1F.h"

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

private:
  TH1F*  ftime;
  TH1F*  fenergyPrimary;
  TH1F*  fenergyGam;
  TH1F*  fangleGam;
  TH1F*  fenergyElec;
  TH1F*  fangleElec;  
  
private:
  void BookHistograms( double maxEnergy );

private:
  TFile* fHistFile;
};

} // end namespace vecphys

#endif
