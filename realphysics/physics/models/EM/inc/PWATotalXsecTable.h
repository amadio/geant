
#ifndef PWATOTALXSECTABLE_H
#define PWATOTALXSECTABLE_H

namespace geantphysics {
//
//  PWATotalXsecZ: sub-class for PWA xsec data that belong to a given Z number
class PWATotalXsecZ {
friend class PWATotalXsecTable;
public:
  //
  // out of energy grid cases
  int    GetLowestEnergyBinIndex() const { return 0; }
  int    GetHighestEnergyBinIndex()const { return gNumTotalXsecBins-1; }
  static double GetLowestEnergy()  { return gPWATotalXsecEnergyGrid[0]; }
  static double GetHighestEnergy() { return gPWATotalXsecEnergyGrid[gNumTotalXsecBins-1]; }

  // see below what is input parameter j
  double GetLowestXsecValue(int j)  const { return fPWAXsecs[j*gNumTotalXsecBins];}
  double GetHighestXsecValue(int j) const { return fPWAXsecs[(j+1)*gNumTotalXsecBins-1];}

  //
  // normal cases i.e. energy is within the grid
  // kinetic energy in GeV ; returns with the index of the lower energy bin edge
  int GetPWATotalXsecEnergyBinIndex(double energy) const;

  //------------------------------------------------------------------------------//
  // The GetPWATotalXsecEnergyBinIndex(energy) will return with the lower energy  //
  // bin edge index = elowindx. Then the following formulas can be used to get the//
  // elastic,  first and second transport mean free path lower bin edge values:   //
  //  index of the lower energy bin edge = j*fgNumTotalXsecBins + elowindex       //
  // where j is                                                                   //
  //  -elastic cross section lower bin edge index:           j = 1.5 + chrage*1.5 //
  //  -first transport cross section lower energy bin index: j = 2.5 + chrage*1.5 //
  //  -first transport cross section lower energy bin index: j = 3.5 + chrage*1.5 //
  // With this, we can avoid to use an IF over particle types (e-/e+)             //
  // Additional note: it's probably a good idea to separate the elowindex comp-   //
  // utation because it depends only on the energy of the particle while the      //
  // cross sections depends on Z and particle type as well                        //
  //------------------------------------------------------------------------------//
  double GetInterpXsec(double energy, int elowindex, int j) const ;
  double GetInterpXsec(double energy, int j) const ;

private:
  // ctr and dtr can be called only by the PWATotalXsecTable friend
  PWATotalXsecZ(int Z);
 ~PWATotalXsecZ() {}

  //  hide assignment operator and cpy ctr.
  PWATotalXsecZ & operator=(const PWATotalXsecZ &right);
  PWATotalXsecZ(const PWATotalXsecZ&);

  void LoadPWATotalXsecZ(int Z);


private:
  //size of the common energy grid //
  static const int gNumTotalXsecBins =  106;

  // common energy grid in [1.e-4;1.e+3] MeV //
  // size is fgNumTotalXsecBins
  static const double gPWATotalXsecEnergyGrid[gNumTotalXsecBins];

  // elastic cross sections, first and second transport cross sections for e-/e+
  // over the common energy grid gPWATotalXsecEnergyGrid in geantV internal length^2
  double fPWAXsecs[gNumTotalXsecBins*6];
  // interpolation parameters if log-log linear interpolation is used
  double fInterpParamA[gNumTotalXsecBins*6];
  double fInterpParamB[gNumTotalXsecBins*6];
};


//
//  PWATotalXsecTable
class PWATotalXsecTable {
public:
  PWATotalXsecTable() {}
 ~PWATotalXsecTable();

  void Initialise();

  const PWATotalXsecZ* GetPWATotalXsecForZet(int Z) const{
        Z = Z>gNumZet ? gNumZet : Z;
        return gPWATotalXsecTable[Z-1];
  }

private:
  //  hide assignment operator and cpy ctr.
  PWATotalXsecTable & operator=(const PWATotalXsecTable &right);
  PWATotalXsecTable(const PWATotalXsecTable&);

 private:
   // size of the table: Z=1-103 //
   static const int gNumZet =  103;

   // PWATotalXsecZ pointers for Z=1-103 //
   static PWATotalXsecZ *gPWATotalXsecTable[gNumZet];
};

}       // namespace geantphysics

#endif  // PWATOTALXSECTABLE_H
