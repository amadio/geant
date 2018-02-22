
#include "NISTElementData.h"

#include "PhysicalConstants.h"

#include <cmath>


namespace geantphysics {
  
void NISTElementData::BuildTable() {
   using geant::units::eV;
   using geant::units::GeV;
   using geant::units::kElectronMassC2;
   using geant::units::kProtonMassC2;
   using geant::units::kAvogadro;
   using geant::units::kCLightSquare;
   using geant::units::kAtomicMassUnit;
   using geant::units::kAtomicMassUnitC2;

   constexpr double u   = kAtomicMassUnit;
   constexpr double kDeuteronMassC2   = 1.875613*GeV; // TODO: update this according to particle properties
   constexpr double kTritonMassC2     = 2.808921*GeV; // TODO: update this according to particle properties
   constexpr double kHe3MassC2        = 2.808391*GeV; // TODO: update this according to particle properties
   constexpr double kAlphaMassC2      = 3.727379*GeV; // TODO: update this according to particle properties

   // Parameterisation of the total binding energy of the electrons from
   // D.Lunney, J.M.Pearson,C .Thibault,
   // Rev.Mod.Phys. 75 (2003) 1021 (APPENDIX A. Eq. (A4))
   for (int i = 0; i < gNumberOfNISTElements; ++i) {
     double Z = (double)(i+1);
     fBindingEnergies[i] = (14.4381*std::pow(Z,2.39) + 1.55468e-6*std::pow(Z,5.35))*eV;
   }

   const double kBindingEnergyH  =  (14.4381 + 1.55468e-6)*eV;
   // (2.0^2.39) = 5.24157361543345  (2.0^5.35) = 40.7859400742164
   const double kBindingEnergyHe =  (14.4381*5.2415736154 + 1.55468e-6*40.7859400742164)*eV;



  // Z = 1-------------------------------------------------------------------------------
   static const std::string   HS[] = {    "H1",     "H2",     "H3",     "H4",     "H5",
                                          "H6",     "H7" };
   static const int      HN[] = {  1,   2,   3,   4,   5,   6,   7 };
   // Garantee consistence with GV Proton, Deuteron and Triton  masses
   static const double   HA[] = {(kProtonMassC2   + kElectronMassC2 - kBindingEnergyH)/kCLightSquare,
                                 (kDeuteronMassC2 + kElectronMassC2 - kBindingEnergyH)/kCLightSquare,
                                 (kTritonMassC2   + kElectronMassC2 - kBindingEnergyH)/kCLightSquare,
                                         4.02643 * u,        5.035311 * u,
                                         6.04496 * u,          7.0527 * u };
   static const double   HW[] = { 0.999885, 0.000115, 0., 0., 0., 0., 0. };
   static double HIsoMass[7];

   fNISTElementDataTable[0].fElemSymbol    = "H";
   fNISTElementDataTable[0].fSymbols       = HS;
   fNISTElementDataTable[0].fZ             = 1;
   fNISTElementDataTable[0].fNumOfIsotopes = 7;
   fNISTElementDataTable[0].fIsNoStableIsotope       = false;
   fNISTElementDataTable[0].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[0].fNIsos         = HN;
   fNISTElementDataTable[0].fAIsos         = HA;
   fNISTElementDataTable[0].fWIsos         = HW;
   fNISTElementDataTable[0].fMassIsos      = HIsoMass;


  // Z = 2-------------------------------------------------------------------------------
   static const std::string  HeS[] = {   "He3",    "He4",    "He5",    "He6",    "He7",
                                         "He8",    "He9",   "He10" };
   static const int     HeN[] = {  3,   4,   5,   6,   7,   8,   9,  10 };
   // Garantee consistence with GV He3 and Alpha  masses
   static const double  HeA[] = {(kHe3MassC2   + 2.0*kElectronMassC2 - kBindingEnergyHe)/kCLightSquare,
                                 (kAlphaMassC2 + 2.0*kElectronMassC2 - kBindingEnergyHe)/kCLightSquare,
                                        5.012057 * u,     6.018885891 * u,
                                       7.0279907 * u,      8.03393439 * u,
                                        9.043946 * u,        10.05279 * u };
   static const double  HeW[] = { 1.34e-06, 0.99999866, 0., 0., 0., 0., 0., 0. };
   static double HeIsoMass[8];

   fNISTElementDataTable[1].fElemSymbol    = "He";
   fNISTElementDataTable[1].fSymbols       = HeS;
   fNISTElementDataTable[1].fZ             = 2;
   fNISTElementDataTable[1].fNumOfIsotopes = 8;
   fNISTElementDataTable[1].fIsNoStableIsotope       = false;
   fNISTElementDataTable[1].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[1].fNIsos         = HeN;
   fNISTElementDataTable[1].fAIsos         = HeA;
   fNISTElementDataTable[1].fWIsos         = HeW;
   fNISTElementDataTable[1].fMassIsos      = HeIsoMass;


  // Z = 3-------------------------------------------------------------------------------
   static const std::string  LiS[] = {   "Li3",    "Li4",    "Li5",    "Li6",    "Li7",
                                         "Li8",    "Li9",   "Li10",   "Li11",   "Li12",
                                        "Li13" };
   static const int     LiN[] = {  3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
                                  13 };
   static const double  LiA[] = {         3.0308 * u,         4.02719 * u,
                                        5.012538 * u,    6.0151228874 * u,
                                    7.0160034366 * u,     8.022486246 * u,
                                      9.02679019 * u,       10.035483 * u,
                                     11.04372358 * u,       12.052517 * u,
                                        13.06263 * u };
   static const double  LiW[] = { 0., 0., 0., 0.0759, 0.9241, 0., 0., 0., 0., 0.,
                                  0. };
   static double LiIsoMass[11];

   fNISTElementDataTable[2].fElemSymbol    = "Li";
   fNISTElementDataTable[2].fSymbols       = LiS;
   fNISTElementDataTable[2].fZ             = 3;
   fNISTElementDataTable[2].fNumOfIsotopes = 11;
   fNISTElementDataTable[2].fIsNoStableIsotope       = false;
   fNISTElementDataTable[2].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[2].fNIsos         = LiN;
   fNISTElementDataTable[2].fAIsos         = LiA;
   fNISTElementDataTable[2].fWIsos         = LiW;
   fNISTElementDataTable[2].fMassIsos      = LiIsoMass;


  // Z = 4-------------------------------------------------------------------------------
   static const std::string  BeS[] = {   "Be5",    "Be6",    "Be7",    "Be8",    "Be9",
                                        "Be10",   "Be11",   "Be12",   "Be13",   "Be14",
                                        "Be15",   "Be16" };
   static const int     BeN[] = {  5,   6,   7,   8,   9,  10,  11,  12,  13,  14,
                                  15,  16 };
   static const double  BeA[] = {         5.0399 * u,       6.0197264 * u,
                                     7.016928717 * u,     8.005305102 * u,
                                     9.012183065 * u,    10.013534695 * u,
                                     11.02166108 * u,      12.0269221 * u,
                                       13.036135 * u,        14.04289 * u,
                                        15.05342 * u,        16.06167 * u };
   static const double  BeW[] = { 0., 0., 0., 0., 1, 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double BeIsoMass[12];

   fNISTElementDataTable[3].fElemSymbol    = "Be";
   fNISTElementDataTable[3].fSymbols       = BeS;
   fNISTElementDataTable[3].fZ             = 4;
   fNISTElementDataTable[3].fNumOfIsotopes = 12;
   fNISTElementDataTable[3].fIsNoStableIsotope       = false;
   fNISTElementDataTable[3].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[3].fNIsos         = BeN;
   fNISTElementDataTable[3].fAIsos         = BeA;
   fNISTElementDataTable[3].fWIsos         = BeW;
   fNISTElementDataTable[3].fMassIsos      = BeIsoMass;


  // Z = 5-------------------------------------------------------------------------------
   static const std::string   BS[] = {    "B6",     "B7",     "B8",     "B9",    "B10",
                                         "B11",    "B12",    "B13",    "B14",    "B15",
                                         "B16",    "B17",    "B18",    "B19",    "B20",
                                         "B21" };
   static const int      BN[] = {  6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
                                  16,  17,  18,  19,  20,  21 };
   static const double   BA[] = {         6.0508 * u,        7.029712 * u,
                                       8.0246073 * u,      9.01332965 * u,
                                     10.01293695 * u,     11.00930536 * u,
                                      12.0143527 * u,      13.0177802 * u,
                                       14.025404 * u,       15.031088 * u,
                                       16.039842 * u,        17.04699 * u,
                                        18.05566 * u,         19.0631 * u,
                                        20.07207 * u,        21.08129 * u };
   static const double   BW[] = { 0., 0., 0., 0., 0.199, 0.801, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double BIsoMass[16];

   fNISTElementDataTable[4].fElemSymbol    = "B";
   fNISTElementDataTable[4].fSymbols       = BS;
   fNISTElementDataTable[4].fZ             = 5;
   fNISTElementDataTable[4].fNumOfIsotopes = 16;
   fNISTElementDataTable[4].fIsNoStableIsotope       = false;
   fNISTElementDataTable[4].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[4].fNIsos         = BN;
   fNISTElementDataTable[4].fAIsos         = BA;
   fNISTElementDataTable[4].fWIsos         = BW;
   fNISTElementDataTable[4].fMassIsos      = BIsoMass;


  // Z = 6-------------------------------------------------------------------------------
   static const std::string   CS[] = {    "C8",     "C9",    "C10",    "C11",    "C12",
                                         "C13",    "C14",    "C15",    "C16",    "C17",
                                         "C18",    "C19",    "C20",    "C21",    "C22",
                                         "C23" };
   static const int      CN[] = {  8,   9,  10,  11,  12,  13,  14,  15,  16,  17,
                                  18,  19,  20,  21,  22,  23 };
   static const double   CA[] = {       8.037643 * u,       9.0310372 * u,
                                     10.01685331 * u,      11.0114336 * u,
                                              12 * u,   13.0033548351 * u,
                                   14.0032419884 * u,     15.01059926 * u,
                                      16.0147013 * u,       17.022577 * u,
                                       18.026751 * u,         19.0348 * u,
                                        20.04032 * u,          21.049 * u,
                                        22.05753 * u,         23.0689 * u };
   static const double   CW[] = { 0., 0., 0., 0., 0.9893, 0.0107, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double CIsoMass[16];

   fNISTElementDataTable[5].fElemSymbol    = "C";
   fNISTElementDataTable[5].fSymbols       = CS;
   fNISTElementDataTable[5].fZ             = 6;
   fNISTElementDataTable[5].fNumOfIsotopes = 16;
   fNISTElementDataTable[5].fIsNoStableIsotope       = false;
   fNISTElementDataTable[5].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[5].fNIsos         = CN;
   fNISTElementDataTable[5].fAIsos         = CA;
   fNISTElementDataTable[5].fWIsos         = CW;
   fNISTElementDataTable[5].fMassIsos      = CIsoMass;


  // Z = 7-------------------------------------------------------------------------------
   static const std::string   NS[] = {   "N10",    "N11",    "N12",    "N13",    "N14",
                                         "N15",    "N16",    "N17",    "N18",    "N19",
                                         "N20",    "N21",    "N22",    "N23",    "N24",
                                         "N25" };
   static const int      NN[] = { 10,  11,  12,  13,  14,  15,  16,  17,  18,  19,
                                  20,  21,  22,  23,  24,  25 };
   static const double   NA[] = {       10.04165 * u,       11.026091 * u,
                                      12.0186132 * u,     13.00573861 * u,
                                   14.0030740044 * u,   15.0001088989 * u,
                                      16.0061019 * u,       17.008449 * u,
                                       18.014078 * u,       19.017022 * u,
                                       20.023366 * u,        21.02711 * u,
                                        22.03439 * u,        23.04114 * u,
                                        24.05039 * u,         25.0601 * u };
   static const double   NW[] = { 0., 0., 0., 0., 0.99636, 0.00364, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double NIsoMass[16];

   fNISTElementDataTable[6].fElemSymbol    = "N";
   fNISTElementDataTable[6].fSymbols       = NS;
   fNISTElementDataTable[6].fZ             = 7;
   fNISTElementDataTable[6].fNumOfIsotopes = 16;
   fNISTElementDataTable[6].fIsNoStableIsotope       = false;
   fNISTElementDataTable[6].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[6].fNIsos         = NN;
   fNISTElementDataTable[6].fAIsos         = NA;
   fNISTElementDataTable[6].fWIsos         = NW;
   fNISTElementDataTable[6].fMassIsos      = NIsoMass;


  // Z = 8-------------------------------------------------------------------------------
   static const std::string   OS[] = {   "O12",    "O13",    "O14",    "O15",    "O16",
                                         "O17",    "O18",    "O19",    "O20",    "O21",
                                         "O22",    "O23",    "O24",    "O25",    "O26",
                                         "O27",    "O28" };
   static const int      ON[] = { 12,  13,  14,  15,  16,  17,  18,  19,  20,  21,
                                  22,  23,  24,  25,  26,  27,  28 };
   static const double   OA[] = {      12.034262 * u,       13.024815 * u,
                                     14.00859636 * u,     15.00306562 * u,
                                   15.9949146196 * u,   16.9991317565 * u,
                                   17.9991596129 * u,       19.003578 * u,
                                     20.00407535 * u,       21.008655 * u,
                                       22.009966 * u,       23.015696 * u,
                                        24.01986 * u,        25.02936 * u,
                                        26.03729 * u,        27.04772 * u,
                                        28.05591 * u };
   static const double   OW[] = { 0., 0., 0., 0., 0.99757, 0.00038, 0.00205, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0. };
   static double OIsoMass[17];

   fNISTElementDataTable[7].fElemSymbol    = "O";
   fNISTElementDataTable[7].fSymbols       = OS;
   fNISTElementDataTable[7].fZ             = 8;
   fNISTElementDataTable[7].fNumOfIsotopes = 17;
   fNISTElementDataTable[7].fIsNoStableIsotope       = false;
   fNISTElementDataTable[7].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[7].fNIsos         = ON;
   fNISTElementDataTable[7].fAIsos         = OA;
   fNISTElementDataTable[7].fWIsos         = OW;
   fNISTElementDataTable[7].fMassIsos      = OIsoMass;


  // Z = 9-------------------------------------------------------------------------------
   static const std::string   FS[] = {   "F14",    "F15",    "F16",    "F17",    "F18",
                                         "F19",    "F20",    "F21",    "F22",    "F23",
                                         "F24",    "F25",    "F26",    "F27",    "F28",
                                         "F29",    "F30",    "F31" };
   static const int      FN[] = { 14,  15,  16,  17,  18,  19,  20,  21,  22,  23,
                                  24,  25,  26,  27,  28,  29,  30,  31 };
   static const double   FA[] = {      14.034315 * u,       15.018043 * u,
                                      16.0114657 * u,     17.00209524 * u,
                                     18.00093733 * u,   18.9984031627 * u,
                                    19.999981252 * u,      20.9999489 * u,
                                       22.002999 * u,       23.003557 * u,
                                       24.008115 * u,       25.012199 * u,
                                       26.020038 * u,        27.02644 * u,
                                        28.03534 * u,        29.04254 * u,
                                        30.05165 * u,        31.05971 * u };
   static const double   FW[] = { 0., 0., 0., 0., 0., 1, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double FIsoMass[18];

   fNISTElementDataTable[8].fElemSymbol    = "F";
   fNISTElementDataTable[8].fSymbols       = FS;
   fNISTElementDataTable[8].fZ             = 9;
   fNISTElementDataTable[8].fNumOfIsotopes = 18;
   fNISTElementDataTable[8].fIsNoStableIsotope       = false;
   fNISTElementDataTable[8].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[8].fNIsos         = FN;
   fNISTElementDataTable[8].fAIsos         = FA;
   fNISTElementDataTable[8].fWIsos         = FW;
   fNISTElementDataTable[8].fMassIsos      = FIsoMass;


  // Z = 10-------------------------------------------------------------------------------
   static const std::string  NeS[] = {  "Ne16",   "Ne17",   "Ne18",   "Ne19",   "Ne20",
                                        "Ne21",   "Ne22",   "Ne23",   "Ne24",   "Ne25",
                                        "Ne26",   "Ne27",   "Ne28",   "Ne29",   "Ne30",
                                        "Ne31",   "Ne32",   "Ne33",   "Ne34" };
   static const int     NeN[] = { 16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
                                  26,  27,  28,  29,  30,  31,  32,  33,  34 };
   static const double  NeA[] = {       16.02575 * u,     17.01771396 * u,
                                      18.0057087 * u,     19.00188091 * u,
                                   19.9924401762 * u,    20.993846685 * u,
                                    21.991385114 * u,     22.99446691 * u,
                                     23.99361065 * u,       24.997789 * u,
                                       26.000515 * u,       27.007553 * u,
                                        28.01212 * u,        29.01975 * u,
                                        30.02473 * u,         31.0331 * u,
                                        32.03972 * u,        33.04938 * u,
                                        34.05673 * u };
   static const double  NeW[] = { 0., 0., 0., 0., 0.9048, 0.0027, 0.0925, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double NeIsoMass[19];

   fNISTElementDataTable[9].fElemSymbol    = "Ne";
   fNISTElementDataTable[9].fSymbols       = NeS;
   fNISTElementDataTable[9].fZ             = 10;
   fNISTElementDataTable[9].fNumOfIsotopes = 19;
   fNISTElementDataTable[9].fIsNoStableIsotope       = false;
   fNISTElementDataTable[9].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[9].fNIsos         = NeN;
   fNISTElementDataTable[9].fAIsos         = NeA;
   fNISTElementDataTable[9].fWIsos         = NeW;
   fNISTElementDataTable[9].fMassIsos      = NeIsoMass;


  // Z = 11-------------------------------------------------------------------------------
   static const std::string  NaS[] = {  "Na18",   "Na19",   "Na20",   "Na21",   "Na22",
                                        "Na23",   "Na24",   "Na25",   "Na26",   "Na27",
                                        "Na28",   "Na29",   "Na30",   "Na31",   "Na32",
                                        "Na33",   "Na34",   "Na35",   "Na36",   "Na37" };
   static const int     NaN[] = { 18,  19,  20,  21,  22,  23,  24,  25,  26,  27,
                                  28,  29,  30,  31,  32,  33,  34,  35,  36,  37 };
   static const double  NaA[] = {       18.02688 * u,        19.01388 * u,
                                      20.0073544 * u,     20.99765469 * u,
                                     21.99443741 * u,    22.989769282 * u,
                                     23.99096295 * u,       24.989954 * u,
                                      25.9926346 * u,      26.9940765 * u,
                                       27.998939 * u,      29.0028771 * u,
                                      30.0090979 * u,       31.013163 * u,
                                        32.02019 * u,        33.02573 * u,
                                        34.03359 * u,        35.04062 * u,
                                        36.04929 * u,        37.05705 * u };
   static const double  NaW[] = { 0., 0., 0., 0., 0., 1, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double NaIsoMass[20];

   fNISTElementDataTable[10].fElemSymbol    = "Na";
   fNISTElementDataTable[10].fSymbols       = NaS;
   fNISTElementDataTable[10].fZ             = 11;
   fNISTElementDataTable[10].fNumOfIsotopes = 20;
   fNISTElementDataTable[10].fIsNoStableIsotope       = false;
   fNISTElementDataTable[10].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[10].fNIsos         = NaN;
   fNISTElementDataTable[10].fAIsos         = NaA;
   fNISTElementDataTable[10].fWIsos         = NaW;
   fNISTElementDataTable[10].fMassIsos      = NaIsoMass;


  // Z = 12-------------------------------------------------------------------------------
   static const std::string  MgS[] = {  "Mg19",   "Mg20",   "Mg21",   "Mg22",   "Mg23",
                                        "Mg24",   "Mg25",   "Mg26",   "Mg27",   "Mg28",
                                        "Mg29",   "Mg30",   "Mg31",   "Mg32",   "Mg33",
                                        "Mg34",   "Mg35",   "Mg36",   "Mg37",   "Mg38",
                                        "Mg39",   "Mg40" };
   static const int     MgN[] = { 19,  20,  21,  22,  23,  24,  25,  26,  27,  28,
                                  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,
                                  39,  40 };
   static const double  MgA[] = {      19.034169 * u,        20.01885 * u,
                                       21.011716 * u,     21.99957065 * u,
                                     22.99412421 * u,    23.985041697 * u,
                                    24.985836976 * u,    25.982592968 * u,
                                    26.984340624 * u,      27.9838767 * u,
                                       28.988617 * u,      29.9904629 * u,
                                       30.996648 * u,      31.9991102 * u,
                                      33.0053271 * u,       34.008935 * u,
                                        35.01679 * u,        36.02188 * u,
                                        37.03037 * u,        38.03658 * u,
                                        39.04538 * u,        40.05218 * u };
   static const double  MgW[] = { 0., 0., 0., 0., 0., 0.7899, 0.1, 0.1101, 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double MgIsoMass[22];

   fNISTElementDataTable[11].fElemSymbol    = "Mg";
   fNISTElementDataTable[11].fSymbols       = MgS;
   fNISTElementDataTable[11].fZ             = 12;
   fNISTElementDataTable[11].fNumOfIsotopes = 22;
   fNISTElementDataTable[11].fIsNoStableIsotope       = false;
   fNISTElementDataTable[11].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[11].fNIsos         = MgN;
   fNISTElementDataTable[11].fAIsos         = MgA;
   fNISTElementDataTable[11].fWIsos         = MgW;
   fNISTElementDataTable[11].fMassIsos      = MgIsoMass;


  // Z = 13-------------------------------------------------------------------------------
   static const std::string  AlS[] = {  "Al21",   "Al22",   "Al23",   "Al24",   "Al25",
                                        "Al26",   "Al27",   "Al28",   "Al29",   "Al30",
                                        "Al31",   "Al32",   "Al33",   "Al34",   "Al35",
                                        "Al36",   "Al37",   "Al38",   "Al39",   "Al40",
                                        "Al41",   "Al42",   "Al43" };
   static const int     AlN[] = { 21,  22,  23,  24,  25,  26,  27,  28,  29,  30,
                                  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
                                  41,  42,  43 };
   static const double  AlA[] = {       21.02897 * u,        22.01954 * u,
                                     23.00724435 * u,      23.9999489 * u,
                                      24.9904281 * u,    25.986891904 * u,
                                     26.98153853 * u,     27.98191021 * u,
                                      28.9804565 * u,        29.98296 * u,
                                       30.983945 * u,       31.988085 * u,
                                       32.990909 * u,       33.996705 * u,
                                       34.999764 * u,        36.00639 * u,
                                        37.01053 * u,         38.0174 * u,
                                        39.02254 * u,        40.03003 * u,
                                        41.03638 * u,        42.04384 * u,
                                        43.05147 * u };
   static const double  AlW[] = { 0., 0., 0., 0., 0., 0., 1, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0. };
   static double AlIsoMass[23];

   fNISTElementDataTable[12].fElemSymbol    = "Al";
   fNISTElementDataTable[12].fSymbols       = AlS;
   fNISTElementDataTable[12].fZ             = 13;
   fNISTElementDataTable[12].fNumOfIsotopes = 23;
   fNISTElementDataTable[12].fIsNoStableIsotope       = false;
   fNISTElementDataTable[12].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[12].fNIsos         = AlN;
   fNISTElementDataTable[12].fAIsos         = AlA;
   fNISTElementDataTable[12].fWIsos         = AlW;
   fNISTElementDataTable[12].fMassIsos      = AlIsoMass;


  // Z = 14-------------------------------------------------------------------------------
   static const std::string  SiS[] = {  "Si22",   "Si23",   "Si24",   "Si25",   "Si26",
                                        "Si27",   "Si28",   "Si29",   "Si30",   "Si31",
                                        "Si32",   "Si33",   "Si34",   "Si35",   "Si36",
                                        "Si37",   "Si38",   "Si39",   "Si40",   "Si41",
                                        "Si42",   "Si43",   "Si44",   "Si45" };
   static const int     SiN[] = { 22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
                                  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
                                  42,  43,  44,  45 };
   static const double  SiA[] = {       22.03579 * u,        23.02544 * u,
                                       24.011535 * u,       25.004109 * u,
                                     25.99233384 * u,     26.98670481 * u,
                                   27.9769265347 * u,   28.9764946649 * u,
                                    29.973770136 * u,    30.975363194 * u,
                                     31.97415154 * u,     32.97797696 * u,
                                       33.978576 * u,       34.984583 * u,
                                       35.986695 * u,       36.992921 * u,
                                       37.995523 * u,       39.002491 * u,
                                        40.00583 * u,        41.01301 * u,
                                        42.01778 * u,         43.0248 * u,
                                        44.03061 * u,        45.03995 * u };
   static const double  SiW[] = { 0., 0., 0., 0., 0., 0., 0.92223, 0.04685, 0.03092, 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0. };
   static double SiIsoMass[24];

   fNISTElementDataTable[13].fElemSymbol    = "Si";
   fNISTElementDataTable[13].fSymbols       = SiS;
   fNISTElementDataTable[13].fZ             = 14;
   fNISTElementDataTable[13].fNumOfIsotopes = 24;
   fNISTElementDataTable[13].fIsNoStableIsotope       = false;
   fNISTElementDataTable[13].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[13].fNIsos         = SiN;
   fNISTElementDataTable[13].fAIsos         = SiA;
   fNISTElementDataTable[13].fWIsos         = SiW;
   fNISTElementDataTable[13].fMassIsos      = SiIsoMass;


  // Z = 15-------------------------------------------------------------------------------
   static const std::string   PS[] = {   "P24",    "P25",    "P26",    "P27",    "P28",
                                         "P29",    "P30",    "P31",    "P32",    "P33",
                                         "P34",    "P35",    "P36",    "P37",    "P38",
                                         "P39",    "P40",    "P41",    "P42",    "P43",
                                         "P44",    "P45",    "P46",    "P47" };
   static const int      PN[] = { 24,  25,  26,  27,  28,  29,  30,  31,  32,  33,
                                  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,
                                  44,  45,  46,  47 };
   static const double   PA[] = {       24.03577 * u,        25.02119 * u,
                                        26.01178 * u,       26.999224 * u,
                                      27.9923266 * u,     28.98180079 * u,
                                     29.97831375 * u,   30.9737619984 * u,
                                    31.973907643 * u,      32.9717257 * u,
                                     33.97364589 * u,      34.9733141 * u,
                                        35.97826 * u,       36.979607 * u,
                                       37.984252 * u,       38.986227 * u,
                                        39.99133 * u,       40.994654 * u,
                                        42.00108 * u,        43.00502 * u,
                                        44.01121 * u,        45.01645 * u,
                                        46.02446 * u,        47.03139 * u };
   static const double   PW[] = { 0., 0., 0., 0., 0., 0., 0., 1, 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0. };
   static double PIsoMass[24];

   fNISTElementDataTable[14].fElemSymbol    = "P";
   fNISTElementDataTable[14].fSymbols       = PS;
   fNISTElementDataTable[14].fZ             = 15;
   fNISTElementDataTable[14].fNumOfIsotopes = 24;
   fNISTElementDataTable[14].fIsNoStableIsotope       = false;
   fNISTElementDataTable[14].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[14].fNIsos         = PN;
   fNISTElementDataTable[14].fAIsos         = PA;
   fNISTElementDataTable[14].fWIsos         = PW;
   fNISTElementDataTable[14].fMassIsos      = PIsoMass;


  // Z = 16-------------------------------------------------------------------------------
   static const std::string   SS[] = {   "S26",    "S27",    "S28",    "S29",    "S30",
                                         "S31",    "S32",    "S33",    "S34",    "S35",
                                         "S36",    "S37",    "S38",    "S39",    "S40",
                                         "S41",    "S42",    "S43",    "S44",    "S45",
                                         "S46",    "S47",    "S48",    "S49" };
   static const int      SN[] = { 26,  27,  28,  29,  30,  31,  32,  33,  34,  35,
                                  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,
                                  46,  47,  48,  49 };
   static const double   SA[] = {       26.02907 * u,        27.01828 * u,
                                        28.00437 * u,       28.996611 * u,
                                     29.98490703 * u,     30.97955701 * u,
                                   31.9720711744 * u,   32.9714589098 * u,
                                    33.967867004 * u,     34.96903231 * u,
                                     35.96708071 * u,     36.97112551 * u,
                                      37.9711633 * u,       38.975134 * u,
                                      39.9754826 * u,      40.9795935 * u,
                                      41.9810651 * u,      42.9869076 * u,
                                      43.9901188 * u,        44.99572 * u,
                                        46.00004 * u,        47.00795 * u,
                                         48.0137 * u,        49.02276 * u };
   static const double   SW[] = { 0., 0., 0., 0., 0., 0., 0.9499, 0.0075, 0.0425, 0.,
                                  0.0001, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0. };
   static double SIsoMass[24];

   fNISTElementDataTable[15].fElemSymbol    = "S";
   fNISTElementDataTable[15].fSymbols       = SS;
   fNISTElementDataTable[15].fZ             = 16;
   fNISTElementDataTable[15].fNumOfIsotopes = 24;
   fNISTElementDataTable[15].fIsNoStableIsotope       = false;
   fNISTElementDataTable[15].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[15].fNIsos         = SN;
   fNISTElementDataTable[15].fAIsos         = SA;
   fNISTElementDataTable[15].fWIsos         = SW;
   fNISTElementDataTable[15].fMassIsos      = SIsoMass;


  // Z = 17-------------------------------------------------------------------------------
   static const std::string  ClS[] = {  "Cl28",   "Cl29",   "Cl30",   "Cl31",   "Cl32",
                                        "Cl33",   "Cl34",   "Cl35",   "Cl36",   "Cl37",
                                        "Cl38",   "Cl39",   "Cl40",   "Cl41",   "Cl42",
                                        "Cl43",   "Cl44",   "Cl45",   "Cl46",   "Cl47",
                                        "Cl48",   "Cl49",   "Cl50",   "Cl51" };
   static const int     ClN[] = { 28,  29,  30,  31,  32,  33,  34,  35,  36,  37,
                                  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
                                  48,  49,  50,  51 };
   static const double  ClA[] = {       28.02954 * u,        29.01478 * u,
                                        30.00477 * u,       30.992414 * u,
                                     31.98568464 * u,     32.97745199 * u,
                                    33.973762485 * u,    34.968852682 * u,
                                    35.968306809 * u,    36.965902602 * u,
                                     37.96801044 * u,      38.9680082 * u,
                                       39.970415 * u,       40.970685 * u,
                                        41.97325 * u,        42.97389 * u,
                                        43.97787 * u,        44.98029 * u,
                                        45.98517 * u,        46.98916 * u,
                                        47.99564 * u,        49.00123 * u,
                                        50.00905 * u,        51.01554 * u };
   static const double  ClW[] = { 0., 0., 0., 0., 0., 0., 0., 0.7576, 0., 0.2424,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0. };
   static double ClIsoMass[24];

   fNISTElementDataTable[16].fElemSymbol    = "Cl";
   fNISTElementDataTable[16].fSymbols       = ClS;
   fNISTElementDataTable[16].fZ             = 17;
   fNISTElementDataTable[16].fNumOfIsotopes = 24;
   fNISTElementDataTable[16].fIsNoStableIsotope       = false;
   fNISTElementDataTable[16].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[16].fNIsos         = ClN;
   fNISTElementDataTable[16].fAIsos         = ClA;
   fNISTElementDataTable[16].fWIsos         = ClW;
   fNISTElementDataTable[16].fMassIsos      = ClIsoMass;


  // Z = 18-------------------------------------------------------------------------------
   static const std::string  ArS[] = {  "Ar30",   "Ar31",   "Ar32",   "Ar33",   "Ar34",
                                        "Ar35",   "Ar36",   "Ar37",   "Ar38",   "Ar39",
                                        "Ar40",   "Ar41",   "Ar42",   "Ar43",   "Ar44",
                                        "Ar45",   "Ar46",   "Ar47",   "Ar48",   "Ar49",
                                        "Ar50",   "Ar51",   "Ar52",   "Ar53" };
   static const int     ArN[] = { 30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
                                  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,
                                  50,  51,  52,  53 };
   static const double  ArA[] = {       30.02307 * u,        31.01212 * u,
                                      31.9976378 * u,     32.98992555 * u,
                                     33.98027009 * u,     34.97525759 * u,
                                    35.967545105 * u,     36.96677633 * u,
                                     37.96273211 * u,       38.964313 * u,
                                   39.9623831237 * u,     40.96450057 * u,
                                      41.9630457 * u,      42.9656361 * u,
                                      43.9649238 * u,     44.96803973 * u,
                                       45.968083 * u,       46.972935 * u,
                                        47.97591 * u,         48.9819 * u,
                                        49.98613 * u,         50.9937 * u,
                                        51.99896 * u,        53.00729 * u };
   static const double  ArW[] = { 0., 0., 0., 0., 0., 0., 0.003336, 0., 0.000629, 0.,
                                  0.996035, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0. };
   static double ArIsoMass[24];

   fNISTElementDataTable[17].fElemSymbol    = "Ar";
   fNISTElementDataTable[17].fSymbols       = ArS;
   fNISTElementDataTable[17].fZ             = 18;
   fNISTElementDataTable[17].fNumOfIsotopes = 24;
   fNISTElementDataTable[17].fIsNoStableIsotope       = false;
   fNISTElementDataTable[17].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[17].fNIsos         = ArN;
   fNISTElementDataTable[17].fAIsos         = ArA;
   fNISTElementDataTable[17].fWIsos         = ArW;
   fNISTElementDataTable[17].fMassIsos      = ArIsoMass;


  // Z = 19-------------------------------------------------------------------------------
   static const std::string   KS[] = {   "K32",    "K33",    "K34",    "K35",    "K36",
                                         "K37",    "K38",    "K39",    "K40",    "K41",
                                         "K42",    "K43",    "K44",    "K45",    "K46",
                                         "K47",    "K48",    "K49",    "K50",    "K51",
                                         "K52",    "K53",    "K54",    "K55",    "K56" };
   static const int      KN[] = { 32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
                                  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,
                                  52,  53,  54,  55,  56 };
   static const double   KA[] = {       32.02265 * u,        33.00756 * u,
                                        33.99869 * u,     34.98800541 * u,
                                     35.98130201 * u,     36.97337589 * u,
                                     37.96908112 * u,   38.9637064864 * u,
                                    39.963998166 * u,   40.9618252579 * u,
                                     41.96240231 * u,      42.9607347 * u,
                                     43.96158699 * u,     44.96069149 * u,
                                     45.96198159 * u,      46.9616616 * u,
                                     47.96534119 * u,     48.96821075 * u,
                                        49.97238 * u,       50.975828 * u,
                                        51.98224 * u,        52.98746 * u,
                                        53.99463 * u,        55.00076 * u,
                                        56.00851 * u };
   static const double   KW[] = { 0., 0., 0., 0., 0., 0., 0., 0.932581, 0.000117, 0.067302,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double KIsoMass[25];

   fNISTElementDataTable[18].fElemSymbol    = "K";
   fNISTElementDataTable[18].fSymbols       = KS;
   fNISTElementDataTable[18].fZ             = 19;
   fNISTElementDataTable[18].fNumOfIsotopes = 25;
   fNISTElementDataTable[18].fIsNoStableIsotope       = false;
   fNISTElementDataTable[18].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[18].fNIsos         = KN;
   fNISTElementDataTable[18].fAIsos         = KA;
   fNISTElementDataTable[18].fWIsos         = KW;
   fNISTElementDataTable[18].fMassIsos      = KIsoMass;


  // Z = 20-------------------------------------------------------------------------------
   static const std::string  CaS[] = {  "Ca34",   "Ca35",   "Ca36",   "Ca37",   "Ca38",
                                        "Ca39",   "Ca40",   "Ca41",   "Ca42",   "Ca43",
                                        "Ca44",   "Ca45",   "Ca46",   "Ca47",   "Ca48",
                                        "Ca49",   "Ca50",   "Ca51",   "Ca52",   "Ca53",
                                        "Ca54",   "Ca55",   "Ca56",   "Ca57",   "Ca58" };
   static const int     CaN[] = { 34,  35,  36,  37,  38,  39,  40,  41,  42,  43,
                                  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,
                                  54,  55,  56,  57,  58 };
   static const double  CaA[] = {       34.01487 * u,        35.00514 * u,
                                       35.993074 * u,     36.98589785 * u,
                                     37.97631922 * u,     38.97071081 * u,
                                    39.962590863 * u,     40.96227792 * u,
                                     41.95861783 * u,     42.95876644 * u,
                                     43.95548156 * u,     44.95618635 * u,
                                       45.953689 * u,      46.9545424 * u,
                                     47.95252276 * u,     48.95566274 * u,
                                      49.9574992 * u,       50.960989 * u,
                                       51.963217 * u,        52.96945 * u,
                                         53.9734 * u,         54.9803 * u,
                                        55.98508 * u,        56.99262 * u,
                                        57.99794 * u };
   static const double  CaW[] = { 0., 0., 0., 0., 0., 0., 0.96941, 0., 0.00647, 0.00135,
                                  0.02086, 0., 4e-05, 0., 0.00187, 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double CaIsoMass[25];

   fNISTElementDataTable[19].fElemSymbol    = "Ca";
   fNISTElementDataTable[19].fSymbols       = CaS;
   fNISTElementDataTable[19].fZ             = 20;
   fNISTElementDataTable[19].fNumOfIsotopes = 25;
   fNISTElementDataTable[19].fIsNoStableIsotope       = false;
   fNISTElementDataTable[19].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[19].fNIsos         = CaN;
   fNISTElementDataTable[19].fAIsos         = CaA;
   fNISTElementDataTable[19].fWIsos         = CaW;
   fNISTElementDataTable[19].fMassIsos      = CaIsoMass;


  // Z = 21-------------------------------------------------------------------------------
   static const std::string  ScS[] = {  "Sc36",   "Sc37",   "Sc38",   "Sc39",   "Sc40",
                                        "Sc41",   "Sc42",   "Sc43",   "Sc44",   "Sc45",
                                        "Sc46",   "Sc47",   "Sc48",   "Sc49",   "Sc50",
                                        "Sc51",   "Sc52",   "Sc53",   "Sc54",   "Sc55",
                                        "Sc56",   "Sc57",   "Sc58",   "Sc59",   "Sc60",
                                        "Sc61" };
   static const int     ScN[] = { 36,  37,  38,  39,  40,  41,  42,  43,  44,  45,
                                  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,
                                  56,  57,  58,  59,  60,  61 };
   static const double  ScA[] = {       36.01648 * u,        37.00374 * u,
                                        37.99512 * u,       38.984785 * u,
                                      39.9779673 * u,    40.969251105 * u,
                                     41.96551653 * u,      42.9611505 * u,
                                      43.9594029 * u,     44.95590828 * u,
                                     45.95516826 * u,      46.9524037 * u,
                                      47.9522236 * u,      48.9500146 * u,
                                       49.952176 * u,       50.953592 * u,
                                        51.95688 * u,        52.95909 * u,
                                        53.96393 * u,        54.96782 * u,
                                        55.97345 * u,        56.97777 * u,
                                        57.98403 * u,        58.98894 * u,
                                        59.99565 * u,          61.001 * u };
   static const double  ScW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 1,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double ScIsoMass[26];

   fNISTElementDataTable[20].fElemSymbol    = "Sc";
   fNISTElementDataTable[20].fSymbols       = ScS;
   fNISTElementDataTable[20].fZ             = 21;
   fNISTElementDataTable[20].fNumOfIsotopes = 26;
   fNISTElementDataTable[20].fIsNoStableIsotope       = false;
   fNISTElementDataTable[20].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[20].fNIsos         = ScN;
   fNISTElementDataTable[20].fAIsos         = ScA;
   fNISTElementDataTable[20].fWIsos         = ScW;
   fNISTElementDataTable[20].fMassIsos      = ScIsoMass;


  // Z = 22-------------------------------------------------------------------------------
   static const std::string  TiS[] = {  "Ti38",   "Ti39",   "Ti40",   "Ti41",   "Ti42",
                                        "Ti43",   "Ti44",   "Ti45",   "Ti46",   "Ti47",
                                        "Ti48",   "Ti49",   "Ti50",   "Ti51",   "Ti52",
                                        "Ti53",   "Ti54",   "Ti55",   "Ti56",   "Ti57",
                                        "Ti58",   "Ti59",   "Ti60",   "Ti61",   "Ti62",
                                        "Ti63" };
   static const int     TiN[] = { 38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
                                  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,
                                  58,  59,  60,  61,  62,  63 };
   static const double  TiA[] = {       38.01145 * u,        39.00236 * u,
                                         39.9905 * u,       40.983148 * u,
                                     41.97304903 * u,      42.9685225 * u,
                                     43.95968995 * u,     44.95812198 * u,
                                     45.95262772 * u,     46.95175879 * u,
                                     47.94794198 * u,     48.94786568 * u,
                                     49.94478689 * u,     50.94661065 * u,
                                       51.946893 * u,        52.94973 * u,
                                        53.95105 * u,        54.95527 * u,
                                        55.95791 * u,        56.96364 * u,
                                         57.9666 * u,        58.97247 * u,
                                        59.97603 * u,        60.98245 * u,
                                        61.98651 * u,        62.99375 * u };
   static const double  TiW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0.0825, 0.0744,
                                  0.7372, 0.0541, 0.0518, 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double TiIsoMass[26];

   fNISTElementDataTable[21].fElemSymbol    = "Ti";
   fNISTElementDataTable[21].fSymbols       = TiS;
   fNISTElementDataTable[21].fZ             = 22;
   fNISTElementDataTable[21].fNumOfIsotopes = 26;
   fNISTElementDataTable[21].fIsNoStableIsotope       = false;
   fNISTElementDataTable[21].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[21].fNIsos         = TiN;
   fNISTElementDataTable[21].fAIsos         = TiA;
   fNISTElementDataTable[21].fWIsos         = TiW;
   fNISTElementDataTable[21].fMassIsos      = TiIsoMass;


  // Z = 23-------------------------------------------------------------------------------
   static const std::string   VS[] = {   "V40",    "V41",    "V42",    "V43",    "V44",
                                         "V45",    "V46",    "V47",    "V48",    "V49",
                                         "V50",    "V51",    "V52",    "V53",    "V54",
                                         "V55",    "V56",    "V57",    "V58",    "V59",
                                         "V60",    "V61",    "V62",    "V63",    "V64",
                                         "V65",    "V66" };
   static const int      VN[] = { 40,  41,  42,  43,  44,  45,  46,  47,  48,  49,
                                  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,
                                  60,  61,  62,  63,  64,  65,  66 };
   static const double   VA[] = {       40.01276 * u,        41.00021 * u,
                                        41.99182 * u,       42.980766 * u,
                                        43.97411 * u,      44.9657748 * u,
                                     45.96019878 * u,     46.95490491 * u,
                                      47.9522522 * u,      48.9485118 * u,
                                     49.94715601 * u,     50.94395704 * u,
                                     51.94477301 * u,      52.9443367 * u,
                                       53.946439 * u,        54.94724 * u,
                                        55.95048 * u,        56.95252 * u,
                                        57.95672 * u,        58.95939 * u,
                                        59.96431 * u,        60.96725 * u,
                                        61.97265 * u,        62.97639 * u,
                                        63.98264 * u,         64.9875 * u,
                                        65.99398 * u };
   static const double   VW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0.0025, 0.9975, 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0. };
   static double VIsoMass[27];

   fNISTElementDataTable[22].fElemSymbol    = "V";
   fNISTElementDataTable[22].fSymbols       = VS;
   fNISTElementDataTable[22].fZ             = 23;
   fNISTElementDataTable[22].fNumOfIsotopes = 27;
   fNISTElementDataTable[22].fIsNoStableIsotope       = false;
   fNISTElementDataTable[22].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[22].fNIsos         = VN;
   fNISTElementDataTable[22].fAIsos         = VA;
   fNISTElementDataTable[22].fWIsos         = VW;
   fNISTElementDataTable[22].fMassIsos      = VIsoMass;


  // Z = 24-------------------------------------------------------------------------------
   static const std::string  CrS[] = {  "Cr42",   "Cr43",   "Cr44",   "Cr45",   "Cr46",
                                        "Cr47",   "Cr48",   "Cr49",   "Cr50",   "Cr51",
                                        "Cr52",   "Cr53",   "Cr54",   "Cr55",   "Cr56",
                                        "Cr57",   "Cr58",   "Cr59",   "Cr60",   "Cr61",
                                        "Cr62",   "Cr63",   "Cr64",   "Cr65",   "Cr66",
                                        "Cr67",   "Cr68" };
   static const int     CrN[] = { 42,  43,  44,  45,  46,  47,  48,  49,  50,  51,
                                  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,
                                  62,  63,  64,  65,  66,  67,  68 };
   static const double  CrA[] = {        42.0067 * u,        42.99753 * u,
                                        43.98536 * u,        44.97905 * u,
                                       45.968359 * u,      46.9628974 * u,
                                      47.9540291 * u,      48.9513333 * u,
                                     49.94604183 * u,     50.94476502 * u,
                                     51.94050623 * u,     52.94064815 * u,
                                     53.93887916 * u,     54.94083843 * u,
                                      55.9406531 * u,       56.943613 * u,
                                        57.94435 * u,        58.94859 * u,
                                        59.95008 * u,        60.95442 * u,
                                         61.9561 * u,        62.96165 * u,
                                        63.96408 * u,        64.96996 * u,
                                        65.97366 * u,        66.98016 * u,
                                        67.98403 * u };
   static const double  CrW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0.04345, 0.,
                                  0.83789, 0.09501, 0.02365, 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0. };
   static double CrIsoMass[27];

   fNISTElementDataTable[23].fElemSymbol    = "Cr";
   fNISTElementDataTable[23].fSymbols       = CrS;
   fNISTElementDataTable[23].fZ             = 24;
   fNISTElementDataTable[23].fNumOfIsotopes = 27;
   fNISTElementDataTable[23].fIsNoStableIsotope       = false;
   fNISTElementDataTable[23].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[23].fNIsos         = CrN;
   fNISTElementDataTable[23].fAIsos         = CrA;
   fNISTElementDataTable[23].fWIsos         = CrW;
   fNISTElementDataTable[23].fMassIsos      = CrIsoMass;


  // Z = 25-------------------------------------------------------------------------------
   static const std::string  MnS[] = {  "Mn44",   "Mn45",   "Mn46",   "Mn47",   "Mn48",
                                        "Mn49",   "Mn50",   "Mn51",   "Mn52",   "Mn53",
                                        "Mn54",   "Mn55",   "Mn56",   "Mn57",   "Mn58",
                                        "Mn59",   "Mn60",   "Mn61",   "Mn62",   "Mn63",
                                        "Mn64",   "Mn65",   "Mn66",   "Mn67",   "Mn68",
                                        "Mn69",   "Mn70",   "Mn71" };
   static const int     MnN[] = { 44,  45,  46,  47,  48,  49,  50,  51,  52,  53,
                                  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
                                  64,  65,  66,  67,  68,  69,  70,  71 };
   static const double  MnA[] = {       44.00715 * u,        44.99449 * u,
                                        45.98609 * u,       46.975775 * u,
                                        47.96852 * u,       48.959595 * u,
                                     49.95423778 * u,     50.94820847 * u,
                                      51.9455639 * u,     52.94128889 * u,
                                      53.9403576 * u,     54.93804391 * u,
                                     55.93890369 * u,      56.9382861 * u,
                                      57.9400666 * u,      58.9403911 * u,
                                      59.9431366 * u,      60.9444525 * u,
                                        61.94795 * u,      62.9496647 * u,
                                      63.9538494 * u,      64.9560198 * u,
                                       65.960547 * u,        66.96424 * u,
                                        67.96962 * u,        68.97366 * u,
                                        69.97937 * u,        70.98368 * u };
   static const double  MnW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 1, 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double MnIsoMass[28];

   fNISTElementDataTable[24].fElemSymbol    = "Mn";
   fNISTElementDataTable[24].fSymbols       = MnS;
   fNISTElementDataTable[24].fZ             = 25;
   fNISTElementDataTable[24].fNumOfIsotopes = 28;
   fNISTElementDataTable[24].fIsNoStableIsotope       = false;
   fNISTElementDataTable[24].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[24].fNIsos         = MnN;
   fNISTElementDataTable[24].fAIsos         = MnA;
   fNISTElementDataTable[24].fWIsos         = MnW;
   fNISTElementDataTable[24].fMassIsos      = MnIsoMass;


  // Z = 26-------------------------------------------------------------------------------
   static const std::string  FeS[] = {  "Fe45",   "Fe46",   "Fe47",   "Fe48",   "Fe49",
                                        "Fe50",   "Fe51",   "Fe52",   "Fe53",   "Fe54",
                                        "Fe55",   "Fe56",   "Fe57",   "Fe58",   "Fe59",
                                        "Fe60",   "Fe61",   "Fe62",   "Fe63",   "Fe64",
                                        "Fe65",   "Fe66",   "Fe67",   "Fe68",   "Fe69",
                                        "Fe70",   "Fe71",   "Fe72",   "Fe73",   "Fe74" };
   static const int     FeN[] = { 45,  46,  47,  48,  49,  50,  51,  52,  53,  54,
                                  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
                                  65,  66,  67,  68,  69,  70,  71,  72,  73,  74 };
   static const double  FeA[] = {       45.01442 * u,        46.00063 * u,
                                        46.99185 * u,        47.98023 * u,
                                       48.973429 * u,       49.962975 * u,
                                       50.956841 * u,      51.9481131 * u,
                                      52.9453064 * u,     53.93960899 * u,
                                     54.93829199 * u,     55.93493633 * u,
                                     56.93539284 * u,     57.93327443 * u,
                                     58.93487434 * u,      59.9340711 * u,
                                      60.9367462 * u,      61.9367918 * u,
                                      62.9402727 * u,      63.9409878 * u,
                                      64.9450115 * u,        65.94625 * u,
                                        66.95054 * u,        67.95295 * u,
                                        68.95807 * u,        69.96102 * u,
                                        70.96672 * u,        71.96983 * u,
                                        72.97572 * u,        73.97935 * u };
   static const double  FeW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.05845,
                                  0., 0.91754, 0.02119, 0.00282, 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double FeIsoMass[30];

   fNISTElementDataTable[25].fElemSymbol    = "Fe";
   fNISTElementDataTable[25].fSymbols       = FeS;
   fNISTElementDataTable[25].fZ             = 26;
   fNISTElementDataTable[25].fNumOfIsotopes = 30;
   fNISTElementDataTable[25].fIsNoStableIsotope       = false;
   fNISTElementDataTable[25].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[25].fNIsos         = FeN;
   fNISTElementDataTable[25].fAIsos         = FeA;
   fNISTElementDataTable[25].fWIsos         = FeW;
   fNISTElementDataTable[25].fMassIsos      = FeIsoMass;


  // Z = 27-------------------------------------------------------------------------------
   static const std::string  CoS[] = {  "Co47",   "Co48",   "Co49",   "Co50",   "Co51",
                                        "Co52",   "Co53",   "Co54",   "Co55",   "Co56",
                                        "Co57",   "Co58",   "Co59",   "Co60",   "Co61",
                                        "Co62",   "Co63",   "Co64",   "Co65",   "Co66",
                                        "Co67",   "Co68",   "Co69",   "Co70",   "Co71",
                                        "Co72",   "Co73",   "Co74",   "Co75",   "Co76" };
   static const int     CoN[] = { 47,  48,  49,  50,  51,  52,  53,  54,  55,  56,
                                  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,
                                  67,  68,  69,  70,  71,  72,  73,  74,  75,  76 };
   static const double  CoA[] = {       47.01057 * u,        48.00093 * u,
                                        48.98891 * u,        49.98091 * u,
                                       50.970647 * u,        51.96351 * u,
                                      52.9542041 * u,     53.94845987 * u,
                                      54.9419972 * u,      55.9398388 * u,
                                     56.93629057 * u,      57.9357521 * u,
                                     58.93319429 * u,      59.9338163 * u,
                                     60.93247662 * u,       61.934059 * u,
                                         62.9336 * u,       63.935811 * u,
                                      64.9364621 * u,       65.939443 * u,
                                      66.9406096 * u,        67.94426 * u,
                                        68.94614 * u,        69.94963 * u,
                                        70.95237 * u,        71.95729 * u,
                                        72.96039 * u,        73.96515 * u,
                                        74.96876 * u,        75.97413 * u };
   static const double  CoW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 1, 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double CoIsoMass[30];

   fNISTElementDataTable[26].fElemSymbol    = "Co";
   fNISTElementDataTable[26].fSymbols       = CoS;
   fNISTElementDataTable[26].fZ             = 27;
   fNISTElementDataTable[26].fNumOfIsotopes = 30;
   fNISTElementDataTable[26].fIsNoStableIsotope       = false;
   fNISTElementDataTable[26].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[26].fNIsos         = CoN;
   fNISTElementDataTable[26].fAIsos         = CoA;
   fNISTElementDataTable[26].fWIsos         = CoW;
   fNISTElementDataTable[26].fMassIsos      = CoIsoMass;


  // Z = 28-------------------------------------------------------------------------------
   static const std::string  NiS[] = {  "Ni48",   "Ni49",   "Ni50",   "Ni51",   "Ni52",
                                        "Ni53",   "Ni54",   "Ni55",   "Ni56",   "Ni57",
                                        "Ni58",   "Ni59",   "Ni60",   "Ni61",   "Ni62",
                                        "Ni63",   "Ni64",   "Ni65",   "Ni66",   "Ni67",
                                        "Ni68",   "Ni69",   "Ni70",   "Ni71",   "Ni72",
                                        "Ni73",   "Ni74",   "Ni75",   "Ni76",   "Ni77",
                                        "Ni78",   "Ni79" };
   static const int     NiN[] = { 48,  49,  50,  51,  52,  53,  54,  55,  56,  57,
                                  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,
                                  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
                                  78,  79 };
   static const double  NiA[] = {       48.01769 * u,         49.0077 * u,
                                        49.99474 * u,        50.98611 * u,
                                         51.9748 * u,        52.96819 * u,
                                       53.957892 * u,     54.95133063 * u,
                                     55.94212855 * u,     56.93979218 * u,
                                     57.93534241 * u,      58.9343462 * u,
                                     59.93078588 * u,     60.93105557 * u,
                                     61.92834537 * u,     62.92966963 * u,
                                     63.92796682 * u,     64.93008517 * u,
                                      65.9291393 * u,      66.9315694 * u,
                                      67.9318688 * u,      68.9356103 * u,
                                      69.9364313 * u,       70.940519 * u,
                                      71.9417859 * u,      72.9462067 * u,
                                        73.94798 * u,         74.9525 * u,
                                        75.95533 * u,        76.96055 * u,
                                        77.96336 * u,        78.97025 * u };
   static const double  NiW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0.68077, 0., 0.26223, 0.011399, 0.036346, 0., 0.009255, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double NiIsoMass[32];

   fNISTElementDataTable[27].fElemSymbol    = "Ni";
   fNISTElementDataTable[27].fSymbols       = NiS;
   fNISTElementDataTable[27].fZ             = 28;
   fNISTElementDataTable[27].fNumOfIsotopes = 32;
   fNISTElementDataTable[27].fIsNoStableIsotope       = false;
   fNISTElementDataTable[27].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[27].fNIsos         = NiN;
   fNISTElementDataTable[27].fAIsos         = NiA;
   fNISTElementDataTable[27].fWIsos         = NiW;
   fNISTElementDataTable[27].fMassIsos      = NiIsoMass;


  // Z = 29-------------------------------------------------------------------------------
   static const std::string  CuS[] = {  "Cu52",   "Cu53",   "Cu54",   "Cu55",   "Cu56",
                                        "Cu57",   "Cu58",   "Cu59",   "Cu60",   "Cu61",
                                        "Cu62",   "Cu63",   "Cu64",   "Cu65",   "Cu66",
                                        "Cu67",   "Cu68",   "Cu69",   "Cu70",   "Cu71",
                                        "Cu72",   "Cu73",   "Cu74",   "Cu75",   "Cu76",
                                        "Cu77",   "Cu78",   "Cu79",   "Cu80",   "Cu81",
                                        "Cu82" };
   static const int     CuN[] = { 52,  53,  54,  55,  56,  57,  58,  59,  60,  61,
                                  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,
                                  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,
                                  82 };
   static const double  CuA[] = {       51.99671 * u,        52.98459 * u,
                                        53.97666 * u,        54.96604 * u,
                                        55.95895 * u,      56.9492125 * u,
                                     57.94453305 * u,     58.93949748 * u,
                                      59.9373645 * u,      60.9334576 * u,
                                     61.93259541 * u,     62.92959772 * u,
                                     63.92976434 * u,      64.9277897 * u,
                                     65.92886903 * u,      66.9277303 * u,
                                      67.9296109 * u,      68.9294293 * u,
                                      69.9323921 * u,      70.9326768 * u,
                                      71.9358203 * u,      72.9366744 * u,
                                      73.9398749 * u,      74.9415226 * u,
                                       75.945275 * u,        76.94792 * u,
                                        77.95223 * u,        78.95502 * u,
                                        79.96089 * u,        80.96587 * u,
                                        81.97244 * u };
   static const double  CuW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0.6915, 0., 0.3085, 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0. };
   static double CuIsoMass[31];

   fNISTElementDataTable[28].fElemSymbol    = "Cu";
   fNISTElementDataTable[28].fSymbols       = CuS;
   fNISTElementDataTable[28].fZ             = 29;
   fNISTElementDataTable[28].fNumOfIsotopes = 31;
   fNISTElementDataTable[28].fIsNoStableIsotope       = false;
   fNISTElementDataTable[28].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[28].fNIsos         = CuN;
   fNISTElementDataTable[28].fAIsos         = CuA;
   fNISTElementDataTable[28].fWIsos         = CuW;
   fNISTElementDataTable[28].fMassIsos      = CuIsoMass;


  // Z = 30-------------------------------------------------------------------------------
   static const std::string  ZnS[] = {  "Zn54",   "Zn55",   "Zn56",   "Zn57",   "Zn58",
                                        "Zn59",   "Zn60",   "Zn61",   "Zn62",   "Zn63",
                                        "Zn64",   "Zn65",   "Zn66",   "Zn67",   "Zn68",
                                        "Zn69",   "Zn70",   "Zn71",   "Zn72",   "Zn73",
                                        "Zn74",   "Zn75",   "Zn76",   "Zn77",   "Zn78",
                                        "Zn79",   "Zn80",   "Zn81",   "Zn82",   "Zn83",
                                        "Zn84",   "Zn85" };
   static const int     ZnN[] = { 54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
                                  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,
                                  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,
                                  84,  85 };
   static const double  ZnA[] = {       53.99204 * u,        54.98398 * u,
                                        55.97254 * u,        56.96506 * u,
                                       57.954591 * u,     58.94931266 * u,
                                      59.9418421 * u,       60.939507 * u,
                                     61.93433397 * u,      62.9332115 * u,
                                     63.92914201 * u,     64.92924077 * u,
                                     65.92603381 * u,     66.92712775 * u,
                                     67.92484455 * u,      68.9265507 * u,
                                      69.9253192 * u,      70.9277196 * u,
                                      71.9268428 * u,      72.9295826 * u,
                                      73.9294073 * u,      74.9328402 * u,
                                       75.933115 * u,      76.9368872 * u,
                                      77.9382892 * u,      78.9426381 * u,
                                      79.9445529 * u,      80.9504026 * u,
                                        81.95426 * u,        82.96056 * u,
                                        83.96521 * u,        84.97226 * u };
   static const double  ZnW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0.4917, 0., 0.2773, 0.0404, 0.1845, 0., 0.0061, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double ZnIsoMass[32];

   fNISTElementDataTable[29].fElemSymbol    = "Zn";
   fNISTElementDataTable[29].fSymbols       = ZnS;
   fNISTElementDataTable[29].fZ             = 30;
   fNISTElementDataTable[29].fNumOfIsotopes = 32;
   fNISTElementDataTable[29].fIsNoStableIsotope       = false;
   fNISTElementDataTable[29].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[29].fNIsos         = ZnN;
   fNISTElementDataTable[29].fAIsos         = ZnA;
   fNISTElementDataTable[29].fWIsos         = ZnW;
   fNISTElementDataTable[29].fMassIsos      = ZnIsoMass;


  // Z = 31-------------------------------------------------------------------------------
   static const std::string  GaS[] = {  "Ga56",   "Ga57",   "Ga58",   "Ga59",   "Ga60",
                                        "Ga61",   "Ga62",   "Ga63",   "Ga64",   "Ga65",
                                        "Ga66",   "Ga67",   "Ga68",   "Ga69",   "Ga70",
                                        "Ga71",   "Ga72",   "Ga73",   "Ga74",   "Ga75",
                                        "Ga76",   "Ga77",   "Ga78",   "Ga79",   "Ga80",
                                        "Ga81",   "Ga82",   "Ga83",   "Ga84",   "Ga85",
                                        "Ga86",   "Ga87" };
   static const int     GaN[] = { 56,  57,  58,  59,  60,  61,  62,  63,  64,  65,
                                  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,
                                  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,
                                  86,  87 };
   static const double  GaA[] = {       55.99536 * u,         56.9832 * u,
                                        57.97478 * u,        58.96353 * u,
                                        59.95729 * u,       60.949399 * u,
                                     61.94419025 * u,      62.9392942 * u,
                                      63.9368404 * u,     64.93273459 * u,
                                      65.9315894 * u,      66.9282025 * u,
                                      67.9279805 * u,      68.9255735 * u,
                                      69.9260219 * u,     70.92470258 * u,
                                     71.92636747 * u,      72.9251747 * u,
                                      73.9269457 * u,      74.9265002 * u,
                                      75.9288276 * u,      76.9291543 * u,
                                      77.9316088 * u,      78.9328523 * u,
                                      79.9364208 * u,      80.9381338 * u,
                                      81.9431765 * u,      82.9471203 * u,
                                        83.95246 * u,        84.95699 * u,
                                        85.96301 * u,        86.96824 * u };
   static const double  GaW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0.60108, 0., 0.39892, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double GaIsoMass[32];

   fNISTElementDataTable[30].fElemSymbol    = "Ga";
   fNISTElementDataTable[30].fSymbols       = GaS;
   fNISTElementDataTable[30].fZ             = 31;
   fNISTElementDataTable[30].fNumOfIsotopes = 32;
   fNISTElementDataTable[30].fIsNoStableIsotope       = false;
   fNISTElementDataTable[30].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[30].fNIsos         = GaN;
   fNISTElementDataTable[30].fAIsos         = GaA;
   fNISTElementDataTable[30].fWIsos         = GaW;
   fNISTElementDataTable[30].fMassIsos      = GaIsoMass;


  // Z = 32-------------------------------------------------------------------------------
   static const std::string  GeS[] = {  "Ge58",   "Ge59",   "Ge60",   "Ge61",   "Ge62",
                                        "Ge63",   "Ge64",   "Ge65",   "Ge66",   "Ge67",
                                        "Ge68",   "Ge69",   "Ge70",   "Ge71",   "Ge72",
                                        "Ge73",   "Ge74",   "Ge75",   "Ge76",   "Ge77",
                                        "Ge78",   "Ge79",   "Ge80",   "Ge81",   "Ge82",
                                        "Ge83",   "Ge84",   "Ge85",   "Ge86",   "Ge87",
                                        "Ge88",   "Ge89",   "Ge90" };
   static const int     GeN[] = { 58,  59,  60,  61,  62,  63,  64,  65,  66,  67,
                                  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
                                  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,
                                  88,  89,  90 };
   static const double  GeA[] = {       57.99172 * u,        58.98249 * u,
                                        59.97036 * u,        60.96379 * u,
                                        61.95502 * u,       62.949628 * u,
                                      63.9416899 * u,      64.9393681 * u,
                                      65.9338621 * u,      66.9327339 * u,
                                      67.9280953 * u,      68.9279645 * u,
                                     69.92424875 * u,     70.92495233 * u,
                                    71.922075826 * u,    72.923458956 * u,
                                    73.921177761 * u,     74.92285837 * u,
                                    75.921402726 * u,    76.923549843 * u,
                                      77.9228529 * u,        78.92536 * u,
                                      79.9253508 * u,      80.9288329 * u,
                                       81.929774 * u,      82.9345391 * u,
                                      83.9375751 * u,      84.9429697 * u,
                                        85.94658 * u,        86.95268 * u,
                                        87.95691 * u,        88.96379 * u,
                                        89.96863 * u };
   static const double  GeW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0.2057, 0., 0.2745, 0.0775, 0.365, 0., 0.0773, 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0. };
   static double GeIsoMass[33];

   fNISTElementDataTable[31].fElemSymbol    = "Ge";
   fNISTElementDataTable[31].fSymbols       = GeS;
   fNISTElementDataTable[31].fZ             = 32;
   fNISTElementDataTable[31].fNumOfIsotopes = 33;
   fNISTElementDataTable[31].fIsNoStableIsotope       = false;
   fNISTElementDataTable[31].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[31].fNIsos         = GeN;
   fNISTElementDataTable[31].fAIsos         = GeA;
   fNISTElementDataTable[31].fWIsos         = GeW;
   fNISTElementDataTable[31].fMassIsos      = GeIsoMass;


  // Z = 33-------------------------------------------------------------------------------
   static const std::string  AsS[] = {  "As60",   "As61",   "As62",   "As63",   "As64",
                                        "As65",   "As66",   "As67",   "As68",   "As69",
                                        "As70",   "As71",   "As72",   "As73",   "As74",
                                        "As75",   "As76",   "As77",   "As78",   "As79",
                                        "As80",   "As81",   "As82",   "As83",   "As84",
                                        "As85",   "As86",   "As87",   "As88",   "As89",
                                        "As90",   "As91",   "As92" };
   static const int     AsN[] = { 60,  61,  62,  63,  64,  65,  66,  67,  68,  69,
                                  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
                                  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,
                                  90,  91,  92 };
   static const double  AsA[] = {       59.99388 * u,        60.98112 * u,
                                        61.97361 * u,         62.9639 * u,
                                        63.95743 * u,       64.949611 * u,
                                      65.9441488 * u,     66.93925111 * u,
                                      67.9367741 * u,       68.932246 * u,
                                       69.930926 * u,      70.9271138 * u,
                                      71.9267523 * u,      72.9238291 * u,
                                      73.9239286 * u,     74.92159457 * u,
                                     75.92239202 * u,      76.9206476 * u,
                                       77.921828 * u,      78.9209484 * u,
                                      79.9224746 * u,      80.9221323 * u,
                                      81.9247412 * u,      82.9252069 * u,
                                      83.9293033 * u,      84.9321637 * u,
                                      85.9367015 * u,      86.9402917 * u,
                                        87.94555 * u,        88.94976 * u,
                                        89.95563 * u,        90.96039 * u,
                                        91.96674 * u };
   static const double  AsW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 1, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0. };
   static double AsIsoMass[33];

   fNISTElementDataTable[32].fElemSymbol    = "As";
   fNISTElementDataTable[32].fSymbols       = AsS;
   fNISTElementDataTable[32].fZ             = 33;
   fNISTElementDataTable[32].fNumOfIsotopes = 33;
   fNISTElementDataTable[32].fIsNoStableIsotope       = false;
   fNISTElementDataTable[32].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[32].fNIsos         = AsN;
   fNISTElementDataTable[32].fAIsos         = AsA;
   fNISTElementDataTable[32].fWIsos         = AsW;
   fNISTElementDataTable[32].fMassIsos      = AsIsoMass;


  // Z = 34-------------------------------------------------------------------------------
   static const std::string  SeS[] = {  "Se64",   "Se65",   "Se66",   "Se67",   "Se68",
                                        "Se69",   "Se70",   "Se71",   "Se72",   "Se73",
                                        "Se74",   "Se75",   "Se76",   "Se77",   "Se78",
                                        "Se79",   "Se80",   "Se81",   "Se82",   "Se83",
                                        "Se84",   "Se85",   "Se86",   "Se87",   "Se88",
                                        "Se89",   "Se90",   "Se91",   "Se92",   "Se93",
                                        "Se94",   "Se95" };
   static const int     SeN[] = { 64,  65,  66,  67,  68,  69,  70,  71,  72,  73,
                                  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,
                                  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,
                                  94,  95 };
   static const double  SeA[] = {       63.97109 * u,         64.9644 * u,
                                        65.95559 * u,       66.949994 * u,
                                     67.94182524 * u,      68.9394148 * u,
                                      69.9335155 * u,      70.9322094 * u,
                                      71.9271405 * u,      72.9267549 * u,
                                    73.922475934 * u,     74.92252287 * u,
                                    75.919213704 * u,    76.919914154 * u,
                                     77.91730928 * u,     78.91849929 * u,
                                      79.9165218 * u,       80.917993 * u,
                                      81.9166995 * u,      82.9191186 * u,
                                      83.9184668 * u,      84.9222608 * u,
                                      85.9243117 * u,      86.9286886 * u,
                                      87.9314175 * u,      88.9366691 * u,
                                         89.9401 * u,        90.94596 * u,
                                        91.94984 * u,        92.95629 * u,
                                        93.96049 * u,         94.9673 * u };
   static const double  SeW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0.0089, 0., 0.0937, 0.0763, 0.2377, 0., 0.4961, 0., 0.0873, 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double SeIsoMass[32];

   fNISTElementDataTable[33].fElemSymbol    = "Se";
   fNISTElementDataTable[33].fSymbols       = SeS;
   fNISTElementDataTable[33].fZ             = 34;
   fNISTElementDataTable[33].fNumOfIsotopes = 32;
   fNISTElementDataTable[33].fIsNoStableIsotope       = false;
   fNISTElementDataTable[33].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[33].fNIsos         = SeN;
   fNISTElementDataTable[33].fAIsos         = SeA;
   fNISTElementDataTable[33].fWIsos         = SeW;
   fNISTElementDataTable[33].fMassIsos      = SeIsoMass;


  // Z = 35-------------------------------------------------------------------------------
   static const std::string  BrS[] = {  "Br67",   "Br68",   "Br69",   "Br70",   "Br71",
                                        "Br72",   "Br73",   "Br74",   "Br75",   "Br76",
                                        "Br77",   "Br78",   "Br79",   "Br80",   "Br81",
                                        "Br82",   "Br83",   "Br84",   "Br85",   "Br86",
                                        "Br87",   "Br88",   "Br89",   "Br90",   "Br91",
                                        "Br92",   "Br93",   "Br94",   "Br95",   "Br96",
                                        "Br97",   "Br98" };
   static const int     BrN[] = { 67,  68,  69,  70,  71,  72,  73,  74,  75,  76,
                                  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,
                                  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,
                                  97,  98 };
   static const double  BrA[] = {       66.96465 * u,        67.95873 * u,
                                       68.950497 * u,       69.944792 * u,
                                      70.9393422 * u,      71.9365886 * u,
                                      72.9316715 * u,      73.9299102 * u,
                                      74.9258105 * u,       75.924542 * u,
                                      76.9213792 * u,      77.9211459 * u,
                                      78.9183376 * u,      79.9185298 * u,
                                      80.9162897 * u,      81.9168032 * u,
                                      82.9151756 * u,       83.916496 * u,
                                      84.9156458 * u,      85.9188054 * u,
                                       86.920674 * u,      87.9240833 * u,
                                      88.9267046 * u,      89.9312928 * u,
                                      90.9343986 * u,      91.9396316 * u,
                                        92.94313 * u,         93.9489 * u,
                                        94.95301 * u,        95.95903 * u,
                                        96.96344 * u,        97.96946 * u };
   static const double  BrW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0.5069, 0., 0.4931, 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double BrIsoMass[32];

   fNISTElementDataTable[34].fElemSymbol    = "Br";
   fNISTElementDataTable[34].fSymbols       = BrS;
   fNISTElementDataTable[34].fZ             = 35;
   fNISTElementDataTable[34].fNumOfIsotopes = 32;
   fNISTElementDataTable[34].fIsNoStableIsotope       = false;
   fNISTElementDataTable[34].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[34].fNIsos         = BrN;
   fNISTElementDataTable[34].fAIsos         = BrA;
   fNISTElementDataTable[34].fWIsos         = BrW;
   fNISTElementDataTable[34].fMassIsos      = BrIsoMass;


  // Z = 36-------------------------------------------------------------------------------
   static const std::string  KrS[] = {  "Kr69",   "Kr70",   "Kr71",   "Kr72",   "Kr73",
                                        "Kr74",   "Kr75",   "Kr76",   "Kr77",   "Kr78",
                                        "Kr79",   "Kr80",   "Kr81",   "Kr82",   "Kr83",
                                        "Kr84",   "Kr85",   "Kr86",   "Kr87",   "Kr88",
                                        "Kr89",   "Kr90",   "Kr91",   "Kr92",   "Kr93",
                                        "Kr94",   "Kr95",   "Kr96",   "Kr97",   "Kr98",
                                        "Kr99",  "Kr100",  "Kr101" };
   static const int     KrN[] = { 69,  70,  71,  72,  73,  74,  75,  76,  77,  78,
                                  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,
                                  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,
                                  99, 100, 101 };
   static const double  KrA[] = {       68.96518 * u,        69.95604 * u,
                                        70.95027 * u,      71.9420924 * u,
                                      72.9392892 * u,       73.933084 * u,
                                      74.9309457 * u,      75.9259103 * u,
                                        76.92467 * u,     77.92036494 * u,
                                      78.9200829 * u,     79.91637808 * u,
                                      80.9165912 * u,     81.91348273 * u,
                                     82.91412716 * u,   83.9114977282 * u,
                                      84.9125273 * u,   85.9106106269 * u,
                                     86.91335476 * u,      87.9144479 * u,
                                      88.9178355 * u,      89.9195279 * u,
                                      90.9238063 * u,      91.9261731 * u,
                                      92.9311472 * u,        93.93414 * u,
                                       94.939711 * u,       95.943017 * u,
                                        96.94909 * u,        97.95243 * u,
                                        98.95839 * u,        99.96237 * u,
                                       100.96873 * u };
   static const double  KrW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.00355,
                                  0., 0.02286, 0., 0.11593, 0.115, 0.56987, 0., 0.17279, 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0. };
   static double KrIsoMass[33];

   fNISTElementDataTable[35].fElemSymbol    = "Kr";
   fNISTElementDataTable[35].fSymbols       = KrS;
   fNISTElementDataTable[35].fZ             = 36;
   fNISTElementDataTable[35].fNumOfIsotopes = 33;
   fNISTElementDataTable[35].fIsNoStableIsotope       = false;
   fNISTElementDataTable[35].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[35].fNIsos         = KrN;
   fNISTElementDataTable[35].fAIsos         = KrA;
   fNISTElementDataTable[35].fWIsos         = KrW;
   fNISTElementDataTable[35].fMassIsos      = KrIsoMass;


  // Z = 37-------------------------------------------------------------------------------
   static const std::string  RbS[] = {  "Rb71",   "Rb72",   "Rb73",   "Rb74",   "Rb75",
                                        "Rb76",   "Rb77",   "Rb78",   "Rb79",   "Rb80",
                                        "Rb81",   "Rb82",   "Rb83",   "Rb84",   "Rb85",
                                        "Rb86",   "Rb87",   "Rb88",   "Rb89",   "Rb90",
                                        "Rb91",   "Rb92",   "Rb93",   "Rb94",   "Rb95",
                                        "Rb96",   "Rb97",   "Rb98",   "Rb99",  "Rb100",
                                       "Rb101",  "Rb102",  "Rb103" };
   static const int     RbN[] = { 71,  72,  73,  74,  75,  76,  77,  78,  79,  80,
                                  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
                                  91,  92,  93,  94,  95,  96,  97,  98,  99, 100,
                                 101, 102, 103 };
   static const double  RbA[] = {       70.96532 * u,        71.95908 * u,
                                        72.95053 * u,      73.9442659 * u,
                                      74.9385732 * u,       75.935073 * u,
                                      76.9304016 * u,      77.9281419 * u,
                                      78.9239899 * u,      79.9225164 * u,
                                      80.9189939 * u,       81.918209 * u,
                                      82.9151142 * u,      83.9143752 * u,
                                   84.9117897379 * u,     85.91116743 * u,
                                    86.909180531 * u,     87.91131559 * u,
                                      88.9122783 * u,      89.9147985 * u,
                                      90.9165372 * u,      91.9197284 * u,
                                      92.9220393 * u,      93.9263948 * u,
                                        94.92926 * u,      95.9341334 * u,
                                      96.9371771 * u,      97.9416869 * u,
                                        98.94503 * u,        99.95003 * u,
                                       100.95404 * u,       101.95952 * u,
                                       102.96392 * u };
   static const double  RbW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0.7217, 0., 0.2783, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0. };
   static double RbIsoMass[33];

   fNISTElementDataTable[36].fElemSymbol    = "Rb";
   fNISTElementDataTable[36].fSymbols       = RbS;
   fNISTElementDataTable[36].fZ             = 37;
   fNISTElementDataTable[36].fNumOfIsotopes = 33;
   fNISTElementDataTable[36].fIsNoStableIsotope       = false;
   fNISTElementDataTable[36].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[36].fNIsos         = RbN;
   fNISTElementDataTable[36].fAIsos         = RbA;
   fNISTElementDataTable[36].fWIsos         = RbW;
   fNISTElementDataTable[36].fMassIsos      = RbIsoMass;


  // Z = 38-------------------------------------------------------------------------------
   static const std::string  SrS[] = {  "Sr73",   "Sr74",   "Sr75",   "Sr76",   "Sr77",
                                        "Sr78",   "Sr79",   "Sr80",   "Sr81",   "Sr82",
                                        "Sr83",   "Sr84",   "Sr85",   "Sr86",   "Sr87",
                                        "Sr88",   "Sr89",   "Sr90",   "Sr91",   "Sr92",
                                        "Sr93",   "Sr94",   "Sr95",   "Sr96",   "Sr97",
                                        "Sr98",   "Sr99",  "Sr100",  "Sr101",  "Sr102",
                                       "Sr103",  "Sr104",  "Sr105",  "Sr106",  "Sr107" };
   static const int     SrN[] = { 73,  74,  75,  76,  77,  78,  79,  80,  81,  82,
                                  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,
                                  93,  94,  95,  96,  97,  98,  99, 100, 101, 102,
                                 103, 104, 105, 106, 107 };
   static const double  SrA[] = {        72.9657 * u,        73.95617 * u,
                                        74.94995 * u,       75.941763 * u,
                                      76.9379455 * u,        77.93218 * u,
                                      78.9297077 * u,      79.9245175 * u,
                                      80.9232114 * u,      81.9183999 * u,
                                      82.9175544 * u,      83.9134191 * u,
                                       84.912932 * u,      85.9092606 * u,
                                      86.9088775 * u,      87.9056125 * u,
                                      88.9074511 * u,        89.90773 * u,
                                      90.9101954 * u,      91.9110382 * u,
                                      92.9140242 * u,      93.9153556 * u,
                                      94.9193529 * u,      95.9217066 * u,
                                       96.926374 * u,      97.9286888 * u,
                                      98.9328907 * u,        99.93577 * u,
                                      100.940352 * u,      101.943791 * u,
                                       102.94909 * u,       103.95265 * u,
                                       104.95855 * u,       105.96265 * u,
                                       106.96897 * u };
   static const double  SrW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0.0056, 0., 0.0986, 0.07, 0.8258, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double SrIsoMass[35];

   fNISTElementDataTable[37].fElemSymbol    = "Sr";
   fNISTElementDataTable[37].fSymbols       = SrS;
   fNISTElementDataTable[37].fZ             = 38;
   fNISTElementDataTable[37].fNumOfIsotopes = 35;
   fNISTElementDataTable[37].fIsNoStableIsotope       = false;
   fNISTElementDataTable[37].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[37].fNIsos         = SrN;
   fNISTElementDataTable[37].fAIsos         = SrA;
   fNISTElementDataTable[37].fWIsos         = SrW;
   fNISTElementDataTable[37].fMassIsos      = SrIsoMass;


  // Z = 39-------------------------------------------------------------------------------
   static const std::string   YS[] = {   "Y76",    "Y77",    "Y78",    "Y79",    "Y80",
                                         "Y81",    "Y82",    "Y83",    "Y84",    "Y85",
                                         "Y86",    "Y87",    "Y88",    "Y89",    "Y90",
                                         "Y91",    "Y92",    "Y93",    "Y94",    "Y95",
                                         "Y96",    "Y97",    "Y98",    "Y99",   "Y100",
                                        "Y101",   "Y102",   "Y103",   "Y104",   "Y105",
                                        "Y106",   "Y107",   "Y108",   "Y109" };
   static const int      YN[] = { 76,  77,  78,  79,  80,  81,  82,  83,  84,  85,
                                  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
                                  96,  97,  98,  99, 100, 101, 102, 103, 104, 105,
                                 106, 107, 108, 109 };
   static const double   YA[] = {       75.95856 * u,       76.949781 * u,
                                        77.94361 * u,        78.93735 * u,
                                      79.9343561 * u,      80.9294556 * u,
                                      81.9269314 * u,       82.922485 * u,
                                      83.9206721 * u,       84.916433 * u,
                                       85.914886 * u,      86.9108761 * u,
                                      87.9095016 * u,      88.9058403 * u,
                                      89.9071439 * u,      90.9072974 * u,
                                      91.9089451 * u,       92.909578 * u,
                                      93.9115906 * u,      94.9128161 * u,
                                      95.9158968 * u,      96.9182741 * u,
                                      97.9223821 * u,       98.924148 * u,
                                       99.927715 * u,     100.9301477 * u,
                                     101.9343277 * u,      102.937243 * u,
                                       103.94196 * u,       104.94544 * u,
                                       105.95056 * u,       106.95452 * u,
                                       107.95996 * u,       108.96436 * u };
   static const double   YW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 1, 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0. };
   static double YIsoMass[34];

   fNISTElementDataTable[38].fElemSymbol    = "Y";
   fNISTElementDataTable[38].fSymbols       = YS;
   fNISTElementDataTable[38].fZ             = 39;
   fNISTElementDataTable[38].fNumOfIsotopes = 34;
   fNISTElementDataTable[38].fIsNoStableIsotope       = false;
   fNISTElementDataTable[38].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[38].fNIsos         = YN;
   fNISTElementDataTable[38].fAIsos         = YA;
   fNISTElementDataTable[38].fWIsos         = YW;
   fNISTElementDataTable[38].fMassIsos      = YIsoMass;


  // Z = 40-------------------------------------------------------------------------------
   static const std::string  ZrS[] = {  "Zr78",   "Zr79",   "Zr80",   "Zr81",   "Zr82",
                                        "Zr83",   "Zr84",   "Zr85",   "Zr86",   "Zr87",
                                        "Zr88",   "Zr89",   "Zr90",   "Zr91",   "Zr92",
                                        "Zr93",   "Zr94",   "Zr95",   "Zr96",   "Zr97",
                                        "Zr98",   "Zr99",  "Zr100",  "Zr101",  "Zr102",
                                       "Zr103",  "Zr104",  "Zr105",  "Zr106",  "Zr107",
                                       "Zr108",  "Zr109",  "Zr110",  "Zr111",  "Zr112" };
   static const int     ZrN[] = { 78,  79,  80,  81,  82,  83,  84,  85,  86,  87,
                                  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,
                                  98,  99, 100, 101, 102, 103, 104, 105, 106, 107,
                                 108, 109, 110, 111, 112 };
   static const double  ZrA[] = {       77.95566 * u,        78.94948 * u,
                                         79.9404 * u,        80.93731 * u,
                                        81.93135 * u,      82.9292421 * u,
                                      83.9233269 * u,      84.9214444 * u,
                                      85.9162972 * u,       86.914818 * u,
                                      87.9102213 * u,      88.9088814 * u,
                                      89.9046977 * u,      90.9056396 * u,
                                      91.9050347 * u,      92.9064699 * u,
                                      93.9063108 * u,      94.9080385 * u,
                                      95.9082714 * u,      96.9109512 * u,
                                      97.9127289 * u,       98.916667 * u,
                                      99.9180006 * u,      100.921448 * u,
                                     101.9231409 * u,      102.927191 * u,
                                      103.929436 * u,      104.934008 * u,
                                       105.93676 * u,       106.94174 * u,
                                       107.94487 * u,       108.95041 * u,
                                       109.95396 * u,       110.95968 * u,
                                        111.9637 * u };
   static const double  ZrW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0.5145, 0.1122, 0.1715, 0., 0.1738, 0., 0.028, 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double ZrIsoMass[35];

   fNISTElementDataTable[39].fElemSymbol    = "Zr";
   fNISTElementDataTable[39].fSymbols       = ZrS;
   fNISTElementDataTable[39].fZ             = 40;
   fNISTElementDataTable[39].fNumOfIsotopes = 35;
   fNISTElementDataTable[39].fIsNoStableIsotope       = false;
   fNISTElementDataTable[39].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[39].fNIsos         = ZrN;
   fNISTElementDataTable[39].fAIsos         = ZrA;
   fNISTElementDataTable[39].fWIsos         = ZrW;
   fNISTElementDataTable[39].fMassIsos      = ZrIsoMass;


  // Z = 41-------------------------------------------------------------------------------
   static const std::string  NbS[] = {  "Nb81",   "Nb82",   "Nb83",   "Nb84",   "Nb85",
                                        "Nb86",   "Nb87",   "Nb88",   "Nb89",   "Nb90",
                                        "Nb91",   "Nb92",   "Nb93",   "Nb94",   "Nb95",
                                        "Nb96",   "Nb97",   "Nb98",   "Nb99",  "Nb100",
                                       "Nb101",  "Nb102",  "Nb103",  "Nb104",  "Nb105",
                                       "Nb106",  "Nb107",  "Nb108",  "Nb109",  "Nb110",
                                       "Nb111",  "Nb112",  "Nb113",  "Nb114",  "Nb115" };
   static const int     NbN[] = { 81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
                                  91,  92,  93,  94,  95,  96,  97,  98,  99, 100,
                                 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
                                 111, 112, 113, 114, 115 };
   static const double  NbA[] = {        80.9496 * u,        81.94396 * u,
                                        82.93729 * u,        83.93449 * u,
                                      84.9288458 * u,      85.9257828 * u,
                                      86.9206937 * u,       87.918222 * u,
                                       88.913445 * u,      89.9112584 * u,
                                      90.9069897 * u,      91.9071881 * u,
                                       92.906373 * u,      93.9072788 * u,
                                      94.9068324 * u,      95.9080973 * u,
                                      96.9080959 * u,      97.9103265 * u,
                                       98.911613 * u,      99.9143276 * u,
                                     100.9153103 * u,     101.9180772 * u,
                                     102.9194572 * u,     103.9228925 * u,
                                     104.9249465 * u,     105.9289317 * u,
                                     106.9315937 * u,     107.9360748 * u,
                                       108.93922 * u,       109.94403 * u,
                                       110.94753 * u,       111.95247 * u,
                                       112.95651 * u,       113.96201 * u,
                                       114.96634 * u };
   static const double  NbW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 1, 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double NbIsoMass[35];

   fNISTElementDataTable[40].fElemSymbol    = "Nb";
   fNISTElementDataTable[40].fSymbols       = NbS;
   fNISTElementDataTable[40].fZ             = 41;
   fNISTElementDataTable[40].fNumOfIsotopes = 35;
   fNISTElementDataTable[40].fIsNoStableIsotope       = false;
   fNISTElementDataTable[40].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[40].fNIsos         = NbN;
   fNISTElementDataTable[40].fAIsos         = NbA;
   fNISTElementDataTable[40].fWIsos         = NbW;
   fNISTElementDataTable[40].fMassIsos      = NbIsoMass;


  // Z = 42-------------------------------------------------------------------------------
   static const std::string  MoS[] = {  "Mo83",   "Mo84",   "Mo85",   "Mo86",   "Mo87",
                                        "Mo88",   "Mo89",   "Mo90",   "Mo91",   "Mo92",
                                        "Mo93",   "Mo94",   "Mo95",   "Mo96",   "Mo97",
                                        "Mo98",   "Mo99",  "Mo100",  "Mo101",  "Mo102",
                                       "Mo103",  "Mo104",  "Mo105",  "Mo106",  "Mo107",
                                       "Mo108",  "Mo109",  "Mo110",  "Mo111",  "Mo112",
                                       "Mo113",  "Mo114",  "Mo115",  "Mo116",  "Mo117" };
   static const int     MoN[] = { 83,  84,  85,  86,  87,  88,  89,  90,  91,  92,
                                  93,  94,  95,  96,  97,  98,  99, 100, 101, 102,
                                 103, 104, 105, 106, 107, 108, 109, 110, 111, 112,
                                 113, 114, 115, 116, 117 };
   static const double  MoA[] = {       82.94988 * u,        83.94149 * u,
                                       84.938261 * u,      85.9311748 * u,
                                      86.9281962 * u,      87.9219678 * u,
                                      88.9194682 * u,      89.9139309 * u,
                                      90.9117453 * u,     91.90680796 * u,
                                     92.90680958 * u,      93.9050849 * u,
                                     94.90583877 * u,     95.90467612 * u,
                                     96.90601812 * u,     97.90540482 * u,
                                     98.90770851 * u,      99.9074718 * u,
                                     100.9103414 * u,     101.9102834 * u,
                                      102.913079 * u,     103.9137344 * u,
                                      104.916969 * u,      105.918259 * u,
                                      106.922106 * u,      107.924033 * u,
                                      108.928424 * u,      109.930704 * u,
                                      110.935654 * u,       111.93831 * u,
                                       112.94335 * u,       113.94653 * u,
                                       114.95196 * u,       115.95545 * u,
                                       116.96117 * u };
   static const double  MoW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.1453,
                                  0., 0.0915, 0.1584, 0.1667, 0.096, 0.2439, 0., 0.0982, 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double MoIsoMass[35];

   fNISTElementDataTable[41].fElemSymbol    = "Mo";
   fNISTElementDataTable[41].fSymbols       = MoS;
   fNISTElementDataTable[41].fZ             = 42;
   fNISTElementDataTable[41].fNumOfIsotopes = 35;
   fNISTElementDataTable[41].fIsNoStableIsotope       = false;
   fNISTElementDataTable[41].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[41].fNIsos         = MoN;
   fNISTElementDataTable[41].fAIsos         = MoA;
   fNISTElementDataTable[41].fWIsos         = MoW;
   fNISTElementDataTable[41].fMassIsos      = MoIsoMass;


  // Z = 43-------------------------------------------------------------------------------
   static const std::string  TcS[] = {  "Tc85",   "Tc86",   "Tc87",   "Tc88",   "Tc89",
                                        "Tc90",   "Tc91",   "Tc92",   "Tc93",   "Tc94",
                                        "Tc95",   "Tc96",   "Tc97",   "Tc98",   "Tc99",
                                       "Tc100",  "Tc101",  "Tc102",  "Tc103",  "Tc104",
                                       "Tc105",  "Tc106",  "Tc107",  "Tc108",  "Tc109",
                                       "Tc110",  "Tc111",  "Tc112",  "Tc113",  "Tc114",
                                       "Tc115",  "Tc116",  "Tc117",  "Tc118",  "Tc119",
                                       "Tc120" };
   static const int     TcN[] = { 85,  86,  87,  88,  89,  90,  91,  92,  93,  94,
                                  95,  96,  97,  98,  99, 100, 101, 102, 103, 104,
                                 105, 106, 107, 108, 109, 110, 111, 112, 113, 114,
                                 115, 116, 117, 118, 119, 120 };
   static const double  TcA[] = {       84.95058 * u,        85.94493 * u,
                                      86.9380672 * u,        87.93378 * u,
                                      88.9276487 * u,      89.9240739 * u,
                                      90.9184254 * u,      91.9152698 * u,
                                       92.910246 * u,      93.9096536 * u,
                                      94.9076536 * u,       95.907868 * u,
                                      96.9063667 * u,      97.9072124 * u,
                                      98.9062508 * u,      99.9076539 * u,
                                      100.907309 * u,     101.9092097 * u,
                                      102.909176 * u,      103.911425 * u,
                                      104.911655 * u,      105.914358 * u,
                                     106.9154606 * u,     107.9184957 * u,
                                      108.920256 * u,      109.923744 * u,
                                      110.925901 * u,     111.9299458 * u,
                                      112.932569 * u,       113.93691 * u,
                                       114.93998 * u,       115.94476 * u,
                                       116.94806 * u,       117.95299 * u,
                                       118.95666 * u,       119.96187 * u };
   static const double  TcW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double TcIsoMass[36];

   fNISTElementDataTable[42].fElemSymbol    = "Tc";
   fNISTElementDataTable[42].fSymbols       = TcS;
   fNISTElementDataTable[42].fZ             = 43;
   fNISTElementDataTable[42].fNumOfIsotopes = 36;
   fNISTElementDataTable[42].fIsNoStableIsotope       = true;
   fNISTElementDataTable[42].fIndxOfMostStableIsotope = 13;
   fNISTElementDataTable[42].fNIsos         = TcN;
   fNISTElementDataTable[42].fAIsos         = TcA;
   fNISTElementDataTable[42].fWIsos         = TcW;
   fNISTElementDataTable[42].fMassIsos      = TcIsoMass;


  // Z = 44-------------------------------------------------------------------------------
   static const std::string  RuS[] = {  "Ru87",   "Ru88",   "Ru89",   "Ru90",   "Ru91",
                                        "Ru92",   "Ru93",   "Ru94",   "Ru95",   "Ru96",
                                        "Ru97",   "Ru98",   "Ru99",  "Ru100",  "Ru101",
                                       "Ru102",  "Ru103",  "Ru104",  "Ru105",  "Ru106",
                                       "Ru107",  "Ru108",  "Ru109",  "Ru110",  "Ru111",
                                       "Ru112",  "Ru113",  "Ru114",  "Ru115",  "Ru116",
                                       "Ru117",  "Ru118",  "Ru119",  "Ru120",  "Ru121",
                                       "Ru122",  "Ru123",  "Ru124" };
   static const int     RuN[] = { 87,  88,  89,  90,  91,  92,  93,  94,  95,  96,
                                  97,  98,  99, 100, 101, 102, 103, 104, 105, 106,
                                 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
                                 117, 118, 119, 120, 121, 122, 123, 124 };
   static const double  RuA[] = {       86.95069 * u,         87.9416 * u,
                                        88.93762 * u,      89.9303444 * u,
                                      90.9267419 * u,      91.9202344 * u,
                                      92.9171044 * u,      93.9113429 * u,
                                       94.910406 * u,     95.90759025 * u,
                                      96.9075471 * u,      97.9052868 * u,
                                      98.9059341 * u,      99.9042143 * u,
                                     100.9055769 * u,     101.9043441 * u,
                                     102.9063186 * u,     103.9054275 * u,
                                     104.9077476 * u,     105.9073291 * u,
                                      106.909972 * u,      107.910188 * u,
                                      108.913326 * u,     109.9140407 * u,
                                       110.91757 * u,      111.918809 * u,
                                      112.922844 * u,     113.9246136 * u,
                                       114.92882 * u,     115.9312192 * u,
                                        116.9361 * u,       117.93853 * u,
                                       118.94357 * u,       119.94631 * u,
                                       120.95164 * u,       121.95447 * u,
                                       122.95989 * u,       123.96305 * u };
   static const double  RuW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0554,
                                  0., 0.0187, 0.1276, 0.126, 0.1706, 0.3155, 0., 0.1862, 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double RuIsoMass[38];

   fNISTElementDataTable[43].fElemSymbol    = "Ru";
   fNISTElementDataTable[43].fSymbols       = RuS;
   fNISTElementDataTable[43].fZ             = 44;
   fNISTElementDataTable[43].fNumOfIsotopes = 38;
   fNISTElementDataTable[43].fIsNoStableIsotope       = false;
   fNISTElementDataTable[43].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[43].fNIsos         = RuN;
   fNISTElementDataTable[43].fAIsos         = RuA;
   fNISTElementDataTable[43].fWIsos         = RuW;
   fNISTElementDataTable[43].fMassIsos      = RuIsoMass;


  // Z = 45-------------------------------------------------------------------------------
   static const std::string  RhS[] = {  "Rh89",   "Rh90",   "Rh91",   "Rh92",   "Rh93",
                                        "Rh94",   "Rh95",   "Rh96",   "Rh97",   "Rh98",
                                        "Rh99",  "Rh100",  "Rh101",  "Rh102",  "Rh103",
                                       "Rh104",  "Rh105",  "Rh106",  "Rh107",  "Rh108",
                                       "Rh109",  "Rh110",  "Rh111",  "Rh112",  "Rh113",
                                       "Rh114",  "Rh115",  "Rh116",  "Rh117",  "Rh118",
                                       "Rh119",  "Rh120",  "Rh121",  "Rh122",  "Rh123",
                                       "Rh124",  "Rh125",  "Rh126" };
   static const int     RhN[] = { 89,  90,  91,  92,  93,  94,  95,  96,  97,  98,
                                  99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
                                 109, 110, 111, 112, 113, 114, 115, 116, 117, 118,
                                 119, 120, 121, 122, 123, 124, 125, 126 };
   static const double  RhA[] = {       88.95058 * u,        89.94422 * u,
                                        90.93688 * u,      91.9323677 * u,
                                      92.9259128 * u,      93.9217305 * u,
                                      94.9158979 * u,       95.914453 * u,
                                       96.911329 * u,       97.910708 * u,
                                      98.9081282 * u,       99.908117 * u,
                                     100.9061606 * u,     101.9068374 * u,
                                      102.905498 * u,     103.9066492 * u,
                                     104.9056885 * u,     105.9072868 * u,
                                      106.906748 * u,      107.908714 * u,
                                     108.9087488 * u,      109.911079 * u,
                                     110.9116423 * u,      111.914403 * u,
                                     112.9154393 * u,      113.918718 * u,
                                     114.9203116 * u,      115.924059 * u,
                                     116.9260354 * u,       117.93034 * u,
                                      118.932557 * u,       119.93686 * u,
                                       120.93942 * u,       121.94399 * u,
                                       122.94685 * u,       123.95151 * u,
                                       124.95469 * u,       125.95946 * u };
   static const double  RhW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 1, 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double RhIsoMass[38];

   fNISTElementDataTable[44].fElemSymbol    = "Rh";
   fNISTElementDataTable[44].fSymbols       = RhS;
   fNISTElementDataTable[44].fZ             = 45;
   fNISTElementDataTable[44].fNumOfIsotopes = 38;
   fNISTElementDataTable[44].fIsNoStableIsotope       = false;
   fNISTElementDataTable[44].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[44].fNIsos         = RhN;
   fNISTElementDataTable[44].fAIsos         = RhA;
   fNISTElementDataTable[44].fWIsos         = RhW;
   fNISTElementDataTable[44].fMassIsos      = RhIsoMass;


  // Z = 46-------------------------------------------------------------------------------
   static const std::string  PdS[] = {  "Pd91",   "Pd92",   "Pd93",   "Pd94",   "Pd95",
                                        "Pd96",   "Pd97",   "Pd98",   "Pd99",  "Pd100",
                                       "Pd101",  "Pd102",  "Pd103",  "Pd104",  "Pd105",
                                       "Pd106",  "Pd107",  "Pd108",  "Pd109",  "Pd110",
                                       "Pd111",  "Pd112",  "Pd113",  "Pd114",  "Pd115",
                                       "Pd116",  "Pd117",  "Pd118",  "Pd119",  "Pd120",
                                       "Pd121",  "Pd122",  "Pd123",  "Pd124",  "Pd125",
                                       "Pd126",  "Pd127",  "Pd128" };
   static const int     PdN[] = { 91,  92,  93,  94,  95,  96,  97,  98,  99, 100,
                                 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
                                 111, 112, 113, 114, 115, 116, 117, 118, 119, 120,
                                 121, 122, 123, 124, 125, 126, 127, 128 };
   static const double  PdA[] = {       90.95032 * u,        91.94088 * u,
                                        92.93651 * u,      93.9290376 * u,
                                      94.9248898 * u,      95.9182151 * u,
                                       96.916472 * u,      97.9126983 * u,
                                      98.9117748 * u,       99.908505 * u,
                                     100.9082864 * u,     101.9056022 * u,
                                     102.9060809 * u,     103.9040305 * u,
                                     104.9050796 * u,     105.9034804 * u,
                                     106.9051282 * u,     107.9038916 * u,
                                     108.9059504 * u,     109.9051722 * u,
                                    110.90768968 * u,     111.9073297 * u,
                                      112.910261 * u,     113.9103686 * u,
                                      114.913659 * u,      115.914297 * u,
                                     116.9179547 * u,     117.9190667 * u,
                                     118.9233402 * u,     119.9245511 * u,
                                     120.9289503 * u,      121.930632 * u,
                                       122.93514 * u,       123.93714 * u,
                                       124.94179 * u,       125.94416 * u,
                                       126.94907 * u,       127.95183 * u };
   static const double  PdW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0.0102, 0., 0.1114, 0.2233, 0.2733, 0., 0.2646, 0., 0.1172,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double PdIsoMass[38];

   fNISTElementDataTable[45].fElemSymbol    = "Pd";
   fNISTElementDataTable[45].fSymbols       = PdS;
   fNISTElementDataTable[45].fZ             = 46;
   fNISTElementDataTable[45].fNumOfIsotopes = 38;
   fNISTElementDataTable[45].fIsNoStableIsotope       = false;
   fNISTElementDataTable[45].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[45].fNIsos         = PdN;
   fNISTElementDataTable[45].fAIsos         = PdA;
   fNISTElementDataTable[45].fWIsos         = PdW;
   fNISTElementDataTable[45].fMassIsos      = PdIsoMass;


  // Z = 47-------------------------------------------------------------------------------
   static const std::string  AgS[] = {  "Ag93",   "Ag94",   "Ag95",   "Ag96",   "Ag97",
                                        "Ag98",   "Ag99",  "Ag100",  "Ag101",  "Ag102",
                                       "Ag103",  "Ag104",  "Ag105",  "Ag106",  "Ag107",
                                       "Ag108",  "Ag109",  "Ag110",  "Ag111",  "Ag112",
                                       "Ag113",  "Ag114",  "Ag115",  "Ag116",  "Ag117",
                                       "Ag118",  "Ag119",  "Ag120",  "Ag121",  "Ag122",
                                       "Ag123",  "Ag124",  "Ag125",  "Ag126",  "Ag127",
                                       "Ag128",  "Ag129",  "Ag130" };
   static const int     AgN[] = { 93,  94,  95,  96,  97,  98,  99, 100, 101, 102,
                                 103, 104, 105, 106, 107, 108, 109, 110, 111, 112,
                                 113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
                                 123, 124, 125, 126, 127, 128, 129, 130 };
   static const double  AgA[] = {       92.95033 * u,        93.94373 * u,
                                        94.93602 * u,       95.930744 * u,
                                        96.92397 * u,        97.92156 * u,
                                      98.9176458 * u,      99.9161154 * u,
                                      100.912684 * u,     101.9117047 * u,
                                     102.9089631 * u,     103.9086239 * u,
                                     104.9065256 * u,     105.9066636 * u,
                                     106.9050916 * u,     107.9059503 * u,
                                     108.9047553 * u,     109.9061102 * u,
                                     110.9052959 * u,     111.9070486 * u,
                                      112.906573 * u,      113.908823 * u,
                                      114.908767 * u,     115.9113868 * u,
                                      116.911774 * u,     117.9145955 * u,
                                       118.91557 * u,     119.9187848 * u,
                                      120.920125 * u,      121.923664 * u,
                                      122.925337 * u,       123.92893 * u,
                                       124.93105 * u,       125.93475 * u,
                                       126.93711 * u,       127.94106 * u,
                                       128.94395 * u,        129.9507 * u };
   static const double  AgW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0.51839, 0., 0.48161, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double AgIsoMass[38];

   fNISTElementDataTable[46].fElemSymbol    = "Ag";
   fNISTElementDataTable[46].fSymbols       = AgS;
   fNISTElementDataTable[46].fZ             = 47;
   fNISTElementDataTable[46].fNumOfIsotopes = 38;
   fNISTElementDataTable[46].fIsNoStableIsotope       = false;
   fNISTElementDataTable[46].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[46].fNIsos         = AgN;
   fNISTElementDataTable[46].fAIsos         = AgA;
   fNISTElementDataTable[46].fWIsos         = AgW;
   fNISTElementDataTable[46].fMassIsos      = AgIsoMass;


  // Z = 48-------------------------------------------------------------------------------
   static const std::string  CdS[] = {  "Cd95",   "Cd96",   "Cd97",   "Cd98",   "Cd99",
                                       "Cd100",  "Cd101",  "Cd102",  "Cd103",  "Cd104",
                                       "Cd105",  "Cd106",  "Cd107",  "Cd108",  "Cd109",
                                       "Cd110",  "Cd111",  "Cd112",  "Cd113",  "Cd114",
                                       "Cd115",  "Cd116",  "Cd117",  "Cd118",  "Cd119",
                                       "Cd120",  "Cd121",  "Cd122",  "Cd123",  "Cd124",
                                       "Cd125",  "Cd126",  "Cd127",  "Cd128",  "Cd129",
                                       "Cd130",  "Cd131",  "Cd132",  "Cd133" };
   static const int     CdN[] = { 95,  96,  97,  98,  99, 100, 101, 102, 103, 104,
                                 105, 106, 107, 108, 109, 110, 111, 112, 113, 114,
                                 115, 116, 117, 118, 119, 120, 121, 122, 123, 124,
                                 125, 126, 127, 128, 129, 130, 131, 132, 133 };
   static const double  CdA[] = {       94.94994 * u,        95.94034 * u,
                                         96.9351 * u,       97.927389 * u,
                                      98.9249258 * u,      99.9203488 * u,
                                     100.9185862 * u,      101.914482 * u,
                                     102.9134165 * u,     103.9098564 * u,
                                     104.9094639 * u,     105.9064599 * u,
                                     106.9066121 * u,     107.9041834 * u,
                                     108.9049867 * u,    109.90300661 * u,
                                    110.90418287 * u,    111.90276287 * u,
                                    112.90440813 * u,    113.90336509 * u,
                                    114.90543751 * u,    115.90476315 * u,
                                      116.907226 * u,      117.906922 * u,
                                      118.909847 * u,     119.9098681 * u,
                                     120.9129637 * u,     121.9134591 * u,
                                     122.9168925 * u,     123.9176574 * u,
                                     124.9212576 * u,     125.9224291 * u,
                                      126.926472 * u,     127.9278129 * u,
                                       128.93182 * u,       129.93394 * u,
                                        130.9406 * u,       131.94604 * u,
                                       132.95285 * u };
   static const double  CdW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0.0125, 0., 0.0089, 0., 0.1249, 0.128, 0.2413, 0.1222, 0.2873,
                                  0., 0.0749, 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double CdIsoMass[39];

   fNISTElementDataTable[47].fElemSymbol    = "Cd";
   fNISTElementDataTable[47].fSymbols       = CdS;
   fNISTElementDataTable[47].fZ             = 48;
   fNISTElementDataTable[47].fNumOfIsotopes = 39;
   fNISTElementDataTable[47].fIsNoStableIsotope       = false;
   fNISTElementDataTable[47].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[47].fNIsos         = CdN;
   fNISTElementDataTable[47].fAIsos         = CdA;
   fNISTElementDataTable[47].fWIsos         = CdW;
   fNISTElementDataTable[47].fMassIsos      = CdIsoMass;


  // Z = 49-------------------------------------------------------------------------------
   static const std::string  InS[] = {  "In97",   "In98",   "In99",  "In100",  "In101",
                                       "In102",  "In103",  "In104",  "In105",  "In106",
                                       "In107",  "In108",  "In109",  "In110",  "In111",
                                       "In112",  "In113",  "In114",  "In115",  "In116",
                                       "In117",  "In118",  "In119",  "In120",  "In121",
                                       "In122",  "In123",  "In124",  "In125",  "In126",
                                       "In127",  "In128",  "In129",  "In130",  "In131",
                                       "In132",  "In133",  "In134",  "In135" };
   static const int     InN[] = { 97,  98,  99, 100, 101, 102, 103, 104, 105, 106,
                                 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
                                 117, 118, 119, 120, 121, 122, 123, 124, 125, 126,
                                 127, 128, 129, 130, 131, 132, 133, 134, 135 };
   static const double  InA[] = {       96.94934 * u,        97.94214 * u,
                                        98.93411 * u,        99.93096 * u,
                                       100.92634 * u,     101.9241071 * u,
                                     102.9198819 * u,     103.9182145 * u,
                                      104.914502 * u,      105.913464 * u,
                                       106.91029 * u,     107.9096935 * u,
                                     108.9071514 * u,       109.90717 * u,
                                     110.9051085 * u,     111.9055377 * u,
                                    112.90406184 * u,    113.90491791 * u,
                                   114.903878776 * u,    115.90525999 * u,
                                     116.9045157 * u,     117.9063566 * u,
                                     118.9058507 * u,      119.907967 * u,
                                      120.907851 * u,      121.910281 * u,
                                      122.910434 * u,      123.913182 * u,
                                      124.913605 * u,      125.916507 * u,
                                      126.917446 * u,        127.9204 * u,
                                     128.9218053 * u,      129.924977 * u,
                                     130.9269715 * u,      131.933001 * u,
                                       132.93831 * u,       133.94454 * u,
                                       134.95005 * u };
   static const double  InW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0.0429, 0., 0.9571, 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double InIsoMass[39];

   fNISTElementDataTable[48].fElemSymbol    = "In";
   fNISTElementDataTable[48].fSymbols       = InS;
   fNISTElementDataTable[48].fZ             = 49;
   fNISTElementDataTable[48].fNumOfIsotopes = 39;
   fNISTElementDataTable[48].fIsNoStableIsotope       = false;
   fNISTElementDataTable[48].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[48].fNIsos         = InN;
   fNISTElementDataTable[48].fAIsos         = InA;
   fNISTElementDataTable[48].fWIsos         = InW;
   fNISTElementDataTable[48].fMassIsos      = InIsoMass;


  // Z = 50-------------------------------------------------------------------------------
   static const std::string  SnS[] = {  "Sn99",  "Sn100",  "Sn101",  "Sn102",  "Sn103",
                                       "Sn104",  "Sn105",  "Sn106",  "Sn107",  "Sn108",
                                       "Sn109",  "Sn110",  "Sn111",  "Sn112",  "Sn113",
                                       "Sn114",  "Sn115",  "Sn116",  "Sn117",  "Sn118",
                                       "Sn119",  "Sn120",  "Sn121",  "Sn122",  "Sn123",
                                       "Sn124",  "Sn125",  "Sn126",  "Sn127",  "Sn128",
                                       "Sn129",  "Sn130",  "Sn131",  "Sn132",  "Sn133",
                                       "Sn134",  "Sn135",  "Sn136",  "Sn137",  "Sn138" };
   static const int     SnN[] = { 99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
                                 109, 110, 111, 112, 113, 114, 115, 116, 117, 118,
                                 119, 120, 121, 122, 123, 124, 125, 126, 127, 128,
                                 129, 130, 131, 132, 133, 134, 135, 136, 137, 138 };
   static const double  SnA[] = {       98.94853 * u,         99.9385 * u,
                                       100.93526 * u,       101.93029 * u,
                                      102.928105 * u,     103.9231052 * u,
                                     104.9212684 * u,     105.9169574 * u,
                                     106.9157137 * u,     107.9118943 * u,
                                     108.9112921 * u,      109.907845 * u,
                                     110.9077401 * u,    111.90482387 * u,
                                     112.9051757 * u,     113.9027827 * u,
                                   114.903344699 * u,     115.9017428 * u,
                                    116.90295398 * u,    117.90160657 * u,
                                    118.90331117 * u,    119.90220163 * u,
                                     120.9042426 * u,     121.9034438 * u,
                                     122.9057252 * u,     123.9052766 * u,
                                     124.9077864 * u,      125.907659 * u,
                                       126.91039 * u,      127.910507 * u,
                                      128.913465 * u,     129.9139738 * u,
                                      130.917045 * u,     131.9178267 * u,
                                     132.9239134 * u,     133.9286821 * u,
                                     134.9349086 * u,       135.93999 * u,
                                       136.94655 * u,       137.95184 * u };
   static const double  SnW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0.0097, 0., 0.0066, 0.0034, 0.1454, 0.0768, 0.2422,
                                  0.0859, 0.3258, 0., 0.0463, 0., 0.0579, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double SnIsoMass[40];

   fNISTElementDataTable[49].fElemSymbol    = "Sn";
   fNISTElementDataTable[49].fSymbols       = SnS;
   fNISTElementDataTable[49].fZ             = 50;
   fNISTElementDataTable[49].fNumOfIsotopes = 40;
   fNISTElementDataTable[49].fIsNoStableIsotope       = false;
   fNISTElementDataTable[49].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[49].fNIsos         = SnN;
   fNISTElementDataTable[49].fAIsos         = SnA;
   fNISTElementDataTable[49].fWIsos         = SnW;
   fNISTElementDataTable[49].fMassIsos      = SnIsoMass;


  // Z = 51-------------------------------------------------------------------------------
   static const std::string  SbS[] = { "Sb103",  "Sb104",  "Sb105",  "Sb106",  "Sb107",
                                       "Sb108",  "Sb109",  "Sb110",  "Sb111",  "Sb112",
                                       "Sb113",  "Sb114",  "Sb115",  "Sb116",  "Sb117",
                                       "Sb118",  "Sb119",  "Sb120",  "Sb121",  "Sb122",
                                       "Sb123",  "Sb124",  "Sb125",  "Sb126",  "Sb127",
                                       "Sb128",  "Sb129",  "Sb130",  "Sb131",  "Sb132",
                                       "Sb133",  "Sb134",  "Sb135",  "Sb136",  "Sb137",
                                       "Sb138",  "Sb139",  "Sb140" };
   static const int     SbN[] = {103, 104, 105, 106, 107, 108, 109, 110, 111, 112,
                                 113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
                                 123, 124, 125, 126, 127, 128, 129, 130, 131, 132,
                                 133, 134, 135, 136, 137, 138, 139, 140 };
   static const double  SbA[] = {      102.93969 * u,       103.93648 * u,
                                      104.931276 * u,      105.928638 * u,
                                     106.9241506 * u,     107.9222267 * u,
                                     108.9181411 * u,     109.9168543 * u,
                                     110.9132182 * u,        111.9124 * u,
                                      112.909375 * u,       113.90929 * u,
                                      114.906598 * u,     115.9067931 * u,
                                     116.9048415 * u,     117.9055321 * u,
                                     118.9039455 * u,     119.9050794 * u,
                                      120.903812 * u,     121.9051699 * u,
                                     122.9042132 * u,      123.905935 * u,
                                      124.905253 * u,      125.907253 * u,
                                     126.9069243 * u,      127.909146 * u,
                                      128.909147 * u,      129.911662 * u,
                                     130.9119888 * u,     131.9145077 * u,
                                     132.9152732 * u,     133.9205357 * u,
                                     134.9251851 * u,     135.9307459 * u,
                                       136.93555 * u,       137.94145 * u,
                                       138.94655 * u,       139.95283 * u };
   static const double  SbW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0.5721, 0.,
                                  0.4279, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double SbIsoMass[38];

   fNISTElementDataTable[50].fElemSymbol    = "Sb";
   fNISTElementDataTable[50].fSymbols       = SbS;
   fNISTElementDataTable[50].fZ             = 51;
   fNISTElementDataTable[50].fNumOfIsotopes = 38;
   fNISTElementDataTable[50].fIsNoStableIsotope       = false;
   fNISTElementDataTable[50].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[50].fNIsos         = SbN;
   fNISTElementDataTable[50].fAIsos         = SbA;
   fNISTElementDataTable[50].fWIsos         = SbW;
   fNISTElementDataTable[50].fMassIsos      = SbIsoMass;


  // Z = 52-------------------------------------------------------------------------------
   static const std::string  TeS[] = { "Te105",  "Te106",  "Te107",  "Te108",  "Te109",
                                       "Te110",  "Te111",  "Te112",  "Te113",  "Te114",
                                       "Te115",  "Te116",  "Te117",  "Te118",  "Te119",
                                       "Te120",  "Te121",  "Te122",  "Te123",  "Te124",
                                       "Te125",  "Te126",  "Te127",  "Te128",  "Te129",
                                       "Te130",  "Te131",  "Te132",  "Te133",  "Te134",
                                       "Te135",  "Te136",  "Te137",  "Te138",  "Te139",
                                       "Te140",  "Te141",  "Te142",  "Te143" };
   static const int     TeN[] = {105, 106, 107, 108, 109, 110, 111, 112, 113, 114,
                                 115, 116, 117, 118, 119, 120, 121, 122, 123, 124,
                                 125, 126, 127, 128, 129, 130, 131, 132, 133, 134,
                                 135, 136, 137, 138, 139, 140, 141, 142, 143 };
   static const double  TeA[] = {       104.9433 * u,        105.9375 * u,
                                      106.935012 * u,     107.9293805 * u,
                                     108.9273045 * u,     109.9224581 * u,
                                     110.9210006 * u,     111.9167279 * u,
                                      112.915891 * u,      113.912089 * u,
                                      114.911902 * u,       115.90846 * u,
                                      116.908646 * u,      117.905854 * u,
                                     118.9064071 * u,     119.9040593 * u,
                                      120.904944 * u,     121.9030435 * u,
                                     122.9042698 * u,     123.9028171 * u,
                                     124.9044299 * u,     125.9033109 * u,
                                     126.9052257 * u,    127.90446128 * u,
                                    128.90659646 * u,   129.906222748 * u,
                                   130.908522213 * u,     131.9085467 * u,
                                     132.9109688 * u,      133.911394 * u,
                                     134.9165557 * u,     135.9201006 * u,
                                     136.9255989 * u,     137.9294722 * u,
                                     138.9353672 * u,      139.939499 * u,
                                        140.9458 * u,       141.95022 * u,
                                       142.95676 * u };
   static const double  TeW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0.0009, 0., 0.0255, 0.0089, 0.0474,
                                  0.0707, 0.1884, 0., 0.3174, 0., 0.3408, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double TeIsoMass[39];

   fNISTElementDataTable[51].fElemSymbol    = "Te";
   fNISTElementDataTable[51].fSymbols       = TeS;
   fNISTElementDataTable[51].fZ             = 52;
   fNISTElementDataTable[51].fNumOfIsotopes = 39;
   fNISTElementDataTable[51].fIsNoStableIsotope       = false;
   fNISTElementDataTable[51].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[51].fNIsos         = TeN;
   fNISTElementDataTable[51].fAIsos         = TeA;
   fNISTElementDataTable[51].fWIsos         = TeW;
   fNISTElementDataTable[51].fMassIsos      = TeIsoMass;


  // Z = 53-------------------------------------------------------------------------------
   static const std::string   IS[] = {  "I107",   "I108",   "I109",   "I110",   "I111",
                                        "I112",   "I113",   "I114",   "I115",   "I116",
                                        "I117",   "I118",   "I119",   "I120",   "I121",
                                        "I122",   "I123",   "I124",   "I125",   "I126",
                                        "I127",   "I128",   "I129",   "I130",   "I131",
                                        "I132",   "I133",   "I134",   "I135",   "I136",
                                        "I137",   "I138",   "I139",   "I140",   "I141",
                                        "I142",   "I143",   "I144",   "I145" };
   static const int      IN[] = {107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
                                 117, 118, 119, 120, 121, 122, 123, 124, 125, 126,
                                 127, 128, 129, 130, 131, 132, 133, 134, 135, 136,
                                 137, 138, 139, 140, 141, 142, 143, 144, 145 };
   static const double   IA[] = {      106.94678 * u,       107.94348 * u,
                                     108.9380853 * u,      109.935089 * u,
                                     110.9302692 * u,      111.928005 * u,
                                     112.9236501 * u,       113.92185 * u,
                                      114.918048 * u,       115.91681 * u,
                                      116.913648 * u,      117.913074 * u,
                                      118.910074 * u,      119.910087 * u,
                                     120.9074051 * u,     121.9075888 * u,
                                     122.9055885 * u,      123.906209 * u,
                                     124.9046294 * u,     125.9056233 * u,
                                     126.9044719 * u,     127.9058086 * u,
                                     128.9049837 * u,     129.9066702 * u,
                                     130.9061263 * u,     131.9079935 * u,
                                      132.907797 * u,     133.9097588 * u,
                                     134.9100488 * u,      135.914604 * u,
                                     136.9180282 * u,     137.9227264 * u,
                                      138.926506 * u,       139.93173 * u,
                                       140.93569 * u,        141.9412 * u,
                                       142.94565 * u,       143.95139 * u,
                                       144.95605 * u };
   static const double   IW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  1, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double IIsoMass[39];

   fNISTElementDataTable[52].fElemSymbol    = "I";
   fNISTElementDataTable[52].fSymbols       = IS;
   fNISTElementDataTable[52].fZ             = 53;
   fNISTElementDataTable[52].fNumOfIsotopes = 39;
   fNISTElementDataTable[52].fIsNoStableIsotope       = false;
   fNISTElementDataTable[52].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[52].fNIsos         = IN;
   fNISTElementDataTable[52].fAIsos         = IA;
   fNISTElementDataTable[52].fWIsos         = IW;
   fNISTElementDataTable[52].fMassIsos      = IIsoMass;


  // Z = 54-------------------------------------------------------------------------------
   static const std::string  XeS[] = { "Xe109",  "Xe110",  "Xe111",  "Xe112",  "Xe113",
                                       "Xe114",  "Xe115",  "Xe116",  "Xe117",  "Xe118",
                                       "Xe119",  "Xe120",  "Xe121",  "Xe122",  "Xe123",
                                       "Xe124",  "Xe125",  "Xe126",  "Xe127",  "Xe128",
                                       "Xe129",  "Xe130",  "Xe131",  "Xe132",  "Xe133",
                                       "Xe134",  "Xe135",  "Xe136",  "Xe137",  "Xe138",
                                       "Xe139",  "Xe140",  "Xe141",  "Xe142",  "Xe143",
                                       "Xe144",  "Xe145",  "Xe146",  "Xe147",  "Xe148" };
   static const int     XeN[] = {109, 110, 111, 112, 113, 114, 115, 116, 117, 118,
                                 119, 120, 121, 122, 123, 124, 125, 126, 127, 128,
                                 129, 130, 131, 132, 133, 134, 135, 136, 137, 138,
                                 139, 140, 141, 142, 143, 144, 145, 146, 147, 148 };
   static const double  XeA[] = {      108.95043 * u,       109.94426 * u,
                                      110.941607 * u,      111.935559 * u,
                                     112.9332217 * u,       113.92798 * u,
                                      114.926294 * u,      115.921581 * u,
                                      116.920359 * u,      117.916179 * u,
                                      118.915411 * u,      119.911784 * u,
                                      120.911453 * u,      121.908368 * u,
                                      122.908482 * u,      123.905892 * u,
                                     124.9063944 * u,     125.9042983 * u,
                                     126.9051829 * u,      127.903531 * u,
                                   128.904780861 * u,   129.903509349 * u,
                                    130.90508406 * u,   131.904155086 * u,
                                     132.9059108 * u,    133.90539466 * u,
                                     134.9072278 * u,   135.907214484 * u,
                                    136.91155778 * u,     137.9141463 * u,
                                     138.9187922 * u,     139.9216458 * u,
                                     140.9267872 * u,     141.9299731 * u,
                                     142.9353696 * u,     143.9389451 * u,
                                       144.94472 * u,      145.948518 * u,
                                       146.95426 * u,       147.95813 * u };
   static const double  XeW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0.000952, 0., 0.00089, 0., 0.019102,
                                  0.264006, 0.04071, 0.212324, 0.269086, 0., 0.104357, 0., 0.088573, 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double XeIsoMass[40];

   fNISTElementDataTable[53].fElemSymbol    = "Xe";
   fNISTElementDataTable[53].fSymbols       = XeS;
   fNISTElementDataTable[53].fZ             = 54;
   fNISTElementDataTable[53].fNumOfIsotopes = 40;
   fNISTElementDataTable[53].fIsNoStableIsotope       = false;
   fNISTElementDataTable[53].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[53].fNIsos         = XeN;
   fNISTElementDataTable[53].fAIsos         = XeA;
   fNISTElementDataTable[53].fWIsos         = XeW;
   fNISTElementDataTable[53].fMassIsos      = XeIsoMass;


  // Z = 55-------------------------------------------------------------------------------
   static const std::string  CsS[] = { "Cs112",  "Cs113",  "Cs114",  "Cs115",  "Cs116",
                                       "Cs117",  "Cs118",  "Cs119",  "Cs120",  "Cs121",
                                       "Cs122",  "Cs123",  "Cs124",  "Cs125",  "Cs126",
                                       "Cs127",  "Cs128",  "Cs129",  "Cs130",  "Cs131",
                                       "Cs132",  "Cs133",  "Cs134",  "Cs135",  "Cs136",
                                       "Cs137",  "Cs138",  "Cs139",  "Cs140",  "Cs141",
                                       "Cs142",  "Cs143",  "Cs144",  "Cs145",  "Cs146",
                                       "Cs147",  "Cs148",  "Cs149",  "Cs150",  "Cs151" };
   static const int     CsN[] = {112, 113, 114, 115, 116, 117, 118, 119, 120, 121,
                                 122, 123, 124, 125, 126, 127, 128, 129, 130, 131,
                                 132, 133, 134, 135, 136, 137, 138, 139, 140, 141,
                                 142, 143, 144, 145, 146, 147, 148, 149, 150, 151 };
   static const double  CsA[] = {     111.950309 * u,     112.9444291 * u,
                                      113.941296 * u,       114.93591 * u,
                                       115.93337 * u,      116.928617 * u,
                                       117.92656 * u,      118.922377 * u,
                                      119.920677 * u,      120.917227 * u,
                                      121.916108 * u,      122.912996 * u,
                                     123.9122578 * u,      124.909728 * u,
                                      125.909446 * u,     126.9074174 * u,
                                     127.9077487 * u,     128.9060657 * u,
                                     129.9067093 * u,     130.9054649 * u,
                                     131.9064339 * u,   132.905451961 * u,
                                   133.906718503 * u,      134.905977 * u,
                                     135.9073114 * u,    136.90708923 * u,
                                     137.9110171 * u,     138.9133638 * u,
                                     139.9172831 * u,     140.9200455 * u,
                                      141.924296 * u,      142.927349 * u,
                                      143.932076 * u,      144.935527 * u,
                                      145.940344 * u,      146.944156 * u,
                                       147.94923 * u,       148.95302 * u,
                                       149.95833 * u,       150.96258 * u };
   static const double  CsW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 1, 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double CsIsoMass[40];

   fNISTElementDataTable[54].fElemSymbol    = "Cs";
   fNISTElementDataTable[54].fSymbols       = CsS;
   fNISTElementDataTable[54].fZ             = 55;
   fNISTElementDataTable[54].fNumOfIsotopes = 40;
   fNISTElementDataTable[54].fIsNoStableIsotope       = false;
   fNISTElementDataTable[54].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[54].fNIsos         = CsN;
   fNISTElementDataTable[54].fAIsos         = CsA;
   fNISTElementDataTable[54].fWIsos         = CsW;
   fNISTElementDataTable[54].fMassIsos      = CsIsoMass;


  // Z = 56-------------------------------------------------------------------------------
   static const std::string  BaS[] = { "Ba114",  "Ba115",  "Ba116",  "Ba117",  "Ba118",
                                       "Ba119",  "Ba120",  "Ba121",  "Ba122",  "Ba123",
                                       "Ba124",  "Ba125",  "Ba126",  "Ba127",  "Ba128",
                                       "Ba129",  "Ba130",  "Ba131",  "Ba132",  "Ba133",
                                       "Ba134",  "Ba135",  "Ba136",  "Ba137",  "Ba138",
                                       "Ba139",  "Ba140",  "Ba141",  "Ba142",  "Ba143",
                                       "Ba144",  "Ba145",  "Ba146",  "Ba147",  "Ba148",
                                       "Ba149",  "Ba150",  "Ba151",  "Ba152",  "Ba153" };
   static const int     BaN[] = {114, 115, 116, 117, 118, 119, 120, 121, 122, 123,
                                 124, 125, 126, 127, 128, 129, 130, 131, 132, 133,
                                 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
                                 144, 145, 146, 147, 148, 149, 150, 151, 152, 153 };
   static const double  BaA[] = {      113.95066 * u,       114.94737 * u,
                                       115.94128 * u,       116.93814 * u,
                                       117.93306 * u,       118.93066 * u,
                                       119.92605 * u,       120.92405 * u,
                                      121.919904 * u,      122.918781 * u,
                                      123.915094 * u,      124.914472 * u,
                                       125.91125 * u,      126.911091 * u,
                                      127.908342 * u,      128.908681 * u,
                                     129.9063207 * u,      130.906941 * u,
                                     131.9050611 * u,     132.9060074 * u,
                                    133.90450818 * u,    134.90568838 * u,
                                    135.90457573 * u,    136.90582714 * u,
                                      137.905247 * u,     138.9088411 * u,
                                     139.9106057 * u,     140.9144033 * u,
                                     141.9164324 * u,     142.9206253 * u,
                                     143.9229549 * u,     144.9275184 * u,
                                      145.930284 * u,      146.935304 * u,
                                      147.938171 * u,       148.94308 * u,
                                       149.94605 * u,       150.95127 * u,
                                       151.95481 * u,       152.96036 * u };
   static const double  BaW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0.00106, 0., 0.00101, 0.,
                                  0.02417, 0.06592, 0.07854, 0.11232, 0.71698, 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double BaIsoMass[40];

   fNISTElementDataTable[55].fElemSymbol    = "Ba";
   fNISTElementDataTable[55].fSymbols       = BaS;
   fNISTElementDataTable[55].fZ             = 56;
   fNISTElementDataTable[55].fNumOfIsotopes = 40;
   fNISTElementDataTable[55].fIsNoStableIsotope       = false;
   fNISTElementDataTable[55].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[55].fNIsos         = BaN;
   fNISTElementDataTable[55].fAIsos         = BaA;
   fNISTElementDataTable[55].fWIsos         = BaW;
   fNISTElementDataTable[55].fMassIsos      = BaIsoMass;


  // Z = 57-------------------------------------------------------------------------------
   static const std::string  LaS[] = { "La116",  "La117",  "La118",  "La119",  "La120",
                                       "La121",  "La122",  "La123",  "La124",  "La125",
                                       "La126",  "La127",  "La128",  "La129",  "La130",
                                       "La131",  "La132",  "La133",  "La134",  "La135",
                                       "La136",  "La137",  "La138",  "La139",  "La140",
                                       "La141",  "La142",  "La143",  "La144",  "La145",
                                       "La146",  "La147",  "La148",  "La149",  "La150",
                                       "La151",  "La152",  "La153",  "La154",  "La155" };
   static const int     LaN[] = {116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
                                 126, 127, 128, 129, 130, 131, 132, 133, 134, 135,
                                 136, 137, 138, 139, 140, 141, 142, 143, 144, 145,
                                 146, 147, 148, 149, 150, 151, 152, 153, 154, 155 };
   static const double  LaA[] = {       115.9563 * u,       116.94999 * u,
                                       117.94673 * u,       118.94099 * u,
                                       119.93807 * u,       120.93315 * u,
                                       121.93071 * u,        122.9263 * u,
                                      123.924574 * u,      124.920816 * u,
                                      125.919513 * u,      126.916375 * u,
                                      127.915592 * u,      128.912694 * u,
                                      129.912369 * u,       130.91007 * u,
                                      131.910119 * u,      132.908218 * u,
                                      133.908514 * u,      134.906984 * u,
                                      135.907635 * u,     136.9064504 * u,
                                     137.9071149 * u,     138.9063563 * u,
                                     139.9094806 * u,      140.910966 * u,
                                     141.9140909 * u,     142.9160795 * u,
                                      143.919646 * u,      144.921808 * u,
                                      145.925875 * u,      146.928418 * u,
                                      147.932679 * u,       148.93535 * u,
                                       149.93947 * u,       150.94232 * u,
                                       151.94682 * u,       152.95036 * u,
                                       153.95517 * u,       154.95901 * u };
   static const double  LaW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0.0008881, 0.9991119, 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double LaIsoMass[40];

   fNISTElementDataTable[56].fElemSymbol    = "La";
   fNISTElementDataTable[56].fSymbols       = LaS;
   fNISTElementDataTable[56].fZ             = 57;
   fNISTElementDataTable[56].fNumOfIsotopes = 40;
   fNISTElementDataTable[56].fIsNoStableIsotope       = false;
   fNISTElementDataTable[56].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[56].fNIsos         = LaN;
   fNISTElementDataTable[56].fAIsos         = LaA;
   fNISTElementDataTable[56].fWIsos         = LaW;
   fNISTElementDataTable[56].fMassIsos      = LaIsoMass;


  // Z = 58-------------------------------------------------------------------------------
   static const std::string  CeS[] = { "Ce119",  "Ce120",  "Ce121",  "Ce122",  "Ce123",
                                       "Ce124",  "Ce125",  "Ce126",  "Ce127",  "Ce128",
                                       "Ce129",  "Ce130",  "Ce131",  "Ce132",  "Ce133",
                                       "Ce134",  "Ce135",  "Ce136",  "Ce137",  "Ce138",
                                       "Ce139",  "Ce140",  "Ce141",  "Ce142",  "Ce143",
                                       "Ce144",  "Ce145",  "Ce146",  "Ce147",  "Ce148",
                                       "Ce149",  "Ce150",  "Ce151",  "Ce152",  "Ce153",
                                       "Ce154",  "Ce155",  "Ce156",  "Ce157" };
   static const int     CeN[] = {119, 120, 121, 122, 123, 124, 125, 126, 127, 128,
                                 129, 130, 131, 132, 133, 134, 135, 136, 137, 138,
                                 139, 140, 141, 142, 143, 144, 145, 146, 147, 148,
                                 149, 150, 151, 152, 153, 154, 155, 156, 157 };
   static const double  CeA[] = {      118.95271 * u,       119.94654 * u,
                                       120.94335 * u,       121.93787 * u,
                                       122.93528 * u,       123.93031 * u,
                                       124.92844 * u,      125.923971 * u,
                                      126.922727 * u,      127.918911 * u,
                                      128.918102 * u,      129.914736 * u,
                                      130.914429 * u,      131.911464 * u,
                                       132.91152 * u,      133.908928 * u,
                                      134.909161 * u,    135.90712921 * u,
                                    136.90776236 * u,      137.905991 * u,
                                     138.9066551 * u,     139.9054431 * u,
                                     140.9082807 * u,     141.9092504 * u,
                                     142.9123921 * u,     143.9136529 * u,
                                      144.917265 * u,      145.918802 * u,
                                     146.9226899 * u,      147.924424 * u,
                                      148.928427 * u,      149.930384 * u,
                                      150.934272 * u,        151.9366 * u,
                                       152.94093 * u,        153.9438 * u,
                                       154.94855 * u,       155.95183 * u,
                                       156.95705 * u };
   static const double  CeW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0.00185, 0., 0.00251,
                                  0., 0.8845, 0., 0.11114, 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double CeIsoMass[39];

   fNISTElementDataTable[57].fElemSymbol    = "Ce";
   fNISTElementDataTable[57].fSymbols       = CeS;
   fNISTElementDataTable[57].fZ             = 58;
   fNISTElementDataTable[57].fNumOfIsotopes = 39;
   fNISTElementDataTable[57].fIsNoStableIsotope       = false;
   fNISTElementDataTable[57].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[57].fNIsos         = CeN;
   fNISTElementDataTable[57].fAIsos         = CeA;
   fNISTElementDataTable[57].fWIsos         = CeW;
   fNISTElementDataTable[57].fMassIsos      = CeIsoMass;


  // Z = 59-------------------------------------------------------------------------------
   static const std::string  PrS[] = { "Pr121",  "Pr122",  "Pr123",  "Pr124",  "Pr125",
                                       "Pr126",  "Pr127",  "Pr128",  "Pr129",  "Pr130",
                                       "Pr131",  "Pr132",  "Pr133",  "Pr134",  "Pr135",
                                       "Pr136",  "Pr137",  "Pr138",  "Pr139",  "Pr140",
                                       "Pr141",  "Pr142",  "Pr143",  "Pr144",  "Pr145",
                                       "Pr146",  "Pr147",  "Pr148",  "Pr149",  "Pr150",
                                       "Pr151",  "Pr152",  "Pr153",  "Pr154",  "Pr155",
                                       "Pr156",  "Pr157",  "Pr158",  "Pr159" };
   static const int     PrN[] = {121, 122, 123, 124, 125, 126, 127, 128, 129, 130,
                                 131, 132, 133, 134, 135, 136, 137, 138, 139, 140,
                                 141, 142, 143, 144, 145, 146, 147, 148, 149, 150,
                                 151, 152, 153, 154, 155, 156, 157, 158, 159 };
   static const double  PrA[] = {      120.95532 * u,       121.95175 * u,
                                       122.94596 * u,       123.94294 * u,
                                        124.9377 * u,       125.93524 * u,
                                       126.93071 * u,      127.928791 * u,
                                      128.925095 * u,       129.92359 * u,
                                      130.920235 * u,      131.919255 * u,
                                      132.916331 * u,      133.915697 * u,
                                      134.913112 * u,      135.912677 * u,
                                     136.9106792 * u,      137.910754 * u,
                                     138.9089408 * u,     139.9090803 * u,
                                     140.9076576 * u,     141.9100496 * u,
                                     142.9108228 * u,     143.9133109 * u,
                                     144.9145182 * u,       145.91768 * u,
                                      146.919008 * u,       147.92213 * u,
                                      148.923736 * u,     149.9266765 * u,
                                      150.928309 * u,      151.931553 * u,
                                      152.933904 * u,       153.93753 * u,
                                      154.940509 * u,       155.94464 * u,
                                       156.94789 * u,       157.95241 * u,
                                       158.95589 * u };
   static const double  PrW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  1, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double PrIsoMass[39];

   fNISTElementDataTable[58].fElemSymbol    = "Pr";
   fNISTElementDataTable[58].fSymbols       = PrS;
   fNISTElementDataTable[58].fZ             = 59;
   fNISTElementDataTable[58].fNumOfIsotopes = 39;
   fNISTElementDataTable[58].fIsNoStableIsotope       = false;
   fNISTElementDataTable[58].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[58].fNIsos         = PrN;
   fNISTElementDataTable[58].fAIsos         = PrA;
   fNISTElementDataTable[58].fWIsos         = PrW;
   fNISTElementDataTable[58].fMassIsos      = PrIsoMass;


  // Z = 60-------------------------------------------------------------------------------
   static const std::string  NdS[] = { "Nd124",  "Nd125",  "Nd126",  "Nd127",  "Nd128",
                                       "Nd129",  "Nd130",  "Nd131",  "Nd132",  "Nd133",
                                       "Nd134",  "Nd135",  "Nd136",  "Nd137",  "Nd138",
                                       "Nd139",  "Nd140",  "Nd141",  "Nd142",  "Nd143",
                                       "Nd144",  "Nd145",  "Nd146",  "Nd147",  "Nd148",
                                       "Nd149",  "Nd150",  "Nd151",  "Nd152",  "Nd153",
                                       "Nd154",  "Nd155",  "Nd156",  "Nd157",  "Nd158",
                                       "Nd159",  "Nd160",  "Nd161" };
   static const int     NdN[] = {124, 125, 126, 127, 128, 129, 130, 131, 132, 133,
                                 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
                                 144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
                                 154, 155, 156, 157, 158, 159, 160, 161 };
   static const double  NdA[] = {       123.9522 * u,        124.9489 * u,
                                       125.94311 * u,       126.94038 * u,
                                       127.93525 * u,        128.9331 * u,
                                      129.928506 * u,      130.927248 * u,
                                      131.923321 * u,      132.922348 * u,
                                       133.91879 * u,      134.918181 * u,
                                      135.914976 * u,      136.914562 * u,
                                       137.91195 * u,      138.911954 * u,
                                       139.90955 * u,     140.9096147 * u,
                                      141.907729 * u,       142.90982 * u,
                                      143.910093 * u,     144.9125793 * u,
                                     145.9131226 * u,     146.9161061 * u,
                                     147.9168993 * u,     148.9201548 * u,
                                     149.9209022 * u,     150.9238403 * u,
                                      151.924692 * u,      152.927718 * u,
                                       153.92948 * u,     154.9331357 * u,
                                       155.93508 * u,      156.939386 * u,
                                       157.94197 * u,       158.94653 * u,
                                        159.9494 * u,       160.95428 * u };
   static const double  NdW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0.27152, 0.12174,
                                  0.23798, 0.08293, 0.17189, 0., 0.05756, 0., 0.05638, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double NdIsoMass[38];

   fNISTElementDataTable[59].fElemSymbol    = "Nd";
   fNISTElementDataTable[59].fSymbols       = NdS;
   fNISTElementDataTable[59].fZ             = 60;
   fNISTElementDataTable[59].fNumOfIsotopes = 38;
   fNISTElementDataTable[59].fIsNoStableIsotope       = false;
   fNISTElementDataTable[59].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[59].fNIsos         = NdN;
   fNISTElementDataTable[59].fAIsos         = NdA;
   fNISTElementDataTable[59].fWIsos         = NdW;
   fNISTElementDataTable[59].fMassIsos      = NdIsoMass;


  // Z = 61-------------------------------------------------------------------------------
   static const std::string  PmS[] = { "Pm126",  "Pm127",  "Pm128",  "Pm129",  "Pm130",
                                       "Pm131",  "Pm132",  "Pm133",  "Pm134",  "Pm135",
                                       "Pm136",  "Pm137",  "Pm138",  "Pm139",  "Pm140",
                                       "Pm141",  "Pm142",  "Pm143",  "Pm144",  "Pm145",
                                       "Pm146",  "Pm147",  "Pm148",  "Pm149",  "Pm150",
                                       "Pm151",  "Pm152",  "Pm153",  "Pm154",  "Pm155",
                                       "Pm156",  "Pm157",  "Pm158",  "Pm159",  "Pm160",
                                       "Pm161",  "Pm162",  "Pm163" };
   static const int     PmN[] = {126, 127, 128, 129, 130, 131, 132, 133, 134, 135,
                                 136, 137, 138, 139, 140, 141, 142, 143, 144, 145,
                                 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
                                 156, 157, 158, 159, 160, 161, 162, 163 };
   static const double  PmA[] = {      125.95792 * u,       126.95192 * u,
                                        127.9487 * u,       128.94323 * u,
                                       129.94053 * u,       130.93567 * u,
                                       131.93384 * u,      132.929782 * u,
                                      133.928353 * u,      134.924823 * u,
                                      135.923585 * u,       136.92048 * u,
                                      137.919548 * u,        138.9168 * u,
                                       139.91604 * u,      140.913555 * u,
                                       141.91289 * u,     142.9109383 * u,
                                     143.9125964 * u,     144.9127559 * u,
                                     145.9147024 * u,      146.915145 * u,
                                     147.9174819 * u,     148.9183423 * u,
                                      149.920991 * u,     150.9212175 * u,
                                      151.923506 * u,     152.9241567 * u,
                                      153.926472 * u,      154.928137 * u,
                                     155.9311175 * u,     156.9331214 * u,
                                      157.936565 * u,      158.939287 * u,
                                        159.9431 * u,       160.94607 * u,
                                       161.95022 * u,       162.95357 * u };
   static const double  PmW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double PmIsoMass[38];

   fNISTElementDataTable[60].fElemSymbol    = "Pm";
   fNISTElementDataTable[60].fSymbols       = PmS;
   fNISTElementDataTable[60].fZ             = 61;
   fNISTElementDataTable[60].fNumOfIsotopes = 38;
   fNISTElementDataTable[60].fIsNoStableIsotope       = true;
   fNISTElementDataTable[60].fIndxOfMostStableIsotope = 19;
   fNISTElementDataTable[60].fNIsos         = PmN;
   fNISTElementDataTable[60].fAIsos         = PmA;
   fNISTElementDataTable[60].fWIsos         = PmW;
   fNISTElementDataTable[60].fMassIsos      = PmIsoMass;


  // Z = 62-------------------------------------------------------------------------------
   static const std::string  SmS[] = { "Sm128",  "Sm129",  "Sm130",  "Sm131",  "Sm132",
                                       "Sm133",  "Sm134",  "Sm135",  "Sm136",  "Sm137",
                                       "Sm138",  "Sm139",  "Sm140",  "Sm141",  "Sm142",
                                       "Sm143",  "Sm144",  "Sm145",  "Sm146",  "Sm147",
                                       "Sm148",  "Sm149",  "Sm150",  "Sm151",  "Sm152",
                                       "Sm153",  "Sm154",  "Sm155",  "Sm156",  "Sm157",
                                       "Sm158",  "Sm159",  "Sm160",  "Sm161",  "Sm162",
                                       "Sm163",  "Sm164",  "Sm165" };
   static const int     SmN[] = {128, 129, 130, 131, 132, 133, 134, 135, 136, 137,
                                 138, 139, 140, 141, 142, 143, 144, 145, 146, 147,
                                 148, 149, 150, 151, 152, 153, 154, 155, 156, 157,
                                 158, 159, 160, 161, 162, 163, 164, 165 };
   static const double  SmA[] = {      127.95842 * u,       128.95476 * u,
                                         129.949 * u,       130.94618 * u,
                                       131.94087 * u,       132.93856 * u,
                                       133.93411 * u,       134.93252 * u,
                                      135.928276 * u,      136.926971 * u,
                                      137.923244 * u,      138.922297 * u,
                                      139.918995 * u,     140.9184816 * u,
                                     141.9152044 * u,     142.9146353 * u,
                                     143.9120065 * u,     144.9134173 * u,
                                      145.913047 * u,     146.9149044 * u,
                                     147.9148292 * u,     148.9171921 * u,
                                     149.9172829 * u,     150.9199398 * u,
                                     151.9197397 * u,     152.9221047 * u,
                                     153.9222169 * u,     154.9246477 * u,
                                      155.925536 * u,     156.9284187 * u,
                                      157.929951 * u,     158.9332172 * u,
                                     159.9353353 * u,     160.9391602 * u,
                                       161.94146 * u,       162.94555 * u,
                                       163.94836 * u,       164.95297 * u };
   static const double  SmW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0.0307, 0., 0., 0.1499,
                                  0.1124, 0.1382, 0.0738, 0., 0.2675, 0., 0.2275, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double SmIsoMass[38];

   fNISTElementDataTable[61].fElemSymbol    = "Sm";
   fNISTElementDataTable[61].fSymbols       = SmS;
   fNISTElementDataTable[61].fZ             = 62;
   fNISTElementDataTable[61].fNumOfIsotopes = 38;
   fNISTElementDataTable[61].fIsNoStableIsotope       = false;
   fNISTElementDataTable[61].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[61].fNIsos         = SmN;
   fNISTElementDataTable[61].fAIsos         = SmA;
   fNISTElementDataTable[61].fWIsos         = SmW;
   fNISTElementDataTable[61].fMassIsos      = SmIsoMass;


  // Z = 63-------------------------------------------------------------------------------
   static const std::string  EuS[] = { "Eu130",  "Eu131",  "Eu132",  "Eu133",  "Eu134",
                                       "Eu135",  "Eu136",  "Eu137",  "Eu138",  "Eu139",
                                       "Eu140",  "Eu141",  "Eu142",  "Eu143",  "Eu144",
                                       "Eu145",  "Eu146",  "Eu147",  "Eu148",  "Eu149",
                                       "Eu150",  "Eu151",  "Eu152",  "Eu153",  "Eu154",
                                       "Eu155",  "Eu156",  "Eu157",  "Eu158",  "Eu159",
                                       "Eu160",  "Eu161",  "Eu162",  "Eu163",  "Eu164",
                                       "Eu165",  "Eu166",  "Eu167" };
   static const int     EuN[] = {130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
                                 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
                                 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
                                 160, 161, 162, 163, 164, 165, 166, 167 };
   static const double  EuA[] = {      129.96369 * u,       130.95784 * u,
                                       131.95467 * u,       132.94929 * u,
                                        133.9464 * u,       134.94187 * u,
                                       135.93962 * u,       136.93546 * u,
                                      137.933709 * u,      138.929792 * u,
                                      139.928088 * u,      140.924932 * u,
                                      141.923442 * u,      142.920299 * u,
                                       143.91882 * u,     144.9162726 * u,
                                      145.917211 * u,     146.9167527 * u,
                                      147.918089 * u,     148.9179378 * u,
                                     149.9197077 * u,     150.9198578 * u,
                                     151.9217522 * u,      152.921238 * u,
                                      153.922987 * u,     154.9229011 * u,
                                     155.9247605 * u,     156.9254334 * u,
                                      157.927799 * u,     158.9291001 * u,
                                      159.931851 * u,      160.933664 * u,
                                      161.936989 * u,      162.939196 * u,
                                       163.94274 * u,       164.94559 * u,
                                       165.94962 * u,       166.95289 * u };
   static const double  EuW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0.4781, 0., 0.5219, 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double EuIsoMass[38];

   fNISTElementDataTable[62].fElemSymbol    = "Eu";
   fNISTElementDataTable[62].fSymbols       = EuS;
   fNISTElementDataTable[62].fZ             = 63;
   fNISTElementDataTable[62].fNumOfIsotopes = 38;
   fNISTElementDataTable[62].fIsNoStableIsotope       = false;
   fNISTElementDataTable[62].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[62].fNIsos         = EuN;
   fNISTElementDataTable[62].fAIsos         = EuA;
   fNISTElementDataTable[62].fWIsos         = EuW;
   fNISTElementDataTable[62].fMassIsos      = EuIsoMass;


  // Z = 64-------------------------------------------------------------------------------
   static const std::string  GdS[] = { "Gd133",  "Gd134",  "Gd135",  "Gd136",  "Gd137",
                                       "Gd138",  "Gd139",  "Gd140",  "Gd141",  "Gd142",
                                       "Gd143",  "Gd144",  "Gd145",  "Gd146",  "Gd147",
                                       "Gd148",  "Gd149",  "Gd150",  "Gd151",  "Gd152",
                                       "Gd153",  "Gd154",  "Gd155",  "Gd156",  "Gd157",
                                       "Gd158",  "Gd159",  "Gd160",  "Gd161",  "Gd162",
                                       "Gd163",  "Gd164",  "Gd165",  "Gd166",  "Gd167",
                                       "Gd168",  "Gd169" };
   static const int     GdN[] = {133, 134, 135, 136, 137, 138, 139, 140, 141, 142,
                                 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,
                                 153, 154, 155, 156, 157, 158, 159, 160, 161, 162,
                                 163, 164, 165, 166, 167, 168, 169 };
   static const double  GdA[] = {      132.96133 * u,       133.95566 * u,
                                       134.95245 * u,        135.9473 * u,
                                       136.94502 * u,       137.94025 * u,
                                       138.93813 * u,      139.933674 * u,
                                      140.932126 * u,      141.928116 * u,
                                       142.92675 * u,      143.922963 * u,
                                      144.921713 * u,     145.9183188 * u,
                                     146.9191014 * u,     147.9181215 * u,
                                     148.9193481 * u,     149.9186644 * u,
                                      150.920356 * u,     151.9197995 * u,
                                      152.921758 * u,     153.9208741 * u,
                                     154.9226305 * u,     155.9221312 * u,
                                     156.9239686 * u,     157.9241123 * u,
                                      158.926397 * u,     159.9270624 * u,
                                     160.9296775 * u,      161.930993 * u,
                                     162.9341769 * u,       163.93583 * u,
                                       164.93936 * u,       165.94146 * u,
                                       166.94545 * u,       167.94808 * u,
                                        168.9526 * u };
   static const double  GdW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.002,
                                  0., 0.0218, 0.148, 0.2047, 0.1565, 0.2484, 0., 0.2186, 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0. };
   static double GdIsoMass[37];

   fNISTElementDataTable[63].fElemSymbol    = "Gd";
   fNISTElementDataTable[63].fSymbols       = GdS;
   fNISTElementDataTable[63].fZ             = 64;
   fNISTElementDataTable[63].fNumOfIsotopes = 37;
   fNISTElementDataTable[63].fIsNoStableIsotope       = false;
   fNISTElementDataTable[63].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[63].fNIsos         = GdN;
   fNISTElementDataTable[63].fAIsos         = GdA;
   fNISTElementDataTable[63].fWIsos         = GdW;
   fNISTElementDataTable[63].fMassIsos      = GdIsoMass;


  // Z = 65-------------------------------------------------------------------------------
   static const std::string  TbS[] = { "Tb135",  "Tb136",  "Tb137",  "Tb138",  "Tb139",
                                       "Tb140",  "Tb141",  "Tb142",  "Tb143",  "Tb144",
                                       "Tb145",  "Tb146",  "Tb147",  "Tb148",  "Tb149",
                                       "Tb150",  "Tb151",  "Tb152",  "Tb153",  "Tb154",
                                       "Tb155",  "Tb156",  "Tb157",  "Tb158",  "Tb159",
                                       "Tb160",  "Tb161",  "Tb162",  "Tb163",  "Tb164",
                                       "Tb165",  "Tb166",  "Tb167",  "Tb168",  "Tb169",
                                       "Tb170",  "Tb171" };
   static const int     TbN[] = {135, 136, 137, 138, 139, 140, 141, 142, 143, 144,
                                 145, 146, 147, 148, 149, 150, 151, 152, 153, 154,
                                 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
                                 165, 166, 167, 168, 169, 170, 171 };
   static const double  TbA[] = {      134.96476 * u,       135.96129 * u,
                                       136.95602 * u,       137.95312 * u,
                                       138.94833 * u,       139.94581 * u,
                                       140.94145 * u,       141.93928 * u,
                                      142.935137 * u,      143.933045 * u,
                                       144.92882 * u,      145.927253 * u,
                                     146.9240548 * u,      147.924282 * u,
                                     148.9232535 * u,     149.9236649 * u,
                                     150.9231096 * u,      151.924083 * u,
                                     152.9234424 * u,      153.924685 * u,
                                      154.923511 * u,     155.9247552 * u,
                                      156.924033 * u,     157.9254209 * u,
                                     158.9253547 * u,     159.9271756 * u,
                                     160.9275778 * u,      161.929495 * u,
                                     162.9306547 * u,       163.93336 * u,
                                       164.93498 * u,       165.93786 * u,
                                       166.93996 * u,        167.9434 * u,
                                       168.94597 * u,       169.94984 * u,
                                       170.95273 * u };
   static const double  TbW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 1, 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0. };
   static double TbIsoMass[37];

   fNISTElementDataTable[64].fElemSymbol    = "Tb";
   fNISTElementDataTable[64].fSymbols       = TbS;
   fNISTElementDataTable[64].fZ             = 65;
   fNISTElementDataTable[64].fNumOfIsotopes = 37;
   fNISTElementDataTable[64].fIsNoStableIsotope       = false;
   fNISTElementDataTable[64].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[64].fNIsos         = TbN;
   fNISTElementDataTable[64].fAIsos         = TbA;
   fNISTElementDataTable[64].fWIsos         = TbW;
   fNISTElementDataTable[64].fMassIsos      = TbIsoMass;


  // Z = 66-------------------------------------------------------------------------------
   static const std::string  DyS[] = { "Dy138",  "Dy139",  "Dy140",  "Dy141",  "Dy142",
                                       "Dy143",  "Dy144",  "Dy145",  "Dy146",  "Dy147",
                                       "Dy148",  "Dy149",  "Dy150",  "Dy151",  "Dy152",
                                       "Dy153",  "Dy154",  "Dy155",  "Dy156",  "Dy157",
                                       "Dy158",  "Dy159",  "Dy160",  "Dy161",  "Dy162",
                                       "Dy163",  "Dy164",  "Dy165",  "Dy166",  "Dy167",
                                       "Dy168",  "Dy169",  "Dy170",  "Dy171",  "Dy172",
                                       "Dy173" };
   static const int     DyN[] = {138, 139, 140, 141, 142, 143, 144, 145, 146, 147,
                                 148, 149, 150, 151, 152, 153, 154, 155, 156, 157,
                                 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
                                 168, 169, 170, 171, 172, 173 };
   static const double  DyA[] = {       137.9625 * u,       138.95959 * u,
                                       139.95402 * u,       140.95128 * u,
                                       141.94619 * u,      142.943994 * u,
                                     143.9392695 * u,      144.937474 * u,
                                     145.9328445 * u,     146.9310827 * u,
                                      147.927157 * u,      148.927322 * u,
                                     149.9255933 * u,     150.9261916 * u,
                                     151.9247253 * u,     152.9257724 * u,
                                     153.9244293 * u,      154.925759 * u,
                                     155.9242847 * u,     156.9254707 * u,
                                     157.9244159 * u,      158.925747 * u,
                                     159.9252046 * u,     160.9269405 * u,
                                     161.9268056 * u,     162.9287383 * u,
                                     163.9291819 * u,     164.9317105 * u,
                                     165.9328139 * u,      166.935661 * u,
                                       167.93713 * u,       168.94031 * u,
                                       169.94239 * u,       170.94612 * u,
                                       171.94846 * u,       172.95283 * u };
   static const double  DyW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0.00056, 0.,
                                  0.00095, 0., 0.02329, 0.18889, 0.25475, 0.24896, 0.2826, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double DyIsoMass[36];

   fNISTElementDataTable[65].fElemSymbol    = "Dy";
   fNISTElementDataTable[65].fSymbols       = DyS;
   fNISTElementDataTable[65].fZ             = 66;
   fNISTElementDataTable[65].fNumOfIsotopes = 36;
   fNISTElementDataTable[65].fIsNoStableIsotope       = false;
   fNISTElementDataTable[65].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[65].fNIsos         = DyN;
   fNISTElementDataTable[65].fAIsos         = DyA;
   fNISTElementDataTable[65].fWIsos         = DyW;
   fNISTElementDataTable[65].fMassIsos      = DyIsoMass;


  // Z = 67-------------------------------------------------------------------------------
   static const std::string  HoS[] = { "Ho140",  "Ho141",  "Ho142",  "Ho143",  "Ho144",
                                       "Ho145",  "Ho146",  "Ho147",  "Ho148",  "Ho149",
                                       "Ho150",  "Ho151",  "Ho152",  "Ho153",  "Ho154",
                                       "Ho155",  "Ho156",  "Ho157",  "Ho158",  "Ho159",
                                       "Ho160",  "Ho161",  "Ho162",  "Ho163",  "Ho164",
                                       "Ho165",  "Ho166",  "Ho167",  "Ho168",  "Ho169",
                                       "Ho170",  "Ho171",  "Ho172",  "Ho173",  "Ho174",
                                       "Ho175" };
   static const int     HoN[] = {140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
                                 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
                                 160, 161, 162, 163, 164, 165, 166, 167, 168, 169,
                                 170, 171, 172, 173, 174, 175 };
   static const double  HoA[] = {      139.96859 * u,       140.96311 * u,
                                       141.96001 * u,       142.95486 * u,
                                     143.9521097 * u,     144.9472674 * u,
                                     145.9449935 * u,     146.9401423 * u,
                                      147.937744 * u,      148.933803 * u,
                                      149.933498 * u,     150.9316983 * u,
                                      151.931724 * u,     152.9302064 * u,
                                     153.9306068 * u,      154.929104 * u,
                                      155.929706 * u,      156.928254 * u,
                                      157.928946 * u,     158.9277197 * u,
                                      159.928737 * u,     160.9278615 * u,
                                     161.9291023 * u,      162.928741 * u,
                                     163.9302403 * u,     164.9303288 * u,
                                     165.9322909 * u,     166.9331385 * u,
                                      167.935522 * u,      168.936878 * u,
                                      169.939625 * u,       170.94147 * u,
                                       171.94473 * u,       172.94702 * u,
                                       173.95095 * u,       174.95362 * u };
   static const double  HoW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 1, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double HoIsoMass[36];

   fNISTElementDataTable[66].fElemSymbol    = "Ho";
   fNISTElementDataTable[66].fSymbols       = HoS;
   fNISTElementDataTable[66].fZ             = 67;
   fNISTElementDataTable[66].fNumOfIsotopes = 36;
   fNISTElementDataTable[66].fIsNoStableIsotope       = false;
   fNISTElementDataTable[66].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[66].fNIsos         = HoN;
   fNISTElementDataTable[66].fAIsos         = HoA;
   fNISTElementDataTable[66].fWIsos         = HoW;
   fNISTElementDataTable[66].fMassIsos      = HoIsoMass;


  // Z = 68-------------------------------------------------------------------------------
   static const std::string  ErS[] = { "Er142",  "Er143",  "Er144",  "Er145",  "Er146",
                                       "Er147",  "Er148",  "Er149",  "Er150",  "Er151",
                                       "Er152",  "Er153",  "Er154",  "Er155",  "Er156",
                                       "Er157",  "Er158",  "Er159",  "Er160",  "Er161",
                                       "Er162",  "Er163",  "Er164",  "Er165",  "Er166",
                                       "Er167",  "Er168",  "Er169",  "Er170",  "Er171",
                                       "Er172",  "Er173",  "Er174",  "Er175",  "Er176",
                                       "Er177" };
   static const int     ErN[] = {142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
                                 152, 153, 154, 155, 156, 157, 158, 159, 160, 161,
                                 162, 163, 164, 165, 166, 167, 168, 169, 170, 171,
                                 172, 173, 174, 175, 176, 177 };
   static const double  ErA[] = {       141.9701 * u,       142.96662 * u,
                                        143.9607 * u,       144.95805 * u,
                                     145.9524184 * u,      146.949964 * u,
                                      147.944735 * u,      148.942306 * u,
                                      149.937916 * u,      150.937449 * u,
                                      151.935057 * u,       152.93508 * u,
                                     153.9327908 * u,     154.9332159 * u,
                                      155.931067 * u,      156.931949 * u,
                                      157.929893 * u,     158.9306918 * u,
                                      159.929077 * u,     160.9300046 * u,
                                     161.9287884 * u,     162.9300408 * u,
                                     163.9292088 * u,     164.9307345 * u,
                                     165.9302995 * u,     166.9320546 * u,
                                     167.9323767 * u,     168.9345968 * u,
                                     169.9354702 * u,     170.9380357 * u,
                                     171.9393619 * u,        172.9424 * u,
                                       173.94423 * u,       174.94777 * u,
                                       175.94994 * u,       176.95399 * u };
   static const double  ErW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0.00139, 0., 0.01601, 0., 0.33503, 0.22869, 0.26978, 0., 0.1491, 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double ErIsoMass[36];

   fNISTElementDataTable[67].fElemSymbol    = "Er";
   fNISTElementDataTable[67].fSymbols       = ErS;
   fNISTElementDataTable[67].fZ             = 68;
   fNISTElementDataTable[67].fNumOfIsotopes = 36;
   fNISTElementDataTable[67].fIsNoStableIsotope       = false;
   fNISTElementDataTable[67].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[67].fNIsos         = ErN;
   fNISTElementDataTable[67].fAIsos         = ErA;
   fNISTElementDataTable[67].fWIsos         = ErW;
   fNISTElementDataTable[67].fMassIsos      = ErIsoMass;


  // Z = 69-------------------------------------------------------------------------------
   static const std::string  TmS[] = { "Tm144",  "Tm145",  "Tm146",  "Tm147",  "Tm148",
                                       "Tm149",  "Tm150",  "Tm151",  "Tm152",  "Tm153",
                                       "Tm154",  "Tm155",  "Tm156",  "Tm157",  "Tm158",
                                       "Tm159",  "Tm160",  "Tm161",  "Tm162",  "Tm163",
                                       "Tm164",  "Tm165",  "Tm166",  "Tm167",  "Tm168",
                                       "Tm169",  "Tm170",  "Tm171",  "Tm172",  "Tm173",
                                       "Tm174",  "Tm175",  "Tm176",  "Tm177",  "Tm178",
                                       "Tm179" };
   static const int     TmN[] = {144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
                                 154, 155, 156, 157, 158, 159, 160, 161, 162, 163,
                                 164, 165, 166, 167, 168, 169, 170, 171, 172, 173,
                                 174, 175, 176, 177, 178, 179 };
   static const double  TmA[] = {      143.97628 * u,       144.97039 * u,
                                       145.96684 * u,     146.9613799 * u,
                                      147.958384 * u,       148.95289 * u,
                                       149.95009 * u,      150.945488 * u,
                                      151.944422 * u,       152.94204 * u,
                                       153.94157 * u,       154.93921 * u,
                                      155.938992 * u,      156.936944 * u,
                                       157.93698 * u,      158.934975 * u,
                                      159.935263 * u,      160.933549 * u,
                                      161.934002 * u,     162.9326592 * u,
                                      163.933544 * u,     164.9324431 * u,
                                      165.933561 * u,     166.9328562 * u,
                                     167.9341774 * u,     168.9342179 * u,
                                      169.935806 * u,     170.9364339 * u,
                                     171.9384055 * u,     172.9396084 * u,
                                      173.942173 * u,      174.943841 * u,
                                         175.947 * u,       176.94904 * u,
                                       177.95264 * u,       178.95534 * u };
   static const double  TmW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 1, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double TmIsoMass[36];

   fNISTElementDataTable[68].fElemSymbol    = "Tm";
   fNISTElementDataTable[68].fSymbols       = TmS;
   fNISTElementDataTable[68].fZ             = 69;
   fNISTElementDataTable[68].fNumOfIsotopes = 36;
   fNISTElementDataTable[68].fIsNoStableIsotope       = false;
   fNISTElementDataTable[68].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[68].fNIsos         = TmN;
   fNISTElementDataTable[68].fAIsos         = TmA;
   fNISTElementDataTable[68].fWIsos         = TmW;
   fNISTElementDataTable[68].fMassIsos      = TmIsoMass;


  // Z = 70-------------------------------------------------------------------------------
   static const std::string  YbS[] = { "Yb148",  "Yb149",  "Yb150",  "Yb151",  "Yb152",
                                       "Yb153",  "Yb154",  "Yb155",  "Yb156",  "Yb157",
                                       "Yb158",  "Yb159",  "Yb160",  "Yb161",  "Yb162",
                                       "Yb163",  "Yb164",  "Yb165",  "Yb166",  "Yb167",
                                       "Yb168",  "Yb169",  "Yb170",  "Yb171",  "Yb172",
                                       "Yb173",  "Yb174",  "Yb175",  "Yb176",  "Yb177",
                                       "Yb178",  "Yb179",  "Yb180",  "Yb181" };
   static const int     YbN[] = {148, 149, 150, 151, 152, 153, 154, 155, 156, 157,
                                 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
                                 168, 169, 170, 171, 172, 173, 174, 175, 176, 177,
                                 178, 179, 180, 181 };
   static const double  YbA[] = {      147.96758 * u,       148.96436 * u,
                                       149.95852 * u,        150.9554 * u,
                                       151.95027 * u,       152.94932 * u,
                                      153.946396 * u,      154.945783 * u,
                                      155.942825 * u,      156.942645 * u,
                                     157.9398705 * u,      158.940055 * u,
                                      159.937557 * u,      160.937907 * u,
                                      161.935774 * u,       162.93634 * u,
                                      163.934495 * u,       164.93527 * u,
                                     165.9338747 * u,      166.934953 * u,
                                     167.9338896 * u,     168.9351825 * u,
                                     169.9347664 * u,     170.9363302 * u,
                                     171.9363859 * u,     172.9382151 * u,
                                     173.9388664 * u,     174.9412808 * u,
                                     175.9425764 * u,     176.9452656 * u,
                                      177.946651 * u,       178.95004 * u,
                                       179.95212 * u,       180.95589 * u };
   static const double  YbW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0.00123, 0., 0.02982, 0.1409, 0.2168, 0.16103, 0.32026, 0., 0.12996, 0.,
                                  0., 0., 0., 0. };
   static double YbIsoMass[34];

   fNISTElementDataTable[69].fElemSymbol    = "Yb";
   fNISTElementDataTable[69].fSymbols       = YbS;
   fNISTElementDataTable[69].fZ             = 70;
   fNISTElementDataTable[69].fNumOfIsotopes = 34;
   fNISTElementDataTable[69].fIsNoStableIsotope       = false;
   fNISTElementDataTable[69].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[69].fNIsos         = YbN;
   fNISTElementDataTable[69].fAIsos         = YbA;
   fNISTElementDataTable[69].fWIsos         = YbW;
   fNISTElementDataTable[69].fMassIsos      = YbIsoMass;


  // Z = 71-------------------------------------------------------------------------------
   static const std::string  LuS[] = { "Lu150",  "Lu151",  "Lu152",  "Lu153",  "Lu154",
                                       "Lu155",  "Lu156",  "Lu157",  "Lu158",  "Lu159",
                                       "Lu160",  "Lu161",  "Lu162",  "Lu163",  "Lu164",
                                       "Lu165",  "Lu166",  "Lu167",  "Lu168",  "Lu169",
                                       "Lu170",  "Lu171",  "Lu172",  "Lu173",  "Lu174",
                                       "Lu175",  "Lu176",  "Lu177",  "Lu178",  "Lu179",
                                       "Lu180",  "Lu181",  "Lu182",  "Lu183",  "Lu184",
                                       "Lu185" };
   static const int     LuN[] = {150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
                                 160, 161, 162, 163, 164, 165, 166, 167, 168, 169,
                                 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
                                 180, 181, 182, 183, 184, 185 };
   static const double  LuA[] = {      149.97355 * u,       150.96768 * u,
                                       151.96412 * u,       152.95875 * u,
                                       153.95736 * u,      154.954321 * u,
                                      155.953033 * u,      156.950127 * u,
                                      157.949316 * u,      158.946636 * u,
                                      159.946033 * u,      160.943572 * u,
                                      161.943283 * u,      162.941179 * u,
                                      163.941339 * u,      164.939407 * u,
                                      165.939859 * u,       166.93827 * u,
                                      167.938736 * u,     168.9376441 * u,
                                      169.938478 * u,      170.937917 * u,
                                     171.9390891 * u,      172.938934 * u,
                                     173.9403409 * u,     174.9407752 * u,
                                     175.9426897 * u,     176.9437615 * u,
                                      177.945958 * u,     178.9473309 * u,
                                      179.949888 * u,       180.95191 * u,
                                       181.95504 * u,      182.957363 * u,
                                       183.96091 * u,       184.96362 * u };
   static const double  LuW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0.97401, 0.02599, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double LuIsoMass[36];

   fNISTElementDataTable[70].fElemSymbol    = "Lu";
   fNISTElementDataTable[70].fSymbols       = LuS;
   fNISTElementDataTable[70].fZ             = 71;
   fNISTElementDataTable[70].fNumOfIsotopes = 36;
   fNISTElementDataTable[70].fIsNoStableIsotope       = false;
   fNISTElementDataTable[70].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[70].fNIsos         = LuN;
   fNISTElementDataTable[70].fAIsos         = LuA;
   fNISTElementDataTable[70].fWIsos         = LuW;
   fNISTElementDataTable[70].fMassIsos      = LuIsoMass;


  // Z = 72-------------------------------------------------------------------------------
   static const std::string  HfS[] = { "Hf153",  "Hf154",  "Hf155",  "Hf156",  "Hf157",
                                       "Hf158",  "Hf159",  "Hf160",  "Hf161",  "Hf162",
                                       "Hf163",  "Hf164",  "Hf165",  "Hf166",  "Hf167",
                                       "Hf168",  "Hf169",  "Hf170",  "Hf171",  "Hf172",
                                       "Hf173",  "Hf174",  "Hf175",  "Hf176",  "Hf177",
                                       "Hf178",  "Hf179",  "Hf180",  "Hf181",  "Hf182",
                                       "Hf183",  "Hf184",  "Hf185",  "Hf186",  "Hf187",
                                       "Hf188",  "Hf189" };
   static const int     HfN[] = {153, 154, 155, 156, 157, 158, 159, 160, 161, 162,
                                 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
                                 173, 174, 175, 176, 177, 178, 179, 180, 181, 182,
                                 183, 184, 185, 186, 187, 188, 189 };
   static const double  HfA[] = {      152.97069 * u,       153.96486 * u,
                                       154.96311 * u,       155.95935 * u,
                                       156.95824 * u,      157.954801 * u,
                                      158.953996 * u,      159.950691 * u,
                                      160.950278 * u,     161.9472148 * u,
                                      162.947113 * u,      163.944371 * u,
                                      164.944567 * u,       165.94218 * u,
                                        166.9426 * u,      167.940568 * u,
                                      168.941259 * u,      169.939609 * u,
                                      170.940492 * u,       171.93945 * u,
                                      172.940513 * u,     173.9400461 * u,
                                     174.9415092 * u,     175.9414076 * u,
                                     176.9432277 * u,     177.9437058 * u,
                                     178.9458232 * u,      179.946557 * u,
                                     180.9491083 * u,     181.9505612 * u,
                                       182.95353 * u,      183.955446 * u,
                                      184.958862 * u,      185.960897 * u,
                                       186.96477 * u,       187.96685 * u,
                                       188.97084 * u };
   static const double  HfW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0.0016, 0., 0.0526, 0.186, 0.2728, 0.1362, 0.3508, 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0. };
   static double HfIsoMass[37];

   fNISTElementDataTable[71].fElemSymbol    = "Hf";
   fNISTElementDataTable[71].fSymbols       = HfS;
   fNISTElementDataTable[71].fZ             = 72;
   fNISTElementDataTable[71].fNumOfIsotopes = 37;
   fNISTElementDataTable[71].fIsNoStableIsotope       = false;
   fNISTElementDataTable[71].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[71].fNIsos         = HfN;
   fNISTElementDataTable[71].fAIsos         = HfA;
   fNISTElementDataTable[71].fWIsos         = HfW;
   fNISTElementDataTable[71].fMassIsos      = HfIsoMass;


  // Z = 73-------------------------------------------------------------------------------
   static const std::string  TaS[] = { "Ta155",  "Ta156",  "Ta157",  "Ta158",  "Ta159",
                                       "Ta160",  "Ta161",  "Ta162",  "Ta163",  "Ta164",
                                       "Ta165",  "Ta166",  "Ta167",  "Ta168",  "Ta169",
                                       "Ta170",  "Ta171",  "Ta172",  "Ta173",  "Ta174",
                                       "Ta175",  "Ta176",  "Ta177",  "Ta178",  "Ta179",
                                       "Ta180",  "Ta181",  "Ta182",  "Ta183",  "Ta184",
                                       "Ta185",  "Ta186",  "Ta187",  "Ta188",  "Ta189",
                                       "Ta190",  "Ta191",  "Ta192" };
   static const int     TaN[] = {155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
                                 165, 166, 167, 168, 169, 170, 171, 172, 173, 174,
                                 175, 176, 177, 178, 179, 180, 181, 182, 183, 184,
                                 185, 186, 187, 188, 189, 190, 191, 192 };
   static const double  TaA[] = {      154.97424 * u,       155.97203 * u,
                                       156.96818 * u,       157.96654 * u,
                                      158.963023 * u,      159.961488 * u,
                                      160.958452 * u,      161.957294 * u,
                                      162.954337 * u,      163.953534 * u,
                                      164.950781 * u,      165.950512 * u,
                                      166.948093 * u,      167.948047 * u,
                                      168.946011 * u,      169.946175 * u,
                                      170.944476 * u,      171.944895 * u,
                                       172.94375 * u,      173.944454 * u,
                                      174.943737 * u,      175.944857 * u,
                                     176.9444795 * u,      177.945678 * u,
                                     178.9459366 * u,     179.9474648 * u,
                                     180.9479958 * u,     181.9501519 * u,
                                     182.9513726 * u,      183.954008 * u,
                                      184.955559 * u,      185.958551 * u,
                                      186.960386 * u,      187.963916 * u,
                                       188.96583 * u,       189.96939 * u,
                                       190.97156 * u,       191.97514 * u };
   static const double  TaW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0.0001201, 0.9998799, 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double TaIsoMass[38];

   fNISTElementDataTable[72].fElemSymbol    = "Ta";
   fNISTElementDataTable[72].fSymbols       = TaS;
   fNISTElementDataTable[72].fZ             = 73;
   fNISTElementDataTable[72].fNumOfIsotopes = 38;
   fNISTElementDataTable[72].fIsNoStableIsotope       = false;
   fNISTElementDataTable[72].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[72].fNIsos         = TaN;
   fNISTElementDataTable[72].fAIsos         = TaA;
   fNISTElementDataTable[72].fWIsos         = TaW;
   fNISTElementDataTable[72].fMassIsos      = TaIsoMass;


  // Z = 74-------------------------------------------------------------------------------
   static const std::string   WS[] = {  "W157",   "W158",   "W159",   "W160",   "W161",
                                        "W162",   "W163",   "W164",   "W165",   "W166",
                                        "W167",   "W168",   "W169",   "W170",   "W171",
                                        "W172",   "W173",   "W174",   "W175",   "W176",
                                        "W177",   "W178",   "W179",   "W180",   "W181",
                                        "W182",   "W183",   "W184",   "W185",   "W186",
                                        "W187",   "W188",   "W189",   "W190",   "W191",
                                        "W192",   "W193",   "W194" };
   static const int      WN[] = {157, 158, 159, 160, 161, 162, 163, 164, 165, 166,
                                 167, 168, 169, 170, 171, 172, 173, 174, 175, 176,
                                 177, 178, 179, 180, 181, 182, 183, 184, 185, 186,
                                 187, 188, 189, 190, 191, 192, 193, 194 };
   static const double   WA[] = {      156.97884 * u,       157.97456 * u,
                                       158.97264 * u,       159.96846 * u,
                                        160.9672 * u,      161.963499 * u,
                                      162.962524 * u,      163.958961 * u,
                                      164.958281 * u,      165.955031 * u,
                                      166.954805 * u,      167.951806 * u,
                                      168.951779 * u,      169.949232 * u,
                                      170.949451 * u,      171.947292 * u,
                                      172.947689 * u,      173.946079 * u,
                                      174.946717 * u,      175.945634 * u,
                                      176.946643 * u,      177.945883 * u,
                                      178.947077 * u,     179.9467108 * u,
                                     180.9481978 * u,    181.94820394 * u,
                                    182.95022275 * u,    183.95093092 * u,
                                    184.95341897 * u,     185.9543628 * u,
                                     186.9571588 * u,     187.9584862 * u,
                                      188.961763 * u,      189.963091 * u,
                                      190.966531 * u,       191.96817 * u,
                                       192.97178 * u,       193.97367 * u };
   static const double   WW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0.0012, 0., 0.265, 0.1431, 0.3064, 0., 0.2843,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double WIsoMass[38];

   fNISTElementDataTable[73].fElemSymbol    = "W";
   fNISTElementDataTable[73].fSymbols       = WS;
   fNISTElementDataTable[73].fZ             = 74;
   fNISTElementDataTable[73].fNumOfIsotopes = 38;
   fNISTElementDataTable[73].fIsNoStableIsotope       = false;
   fNISTElementDataTable[73].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[73].fNIsos         = WN;
   fNISTElementDataTable[73].fAIsos         = WA;
   fNISTElementDataTable[73].fWIsos         = WW;
   fNISTElementDataTable[73].fMassIsos      = WIsoMass;


  // Z = 75-------------------------------------------------------------------------------
   static const std::string  ReS[] = { "Re159",  "Re160",  "Re161",  "Re162",  "Re163",
                                       "Re164",  "Re165",  "Re166",  "Re167",  "Re168",
                                       "Re169",  "Re170",  "Re171",  "Re172",  "Re173",
                                       "Re174",  "Re175",  "Re176",  "Re177",  "Re178",
                                       "Re179",  "Re180",  "Re181",  "Re182",  "Re183",
                                       "Re184",  "Re185",  "Re186",  "Re187",  "Re188",
                                       "Re189",  "Re190",  "Re191",  "Re192",  "Re193",
                                       "Re194",  "Re195",  "Re196",  "Re197",  "Re198" };
   static const int     ReN[] = {159, 160, 161, 162, 163, 164, 165, 166, 167, 168,
                                 169, 170, 171, 172, 173, 174, 175, 176, 177, 178,
                                 179, 180, 181, 182, 183, 184, 185, 186, 187, 188,
                                 189, 190, 191, 192, 193, 194, 195, 196, 197, 198 };
   static const double  ReA[] = {      158.98418 * u,       159.98182 * u,
                                       160.97757 * u,       161.97584 * u,
                                       162.97208 * u,      163.970453 * u,
                                      164.967103 * u,      165.965761 * u,
                                      166.962595 * u,      167.961573 * u,
                                      168.958766 * u,       169.95822 * u,
                                      170.955716 * u,       171.95542 * u,
                                      172.953243 * u,      173.953115 * u,
                                      174.951381 * u,      175.951623 * u,
                                      176.950328 * u,      177.950989 * u,
                                      178.949989 * u,      179.950792 * u,
                                      180.950058 * u,       181.95121 * u,
                                     182.9508196 * u,     183.9525228 * u,
                                     184.9529545 * u,     185.9549856 * u,
                                     186.9557501 * u,     187.9581115 * u,
                                      188.959226 * u,      189.961744 * u,
                                      190.963122 * u,      191.966088 * u,
                                      192.967541 * u,       193.97076 * u,
                                       194.97254 * u,        195.9758 * u,
                                       196.97799 * u,        197.9816 * u };
   static const double  ReW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0.374, 0., 0.626, 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double ReIsoMass[40];

   fNISTElementDataTable[74].fElemSymbol    = "Re";
   fNISTElementDataTable[74].fSymbols       = ReS;
   fNISTElementDataTable[74].fZ             = 75;
   fNISTElementDataTable[74].fNumOfIsotopes = 40;
   fNISTElementDataTable[74].fIsNoStableIsotope       = false;
   fNISTElementDataTable[74].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[74].fNIsos         = ReN;
   fNISTElementDataTable[74].fAIsos         = ReA;
   fNISTElementDataTable[74].fWIsos         = ReW;
   fNISTElementDataTable[74].fMassIsos      = ReIsoMass;


  // Z = 76-------------------------------------------------------------------------------
   static const std::string  OsS[] = { "Os161",  "Os162",  "Os163",  "Os164",  "Os165",
                                       "Os166",  "Os167",  "Os168",  "Os169",  "Os170",
                                       "Os171",  "Os172",  "Os173",  "Os174",  "Os175",
                                       "Os176",  "Os177",  "Os178",  "Os179",  "Os180",
                                       "Os181",  "Os182",  "Os183",  "Os184",  "Os185",
                                       "Os186",  "Os187",  "Os188",  "Os189",  "Os190",
                                       "Os191",  "Os192",  "Os193",  "Os194",  "Os195",
                                       "Os196",  "Os197",  "Os198",  "Os199",  "Os200",
                                       "Os201",  "Os202" };
   static const int     OsN[] = {161, 162, 163, 164, 165, 166, 167, 168, 169, 170,
                                 171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
                                 181, 182, 183, 184, 185, 186, 187, 188, 189, 190,
                                 191, 192, 193, 194, 195, 196, 197, 198, 199, 200,
                                 201, 202 };
   static const double  OsA[] = {      160.98903 * u,       161.98443 * u,
                                       162.98241 * u,       163.97802 * u,
                                        164.9766 * u,      165.972692 * u,
                                      166.971549 * u,      167.967808 * u,
                                      168.967018 * u,      169.963578 * u,
                                      170.963174 * u,      171.960017 * u,
                                      172.959808 * u,      173.957064 * u,
                                      174.956945 * u,      175.954806 * u,
                                      176.954966 * u,      177.953254 * u,
                                      178.953817 * u,      179.952375 * u,
                                      180.953247 * u,       181.95211 * u,
                                      182.953125 * u,     183.9524885 * u,
                                     184.9540417 * u,      185.953835 * u,
                                     186.9557474 * u,     187.9558352 * u,
                                     188.9581442 * u,     189.9584437 * u,
                                     190.9609264 * u,      191.961477 * u,
                                     192.9641479 * u,     193.9651772 * u,
                                      194.968318 * u,      195.969641 * u,
                                       196.97283 * u,       197.97441 * u,
                                       198.97801 * u,       199.97984 * u,
                                       200.98364 * u,       201.98595 * u };
   static const double  OsW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0.0002, 0., 0.0159, 0.0196, 0.1324, 0.1615, 0.2626,
                                  0., 0.4078, 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double OsIsoMass[42];

   fNISTElementDataTable[75].fElemSymbol    = "Os";
   fNISTElementDataTable[75].fSymbols       = OsS;
   fNISTElementDataTable[75].fZ             = 76;
   fNISTElementDataTable[75].fNumOfIsotopes = 42;
   fNISTElementDataTable[75].fIsNoStableIsotope       = false;
   fNISTElementDataTable[75].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[75].fNIsos         = OsN;
   fNISTElementDataTable[75].fAIsos         = OsA;
   fNISTElementDataTable[75].fWIsos         = OsW;
   fNISTElementDataTable[75].fMassIsos      = OsIsoMass;


  // Z = 77-------------------------------------------------------------------------------
   static const std::string  IrS[] = { "Ir164",  "Ir165",  "Ir166",  "Ir167",  "Ir168",
                                       "Ir169",  "Ir170",  "Ir171",  "Ir172",  "Ir173",
                                       "Ir174",  "Ir175",  "Ir176",  "Ir177",  "Ir178",
                                       "Ir179",  "Ir180",  "Ir181",  "Ir182",  "Ir183",
                                       "Ir184",  "Ir185",  "Ir186",  "Ir187",  "Ir188",
                                       "Ir189",  "Ir190",  "Ir191",  "Ir192",  "Ir193",
                                       "Ir194",  "Ir195",  "Ir196",  "Ir197",  "Ir198",
                                       "Ir199",  "Ir200",  "Ir201",  "Ir202",  "Ir203",
                                       "Ir204" };
   static const int     IrN[] = {164, 165, 166, 167, 168, 169, 170, 171, 172, 173,
                                 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
                                 184, 185, 186, 187, 188, 189, 190, 191, 192, 193,
                                 194, 195, 196, 197, 198, 199, 200, 201, 202, 203,
                                 204 };
   static const double  IrA[] = {      163.99191 * u,        164.9875 * u,
                                       165.98566 * u,      166.981666 * u,
                                      167.979907 * u,      168.976298 * u,
                                      169.974922 * u,       170.97164 * u,
                                      171.970607 * u,      172.967506 * u,
                                      173.966861 * u,       174.96415 * u,
                                       175.96365 * u,      176.961301 * u,
                                      177.961082 * u,       178.95912 * u,
                                      179.959229 * u,      180.957625 * u,
                                      181.958076 * u,       182.95684 * u,
                                      183.957476 * u,      184.956698 * u,
                                      185.957944 * u,      186.957542 * u,
                                      187.958828 * u,      188.958715 * u,
                                     189.9605412 * u,     190.9605893 * u,
                                     191.9626002 * u,     192.9629216 * u,
                                     193.9650735 * u,     194.9659747 * u,
                                      195.968397 * u,      196.969655 * u,
                                       197.97228 * u,      198.973805 * u,
                                        199.9768 * u,       200.97864 * u,
                                       201.98199 * u,       202.98423 * u,
                                        203.9896 * u };
   static const double  IrW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0.373, 0., 0.627,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0. };
   static double IrIsoMass[41];

   fNISTElementDataTable[76].fElemSymbol    = "Ir";
   fNISTElementDataTable[76].fSymbols       = IrS;
   fNISTElementDataTable[76].fZ             = 77;
   fNISTElementDataTable[76].fNumOfIsotopes = 41;
   fNISTElementDataTable[76].fIsNoStableIsotope       = false;
   fNISTElementDataTable[76].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[76].fNIsos         = IrN;
   fNISTElementDataTable[76].fAIsos         = IrA;
   fNISTElementDataTable[76].fWIsos         = IrW;
   fNISTElementDataTable[76].fMassIsos      = IrIsoMass;


  // Z = 78-------------------------------------------------------------------------------
   static const std::string  PtS[] = { "Pt166",  "Pt167",  "Pt168",  "Pt169",  "Pt170",
                                       "Pt171",  "Pt172",  "Pt173",  "Pt174",  "Pt175",
                                       "Pt176",  "Pt177",  "Pt178",  "Pt179",  "Pt180",
                                       "Pt181",  "Pt182",  "Pt183",  "Pt184",  "Pt185",
                                       "Pt186",  "Pt187",  "Pt188",  "Pt189",  "Pt190",
                                       "Pt191",  "Pt192",  "Pt193",  "Pt194",  "Pt195",
                                       "Pt196",  "Pt197",  "Pt198",  "Pt199",  "Pt200",
                                       "Pt201",  "Pt202",  "Pt203",  "Pt204",  "Pt205",
                                       "Pt206" };
   static const int     PtN[] = {166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
                                 176, 177, 178, 179, 180, 181, 182, 183, 184, 185,
                                 186, 187, 188, 189, 190, 191, 192, 193, 194, 195,
                                 196, 197, 198, 199, 200, 201, 202, 203, 204, 205,
                                 206 };
   static const double  PtA[] = {      165.99486 * u,       166.99269 * u,
                                       167.98813 * u,       168.98657 * u,
                                      169.982496 * u,      170.981245 * u,
                                      171.977351 * u,      172.976443 * u,
                                       173.97282 * u,       174.97241 * u,
                                      175.968938 * u,       176.96847 * u,
                                       177.96565 * u,      178.965359 * u,
                                      179.963032 * u,      180.963098 * u,
                                      181.961172 * u,      182.961597 * u,
                                      183.959915 * u,      184.960614 * u,
                                      185.959351 * u,      186.960617 * u,
                                     187.9593889 * u,      188.960831 * u,
                                     189.9599297 * u,     190.9616729 * u,
                                     191.9610387 * u,     192.9629824 * u,
                                     193.9626809 * u,     194.9647917 * u,
                                    195.96495209 * u,    196.96734069 * u,
                                     197.9678949 * u,     198.9705952 * u,
                                      199.971443 * u,      200.974513 * u,
                                      201.975639 * u,       202.97893 * u,
                                       203.98076 * u,       204.98608 * u,
                                       205.98966 * u };
   static const double  PtW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0.00012, 0., 0.00782, 0., 0.3286, 0.3378,
                                  0.2521, 0., 0.07356, 0., 0., 0., 0., 0., 0., 0.,
                                  0. };
   static double PtIsoMass[41];

   fNISTElementDataTable[77].fElemSymbol    = "Pt";
   fNISTElementDataTable[77].fSymbols       = PtS;
   fNISTElementDataTable[77].fZ             = 78;
   fNISTElementDataTable[77].fNumOfIsotopes = 41;
   fNISTElementDataTable[77].fIsNoStableIsotope       = false;
   fNISTElementDataTable[77].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[77].fNIsos         = PtN;
   fNISTElementDataTable[77].fAIsos         = PtA;
   fNISTElementDataTable[77].fWIsos         = PtW;
   fNISTElementDataTable[77].fMassIsos      = PtIsoMass;


  // Z = 79-------------------------------------------------------------------------------
   static const std::string  AuS[] = { "Au169",  "Au170",  "Au171",  "Au172",  "Au173",
                                       "Au174",  "Au175",  "Au176",  "Au177",  "Au178",
                                       "Au179",  "Au180",  "Au181",  "Au182",  "Au183",
                                       "Au184",  "Au185",  "Au186",  "Au187",  "Au188",
                                       "Au189",  "Au190",  "Au191",  "Au192",  "Au193",
                                       "Au194",  "Au195",  "Au196",  "Au197",  "Au198",
                                       "Au199",  "Au200",  "Au201",  "Au202",  "Au203",
                                       "Au204",  "Au205",  "Au206",  "Au207",  "Au208",
                                       "Au209",  "Au210" };
   static const int     AuN[] = {169, 170, 171, 172, 173, 174, 175, 176, 177, 178,
                                 179, 180, 181, 182, 183, 184, 185, 186, 187, 188,
                                 189, 190, 191, 192, 193, 194, 195, 196, 197, 198,
                                 199, 200, 201, 202, 203, 204, 205, 206, 207, 208,
                                 209, 210 };
   static const double  AuA[] = {      168.99808 * u,       169.99597 * u,
                                      170.991876 * u,      171.989942 * u,
                                      172.986241 * u,      173.984717 * u,
                                      174.981304 * u,       175.98025 * u,
                                       176.97687 * u,      177.976032 * u,
                                      178.973174 * u,      179.972523 * u,
                                      180.970079 * u,      181.969618 * u,
                                      182.967591 * u,      183.967452 * u,
                                       184.96579 * u,      185.965953 * u,
                                      186.964543 * u,      187.965349 * u,
                                      188.963948 * u,      189.964698 * u,
                                      190.963702 * u,      191.964814 * u,
                                     192.9641373 * u,     193.9654178 * u,
                                     194.9650352 * u,     195.9665699 * u,
                                    196.96656879 * u,    197.96824242 * u,
                                    198.96876528 * u,      199.970756 * u,
                                     200.9716575 * u,      201.973856 * u,
                                     202.9751544 * u,       203.97783 * u,
                                       204.97985 * u,       205.98474 * u,
                                        206.9884 * u,       207.99345 * u,
                                       208.99735 * u,        210.0025 * u };
   static const double  AuW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 1, 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double AuIsoMass[42];

   fNISTElementDataTable[78].fElemSymbol    = "Au";
   fNISTElementDataTable[78].fSymbols       = AuS;
   fNISTElementDataTable[78].fZ             = 79;
   fNISTElementDataTable[78].fNumOfIsotopes = 42;
   fNISTElementDataTable[78].fIsNoStableIsotope       = false;
   fNISTElementDataTable[78].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[78].fNIsos         = AuN;
   fNISTElementDataTable[78].fAIsos         = AuA;
   fNISTElementDataTable[78].fWIsos         = AuW;
   fNISTElementDataTable[78].fMassIsos      = AuIsoMass;


  // Z = 80-------------------------------------------------------------------------------
   static const std::string  HgS[] = { "Hg171",  "Hg172",  "Hg173",  "Hg174",  "Hg175",
                                       "Hg176",  "Hg177",  "Hg178",  "Hg179",  "Hg180",
                                       "Hg181",  "Hg182",  "Hg183",  "Hg184",  "Hg185",
                                       "Hg186",  "Hg187",  "Hg188",  "Hg189",  "Hg190",
                                       "Hg191",  "Hg192",  "Hg193",  "Hg194",  "Hg195",
                                       "Hg196",  "Hg197",  "Hg198",  "Hg199",  "Hg200",
                                       "Hg201",  "Hg202",  "Hg203",  "Hg204",  "Hg205",
                                       "Hg206",  "Hg207",  "Hg208",  "Hg209",  "Hg210",
                                       "Hg211",  "Hg212",  "Hg213",  "Hg214",  "Hg215",
                                       "Hg216" };
   static const int     HgN[] = {171, 172, 173, 174, 175, 176, 177, 178, 179, 180,
                                 181, 182, 183, 184, 185, 186, 187, 188, 189, 190,
                                 191, 192, 193, 194, 195, 196, 197, 198, 199, 200,
                                 201, 202, 203, 204, 205, 206, 207, 208, 209, 210,
                                 211, 212, 213, 214, 215, 216 };
   static const double  HgA[] = {      171.00353 * u,       171.99881 * u,
                                       172.99709 * u,      173.992865 * u,
                                      174.991441 * u,      175.987361 * u,
                                      176.986277 * u,      177.982484 * u,
                                      178.981831 * u,       179.97826 * u,
                                      180.977819 * u,      181.974689 * u,
                                     182.9744448 * u,      183.971714 * u,
                                      184.971899 * u,      185.969362 * u,
                                      186.969814 * u,      187.967567 * u,
                                      188.968195 * u,      189.966323 * u,
                                      190.967157 * u,      191.965635 * u,
                                      192.966653 * u,     193.9654491 * u,
                                      194.966721 * u,     195.9658326 * u,
                                     196.9672128 * u,     197.9667686 * u,
                                    198.96828064 * u,    199.96832659 * u,
                                    200.97030284 * u,     201.9706434 * u,
                                     202.9728728 * u,    203.97349398 * u,
                                     204.9760734 * u,      205.977514 * u,
                                        206.9823 * u,      207.985759 * u,
                                       208.99072 * u,       209.99424 * u,
                                       210.99933 * u,       212.00296 * u,
                                       213.00823 * u,         214.012 * u,
                                        215.0174 * u,       216.02132 * u };
   static const double  HgW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0.0015, 0., 0.0997, 0.1687, 0.231,
                                  0.1318, 0.2986, 0., 0.0687, 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double HgIsoMass[46];

   fNISTElementDataTable[79].fElemSymbol    = "Hg";
   fNISTElementDataTable[79].fSymbols       = HgS;
   fNISTElementDataTable[79].fZ             = 80;
   fNISTElementDataTable[79].fNumOfIsotopes = 46;
   fNISTElementDataTable[79].fIsNoStableIsotope       = false;
   fNISTElementDataTable[79].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[79].fNIsos         = HgN;
   fNISTElementDataTable[79].fAIsos         = HgA;
   fNISTElementDataTable[79].fWIsos         = HgW;
   fNISTElementDataTable[79].fMassIsos      = HgIsoMass;


  // Z = 81-------------------------------------------------------------------------------
   static const std::string  TlS[] = { "Tl176",  "Tl177",  "Tl178",  "Tl179",  "Tl180",
                                       "Tl181",  "Tl182",  "Tl183",  "Tl184",  "Tl185",
                                       "Tl186",  "Tl187",  "Tl188",  "Tl189",  "Tl190",
                                       "Tl191",  "Tl192",  "Tl193",  "Tl194",  "Tl195",
                                       "Tl196",  "Tl197",  "Tl198",  "Tl199",  "Tl200",
                                       "Tl201",  "Tl202",  "Tl203",  "Tl204",  "Tl205",
                                       "Tl206",  "Tl207",  "Tl208",  "Tl209",  "Tl210",
                                       "Tl211",  "Tl212",  "Tl213",  "Tl214",  "Tl215",
                                       "Tl216",  "Tl217",  "Tl218" };
   static const int     TlN[] = {176, 177, 178, 179, 180, 181, 182, 183, 184, 185,
                                 186, 187, 188, 189, 190, 191, 192, 193, 194, 195,
                                 196, 197, 198, 199, 200, 201, 202, 203, 204, 205,
                                 206, 207, 208, 209, 210, 211, 212, 213, 214, 215,
                                 216, 217, 218 };
   static const double  TlA[] = {     176.000624 * u,      176.996431 * u,
                                       177.99485 * u,      178.991111 * u,
                                      179.990057 * u,       180.98626 * u,
                                      181.985713 * u,      182.982193 * u,
                                      183.981886 * u,      184.978789 * u,
                                      185.978651 * u,     186.9759063 * u,
                                      187.976021 * u,      188.973588 * u,
                                      189.973828 * u,     190.9717842 * u,
                                      191.972225 * u,      192.970502 * u,
                                      193.971081 * u,      194.969774 * u,
                                      195.970481 * u,      196.969576 * u,
                                      197.970483 * u,      198.969877 * u,
                                     199.9709633 * u,      200.970822 * u,
                                      201.972102 * u,     202.9723446 * u,
                                     203.9738639 * u,     204.9744278 * u,
                                     205.9761106 * u,     206.9774197 * u,
                                      207.982019 * u,     208.9853594 * u,
                                      209.990074 * u,      210.993475 * u,
                                       211.99834 * u,      213.001915 * u,
                                       214.00694 * u,       215.01064 * u,
                                        216.0158 * u,       217.01966 * u,
                                       218.02479 * u };
   static const double  TlW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0.2952, 0., 0.7048,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0. };
   static double TlIsoMass[43];

   fNISTElementDataTable[80].fElemSymbol    = "Tl";
   fNISTElementDataTable[80].fSymbols       = TlS;
   fNISTElementDataTable[80].fZ             = 81;
   fNISTElementDataTable[80].fNumOfIsotopes = 43;
   fNISTElementDataTable[80].fIsNoStableIsotope       = false;
   fNISTElementDataTable[80].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[80].fNIsos         = TlN;
   fNISTElementDataTable[80].fAIsos         = TlA;
   fNISTElementDataTable[80].fWIsos         = TlW;
   fNISTElementDataTable[80].fMassIsos      = TlIsoMass;


  // Z = 82-------------------------------------------------------------------------------
   static const std::string  PbS[] = { "Pb178",  "Pb179",  "Pb180",  "Pb181",  "Pb182",
                                       "Pb183",  "Pb184",  "Pb185",  "Pb186",  "Pb187",
                                       "Pb188",  "Pb189",  "Pb190",  "Pb191",  "Pb192",
                                       "Pb193",  "Pb194",  "Pb195",  "Pb196",  "Pb197",
                                       "Pb198",  "Pb199",  "Pb200",  "Pb201",  "Pb202",
                                       "Pb203",  "Pb204",  "Pb205",  "Pb206",  "Pb207",
                                       "Pb208",  "Pb209",  "Pb210",  "Pb211",  "Pb212",
                                       "Pb213",  "Pb214",  "Pb215",  "Pb216",  "Pb217",
                                       "Pb218",  "Pb219",  "Pb220" };
   static const int     PbN[] = {178, 179, 180, 181, 182, 183, 184, 185, 186, 187,
                                 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
                                 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
                                 208, 209, 210, 211, 212, 213, 214, 215, 216, 217,
                                 218, 219, 220 };
   static const double  PbA[] = {     178.003831 * u,      179.002201 * u,
                                      179.997928 * u,      180.996653 * u,
                                      181.992672 * u,      182.991872 * u,
                                      183.988136 * u,       184.98761 * u,
                                      185.984238 * u,     186.9839109 * u,
                                      187.980875 * u,      188.980807 * u,
                                      189.978082 * u,      190.978276 * u,
                                      191.975775 * u,      192.976173 * u,
                                      193.974012 * u,      194.974543 * u,
                                      195.972774 * u,     196.9734312 * u,
                                      197.972034 * u,      198.972913 * u,
                                      199.971819 * u,      200.972883 * u,
                                      201.972152 * u,     202.9733911 * u,
                                      203.973044 * u,     204.9744822 * u,
                                     205.9744657 * u,     206.9758973 * u,
                                     207.9766525 * u,     208.9810905 * u,
                                     209.9841889 * u,     210.9887371 * u,
                                     211.9918977 * u,     212.9965629 * u,
                                     213.9998059 * u,       215.00474 * u,
                                       216.00803 * u,       217.01314 * u,
                                       218.01659 * u,       219.02177 * u,
                                       220.02541 * u };
   static const double  PbW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0.014, 0., 0.241, 0.221,
                                  0.524, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0. };
   static double PbIsoMass[43];

   fNISTElementDataTable[81].fElemSymbol    = "Pb";
   fNISTElementDataTable[81].fSymbols       = PbS;
   fNISTElementDataTable[81].fZ             = 82;
   fNISTElementDataTable[81].fNumOfIsotopes = 43;
   fNISTElementDataTable[81].fIsNoStableIsotope       = false;
   fNISTElementDataTable[81].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[81].fNIsos         = PbN;
   fNISTElementDataTable[81].fAIsos         = PbA;
   fNISTElementDataTable[81].fWIsos         = PbW;
   fNISTElementDataTable[81].fMassIsos      = PbIsoMass;


  // Z = 83-------------------------------------------------------------------------------
   static const std::string  BiS[] = { "Bi184",  "Bi185",  "Bi186",  "Bi187",  "Bi188",
                                       "Bi189",  "Bi190",  "Bi191",  "Bi192",  "Bi193",
                                       "Bi194",  "Bi195",  "Bi196",  "Bi197",  "Bi198",
                                       "Bi199",  "Bi200",  "Bi201",  "Bi202",  "Bi203",
                                       "Bi204",  "Bi205",  "Bi206",  "Bi207",  "Bi208",
                                       "Bi209",  "Bi210",  "Bi211",  "Bi212",  "Bi213",
                                       "Bi214",  "Bi215",  "Bi216",  "Bi217",  "Bi218",
                                       "Bi219",  "Bi220",  "Bi221",  "Bi222",  "Bi223",
                                       "Bi224" };
   static const int     BiN[] = {184, 185, 186, 187, 188, 189, 190, 191, 192, 193,
                                 194, 195, 196, 197, 198, 199, 200, 201, 202, 203,
                                 204, 205, 206, 207, 208, 209, 210, 211, 212, 213,
                                 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
                                 224 };
   static const double  BiA[] = {     184.001275 * u,        184.9976 * u,
                                      185.996644 * u,      186.993147 * u,
                                      187.992287 * u,      188.989195 * u,
                                      189.988622 * u,     190.9857866 * u,
                                      191.985469 * u,       192.98296 * u,
                                      193.982785 * u,     194.9806488 * u,
                                      195.980667 * u,     196.9788651 * u,
                                      197.979206 * u,      198.977673 * u,
                                      199.978131 * u,       200.97701 * u,
                                      201.977734 * u,      202.976893 * u,
                                     203.9778361 * u,     204.9773867 * u,
                                     205.9784993 * u,      206.978471 * u,
                                     207.9797425 * u,     208.9803991 * u,
                                     209.9841207 * u,     210.9872697 * u,
                                      211.991286 * u,     212.9943851 * u,
                                      213.998712 * u,       215.00177 * u,
                                      216.006306 * u,      217.009372 * u,
                                      218.014188 * u,       219.01748 * u,
                                       220.02235 * u,       221.02587 * u,
                                       222.03078 * u,        223.0345 * u,
                                       224.03947 * u };
   static const double  BiW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 1, 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0. };
   static double BiIsoMass[41];

   fNISTElementDataTable[82].fElemSymbol    = "Bi";
   fNISTElementDataTable[82].fSymbols       = BiS;
   fNISTElementDataTable[82].fZ             = 83;
   fNISTElementDataTable[82].fNumOfIsotopes = 41;
   fNISTElementDataTable[82].fIsNoStableIsotope       = false;
   fNISTElementDataTable[82].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[82].fNIsos         = BiN;
   fNISTElementDataTable[82].fAIsos         = BiA;
   fNISTElementDataTable[82].fWIsos         = BiW;
   fNISTElementDataTable[82].fMassIsos      = BiIsoMass;


  // Z = 84-------------------------------------------------------------------------------
   static const std::string  PoS[] = { "Po186",  "Po187",  "Po188",  "Po189",  "Po190",
                                       "Po191",  "Po192",  "Po193",  "Po194",  "Po195",
                                       "Po196",  "Po197",  "Po198",  "Po199",  "Po200",
                                       "Po201",  "Po202",  "Po203",  "Po204",  "Po205",
                                       "Po206",  "Po207",  "Po208",  "Po209",  "Po210",
                                       "Po211",  "Po212",  "Po213",  "Po214",  "Po215",
                                       "Po216",  "Po217",  "Po218",  "Po219",  "Po220",
                                       "Po221",  "Po222",  "Po223",  "Po224",  "Po225",
                                       "Po226",  "Po227" };
   static const int     PoN[] = {186, 187, 188, 189, 190, 191, 192, 193, 194, 195,
                                 196, 197, 198, 199, 200, 201, 202, 203, 204, 205,
                                 206, 207, 208, 209, 210, 211, 212, 213, 214, 215,
                                 216, 217, 218, 219, 220, 221, 222, 223, 224, 225,
                                 226, 227 };
   static const double  PoA[] = {     186.004393 * u,      187.003041 * u,
                                      187.999416 * u,      188.998473 * u,
                                      189.995101 * u,     190.9945585 * u,
                                      191.991336 * u,      192.991026 * u,
                                      193.988186 * u,      194.988126 * u,
                                      195.985526 * u,       196.98566 * u,
                                      197.983389 * u,      198.983667 * u,
                                      199.981799 * u,     200.9822598 * u,
                                      201.980758 * u,     202.9814161 * u,
                                       203.98031 * u,      204.981203 * u,
                                      205.980474 * u,     206.9815938 * u,
                                     207.9812461 * u,     208.9824308 * u,
                                     209.9828741 * u,     210.9866536 * u,
                                     211.9888684 * u,     212.9928576 * u,
                                     213.9952017 * u,     214.9994201 * u,
                                     216.0019152 * u,     217.0063182 * u,
                                     218.0089735 * u,      219.013614 * u,
                                      220.016386 * u,      221.021228 * u,
                                       222.02414 * u,       223.02907 * u,
                                       224.03211 * u,       225.03707 * u,
                                       226.04031 * u,       227.04539 * u };
   static const double  PoW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double PoIsoMass[42];

   fNISTElementDataTable[83].fElemSymbol    = "Po";
   fNISTElementDataTable[83].fSymbols       = PoS;
   fNISTElementDataTable[83].fZ             = 84;
   fNISTElementDataTable[83].fNumOfIsotopes = 42;
   fNISTElementDataTable[83].fIsNoStableIsotope       = true;
   fNISTElementDataTable[83].fIndxOfMostStableIsotope = 23;
   fNISTElementDataTable[83].fNIsos         = PoN;
   fNISTElementDataTable[83].fAIsos         = PoA;
   fNISTElementDataTable[83].fWIsos         = PoW;
   fNISTElementDataTable[83].fMassIsos      = PoIsoMass;


  // Z = 85-------------------------------------------------------------------------------
   static const std::string  AtS[] = { "At191",  "At192",  "At193",  "At194",  "At195",
                                       "At196",  "At197",  "At198",  "At199",  "At200",
                                       "At201",  "At202",  "At203",  "At204",  "At205",
                                       "At206",  "At207",  "At208",  "At209",  "At210",
                                       "At211",  "At212",  "At213",  "At214",  "At215",
                                       "At216",  "At217",  "At218",  "At219",  "At220",
                                       "At221",  "At222",  "At223",  "At224",  "At225",
                                       "At226",  "At227",  "At228",  "At229" };
   static const int     AtN[] = {191, 192, 193, 194, 195, 196, 197, 198, 199, 200,
                                 201, 202, 203, 204, 205, 206, 207, 208, 209, 210,
                                 211, 212, 213, 214, 215, 216, 217, 218, 219, 220,
                                 221, 222, 223, 224, 225, 226, 227, 228, 229 };
   static const double  AtA[] = {     191.004148 * u,      192.003152 * u,
                                      192.999927 * u,      193.999236 * u,
                                     194.9962685 * u,        195.9958 * u,
                                      196.993189 * u,      197.992784 * u,
                                     198.9905277 * u,      199.990351 * u,
                                     200.9884171 * u,       201.98863 * u,
                                      202.986943 * u,      203.987251 * u,
                                      204.986076 * u,      205.986657 * u,
                                        206.9858 * u,     207.9866133 * u,
                                     208.9861702 * u,     209.9871479 * u,
                                     210.9874966 * u,     211.9907377 * u,
                                      212.992937 * u,     213.9963721 * u,
                                     214.9986528 * u,     216.0024236 * u,
                                     217.0047192 * u,      218.008695 * u,
                                     219.0111618 * u,      220.015433 * u,
                                      221.018017 * u,      222.022494 * u,
                                      223.025151 * u,      224.029749 * u,
                                       225.03263 * u,       226.03716 * u,
                                       227.04024 * u,       228.04475 * u,
                                       229.04812 * u };
   static const double  AtW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double AtIsoMass[39];

   fNISTElementDataTable[84].fElemSymbol    = "At";
   fNISTElementDataTable[84].fSymbols       = AtS;
   fNISTElementDataTable[84].fZ             = 85;
   fNISTElementDataTable[84].fNumOfIsotopes = 39;
   fNISTElementDataTable[84].fIsNoStableIsotope       = true;
   fNISTElementDataTable[84].fIndxOfMostStableIsotope = 19;
   fNISTElementDataTable[84].fNIsos         = AtN;
   fNISTElementDataTable[84].fAIsos         = AtA;
   fNISTElementDataTable[84].fWIsos         = AtW;
   fNISTElementDataTable[84].fMassIsos      = AtIsoMass;


  // Z = 86-------------------------------------------------------------------------------
   static const std::string  RnS[] = { "Rn193",  "Rn194",  "Rn195",  "Rn196",  "Rn197",
                                       "Rn198",  "Rn199",  "Rn200",  "Rn201",  "Rn202",
                                       "Rn203",  "Rn204",  "Rn205",  "Rn206",  "Rn207",
                                       "Rn208",  "Rn209",  "Rn210",  "Rn211",  "Rn212",
                                       "Rn213",  "Rn214",  "Rn215",  "Rn216",  "Rn217",
                                       "Rn218",  "Rn219",  "Rn220",  "Rn221",  "Rn222",
                                       "Rn223",  "Rn224",  "Rn225",  "Rn226",  "Rn227",
                                       "Rn228",  "Rn229",  "Rn230",  "Rn231" };
   static const int     RnN[] = {193, 194, 195, 196, 197, 198, 199, 200, 201, 202,
                                 203, 204, 205, 206, 207, 208, 209, 210, 211, 212,
                                 213, 214, 215, 216, 217, 218, 219, 220, 221, 222,
                                 223, 224, 225, 226, 227, 228, 229, 230, 231 };
   static const double  RnA[] = {     193.009708 * u,      194.006144 * u,
                                      195.005422 * u,      196.002116 * u,
                                      197.001585 * u,      197.998679 * u,
                                       198.99839 * u,       199.99569 * u,
                                      200.995628 * u,      201.993264 * u,
                                      202.993388 * u,       203.99143 * u,
                                      204.991719 * u,      205.990214 * u,
                                     206.9907303 * u,      207.989635 * u,
                                      208.990415 * u,     209.9896891 * u,
                                     210.9906011 * u,     211.9907039 * u,
                                     212.9938831 * u,      213.995363 * u,
                                     214.9987459 * u,     216.0002719 * u,
                                      217.003928 * u,     218.0056016 * u,
                                     219.0094804 * u,     220.0113941 * u,
                                     221.0155371 * u,     222.0175782 * u,
                                     223.0218893 * u,      224.024096 * u,
                                      225.028486 * u,      226.030861 * u,
                                      227.035304 * u,      228.037835 * u,
                                      229.042257 * u,       230.04514 * u,
                                       231.04987 * u };
   static const double  RnW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double RnIsoMass[39];

   fNISTElementDataTable[85].fElemSymbol    = "Rn";
   fNISTElementDataTable[85].fSymbols       = RnS;
   fNISTElementDataTable[85].fZ             = 86;
   fNISTElementDataTable[85].fNumOfIsotopes = 39;
   fNISTElementDataTable[85].fIsNoStableIsotope       = true;
   fNISTElementDataTable[85].fIndxOfMostStableIsotope = 29;
   fNISTElementDataTable[85].fNIsos         = RnN;
   fNISTElementDataTable[85].fAIsos         = RnA;
   fNISTElementDataTable[85].fWIsos         = RnW;
   fNISTElementDataTable[85].fMassIsos      = RnIsoMass;


  // Z = 87-------------------------------------------------------------------------------
   static const std::string  FrS[] = { "Fr199",  "Fr200",  "Fr201",  "Fr202",  "Fr203",
                                       "Fr204",  "Fr205",  "Fr206",  "Fr207",  "Fr208",
                                       "Fr209",  "Fr210",  "Fr211",  "Fr212",  "Fr213",
                                       "Fr214",  "Fr215",  "Fr216",  "Fr217",  "Fr218",
                                       "Fr219",  "Fr220",  "Fr221",  "Fr222",  "Fr223",
                                       "Fr224",  "Fr225",  "Fr226",  "Fr227",  "Fr228",
                                       "Fr229",  "Fr230",  "Fr231",  "Fr232",  "Fr233" };
   static const int     FrN[] = {199, 200, 201, 202, 203, 204, 205, 206, 207, 208,
                                 209, 210, 211, 212, 213, 214, 215, 216, 217, 218,
                                 219, 220, 221, 222, 223, 224, 225, 226, 227, 228,
                                 229, 230, 231, 232, 233 };
   static const double  FrA[] = {     199.007259 * u,      200.006586 * u,
                                      201.003867 * u,       202.00332 * u,
                                     203.0009407 * u,      204.000652 * u,
                                     204.9985939 * u,      205.998666 * u,
                                      206.996946 * u,      207.997138 * u,
                                      208.995955 * u,      209.996422 * u,
                                      210.995556 * u,     211.9962257 * u,
                                      212.996186 * u,     213.9989713 * u,
                                     215.0003418 * u,     216.0031899 * u,
                                     217.0046323 * u,     218.0075787 * u,
                                     219.0092524 * u,     220.0123277 * u,
                                     221.0142552 * u,      222.017552 * u,
                                      223.019736 * u,      224.023398 * u,
                                      225.025573 * u,      226.029566 * u,
                                      227.031869 * u,      228.035823 * u,
                                      229.038298 * u,      230.042416 * u,
                                      231.045158 * u,       232.04937 * u,
                                       233.05264 * u };
   static const double  FrW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double FrIsoMass[35];

   fNISTElementDataTable[86].fElemSymbol    = "Fr";
   fNISTElementDataTable[86].fSymbols       = FrS;
   fNISTElementDataTable[86].fZ             = 87;
   fNISTElementDataTable[86].fNumOfIsotopes = 35;
   fNISTElementDataTable[86].fIsNoStableIsotope       = true;
   fNISTElementDataTable[86].fIndxOfMostStableIsotope = 24;
   fNISTElementDataTable[86].fNIsos         = FrN;
   fNISTElementDataTable[86].fAIsos         = FrA;
   fNISTElementDataTable[86].fWIsos         = FrW;
   fNISTElementDataTable[86].fMassIsos      = FrIsoMass;


  // Z = 88-------------------------------------------------------------------------------
   static const std::string  RaS[] = { "Ra201",  "Ra202",  "Ra203",  "Ra204",  "Ra205",
                                       "Ra206",  "Ra207",  "Ra208",  "Ra209",  "Ra210",
                                       "Ra211",  "Ra212",  "Ra213",  "Ra214",  "Ra215",
                                       "Ra216",  "Ra217",  "Ra218",  "Ra219",  "Ra220",
                                       "Ra221",  "Ra222",  "Ra223",  "Ra224",  "Ra225",
                                       "Ra226",  "Ra227",  "Ra228",  "Ra229",  "Ra230",
                                       "Ra231",  "Ra232",  "Ra233",  "Ra234",  "Ra235" };
   static const int     RaN[] = {201, 202, 203, 204, 205, 206, 207, 208, 209, 210,
                                 211, 212, 213, 214, 215, 216, 217, 218, 219, 220,
                                 221, 222, 223, 224, 225, 226, 227, 228, 229, 230,
                                 231, 232, 233, 234, 235 };
   static const double  RaA[] = {      201.01271 * u,       202.00976 * u,
                                      203.009304 * u,      204.006492 * u,
                                      205.006268 * u,      206.003828 * u,
                                      207.003799 * u,      208.001841 * u,
                                       209.00199 * u,      210.000494 * u,
                                     211.0008932 * u,      211.999787 * u,
                                      213.000384 * u,     214.0000997 * u,
                                     215.0027204 * u,     216.0035334 * u,
                                     217.0063207 * u,      218.007141 * u,
                                     219.0100855 * u,     220.0110259 * u,
                                     221.0139177 * u,     222.0153748 * u,
                                     223.0185023 * u,      224.020212 * u,
                                     225.0236119 * u,     226.0254103 * u,
                                     227.0291783 * u,     228.0310707 * u,
                                      229.034942 * u,      230.037055 * u,
                                      231.041027 * u,     232.0434753 * u,
                                      233.047582 * u,      234.050342 * u,
                                       235.05497 * u };
   static const double  RaW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double RaIsoMass[35];

   fNISTElementDataTable[87].fElemSymbol    = "Ra";
   fNISTElementDataTable[87].fSymbols       = RaS;
   fNISTElementDataTable[87].fZ             = 88;
   fNISTElementDataTable[87].fNumOfIsotopes = 35;
   fNISTElementDataTable[87].fIsNoStableIsotope       = true;
   fNISTElementDataTable[87].fIndxOfMostStableIsotope = 25;
   fNISTElementDataTable[87].fNIsos         = RaN;
   fNISTElementDataTable[87].fAIsos         = RaA;
   fNISTElementDataTable[87].fWIsos         = RaW;
   fNISTElementDataTable[87].fMassIsos      = RaIsoMass;


  // Z = 89-------------------------------------------------------------------------------
   static const std::string  AcS[] = { "Ac206",  "Ac207",  "Ac208",  "Ac209",  "Ac210",
                                       "Ac211",  "Ac212",  "Ac213",  "Ac214",  "Ac215",
                                       "Ac216",  "Ac217",  "Ac218",  "Ac219",  "Ac220",
                                       "Ac221",  "Ac222",  "Ac223",  "Ac224",  "Ac225",
                                       "Ac226",  "Ac227",  "Ac228",  "Ac229",  "Ac230",
                                       "Ac231",  "Ac232",  "Ac233",  "Ac234",  "Ac235",
                                       "Ac236",  "Ac237" };
   static const int     AcN[] = {206, 207, 208, 209, 210, 211, 212, 213, 214, 215,
                                 216, 217, 218, 219, 220, 221, 222, 223, 224, 225,
                                 226, 227, 228, 229, 230, 231, 232, 233, 234, 235,
                                 236, 237 };
   static const double  AcA[] = {     206.014452 * u,      207.011966 * u,
                                       208.01155 * u,      209.009495 * u,
                                      210.009436 * u,      211.007732 * u,
                                      212.007813 * u,      213.006609 * u,
                                      214.006918 * u,      215.006475 * u,
                                      216.008743 * u,      217.009344 * u,
                                      218.011642 * u,      219.012421 * u,
                                     220.0147549 * u,      221.015592 * u,
                                     222.0178442 * u,     223.0191377 * u,
                                     224.0217232 * u,       225.02323 * u,
                                     226.0260984 * u,     227.0277523 * u,
                                     228.0310215 * u,      229.032956 * u,
                                      230.036327 * u,      231.038393 * u,
                                      232.042034 * u,      233.044346 * u,
                                      234.048139 * u,       235.05084 * u,
                                      236.054988 * u,       237.05827 * u };
   static const double  AcW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double AcIsoMass[32];

   fNISTElementDataTable[88].fElemSymbol    = "Ac";
   fNISTElementDataTable[88].fSymbols       = AcS;
   fNISTElementDataTable[88].fZ             = 89;
   fNISTElementDataTable[88].fNumOfIsotopes = 32;
   fNISTElementDataTable[88].fIsNoStableIsotope       = true;
   fNISTElementDataTable[88].fIndxOfMostStableIsotope = 21;
   fNISTElementDataTable[88].fNIsos         = AcN;
   fNISTElementDataTable[88].fAIsos         = AcA;
   fNISTElementDataTable[88].fWIsos         = AcW;
   fNISTElementDataTable[88].fMassIsos      = AcIsoMass;


  // Z = 90-------------------------------------------------------------------------------
   static const std::string  ThS[] = { "Th208",  "Th209",  "Th210",  "Th211",  "Th212",
                                       "Th213",  "Th214",  "Th215",  "Th216",  "Th217",
                                       "Th218",  "Th219",  "Th220",  "Th221",  "Th222",
                                       "Th223",  "Th224",  "Th225",  "Th226",  "Th227",
                                       "Th228",  "Th229",  "Th230",  "Th231",  "Th232",
                                       "Th233",  "Th234",  "Th235",  "Th236",  "Th237",
                                       "Th238",  "Th239" };
   static const int     ThN[] = {208, 209, 210, 211, 212, 213, 214, 215, 216, 217,
                                 218, 219, 220, 221, 222, 223, 224, 225, 226, 227,
                                 228, 229, 230, 231, 232, 233, 234, 235, 236, 237,
                                 238, 239 };
   static const double  ThA[] = {       208.0179 * u,      209.017753 * u,
                                      210.015094 * u,      211.014929 * u,
                                      212.012988 * u,      213.013009 * u,
                                        214.0115 * u,     215.0117248 * u,
                                      216.011056 * u,      217.013117 * u,
                                      218.013276 * u,      219.015537 * u,
                                      220.015748 * u,      221.018184 * u,
                                      222.018469 * u,     223.0208119 * u,
                                      224.021464 * u,     225.0239514 * u,
                                     226.0249034 * u,     227.0277042 * u,
                                     228.0287413 * u,     229.0317627 * u,
                                     230.0331341 * u,     231.0363046 * u,
                                     232.0380558 * u,     233.0415823 * u,
                                     234.0436014 * u,      235.047255 * u,
                                      236.049657 * u,      237.053629 * u,
                                        238.0565 * u,       239.06077 * u };
   static const double  ThW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 1, 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double ThIsoMass[32];

   fNISTElementDataTable[89].fElemSymbol    = "Th";
   fNISTElementDataTable[89].fSymbols       = ThS;
   fNISTElementDataTable[89].fZ             = 90;
   fNISTElementDataTable[89].fNumOfIsotopes = 32;
   fNISTElementDataTable[89].fIsNoStableIsotope       = false;
   fNISTElementDataTable[89].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[89].fNIsos         = ThN;
   fNISTElementDataTable[89].fAIsos         = ThA;
   fNISTElementDataTable[89].fWIsos         = ThW;
   fNISTElementDataTable[89].fMassIsos      = ThIsoMass;


  // Z = 91-------------------------------------------------------------------------------
   static const std::string  PaS[] = { "Pa212",  "Pa213",  "Pa214",  "Pa215",  "Pa216",
                                       "Pa217",  "Pa218",  "Pa219",  "Pa220",  "Pa221",
                                       "Pa222",  "Pa223",  "Pa224",  "Pa225",  "Pa226",
                                       "Pa227",  "Pa228",  "Pa229",  "Pa230",  "Pa231",
                                       "Pa232",  "Pa233",  "Pa234",  "Pa235",  "Pa236",
                                       "Pa237",  "Pa238",  "Pa239",  "Pa240",  "Pa241" };
   static const int     PaN[] = {212, 213, 214, 215, 216, 217, 218, 219, 220, 221,
                                 222, 223, 224, 225, 226, 227, 228, 229, 230, 231,
                                 232, 233, 234, 235, 236, 237, 238, 239, 240, 241 };
   static const double  PaA[] = {     212.023203 * u,      213.021109 * u,
                                      214.020918 * u,      215.019183 * u,
                                      216.019109 * u,      217.018325 * u,
                                      218.020059 * u,      219.019904 * u,
                                      220.021705 * u,      221.021875 * u,
                                      222.023784 * u,      223.023963 * u,
                                     224.0256176 * u,      225.026131 * u,
                                      226.027948 * u,     227.0288054 * u,
                                     228.0310517 * u,     229.0320972 * u,
                                      230.034541 * u,     231.0358842 * u,
                                     232.0385917 * u,     233.0402472 * u,
                                     234.0433072 * u,      235.045399 * u,
                                      236.048668 * u,      237.051023 * u,
                                      238.054637 * u,       239.05726 * u,
                                       240.06098 * u,       241.06408 * u };
   static const double  PaW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 1,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double PaIsoMass[30];

   fNISTElementDataTable[90].fElemSymbol    = "Pa";
   fNISTElementDataTable[90].fSymbols       = PaS;
   fNISTElementDataTable[90].fZ             = 91;
   fNISTElementDataTable[90].fNumOfIsotopes = 30;
   fNISTElementDataTable[90].fIsNoStableIsotope       = false;
   fNISTElementDataTable[90].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[90].fNIsos         = PaN;
   fNISTElementDataTable[90].fAIsos         = PaA;
   fNISTElementDataTable[90].fWIsos         = PaW;
   fNISTElementDataTable[90].fMassIsos      = PaIsoMass;


  // Z = 92-------------------------------------------------------------------------------
   static const std::string   US[] = {  "U217",   "U218",   "U219",   "U220",   "U221",
                                        "U222",   "U223",   "U224",   "U225",   "U226",
                                        "U227",   "U228",   "U229",   "U230",   "U231",
                                        "U232",   "U233",   "U234",   "U235",   "U236",
                                        "U237",   "U238",   "U239",   "U240",   "U241",
                                        "U242",   "U243" };
   static const int      UN[] = {217, 218, 219, 220, 221, 222, 223, 224, 225, 226,
                                 227, 228, 229, 230, 231, 232, 233, 234, 235, 236,
                                 237, 238, 239, 240, 241, 242, 243 };
   static const double   UA[] = {      217.02466 * u,      218.023523 * u,
                                      219.024999 * u,       220.02462 * u,
                                       221.02628 * u,         222.026 * u,
                                      223.027739 * u,      224.027605 * u,
                                      225.029391 * u,      226.029339 * u,
                                      227.031157 * u,      228.031371 * u,
                                     229.0335063 * u,     230.0339401 * u,
                                     231.0362939 * u,     232.0371563 * u,
                                     233.0396355 * u,     234.0409523 * u,
                                     235.0439301 * u,     236.0455682 * u,
                                     237.0487304 * u,     238.0507884 * u,
                                     239.0542935 * u,     240.0565934 * u,
                                       241.06033 * u,       242.06293 * u,
                                       243.06699 * u };
   static const double   UW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 5.4e-05, 0.007204, 0.,
                                  0., 0.992742, 0., 0., 0., 0., 0. };
   static double UIsoMass[27];

   fNISTElementDataTable[91].fElemSymbol    = "U";
   fNISTElementDataTable[91].fSymbols       = US;
   fNISTElementDataTable[91].fZ             = 92;
   fNISTElementDataTable[91].fNumOfIsotopes = 27;
   fNISTElementDataTable[91].fIsNoStableIsotope       = false;
   fNISTElementDataTable[91].fIndxOfMostStableIsotope = -1;
   fNISTElementDataTable[91].fNIsos         = UN;
   fNISTElementDataTable[91].fAIsos         = UA;
   fNISTElementDataTable[91].fWIsos         = UW;
   fNISTElementDataTable[91].fMassIsos      = UIsoMass;


  // Z = 93-------------------------------------------------------------------------------
   static const std::string  NpS[] = { "Np219",  "Np220",  "Np221",  "Np222",  "Np223",
                                       "Np224",  "Np225",  "Np226",  "Np227",  "Np228",
                                       "Np229",  "Np230",  "Np231",  "Np232",  "Np233",
                                       "Np234",  "Np235",  "Np236",  "Np237",  "Np238",
                                       "Np239",  "Np240",  "Np241",  "Np242",  "Np243",
                                       "Np244",  "Np245" };
   static const int     NpN[] = {219, 220, 221, 222, 223, 224, 225, 226, 227, 228,
                                 229, 230, 231, 232, 233, 234, 235, 236, 237, 238,
                                 239, 240, 241, 242, 243, 244, 245 };
   static const double  NpA[] = {      219.03143 * u,       220.03254 * u,
                                       221.03204 * u,        222.0333 * u,
                                       223.03285 * u,       224.03422 * u,
                                      225.033911 * u,      226.035188 * u,
                                      227.034957 * u,      228.036067 * u,
                                      229.036264 * u,      230.037828 * u,
                                      231.038245 * u,       232.04011 * u,
                                      233.040741 * u,     234.0428953 * u,
                                     235.0440635 * u,       236.04657 * u,
                                     237.0481736 * u,     238.0509466 * u,
                                     239.0529392 * u,      240.056165 * u,
                                      241.058253 * u,       242.06164 * u,
                                       243.06428 * u,       244.06785 * u,
                                        245.0708 * u };
   static const double  NpW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0. };
   static double NpIsoMass[27];

   fNISTElementDataTable[92].fElemSymbol    = "Np";
   fNISTElementDataTable[92].fSymbols       = NpS;
   fNISTElementDataTable[92].fZ             = 93;
   fNISTElementDataTable[92].fNumOfIsotopes = 27;
   fNISTElementDataTable[92].fIsNoStableIsotope       = true;
   fNISTElementDataTable[92].fIndxOfMostStableIsotope = 18;
   fNISTElementDataTable[92].fNIsos         = NpN;
   fNISTElementDataTable[92].fAIsos         = NpA;
   fNISTElementDataTable[92].fWIsos         = NpW;
   fNISTElementDataTable[92].fMassIsos      = NpIsoMass;


  // Z = 94-------------------------------------------------------------------------------
   static const std::string  PuS[] = { "Pu228",  "Pu229",  "Pu230",  "Pu231",  "Pu232",
                                       "Pu233",  "Pu234",  "Pu235",  "Pu236",  "Pu237",
                                       "Pu238",  "Pu239",  "Pu240",  "Pu241",  "Pu242",
                                       "Pu243",  "Pu244",  "Pu245",  "Pu246",  "Pu247" };
   static const int     PuN[] = {228, 229, 230, 231, 232, 233, 234, 235, 236, 237,
                                 238, 239, 240, 241, 242, 243, 244, 245, 246, 247 };
   static const double  PuA[] = {     228.038732 * u,      229.040144 * u,
                                       230.03965 * u,      231.041102 * u,
                                      232.041185 * u,      233.042998 * u,
                                     234.0433174 * u,      235.045286 * u,
                                     236.0460581 * u,     237.0484098 * u,
                                     238.0495601 * u,     239.0521636 * u,
                                     240.0538138 * u,     241.0568517 * u,
                                     242.0587428 * u,     243.0620036 * u,
                                     244.0642053 * u,      245.067826 * u,
                                      246.070205 * u,       247.07419 * u };
   static const double  PuW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double PuIsoMass[20];

   fNISTElementDataTable[93].fElemSymbol    = "Pu";
   fNISTElementDataTable[93].fSymbols       = PuS;
   fNISTElementDataTable[93].fZ             = 94;
   fNISTElementDataTable[93].fNumOfIsotopes = 20;
   fNISTElementDataTable[93].fIsNoStableIsotope       = true;
   fNISTElementDataTable[93].fIndxOfMostStableIsotope = 16;
   fNISTElementDataTable[93].fNIsos         = PuN;
   fNISTElementDataTable[93].fAIsos         = PuA;
   fNISTElementDataTable[93].fWIsos         = PuW;
   fNISTElementDataTable[93].fMassIsos      = PuIsoMass;


  // Z = 95-------------------------------------------------------------------------------
   static const std::string  AmS[] = { "Am230",  "Am231",  "Am232",  "Am233",  "Am234",
                                       "Am235",  "Am236",  "Am237",  "Am238",  "Am239",
                                       "Am240",  "Am241",  "Am242",  "Am243",  "Am244",
                                       "Am245",  "Am246",  "Am247",  "Am248",  "Am249" };
   static const int     AmN[] = {230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
                                 240, 241, 242, 243, 244, 245, 246, 247, 248, 249 };
   static const double  AmA[] = {      230.04609 * u,       231.04556 * u,
                                       232.04645 * u,       233.04644 * u,
                                       234.04773 * u,      235.047908 * u,
                                       236.04943 * u,      237.049996 * u,
                                      238.051985 * u,     239.0530247 * u,
                                        240.0553 * u,     241.0568293 * u,
                                     242.0595494 * u,     243.0613813 * u,
                                     244.0642851 * u,     245.0664548 * u,
                                      246.069775 * u,       247.07209 * u,
                                       248.07575 * u,       249.07848 * u };
   static const double  AmW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double AmIsoMass[20];

   fNISTElementDataTable[94].fElemSymbol    = "Am";
   fNISTElementDataTable[94].fSymbols       = AmS;
   fNISTElementDataTable[94].fZ             = 95;
   fNISTElementDataTable[94].fNumOfIsotopes = 20;
   fNISTElementDataTable[94].fIsNoStableIsotope       = true;
   fNISTElementDataTable[94].fIndxOfMostStableIsotope = 13;
   fNISTElementDataTable[94].fNIsos         = AmN;
   fNISTElementDataTable[94].fAIsos         = AmA;
   fNISTElementDataTable[94].fWIsos         = AmW;
   fNISTElementDataTable[94].fMassIsos      = AmIsoMass;


  // Z = 96-------------------------------------------------------------------------------
   static const std::string  CmS[] = { "Cm232",  "Cm233",  "Cm234",  "Cm235",  "Cm236",
                                       "Cm237",  "Cm238",  "Cm239",  "Cm240",  "Cm241",
                                       "Cm242",  "Cm243",  "Cm244",  "Cm245",  "Cm246",
                                       "Cm247",  "Cm248",  "Cm249",  "Cm250",  "Cm251",
                                       "Cm252" };
   static const int     CmN[] = {232, 233, 234, 235, 236, 237, 238, 239, 240, 241,
                                 242, 243, 244, 245, 246, 247, 248, 249, 250, 251,
                                 252 };
   static const double  CmA[] = {      232.04982 * u,       233.05077 * u,
                                       234.05016 * u,       235.05154 * u,
                                      236.051374 * u,      237.052869 * u,
                                      238.053081 * u,       239.05491 * u,
                                     240.0555297 * u,     241.0576532 * u,
                                      242.058836 * u,     243.0613893 * u,
                                     244.0627528 * u,     245.0654915 * u,
                                     246.0672238 * u,     247.0703541 * u,
                                     248.0723499 * u,     249.0759548 * u,
                                      250.078358 * u,      251.082286 * u,
                                       252.08487 * u };
   static const double  CmW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0. };
   static double CmIsoMass[21];

   fNISTElementDataTable[95].fElemSymbol    = "Cm";
   fNISTElementDataTable[95].fSymbols       = CmS;
   fNISTElementDataTable[95].fZ             = 96;
   fNISTElementDataTable[95].fNumOfIsotopes = 21;
   fNISTElementDataTable[95].fIsNoStableIsotope       = true;
   fNISTElementDataTable[95].fIndxOfMostStableIsotope = 15;
   fNISTElementDataTable[95].fNIsos         = CmN;
   fNISTElementDataTable[95].fAIsos         = CmA;
   fNISTElementDataTable[95].fWIsos         = CmW;
   fNISTElementDataTable[95].fMassIsos      = CmIsoMass;


  // Z = 97-------------------------------------------------------------------------------
   static const std::string  BkS[] = { "Bk234",  "Bk235",  "Bk236",  "Bk237",  "Bk238",
                                       "Bk239",  "Bk240",  "Bk241",  "Bk242",  "Bk243",
                                       "Bk244",  "Bk245",  "Bk246",  "Bk247",  "Bk248",
                                       "Bk249",  "Bk250",  "Bk251",  "Bk252",  "Bk253",
                                       "Bk254" };
   static const int     BkN[] = {234, 235, 236, 237, 238, 239, 240, 241, 242, 243,
                                 244, 245, 246, 247, 248, 249, 250, 251, 252, 253,
                                 254 };
   static const double  BkA[] = {      234.05727 * u,       235.05658 * u,
                                       236.05748 * u,        237.0571 * u,
                                        238.0582 * u,       239.05824 * u,
                                       240.05976 * u,       241.06016 * u,
                                       242.06198 * u,     243.0630078 * u,
                                      244.065181 * u,     245.0663618 * u,
                                      246.068673 * u,     247.0703073 * u,
                                      248.073088 * u,     249.0749877 * u,
                                     250.0783167 * u,      251.080762 * u,
                                       252.08431 * u,       253.08688 * u,
                                        254.0906 * u };
   static const double  BkW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0. };
   static double BkIsoMass[21];

   fNISTElementDataTable[96].fElemSymbol    = "Bk";
   fNISTElementDataTable[96].fSymbols       = BkS;
   fNISTElementDataTable[96].fZ             = 97;
   fNISTElementDataTable[96].fNumOfIsotopes = 21;
   fNISTElementDataTable[96].fIsNoStableIsotope       = true;
   fNISTElementDataTable[96].fIndxOfMostStableIsotope = 13;
   fNISTElementDataTable[96].fNIsos         = BkN;
   fNISTElementDataTable[96].fAIsos         = BkA;
   fNISTElementDataTable[96].fWIsos         = BkW;
   fNISTElementDataTable[96].fMassIsos      = BkIsoMass;


  // Z = 98-------------------------------------------------------------------------------
   static const std::string  CfS[] = { "Cf237",  "Cf238",  "Cf239",  "Cf240",  "Cf241",
                                       "Cf242",  "Cf243",  "Cf244",  "Cf245",  "Cf246",
                                       "Cf247",  "Cf248",  "Cf249",  "Cf250",  "Cf251",
                                       "Cf252",  "Cf253",  "Cf254",  "Cf255",  "Cf256" };
   static const int     CfN[] = {237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
                                 247, 248, 249, 250, 251, 252, 253, 254, 255, 256 };
   static const double  CfA[] = {     237.062198 * u,       238.06149 * u,
                                       239.06253 * u,      240.062256 * u,
                                       241.06369 * u,      242.063754 * u,
                                       243.06548 * u,     244.0660008 * u,
                                     245.0680487 * u,     246.0688055 * u,
                                      247.070965 * u,     248.0721851 * u,
                                     249.0748539 * u,     250.0764062 * u,
                                     251.0795886 * u,     252.0816272 * u,
                                     253.0851345 * u,      254.087324 * u,
                                       255.09105 * u,       256.09344 * u };
   static const double  CfW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double CfIsoMass[20];

   fNISTElementDataTable[97].fElemSymbol    = "Cf";
   fNISTElementDataTable[97].fSymbols       = CfS;
   fNISTElementDataTable[97].fZ             = 98;
   fNISTElementDataTable[97].fNumOfIsotopes = 20;
   fNISTElementDataTable[97].fIsNoStableIsotope       = true;
   fNISTElementDataTable[97].fIndxOfMostStableIsotope = 14;
   fNISTElementDataTable[97].fNIsos         = CfN;
   fNISTElementDataTable[97].fAIsos         = CfA;
   fNISTElementDataTable[97].fWIsos         = CfW;
   fNISTElementDataTable[97].fMassIsos      = CfIsoMass;


  // Z = 99-------------------------------------------------------------------------------
   static const std::string  EsS[] = { "Es239",  "Es240",  "Es241",  "Es242",  "Es243",
                                       "Es244",  "Es245",  "Es246",  "Es247",  "Es248",
                                       "Es249",  "Es250",  "Es251",  "Es252",  "Es253",
                                       "Es254",  "Es255",  "Es256",  "Es257",  "Es258" };
   static const int     EsN[] = {239, 240, 241, 242, 243, 244, 245, 246, 247, 248,
                                 249, 250, 251, 252, 253, 254, 255, 256, 257, 258 };
   static const double  EsA[] = {      239.06823 * u,       240.06892 * u,
                                       241.06856 * u,       242.06957 * u,
                                       243.06951 * u,       244.07088 * u,
                                       245.07125 * u,        246.0729 * u,
                                      247.073622 * u,      248.075471 * u,
                                      249.076411 * u,       250.07861 * u,
                                     251.0799936 * u,       252.08298 * u,
                                     253.0848257 * u,     254.0880222 * u,
                                      255.090275 * u,        256.0936 * u,
                                       257.09598 * u,       258.09952 * u };
   static const double  EsW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double EsIsoMass[20];

   fNISTElementDataTable[98].fElemSymbol    = "Es";
   fNISTElementDataTable[98].fSymbols       = EsS;
   fNISTElementDataTable[98].fZ             = 99;
   fNISTElementDataTable[98].fNumOfIsotopes = 20;
   fNISTElementDataTable[98].fIsNoStableIsotope       = true;
   fNISTElementDataTable[98].fIndxOfMostStableIsotope = 13;
   fNISTElementDataTable[98].fNIsos         = EsN;
   fNISTElementDataTable[98].fAIsos         = EsA;
   fNISTElementDataTable[98].fWIsos         = EsW;
   fNISTElementDataTable[98].fMassIsos      = EsIsoMass;


  // Z = 100-------------------------------------------------------------------------------
   static const std::string  FmS[] = { "Fm241",  "Fm242",  "Fm243",  "Fm244",  "Fm245",
                                       "Fm246",  "Fm247",  "Fm248",  "Fm249",  "Fm250",
                                       "Fm251",  "Fm252",  "Fm253",  "Fm254",  "Fm255",
                                       "Fm256",  "Fm257",  "Fm258",  "Fm259",  "Fm260" };
   static const int     FmN[] = {241, 242, 243, 244, 245, 246, 247, 248, 249, 250,
                                 251, 252, 253, 254, 255, 256, 257, 258, 259, 260 };
   static const double  FmA[] = {      241.07421 * u,       242.07343 * u,
                                       243.07446 * u,       244.07404 * u,
                                       245.07535 * u,       246.07535 * u,
                                       247.07694 * u,     248.0771865 * u,
                                     249.0789275 * u,      250.079521 * u,
                                       251.08154 * u,     252.0824671 * u,
                                     253.0851846 * u,     254.0868544 * u,
                                      255.089964 * u,     256.0917745 * u,
                                     257.0951061 * u,       258.09708 * u,
                                        259.1006 * u,       260.10281 * u };
   static const double  FmW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double FmIsoMass[20];

   fNISTElementDataTable[99].fElemSymbol    = "Fm";
   fNISTElementDataTable[99].fSymbols       = FmS;
   fNISTElementDataTable[99].fZ             = 100;
   fNISTElementDataTable[99].fNumOfIsotopes = 20;
   fNISTElementDataTable[99].fIsNoStableIsotope       = true;
   fNISTElementDataTable[99].fIndxOfMostStableIsotope = 16;
   fNISTElementDataTable[99].fNIsos         = FmN;
   fNISTElementDataTable[99].fAIsos         = FmA;
   fNISTElementDataTable[99].fWIsos         = FmW;
   fNISTElementDataTable[99].fMassIsos      = FmIsoMass;


  // Z = 101-------------------------------------------------------------------------------
   static const std::string  MdS[] = { "Md245",  "Md246",  "Md247",  "Md248",  "Md249",
                                       "Md250",  "Md251",  "Md252",  "Md253",  "Md254",
                                       "Md255",  "Md256",  "Md257",  "Md258",  "Md259",
                                       "Md260",  "Md261",  "Md262" };
   static const int     MdN[] = {245, 246, 247, 248, 249, 250, 251, 252, 253, 254,
                                 255, 256, 257, 258, 259, 260, 261, 262 };
   static const double  MdA[] = {      245.08081 * u,       246.08171 * u,
                                       247.08152 * u,       248.08282 * u,
                                       249.08291 * u,       250.08441 * u,
                                      251.084774 * u,       252.08643 * u,
                                      253.087144 * u,       254.08959 * u,
                                     255.0910841 * u,       256.09389 * u,
                                     257.0955424 * u,     258.0984315 * u,
                                       259.10051 * u,       260.10365 * u,
                                       261.10583 * u,        262.1091 * u };
   static const double  MdW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0., 0. };
   static double MdIsoMass[18];

   fNISTElementDataTable[100].fElemSymbol    = "Md";
   fNISTElementDataTable[100].fSymbols       = MdS;
   fNISTElementDataTable[100].fZ             = 101;
   fNISTElementDataTable[100].fNumOfIsotopes = 18;
   fNISTElementDataTable[100].fIsNoStableIsotope       = true;
   fNISTElementDataTable[100].fIndxOfMostStableIsotope = 13;
   fNISTElementDataTable[100].fNIsos         = MdN;
   fNISTElementDataTable[100].fAIsos         = MdA;
   fNISTElementDataTable[100].fWIsos         = MdW;
   fNISTElementDataTable[100].fMassIsos      = MdIsoMass;


  // Z = 102-------------------------------------------------------------------------------
   static const std::string  NoS[] = { "No248",  "No249",  "No250",  "No251",  "No252",
                                       "No253",  "No254",  "No255",  "No256",  "No257",
                                       "No258",  "No259",  "No260",  "No261",  "No262",
                                       "No263",  "No264" };
   static const int     NoN[] = {248, 249, 250, 251, 252, 253, 254, 255, 256, 257,
                                 258, 259, 260, 261, 262, 263, 264 };
   static const double  NoA[] = {      248.08655 * u,        249.0878 * u,
                                       250.08756 * u,       251.08894 * u,
                                      252.088967 * u,     253.0905641 * u,
                                      254.090956 * u,      255.093191 * u,
                                     256.0942829 * u,     257.0968878 * u,
                                       258.09821 * u,       259.10103 * u,
                                       260.10264 * u,        261.1057 * u,
                                       262.10746 * u,       263.11071 * u,
                                       264.11273 * u };
   static const double  NoW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0., 0. };
   static double NoIsoMass[17];

   fNISTElementDataTable[101].fElemSymbol    = "No";
   fNISTElementDataTable[101].fSymbols       = NoS;
   fNISTElementDataTable[101].fZ             = 102;
   fNISTElementDataTable[101].fNumOfIsotopes = 17;
   fNISTElementDataTable[101].fIsNoStableIsotope       = true;
   fNISTElementDataTable[101].fIndxOfMostStableIsotope = 11;
   fNISTElementDataTable[101].fNIsos         = NoN;
   fNISTElementDataTable[101].fAIsos         = NoA;
   fNISTElementDataTable[101].fWIsos         = NoW;
   fNISTElementDataTable[101].fMassIsos      = NoIsoMass;


  // Z = 103-------------------------------------------------------------------------------
   static const std::string  LrS[] = { "Lr251",  "Lr252",  "Lr253",  "Lr254",  "Lr255",
                                       "Lr256",  "Lr257",  "Lr258",  "Lr259",  "Lr260",
                                       "Lr261",  "Lr262",  "Lr263",  "Lr264",  "Lr265",
                                       "Lr266" };
   static const int     LrN[] = {251, 252, 253, 254, 255, 256, 257, 258, 259, 260,
                                 261, 262, 263, 264, 265, 266 };
   static const double  LrA[] = {      251.09418 * u,       252.09526 * u,
                                       253.09509 * u,       254.09648 * u,
                                      255.096562 * u,      256.098494 * u,
                                      257.099418 * u,       258.10176 * u,
                                      259.102902 * u,        260.1055 * u,
                                       261.10688 * u,       262.10961 * u,
                                       263.11136 * u,        264.1142 * u,
                                       265.11619 * u,       266.11983 * u };
   static const double  LrW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double LrIsoMass[16];

   fNISTElementDataTable[102].fElemSymbol    = "Lr";
   fNISTElementDataTable[102].fSymbols       = LrS;
   fNISTElementDataTable[102].fZ             = 103;
   fNISTElementDataTable[102].fNumOfIsotopes = 16;
   fNISTElementDataTable[102].fIsNoStableIsotope       = true;
   fNISTElementDataTable[102].fIndxOfMostStableIsotope = 15;
   fNISTElementDataTable[102].fNIsos         = LrN;
   fNISTElementDataTable[102].fAIsos         = LrA;
   fNISTElementDataTable[102].fWIsos         = LrW;
   fNISTElementDataTable[102].fMassIsos      = LrIsoMass;


  // Z = 104-------------------------------------------------------------------------------
   static const std::string  RfS[] = { "Rf253",  "Rf254",  "Rf255",  "Rf256",  "Rf257",
                                       "Rf258",  "Rf259",  "Rf260",  "Rf261",  "Rf262",
                                       "Rf263",  "Rf264",  "Rf265",  "Rf266",  "Rf267",
                                       "Rf268" };
   static const int     RfN[] = {253, 254, 255, 256, 257, 258, 259, 260, 261, 262,
                                 263, 264, 265, 266, 267, 268 };
   static const double  RfA[] = {      253.10044 * u,       254.10005 * u,
                                       255.10127 * u,      256.101152 * u,
                                      257.102918 * u,      258.103428 * u,
                                      259.105596 * u,       260.10644 * u,
                                      261.108773 * u,       262.10992 * u,
                                       263.11249 * u,       264.11388 * u,
                                       265.11668 * u,       266.11817 * u,
                                       267.12179 * u,       268.12397 * u };
   static const double  RfW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double RfIsoMass[16];

   fNISTElementDataTable[103].fElemSymbol    = "Rf";
   fNISTElementDataTable[103].fSymbols       = RfS;
   fNISTElementDataTable[103].fZ             = 104;
   fNISTElementDataTable[103].fNumOfIsotopes = 16;
   fNISTElementDataTable[103].fIsNoStableIsotope       = true;
   fNISTElementDataTable[103].fIndxOfMostStableIsotope = 14;
   fNISTElementDataTable[103].fNIsos         = RfN;
   fNISTElementDataTable[103].fAIsos         = RfA;
   fNISTElementDataTable[103].fWIsos         = RfW;
   fNISTElementDataTable[103].fMassIsos      = RfIsoMass;


  // Z = 105-------------------------------------------------------------------------------
   static const std::string  DbS[] = { "Db255",  "Db256",  "Db257",  "Db258",  "Db259",
                                       "Db260",  "Db261",  "Db262",  "Db263",  "Db264",
                                       "Db265",  "Db266",  "Db267",  "Db268",  "Db269",
                                       "Db270" };
   static const int     DbN[] = {255, 256, 257, 258, 259, 260, 261, 262, 263, 264,
                                 265, 266, 267, 268, 269, 270 };
   static const double  DbA[] = {      255.10707 * u,       256.10789 * u,
                                       257.10758 * u,       258.10928 * u,
                                      259.109492 * u,        260.1113 * u,
                                       261.11192 * u,       262.11407 * u,
                                       263.11499 * u,       264.11741 * u,
                                       265.11861 * u,       266.12103 * u,
                                       267.12247 * u,       268.12567 * u,
                                       269.12791 * u,       270.13136 * u };
   static const double  DbW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double DbIsoMass[16];

   fNISTElementDataTable[104].fElemSymbol    = "Db";
   fNISTElementDataTable[104].fSymbols       = DbS;
   fNISTElementDataTable[104].fZ             = 105;
   fNISTElementDataTable[104].fNumOfIsotopes = 16;
   fNISTElementDataTable[104].fIsNoStableIsotope       = true;
   fNISTElementDataTable[104].fIndxOfMostStableIsotope = 13;
   fNISTElementDataTable[104].fNIsos         = DbN;
   fNISTElementDataTable[104].fAIsos         = DbA;
   fNISTElementDataTable[104].fWIsos         = DbW;
   fNISTElementDataTable[104].fMassIsos      = DbIsoMass;


  // Z = 106-------------------------------------------------------------------------------
   static const std::string  SgS[] = { "Sg258",  "Sg259",  "Sg260",  "Sg261",  "Sg262",
                                       "Sg263",  "Sg264",  "Sg265",  "Sg266",  "Sg267",
                                       "Sg268",  "Sg269",  "Sg270",  "Sg271",  "Sg272",
                                       "Sg273" };
   static const int     SgN[] = {258, 259, 260, 261, 262, 263, 264, 265, 266, 267,
                                 268, 269, 270, 271, 272, 273 };
   static const double  SgA[] = {      258.11298 * u,        259.1144 * u,
                                      260.114384 * u,      261.115949 * u,
                                      262.116337 * u,       263.11829 * u,
                                       264.11893 * u,       265.12109 * u,
                                       266.12198 * u,       267.12436 * u,
                                       268.12539 * u,       269.12863 * u,
                                       270.13043 * u,       271.13393 * u,
                                       272.13589 * u,       273.13958 * u };
   static const double  SgW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double SgIsoMass[16];

   fNISTElementDataTable[105].fElemSymbol    = "Sg";
   fNISTElementDataTable[105].fSymbols       = SgS;
   fNISTElementDataTable[105].fZ             = 106;
   fNISTElementDataTable[105].fNumOfIsotopes = 16;
   fNISTElementDataTable[105].fIsNoStableIsotope       = true;
   fNISTElementDataTable[105].fIndxOfMostStableIsotope = 11;
   fNISTElementDataTable[105].fNIsos         = SgN;
   fNISTElementDataTable[105].fAIsos         = SgA;
   fNISTElementDataTable[105].fWIsos         = SgW;
   fNISTElementDataTable[105].fMassIsos      = SgIsoMass;


  // Z = 107-------------------------------------------------------------------------------
   static const std::string  BhS[] = { "Bh260",  "Bh261",  "Bh262",  "Bh263",  "Bh264",
                                       "Bh265",  "Bh266",  "Bh267",  "Bh268",  "Bh269",
                                       "Bh270",  "Bh271",  "Bh272",  "Bh273",  "Bh274",
                                       "Bh275" };
   static const int     BhN[] = {260, 261, 262, 263, 264, 265, 266, 267, 268, 269,
                                 270, 271, 272, 273, 274, 275 };
   static const double  BhA[] = {      260.12166 * u,       261.12145 * u,
                                       262.12297 * u,       263.12292 * u,
                                       264.12459 * u,       265.12491 * u,
                                       266.12679 * u,        267.1275 * u,
                                       268.12969 * u,       269.13042 * u,
                                       270.13336 * u,       271.13526 * u,
                                       272.13826 * u,       273.14024 * u,
                                       274.14355 * u,       275.14567 * u };
   static const double  BhW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0. };
   static double BhIsoMass[16];

   fNISTElementDataTable[106].fElemSymbol    = "Bh";
   fNISTElementDataTable[106].fSymbols       = BhS;
   fNISTElementDataTable[106].fZ             = 107;
   fNISTElementDataTable[106].fNumOfIsotopes = 16;
   fNISTElementDataTable[106].fIsNoStableIsotope       = true;
   fNISTElementDataTable[106].fIndxOfMostStableIsotope = 10;
   fNISTElementDataTable[106].fNIsos         = BhN;
   fNISTElementDataTable[106].fAIsos         = BhA;
   fNISTElementDataTable[106].fWIsos         = BhW;
   fNISTElementDataTable[106].fMassIsos      = BhIsoMass;


  // Z = 108-------------------------------------------------------------------------------
   static const std::string  HsS[] = { "Hs263",  "Hs264",  "Hs265",  "Hs266",  "Hs267",
                                       "Hs268",  "Hs269",  "Hs270",  "Hs271",  "Hs272",
                                       "Hs273",  "Hs274",  "Hs275",  "Hs276",  "Hs277" };
   static const int     HsN[] = {263, 264, 265, 266, 267, 268, 269, 270, 271, 272,
                                 273, 274, 275, 276, 277 };
   static const double  HsA[] = {      263.12852 * u,      264.128357 * u,
                                      265.129793 * u,      266.130046 * u,
                                       267.13167 * u,       268.13186 * u,
                                       269.13375 * u,       270.13429 * u,
                                       271.13717 * u,        272.1385 * u,
                                       273.14168 * u,        274.1433 * u,
                                       275.14667 * u,       276.14846 * u,
                                        277.1519 * u };
   static const double  HsW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double HsIsoMass[15];

   fNISTElementDataTable[107].fElemSymbol    = "Hs";
   fNISTElementDataTable[107].fSymbols       = HsS;
   fNISTElementDataTable[107].fZ             = 108;
   fNISTElementDataTable[107].fNumOfIsotopes = 15;
   fNISTElementDataTable[107].fIsNoStableIsotope       = true;
   fNISTElementDataTable[107].fIndxOfMostStableIsotope = 6;
   fNISTElementDataTable[107].fNIsos         = HsN;
   fNISTElementDataTable[107].fAIsos         = HsA;
   fNISTElementDataTable[107].fWIsos         = HsW;
   fNISTElementDataTable[107].fMassIsos      = HsIsoMass;


  // Z = 109-------------------------------------------------------------------------------
   static const std::string  MtS[] = { "Mt265",  "Mt266",  "Mt267",  "Mt268",  "Mt269",
                                       "Mt270",  "Mt271",  "Mt272",  "Mt273",  "Mt274",
                                       "Mt275",  "Mt276",  "Mt277",  "Mt278",  "Mt279" };
   static const int     MtN[] = {265, 266, 267, 268, 269, 270, 271, 272, 273, 274,
                                 275, 276, 277, 278, 279 };
   static const double  MtA[] = {        265.136 * u,       266.13737 * u,
                                       267.13719 * u,       268.13865 * u,
                                       269.13882 * u,       270.14033 * u,
                                       271.14074 * u,       272.14341 * u,
                                        273.1444 * u,       274.14724 * u,
                                       275.14882 * u,       276.15159 * u,
                                       277.15327 * u,       278.15631 * u,
                                       279.15808 * u };
   static const double  MtW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double MtIsoMass[15];

   fNISTElementDataTable[108].fElemSymbol    = "Mt";
   fNISTElementDataTable[108].fSymbols       = MtS;
   fNISTElementDataTable[108].fZ             = 109;
   fNISTElementDataTable[108].fNumOfIsotopes = 15;
   fNISTElementDataTable[108].fIsNoStableIsotope       = true;
   fNISTElementDataTable[108].fIndxOfMostStableIsotope = 13;
   fNISTElementDataTable[108].fNIsos         = MtN;
   fNISTElementDataTable[108].fAIsos         = MtA;
   fNISTElementDataTable[108].fWIsos         = MtW;
   fNISTElementDataTable[108].fMassIsos      = MtIsoMass;


  // Z = 110-------------------------------------------------------------------------------
   static const std::string  DsS[] = { "Ds267",  "Ds268",  "Ds269",  "Ds270",  "Ds271",
                                       "Ds272",  "Ds273",  "Ds274",  "Ds275",  "Ds276",
                                       "Ds277",  "Ds278",  "Ds279",  "Ds280",  "Ds281" };
   static const int     DsN[] = {267, 268, 269, 270, 271, 272, 273, 274, 275, 276,
                                 277, 278, 279, 280, 281 };
   static const double  DsA[] = {      267.14377 * u,       268.14348 * u,
                                      269.144752 * u,      270.144584 * u,
                                       271.14595 * u,       272.14602 * u,
                                       273.14856 * u,       274.14941 * u,
                                       275.15203 * u,       276.15303 * u,
                                       277.15591 * u,       278.15704 * u,
                                        279.1601 * u,       280.16131 * u,
                                       281.16451 * u };
   static const double  DsW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0. };
   static double DsIsoMass[15];

   fNISTElementDataTable[109].fElemSymbol    = "Ds";
   fNISTElementDataTable[109].fSymbols       = DsS;
   fNISTElementDataTable[109].fZ             = 110;
   fNISTElementDataTable[109].fNumOfIsotopes = 15;
   fNISTElementDataTable[109].fIsNoStableIsotope       = true;
   fNISTElementDataTable[109].fIndxOfMostStableIsotope = 14;
   fNISTElementDataTable[109].fNIsos         = DsN;
   fNISTElementDataTable[109].fAIsos         = DsA;
   fNISTElementDataTable[109].fWIsos         = DsW;
   fNISTElementDataTable[109].fMassIsos      = DsIsoMass;


  // Z = 111-------------------------------------------------------------------------------
   static const std::string  RgS[] = { "Rg272",  "Rg273",  "Rg274",  "Rg275",  "Rg276",
                                       "Rg277",  "Rg278",  "Rg279",  "Rg280",  "Rg281",
                                       "Rg282",  "Rg283" };
   static const int     RgN[] = {272, 273, 274, 275, 276, 277, 278, 279, 280, 281,
                                 282, 283 };
   static const double  RgA[] = {      272.15327 * u,       273.15313 * u,
                                       274.15525 * u,       275.15594 * u,
                                       276.15833 * u,       277.15907 * u,
                                       278.16149 * u,       279.16272 * u,
                                       280.16514 * u,       281.16636 * u,
                                       282.16912 * u,       283.17054 * u };
   static const double  RgW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                  0., 0. };
   static double RgIsoMass[12];

   fNISTElementDataTable[110].fElemSymbol    = "Rg";
   fNISTElementDataTable[110].fSymbols       = RgS;
   fNISTElementDataTable[110].fZ             = 111;
   fNISTElementDataTable[110].fNumOfIsotopes = 12;
   fNISTElementDataTable[110].fIsNoStableIsotope       = true;
   fNISTElementDataTable[110].fIndxOfMostStableIsotope = 10;
   fNISTElementDataTable[110].fNIsos         = RgN;
   fNISTElementDataTable[110].fAIsos         = RgA;
   fNISTElementDataTable[110].fWIsos         = RgW;
   fNISTElementDataTable[110].fMassIsos      = RgIsoMass;


  // Z = 112-------------------------------------------------------------------------------
   static const std::string  CnS[] = { "Cn276",  "Cn277",  "Cn278",  "Cn279",  "Cn280",
                                       "Cn281",  "Cn282",  "Cn283",  "Cn284",  "Cn285" };
   static const int     CnN[] = {276, 277, 278, 279, 280, 281, 282, 283, 284, 285 };
   static const double  CnA[] = {      276.16141 * u,       277.16364 * u,
                                       278.16416 * u,       279.16654 * u,
                                       280.16715 * u,       281.16975 * u,
                                        282.1705 * u,       283.17327 * u,
                                       284.17416 * u,       285.17712 * u };
   static const double  CnW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double CnIsoMass[10];

   fNISTElementDataTable[111].fElemSymbol    = "Cn";
   fNISTElementDataTable[111].fSymbols       = CnS;
   fNISTElementDataTable[111].fZ             = 112;
   fNISTElementDataTable[111].fNumOfIsotopes = 10;
   fNISTElementDataTable[111].fIsNoStableIsotope       = true;
   fNISTElementDataTable[111].fIndxOfMostStableIsotope = 9;
   fNISTElementDataTable[111].fNIsos         = CnN;
   fNISTElementDataTable[111].fAIsos         = CnA;
   fNISTElementDataTable[111].fWIsos         = CnW;
   fNISTElementDataTable[111].fMassIsos      = CnIsoMass;


  // Z = 113-------------------------------------------------------------------------------
   static const std::string UutS[] = {"Uut278", "Uut279", "Uut280", "Uut281", "Uut282",
                                      "Uut283", "Uut284", "Uut285", "Uut286", "Uut287" };
   static const int    UutN[] = {278, 279, 280, 281, 282, 283, 284, 285, 286, 287 };
   static const double UutA[] = {      278.17058 * u,       279.17095 * u,
                                       280.17293 * u,       281.17348 * u,
                                       282.17567 * u,       283.17657 * u,
                                       284.17873 * u,       285.17973 * u,
                                       286.18221 * u,       287.18339 * u };
   static const double UutW[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
   static double UutIsoMass[10];

   fNISTElementDataTable[112].fElemSymbol    = "Uut";
   fNISTElementDataTable[112].fSymbols       = UutS;
   fNISTElementDataTable[112].fZ             = 113;
   fNISTElementDataTable[112].fNumOfIsotopes = 10;
   fNISTElementDataTable[112].fIsNoStableIsotope       = true;
   fNISTElementDataTable[112].fIndxOfMostStableIsotope = 8;
   fNISTElementDataTable[112].fNIsos         = UutN;
   fNISTElementDataTable[112].fAIsos         = UutA;
   fNISTElementDataTable[112].fWIsos         = UutW;
   fNISTElementDataTable[112].fMassIsos      = UutIsoMass;


  // Z = 114-------------------------------------------------------------------------------
   static const std::string  FlS[] = { "Fl285",  "Fl286",  "Fl287",  "Fl288",  "Fl289" };
   static const int     FlN[] = {285, 286, 287, 288, 289 };
   static const double  FlA[] = {      285.18364 * u,       286.18423 * u,
                                       287.18678 * u,       288.18757 * u,
                                       289.19042 * u };
   static const double  FlW[] = { 0., 0., 0., 0., 0. };
   static double FlIsoMass[5];

   fNISTElementDataTable[113].fElemSymbol    = "Fl";
   fNISTElementDataTable[113].fSymbols       = FlS;
   fNISTElementDataTable[113].fZ             = 114;
   fNISTElementDataTable[113].fNumOfIsotopes = 5;
   fNISTElementDataTable[113].fIsNoStableIsotope       = true;
   fNISTElementDataTable[113].fIndxOfMostStableIsotope = 4;
   fNISTElementDataTable[113].fNIsos         = FlN;
   fNISTElementDataTable[113].fAIsos         = FlA;
   fNISTElementDataTable[113].fWIsos         = FlW;
   fNISTElementDataTable[113].fMassIsos      = FlIsoMass;


  // Z = 115-------------------------------------------------------------------------------
   static const std::string UupS[] = {"Uup287", "Uup288", "Uup289", "Uup290", "Uup291" };
   static const int    UupN[] = {287, 288, 289, 290, 291 };
   static const double UupA[] = {       287.1907 * u,       288.19274 * u,
                                       289.19363 * u,       290.19598 * u,
                                       291.19707 * u };
   static const double UupW[] = { 0., 0., 0., 0., 0. };
   static double UupIsoMass[5];

   fNISTElementDataTable[114].fElemSymbol    = "Uup";
   fNISTElementDataTable[114].fSymbols       = UupS;
   fNISTElementDataTable[114].fZ             = 115;
   fNISTElementDataTable[114].fNumOfIsotopes = 5;
   fNISTElementDataTable[114].fIsNoStableIsotope       = true;
   fNISTElementDataTable[114].fIndxOfMostStableIsotope = 2;
   fNISTElementDataTable[114].fNIsos         = UupN;
   fNISTElementDataTable[114].fAIsos         = UupA;
   fNISTElementDataTable[114].fWIsos         = UupW;
   fNISTElementDataTable[114].fMassIsos      = UupIsoMass;


  // Z = 116-------------------------------------------------------------------------------
   static const std::string  LvS[] = { "Lv289",  "Lv290",  "Lv291",  "Lv292",  "Lv293" };
   static const int     LvN[] = {289, 290, 291, 292, 293 };
   static const double  LvA[] = {      289.19816 * u,       290.19864 * u,
                                       291.20108 * u,       292.20174 * u,
                                       293.20449 * u };
   static const double  LvW[] = { 0., 0., 0., 0., 0. };
   static double LvIsoMass[5];

   fNISTElementDataTable[115].fElemSymbol    = "Lv";
   fNISTElementDataTable[115].fSymbols       = LvS;
   fNISTElementDataTable[115].fZ             = 116;
   fNISTElementDataTable[115].fNumOfIsotopes = 5;
   fNISTElementDataTable[115].fIsNoStableIsotope       = true;
   fNISTElementDataTable[115].fIndxOfMostStableIsotope = 4;
   fNISTElementDataTable[115].fNIsos         = LvN;
   fNISTElementDataTable[115].fAIsos         = LvA;
   fNISTElementDataTable[115].fWIsos         = LvW;
   fNISTElementDataTable[115].fMassIsos      = LvIsoMass;


  // Z = 117-------------------------------------------------------------------------------
   static const std::string UusS[] = {"Uus291", "Uus292", "Uus293", "Uus294" };
   static const int    UusN[] = {291, 292, 293, 294 };
   static const double UusA[] = {      291.20553 * u,       292.20746 * u,
                                       293.20824 * u,       294.21046 * u };
   static const double UusW[] = { 0., 0., 0., 0. };
   static double UusIsoMass[4];

   fNISTElementDataTable[116].fElemSymbol    = "Uus";
   fNISTElementDataTable[116].fSymbols       = UusS;
   fNISTElementDataTable[116].fZ             = 117;
   fNISTElementDataTable[116].fNumOfIsotopes = 4;
   fNISTElementDataTable[116].fIsNoStableIsotope       = true;
   fNISTElementDataTable[116].fIndxOfMostStableIsotope = 3;
   fNISTElementDataTable[116].fNIsos         = UusN;
   fNISTElementDataTable[116].fAIsos         = UusA;
   fNISTElementDataTable[116].fWIsos         = UusW;
   fNISTElementDataTable[116].fMassIsos      = UusIsoMass;


  // Z = 118-------------------------------------------------------------------------------
   static const std::string UuoS[] = {"Uuo293", "Uuo294", "Uuo295" };
   static const int    UuoN[] = {293, 294, 295 };
   static const double UuoA[] = {      293.21356 * u,       294.21392 * u,
                                       295.21624 * u };
   static const double UuoW[] = { 0., 0., 0. };
   static double UuoIsoMass[3];

   fNISTElementDataTable[117].fElemSymbol    = "Uuo";
   fNISTElementDataTable[117].fSymbols       = UuoS;
   fNISTElementDataTable[117].fZ             = 118;
   fNISTElementDataTable[117].fNumOfIsotopes = 3;
   fNISTElementDataTable[117].fIsNoStableIsotope       = true;
   fNISTElementDataTable[117].fIndxOfMostStableIsotope = 1;
   fNISTElementDataTable[117].fNIsos         = UuoN;
   fNISTElementDataTable[117].fAIsos         = UuoA;
   fNISTElementDataTable[117].fWIsos         = UuoW;
   fNISTElementDataTable[117].fMassIsos      = UuoIsoMass;

   // Compute mass of the isotopes from the atomic masses
   // i.e. without the electrons in internal [Energy] units
   for (int i = 0; i < gNumberOfNISTElements; i++) {
     int numisos = fNISTElementDataTable[i].fNumOfIsotopes;
     int zet     = fNISTElementDataTable[i].fZ;
     double be   = fBindingEnergies[zet-1];
     for (int j = 0; j < numisos; ++j) {
       fNISTElementDataTable[i].fMassIsos[j] =
       fNISTElementDataTable[i].fAIsos[j]*kCLightSquare - zet*kElectronMassC2 + be;
     }
   }

   // Compute mean atomic masses of the elements in internal [weight/mole] units
   for (int i = 0; i < gNumberOfNISTElements; i++) {
     int numisos        = fNISTElementDataTable[i].fNumOfIsotopes;
     fNISTElementDataTable[i].fMeanAtomicMass = 0.0;
     //double amass = 0.0;
     // if the element has stable isotope
     if (!fNISTElementDataTable[i].fIsNoStableIsotope) {
       for (int j = 0; j < numisos; ++j) {
         fNISTElementDataTable[i].fMeanAtomicMass +=
         fNISTElementDataTable[i].fWIsos[j]*fNISTElementDataTable[i].fAIsos[j];
       }
     } else {
       fNISTElementDataTable[i].fMeanAtomicMass =
       fNISTElementDataTable[i].fAIsos[fNISTElementDataTable[i].fIndxOfMostStableIsotope];
     }     fNISTElementDataTable[i].fMeanAtomicMass *= kAvogadro;
   }


} // namespace geantphysics

}
