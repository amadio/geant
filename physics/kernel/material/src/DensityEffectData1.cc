
#include "DensityEffectData.h"
#include "Geant/PhysicalConstants.h"
#include "MaterialState.h"

namespace geantphysics {

void DensityEffectData::BuildTable() {
  using geant::units::eV;

  // NIST_MAT_lH2 --------------------------------------------------------------------------
  fDensityEffectDataTable[  0].fName              = "NIST_MAT_lH2";
  fDensityEffectDataTable[  0].fPlasmaEnergy      =         7.031 * eV;
  fDensityEffectDataTable[  0].fSternheimerFactor =         1.546;
  fDensityEffectDataTable[  0].fParameterC        =        3.2632;
  fDensityEffectDataTable[  0].fParameterFitX0    =        0.4759;
  fDensityEffectDataTable[  0].fParameterFitX1    =        1.9215;
  fDensityEffectDataTable[  0].fParameterFitA     =       0.13483;
  fDensityEffectDataTable[  0].fParameterFitM     =        5.6249;
  fDensityEffectDataTable[  0].fParameterDelta0   =             0;
  fDensityEffectDataTable[  0].fDeltaErrorMax     =         0.021;
  fDensityEffectDataTable[  0].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_H ----------------------------------------------------------------------------
  fDensityEffectDataTable[  1].fName              = "NIST_MAT_H";
  fDensityEffectDataTable[  1].fPlasmaEnergy      =         0.263 * eV;
  fDensityEffectDataTable[  1].fSternheimerFactor =         1.412;
  fDensityEffectDataTable[  1].fParameterC        =        9.5835;
  fDensityEffectDataTable[  1].fParameterFitX0    =        1.8639;
  fDensityEffectDataTable[  1].fParameterFitX1    =        3.2718;
  fDensityEffectDataTable[  1].fParameterFitA     =       0.14092;
  fDensityEffectDataTable[  1].fParameterFitM     =        5.7273;
  fDensityEffectDataTable[  1].fParameterDelta0   =             0;
  fDensityEffectDataTable[  1].fDeltaErrorMax     =         0.024;
  fDensityEffectDataTable[  1].fState             = MaterialState::kStateGas;


  // NIST_MAT_He ---------------------------------------------------------------------------
  fDensityEffectDataTable[  2].fName              = "NIST_MAT_He";
  fDensityEffectDataTable[  2].fPlasmaEnergy      =         0.263 * eV;
  fDensityEffectDataTable[  2].fSternheimerFactor =           1.7;
  fDensityEffectDataTable[  2].fParameterC        =       11.1393;
  fDensityEffectDataTable[  2].fParameterFitX0    =        2.2017;
  fDensityEffectDataTable[  2].fParameterFitX1    =        3.6122;
  fDensityEffectDataTable[  2].fParameterFitA     =       0.13443;
  fDensityEffectDataTable[  2].fParameterFitM     =        5.8347;
  fDensityEffectDataTable[  2].fParameterDelta0   =             0;
  fDensityEffectDataTable[  2].fDeltaErrorMax     =         0.024;
  fDensityEffectDataTable[  2].fState             = MaterialState::kStateGas;


  // NIST_MAT_Li ---------------------------------------------------------------------------
  fDensityEffectDataTable[  3].fName              = "NIST_MAT_Li";
  fDensityEffectDataTable[  3].fPlasmaEnergy      =        13.844 * eV;
  fDensityEffectDataTable[  3].fSternheimerFactor =         1.535;
  fDensityEffectDataTable[  3].fParameterC        =        3.1221;
  fDensityEffectDataTable[  3].fParameterFitX0    =        0.1304;
  fDensityEffectDataTable[  3].fParameterFitX1    =        1.6397;
  fDensityEffectDataTable[  3].fParameterFitA     =       0.95136;
  fDensityEffectDataTable[  3].fParameterFitM     =        2.4993;
  fDensityEffectDataTable[  3].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[  3].fDeltaErrorMax     =         0.062;
  fDensityEffectDataTable[  3].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Be ---------------------------------------------------------------------------
  fDensityEffectDataTable[  4].fName              = "NIST_MAT_Be";
  fDensityEffectDataTable[  4].fPlasmaEnergy      =        26.096 * eV;
  fDensityEffectDataTable[  4].fSternheimerFactor =         1.908;
  fDensityEffectDataTable[  4].fParameterC        =        2.7847;
  fDensityEffectDataTable[  4].fParameterFitX0    =        0.0392;
  fDensityEffectDataTable[  4].fParameterFitX1    =        1.6922;
  fDensityEffectDataTable[  4].fParameterFitA     =       0.80392;
  fDensityEffectDataTable[  4].fParameterFitM     =        2.4339;
  fDensityEffectDataTable[  4].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[  4].fDeltaErrorMax     =         0.029;
  fDensityEffectDataTable[  4].fState             = MaterialState::kStateSolid;


  // NIST_MAT_B ----------------------------------------------------------------------------
  fDensityEffectDataTable[  5].fName              = "NIST_MAT_B";
  fDensityEffectDataTable[  5].fPlasmaEnergy      =         30.17 * eV;
  fDensityEffectDataTable[  5].fSternheimerFactor =          2.32;
  fDensityEffectDataTable[  5].fParameterC        =        2.8477;
  fDensityEffectDataTable[  5].fParameterFitX0    =        0.0305;
  fDensityEffectDataTable[  5].fParameterFitX1    =        1.9688;
  fDensityEffectDataTable[  5].fParameterFitA     =       0.56224;
  fDensityEffectDataTable[  5].fParameterFitM     =        2.4512;
  fDensityEffectDataTable[  5].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[  5].fDeltaErrorMax     =         0.024;
  fDensityEffectDataTable[  5].fState             = MaterialState::kStateSolid;


  // NIST_MAT_C ----------------------------------------------------------------------------
  fDensityEffectDataTable[  6].fName              = "NIST_MAT_C";
  fDensityEffectDataTable[  6].fPlasmaEnergy      =        28.803 * eV;
  fDensityEffectDataTable[  6].fSternheimerFactor =         2.376;
  fDensityEffectDataTable[  6].fParameterC        =        2.9925;
  fDensityEffectDataTable[  6].fParameterFitX0    =       -0.0351;
  fDensityEffectDataTable[  6].fParameterFitX1    =         2.486;
  fDensityEffectDataTable[  6].fParameterFitA     =        0.2024;
  fDensityEffectDataTable[  6].fParameterFitM     =        3.0036;
  fDensityEffectDataTable[  6].fParameterDelta0   =           0.1;
  fDensityEffectDataTable[  6].fDeltaErrorMax     =         0.038;
  fDensityEffectDataTable[  6].fState             = MaterialState::kStateSolid;


  // NIST_MAT_N ----------------------------------------------------------------------------
  fDensityEffectDataTable[  7].fName              = "NIST_MAT_N";
  fDensityEffectDataTable[  7].fPlasmaEnergy      =         0.695 * eV;
  fDensityEffectDataTable[  7].fSternheimerFactor =         1.984;
  fDensityEffectDataTable[  7].fParameterC        =         10.54;
  fDensityEffectDataTable[  7].fParameterFitX0    =        1.7378;
  fDensityEffectDataTable[  7].fParameterFitX1    =        4.1323;
  fDensityEffectDataTable[  7].fParameterFitA     =       0.15349;
  fDensityEffectDataTable[  7].fParameterFitM     =        3.2125;
  fDensityEffectDataTable[  7].fParameterDelta0   =             0;
  fDensityEffectDataTable[  7].fDeltaErrorMax     =         0.086;
  fDensityEffectDataTable[  7].fState             = MaterialState::kStateGas;


  // NIST_MAT_O ----------------------------------------------------------------------------
  fDensityEffectDataTable[  8].fName              = "NIST_MAT_O";
  fDensityEffectDataTable[  8].fPlasmaEnergy      =         0.744 * eV;
  fDensityEffectDataTable[  8].fSternheimerFactor =         2.314;
  fDensityEffectDataTable[  8].fParameterC        =       10.7004;
  fDensityEffectDataTable[  8].fParameterFitX0    =        1.7541;
  fDensityEffectDataTable[  8].fParameterFitX1    =        4.3213;
  fDensityEffectDataTable[  8].fParameterFitA     =       0.11778;
  fDensityEffectDataTable[  8].fParameterFitM     =        3.2913;
  fDensityEffectDataTable[  8].fParameterDelta0   =             0;
  fDensityEffectDataTable[  8].fDeltaErrorMax     =         0.101;
  fDensityEffectDataTable[  8].fState             = MaterialState::kStateGas;


  // NIST_MAT_F ----------------------------------------------------------------------------
  fDensityEffectDataTable[  9].fName              = "NIST_MAT_F";
  fDensityEffectDataTable[  9].fPlasmaEnergy      =         0.788 * eV;
  fDensityEffectDataTable[  9].fSternheimerFactor =          2.45;
  fDensityEffectDataTable[  9].fParameterC        =       10.9653;
  fDensityEffectDataTable[  9].fParameterFitX0    =        1.8433;
  fDensityEffectDataTable[  9].fParameterFitX1    =        4.4096;
  fDensityEffectDataTable[  9].fParameterFitA     =       0.11083;
  fDensityEffectDataTable[  9].fParameterFitM     =        3.2962;
  fDensityEffectDataTable[  9].fParameterDelta0   =             0;
  fDensityEffectDataTable[  9].fDeltaErrorMax     =         0.121;
  fDensityEffectDataTable[  9].fState             = MaterialState::kStateGas;


  // NIST_MAT_Ne ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 10].fName              = "NIST_MAT_Ne";
  fDensityEffectDataTable[ 10].fPlasmaEnergy      =         0.587 * eV;
  fDensityEffectDataTable[ 10].fSternheimerFactor =         2.577;
  fDensityEffectDataTable[ 10].fParameterC        =       11.9041;
  fDensityEffectDataTable[ 10].fParameterFitX0    =        2.0735;
  fDensityEffectDataTable[ 10].fParameterFitX1    =        4.6421;
  fDensityEffectDataTable[ 10].fParameterFitA     =       0.08064;
  fDensityEffectDataTable[ 10].fParameterFitM     =        3.5771;
  fDensityEffectDataTable[ 10].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 10].fDeltaErrorMax     =          0.11;
  fDensityEffectDataTable[ 10].fState             = MaterialState::kStateGas;


  // NIST_MAT_Na ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 11].fName              = "NIST_MAT_Na";
  fDensityEffectDataTable[ 11].fPlasmaEnergy      =        19.641 * eV;
  fDensityEffectDataTable[ 11].fSternheimerFactor =         2.648;
  fDensityEffectDataTable[ 11].fParameterC        =        5.0526;
  fDensityEffectDataTable[ 11].fParameterFitX0    =         0.288;
  fDensityEffectDataTable[ 11].fParameterFitX1    =        3.1962;
  fDensityEffectDataTable[ 11].fParameterFitA     =       0.07772;
  fDensityEffectDataTable[ 11].fParameterFitM     =        3.6452;
  fDensityEffectDataTable[ 11].fParameterDelta0   =          0.08;
  fDensityEffectDataTable[ 11].fDeltaErrorMax     =         0.098;
  fDensityEffectDataTable[ 11].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Mg ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 12].fName              = "NIST_MAT_Mg";
  fDensityEffectDataTable[ 12].fPlasmaEnergy      =        26.708 * eV;
  fDensityEffectDataTable[ 12].fSternheimerFactor =         2.331;
  fDensityEffectDataTable[ 12].fParameterC        =        4.5297;
  fDensityEffectDataTable[ 12].fParameterFitX0    =        0.1499;
  fDensityEffectDataTable[ 12].fParameterFitX1    =        3.0668;
  fDensityEffectDataTable[ 12].fParameterFitA     =      0.081163;
  fDensityEffectDataTable[ 12].fParameterFitM     =        3.6166;
  fDensityEffectDataTable[ 12].fParameterDelta0   =          0.08;
  fDensityEffectDataTable[ 12].fDeltaErrorMax     =         0.073;
  fDensityEffectDataTable[ 12].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Al ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 13].fName              = "NIST_MAT_Al";
  fDensityEffectDataTable[ 13].fPlasmaEnergy      =         32.86 * eV;
  fDensityEffectDataTable[ 13].fSternheimerFactor =          2.18;
  fDensityEffectDataTable[ 13].fParameterC        =        4.2395;
  fDensityEffectDataTable[ 13].fParameterFitX0    =        0.1708;
  fDensityEffectDataTable[ 13].fParameterFitX1    =        3.0127;
  fDensityEffectDataTable[ 13].fParameterFitA     =       0.08024;
  fDensityEffectDataTable[ 13].fParameterFitM     =        3.6345;
  fDensityEffectDataTable[ 13].fParameterDelta0   =          0.12;
  fDensityEffectDataTable[ 13].fDeltaErrorMax     =         0.061;
  fDensityEffectDataTable[ 13].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Si ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 14].fName              = "NIST_MAT_Si";
  fDensityEffectDataTable[ 14].fPlasmaEnergy      =        31.055 * eV;
  fDensityEffectDataTable[ 14].fSternheimerFactor =         2.103;
  fDensityEffectDataTable[ 14].fParameterC        =        4.4351;
  fDensityEffectDataTable[ 14].fParameterFitX0    =        0.2014;
  fDensityEffectDataTable[ 14].fParameterFitX1    =        2.8715;
  fDensityEffectDataTable[ 14].fParameterFitA     =       0.14921;
  fDensityEffectDataTable[ 14].fParameterFitM     =        3.2546;
  fDensityEffectDataTable[ 14].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 14].fDeltaErrorMax     =         0.059;
  fDensityEffectDataTable[ 14].fState             = MaterialState::kStateSolid;


  // NIST_MAT_P ----------------------------------------------------------------------------
  fDensityEffectDataTable[ 15].fName              = "NIST_MAT_P";
  fDensityEffectDataTable[ 15].fPlasmaEnergy      =        29.743 * eV;
  fDensityEffectDataTable[ 15].fSternheimerFactor =         2.056;
  fDensityEffectDataTable[ 15].fParameterC        =        4.5214;
  fDensityEffectDataTable[ 15].fParameterFitX0    =        0.1696;
  fDensityEffectDataTable[ 15].fParameterFitX1    =        2.7815;
  fDensityEffectDataTable[ 15].fParameterFitA     =        0.2361;
  fDensityEffectDataTable[ 15].fParameterFitM     =        2.9158;
  fDensityEffectDataTable[ 15].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 15].fDeltaErrorMax     =         0.057;
  fDensityEffectDataTable[ 15].fState             = MaterialState::kStateSolid;


  // NIST_MAT_S ----------------------------------------------------------------------------
  fDensityEffectDataTable[ 16].fName              = "NIST_MAT_S";
  fDensityEffectDataTable[ 16].fPlasmaEnergy      =        28.789 * eV;
  fDensityEffectDataTable[ 16].fSternheimerFactor =         2.131;
  fDensityEffectDataTable[ 16].fParameterC        =        4.6659;
  fDensityEffectDataTable[ 16].fParameterFitX0    =         0.158;
  fDensityEffectDataTable[ 16].fParameterFitX1    =        2.7159;
  fDensityEffectDataTable[ 16].fParameterFitA     =       0.33992;
  fDensityEffectDataTable[ 16].fParameterFitM     =        2.6456;
  fDensityEffectDataTable[ 16].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 16].fDeltaErrorMax     =         0.059;
  fDensityEffectDataTable[ 16].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Cl ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 17].fName              = "NIST_MAT_Cl";
  fDensityEffectDataTable[ 17].fPlasmaEnergy      =         1.092 * eV;
  fDensityEffectDataTable[ 17].fSternheimerFactor =         1.734;
  fDensityEffectDataTable[ 17].fParameterC        =       11.1421;
  fDensityEffectDataTable[ 17].fParameterFitX0    =        1.5555;
  fDensityEffectDataTable[ 17].fParameterFitX1    =        4.2994;
  fDensityEffectDataTable[ 17].fParameterFitA     =       0.19849;
  fDensityEffectDataTable[ 17].fParameterFitM     =        2.9702;
  fDensityEffectDataTable[ 17].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 17].fDeltaErrorMax     =         0.041;
  fDensityEffectDataTable[ 17].fState             = MaterialState::kStateGas;


  // NIST_MAT_Ar ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 18].fName              = "NIST_MAT_Ar";
  fDensityEffectDataTable[ 18].fPlasmaEnergy      =         0.789 * eV;
  fDensityEffectDataTable[ 18].fSternheimerFactor =         1.753;
  fDensityEffectDataTable[ 18].fParameterC        =        11.948;
  fDensityEffectDataTable[ 18].fParameterFitX0    =        1.7635;
  fDensityEffectDataTable[ 18].fParameterFitX1    =        4.4855;
  fDensityEffectDataTable[ 18].fParameterFitA     =       0.19714;
  fDensityEffectDataTable[ 18].fParameterFitM     =        2.9618;
  fDensityEffectDataTable[ 18].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 18].fDeltaErrorMax     =         0.037;
  fDensityEffectDataTable[ 18].fState             = MaterialState::kStateGas;


  // NIST_MAT_K ----------------------------------------------------------------------------
  fDensityEffectDataTable[ 19].fName              = "NIST_MAT_K";
  fDensityEffectDataTable[ 19].fPlasmaEnergy      =         18.65 * eV;
  fDensityEffectDataTable[ 19].fSternheimerFactor =          1.83;
  fDensityEffectDataTable[ 19].fParameterC        =        5.6423;
  fDensityEffectDataTable[ 19].fParameterFitX0    =        0.3851;
  fDensityEffectDataTable[ 19].fParameterFitX1    =        3.1724;
  fDensityEffectDataTable[ 19].fParameterFitA     =       0.19827;
  fDensityEffectDataTable[ 19].fParameterFitM     =        2.9233;
  fDensityEffectDataTable[ 19].fParameterDelta0   =           0.1;
  fDensityEffectDataTable[ 19].fDeltaErrorMax     =         0.035;
  fDensityEffectDataTable[ 19].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ca ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 20].fName              = "NIST_MAT_Ca";
  fDensityEffectDataTable[ 20].fPlasmaEnergy      =        25.342 * eV;
  fDensityEffectDataTable[ 20].fSternheimerFactor =         1.666;
  fDensityEffectDataTable[ 20].fParameterC        =        5.0396;
  fDensityEffectDataTable[ 20].fParameterFitX0    =        0.3228;
  fDensityEffectDataTable[ 20].fParameterFitX1    =        3.1191;
  fDensityEffectDataTable[ 20].fParameterFitA     =       0.15643;
  fDensityEffectDataTable[ 20].fParameterFitM     =        3.0745;
  fDensityEffectDataTable[ 20].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 20].fDeltaErrorMax     =         0.031;
  fDensityEffectDataTable[ 20].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Sc ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 21].fName              = "NIST_MAT_Sc";
  fDensityEffectDataTable[ 21].fPlasmaEnergy      =         34.05 * eV;
  fDensityEffectDataTable[ 21].fSternheimerFactor =         1.826;
  fDensityEffectDataTable[ 21].fParameterC        =        4.6949;
  fDensityEffectDataTable[ 21].fParameterFitX0    =         0.164;
  fDensityEffectDataTable[ 21].fParameterFitX1    =        3.0593;
  fDensityEffectDataTable[ 21].fParameterFitA     =       0.15754;
  fDensityEffectDataTable[ 21].fParameterFitM     =        3.0517;
  fDensityEffectDataTable[ 21].fParameterDelta0   =           0.1;
  fDensityEffectDataTable[ 21].fDeltaErrorMax     =         0.027;
  fDensityEffectDataTable[ 21].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ti ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 22].fName              = "NIST_MAT_Ti";
  fDensityEffectDataTable[ 22].fPlasmaEnergy      =        41.619 * eV;
  fDensityEffectDataTable[ 22].fSternheimerFactor =         1.969;
  fDensityEffectDataTable[ 22].fParameterC        =         4.445;
  fDensityEffectDataTable[ 22].fParameterFitX0    =        0.0957;
  fDensityEffectDataTable[ 22].fParameterFitX1    =        3.0386;
  fDensityEffectDataTable[ 22].fParameterFitA     =       0.15662;
  fDensityEffectDataTable[ 22].fParameterFitM     =        3.0302;
  fDensityEffectDataTable[ 22].fParameterDelta0   =          0.12;
  fDensityEffectDataTable[ 22].fDeltaErrorMax     =         0.025;
  fDensityEffectDataTable[ 22].fState             = MaterialState::kStateSolid;


  // NIST_MAT_V ----------------------------------------------------------------------------
  fDensityEffectDataTable[ 23].fName              = "NIST_MAT_V";
  fDensityEffectDataTable[ 23].fPlasmaEnergy      =        47.861 * eV;
  fDensityEffectDataTable[ 23].fSternheimerFactor =          2.07;
  fDensityEffectDataTable[ 23].fParameterC        =        4.2659;
  fDensityEffectDataTable[ 23].fParameterFitX0    =        0.0691;
  fDensityEffectDataTable[ 23].fParameterFitX1    =        3.0322;
  fDensityEffectDataTable[ 23].fParameterFitA     =       0.15436;
  fDensityEffectDataTable[ 23].fParameterFitM     =        3.0163;
  fDensityEffectDataTable[ 23].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 23].fDeltaErrorMax     =         0.024;
  fDensityEffectDataTable[ 23].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Cr ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 24].fName              = "NIST_MAT_Cr";
  fDensityEffectDataTable[ 24].fPlasmaEnergy      =        52.458 * eV;
  fDensityEffectDataTable[ 24].fSternheimerFactor =         2.181;
  fDensityEffectDataTable[ 24].fParameterC        =        4.1781;
  fDensityEffectDataTable[ 24].fParameterFitX0    =         0.034;
  fDensityEffectDataTable[ 24].fParameterFitX1    =        3.0451;
  fDensityEffectDataTable[ 24].fParameterFitA     =       0.15419;
  fDensityEffectDataTable[ 24].fParameterFitM     =        2.9896;
  fDensityEffectDataTable[ 24].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 24].fDeltaErrorMax     =         0.023;
  fDensityEffectDataTable[ 24].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Mn ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 25].fName              = "NIST_MAT_Mn";
  fDensityEffectDataTable[ 25].fPlasmaEnergy      =        53.022 * eV;
  fDensityEffectDataTable[ 25].fSternheimerFactor =         2.347;
  fDensityEffectDataTable[ 25].fParameterC        =        4.2702;
  fDensityEffectDataTable[ 25].fParameterFitX0    =        0.0447;
  fDensityEffectDataTable[ 25].fParameterFitX1    =        3.1074;
  fDensityEffectDataTable[ 25].fParameterFitA     =       0.14973;
  fDensityEffectDataTable[ 25].fParameterFitM     =        2.9796;
  fDensityEffectDataTable[ 25].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 25].fDeltaErrorMax     =         0.021;
  fDensityEffectDataTable[ 25].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Fe ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 26].fName              = "NIST_MAT_Fe";
  fDensityEffectDataTable[ 26].fPlasmaEnergy      =        55.172 * eV;
  fDensityEffectDataTable[ 26].fSternheimerFactor =         2.504;
  fDensityEffectDataTable[ 26].fParameterC        =        4.2911;
  fDensityEffectDataTable[ 26].fParameterFitX0    =       -0.0012;
  fDensityEffectDataTable[ 26].fParameterFitX1    =        3.1531;
  fDensityEffectDataTable[ 26].fParameterFitA     =         0.146;
  fDensityEffectDataTable[ 26].fParameterFitM     =        2.9632;
  fDensityEffectDataTable[ 26].fParameterDelta0   =          0.12;
  fDensityEffectDataTable[ 26].fDeltaErrorMax     =         0.021;
  fDensityEffectDataTable[ 26].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Co ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 27].fName              = "NIST_MAT_Co";
  fDensityEffectDataTable[ 27].fPlasmaEnergy      =        58.188 * eV;
  fDensityEffectDataTable[ 27].fSternheimerFactor =         2.626;
  fDensityEffectDataTable[ 27].fParameterC        =        4.2601;
  fDensityEffectDataTable[ 27].fParameterFitX0    =       -0.0187;
  fDensityEffectDataTable[ 27].fParameterFitX1    =         3.179;
  fDensityEffectDataTable[ 27].fParameterFitA     =       0.14474;
  fDensityEffectDataTable[ 27].fParameterFitM     =        2.9502;
  fDensityEffectDataTable[ 27].fParameterDelta0   =          0.12;
  fDensityEffectDataTable[ 27].fDeltaErrorMax     =         0.019;
  fDensityEffectDataTable[ 27].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ni ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 28].fName              = "NIST_MAT_Ni";
  fDensityEffectDataTable[ 28].fPlasmaEnergy      =        59.385 * eV;
  fDensityEffectDataTable[ 28].fSternheimerFactor =         2.889;
  fDensityEffectDataTable[ 28].fParameterC        =        4.3115;
  fDensityEffectDataTable[ 28].fParameterFitX0    =       -0.0566;
  fDensityEffectDataTable[ 28].fParameterFitX1    =        3.1851;
  fDensityEffectDataTable[ 28].fParameterFitA     =       0.16496;
  fDensityEffectDataTable[ 28].fParameterFitM     =         2.843;
  fDensityEffectDataTable[ 28].fParameterDelta0   =           0.1;
  fDensityEffectDataTable[ 28].fDeltaErrorMax     =          0.02;
  fDensityEffectDataTable[ 28].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Cu ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 29].fName              = "NIST_MAT_Cu";
  fDensityEffectDataTable[ 29].fPlasmaEnergy      =         58.27 * eV;
  fDensityEffectDataTable[ 29].fSternheimerFactor =         2.956;
  fDensityEffectDataTable[ 29].fParameterC        =         4.419;
  fDensityEffectDataTable[ 29].fParameterFitX0    =       -0.0254;
  fDensityEffectDataTable[ 29].fParameterFitX1    =        3.2792;
  fDensityEffectDataTable[ 29].fParameterFitA     =       0.14339;
  fDensityEffectDataTable[ 29].fParameterFitM     =        2.9044;
  fDensityEffectDataTable[ 29].fParameterDelta0   =          0.08;
  fDensityEffectDataTable[ 29].fDeltaErrorMax     =         0.019;
  fDensityEffectDataTable[ 29].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Zn ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 30].fName              = "NIST_MAT_Zn";
  fDensityEffectDataTable[ 30].fPlasmaEnergy      =        52.132 * eV;
  fDensityEffectDataTable[ 30].fSternheimerFactor =         3.142;
  fDensityEffectDataTable[ 30].fParameterC        =        4.6906;
  fDensityEffectDataTable[ 30].fParameterFitX0    =        0.0049;
  fDensityEffectDataTable[ 30].fParameterFitX1    =        3.3668;
  fDensityEffectDataTable[ 30].fParameterFitA     =       0.14714;
  fDensityEffectDataTable[ 30].fParameterFitM     =        2.8652;
  fDensityEffectDataTable[ 30].fParameterDelta0   =          0.08;
  fDensityEffectDataTable[ 30].fDeltaErrorMax     =         0.019;
  fDensityEffectDataTable[ 30].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ga ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 31].fName              = "NIST_MAT_Ga";
  fDensityEffectDataTable[ 31].fPlasmaEnergy      =        46.688 * eV;
  fDensityEffectDataTable[ 31].fSternheimerFactor =         2.747;
  fDensityEffectDataTable[ 31].fParameterC        =        4.9353;
  fDensityEffectDataTable[ 31].fParameterFitX0    =        0.2267;
  fDensityEffectDataTable[ 31].fParameterFitX1    =        3.5434;
  fDensityEffectDataTable[ 31].fParameterFitA     =        0.0944;
  fDensityEffectDataTable[ 31].fParameterFitM     =        3.1314;
  fDensityEffectDataTable[ 31].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 31].fDeltaErrorMax     =         0.019;
  fDensityEffectDataTable[ 31].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ge ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 32].fName              = "NIST_MAT_Ge";
  fDensityEffectDataTable[ 32].fPlasmaEnergy      =        44.141 * eV;
  fDensityEffectDataTable[ 32].fSternheimerFactor =         2.461;
  fDensityEffectDataTable[ 32].fParameterC        =        5.1411;
  fDensityEffectDataTable[ 32].fParameterFitX0    =        0.3376;
  fDensityEffectDataTable[ 32].fParameterFitX1    =        3.6096;
  fDensityEffectDataTable[ 32].fParameterFitA     =       0.07188;
  fDensityEffectDataTable[ 32].fParameterFitM     =        3.3306;
  fDensityEffectDataTable[ 32].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 32].fDeltaErrorMax     =         0.025;
  fDensityEffectDataTable[ 32].fState             = MaterialState::kStateSolid;


  // NIST_MAT_As ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 33].fName              = "NIST_MAT_As";
  fDensityEffectDataTable[ 33].fPlasmaEnergy      =        45.779 * eV;
  fDensityEffectDataTable[ 33].fSternheimerFactor =         2.219;
  fDensityEffectDataTable[ 33].fParameterC        =         5.051;
  fDensityEffectDataTable[ 33].fParameterFitX0    =        0.1767;
  fDensityEffectDataTable[ 33].fParameterFitX1    =        3.5702;
  fDensityEffectDataTable[ 33].fParameterFitA     =       0.06633;
  fDensityEffectDataTable[ 33].fParameterFitM     =        3.4176;
  fDensityEffectDataTable[ 33].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 33].fDeltaErrorMax     =          0.03;
  fDensityEffectDataTable[ 33].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Se ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 34].fName              = "NIST_MAT_Se";
  fDensityEffectDataTable[ 34].fPlasmaEnergy      =        40.112 * eV;
  fDensityEffectDataTable[ 34].fSternheimerFactor =         2.104;
  fDensityEffectDataTable[ 34].fParameterC        =         5.321;
  fDensityEffectDataTable[ 34].fParameterFitX0    =        0.2258;
  fDensityEffectDataTable[ 34].fParameterFitX1    =        3.6264;
  fDensityEffectDataTable[ 34].fParameterFitA     =       0.06568;
  fDensityEffectDataTable[ 34].fParameterFitM     =        3.4317;
  fDensityEffectDataTable[ 34].fParameterDelta0   =           0.1;
  fDensityEffectDataTable[ 34].fDeltaErrorMax     =         0.024;
  fDensityEffectDataTable[ 34].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Br ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 35].fName              = "NIST_MAT_Br";
  fDensityEffectDataTable[ 35].fPlasmaEnergy      =         1.604 * eV;
  fDensityEffectDataTable[ 35].fSternheimerFactor =         1.845;
  fDensityEffectDataTable[ 35].fParameterC        =       11.7307;
  fDensityEffectDataTable[ 35].fParameterFitX0    =        1.5262;
  fDensityEffectDataTable[ 35].fParameterFitX1    =        4.9899;
  fDensityEffectDataTable[ 35].fParameterFitA     =       0.06335;
  fDensityEffectDataTable[ 35].fParameterFitM     =         3.467;
  fDensityEffectDataTable[ 35].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 35].fDeltaErrorMax     =         0.022;
  fDensityEffectDataTable[ 35].fState             = MaterialState::kStateGas;


  // NIST_MAT_Kr ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 36].fName              = "NIST_MAT_Kr";
  fDensityEffectDataTable[ 36].fPlasmaEnergy      =         1.114 * eV;
  fDensityEffectDataTable[ 36].fSternheimerFactor =          1.77;
  fDensityEffectDataTable[ 36].fParameterC        =       12.5115;
  fDensityEffectDataTable[ 36].fParameterFitX0    =        1.7158;
  fDensityEffectDataTable[ 36].fParameterFitX1    =        5.0748;
  fDensityEffectDataTable[ 36].fParameterFitA     =       0.07446;
  fDensityEffectDataTable[ 36].fParameterFitM     =        3.4051;
  fDensityEffectDataTable[ 36].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 36].fDeltaErrorMax     =         0.025;
  fDensityEffectDataTable[ 36].fState             = MaterialState::kStateGas;


  // NIST_MAT_Ru ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 37].fName              = "NIST_MAT_Ru";
  fDensityEffectDataTable[ 37].fPlasmaEnergy      =        23.467 * eV;
  fDensityEffectDataTable[ 37].fSternheimerFactor =         1.823;
  fDensityEffectDataTable[ 37].fParameterC        =        6.4776;
  fDensityEffectDataTable[ 37].fParameterFitX0    =        0.5737;
  fDensityEffectDataTable[ 37].fParameterFitX1    =        3.7995;
  fDensityEffectDataTable[ 37].fParameterFitA     =       0.07261;
  fDensityEffectDataTable[ 37].fParameterFitM     =        3.4177;
  fDensityEffectDataTable[ 37].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 37].fDeltaErrorMax     =         0.026;
  fDensityEffectDataTable[ 37].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Sr ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 38].fName              = "NIST_MAT_Sr";
  fDensityEffectDataTable[ 38].fPlasmaEnergy      =        30.244 * eV;
  fDensityEffectDataTable[ 38].fSternheimerFactor =         1.707;
  fDensityEffectDataTable[ 38].fParameterC        =        5.9867;
  fDensityEffectDataTable[ 38].fParameterFitX0    =        0.4585;
  fDensityEffectDataTable[ 38].fParameterFitX1    =        3.6778;
  fDensityEffectDataTable[ 38].fParameterFitA     =       0.07165;
  fDensityEffectDataTable[ 38].fParameterFitM     =        3.4435;
  fDensityEffectDataTable[ 38].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 38].fDeltaErrorMax     =         0.026;
  fDensityEffectDataTable[ 38].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Y ----------------------------------------------------------------------------
  fDensityEffectDataTable[ 39].fName              = "NIST_MAT_Y";
  fDensityEffectDataTable[ 39].fPlasmaEnergy      =        40.346 * eV;
  fDensityEffectDataTable[ 39].fSternheimerFactor =         1.649;
  fDensityEffectDataTable[ 39].fParameterC        =        5.4801;
  fDensityEffectDataTable[ 39].fParameterFitX0    =        0.3608;
  fDensityEffectDataTable[ 39].fParameterFitX1    =        3.5542;
  fDensityEffectDataTable[ 39].fParameterFitA     =       0.07138;
  fDensityEffectDataTable[ 39].fParameterFitM     =        3.4565;
  fDensityEffectDataTable[ 39].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 39].fDeltaErrorMax     =         0.027;
  fDensityEffectDataTable[ 39].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Zr ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 40].fName              = "NIST_MAT_Zr";
  fDensityEffectDataTable[ 40].fPlasmaEnergy      =        48.671 * eV;
  fDensityEffectDataTable[ 40].fSternheimerFactor =         1.638;
  fDensityEffectDataTable[ 40].fParameterC        =        5.1774;
  fDensityEffectDataTable[ 40].fParameterFitX0    =        0.2957;
  fDensityEffectDataTable[ 40].fParameterFitX1    =         3.489;
  fDensityEffectDataTable[ 40].fParameterFitA     =       0.07177;
  fDensityEffectDataTable[ 40].fParameterFitM     =        3.4533;
  fDensityEffectDataTable[ 40].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 40].fDeltaErrorMax     =         0.028;
  fDensityEffectDataTable[ 40].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Nb ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 41].fName              = "NIST_MAT_Nb";
  fDensityEffectDataTable[ 41].fPlasmaEnergy      =        56.039 * eV;
  fDensityEffectDataTable[ 41].fSternheimerFactor =         1.734;
  fDensityEffectDataTable[ 41].fParameterC        =        5.0141;
  fDensityEffectDataTable[ 41].fParameterFitX0    =        0.1785;
  fDensityEffectDataTable[ 41].fParameterFitX1    =        3.2201;
  fDensityEffectDataTable[ 41].fParameterFitA     =       0.13883;
  fDensityEffectDataTable[ 41].fParameterFitM     =         3.093;
  fDensityEffectDataTable[ 41].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 41].fDeltaErrorMax     =         0.036;
  fDensityEffectDataTable[ 41].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Mo ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 42].fName              = "NIST_MAT_Mo";
  fDensityEffectDataTable[ 42].fPlasmaEnergy      =        60.951 * eV;
  fDensityEffectDataTable[ 42].fSternheimerFactor =         1.658;
  fDensityEffectDataTable[ 42].fParameterC        =        4.8793;
  fDensityEffectDataTable[ 42].fParameterFitX0    =        0.2267;
  fDensityEffectDataTable[ 42].fParameterFitX1    =        3.2784;
  fDensityEffectDataTable[ 42].fParameterFitA     =       0.10525;
  fDensityEffectDataTable[ 42].fParameterFitM     =        3.2549;
  fDensityEffectDataTable[ 42].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 42].fDeltaErrorMax     =          0.03;
  fDensityEffectDataTable[ 42].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Tc ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 43].fName              = "NIST_MAT_Tc";
  fDensityEffectDataTable[ 43].fPlasmaEnergy      =         64.76 * eV;
  fDensityEffectDataTable[ 43].fSternheimerFactor =         1.727;
  fDensityEffectDataTable[ 43].fParameterC        =        4.7769;
  fDensityEffectDataTable[ 43].fParameterFitX0    =        0.0949;
  fDensityEffectDataTable[ 43].fParameterFitX1    =        3.1253;
  fDensityEffectDataTable[ 43].fParameterFitA     =       0.16572;
  fDensityEffectDataTable[ 43].fParameterFitM     =        2.9738;
  fDensityEffectDataTable[ 43].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 43].fDeltaErrorMax     =          0.04;
  fDensityEffectDataTable[ 43].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ru ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 44].fName              = "NIST_MAT_Ru";
  fDensityEffectDataTable[ 44].fPlasmaEnergy      =        66.978 * eV;
  fDensityEffectDataTable[ 44].fSternheimerFactor =          1.78;
  fDensityEffectDataTable[ 44].fParameterC        =        4.7694;
  fDensityEffectDataTable[ 44].fParameterFitX0    =        0.0599;
  fDensityEffectDataTable[ 44].fParameterFitX1    =        3.0834;
  fDensityEffectDataTable[ 44].fParameterFitA     =       0.19342;
  fDensityEffectDataTable[ 44].fParameterFitM     =        2.8707;
  fDensityEffectDataTable[ 44].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 44].fDeltaErrorMax     =         0.046;
  fDensityEffectDataTable[ 44].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Rh ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 45].fName              = "NIST_MAT_Rh";
  fDensityEffectDataTable[ 45].fPlasmaEnergy      =        67.128 * eV;
  fDensityEffectDataTable[ 45].fSternheimerFactor =         1.804;
  fDensityEffectDataTable[ 45].fParameterC        =        4.8008;
  fDensityEffectDataTable[ 45].fParameterFitX0    =        0.0576;
  fDensityEffectDataTable[ 45].fParameterFitX1    =        3.1069;
  fDensityEffectDataTable[ 45].fParameterFitA     =       0.19205;
  fDensityEffectDataTable[ 45].fParameterFitM     =        2.8633;
  fDensityEffectDataTable[ 45].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 45].fDeltaErrorMax     =         0.046;
  fDensityEffectDataTable[ 45].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Pd ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 46].fName              = "NIST_MAT_Pd";
  fDensityEffectDataTable[ 46].fPlasmaEnergy      =        65.683 * eV;
  fDensityEffectDataTable[ 46].fSternheimerFactor =         1.911;
  fDensityEffectDataTable[ 46].fParameterC        =        4.9358;
  fDensityEffectDataTable[ 46].fParameterFitX0    =        0.0563;
  fDensityEffectDataTable[ 46].fParameterFitX1    =        3.0555;
  fDensityEffectDataTable[ 46].fParameterFitA     =       0.24178;
  fDensityEffectDataTable[ 46].fParameterFitM     =        2.7239;
  fDensityEffectDataTable[ 46].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 46].fDeltaErrorMax     =         0.047;
  fDensityEffectDataTable[ 46].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ag ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 47].fName              = "NIST_MAT_Ag";
  fDensityEffectDataTable[ 47].fPlasmaEnergy      =        61.635 * eV;
  fDensityEffectDataTable[ 47].fSternheimerFactor =         1.933;
  fDensityEffectDataTable[ 47].fParameterC        =         5.063;
  fDensityEffectDataTable[ 47].fParameterFitX0    =        0.0657;
  fDensityEffectDataTable[ 47].fParameterFitX1    =        3.1074;
  fDensityEffectDataTable[ 47].fParameterFitA     =       0.24585;
  fDensityEffectDataTable[ 47].fParameterFitM     =        2.6899;
  fDensityEffectDataTable[ 47].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 47].fDeltaErrorMax     =         0.052;
  fDensityEffectDataTable[ 47].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Cd ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 48].fName              = "NIST_MAT_Cd";
  fDensityEffectDataTable[ 48].fPlasmaEnergy      =        55.381 * eV;
  fDensityEffectDataTable[ 48].fSternheimerFactor =         1.895;
  fDensityEffectDataTable[ 48].fParameterC        =        5.2727;
  fDensityEffectDataTable[ 48].fParameterFitX0    =        0.1281;
  fDensityEffectDataTable[ 48].fParameterFitX1    =        3.1667;
  fDensityEffectDataTable[ 48].fParameterFitA     =       0.24609;
  fDensityEffectDataTable[ 48].fParameterFitM     =        2.6772;
  fDensityEffectDataTable[ 48].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 48].fDeltaErrorMax     =         0.051;
  fDensityEffectDataTable[ 48].fState             = MaterialState::kStateSolid;


  // NIST_MAT_In ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 49].fName              = "NIST_MAT_In";
  fDensityEffectDataTable[ 49].fPlasmaEnergy      =        50.896 * eV;
  fDensityEffectDataTable[ 49].fSternheimerFactor =         1.851;
  fDensityEffectDataTable[ 49].fParameterC        =        5.5211;
  fDensityEffectDataTable[ 49].fParameterFitX0    =        0.2406;
  fDensityEffectDataTable[ 49].fParameterFitX1    =        3.2032;
  fDensityEffectDataTable[ 49].fParameterFitA     =       0.23879;
  fDensityEffectDataTable[ 49].fParameterFitM     =        2.7144;
  fDensityEffectDataTable[ 49].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 49].fDeltaErrorMax     =         0.044;
  fDensityEffectDataTable[ 49].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Sn ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 50].fName              = "NIST_MAT_Sn";
  fDensityEffectDataTable[ 50].fPlasmaEnergy      =        50.567 * eV;
  fDensityEffectDataTable[ 50].fSternheimerFactor =         1.732;
  fDensityEffectDataTable[ 50].fParameterC        =         5.534;
  fDensityEffectDataTable[ 50].fParameterFitX0    =        0.2879;
  fDensityEffectDataTable[ 50].fParameterFitX1    =        3.2959;
  fDensityEffectDataTable[ 50].fParameterFitA     =       0.18689;
  fDensityEffectDataTable[ 50].fParameterFitM     =        2.8576;
  fDensityEffectDataTable[ 50].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 50].fDeltaErrorMax     =         0.037;
  fDensityEffectDataTable[ 50].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Sb ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 51].fName              = "NIST_MAT_Sb";
  fDensityEffectDataTable[ 51].fPlasmaEnergy      =        48.242 * eV;
  fDensityEffectDataTable[ 51].fSternheimerFactor =         1.645;
  fDensityEffectDataTable[ 51].fParameterC        =        5.6241;
  fDensityEffectDataTable[ 51].fParameterFitX0    =        0.3189;
  fDensityEffectDataTable[ 51].fParameterFitX1    =        3.3489;
  fDensityEffectDataTable[ 51].fParameterFitA     =       0.16652;
  fDensityEffectDataTable[ 51].fParameterFitM     =        2.9319;
  fDensityEffectDataTable[ 51].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 51].fDeltaErrorMax     =         0.034;
  fDensityEffectDataTable[ 51].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Te ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 52].fName              = "NIST_MAT_Te";
  fDensityEffectDataTable[ 52].fPlasmaEnergy      =        45.952 * eV;
  fDensityEffectDataTable[ 52].fSternheimerFactor =         1.577;
  fDensityEffectDataTable[ 52].fParameterC        =        5.7131;
  fDensityEffectDataTable[ 52].fParameterFitX0    =        0.3296;
  fDensityEffectDataTable[ 52].fParameterFitX1    =        3.4418;
  fDensityEffectDataTable[ 52].fParameterFitA     =       0.13815;
  fDensityEffectDataTable[ 52].fParameterFitM     =        3.0354;
  fDensityEffectDataTable[ 52].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 52].fDeltaErrorMax     =         0.033;
  fDensityEffectDataTable[ 52].fState             = MaterialState::kStateSolid;


  // NIST_MAT_I ----------------------------------------------------------------------------
  fDensityEffectDataTable[ 53].fName              = "NIST_MAT_I";
  fDensityEffectDataTable[ 53].fPlasmaEnergy      =        41.348 * eV;
  fDensityEffectDataTable[ 53].fSternheimerFactor =         1.498;
  fDensityEffectDataTable[ 53].fParameterC        =        5.9488;
  fDensityEffectDataTable[ 53].fParameterFitX0    =        0.0549;
  fDensityEffectDataTable[ 53].fParameterFitX1    =        3.2596;
  fDensityEffectDataTable[ 53].fParameterFitA     =       0.23766;
  fDensityEffectDataTable[ 53].fParameterFitM     =        2.7276;
  fDensityEffectDataTable[ 53].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 53].fDeltaErrorMax     =         0.045;
  fDensityEffectDataTable[ 53].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Xe ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 54].fName              = "NIST_MAT_Xe";
  fDensityEffectDataTable[ 54].fPlasmaEnergy      =         1.369 * eV;
  fDensityEffectDataTable[ 54].fSternheimerFactor =         1.435;
  fDensityEffectDataTable[ 54].fParameterC        =       12.7281;
  fDensityEffectDataTable[ 54].fParameterFitX0    =         1.563;
  fDensityEffectDataTable[ 54].fParameterFitX1    =        4.7371;
  fDensityEffectDataTable[ 54].fParameterFitA     =       0.23314;
  fDensityEffectDataTable[ 54].fParameterFitM     =        2.7414;
  fDensityEffectDataTable[ 54].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 54].fDeltaErrorMax     =         0.043;
  fDensityEffectDataTable[ 54].fState             = MaterialState::kStateGas;


  // NIST_MAT_Cs ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 55].fName              = "NIST_MAT_Cs";
  fDensityEffectDataTable[ 55].fPlasmaEnergy      =         25.37 * eV;
  fDensityEffectDataTable[ 55].fSternheimerFactor =         1.462;
  fDensityEffectDataTable[ 55].fParameterC        =        6.9135;
  fDensityEffectDataTable[ 55].fParameterFitX0    =        0.5473;
  fDensityEffectDataTable[ 55].fParameterFitX1    =        3.5914;
  fDensityEffectDataTable[ 55].fParameterFitA     =       0.18233;
  fDensityEffectDataTable[ 55].fParameterFitM     =        2.8866;
  fDensityEffectDataTable[ 55].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 55].fDeltaErrorMax     =         0.035;
  fDensityEffectDataTable[ 55].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ba ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 56].fName              = "NIST_MAT_Ba";
  fDensityEffectDataTable[ 56].fPlasmaEnergy      =        34.425 * eV;
  fDensityEffectDataTable[ 56].fSternheimerFactor =          1.41;
  fDensityEffectDataTable[ 56].fParameterC        =        6.3153;
  fDensityEffectDataTable[ 56].fParameterFitX0    =         0.419;
  fDensityEffectDataTable[ 56].fParameterFitX1    =        3.4547;
  fDensityEffectDataTable[ 56].fParameterFitA     =       0.18268;
  fDensityEffectDataTable[ 56].fParameterFitM     =        2.8906;
  fDensityEffectDataTable[ 56].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 56].fDeltaErrorMax     =         0.035;
  fDensityEffectDataTable[ 56].fState             = MaterialState::kStateSolid;


  // NIST_MAT_La ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 57].fName              = "NIST_MAT_La";
  fDensityEffectDataTable[ 57].fPlasmaEnergy      =        45.792 * eV;
  fDensityEffectDataTable[ 57].fSternheimerFactor =         1.392;
  fDensityEffectDataTable[ 57].fParameterC        =         5.785;
  fDensityEffectDataTable[ 57].fParameterFitX0    =        0.3161;
  fDensityEffectDataTable[ 57].fParameterFitX1    =        3.3293;
  fDensityEffectDataTable[ 57].fParameterFitA     =       0.18591;
  fDensityEffectDataTable[ 57].fParameterFitM     =        2.8828;
  fDensityEffectDataTable[ 57].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 57].fDeltaErrorMax     =         0.036;
  fDensityEffectDataTable[ 57].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ce ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 58].fName              = "NIST_MAT_Ce";
  fDensityEffectDataTable[ 58].fPlasmaEnergy      =        47.834 * eV;
  fDensityEffectDataTable[ 58].fSternheimerFactor =         1.461;
  fDensityEffectDataTable[ 58].fParameterC        =        5.7837;
  fDensityEffectDataTable[ 58].fParameterFitX0    =        0.2713;
  fDensityEffectDataTable[ 58].fParameterFitX1    =        3.3432;
  fDensityEffectDataTable[ 58].fParameterFitA     =       0.18885;
  fDensityEffectDataTable[ 58].fParameterFitM     =        2.8592;
  fDensityEffectDataTable[ 58].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 58].fDeltaErrorMax     =          0.04;
  fDensityEffectDataTable[ 58].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Pr ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 59].fName              = "NIST_MAT_Pr";
  fDensityEffectDataTable[ 59].fPlasmaEnergy      =        48.301 * eV;
  fDensityEffectDataTable[ 59].fSternheimerFactor =          1.52;
  fDensityEffectDataTable[ 59].fParameterC        =        5.8096;
  fDensityEffectDataTable[ 59].fParameterFitX0    =        0.2333;
  fDensityEffectDataTable[ 59].fParameterFitX1    =        3.2773;
  fDensityEffectDataTable[ 59].fParameterFitA     =       0.23265;
  fDensityEffectDataTable[ 59].fParameterFitM     =        2.7331;
  fDensityEffectDataTable[ 59].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 59].fDeltaErrorMax     =         0.041;
  fDensityEffectDataTable[ 59].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ne ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 60].fName              = "NIST_MAT_Ne";
  fDensityEffectDataTable[ 60].fPlasmaEnergy      =        48.819 * eV;
  fDensityEffectDataTable[ 60].fSternheimerFactor =         1.588;
  fDensityEffectDataTable[ 60].fParameterC        =         5.829;
  fDensityEffectDataTable[ 60].fParameterFitX0    =        0.1984;
  fDensityEffectDataTable[ 60].fParameterFitX1    =        3.3063;
  fDensityEffectDataTable[ 60].fParameterFitA     =        0.2353;
  fDensityEffectDataTable[ 60].fParameterFitM     =         2.705;
  fDensityEffectDataTable[ 60].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 60].fDeltaErrorMax     =         0.044;
  fDensityEffectDataTable[ 60].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Pr ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 61].fName              = "NIST_MAT_Pr";
  fDensityEffectDataTable[ 61].fPlasmaEnergy      =        50.236 * eV;
  fDensityEffectDataTable[ 61].fSternheimerFactor =         1.672;
  fDensityEffectDataTable[ 61].fParameterC        =        5.8224;
  fDensityEffectDataTable[ 61].fParameterFitX0    =        0.1627;
  fDensityEffectDataTable[ 61].fParameterFitX1    =        3.3199;
  fDensityEffectDataTable[ 61].fParameterFitA     =        0.2428;
  fDensityEffectDataTable[ 61].fParameterFitM     =        2.6674;
  fDensityEffectDataTable[ 61].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 61].fDeltaErrorMax     =         0.048;
  fDensityEffectDataTable[ 61].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Sa ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 62].fName              = "NIST_MAT_Sa";
  fDensityEffectDataTable[ 62].fPlasmaEnergy      =         50.54 * eV;
  fDensityEffectDataTable[ 62].fSternheimerFactor =         1.749;
  fDensityEffectDataTable[ 62].fParameterC        =        5.8597;
  fDensityEffectDataTable[ 62].fParameterFitX0    =         0.152;
  fDensityEffectDataTable[ 62].fParameterFitX1    =         3.346;
  fDensityEffectDataTable[ 62].fParameterFitA     =       0.24698;
  fDensityEffectDataTable[ 62].fParameterFitM     =        2.6403;
  fDensityEffectDataTable[ 62].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 62].fDeltaErrorMax     =         0.053;
  fDensityEffectDataTable[ 62].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Eu ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 63].fName              = "NIST_MAT_Eu";
  fDensityEffectDataTable[ 63].fPlasmaEnergy      =        42.484 * eV;
  fDensityEffectDataTable[ 63].fSternheimerFactor =         1.838;
  fDensityEffectDataTable[ 63].fParameterC        =        6.2278;
  fDensityEffectDataTable[ 63].fParameterFitX0    =        0.1888;
  fDensityEffectDataTable[ 63].fParameterFitX1    =        3.4633;
  fDensityEffectDataTable[ 63].fParameterFitA     =       0.24448;
  fDensityEffectDataTable[ 63].fParameterFitM     =        2.6245;
  fDensityEffectDataTable[ 63].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 63].fDeltaErrorMax     =          0.06;
  fDensityEffectDataTable[ 63].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Gd ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 64].fName              = "NIST_MAT_Gd";
  fDensityEffectDataTable[ 64].fPlasmaEnergy      =        51.672 * eV;
  fDensityEffectDataTable[ 64].fSternheimerFactor =         1.882;
  fDensityEffectDataTable[ 64].fParameterC        =        5.8738;
  fDensityEffectDataTable[ 64].fParameterFitX0    =        0.1058;
  fDensityEffectDataTable[ 64].fParameterFitX1    =        3.3932;
  fDensityEffectDataTable[ 64].fParameterFitA     =       0.25109;
  fDensityEffectDataTable[ 64].fParameterFitM     =        2.5977;
  fDensityEffectDataTable[ 64].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 64].fDeltaErrorMax     =         0.061;
  fDensityEffectDataTable[ 64].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Tb ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 65].fName              = "NIST_MAT_Tb";
  fDensityEffectDataTable[ 65].fPlasmaEnergy      =        52.865 * eV;
  fDensityEffectDataTable[ 65].fSternheimerFactor =         1.993;
  fDensityEffectDataTable[ 65].fParameterC        =        5.9045;
  fDensityEffectDataTable[ 65].fParameterFitX0    =        0.0947;
  fDensityEffectDataTable[ 65].fParameterFitX1    =        3.4224;
  fDensityEffectDataTable[ 65].fParameterFitA     =       0.24453;
  fDensityEffectDataTable[ 65].fParameterFitM     =        2.6056;
  fDensityEffectDataTable[ 65].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 65].fDeltaErrorMax     =         0.063;
  fDensityEffectDataTable[ 65].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Dy ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 66].fName              = "NIST_MAT_Dy";
  fDensityEffectDataTable[ 66].fPlasmaEnergy      =        53.698 * eV;
  fDensityEffectDataTable[ 66].fSternheimerFactor =         2.081;
  fDensityEffectDataTable[ 66].fParameterC        =        5.9183;
  fDensityEffectDataTable[ 66].fParameterFitX0    =        0.0822;
  fDensityEffectDataTable[ 66].fParameterFitX1    =        3.4474;
  fDensityEffectDataTable[ 66].fParameterFitA     =       0.24665;
  fDensityEffectDataTable[ 66].fParameterFitM     =        2.5849;
  fDensityEffectDataTable[ 66].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 66].fDeltaErrorMax     =         0.061;
  fDensityEffectDataTable[ 66].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ho ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 67].fName              = "NIST_MAT_Ho";
  fDensityEffectDataTable[ 67].fPlasmaEnergy      =        54.467 * eV;
  fDensityEffectDataTable[ 67].fSternheimerFactor =         2.197;
  fDensityEffectDataTable[ 67].fParameterC        =        5.9587;
  fDensityEffectDataTable[ 67].fParameterFitX0    =        0.0761;
  fDensityEffectDataTable[ 67].fParameterFitX1    =        3.4782;
  fDensityEffectDataTable[ 67].fParameterFitA     =       0.24638;
  fDensityEffectDataTable[ 67].fParameterFitM     =        2.5726;
  fDensityEffectDataTable[ 67].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 67].fDeltaErrorMax     =         0.062;
  fDensityEffectDataTable[ 67].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Er ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 68].fName              = "NIST_MAT_Er";
  fDensityEffectDataTable[ 68].fPlasmaEnergy      =        55.322 * eV;
  fDensityEffectDataTable[ 68].fSternheimerFactor =          2.26;
  fDensityEffectDataTable[ 68].fParameterC        =        5.9521;
  fDensityEffectDataTable[ 68].fParameterFitX0    =        0.0648;
  fDensityEffectDataTable[ 68].fParameterFitX1    =        3.4922;
  fDensityEffectDataTable[ 68].fParameterFitA     =       0.24823;
  fDensityEffectDataTable[ 68].fParameterFitM     =        2.5573;
  fDensityEffectDataTable[ 68].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 68].fDeltaErrorMax     =         0.061;
  fDensityEffectDataTable[ 68].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Tm ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 69].fName              = "NIST_MAT_Tm";
  fDensityEffectDataTable[ 69].fPlasmaEnergy      =        56.225 * eV;
  fDensityEffectDataTable[ 69].fSternheimerFactor =         2.333;
  fDensityEffectDataTable[ 69].fParameterC        =        5.9677;
  fDensityEffectDataTable[ 69].fParameterFitX0    =        0.0812;
  fDensityEffectDataTable[ 69].fParameterFitX1    =        3.5085;
  fDensityEffectDataTable[ 69].fParameterFitA     =       0.24189;
  fDensityEffectDataTable[ 69].fParameterFitM     =        2.5469;
  fDensityEffectDataTable[ 69].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 69].fDeltaErrorMax     =         0.062;
  fDensityEffectDataTable[ 69].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Yb ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 70].fName              = "NIST_MAT_Yb";
  fDensityEffectDataTable[ 70].fPlasmaEnergy      =        47.546 * eV;
  fDensityEffectDataTable[ 70].fSternheimerFactor =         2.505;
  fDensityEffectDataTable[ 70].fParameterC        =        6.3325;
  fDensityEffectDataTable[ 70].fParameterFitX0    =        0.1199;
  fDensityEffectDataTable[ 70].fParameterFitX1    =        3.6246;
  fDensityEffectDataTable[ 70].fParameterFitA     =       0.25295;
  fDensityEffectDataTable[ 70].fParameterFitM     =        2.5141;
  fDensityEffectDataTable[ 70].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 70].fDeltaErrorMax     =         0.071;
  fDensityEffectDataTable[ 70].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Lu ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 71].fName              = "NIST_MAT_Lu";
  fDensityEffectDataTable[ 71].fPlasmaEnergy      =        57.581 * eV;
  fDensityEffectDataTable[ 71].fSternheimerFactor =         2.348;
  fDensityEffectDataTable[ 71].fParameterC        =        5.9785;
  fDensityEffectDataTable[ 71].fParameterFitX0    =         0.156;
  fDensityEffectDataTable[ 71].fParameterFitX1    =        3.5218;
  fDensityEffectDataTable[ 71].fParameterFitA     =       0.24033;
  fDensityEffectDataTable[ 71].fParameterFitM     =        2.5643;
  fDensityEffectDataTable[ 71].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 71].fDeltaErrorMax     =         0.054;
  fDensityEffectDataTable[ 71].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Hf ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 72].fName              = "NIST_MAT_Hf";
  fDensityEffectDataTable[ 72].fPlasmaEnergy      =         66.77 * eV;
  fDensityEffectDataTable[ 72].fSternheimerFactor =         2.174;
  fDensityEffectDataTable[ 72].fParameterC        =        5.7139;
  fDensityEffectDataTable[ 72].fParameterFitX0    =        0.1965;
  fDensityEffectDataTable[ 72].fParameterFitX1    =        3.4337;
  fDensityEffectDataTable[ 72].fParameterFitA     =       0.22918;
  fDensityEffectDataTable[ 72].fParameterFitM     =        2.6155;
  fDensityEffectDataTable[ 72].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 72].fDeltaErrorMax     =         0.035;
  fDensityEffectDataTable[ 72].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ta ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 73].fName              = "NIST_MAT_Ta";
  fDensityEffectDataTable[ 73].fPlasmaEnergy      =        74.692 * eV;
  fDensityEffectDataTable[ 73].fSternheimerFactor =          2.07;
  fDensityEffectDataTable[ 73].fParameterC        =        5.5262;
  fDensityEffectDataTable[ 73].fParameterFitX0    =        0.2117;
  fDensityEffectDataTable[ 73].fParameterFitX1    =        3.4805;
  fDensityEffectDataTable[ 73].fParameterFitA     =       0.17798;
  fDensityEffectDataTable[ 73].fParameterFitM     =        2.7623;
  fDensityEffectDataTable[ 73].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 73].fDeltaErrorMax     =          0.03;
  fDensityEffectDataTable[ 73].fState             = MaterialState::kStateSolid;


  // NIST_MAT_W ----------------------------------------------------------------------------
  fDensityEffectDataTable[ 74].fName              = "NIST_MAT_W";
  fDensityEffectDataTable[ 74].fPlasmaEnergy      =        80.315 * eV;
  fDensityEffectDataTable[ 74].fSternheimerFactor =         1.997;
  fDensityEffectDataTable[ 74].fParameterC        =        5.4059;
  fDensityEffectDataTable[ 74].fParameterFitX0    =        0.2167;
  fDensityEffectDataTable[ 74].fParameterFitX1    =         3.496;
  fDensityEffectDataTable[ 74].fParameterFitA     =       0.15509;
  fDensityEffectDataTable[ 74].fParameterFitM     =        2.8447;
  fDensityEffectDataTable[ 74].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 74].fDeltaErrorMax     =         0.027;
  fDensityEffectDataTable[ 74].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Re ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 75].fName              = "NIST_MAT_Re";
  fDensityEffectDataTable[ 75].fPlasmaEnergy      =        83.846 * eV;
  fDensityEffectDataTable[ 75].fSternheimerFactor =         1.976;
  fDensityEffectDataTable[ 75].fParameterC        =        5.3445;
  fDensityEffectDataTable[ 75].fParameterFitX0    =        0.0559;
  fDensityEffectDataTable[ 75].fParameterFitX1    =        3.4845;
  fDensityEffectDataTable[ 75].fParameterFitA     =       0.15184;
  fDensityEffectDataTable[ 75].fParameterFitM     =        2.8627;
  fDensityEffectDataTable[ 75].fParameterDelta0   =          0.08;
  fDensityEffectDataTable[ 75].fDeltaErrorMax     =         0.026;
  fDensityEffectDataTable[ 75].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Os ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 76].fName              = "NIST_MAT_Os";
  fDensityEffectDataTable[ 76].fPlasmaEnergy      =        86.537 * eV;
  fDensityEffectDataTable[ 76].fSternheimerFactor =         1.947;
  fDensityEffectDataTable[ 76].fParameterC        =        5.3083;
  fDensityEffectDataTable[ 76].fParameterFitX0    =        0.0891;
  fDensityEffectDataTable[ 76].fParameterFitX1    =        3.5414;
  fDensityEffectDataTable[ 76].fParameterFitA     =       0.12751;
  fDensityEffectDataTable[ 76].fParameterFitM     =        2.9608;
  fDensityEffectDataTable[ 76].fParameterDelta0   =           0.1;
  fDensityEffectDataTable[ 76].fDeltaErrorMax     =         0.023;
  fDensityEffectDataTable[ 76].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ir ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 77].fName              = "NIST_MAT_Ir";
  fDensityEffectDataTable[ 77].fPlasmaEnergy      =        86.357 * eV;
  fDensityEffectDataTable[ 77].fSternheimerFactor =         1.927;
  fDensityEffectDataTable[ 77].fParameterC        =        5.3418;
  fDensityEffectDataTable[ 77].fParameterFitX0    =        0.0819;
  fDensityEffectDataTable[ 77].fParameterFitX1    =         3.548;
  fDensityEffectDataTable[ 77].fParameterFitA     =        0.1269;
  fDensityEffectDataTable[ 77].fParameterFitM     =        2.9658;
  fDensityEffectDataTable[ 77].fParameterDelta0   =           0.1;
  fDensityEffectDataTable[ 77].fDeltaErrorMax     =         0.023;
  fDensityEffectDataTable[ 77].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Pt ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 78].fName              = "NIST_MAT_Pt";
  fDensityEffectDataTable[ 78].fPlasmaEnergy      =        84.389 * eV;
  fDensityEffectDataTable[ 78].fSternheimerFactor =         1.965;
  fDensityEffectDataTable[ 78].fParameterC        =        5.4732;
  fDensityEffectDataTable[ 78].fParameterFitX0    =        0.1484;
  fDensityEffectDataTable[ 78].fParameterFitX1    =        3.6212;
  fDensityEffectDataTable[ 78].fParameterFitA     =       0.11128;
  fDensityEffectDataTable[ 78].fParameterFitM     =        3.0417;
  fDensityEffectDataTable[ 78].fParameterDelta0   =          0.12;
  fDensityEffectDataTable[ 78].fDeltaErrorMax     =         0.021;
  fDensityEffectDataTable[ 78].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Au ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 79].fName              = "NIST_MAT_Au";
  fDensityEffectDataTable[ 79].fPlasmaEnergy      =        80.215 * eV;
  fDensityEffectDataTable[ 79].fSternheimerFactor =         1.926;
  fDensityEffectDataTable[ 79].fParameterC        =        5.5747;
  fDensityEffectDataTable[ 79].fParameterFitX0    =        0.2021;
  fDensityEffectDataTable[ 79].fParameterFitX1    =        3.6979;
  fDensityEffectDataTable[ 79].fParameterFitA     =       0.09756;
  fDensityEffectDataTable[ 79].fParameterFitM     =        3.1101;
  fDensityEffectDataTable[ 79].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 79].fDeltaErrorMax     =          0.02;
  fDensityEffectDataTable[ 79].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Hg ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 80].fName              = "NIST_MAT_Hg";
  fDensityEffectDataTable[ 80].fPlasmaEnergy      =        66.977 * eV;
  fDensityEffectDataTable[ 80].fSternheimerFactor =         1.904;
  fDensityEffectDataTable[ 80].fParameterC        =        5.9605;
  fDensityEffectDataTable[ 80].fParameterFitX0    =        0.2756;
  fDensityEffectDataTable[ 80].fParameterFitX1    =        3.7275;
  fDensityEffectDataTable[ 80].fParameterFitA     =       0.11014;
  fDensityEffectDataTable[ 80].fParameterFitM     =        3.0519;
  fDensityEffectDataTable[ 80].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 80].fDeltaErrorMax     =         0.021;
  fDensityEffectDataTable[ 80].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Tl ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 81].fName              = "NIST_MAT_Tl";
  fDensityEffectDataTable[ 81].fPlasmaEnergy      =        62.104 * eV;
  fDensityEffectDataTable[ 81].fSternheimerFactor =         1.814;
  fDensityEffectDataTable[ 81].fParameterC        =        6.1365;
  fDensityEffectDataTable[ 81].fParameterFitX0    =        0.3491;
  fDensityEffectDataTable[ 81].fParameterFitX1    =        3.8044;
  fDensityEffectDataTable[ 81].fParameterFitA     =       0.09455;
  fDensityEffectDataTable[ 81].fParameterFitM     =         3.145;
  fDensityEffectDataTable[ 81].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 81].fDeltaErrorMax     =         0.019;
  fDensityEffectDataTable[ 81].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Pb ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 82].fName              = "NIST_MAT_Pb";
  fDensityEffectDataTable[ 82].fPlasmaEnergy      =        61.072 * eV;
  fDensityEffectDataTable[ 82].fSternheimerFactor =         1.755;
  fDensityEffectDataTable[ 82].fParameterC        =        6.2018;
  fDensityEffectDataTable[ 82].fParameterFitX0    =        0.3776;
  fDensityEffectDataTable[ 82].fParameterFitX1    =        3.8073;
  fDensityEffectDataTable[ 82].fParameterFitA     =       0.09359;
  fDensityEffectDataTable[ 82].fParameterFitM     =        3.1608;
  fDensityEffectDataTable[ 82].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 82].fDeltaErrorMax     =         0.019;
  fDensityEffectDataTable[ 82].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Bi ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 83].fName              = "NIST_MAT_Bi";
  fDensityEffectDataTable[ 83].fPlasmaEnergy      =        56.696 * eV;
  fDensityEffectDataTable[ 83].fSternheimerFactor =         1.684;
  fDensityEffectDataTable[ 83].fParameterC        =        6.3505;
  fDensityEffectDataTable[ 83].fParameterFitX0    =        0.4152;
  fDensityEffectDataTable[ 83].fParameterFitX1    =        3.8248;
  fDensityEffectDataTable[ 83].fParameterFitA     =        0.0941;
  fDensityEffectDataTable[ 83].fParameterFitM     =        3.1671;
  fDensityEffectDataTable[ 83].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 83].fDeltaErrorMax     =          0.02;
  fDensityEffectDataTable[ 83].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Po ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 84].fName              = "NIST_MAT_Po";
  fDensityEffectDataTable[ 84].fPlasmaEnergy      =        55.773 * eV;
  fDensityEffectDataTable[ 84].fSternheimerFactor =         1.637;
  fDensityEffectDataTable[ 84].fParameterC        =        6.4003;
  fDensityEffectDataTable[ 84].fParameterFitX0    =        0.4267;
  fDensityEffectDataTable[ 84].fParameterFitX1    =        3.8293;
  fDensityEffectDataTable[ 84].fParameterFitA     =       0.09282;
  fDensityEffectDataTable[ 84].fParameterFitM     =         3.183;
  fDensityEffectDataTable[ 84].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 84].fDeltaErrorMax     =          0.02;
  fDensityEffectDataTable[ 84].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Rn ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 85].fName              = "NIST_MAT_Rn";
  fDensityEffectDataTable[ 85].fPlasmaEnergy      =         1.708 * eV;
  fDensityEffectDataTable[ 85].fSternheimerFactor =         1.458;
  fDensityEffectDataTable[ 85].fParameterC        =       13.2839;
  fDensityEffectDataTable[ 85].fParameterFitX0    =        1.5368;
  fDensityEffectDataTable[ 85].fParameterFitX1    =        4.9889;
  fDensityEffectDataTable[ 85].fParameterFitA     =       0.20798;
  fDensityEffectDataTable[ 85].fParameterFitM     =        2.7409;
  fDensityEffectDataTable[ 85].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 85].fDeltaErrorMax     =         0.057;
  fDensityEffectDataTable[ 85].fState             = MaterialState::kStateGas;


  // NIST_MAT_Ra ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 86].fName              = "NIST_MAT_Ra";
  fDensityEffectDataTable[ 86].fPlasmaEnergy      =        40.205 * eV;
  fDensityEffectDataTable[ 86].fSternheimerFactor =         1.403;
  fDensityEffectDataTable[ 86].fParameterC        =        7.0452;
  fDensityEffectDataTable[ 86].fParameterFitX0    =        0.5991;
  fDensityEffectDataTable[ 86].fParameterFitX1    =        3.9428;
  fDensityEffectDataTable[ 86].fParameterFitA     =       0.08804;
  fDensityEffectDataTable[ 86].fParameterFitM     =        3.2454;
  fDensityEffectDataTable[ 86].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 86].fDeltaErrorMax     =         0.022;
  fDensityEffectDataTable[ 86].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Ac ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 87].fName              = "NIST_MAT_Ac";
  fDensityEffectDataTable[ 87].fPlasmaEnergy      =        57.254 * eV;
  fDensityEffectDataTable[ 87].fSternheimerFactor =          1.38;
  fDensityEffectDataTable[ 87].fParameterC        =        6.3742;
  fDensityEffectDataTable[ 87].fParameterFitX0    =        0.4559;
  fDensityEffectDataTable[ 87].fParameterFitX1    =        3.7966;
  fDensityEffectDataTable[ 87].fParameterFitA     =       0.08567;
  fDensityEffectDataTable[ 87].fParameterFitM     =        3.2683;
  fDensityEffectDataTable[ 87].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 87].fDeltaErrorMax     =         0.023;
  fDensityEffectDataTable[ 87].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Th ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 88].fName              = "NIST_MAT_Th";
  fDensityEffectDataTable[ 88].fPlasmaEnergy      =        61.438 * eV;
  fDensityEffectDataTable[ 88].fSternheimerFactor =         1.363;
  fDensityEffectDataTable[ 88].fParameterC        =        6.2473;
  fDensityEffectDataTable[ 88].fParameterFitX0    =        0.4202;
  fDensityEffectDataTable[ 88].fParameterFitX1    =        3.7681;
  fDensityEffectDataTable[ 88].fParameterFitA     =       0.08655;
  fDensityEffectDataTable[ 88].fParameterFitM     =         3.261;
  fDensityEffectDataTable[ 88].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 88].fDeltaErrorMax     =         0.025;
  fDensityEffectDataTable[ 88].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Pa ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 89].fName              = "NIST_MAT_Pa";
  fDensityEffectDataTable[ 89].fPlasmaEnergy      =        70.901 * eV;
  fDensityEffectDataTable[ 89].fSternheimerFactor =          1.42;
  fDensityEffectDataTable[ 89].fParameterC        =        6.0327;
  fDensityEffectDataTable[ 89].fParameterFitX0    =        0.3144;
  fDensityEffectDataTable[ 89].fParameterFitX1    =        3.5079;
  fDensityEffectDataTable[ 89].fParameterFitA     =        0.1477;
  fDensityEffectDataTable[ 89].fParameterFitM     =        2.9845;
  fDensityEffectDataTable[ 89].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 89].fDeltaErrorMax     =         0.036;
  fDensityEffectDataTable[ 89].fState             = MaterialState::kStateSolid;


  // NIST_MAT_U ----------------------------------------------------------------------------
  fDensityEffectDataTable[ 90].fName              = "NIST_MAT_U";
  fDensityEffectDataTable[ 90].fPlasmaEnergy      =        77.986 * eV;
  fDensityEffectDataTable[ 90].fSternheimerFactor =         1.447;
  fDensityEffectDataTable[ 90].fParameterC        =        5.8694;
  fDensityEffectDataTable[ 90].fParameterFitX0    =         0.226;
  fDensityEffectDataTable[ 90].fParameterFitX1    =        3.3721;
  fDensityEffectDataTable[ 90].fParameterFitA     =       0.19677;
  fDensityEffectDataTable[ 90].fParameterFitM     =        2.8171;
  fDensityEffectDataTable[ 90].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 90].fDeltaErrorMax     =         0.043;
  fDensityEffectDataTable[ 90].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Np ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 91].fName              = "NIST_MAT_Np";
  fDensityEffectDataTable[ 91].fPlasmaEnergy      =        81.221 * eV;
  fDensityEffectDataTable[ 91].fSternheimerFactor =         1.468;
  fDensityEffectDataTable[ 91].fParameterC        =        5.8149;
  fDensityEffectDataTable[ 91].fParameterFitX0    =        0.1869;
  fDensityEffectDataTable[ 91].fParameterFitX1    =         3.369;
  fDensityEffectDataTable[ 91].fParameterFitA     =       0.19741;
  fDensityEffectDataTable[ 91].fParameterFitM     =        2.8082;
  fDensityEffectDataTable[ 91].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 91].fDeltaErrorMax     =         0.043;
  fDensityEffectDataTable[ 91].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Pu ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 92].fName              = "NIST_MAT_Pu";
  fDensityEffectDataTable[ 92].fPlasmaEnergy      =        80.486 * eV;
  fDensityEffectDataTable[ 92].fSternheimerFactor =         1.519;
  fDensityEffectDataTable[ 92].fParameterC        =        5.8748;
  fDensityEffectDataTable[ 92].fParameterFitX0    =        0.1557;
  fDensityEffectDataTable[ 92].fParameterFitX1    =        3.3981;
  fDensityEffectDataTable[ 92].fParameterFitA     =       0.20419;
  fDensityEffectDataTable[ 92].fParameterFitM     =        2.7679;
  fDensityEffectDataTable[ 92].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 92].fDeltaErrorMax     =         0.057;
  fDensityEffectDataTable[ 92].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Am ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 93].fName              = "NIST_MAT_Am";
  fDensityEffectDataTable[ 93].fPlasmaEnergy      =        66.607 * eV;
  fDensityEffectDataTable[ 93].fSternheimerFactor =         1.552;
  fDensityEffectDataTable[ 93].fParameterC        =        6.2813;
  fDensityEffectDataTable[ 93].fParameterFitX0    =        0.2274;
  fDensityEffectDataTable[ 93].fParameterFitX1    =        3.5021;
  fDensityEffectDataTable[ 93].fParameterFitA     =       0.20308;
  fDensityEffectDataTable[ 93].fParameterFitM     =        2.7615;
  fDensityEffectDataTable[ 93].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 93].fDeltaErrorMax     =         0.056;
  fDensityEffectDataTable[ 93].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Cm ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 94].fName              = "NIST_MAT_Cm";
  fDensityEffectDataTable[ 94].fPlasmaEnergy      =        66.022 * eV;
  fDensityEffectDataTable[ 94].fSternheimerFactor =         1.559;
  fDensityEffectDataTable[ 94].fParameterC        =        6.3097;
  fDensityEffectDataTable[ 94].fParameterFitX0    =        0.2484;
  fDensityEffectDataTable[ 94].fParameterFitX1    =         3.516;
  fDensityEffectDataTable[ 94].fParameterFitA     =       0.20257;
  fDensityEffectDataTable[ 94].fParameterFitM     =        2.7579;
  fDensityEffectDataTable[ 94].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 94].fDeltaErrorMax     =         0.056;
  fDensityEffectDataTable[ 94].fState             = MaterialState::kStateSolid;


  // NIST_MAT_Bk ---------------------------------------------------------------------------
  fDensityEffectDataTable[ 95].fName              = "NIST_MAT_Bk";
  fDensityEffectDataTable[ 95].fPlasmaEnergy      =        67.557 * eV;
  fDensityEffectDataTable[ 95].fSternheimerFactor =         1.574;
  fDensityEffectDataTable[ 95].fParameterC        =        6.2912;
  fDensityEffectDataTable[ 95].fParameterFitX0    =        0.2378;
  fDensityEffectDataTable[ 95].fParameterFitX1    =        3.5186;
  fDensityEffectDataTable[ 95].fParameterFitA     =       0.20192;
  fDensityEffectDataTable[ 95].fParameterFitM     =         2.756;
  fDensityEffectDataTable[ 95].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[ 95].fDeltaErrorMax     =         0.062;
  fDensityEffectDataTable[ 95].fState             = MaterialState::kStateSolid;


  // NIST_MAT_A-150_TISSUE -----------------------------------------------------------------
  fDensityEffectDataTable[ 96].fName              = "NIST_MAT_A-150_TISSUE";
  fDensityEffectDataTable[ 96].fPlasmaEnergy      =        22.667 * eV;
  fDensityEffectDataTable[ 96].fSternheimerFactor =          1.95;
  fDensityEffectDataTable[ 96].fParameterC        =          3.11;
  fDensityEffectDataTable[ 96].fParameterFitX0    =        0.1329;
  fDensityEffectDataTable[ 96].fParameterFitX1    =        2.6234;
  fDensityEffectDataTable[ 96].fParameterFitA     =       0.10783;
  fDensityEffectDataTable[ 96].fParameterFitM     =        3.4442;
  fDensityEffectDataTable[ 96].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 96].fDeltaErrorMax     =         0.048;
  fDensityEffectDataTable[ 96].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ACETONE ----------------------------------------------------------------------
  fDensityEffectDataTable[ 97].fName              = "NIST_MAT_ACETONE";
  fDensityEffectDataTable[ 97].fPlasmaEnergy      =         19.01 * eV;
  fDensityEffectDataTable[ 97].fSternheimerFactor =         1.976;
  fDensityEffectDataTable[ 97].fParameterC        =        3.4341;
  fDensityEffectDataTable[ 97].fParameterFitX0    =        0.2197;
  fDensityEffectDataTable[ 97].fParameterFitX1    =        2.6928;
  fDensityEffectDataTable[ 97].fParameterFitA     =         0.111;
  fDensityEffectDataTable[ 97].fParameterFitM     =        3.4047;
  fDensityEffectDataTable[ 97].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 97].fDeltaErrorMax     =         0.069;
  fDensityEffectDataTable[ 97].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ACETYLENE --------------------------------------------------------------------
  fDensityEffectDataTable[ 98].fName              = "NIST_MAT_ACETYLENE";
  fDensityEffectDataTable[ 98].fPlasmaEnergy      =           0.7 * eV;
  fDensityEffectDataTable[ 98].fSternheimerFactor =         1.784;
  fDensityEffectDataTable[ 98].fParameterC        =        9.8419;
  fDensityEffectDataTable[ 98].fParameterFitX0    =        1.6017;
  fDensityEffectDataTable[ 98].fParameterFitX1    =        4.0074;
  fDensityEffectDataTable[ 98].fParameterFitA     =       0.12167;
  fDensityEffectDataTable[ 98].fParameterFitM     =        3.4277;
  fDensityEffectDataTable[ 98].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 98].fDeltaErrorMax     =          0.08;
  fDensityEffectDataTable[ 98].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ADENINE ----------------------------------------------------------------------
  fDensityEffectDataTable[ 99].fName              = "NIST_MAT_ADENINE";
  fDensityEffectDataTable[ 99].fPlasmaEnergy      =        24.098 * eV;
  fDensityEffectDataTable[ 99].fSternheimerFactor =         1.892;
  fDensityEffectDataTable[ 99].fParameterC        =        3.1724;
  fDensityEffectDataTable[ 99].fParameterFitX0    =        0.1295;
  fDensityEffectDataTable[ 99].fParameterFitX1    =        2.4219;
  fDensityEffectDataTable[ 99].fParameterFitA     =       0.20908;
  fDensityEffectDataTable[ 99].fParameterFitM     =        3.0271;
  fDensityEffectDataTable[ 99].fParameterDelta0   =             0;
  fDensityEffectDataTable[ 99].fDeltaErrorMax     =         0.052;
  fDensityEffectDataTable[ 99].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ADIPOSE_TISSUE_ICRP ----------------------------------------------------------
  fDensityEffectDataTable[100].fName              = "NIST_MAT_ADIPOSE_TISSUE_ICRP";
  fDensityEffectDataTable[100].fPlasmaEnergy      =        20.655 * eV;
  fDensityEffectDataTable[100].fSternheimerFactor =         1.987;
  fDensityEffectDataTable[100].fParameterC        =        3.2367;
  fDensityEffectDataTable[100].fParameterFitX0    =        0.1827;
  fDensityEffectDataTable[100].fParameterFitX1    =         2.653;
  fDensityEffectDataTable[100].fParameterFitA     =       0.10278;
  fDensityEffectDataTable[100].fParameterFitM     =        3.4817;
  fDensityEffectDataTable[100].fParameterDelta0   =             0;
  fDensityEffectDataTable[100].fDeltaErrorMax     =          0.06;
  fDensityEffectDataTable[100].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_AIR --------------------------------------------------------------------------
  fDensityEffectDataTable[101].fName              = "NIST_MAT_AIR";
  fDensityEffectDataTable[101].fPlasmaEnergy      =         0.707 * eV;
  fDensityEffectDataTable[101].fSternheimerFactor =         2.054;
  fDensityEffectDataTable[101].fParameterC        =       10.5961;
  fDensityEffectDataTable[101].fParameterFitX0    =        1.7418;
  fDensityEffectDataTable[101].fParameterFitX1    =        4.2759;
  fDensityEffectDataTable[101].fParameterFitA     =       0.10914;
  fDensityEffectDataTable[101].fParameterFitM     =        3.3994;
  fDensityEffectDataTable[101].fParameterDelta0   =             0;
  fDensityEffectDataTable[101].fDeltaErrorMax     =          0.09;
  fDensityEffectDataTable[101].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ALANINE ----------------------------------------------------------------------
  fDensityEffectDataTable[102].fName              = "NIST_MAT_ALANINE";
  fDensityEffectDataTable[102].fPlasmaEnergy      =        25.204 * eV;
  fDensityEffectDataTable[102].fSternheimerFactor =         2.074;
  fDensityEffectDataTable[102].fParameterC        =        3.0965;
  fDensityEffectDataTable[102].fParameterFitX0    =        0.1354;
  fDensityEffectDataTable[102].fParameterFitX1    =        2.6336;
  fDensityEffectDataTable[102].fParameterFitA     =       0.11484;
  fDensityEffectDataTable[102].fParameterFitM     =        3.3526;
  fDensityEffectDataTable[102].fParameterDelta0   =             0;
  fDensityEffectDataTable[102].fDeltaErrorMax     =         0.056;
  fDensityEffectDataTable[102].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ALUMINIM_OXIDE ---------------------------------------------------------------
  fDensityEffectDataTable[103].fName              = "NIST_MAT_ALUMINIM_OXIDE";
  fDensityEffectDataTable[103].fPlasmaEnergy      =        40.206 * eV;
  fDensityEffectDataTable[103].fSternheimerFactor =         2.394;
  fDensityEffectDataTable[103].fParameterC        =        3.5682;
  fDensityEffectDataTable[103].fParameterFitX0    =        0.0402;
  fDensityEffectDataTable[103].fParameterFitX1    =        2.8665;
  fDensityEffectDataTable[103].fParameterFitA     =         0.085;
  fDensityEffectDataTable[103].fParameterFitM     =        3.5458;
  fDensityEffectDataTable[103].fParameterDelta0   =             0;
  fDensityEffectDataTable[103].fDeltaErrorMax     =         0.031;
  fDensityEffectDataTable[103].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_AMBER ------------------------------------------------------------------------
  fDensityEffectDataTable[104].fName              = "NIST_MAT_AMBER";
  fDensityEffectDataTable[104].fPlasmaEnergy      =         22.45 * eV;
  fDensityEffectDataTable[104].fSternheimerFactor =         1.946;
  fDensityEffectDataTable[104].fParameterC        =        3.0701;
  fDensityEffectDataTable[104].fParameterFitX0    =        0.1335;
  fDensityEffectDataTable[104].fParameterFitX1    =         2.561;
  fDensityEffectDataTable[104].fParameterFitA     =       0.11934;
  fDensityEffectDataTable[104].fParameterFitM     =        3.4098;
  fDensityEffectDataTable[104].fParameterDelta0   =             0;
  fDensityEffectDataTable[104].fDeltaErrorMax     =         0.053;
  fDensityEffectDataTable[104].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_AMMONIA ----------------------------------------------------------------------
  fDensityEffectDataTable[105].fName              = "NIST_MAT_AMMONIA";
  fDensityEffectDataTable[105].fPlasmaEnergy      =         0.635 * eV;
  fDensityEffectDataTable[105].fSternheimerFactor =         1.814;
  fDensityEffectDataTable[105].fParameterC        =        9.8763;
  fDensityEffectDataTable[105].fParameterFitX0    =        1.6822;
  fDensityEffectDataTable[105].fParameterFitX1    =        4.1158;
  fDensityEffectDataTable[105].fParameterFitA     =       0.08315;
  fDensityEffectDataTable[105].fParameterFitM     =        3.6464;
  fDensityEffectDataTable[105].fParameterDelta0   =             0;
  fDensityEffectDataTable[105].fDeltaErrorMax     =         0.102;
  fDensityEffectDataTable[105].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ANILINE ----------------------------------------------------------------------
  fDensityEffectDataTable[106].fName              = "NIST_MAT_ANILINE";
  fDensityEffectDataTable[106].fPlasmaEnergy      =        21.361 * eV;
  fDensityEffectDataTable[106].fSternheimerFactor =         1.938;
  fDensityEffectDataTable[106].fParameterC        =        3.2622;
  fDensityEffectDataTable[106].fParameterFitX0    =        0.1618;
  fDensityEffectDataTable[106].fParameterFitX1    =        2.5805;
  fDensityEffectDataTable[106].fParameterFitA     =       0.13134;
  fDensityEffectDataTable[106].fParameterFitM     =        3.3434;
  fDensityEffectDataTable[106].fParameterDelta0   =             0;
  fDensityEffectDataTable[106].fDeltaErrorMax     =         0.052;
  fDensityEffectDataTable[106].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ANTHRACENE -------------------------------------------------------------------
  fDensityEffectDataTable[107].fName              = "NIST_MAT_ANTHRACENE";
  fDensityEffectDataTable[107].fPlasmaEnergy      =        23.704 * eV;
  fDensityEffectDataTable[107].fSternheimerFactor =         1.954;
  fDensityEffectDataTable[107].fParameterC        =        3.1514;
  fDensityEffectDataTable[107].fParameterFitX0    =        0.1146;
  fDensityEffectDataTable[107].fParameterFitX1    =        2.5213;
  fDensityEffectDataTable[107].fParameterFitA     =       0.14677;
  fDensityEffectDataTable[107].fParameterFitM     =        3.2831;
  fDensityEffectDataTable[107].fParameterDelta0   =             0;
  fDensityEffectDataTable[107].fDeltaErrorMax     =         0.042;
  fDensityEffectDataTable[107].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_B-100_BONE -------------------------------------------------------------------
  fDensityEffectDataTable[108].fName              = "NIST_MAT_B-100_BONE";
  fDensityEffectDataTable[108].fPlasmaEnergy      =        25.199 * eV;
  fDensityEffectDataTable[108].fSternheimerFactor =         2.013;
  fDensityEffectDataTable[108].fParameterC        =        3.4528;
  fDensityEffectDataTable[108].fParameterFitX0    =        0.1252;
  fDensityEffectDataTable[108].fParameterFitX1    =         3.042;
  fDensityEffectDataTable[108].fParameterFitA     =       0.05268;
  fDensityEffectDataTable[108].fParameterFitM     =        3.7365;
  fDensityEffectDataTable[108].fParameterDelta0   =             0;
  fDensityEffectDataTable[108].fDeltaErrorMax     =         0.043;
  fDensityEffectDataTable[108].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BAKELITE ---------------------------------------------------------------------
  fDensityEffectDataTable[109].fName              = "NIST_MAT_BAKELITE";
  fDensityEffectDataTable[109].fPlasmaEnergy      =        23.408 * eV;
  fDensityEffectDataTable[109].fSternheimerFactor =         2.046;
  fDensityEffectDataTable[109].fParameterC        =        3.2582;
  fDensityEffectDataTable[109].fParameterFitX0    =        0.1471;
  fDensityEffectDataTable[109].fParameterFitX1    =        2.6055;
  fDensityEffectDataTable[109].fParameterFitA     =       0.12713;
  fDensityEffectDataTable[109].fParameterFitM     =         3.347;
  fDensityEffectDataTable[109].fParameterDelta0   =             0;
  fDensityEffectDataTable[109].fDeltaErrorMax     =         0.052;
  fDensityEffectDataTable[109].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BARIUM_FLUORIDE --------------------------------------------------------------
  fDensityEffectDataTable[110].fName              = "NIST_MAT_BARIUM_FLUORIDE";
  fDensityEffectDataTable[110].fPlasmaEnergy      =        41.398 * eV;
  fDensityEffectDataTable[110].fSternheimerFactor =         1.727;
  fDensityEffectDataTable[110].fParameterC        =        5.4122;
  fDensityEffectDataTable[110].fParameterFitX0    =       -0.0098;
  fDensityEffectDataTable[110].fParameterFitX1    =        3.3871;
  fDensityEffectDataTable[110].fParameterFitA     =       0.15991;
  fDensityEffectDataTable[110].fParameterFitM     =        2.8867;
  fDensityEffectDataTable[110].fParameterDelta0   =             0;
  fDensityEffectDataTable[110].fDeltaErrorMax     =         0.034;
  fDensityEffectDataTable[110].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BARIUM_SULFATE ---------------------------------------------------------------
  fDensityEffectDataTable[111].fName              = "NIST_MAT_BARIUM_SULFATE";
  fDensityEffectDataTable[111].fPlasmaEnergy      =        40.805 * eV;
  fDensityEffectDataTable[111].fSternheimerFactor =         1.893;
  fDensityEffectDataTable[111].fParameterC        =        4.8923;
  fDensityEffectDataTable[111].fParameterFitX0    =       -0.0128;
  fDensityEffectDataTable[111].fParameterFitX1    =        3.4069;
  fDensityEffectDataTable[111].fParameterFitA     =       0.11747;
  fDensityEffectDataTable[111].fParameterFitM     =        3.0427;
  fDensityEffectDataTable[111].fParameterDelta0   =             0;
  fDensityEffectDataTable[111].fDeltaErrorMax     =          0.03;
  fDensityEffectDataTable[111].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BENZENE ----------------------------------------------------------------------
  fDensityEffectDataTable[112].fName              = "NIST_MAT_BENZENE";
  fDensityEffectDataTable[112].fPlasmaEnergy      =        19.806 * eV;
  fDensityEffectDataTable[112].fSternheimerFactor =         1.873;
  fDensityEffectDataTable[112].fParameterC        =        3.3269;
  fDensityEffectDataTable[112].fParameterFitX0    =         0.171;
  fDensityEffectDataTable[112].fParameterFitX1    =        2.5091;
  fDensityEffectDataTable[112].fParameterFitA     =       0.16519;
  fDensityEffectDataTable[112].fParameterFitM     =        3.2174;
  fDensityEffectDataTable[112].fParameterDelta0   =             0;
  fDensityEffectDataTable[112].fDeltaErrorMax     =         0.052;
  fDensityEffectDataTable[112].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BERYLLIUM_OXIDE --------------------------------------------------------------
  fDensityEffectDataTable[113].fName              = "NIST_MAT_BERYLLIUM_OXIDE";
  fDensityEffectDataTable[113].fPlasmaEnergy      =        34.629 * eV;
  fDensityEffectDataTable[113].fSternheimerFactor =         2.296;
  fDensityEffectDataTable[113].fParameterC        =        2.9801;
  fDensityEffectDataTable[113].fParameterFitX0    =        0.0241;
  fDensityEffectDataTable[113].fParameterFitX1    =        2.5846;
  fDensityEffectDataTable[113].fParameterFitA     =       0.10755;
  fDensityEffectDataTable[113].fParameterFitM     =        3.4927;
  fDensityEffectDataTable[113].fParameterDelta0   =             0;
  fDensityEffectDataTable[113].fDeltaErrorMax     =         0.031;
  fDensityEffectDataTable[113].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BGO --------------------------------------------------------------------------
  fDensityEffectDataTable[114].fName              = "NIST_MAT_BGO";
  fDensityEffectDataTable[114].fPlasmaEnergy      =        49.904 * eV;
  fDensityEffectDataTable[114].fSternheimerFactor =         2.121;
  fDensityEffectDataTable[114].fParameterC        =        5.7409;
  fDensityEffectDataTable[114].fParameterFitX0    =        0.0456;
  fDensityEffectDataTable[114].fParameterFitX1    =        3.7816;
  fDensityEffectDataTable[114].fParameterFitA     =       0.09569;
  fDensityEffectDataTable[114].fParameterFitM     =        3.0781;
  fDensityEffectDataTable[114].fParameterDelta0   =             0;
  fDensityEffectDataTable[114].fDeltaErrorMax     =         0.023;
  fDensityEffectDataTable[114].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BLOOD_ICRP -------------------------------------------------------------------
  fDensityEffectDataTable[115].fName              = "NIST_MAT_BLOOD_ICRP";
  fDensityEffectDataTable[115].fPlasmaEnergy      =        22.001 * eV;
  fDensityEffectDataTable[115].fSternheimerFactor =         2.184;
  fDensityEffectDataTable[115].fParameterC        =        3.4581;
  fDensityEffectDataTable[115].fParameterFitX0    =        0.2239;
  fDensityEffectDataTable[115].fParameterFitX1    =        2.8017;
  fDensityEffectDataTable[115].fParameterFitA     =       0.08492;
  fDensityEffectDataTable[115].fParameterFitM     =        3.5406;
  fDensityEffectDataTable[115].fParameterDelta0   =             0;
  fDensityEffectDataTable[115].fDeltaErrorMax     =         0.088;
  fDensityEffectDataTable[115].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BONE_COMPACT_ICRU ------------------------------------------------------------
  fDensityEffectDataTable[116].fName              = "NIST_MAT_BONE_COMPACT_ICRU";
  fDensityEffectDataTable[116].fPlasmaEnergy      =        28.536 * eV;
  fDensityEffectDataTable[116].fSternheimerFactor =         2.091;
  fDensityEffectDataTable[116].fParameterC        =         3.339;
  fDensityEffectDataTable[116].fParameterFitX0    =        0.0944;
  fDensityEffectDataTable[116].fParameterFitX1    =        3.0201;
  fDensityEffectDataTable[116].fParameterFitA     =       0.05822;
  fDensityEffectDataTable[116].fParameterFitM     =        3.6419;
  fDensityEffectDataTable[116].fParameterDelta0   =             0;
  fDensityEffectDataTable[116].fDeltaErrorMax     =         0.042;
  fDensityEffectDataTable[116].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BONE_CORTICAL_ICRP -----------------------------------------------------------
  fDensityEffectDataTable[117].fName              = "NIST_MAT_BONE_CORTICAL_ICRP";
  fDensityEffectDataTable[117].fPlasmaEnergy      =        28.298 * eV;
  fDensityEffectDataTable[117].fSternheimerFactor =         2.118;
  fDensityEffectDataTable[117].fParameterC        =        3.6488;
  fDensityEffectDataTable[117].fParameterFitX0    =        0.1161;
  fDensityEffectDataTable[117].fParameterFitX1    =        3.0919;
  fDensityEffectDataTable[117].fParameterFitA     =       0.06198;
  fDensityEffectDataTable[117].fParameterFitM     =        3.5919;
  fDensityEffectDataTable[117].fParameterDelta0   =             0;
  fDensityEffectDataTable[117].fDeltaErrorMax     =          0.04;
  fDensityEffectDataTable[117].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BORON_CARBIDE ----------------------------------------------------------------
  fDensityEffectDataTable[118].fName              = "NIST_MAT_BORON_CARBIDE";
  fDensityEffectDataTable[118].fPlasmaEnergy      =         31.38 * eV;
  fDensityEffectDataTable[118].fSternheimerFactor =          2.14;
  fDensityEffectDataTable[118].fParameterC        =        2.9859;
  fDensityEffectDataTable[118].fParameterFitX0    =        0.0093;
  fDensityEffectDataTable[118].fParameterFitX1    =        2.1006;
  fDensityEffectDataTable[118].fParameterFitA     =       0.37087;
  fDensityEffectDataTable[118].fParameterFitM     =        2.8076;
  fDensityEffectDataTable[118].fParameterDelta0   =             0;
  fDensityEffectDataTable[118].fDeltaErrorMax     =         0.022;
  fDensityEffectDataTable[118].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BORON_OXIDE ------------------------------------------------------------------
  fDensityEffectDataTable[119].fName              = "NIST_MAT_BORON_OXIDE";
  fDensityEffectDataTable[119].fPlasmaEnergy      =        27.107 * eV;
  fDensityEffectDataTable[119].fSternheimerFactor =         2.446;
  fDensityEffectDataTable[119].fParameterC        =        3.6027;
  fDensityEffectDataTable[119].fParameterFitX0    =        0.1843;
  fDensityEffectDataTable[119].fParameterFitX1    =        2.7379;
  fDensityEffectDataTable[119].fParameterFitA     =       0.11548;
  fDensityEffectDataTable[119].fParameterFitM     =        3.3832;
  fDensityEffectDataTable[119].fParameterDelta0   =             0;
  fDensityEffectDataTable[119].fDeltaErrorMax     =         0.053;
  fDensityEffectDataTable[119].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BRAIN_ICRP -------------------------------------------------------------------
  fDensityEffectDataTable[120].fName              = "NIST_MAT_BRAIN_ICRP";
  fDensityEffectDataTable[120].fPlasmaEnergy      =        21.772 * eV;
  fDensityEffectDataTable[120].fSternheimerFactor =         2.162;
  fDensityEffectDataTable[120].fParameterC        =        3.4279;
  fDensityEffectDataTable[120].fParameterFitX0    =        0.2206;
  fDensityEffectDataTable[120].fParameterFitX1    =        2.8021;
  fDensityEffectDataTable[120].fParameterFitA     =       0.08255;
  fDensityEffectDataTable[120].fParameterFitM     =        3.5585;
  fDensityEffectDataTable[120].fParameterDelta0   =             0;
  fDensityEffectDataTable[120].fDeltaErrorMax     =         0.086;
  fDensityEffectDataTable[120].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_BUTANE -----------------------------------------------------------------------
  fDensityEffectDataTable[121].fName              = "NIST_MAT_BUTANE";
  fDensityEffectDataTable[121].fPlasmaEnergy      =         1.101 * eV;
  fDensityEffectDataTable[121].fSternheimerFactor =         1.727;
  fDensityEffectDataTable[121].fParameterC        =        8.5633;
  fDensityEffectDataTable[121].fParameterFitX0    =        1.3788;
  fDensityEffectDataTable[121].fParameterFitX1    =        3.7524;
  fDensityEffectDataTable[121].fParameterFitA     =       0.10852;
  fDensityEffectDataTable[121].fParameterFitM     =        3.4184;
  fDensityEffectDataTable[121].fParameterDelta0   =             0;
  fDensityEffectDataTable[121].fDeltaErrorMax     =           0.1;
  fDensityEffectDataTable[121].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_N-BUTYL_ALCOHOL --------------------------------------------------------------
  fDensityEffectDataTable[122].fName              = "NIST_MAT_N-BUTYL_ALCOHOL";
  fDensityEffectDataTable[122].fPlasmaEnergy      =         19.52 * eV;
  fDensityEffectDataTable[122].fSternheimerFactor =         1.942;
  fDensityEffectDataTable[122].fParameterC        =        3.2425;
  fDensityEffectDataTable[122].fParameterFitX0    =        0.1937;
  fDensityEffectDataTable[122].fParameterFitX1    =        2.6439;
  fDensityEffectDataTable[122].fParameterFitA     =       0.10081;
  fDensityEffectDataTable[122].fParameterFitM     =        3.5139;
  fDensityEffectDataTable[122].fParameterDelta0   =             0;
  fDensityEffectDataTable[122].fDeltaErrorMax     =         0.065;
  fDensityEffectDataTable[122].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_C-552 ------------------------------------------------------------------------
  fDensityEffectDataTable[123].fName              = "NIST_MAT_C-552";
  fDensityEffectDataTable[123].fPlasmaEnergy      =        27.023 * eV;
  fDensityEffectDataTable[123].fSternheimerFactor =         2.128;
  fDensityEffectDataTable[123].fParameterC        =        3.3338;
  fDensityEffectDataTable[123].fParameterFitX0    =         0.151;
  fDensityEffectDataTable[123].fParameterFitX1    =        2.7083;
  fDensityEffectDataTable[123].fParameterFitA     =       0.10492;
  fDensityEffectDataTable[123].fParameterFitM     =        3.4344;
  fDensityEffectDataTable[123].fParameterDelta0   =             0;
  fDensityEffectDataTable[123].fDeltaErrorMax     =         0.053;
  fDensityEffectDataTable[123].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CADMIUM_TELLURIDE ------------------------------------------------------------
  fDensityEffectDataTable[124].fName              = "NIST_MAT_CADMIUM_TELLURIDE";
  fDensityEffectDataTable[124].fPlasmaEnergy      =        46.314 * eV;
  fDensityEffectDataTable[124].fSternheimerFactor =         1.935;
  fDensityEffectDataTable[124].fParameterC        =        5.9096;
  fDensityEffectDataTable[124].fParameterFitX0    =        0.0438;
  fDensityEffectDataTable[124].fParameterFitX1    =        3.2836;
  fDensityEffectDataTable[124].fParameterFitA     =        0.2484;
  fDensityEffectDataTable[124].fParameterFitM     =        2.6665;
  fDensityEffectDataTable[124].fParameterDelta0   =             0;
  fDensityEffectDataTable[124].fDeltaErrorMax     =         0.057;
  fDensityEffectDataTable[124].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CADMIUM_TUNGSTATE ------------------------------------------------------------
  fDensityEffectDataTable[125].fName              = "NIST_MAT_CADMIUM_TUNGSTATE";
  fDensityEffectDataTable[125].fPlasmaEnergy      =        52.954 * eV;
  fDensityEffectDataTable[125].fSternheimerFactor =         2.289;
  fDensityEffectDataTable[125].fParameterC        =        5.3594;
  fDensityEffectDataTable[125].fParameterFitX0    =        0.0123;
  fDensityEffectDataTable[125].fParameterFitX1    =        3.5941;
  fDensityEffectDataTable[125].fParameterFitA     =       0.12861;
  fDensityEffectDataTable[125].fParameterFitM     =         2.915;
  fDensityEffectDataTable[125].fParameterDelta0   =             0;
  fDensityEffectDataTable[125].fDeltaErrorMax     =         0.027;
  fDensityEffectDataTable[125].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CALCIUM_CARBONATE ------------------------------------------------------------
  fDensityEffectDataTable[126].fName              = "NIST_MAT_CALCIUM_CARBONATE";
  fDensityEffectDataTable[126].fPlasmaEnergy      =         34.08 * eV;
  fDensityEffectDataTable[126].fSternheimerFactor =         2.141;
  fDensityEffectDataTable[126].fParameterC        =        3.7738;
  fDensityEffectDataTable[126].fParameterFitX0    =        0.0492;
  fDensityEffectDataTable[126].fParameterFitX1    =        3.0549;
  fDensityEffectDataTable[126].fParameterFitA     =       0.08301;
  fDensityEffectDataTable[126].fParameterFitM     =         3.412;
  fDensityEffectDataTable[126].fParameterDelta0   =             0;
  fDensityEffectDataTable[126].fDeltaErrorMax     =         0.037;
  fDensityEffectDataTable[126].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CALCIUM_FLUORIDE -------------------------------------------------------------
  fDensityEffectDataTable[127].fName              = "NIST_MAT_CALCIUM_FLUORIDE";
  fDensityEffectDataTable[127].fPlasmaEnergy      =        35.849 * eV;
  fDensityEffectDataTable[127].fSternheimerFactor =         2.127;
  fDensityEffectDataTable[127].fParameterC        =        4.0653;
  fDensityEffectDataTable[127].fParameterFitX0    =        0.0676;
  fDensityEffectDataTable[127].fParameterFitX1    =        3.1683;
  fDensityEffectDataTable[127].fParameterFitA     =       0.06942;
  fDensityEffectDataTable[127].fParameterFitM     =        3.5263;
  fDensityEffectDataTable[127].fParameterDelta0   =             0;
  fDensityEffectDataTable[127].fDeltaErrorMax     =         0.044;
  fDensityEffectDataTable[127].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CALCIUM_OXIDE ----------------------------------------------------------------
  fDensityEffectDataTable[128].fName              = "NIST_MAT_CALCIUM_OXIDE";
  fDensityEffectDataTable[128].fPlasmaEnergy      =        36.988 * eV;
  fDensityEffectDataTable[128].fSternheimerFactor =         1.973;
  fDensityEffectDataTable[128].fParameterC        =        4.1209;
  fDensityEffectDataTable[128].fParameterFitX0    =       -0.0172;
  fDensityEffectDataTable[128].fParameterFitX1    =        3.0171;
  fDensityEffectDataTable[128].fParameterFitA     =       0.12128;
  fDensityEffectDataTable[128].fParameterFitM     =        3.1936;
  fDensityEffectDataTable[128].fParameterDelta0   =             0;
  fDensityEffectDataTable[128].fDeltaErrorMax     =         0.024;
  fDensityEffectDataTable[128].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CALCIUM_SULFATE --------------------------------------------------------------
  fDensityEffectDataTable[129].fName              = "NIST_MAT_CALCIUM_SULFATE";
  fDensityEffectDataTable[129].fPlasmaEnergy      =        35.038 * eV;
  fDensityEffectDataTable[129].fSternheimerFactor =         2.179;
  fDensityEffectDataTable[129].fParameterC        =        3.9388;
  fDensityEffectDataTable[129].fParameterFitX0    =        0.0587;
  fDensityEffectDataTable[129].fParameterFitX1    =        3.1229;
  fDensityEffectDataTable[129].fParameterFitA     =       0.07708;
  fDensityEffectDataTable[129].fParameterFitM     =        3.4495;
  fDensityEffectDataTable[129].fParameterDelta0   =             0;
  fDensityEffectDataTable[129].fDeltaErrorMax     =         0.021;
  fDensityEffectDataTable[129].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CALCIUM_TUNGSTATE ------------------------------------------------------------
  fDensityEffectDataTable[130].fName              = "NIST_MAT_CALCIUM_TUNGSTATE";
  fDensityEffectDataTable[130].fPlasmaEnergy      =        46.934 * eV;
  fDensityEffectDataTable[130].fSternheimerFactor =         2.262;
  fDensityEffectDataTable[130].fParameterC        =        5.2603;
  fDensityEffectDataTable[130].fParameterFitX0    =        0.0323;
  fDensityEffectDataTable[130].fParameterFitX1    =        3.8932;
  fDensityEffectDataTable[130].fParameterFitA     =        0.0621;
  fDensityEffectDataTable[130].fParameterFitM     =        3.2649;
  fDensityEffectDataTable[130].fParameterDelta0   =             0;
  fDensityEffectDataTable[130].fDeltaErrorMax     =         0.021;
  fDensityEffectDataTable[130].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CARBON_DIOXIDE ---------------------------------------------------------------
  fDensityEffectDataTable[131].fName              = "NIST_MAT_CARBON_DIOXIDE";
  fDensityEffectDataTable[131].fPlasmaEnergy      =         0.874 * eV;
  fDensityEffectDataTable[131].fSternheimerFactor =         2.118;
  fDensityEffectDataTable[131].fParameterC        =       10.1537;
  fDensityEffectDataTable[131].fParameterFitX0    =        1.6294;
  fDensityEffectDataTable[131].fParameterFitX1    =        4.1825;
  fDensityEffectDataTable[131].fParameterFitA     =       0.11768;
  fDensityEffectDataTable[131].fParameterFitM     =        3.3227;
  fDensityEffectDataTable[131].fParameterDelta0   =             0;
  fDensityEffectDataTable[131].fDeltaErrorMax     =         0.091;
  fDensityEffectDataTable[131].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CARBON_TETRACHLORIDE ---------------------------------------------------------
  fDensityEffectDataTable[132].fName              = "NIST_MAT_CARBON_TETRACHLORIDE";
  fDensityEffectDataTable[132].fPlasmaEnergy      =        25.234 * eV;
  fDensityEffectDataTable[132].fSternheimerFactor =         1.742;
  fDensityEffectDataTable[132].fParameterC        =        4.7712;
  fDensityEffectDataTable[132].fParameterFitX0    =        0.1773;
  fDensityEffectDataTable[132].fParameterFitX1    =        2.9165;
  fDensityEffectDataTable[132].fParameterFitA     =       0.19018;
  fDensityEffectDataTable[132].fParameterFitM     =        3.0116;
  fDensityEffectDataTable[132].fParameterDelta0   =             0;
  fDensityEffectDataTable[132].fDeltaErrorMax     =         0.041;
  fDensityEffectDataTable[132].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CELLULOSE_CELLOPHANE ---------------------------------------------------------
  fDensityEffectDataTable[133].fName              = "NIST_MAT_CELLULOSE_CELLOPHANE";
  fDensityEffectDataTable[133].fPlasmaEnergy      =        25.008 * eV;
  fDensityEffectDataTable[133].fSternheimerFactor =          2.17;
  fDensityEffectDataTable[133].fParameterC        =        3.2647;
  fDensityEffectDataTable[133].fParameterFitX0    =         0.158;
  fDensityEffectDataTable[133].fParameterFitX1    =        2.6778;
  fDensityEffectDataTable[133].fParameterFitA     =       0.11151;
  fDensityEffectDataTable[133].fParameterFitM     =         3.381;
  fDensityEffectDataTable[133].fParameterDelta0   =             0;
  fDensityEffectDataTable[133].fDeltaErrorMax     =          0.06;
  fDensityEffectDataTable[133].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CELLULOSE_BUTYRATE -----------------------------------------------------------
  fDensityEffectDataTable[134].fName              = "NIST_MAT_CELLULOSE_BUTYRATE";
  fDensityEffectDataTable[134].fPlasmaEnergy      =        23.041 * eV;
  fDensityEffectDataTable[134].fSternheimerFactor =         2.128;
  fDensityEffectDataTable[134].fParameterC        =        3.3497;
  fDensityEffectDataTable[134].fParameterFitX0    =        0.1794;
  fDensityEffectDataTable[134].fParameterFitX1    =        2.6809;
  fDensityEffectDataTable[134].fParameterFitA     =       0.11444;
  fDensityEffectDataTable[134].fParameterFitM     =        3.3738;
  fDensityEffectDataTable[134].fParameterDelta0   =             0;
  fDensityEffectDataTable[134].fDeltaErrorMax     =         0.056;
  fDensityEffectDataTable[134].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CELLULOSE_NITRATE ------------------------------------------------------------
  fDensityEffectDataTable[135].fName              = "NIST_MAT_CELLULOSE_NITRATE";
  fDensityEffectDataTable[135].fPlasmaEnergy      =        25.224 * eV;
  fDensityEffectDataTable[135].fSternheimerFactor =         2.252;
  fDensityEffectDataTable[135].fParameterC        =        3.4762;
  fDensityEffectDataTable[135].fParameterFitX0    =        0.1897;
  fDensityEffectDataTable[135].fParameterFitX1    =        2.7253;
  fDensityEffectDataTable[135].fParameterFitA     =       0.11813;
  fDensityEffectDataTable[135].fParameterFitM     =        3.3237;
  fDensityEffectDataTable[135].fParameterDelta0   =             0;
  fDensityEffectDataTable[135].fDeltaErrorMax     =         0.063;
  fDensityEffectDataTable[135].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CERIC_SULFATE ----------------------------------------------------------------
  fDensityEffectDataTable[136].fName              = "NIST_MAT_CERIC_SULFATE";
  fDensityEffectDataTable[136].fPlasmaEnergy      =        21.743 * eV;
  fDensityEffectDataTable[136].fSternheimerFactor =         2.205;
  fDensityEffectDataTable[136].fParameterC        =        3.5212;
  fDensityEffectDataTable[136].fParameterFitX0    =        0.2363;
  fDensityEffectDataTable[136].fParameterFitX1    =        2.8769;
  fDensityEffectDataTable[136].fParameterFitA     =       0.07666;
  fDensityEffectDataTable[136].fParameterFitM     =        3.5607;
  fDensityEffectDataTable[136].fParameterDelta0   =             0;
  fDensityEffectDataTable[136].fDeltaErrorMax     =         0.095;
  fDensityEffectDataTable[136].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CESIUM_FLUORIDE --------------------------------------------------------------
  fDensityEffectDataTable[137].fName              = "NIST_MAT_CESIUM_FLUORIDE";
  fDensityEffectDataTable[137].fPlasmaEnergy      =        37.942 * eV;
  fDensityEffectDataTable[137].fSternheimerFactor =         1.714;
  fDensityEffectDataTable[137].fParameterC        =        5.9046;
  fDensityEffectDataTable[137].fParameterFitX0    =        0.0084;
  fDensityEffectDataTable[137].fParameterFitX1    =        3.3374;
  fDensityEffectDataTable[137].fParameterFitA     =       0.22052;
  fDensityEffectDataTable[137].fParameterFitM     =         2.728;
  fDensityEffectDataTable[137].fParameterDelta0   =             0;
  fDensityEffectDataTable[137].fDeltaErrorMax     =         0.044;
  fDensityEffectDataTable[137].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CESIUM_IODIDE ----------------------------------------------------------------
  fDensityEffectDataTable[138].fName              = "NIST_MAT_CESIUM_IODIDE";
  fDensityEffectDataTable[138].fPlasmaEnergy      =        39.455 * eV;
  fDensityEffectDataTable[138].fSternheimerFactor =         1.672;
  fDensityEffectDataTable[138].fParameterC        =        6.2807;
  fDensityEffectDataTable[138].fParameterFitX0    =        0.0395;
  fDensityEffectDataTable[138].fParameterFitX1    =        3.3353;
  fDensityEffectDataTable[138].fParameterFitA     =       0.25381;
  fDensityEffectDataTable[138].fParameterFitM     =        2.6657;
  fDensityEffectDataTable[138].fParameterDelta0   =             0;
  fDensityEffectDataTable[138].fDeltaErrorMax     =         0.067;
  fDensityEffectDataTable[138].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CHLOROBENZENE ----------------------------------------------------------------
  fDensityEffectDataTable[139].fName              = "NIST_MAT_CHLOROBENZENE";
  fDensityEffectDataTable[139].fPlasmaEnergy      =        21.752 * eV;
  fDensityEffectDataTable[139].fSternheimerFactor =         1.889;
  fDensityEffectDataTable[139].fParameterC        =        3.8201;
  fDensityEffectDataTable[139].fParameterFitX0    =        0.1714;
  fDensityEffectDataTable[139].fParameterFitX1    =        2.9272;
  fDensityEffectDataTable[139].fParameterFitA     =       0.09856;
  fDensityEffectDataTable[139].fParameterFitM     =        3.3797;
  fDensityEffectDataTable[139].fParameterDelta0   =             0;
  fDensityEffectDataTable[139].fDeltaErrorMax     =         0.031;
  fDensityEffectDataTable[139].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CHLOROFORM -------------------------------------------------------------------
  fDensityEffectDataTable[140].fName              = "NIST_MAT_CHLOROFORM";
  fDensityEffectDataTable[140].fPlasmaEnergy      =        24.462 * eV;
  fDensityEffectDataTable[140].fSternheimerFactor =         1.734;
  fDensityEffectDataTable[140].fParameterC        =        4.7055;
  fDensityEffectDataTable[140].fParameterFitX0    =        0.1786;
  fDensityEffectDataTable[140].fParameterFitX1    =        2.9581;
  fDensityEffectDataTable[140].fParameterFitA     =       0.16959;
  fDensityEffectDataTable[140].fParameterFitM     =        3.0627;
  fDensityEffectDataTable[140].fParameterDelta0   =             0;
  fDensityEffectDataTable[140].fDeltaErrorMax     =         0.038;
  fDensityEffectDataTable[140].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CONCRETE ---------------------------------------------------------------------
  fDensityEffectDataTable[141].fName              = "NIST_MAT_CONCRETE";
  fDensityEffectDataTable[141].fPlasmaEnergy      =        30.986 * eV;
  fDensityEffectDataTable[141].fSternheimerFactor =         2.322;
  fDensityEffectDataTable[141].fParameterC        =        3.9464;
  fDensityEffectDataTable[141].fParameterFitX0    =        0.1301;
  fDensityEffectDataTable[141].fParameterFitX1    =        3.0466;
  fDensityEffectDataTable[141].fParameterFitA     =       0.07515;
  fDensityEffectDataTable[141].fParameterFitM     =        3.5467;
  fDensityEffectDataTable[141].fParameterDelta0   =             0;
  fDensityEffectDataTable[141].fDeltaErrorMax     =         0.024;
  fDensityEffectDataTable[141].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_CYCLOHEXANE ------------------------------------------------------------------
  fDensityEffectDataTable[142].fName              = "NIST_MAT_CYCLOHEXANE";
  fDensityEffectDataTable[142].fPlasmaEnergy      =        19.207 * eV;
  fDensityEffectDataTable[142].fSternheimerFactor =         1.861;
  fDensityEffectDataTable[142].fParameterC        =        3.1544;
  fDensityEffectDataTable[142].fParameterFitX0    =        0.1728;
  fDensityEffectDataTable[142].fParameterFitX1    =        2.5549;
  fDensityEffectDataTable[142].fParameterFitA     =       0.12035;
  fDensityEffectDataTable[142].fParameterFitM     =        3.4278;
  fDensityEffectDataTable[142].fParameterDelta0   =             0;
  fDensityEffectDataTable[142].fDeltaErrorMax     =         0.057;
  fDensityEffectDataTable[142].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_1,2-DICHLOROBENZENE ----------------------------------------------------------
  fDensityEffectDataTable[143].fName              = "NIST_MAT_1,2-DICHLOROBENZENE";
  fDensityEffectDataTable[143].fPlasmaEnergy      =        23.354 * eV;
  fDensityEffectDataTable[143].fSternheimerFactor =         1.862;
  fDensityEffectDataTable[143].fParameterC        =        4.0348;
  fDensityEffectDataTable[143].fParameterFitX0    =        0.1587;
  fDensityEffectDataTable[143].fParameterFitX1    =        2.8276;
  fDensityEffectDataTable[143].fParameterFitA     =        0.1601;
  fDensityEffectDataTable[143].fParameterFitM     =        3.0836;
  fDensityEffectDataTable[143].fParameterDelta0   =             0;
  fDensityEffectDataTable[143].fDeltaErrorMax     =         0.029;
  fDensityEffectDataTable[143].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_DICHLORODIETHYL_ETHER --------------------------------------------------------
  fDensityEffectDataTable[144].fName              = "NIST_MAT_DICHLORODIETHYL_ETHER";
  fDensityEffectDataTable[144].fPlasmaEnergy      =        22.894 * eV;
  fDensityEffectDataTable[144].fSternheimerFactor =         1.903;
  fDensityEffectDataTable[144].fParameterC        =        4.0135;
  fDensityEffectDataTable[144].fParameterFitX0    =        0.1773;
  fDensityEffectDataTable[144].fParameterFitX1    =        3.1586;
  fDensityEffectDataTable[144].fParameterFitA     =       0.06799;
  fDensityEffectDataTable[144].fParameterFitM     =         3.525;
  fDensityEffectDataTable[144].fParameterDelta0   =             0;
  fDensityEffectDataTable[144].fDeltaErrorMax     =         0.026;
  fDensityEffectDataTable[144].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_1,2-DICHLOROETHANE -----------------------------------------------------------
  fDensityEffectDataTable[145].fName              = "NIST_MAT_1,2-DICHLOROETHANE";
  fDensityEffectDataTable[145].fPlasmaEnergy      =        22.764 * eV;
  fDensityEffectDataTable[145].fSternheimerFactor =         1.618;
  fDensityEffectDataTable[145].fParameterC        =        4.1849;
  fDensityEffectDataTable[145].fParameterFitX0    =        0.1375;
  fDensityEffectDataTable[145].fParameterFitX1    =        2.9529;
  fDensityEffectDataTable[145].fParameterFitA     =       0.13383;
  fDensityEffectDataTable[145].fParameterFitM     =        3.1675;
  fDensityEffectDataTable[145].fParameterDelta0   =             0;
  fDensityEffectDataTable[145].fDeltaErrorMax     =          0.03;
  fDensityEffectDataTable[145].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_DIETHYL_ETHER ----------------------------------------------------------------
  fDensityEffectDataTable[146].fName              = "NIST_MAT_DIETHYL_ETHER";
  fDensityEffectDataTable[146].fPlasmaEnergy      =        18.326 * eV;
  fDensityEffectDataTable[146].fSternheimerFactor =         1.951;
  fDensityEffectDataTable[146].fParameterC        =        3.3721;
  fDensityEffectDataTable[146].fParameterFitX0    =        0.2231;
  fDensityEffectDataTable[146].fParameterFitX1    =        2.6745;
  fDensityEffectDataTable[146].fParameterFitA     =        0.1055;
  fDensityEffectDataTable[146].fParameterFitM     =        3.4586;
  fDensityEffectDataTable[146].fParameterDelta0   =             0;
  fDensityEffectDataTable[146].fDeltaErrorMax     =          0.07;
  fDensityEffectDataTable[146].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_N,N-DIMETHYL_FORMAMIDE -------------------------------------------------------
  fDensityEffectDataTable[147].fName              = "NIST_MAT_N,N-DIMETHYL_FORMAMIDE";
  fDensityEffectDataTable[147].fPlasmaEnergy      =        20.763 * eV;
  fDensityEffectDataTable[147].fSternheimerFactor =         2.005;
  fDensityEffectDataTable[147].fParameterC        =        3.3311;
  fDensityEffectDataTable[147].fParameterFitX0    =        0.1977;
  fDensityEffectDataTable[147].fParameterFitX1    =        2.6686;
  fDensityEffectDataTable[147].fParameterFitA     =        0.1147;
  fDensityEffectDataTable[147].fParameterFitM     =         3.371;
  fDensityEffectDataTable[147].fParameterDelta0   =             0;
  fDensityEffectDataTable[147].fDeltaErrorMax     =         0.065;
  fDensityEffectDataTable[147].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_DIMETHYL_SULFOXIDE -----------------------------------------------------------
  fDensityEffectDataTable[148].fName              = "NIST_MAT_DIMETHYL_SULFOXIDE";
  fDensityEffectDataTable[148].fPlasmaEnergy      =        22.173 * eV;
  fDensityEffectDataTable[148].fSternheimerFactor =         2.075;
  fDensityEffectDataTable[148].fParameterC        =        3.9844;
  fDensityEffectDataTable[148].fParameterFitX0    =        0.2021;
  fDensityEffectDataTable[148].fParameterFitX1    =        3.1263;
  fDensityEffectDataTable[148].fParameterFitA     =       0.06619;
  fDensityEffectDataTable[148].fParameterFitM     =        3.5708;
  fDensityEffectDataTable[148].fParameterDelta0   =             0;
  fDensityEffectDataTable[148].fDeltaErrorMax     =          0.03;
  fDensityEffectDataTable[148].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ETHANE -----------------------------------------------------------------------
  fDensityEffectDataTable[149].fName              = "NIST_MAT_ETHANE";
  fDensityEffectDataTable[149].fPlasmaEnergy      =         0.789 * eV;
  fDensityEffectDataTable[149].fSternheimerFactor =          1.69;
  fDensityEffectDataTable[149].fParameterC        =        9.1043;
  fDensityEffectDataTable[149].fParameterFitX0    =        1.5107;
  fDensityEffectDataTable[149].fParameterFitX1    =        3.8743;
  fDensityEffectDataTable[149].fParameterFitA     =       0.09627;
  fDensityEffectDataTable[149].fParameterFitM     =        3.6095;
  fDensityEffectDataTable[149].fParameterDelta0   =             0;
  fDensityEffectDataTable[149].fDeltaErrorMax     =         0.097;
  fDensityEffectDataTable[149].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ETHYL_ALCOHOL ----------------------------------------------------------------
  fDensityEffectDataTable[150].fName              = "NIST_MAT_ETHYL_ALCOHOL";
  fDensityEffectDataTable[150].fPlasmaEnergy      =        19.232 * eV;
  fDensityEffectDataTable[150].fSternheimerFactor =         2.013;
  fDensityEffectDataTable[150].fParameterC        =        3.3699;
  fDensityEffectDataTable[150].fParameterFitX0    =        0.2218;
  fDensityEffectDataTable[150].fParameterFitX1    =        2.7052;
  fDensityEffectDataTable[150].fParameterFitA     =       0.09878;
  fDensityEffectDataTable[150].fParameterFitM     =        3.4834;
  fDensityEffectDataTable[150].fParameterDelta0   =             0;
  fDensityEffectDataTable[150].fDeltaErrorMax     =         0.071;
  fDensityEffectDataTable[150].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ETHYL_CELLULOSE --------------------------------------------------------------
  fDensityEffectDataTable[151].fName              = "NIST_MAT_ETHYL_CELLULOSE";
  fDensityEffectDataTable[151].fPlasmaEnergy      =        22.594 * eV;
  fDensityEffectDataTable[151].fSternheimerFactor =         2.065;
  fDensityEffectDataTable[151].fParameterC        =        3.2415;
  fDensityEffectDataTable[151].fParameterFitX0    =        0.1683;
  fDensityEffectDataTable[151].fParameterFitX1    =        2.6527;
  fDensityEffectDataTable[151].fParameterFitA     =       0.11077;
  fDensityEffectDataTable[151].fParameterFitM     =        3.4098;
  fDensityEffectDataTable[151].fParameterDelta0   =             0;
  fDensityEffectDataTable[151].fDeltaErrorMax     =         0.057;
  fDensityEffectDataTable[151].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ETHYLENE ---------------------------------------------------------------------
  fDensityEffectDataTable[152].fName              = "NIST_MAT_ETHYLENE";
  fDensityEffectDataTable[152].fPlasmaEnergy      =         0.746 * eV;
  fDensityEffectDataTable[152].fSternheimerFactor =         1.733;
  fDensityEffectDataTable[152].fParameterC        =         9.438;
  fDensityEffectDataTable[152].fParameterFitX0    =        1.5528;
  fDensityEffectDataTable[152].fParameterFitX1    =        3.9327;
  fDensityEffectDataTable[152].fParameterFitA     =       0.10636;
  fDensityEffectDataTable[152].fParameterFitM     =        3.5387;
  fDensityEffectDataTable[152].fParameterDelta0   =             0;
  fDensityEffectDataTable[152].fDeltaErrorMax     =         0.085;
  fDensityEffectDataTable[152].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_EYE_LENS_ICRP ----------------------------------------------------------------
  fDensityEffectDataTable[153].fName              = "NIST_MAT_EYE_LENS_ICRP";
  fDensityEffectDataTable[153].fPlasmaEnergy      =        22.388 * eV;
  fDensityEffectDataTable[153].fSternheimerFactor =         2.154;
  fDensityEffectDataTable[153].fParameterC        =         3.372;
  fDensityEffectDataTable[153].fParameterFitX0    =         0.207;
  fDensityEffectDataTable[153].fParameterFitX1    =        2.7446;
  fDensityEffectDataTable[153].fParameterFitA     =        0.0969;
  fDensityEffectDataTable[153].fParameterFitM     =         3.455;
  fDensityEffectDataTable[153].fParameterDelta0   =             0;
  fDensityEffectDataTable[153].fDeltaErrorMax     =         0.077;
  fDensityEffectDataTable[153].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_FERRIC_OXIDE -----------------------------------------------------------------
  fDensityEffectDataTable[154].fName              = "NIST_MAT_FERRIC_OXIDE";
  fDensityEffectDataTable[154].fPlasmaEnergy      =        45.331 * eV;
  fDensityEffectDataTable[154].fSternheimerFactor =         2.747;
  fDensityEffectDataTable[154].fParameterC        =        4.2245;
  fDensityEffectDataTable[154].fParameterFitX0    =       -0.0074;
  fDensityEffectDataTable[154].fParameterFitX1    =        3.2573;
  fDensityEffectDataTable[154].fParameterFitA     =       0.10478;
  fDensityEffectDataTable[154].fParameterFitM     =        3.1313;
  fDensityEffectDataTable[154].fParameterDelta0   =             0;
  fDensityEffectDataTable[154].fDeltaErrorMax     =         0.026;
  fDensityEffectDataTable[154].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_FERROBORIDE ------------------------------------------------------------------
  fDensityEffectDataTable[155].fName              = "NIST_MAT_FERROBORIDE";
  fDensityEffectDataTable[155].fPlasmaEnergy      =        52.546 * eV;
  fDensityEffectDataTable[155].fSternheimerFactor =         2.726;
  fDensityEffectDataTable[155].fParameterC        =        4.2057;
  fDensityEffectDataTable[155].fParameterFitX0    =       -0.0988;
  fDensityEffectDataTable[155].fParameterFitX1    =        3.1749;
  fDensityEffectDataTable[155].fParameterFitA     =       0.12911;
  fDensityEffectDataTable[155].fParameterFitM     =         3.024;
  fDensityEffectDataTable[155].fParameterDelta0   =             0;
  fDensityEffectDataTable[155].fDeltaErrorMax     =         0.022;
  fDensityEffectDataTable[155].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_FERROUS_OXIDE ----------------------------------------------------------------
  fDensityEffectDataTable[156].fName              = "NIST_MAT_FERROUS_OXIDE";
  fDensityEffectDataTable[156].fPlasmaEnergy      =        47.327 * eV;
  fDensityEffectDataTable[156].fSternheimerFactor =         2.769;
  fDensityEffectDataTable[156].fParameterC        =        4.3175;
  fDensityEffectDataTable[156].fParameterFitX0    =       -0.0279;
  fDensityEffectDataTable[156].fParameterFitX1    =        3.2002;
  fDensityEffectDataTable[156].fParameterFitA     =       0.12959;
  fDensityEffectDataTable[156].fParameterFitM     =        3.0168;
  fDensityEffectDataTable[156].fParameterDelta0   =             0;
  fDensityEffectDataTable[156].fDeltaErrorMax     =         0.022;
  fDensityEffectDataTable[156].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_FERROUS_SULFATE --------------------------------------------------------------
  fDensityEffectDataTable[157].fName              = "NIST_MAT_FERROUS_SULFATE";
  fDensityEffectDataTable[157].fPlasmaEnergy      =         21.69 * eV;
  fDensityEffectDataTable[157].fSternheimerFactor =         2.208;
  fDensityEffectDataTable[157].fParameterC        =        3.5183;
  fDensityEffectDataTable[157].fParameterFitX0    =        0.2378;
  fDensityEffectDataTable[157].fParameterFitX1    =        2.8254;
  fDensityEffectDataTable[157].fParameterFitA     =       0.08759;
  fDensityEffectDataTable[157].fParameterFitM     =        3.4923;
  fDensityEffectDataTable[157].fParameterDelta0   =             0;
  fDensityEffectDataTable[157].fDeltaErrorMax     =         0.096;
  fDensityEffectDataTable[157].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_FREON-12 ---------------------------------------------------------------------
  fDensityEffectDataTable[158].fName              = "NIST_MAT_FREON-12";
  fDensityEffectDataTable[158].fPlasmaEnergy      =        21.121 * eV;
  fDensityEffectDataTable[158].fSternheimerFactor =         1.974;
  fDensityEffectDataTable[158].fParameterC        =        4.8251;
  fDensityEffectDataTable[158].fParameterFitX0    =        0.3035;
  fDensityEffectDataTable[158].fParameterFitX1    =        3.2659;
  fDensityEffectDataTable[158].fParameterFitA     =       0.07978;
  fDensityEffectDataTable[158].fParameterFitM     =        3.4626;
  fDensityEffectDataTable[158].fParameterDelta0   =             0;
  fDensityEffectDataTable[158].fDeltaErrorMax     =         0.025;
  fDensityEffectDataTable[158].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_FREON-12B2 -------------------------------------------------------------------
  fDensityEffectDataTable[159].fName              = "NIST_MAT_FREON-12B2";
  fDensityEffectDataTable[159].fPlasmaEnergy      =        25.877 * eV;
  fDensityEffectDataTable[159].fSternheimerFactor =         2.195;
  fDensityEffectDataTable[159].fParameterC        =        5.7976;
  fDensityEffectDataTable[159].fParameterFitX0    =        0.3406;
  fDensityEffectDataTable[159].fParameterFitX1    =        3.7956;
  fDensityEffectDataTable[159].fParameterFitA     =       0.05144;
  fDensityEffectDataTable[159].fParameterFitM     =        3.5565;
  fDensityEffectDataTable[159].fParameterDelta0   =             0;
  fDensityEffectDataTable[159].fDeltaErrorMax     =         0.021;
  fDensityEffectDataTable[159].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_FREON-13 ---------------------------------------------------------------------
  fDensityEffectDataTable[160].fName              = "NIST_MAT_FREON-13";
  fDensityEffectDataTable[160].fPlasmaEnergy      =        19.432 * eV;
  fDensityEffectDataTable[160].fSternheimerFactor =         2.116;
  fDensityEffectDataTable[160].fParameterC        =        4.7483;
  fDensityEffectDataTable[160].fParameterFitX0    =        0.3659;
  fDensityEffectDataTable[160].fParameterFitX1    =        3.2337;
  fDensityEffectDataTable[160].fParameterFitA     =       0.07238;
  fDensityEffectDataTable[160].fParameterFitM     =        3.5551;
  fDensityEffectDataTable[160].fParameterDelta0   =             0;
  fDensityEffectDataTable[160].fDeltaErrorMax     =          0.05;
  fDensityEffectDataTable[160].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_FREON-13B1 -------------------------------------------------------------------
  fDensityEffectDataTable[161].fName              = "NIST_MAT_FREON-13B1";
  fDensityEffectDataTable[161].fPlasmaEnergy      =        23.849 * eV;
  fDensityEffectDataTable[161].fSternheimerFactor =         2.233;
  fDensityEffectDataTable[161].fParameterC        =        5.3555;
  fDensityEffectDataTable[161].fParameterFitX0    =        0.3522;
  fDensityEffectDataTable[161].fParameterFitX1    =        3.7554;
  fDensityEffectDataTable[161].fParameterFitA     =       0.03925;
  fDensityEffectDataTable[161].fParameterFitM     =        3.7194;
  fDensityEffectDataTable[161].fParameterDelta0   =             0;
  fDensityEffectDataTable[161].fDeltaErrorMax     =         0.036;
  fDensityEffectDataTable[161].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_FREON-13I1 -------------------------------------------------------------------
  fDensityEffectDataTable[162].fName              = "NIST_MAT_FREON-13I1";
  fDensityEffectDataTable[162].fPlasmaEnergy      =        25.615 * eV;
  fDensityEffectDataTable[162].fSternheimerFactor =         1.924;
  fDensityEffectDataTable[162].fParameterC        =        5.8774;
  fDensityEffectDataTable[162].fParameterFitX0    =        0.2847;
  fDensityEffectDataTable[162].fParameterFitX1    =         3.728;
  fDensityEffectDataTable[162].fParameterFitA     =       0.09112;
  fDensityEffectDataTable[162].fParameterFitM     =        3.1658;
  fDensityEffectDataTable[162].fParameterDelta0   =             0;
  fDensityEffectDataTable[162].fDeltaErrorMax     =         0.025;
  fDensityEffectDataTable[162].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GADOLINIUM_OXYSULFIDE --------------------------------------------------------
  fDensityEffectDataTable[163].fName              = "NIST_MAT_GADOLINIUM_OXYSULFIDE";
  fDensityEffectDataTable[163].fPlasmaEnergy      =        51.099 * eV;
  fDensityEffectDataTable[163].fSternheimerFactor =         2.179;
  fDensityEffectDataTable[163].fParameterC        =        5.5347;
  fDensityEffectDataTable[163].fParameterFitX0    =       -0.1774;
  fDensityEffectDataTable[163].fParameterFitX1    =        3.4045;
  fDensityEffectDataTable[163].fParameterFitA     =       0.22161;
  fDensityEffectDataTable[163].fParameterFitM     =          2.63;
  fDensityEffectDataTable[163].fParameterDelta0   =             0;
  fDensityEffectDataTable[163].fDeltaErrorMax     =         0.056;
  fDensityEffectDataTable[163].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GALLIUM_ARSENIDE -------------------------------------------------------------
  fDensityEffectDataTable[164].fName              = "NIST_MAT_GALLIUM_ARSENIDE";
  fDensityEffectDataTable[164].fPlasmaEnergy      =         44.17 * eV;
  fDensityEffectDataTable[164].fSternheimerFactor =         2.652;
  fDensityEffectDataTable[164].fParameterC        =        5.3299;
  fDensityEffectDataTable[164].fParameterFitX0    =        0.1764;
  fDensityEffectDataTable[164].fParameterFitX1    =         3.642;
  fDensityEffectDataTable[164].fParameterFitA     =       0.07152;
  fDensityEffectDataTable[164].fParameterFitM     =        3.3356;
  fDensityEffectDataTable[164].fParameterDelta0   =             0;
  fDensityEffectDataTable[164].fDeltaErrorMax     =         0.027;
  fDensityEffectDataTable[164].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GEL_PHOTO_EMULSION -----------------------------------------------------------
  fDensityEffectDataTable[165].fName              = "NIST_MAT_GEL_PHOTO_EMULSION";
  fDensityEffectDataTable[165].fPlasmaEnergy      =        24.058 * eV;
  fDensityEffectDataTable[165].fSternheimerFactor =         2.156;
  fDensityEffectDataTable[165].fParameterC        =        3.2687;
  fDensityEffectDataTable[165].fParameterFitX0    =        0.1709;
  fDensityEffectDataTable[165].fParameterFitX1    =        2.7058;
  fDensityEffectDataTable[165].fParameterFitA     =       0.10102;
  fDensityEffectDataTable[165].fParameterFitM     =        3.4418;
  fDensityEffectDataTable[165].fParameterDelta0   =             0;
  fDensityEffectDataTable[165].fDeltaErrorMax     =          0.06;
  fDensityEffectDataTable[165].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_Pyrex_Glass ------------------------------------------------------------------
  fDensityEffectDataTable[166].fName              = "NIST_MAT_Pyrex_Glass";
  fDensityEffectDataTable[166].fPlasmaEnergy      =        30.339 * eV;
  fDensityEffectDataTable[166].fSternheimerFactor =         2.369;
  fDensityEffectDataTable[166].fParameterC        =        3.9708;
  fDensityEffectDataTable[166].fParameterFitX0    =        0.1479;
  fDensityEffectDataTable[166].fParameterFitX1    =        2.9933;
  fDensityEffectDataTable[166].fParameterFitA     =        0.0827;
  fDensityEffectDataTable[166].fParameterFitM     =        3.5224;
  fDensityEffectDataTable[166].fParameterDelta0   =             0;
  fDensityEffectDataTable[166].fDeltaErrorMax     =         0.022;
  fDensityEffectDataTable[166].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GLASS_LEAD -------------------------------------------------------------------
  fDensityEffectDataTable[167].fName              = "NIST_MAT_GLASS_LEAD";
  fDensityEffectDataTable[167].fPlasmaEnergy      =        46.631 * eV;
  fDensityEffectDataTable[167].fSternheimerFactor =         2.085;
  fDensityEffectDataTable[167].fParameterC        =        5.8476;
  fDensityEffectDataTable[167].fParameterFitX0    =        0.0614;
  fDensityEffectDataTable[167].fParameterFitX1    =        3.8146;
  fDensityEffectDataTable[167].fParameterFitA     =       0.09544;
  fDensityEffectDataTable[167].fParameterFitM     =         3.074;
  fDensityEffectDataTable[167].fParameterDelta0   =             0;
  fDensityEffectDataTable[167].fDeltaErrorMax     =         0.025;
  fDensityEffectDataTable[167].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GLASS_PLATE ------------------------------------------------------------------
  fDensityEffectDataTable[168].fName              = "NIST_MAT_GLASS_PLATE";
  fDensityEffectDataTable[168].fPlasmaEnergy      =        31.481 * eV;
  fDensityEffectDataTable[168].fSternheimerFactor =         2.329;
  fDensityEffectDataTable[168].fParameterC        =        4.0602;
  fDensityEffectDataTable[168].fParameterFitX0    =        0.1237;
  fDensityEffectDataTable[168].fParameterFitX1    =        3.0649;
  fDensityEffectDataTable[168].fParameterFitA     =       0.07678;
  fDensityEffectDataTable[168].fParameterFitM     =        3.5381;
  fDensityEffectDataTable[168].fParameterDelta0   =             0;
  fDensityEffectDataTable[168].fDeltaErrorMax     =         0.025;
  fDensityEffectDataTable[168].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GLUCOSE ----------------------------------------------------------------------
  fDensityEffectDataTable[169].fName              = "NIST_MAT_GLUCOSE";
  fDensityEffectDataTable[169].fPlasmaEnergy      =        26.153 * eV;
  fDensityEffectDataTable[169].fSternheimerFactor =         2.174;
  fDensityEffectDataTable[169].fParameterC        =        3.1649;
  fDensityEffectDataTable[169].fParameterFitX0    =        0.1411;
  fDensityEffectDataTable[169].fParameterFitX1    =          2.67;
  fDensityEffectDataTable[169].fParameterFitA     =       0.10783;
  fDensityEffectDataTable[169].fParameterFitM     =        3.3946;
  fDensityEffectDataTable[169].fParameterDelta0   =             0;
  fDensityEffectDataTable[169].fDeltaErrorMax     =         0.061;
  fDensityEffectDataTable[169].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GLUTAMINE --------------------------------------------------------------------
  fDensityEffectDataTable[170].fName              = "NIST_MAT_GLUTAMINE";
  fDensityEffectDataTable[170].fPlasmaEnergy      =        25.437 * eV;
  fDensityEffectDataTable[170].fSternheimerFactor =         2.077;
  fDensityEffectDataTable[170].fParameterC        =        3.1167;
  fDensityEffectDataTable[170].fParameterFitX0    =        0.1347;
  fDensityEffectDataTable[170].fParameterFitX1    =        2.6301;
  fDensityEffectDataTable[170].fParameterFitA     =       0.11931;
  fDensityEffectDataTable[170].fParameterFitM     =        3.3254;
  fDensityEffectDataTable[170].fParameterDelta0   =             0;
  fDensityEffectDataTable[170].fDeltaErrorMax     =         0.055;
  fDensityEffectDataTable[170].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GLYCEROL ---------------------------------------------------------------------
  fDensityEffectDataTable[171].fName              = "NIST_MAT_GLYCEROL";
  fDensityEffectDataTable[171].fPlasmaEnergy      =        23.846 * eV;
  fDensityEffectDataTable[171].fSternheimerFactor =          2.12;
  fDensityEffectDataTable[171].fParameterC        =        3.2267;
  fDensityEffectDataTable[171].fParameterFitX0    =        0.1653;
  fDensityEffectDataTable[171].fParameterFitX1    =        2.6862;
  fDensityEffectDataTable[171].fParameterFitA     =       0.10168;
  fDensityEffectDataTable[171].fParameterFitM     =        3.4481;
  fDensityEffectDataTable[171].fParameterDelta0   =             0;
  fDensityEffectDataTable[171].fDeltaErrorMax     =         0.067;
  fDensityEffectDataTable[171].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GUANINE ----------------------------------------------------------------------
  fDensityEffectDataTable[172].fName              = "NIST_MAT_GUANINE";
  fDensityEffectDataTable[172].fPlasmaEnergy      =        26.022 * eV;
  fDensityEffectDataTable[172].fSternheimerFactor =          1.97;
  fDensityEffectDataTable[172].fParameterC        =        3.1171;
  fDensityEffectDataTable[172].fParameterFitX0    =        0.1163;
  fDensityEffectDataTable[172].fParameterFitX1    =        2.4296;
  fDensityEffectDataTable[172].fParameterFitA     =        0.2053;
  fDensityEffectDataTable[172].fParameterFitM     =        3.0186;
  fDensityEffectDataTable[172].fParameterDelta0   =             0;
  fDensityEffectDataTable[172].fDeltaErrorMax     =         0.069;
  fDensityEffectDataTable[172].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GYPSUM -----------------------------------------------------------------------
  fDensityEffectDataTable[173].fName              = "NIST_MAT_GYPSUM";
  fDensityEffectDataTable[173].fPlasmaEnergy      =        31.379 * eV;
  fDensityEffectDataTable[173].fSternheimerFactor =         2.187;
  fDensityEffectDataTable[173].fParameterC        =        3.8382;
  fDensityEffectDataTable[173].fParameterFitX0    =        0.0995;
  fDensityEffectDataTable[173].fParameterFitX1    =        3.1206;
  fDensityEffectDataTable[173].fParameterFitA     =       0.06949;
  fDensityEffectDataTable[173].fParameterFitM     =        3.5134;
  fDensityEffectDataTable[173].fParameterDelta0   =             0;
  fDensityEffectDataTable[173].fDeltaErrorMax     =         0.038;
  fDensityEffectDataTable[173].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_N-HEPTANE --------------------------------------------------------------------
  fDensityEffectDataTable[174].fName              = "NIST_MAT_N-HEPTANE";
  fDensityEffectDataTable[174].fPlasmaEnergy      =        18.128 * eV;
  fDensityEffectDataTable[174].fSternheimerFactor =         1.848;
  fDensityEffectDataTable[174].fParameterC        =        3.1978;
  fDensityEffectDataTable[174].fParameterFitX0    =        0.1928;
  fDensityEffectDataTable[174].fParameterFitX1    =        2.5706;
  fDensityEffectDataTable[174].fParameterFitA     =       0.11255;
  fDensityEffectDataTable[174].fParameterFitM     =        3.4885;
  fDensityEffectDataTable[174].fParameterDelta0   =             0;
  fDensityEffectDataTable[174].fDeltaErrorMax     =         0.059;
  fDensityEffectDataTable[174].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_N-HEXANE ---------------------------------------------------------------------
  fDensityEffectDataTable[175].fName              = "NIST_MAT_N-HEXANE";
  fDensityEffectDataTable[175].fPlasmaEnergy      =        17.836 * eV;
  fDensityEffectDataTable[175].fSternheimerFactor =         1.843;
  fDensityEffectDataTable[175].fParameterC        =        3.2156;
  fDensityEffectDataTable[175].fParameterFitX0    =        0.1984;
  fDensityEffectDataTable[175].fParameterFitX1    =        2.5757;
  fDensityEffectDataTable[175].fParameterFitA     =       0.11085;
  fDensityEffectDataTable[175].fParameterFitM     =        3.5027;
  fDensityEffectDataTable[175].fParameterDelta0   =             0;
  fDensityEffectDataTable[175].fDeltaErrorMax     =         0.061;
  fDensityEffectDataTable[175].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_KAPTON -----------------------------------------------------------------------
  fDensityEffectDataTable[176].fName              = "NIST_MAT_KAPTON";
  fDensityEffectDataTable[176].fPlasmaEnergy      =        24.586 * eV;
  fDensityEffectDataTable[176].fSternheimerFactor =         2.109;
  fDensityEffectDataTable[176].fParameterC        =        3.3497;
  fDensityEffectDataTable[176].fParameterFitX0    =        0.1509;
  fDensityEffectDataTable[176].fParameterFitX1    =        2.5631;
  fDensityEffectDataTable[176].fParameterFitA     =       0.15972;
  fDensityEffectDataTable[176].fParameterFitM     =        3.1921;
  fDensityEffectDataTable[176].fParameterDelta0   =             0;
  fDensityEffectDataTable[176].fDeltaErrorMax     =          0.05;
  fDensityEffectDataTable[176].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LANTHANUM_OXYBROMIDE ---------------------------------------------------------
  fDensityEffectDataTable[177].fName              = "NIST_MAT_LANTHANUM_OXYBROMIDE";
  fDensityEffectDataTable[177].fPlasmaEnergy      =        47.125 * eV;
  fDensityEffectDataTable[177].fSternheimerFactor =         1.831;
  fDensityEffectDataTable[177].fParameterC        =        5.4666;
  fDensityEffectDataTable[177].fParameterFitX0    =        -0.035;
  fDensityEffectDataTable[177].fParameterFitX1    =        3.3288;
  fDensityEffectDataTable[177].fParameterFitA     =        0.1783;
  fDensityEffectDataTable[177].fParameterFitM     =        2.8457;
  fDensityEffectDataTable[177].fParameterDelta0   =             0;
  fDensityEffectDataTable[177].fDeltaErrorMax     =          0.04;
  fDensityEffectDataTable[177].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LANTHANUM_OXYSULFIDE ---------------------------------------------------------
  fDensityEffectDataTable[178].fName              = "NIST_MAT_LANTHANUM_OXYSULFIDE";
  fDensityEffectDataTable[178].fPlasmaEnergy      =        45.394 * eV;
  fDensityEffectDataTable[178].fSternheimerFactor =         1.681;
  fDensityEffectDataTable[178].fParameterC        =        5.6151;
  fDensityEffectDataTable[178].fParameterFitX0    =       -0.0934;
  fDensityEffectDataTable[178].fParameterFitX1    =        3.2741;
  fDensityEffectDataTable[178].fParameterFitA     =       0.22579;
  fDensityEffectDataTable[178].fParameterFitM     =        2.7075;
  fDensityEffectDataTable[178].fParameterDelta0   =             0;
  fDensityEffectDataTable[178].fDeltaErrorMax     =         0.065;
  fDensityEffectDataTable[178].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LEAD_OXIDE -------------------------------------------------------------------
  fDensityEffectDataTable[179].fName              = "NIST_MAT_LEAD_OXIDE";
  fDensityEffectDataTable[179].fPlasmaEnergy      =        56.488 * eV;
  fDensityEffectDataTable[179].fSternheimerFactor =         2.012;
  fDensityEffectDataTable[179].fParameterC        =        6.2162;
  fDensityEffectDataTable[179].fParameterFitX0    =        0.0356;
  fDensityEffectDataTable[179].fParameterFitX1    =        3.5456;
  fDensityEffectDataTable[179].fParameterFitA     =       0.19645;
  fDensityEffectDataTable[179].fParameterFitM     =        2.7299;
  fDensityEffectDataTable[179].fParameterDelta0   =             0;
  fDensityEffectDataTable[179].fDeltaErrorMax     =         0.039;
  fDensityEffectDataTable[179].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LITHIUM_AMIDE ----------------------------------------------------------------
  fDensityEffectDataTable[180].fName              = "NIST_MAT_LITHIUM_AMIDE";
  fDensityEffectDataTable[180].fPlasmaEnergy      =        22.609 * eV;
  fDensityEffectDataTable[180].fSternheimerFactor =          1.74;
  fDensityEffectDataTable[180].fParameterC        =        2.7961;
  fDensityEffectDataTable[180].fParameterFitX0    =        0.0198;
  fDensityEffectDataTable[180].fParameterFitX1    =        2.5152;
  fDensityEffectDataTable[180].fParameterFitA     =        0.0874;
  fDensityEffectDataTable[180].fParameterFitM     =        3.7534;
  fDensityEffectDataTable[180].fParameterDelta0   =             0;
  fDensityEffectDataTable[180].fDeltaErrorMax     =          0.05;
  fDensityEffectDataTable[180].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LITHIUM_CARBONATE ------------------------------------------------------------
  fDensityEffectDataTable[181].fName              = "NIST_MAT_LITHIUM_CARBONATE";
  fDensityEffectDataTable[181].fPlasmaEnergy      =        29.217 * eV;
  fDensityEffectDataTable[181].fSternheimerFactor =         2.246;
  fDensityEffectDataTable[181].fParameterC        =        3.2029;
  fDensityEffectDataTable[181].fParameterFitX0    =        0.0551;
  fDensityEffectDataTable[181].fParameterFitX1    =        2.6598;
  fDensityEffectDataTable[181].fParameterFitA     =       0.09936;
  fDensityEffectDataTable[181].fParameterFitM     =        3.5417;
  fDensityEffectDataTable[181].fParameterDelta0   =             0;
  fDensityEffectDataTable[181].fDeltaErrorMax     =         0.062;
  fDensityEffectDataTable[181].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LITHIUM_FLUORIDE -------------------------------------------------------------
  fDensityEffectDataTable[182].fName              = "NIST_MAT_LITHIUM_FLUORIDE";
  fDensityEffectDataTable[182].fPlasmaEnergy      =        31.815 * eV;
  fDensityEffectDataTable[182].fSternheimerFactor =         2.197;
  fDensityEffectDataTable[182].fParameterC        =        3.1667;
  fDensityEffectDataTable[182].fParameterFitX0    =        0.0171;
  fDensityEffectDataTable[182].fParameterFitX1    =        2.7049;
  fDensityEffectDataTable[182].fParameterFitA     =       0.07593;
  fDensityEffectDataTable[182].fParameterFitM     =        3.7478;
  fDensityEffectDataTable[182].fParameterDelta0   =             0;
  fDensityEffectDataTable[182].fDeltaErrorMax     =         0.084;
  fDensityEffectDataTable[182].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LITHIUM_HYDRIDE --------------------------------------------------------------
  fDensityEffectDataTable[183].fName              = "NIST_MAT_LITHIUM_HYDRIDE";
  fDensityEffectDataTable[183].fPlasmaEnergy      =         18.51 * eV;
  fDensityEffectDataTable[183].fSternheimerFactor =         1.482;
  fDensityEffectDataTable[183].fParameterC        =         2.358;
  fDensityEffectDataTable[183].fParameterFitX0    =       -0.0988;
  fDensityEffectDataTable[183].fParameterFitX1    =        1.4515;
  fDensityEffectDataTable[183].fParameterFitA     =       0.90567;
  fDensityEffectDataTable[183].fParameterFitM     =        2.5849;
  fDensityEffectDataTable[183].fParameterDelta0   =             0;
  fDensityEffectDataTable[183].fDeltaErrorMax     =         0.035;
  fDensityEffectDataTable[183].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LITHIUM_IODIDE ---------------------------------------------------------------
  fDensityEffectDataTable[184].fName              = "NIST_MAT_LITHIUM_IODIDE";
  fDensityEffectDataTable[184].fPlasmaEnergy      =        34.841 * eV;
  fDensityEffectDataTable[184].fSternheimerFactor =         1.706;
  fDensityEffectDataTable[184].fParameterC        =        6.2671;
  fDensityEffectDataTable[184].fParameterFitX0    =        0.0892;
  fDensityEffectDataTable[184].fParameterFitX1    =        3.3702;
  fDensityEffectDataTable[184].fParameterFitA     =       0.23274;
  fDensityEffectDataTable[184].fParameterFitM     =        2.7146;
  fDensityEffectDataTable[184].fParameterDelta0   =             0;
  fDensityEffectDataTable[184].fDeltaErrorMax     =         0.043;
  fDensityEffectDataTable[184].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LITHIUM_OXIDE ----------------------------------------------------------------
  fDensityEffectDataTable[185].fName              = "NIST_MAT_LITHIUM_OXIDE";
  fDensityEffectDataTable[185].fPlasmaEnergy      =        27.984 * eV;
  fDensityEffectDataTable[185].fSternheimerFactor =         2.039;
  fDensityEffectDataTable[185].fParameterC        =         2.934;
  fDensityEffectDataTable[185].fParameterFitX0    =       -0.0511;
  fDensityEffectDataTable[185].fParameterFitX1    =        2.5874;
  fDensityEffectDataTable[185].fParameterFitA     =       0.08035;
  fDensityEffectDataTable[185].fParameterFitM     =        3.7878;
  fDensityEffectDataTable[185].fParameterDelta0   =             0;
  fDensityEffectDataTable[185].fDeltaErrorMax     =         0.043;
  fDensityEffectDataTable[185].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LITHIUM_TETRABORATE ----------------------------------------------------------
  fDensityEffectDataTable[186].fName              = "NIST_MAT_LITHIUM_TETRABORATE";
  fDensityEffectDataTable[186].fPlasmaEnergy      =        31.343 * eV;
  fDensityEffectDataTable[186].fSternheimerFactor =          2.36;
  fDensityEffectDataTable[186].fParameterC        =        3.2093;
  fDensityEffectDataTable[186].fParameterFitX0    =        0.0737;
  fDensityEffectDataTable[186].fParameterFitX1    =        2.6502;
  fDensityEffectDataTable[186].fParameterFitA     =       0.11075;
  fDensityEffectDataTable[186].fParameterFitM     =        3.4389;
  fDensityEffectDataTable[186].fParameterDelta0   =             0;
  fDensityEffectDataTable[186].fDeltaErrorMax     =         0.048;
  fDensityEffectDataTable[186].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LUNG_ICRP --------------------------------------------------------------------
  fDensityEffectDataTable[187].fName              = "NIST_MAT_LUNG_ICRP";
  fDensityEffectDataTable[187].fPlasmaEnergy      =        21.891 * eV;
  fDensityEffectDataTable[187].fSternheimerFactor =         2.184;
  fDensityEffectDataTable[187].fParameterC        =        3.4708;
  fDensityEffectDataTable[187].fParameterFitX0    =        0.2261;
  fDensityEffectDataTable[187].fParameterFitX1    =        2.8001;
  fDensityEffectDataTable[187].fParameterFitA     =       0.08588;
  fDensityEffectDataTable[187].fParameterFitM     =        3.5353;
  fDensityEffectDataTable[187].fParameterDelta0   =             0;
  fDensityEffectDataTable[187].fDeltaErrorMax     =         0.089;
  fDensityEffectDataTable[187].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_M3_WAX -----------------------------------------------------------------------
  fDensityEffectDataTable[188].fName              = "NIST_MAT_M3_WAX";
  fDensityEffectDataTable[188].fPlasmaEnergy      =            22 * eV;
  fDensityEffectDataTable[188].fSternheimerFactor =         1.975;
  fDensityEffectDataTable[188].fParameterC        =         3.254;
  fDensityEffectDataTable[188].fParameterFitX0    =        0.1523;
  fDensityEffectDataTable[188].fParameterFitX1    =        2.7529;
  fDensityEffectDataTable[188].fParameterFitA     =       0.07864;
  fDensityEffectDataTable[188].fParameterFitM     =        3.6412;
  fDensityEffectDataTable[188].fParameterDelta0   =             0;
  fDensityEffectDataTable[188].fDeltaErrorMax     =         0.044;
  fDensityEffectDataTable[188].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MAGNESIUM_CARBONATE ----------------------------------------------------------
  fDensityEffectDataTable[189].fName              = "NIST_MAT_MAGNESIUM_CARBONATE";
  fDensityEffectDataTable[189].fPlasmaEnergy      =        34.979 * eV;
  fDensityEffectDataTable[189].fSternheimerFactor =         2.388;
  fDensityEffectDataTable[189].fParameterC        =        3.4319;
  fDensityEffectDataTable[189].fParameterFitX0    =         0.086;
  fDensityEffectDataTable[189].fParameterFitX1    =        2.7997;
  fDensityEffectDataTable[189].fParameterFitA     =       0.09219;
  fDensityEffectDataTable[189].fParameterFitM     =        3.5003;
  fDensityEffectDataTable[189].fParameterDelta0   =             0;
  fDensityEffectDataTable[189].fDeltaErrorMax     =         0.045;
  fDensityEffectDataTable[189].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MAGNESIUM_FLUORIDE -----------------------------------------------------------
  fDensityEffectDataTable[190].fName              = "NIST_MAT_MAGNESIUM_FLUORIDE";
  fDensityEffectDataTable[190].fPlasmaEnergy      =        34.634 * eV;
  fDensityEffectDataTable[190].fSternheimerFactor =          2.33;
  fDensityEffectDataTable[190].fParameterC        =        3.7105;
  fDensityEffectDataTable[190].fParameterFitX0    =        0.1369;
  fDensityEffectDataTable[190].fParameterFitX1    =         2.863;
  fDensityEffectDataTable[190].fParameterFitA     =       0.07934;
  fDensityEffectDataTable[190].fParameterFitM     =        3.6485;
  fDensityEffectDataTable[190].fParameterDelta0   =             0;
  fDensityEffectDataTable[190].fDeltaErrorMax     =         0.085;
  fDensityEffectDataTable[190].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MAGNESIUM_OXIDE --------------------------------------------------------------
  fDensityEffectDataTable[191].fName              = "NIST_MAT_MAGNESIUM_OXIDE";
  fDensityEffectDataTable[191].fPlasmaEnergy      =        38.407 * eV;
  fDensityEffectDataTable[191].fSternheimerFactor =         2.412;
  fDensityEffectDataTable[191].fParameterC        =        3.6404;
  fDensityEffectDataTable[191].fParameterFitX0    =        0.0575;
  fDensityEffectDataTable[191].fParameterFitX1    =         2.858;
  fDensityEffectDataTable[191].fParameterFitA     =       0.08313;
  fDensityEffectDataTable[191].fParameterFitM     =        3.5968;
  fDensityEffectDataTable[191].fParameterDelta0   =             0;
  fDensityEffectDataTable[191].fDeltaErrorMax     =         0.055;
  fDensityEffectDataTable[191].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MAGNESIUM_TETRABORATE --------------------------------------------------------
  fDensityEffectDataTable[192].fName              = "NIST_MAT_MAGNESIUM_TETRABORATE";
  fDensityEffectDataTable[192].fPlasmaEnergy      =        32.089 * eV;
  fDensityEffectDataTable[192].fSternheimerFactor =          2.43;
  fDensityEffectDataTable[192].fParameterC        =        3.4328;
  fDensityEffectDataTable[192].fParameterFitX0    =        0.1147;
  fDensityEffectDataTable[192].fParameterFitX1    =        2.7635;
  fDensityEffectDataTable[192].fParameterFitA     =       0.09703;
  fDensityEffectDataTable[192].fParameterFitM     =        3.4893;
  fDensityEffectDataTable[192].fParameterDelta0   =             0;
  fDensityEffectDataTable[192].fDeltaErrorMax     =         0.044;
  fDensityEffectDataTable[192].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MERCURIC_IODIDE --------------------------------------------------------------
  fDensityEffectDataTable[193].fName              = "NIST_MAT_MERCURIC_IODIDE";
  fDensityEffectDataTable[193].fPlasmaEnergy      =        46.494 * eV;
  fDensityEffectDataTable[193].fSternheimerFactor =         1.892;
  fDensityEffectDataTable[193].fParameterC        =        6.3787;
  fDensityEffectDataTable[193].fParameterFitX0    =         0.104;
  fDensityEffectDataTable[193].fParameterFitX1    =        3.4728;
  fDensityEffectDataTable[193].fParameterFitA     =       0.21513;
  fDensityEffectDataTable[193].fParameterFitM     =        2.7264;
  fDensityEffectDataTable[193].fParameterDelta0   =             0;
  fDensityEffectDataTable[193].fDeltaErrorMax     =         0.047;
  fDensityEffectDataTable[193].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_METHANE ----------------------------------------------------------------------
  fDensityEffectDataTable[194].fName              = "NIST_MAT_METHANE";
  fDensityEffectDataTable[194].fPlasmaEnergy      =         0.588 * eV;
  fDensityEffectDataTable[194].fSternheimerFactor =         1.662;
  fDensityEffectDataTable[194].fParameterC        =        9.5243;
  fDensityEffectDataTable[194].fParameterFitX0    =        1.6263;
  fDensityEffectDataTable[194].fParameterFitX1    =        3.9716;
  fDensityEffectDataTable[194].fParameterFitA     =       0.09253;
  fDensityEffectDataTable[194].fParameterFitM     =        3.6257;
  fDensityEffectDataTable[194].fParameterDelta0   =             0;
  fDensityEffectDataTable[194].fDeltaErrorMax     =         0.112;
  fDensityEffectDataTable[194].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_METHANOL ---------------------------------------------------------------------
  fDensityEffectDataTable[195].fName              = "NIST_MAT_METHANOL";
  fDensityEffectDataTable[195].fPlasmaEnergy      =        19.214 * eV;
  fDensityEffectDataTable[195].fSternheimerFactor =         2.125;
  fDensityEffectDataTable[195].fParameterC        =         3.516;
  fDensityEffectDataTable[195].fParameterFitX0    =        0.2529;
  fDensityEffectDataTable[195].fParameterFitX1    =        2.7639;
  fDensityEffectDataTable[195].fParameterFitA     =        0.0897;
  fDensityEffectDataTable[195].fParameterFitM     =        3.5477;
  fDensityEffectDataTable[195].fParameterDelta0   =             0;
  fDensityEffectDataTable[195].fDeltaErrorMax     =          0.08;
  fDensityEffectDataTable[195].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MIX_D_WAX --------------------------------------------------------------------
  fDensityEffectDataTable[196].fName              = "NIST_MAT_MIX_D_WAX";
  fDensityEffectDataTable[196].fPlasmaEnergy      =        21.547 * eV;
  fDensityEffectDataTable[196].fSternheimerFactor =         1.905;
  fDensityEffectDataTable[196].fParameterC        =         3.078;
  fDensityEffectDataTable[196].fParameterFitX0    =        0.1371;
  fDensityEffectDataTable[196].fParameterFitX1    =        2.7145;
  fDensityEffectDataTable[196].fParameterFitA     =        0.0749;
  fDensityEffectDataTable[196].fParameterFitM     =        3.6823;
  fDensityEffectDataTable[196].fParameterDelta0   =             0;
  fDensityEffectDataTable[196].fDeltaErrorMax     =         0.047;
  fDensityEffectDataTable[196].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MS20_TISSUE ------------------------------------------------------------------
  fDensityEffectDataTable[197].fName              = "NIST_MAT_MS20_TISSUE";
  fDensityEffectDataTable[197].fPlasmaEnergy      =        21.153 * eV;
  fDensityEffectDataTable[197].fSternheimerFactor =          2.07;
  fDensityEffectDataTable[197].fParameterC        =        3.5241;
  fDensityEffectDataTable[197].fParameterFitX0    =        0.1997;
  fDensityEffectDataTable[197].fParameterFitX1    =        2.8033;
  fDensityEffectDataTable[197].fParameterFitA     =       0.08294;
  fDensityEffectDataTable[197].fParameterFitM     =        3.6061;
  fDensityEffectDataTable[197].fParameterDelta0   =             0;
  fDensityEffectDataTable[197].fDeltaErrorMax     =         0.053;
  fDensityEffectDataTable[197].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MUSCLE_SCELETAL_ICRP ---------------------------------------------------------
  fDensityEffectDataTable[198].fName              = "NIST_MAT_MUSCLE_SCELETAL_ICRP";
  fDensityEffectDataTable[198].fPlasmaEnergy      =        21.781 * eV;
  fDensityEffectDataTable[198].fSternheimerFactor =         2.185;
  fDensityEffectDataTable[198].fParameterC        =        3.4809;
  fDensityEffectDataTable[198].fParameterFitX0    =        0.2282;
  fDensityEffectDataTable[198].fParameterFitX1    =        2.7999;
  fDensityEffectDataTable[198].fParameterFitA     =       0.08636;
  fDensityEffectDataTable[198].fParameterFitM     =         3.533;
  fDensityEffectDataTable[198].fParameterDelta0   =             0;
  fDensityEffectDataTable[198].fDeltaErrorMax     =         0.089;
  fDensityEffectDataTable[198].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MUSCLE_STRIATED_ICRU ---------------------------------------------------------
  fDensityEffectDataTable[199].fName              = "NIST_MAT_MUSCLE_STRIATED_ICRU";
  fDensityEffectDataTable[199].fPlasmaEnergy      =        21.795 * eV;
  fDensityEffectDataTable[199].fSternheimerFactor =         2.174;
  fDensityEffectDataTable[199].fParameterC        =        3.4636;
  fDensityEffectDataTable[199].fParameterFitX0    =        0.2249;
  fDensityEffectDataTable[199].fParameterFitX1    =        2.8032;
  fDensityEffectDataTable[199].fParameterFitA     =       0.08507;
  fDensityEffectDataTable[199].fParameterFitM     =        3.5383;
  fDensityEffectDataTable[199].fParameterDelta0   =             0;
  fDensityEffectDataTable[199].fDeltaErrorMax     =         0.086;
  fDensityEffectDataTable[199].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MUSCLE_WITH_SUCROSE ----------------------------------------------------------
  fDensityEffectDataTable[200].fName              = "NIST_MAT_MUSCLE_WITH_SUCROSE";
  fDensityEffectDataTable[200].fPlasmaEnergy      =         22.48 * eV;
  fDensityEffectDataTable[200].fSternheimerFactor =         2.169;
  fDensityEffectDataTable[200].fParameterC        =         3.391;
  fDensityEffectDataTable[200].fParameterFitX0    =        0.2098;
  fDensityEffectDataTable[200].fParameterFitX1    =         2.755;
  fDensityEffectDataTable[200].fParameterFitA     =       0.09481;
  fDensityEffectDataTable[200].fParameterFitM     =        3.4699;
  fDensityEffectDataTable[200].fParameterDelta0   =             0;
  fDensityEffectDataTable[200].fDeltaErrorMax     =          0.08;
  fDensityEffectDataTable[200].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MUSCLE_WITHOUT_SUCROSE -------------------------------------------------------
  fDensityEffectDataTable[201].fName              = "NIST_MAT_MUSCLE_WITHOUT_SUCROSE";
  fDensityEffectDataTable[201].fPlasmaEnergy      =        22.109 * eV;
  fDensityEffectDataTable[201].fSternheimerFactor =         2.173;
  fDensityEffectDataTable[201].fParameterC        =        3.4216;
  fDensityEffectDataTable[201].fParameterFitX0    =        0.2187;
  fDensityEffectDataTable[201].fParameterFitX1    =         2.768;
  fDensityEffectDataTable[201].fParameterFitA     =       0.09143;
  fDensityEffectDataTable[201].fParameterFitM     =        3.4982;
  fDensityEffectDataTable[201].fParameterDelta0   =             0;
  fDensityEffectDataTable[201].fDeltaErrorMax     =         0.056;
  fDensityEffectDataTable[201].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_NAPHTALENE -------------------------------------------------------------------
  fDensityEffectDataTable[202].fName              = "NIST_MAT_NAPHTALENE";
  fDensityEffectDataTable[202].fPlasmaEnergy      =        22.459 * eV;
  fDensityEffectDataTable[202].fSternheimerFactor =         1.956;
  fDensityEffectDataTable[202].fParameterC        =        3.2274;
  fDensityEffectDataTable[202].fParameterFitX0    =        0.1374;
  fDensityEffectDataTable[202].fParameterFitX1    =        2.5429;
  fDensityEffectDataTable[202].fParameterFitA     =       0.14766;
  fDensityEffectDataTable[202].fParameterFitM     =        3.2654;
  fDensityEffectDataTable[202].fParameterDelta0   =             0;
  fDensityEffectDataTable[202].fDeltaErrorMax     =         0.051;
  fDensityEffectDataTable[202].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_NITROBENZENE -----------------------------------------------------------------
  fDensityEffectDataTable[203].fName              = "NIST_MAT_NITROBENZENE";
  fDensityEffectDataTable[203].fPlasmaEnergy      =        22.747 * eV;
  fDensityEffectDataTable[203].fSternheimerFactor =         2.065;
  fDensityEffectDataTable[203].fParameterC        =        3.4073;
  fDensityEffectDataTable[203].fParameterFitX0    =        0.1777;
  fDensityEffectDataTable[203].fParameterFitX1    =         2.663;
  fDensityEffectDataTable[203].fParameterFitA     =       0.12727;
  fDensityEffectDataTable[203].fParameterFitM     =        3.3091;
  fDensityEffectDataTable[203].fParameterDelta0   =             0;
  fDensityEffectDataTable[203].fDeltaErrorMax     =         0.051;
  fDensityEffectDataTable[203].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_NITROUS_OXIDE ----------------------------------------------------------------
  fDensityEffectDataTable[204].fName              = "NIST_MAT_NITROUS_OXIDE";
  fDensityEffectDataTable[204].fPlasmaEnergy      =         0.872 * eV;
  fDensityEffectDataTable[204].fSternheimerFactor =         2.059;
  fDensityEffectDataTable[204].fParameterC        =       10.1575;
  fDensityEffectDataTable[204].fParameterFitX0    =        1.6477;
  fDensityEffectDataTable[204].fParameterFitX1    =        4.1565;
  fDensityEffectDataTable[204].fParameterFitA     =       0.11992;
  fDensityEffectDataTable[204].fParameterFitM     =        3.3318;
  fDensityEffectDataTable[204].fParameterDelta0   =             0;
  fDensityEffectDataTable[204].fDeltaErrorMax     =         0.086;
  fDensityEffectDataTable[204].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_NYLON-8062 -------------------------------------------------------------------
  fDensityEffectDataTable[205].fName              = "NIST_MAT_NYLON-8062";
  fDensityEffectDataTable[205].fPlasmaEnergy      =        22.221 * eV;
  fDensityEffectDataTable[205].fSternheimerFactor =         1.967;
  fDensityEffectDataTable[205].fParameterC        =         3.125;
  fDensityEffectDataTable[205].fParameterFitX0    =        0.1503;
  fDensityEffectDataTable[205].fParameterFitX1    =        2.6004;
  fDensityEffectDataTable[205].fParameterFitA     =       0.11513;
  fDensityEffectDataTable[205].fParameterFitM     =        3.4044;
  fDensityEffectDataTable[205].fParameterDelta0   =             0;
  fDensityEffectDataTable[205].fDeltaErrorMax     =         0.054;
  fDensityEffectDataTable[205].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_NYLON-6-6 --------------------------------------------------------------------
  fDensityEffectDataTable[206].fName              = "NIST_MAT_NYLON-6-6";
  fDensityEffectDataTable[206].fPlasmaEnergy      =        22.774 * eV;
  fDensityEffectDataTable[206].fSternheimerFactor =         1.931;
  fDensityEffectDataTable[206].fParameterC        =        3.0634;
  fDensityEffectDataTable[206].fParameterFitX0    =        0.1336;
  fDensityEffectDataTable[206].fParameterFitX1    =        2.5834;
  fDensityEffectDataTable[206].fParameterFitA     =       0.11818;
  fDensityEffectDataTable[206].fParameterFitM     =        3.3826;
  fDensityEffectDataTable[206].fParameterDelta0   =             0;
  fDensityEffectDataTable[206].fDeltaErrorMax     =         0.051;
  fDensityEffectDataTable[206].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_NYLON-6-10 -------------------------------------------------------------------
  fDensityEffectDataTable[207].fName              = "NIST_MAT_NYLON-6-10";
  fDensityEffectDataTable[207].fPlasmaEnergy      =        22.866 * eV;
  fDensityEffectDataTable[207].fSternheimerFactor =         1.942;
  fDensityEffectDataTable[207].fParameterC        =        3.0333;
  fDensityEffectDataTable[207].fParameterFitX0    =        0.1304;
  fDensityEffectDataTable[207].fParameterFitX1    =        2.5681;
  fDensityEffectDataTable[207].fParameterFitA     =       0.11852;
  fDensityEffectDataTable[207].fParameterFitM     =        3.3912;
  fDensityEffectDataTable[207].fParameterDelta0   =             0;
  fDensityEffectDataTable[207].fDeltaErrorMax     =          0.05;
  fDensityEffectDataTable[207].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_NYLON-11_RILSAN --------------------------------------------------------------
  fDensityEffectDataTable[208].fName              = "NIST_MAT_NYLON-11_RILSAN";
  fDensityEffectDataTable[208].fPlasmaEnergy      =        25.661 * eV;
  fDensityEffectDataTable[208].fSternheimerFactor =         1.902;
  fDensityEffectDataTable[208].fParameterC        =        2.7514;
  fDensityEffectDataTable[208].fParameterFitX0    =        0.0678;
  fDensityEffectDataTable[208].fParameterFitX1    =        2.4281;
  fDensityEffectDataTable[208].fParameterFitA     =       0.14868;
  fDensityEffectDataTable[208].fParameterFitM     =        3.2576;
  fDensityEffectDataTable[208].fParameterDelta0   =             0;
  fDensityEffectDataTable[208].fDeltaErrorMax     =         0.044;
  fDensityEffectDataTable[208].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_OCTANE -----------------------------------------------------------------------
  fDensityEffectDataTable[209].fName              = "NIST_MAT_OCTANE";
  fDensityEffectDataTable[209].fPlasmaEnergy      =         18.36 * eV;
  fDensityEffectDataTable[209].fSternheimerFactor =         1.851;
  fDensityEffectDataTable[209].fParameterC        =        3.1834;
  fDensityEffectDataTable[209].fParameterFitX0    =        0.1882;
  fDensityEffectDataTable[209].fParameterFitX1    =        2.5664;
  fDensityEffectDataTable[209].fParameterFitA     =       0.11387;
  fDensityEffectDataTable[209].fParameterFitM     =        3.4776;
  fDensityEffectDataTable[209].fParameterDelta0   =             0;
  fDensityEffectDataTable[209].fDeltaErrorMax     =         0.057;
  fDensityEffectDataTable[209].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_PARAFFIN ---------------------------------------------------------------------
  fDensityEffectDataTable[210].fName              = "NIST_MAT_PARAFFIN";
  fDensityEffectDataTable[210].fPlasmaEnergy      =        21.031 * eV;
  fDensityEffectDataTable[210].fSternheimerFactor =         1.844;
  fDensityEffectDataTable[210].fParameterC        =        2.9551;
  fDensityEffectDataTable[210].fParameterFitX0    =        0.1289;
  fDensityEffectDataTable[210].fParameterFitX1    =        2.5084;
  fDensityEffectDataTable[210].fParameterFitA     =       0.12087;
  fDensityEffectDataTable[210].fParameterFitM     =        3.4288;
  fDensityEffectDataTable[210].fParameterDelta0   =             0;
  fDensityEffectDataTable[210].fDeltaErrorMax     =         0.052;
  fDensityEffectDataTable[210].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_N-PENTANE --------------------------------------------------------------------
  fDensityEffectDataTable[211].fName              = "NIST_MAT_N-PENTANE";
  fDensityEffectDataTable[211].fPlasmaEnergy      =        17.398 * eV;
  fDensityEffectDataTable[211].fSternheimerFactor =         1.842;
  fDensityEffectDataTable[211].fParameterC        =        3.2504;
  fDensityEffectDataTable[211].fParameterFitX0    =        0.2086;
  fDensityEffectDataTable[211].fParameterFitX1    =        2.5855;
  fDensityEffectDataTable[211].fParameterFitA     =       0.10809;
  fDensityEffectDataTable[211].fParameterFitM     =        3.5265;
  fDensityEffectDataTable[211].fParameterDelta0   =             0;
  fDensityEffectDataTable[211].fDeltaErrorMax     =         0.064;
  fDensityEffectDataTable[211].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_PHOTO_EMULSION ---------------------------------------------------------------
  fDensityEffectDataTable[212].fName              = "NIST_MAT_PHOTO_EMULSION";
  fDensityEffectDataTable[212].fPlasmaEnergy      =        37.946 * eV;
  fDensityEffectDataTable[212].fSternheimerFactor =         2.264;
  fDensityEffectDataTable[212].fParameterC        =        5.3319;
  fDensityEffectDataTable[212].fParameterFitX0    =        0.1009;
  fDensityEffectDataTable[212].fParameterFitX1    =        3.4866;
  fDensityEffectDataTable[212].fParameterFitA     =       0.12399;
  fDensityEffectDataTable[212].fParameterFitM     =        3.0094;
  fDensityEffectDataTable[212].fParameterDelta0   =             0;
  fDensityEffectDataTable[212].fDeltaErrorMax     =         0.028;
  fDensityEffectDataTable[212].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_PLASTIC_SC_VINYLTOLUENE ------------------------------------------------------
  fDensityEffectDataTable[213].fName              = "NIST_MAT_PLASTIC_SC_VINYLTOLUENE";
  fDensityEffectDataTable[213].fPlasmaEnergy      =         21.54 * eV;
  fDensityEffectDataTable[213].fSternheimerFactor =         1.929;
  fDensityEffectDataTable[213].fParameterC        =        3.1997;
  fDensityEffectDataTable[213].fParameterFitX0    =        0.1464;
  fDensityEffectDataTable[213].fParameterFitX1    =        2.4855;
  fDensityEffectDataTable[213].fParameterFitA     =       0.16101;
  fDensityEffectDataTable[213].fParameterFitM     =        3.2393;
  fDensityEffectDataTable[213].fParameterDelta0   =             0;
  fDensityEffectDataTable[213].fDeltaErrorMax     =          0.05;
  fDensityEffectDataTable[213].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_PLUTONIUM_DIOXIDE ------------------------------------------------------------
  fDensityEffectDataTable[214].fName              = "NIST_MAT_PLUTONIUM_DIOXIDE";
  fDensityEffectDataTable[214].fPlasmaEnergy      =        62.143 * eV;
  fDensityEffectDataTable[214].fSternheimerFactor =         1.846;
  fDensityEffectDataTable[214].fParameterC        =        5.9719;
  fDensityEffectDataTable[214].fParameterFitX0    =       -0.2311;
  fDensityEffectDataTable[214].fParameterFitX1    =        3.5554;
  fDensityEffectDataTable[214].fParameterFitA     =       0.20594;
  fDensityEffectDataTable[214].fParameterFitM     =        2.6522;
  fDensityEffectDataTable[214].fParameterDelta0   =             0;
  fDensityEffectDataTable[214].fDeltaErrorMax     =         0.111;
  fDensityEffectDataTable[214].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYACRYLONITRILE ------------------------------------------------------------
  fDensityEffectDataTable[215].fName              = "NIST_MAT_POLYACRYLONITRILE";
  fDensityEffectDataTable[215].fPlasmaEnergy      =        22.642 * eV;
  fDensityEffectDataTable[215].fSternheimerFactor =         1.955;
  fDensityEffectDataTable[215].fParameterC        =        3.2459;
  fDensityEffectDataTable[215].fParameterFitX0    =        0.1504;
  fDensityEffectDataTable[215].fParameterFitX1    =        2.5159;
  fDensityEffectDataTable[215].fParameterFitA     =       0.16275;
  fDensityEffectDataTable[215].fParameterFitM     =        3.1975;
  fDensityEffectDataTable[215].fParameterDelta0   =             0;
  fDensityEffectDataTable[215].fDeltaErrorMax     =          0.05;
  fDensityEffectDataTable[215].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYCARBONATE ----------------------------------------------------------------
  fDensityEffectDataTable[216].fName              = "NIST_MAT_POLYCARBONATE";
  fDensityEffectDataTable[216].fPlasmaEnergy      =        22.915 * eV;
  fDensityEffectDataTable[216].fSternheimerFactor =          2.06;
  fDensityEffectDataTable[216].fParameterC        =        3.3201;
  fDensityEffectDataTable[216].fParameterFitX0    =        0.1606;
  fDensityEffectDataTable[216].fParameterFitX1    =        2.6225;
  fDensityEffectDataTable[216].fParameterFitA     =        0.1286;
  fDensityEffectDataTable[216].fParameterFitM     =        3.3288;
  fDensityEffectDataTable[216].fParameterDelta0   =             0;
  fDensityEffectDataTable[216].fDeltaErrorMax     =         0.049;
  fDensityEffectDataTable[216].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYCHLOROSTYRENE ------------------------------------------------------------
  fDensityEffectDataTable[217].fName              = "NIST_MAT_POLYCHLOROSTYRENE";
  fDensityEffectDataTable[217].fPlasmaEnergy      =         23.81 * eV;
  fDensityEffectDataTable[217].fSternheimerFactor =         1.902;
  fDensityEffectDataTable[217].fParameterC        =        3.4659;
  fDensityEffectDataTable[217].fParameterFitX0    =        0.1238;
  fDensityEffectDataTable[217].fParameterFitX1    =        2.9241;
  fDensityEffectDataTable[217].fParameterFitA     =        0.0753;
  fDensityEffectDataTable[217].fParameterFitM     =        3.5441;
  fDensityEffectDataTable[217].fParameterDelta0   =             0;
  fDensityEffectDataTable[217].fDeltaErrorMax     =         0.029;
  fDensityEffectDataTable[217].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYETHYLENE -----------------------------------------------------------------
  fDensityEffectDataTable[218].fName              = "NIST_MAT_POLYETHYLENE";
  fDensityEffectDataTable[218].fPlasmaEnergy      =        21.099 * eV;
  fDensityEffectDataTable[218].fSternheimerFactor =         1.882;
  fDensityEffectDataTable[218].fParameterC        =        3.0016;
  fDensityEffectDataTable[218].fParameterFitX0    =         0.137;
  fDensityEffectDataTable[218].fParameterFitX1    =        2.5177;
  fDensityEffectDataTable[218].fParameterFitA     =       0.12108;
  fDensityEffectDataTable[218].fParameterFitM     =        3.4292;
  fDensityEffectDataTable[218].fParameterDelta0   =             0;
  fDensityEffectDataTable[218].fDeltaErrorMax     =         0.051;
  fDensityEffectDataTable[218].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_MYLAR ------------------------------------------------------------------------
  fDensityEffectDataTable[219].fName              = "NIST_MAT_MYLAR";
  fDensityEffectDataTable[219].fPlasmaEnergy      =        24.595 * eV;
  fDensityEffectDataTable[219].fSternheimerFactor =         2.144;
  fDensityEffectDataTable[219].fParameterC        =        3.3262;
  fDensityEffectDataTable[219].fParameterFitX0    =        0.1562;
  fDensityEffectDataTable[219].fParameterFitX1    =        2.6507;
  fDensityEffectDataTable[219].fParameterFitA     =       0.12679;
  fDensityEffectDataTable[219].fParameterFitM     =        3.3076;
  fDensityEffectDataTable[219].fParameterDelta0   =             0;
  fDensityEffectDataTable[219].fDeltaErrorMax     =         0.052;
  fDensityEffectDataTable[219].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_LUCITE -----------------------------------------------------------------------
  fDensityEffectDataTable[220].fName              = "NIST_MAT_LUCITE";
  fDensityEffectDataTable[220].fPlasmaEnergy      =        23.086 * eV;
  fDensityEffectDataTable[220].fSternheimerFactor =         2.173;
  fDensityEffectDataTable[220].fParameterC        =        3.3297;
  fDensityEffectDataTable[220].fParameterFitX0    =        0.1824;
  fDensityEffectDataTable[220].fParameterFitX1    =        2.6681;
  fDensityEffectDataTable[220].fParameterFitA     =       0.11433;
  fDensityEffectDataTable[220].fParameterFitM     =        3.3836;
  fDensityEffectDataTable[220].fParameterDelta0   =             0;
  fDensityEffectDataTable[220].fDeltaErrorMax     =         0.056;
  fDensityEffectDataTable[220].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYOXOMETHYLENE -------------------------------------------------------------
  fDensityEffectDataTable[221].fName              = "NIST_MAT_POLYOXOMETHYLENE";
  fDensityEffectDataTable[221].fPlasmaEnergy      =         25.11 * eV;
  fDensityEffectDataTable[221].fSternheimerFactor =         2.175;
  fDensityEffectDataTable[221].fParameterC        =        3.2514;
  fDensityEffectDataTable[221].fParameterFitX0    =        0.1584;
  fDensityEffectDataTable[221].fParameterFitX1    =        2.6838;
  fDensityEffectDataTable[221].fParameterFitA     =       0.10808;
  fDensityEffectDataTable[221].fParameterFitM     =        3.4002;
  fDensityEffectDataTable[221].fParameterDelta0   =             0;
  fDensityEffectDataTable[221].fDeltaErrorMax     =         0.063;
  fDensityEffectDataTable[221].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYPROPYLENE ----------------------------------------------------------------
  fDensityEffectDataTable[222].fName              = "NIST_MAT_POLYPROPYLENE";
  fDensityEffectDataTable[222].fPlasmaEnergy      =        20.457 * eV;
  fDensityEffectDataTable[222].fSternheimerFactor =         1.884;
  fDensityEffectDataTable[222].fParameterC        =        3.1252;
  fDensityEffectDataTable[222].fParameterFitX0    =        0.1534;
  fDensityEffectDataTable[222].fParameterFitX1    =        2.4822;
  fDensityEffectDataTable[222].fParameterFitA     =       0.15045;
  fDensityEffectDataTable[222].fParameterFitM     =        3.2855;
  fDensityEffectDataTable[222].fParameterDelta0   =             0;
  fDensityEffectDataTable[222].fDeltaErrorMax     =         0.055;
  fDensityEffectDataTable[222].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYSTYRENE ------------------------------------------------------------------
  fDensityEffectDataTable[223].fName              = "NIST_MAT_POLYSTYRENE";
  fDensityEffectDataTable[223].fPlasmaEnergy      =        21.754 * eV;
  fDensityEffectDataTable[223].fSternheimerFactor =         2.027;
  fDensityEffectDataTable[223].fParameterC        =        3.2999;
  fDensityEffectDataTable[223].fParameterFitX0    =        0.1647;
  fDensityEffectDataTable[223].fParameterFitX1    =        2.5031;
  fDensityEffectDataTable[223].fParameterFitA     =       0.16454;
  fDensityEffectDataTable[223].fParameterFitM     =        3.2224;
  fDensityEffectDataTable[223].fParameterDelta0   =             0;
  fDensityEffectDataTable[223].fDeltaErrorMax     =         0.051;
  fDensityEffectDataTable[223].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TEFLON -----------------------------------------------------------------------
  fDensityEffectDataTable[224].fName              = "NIST_MAT_TEFLON";
  fDensityEffectDataTable[224].fPlasmaEnergy      =        29.609 * eV;
  fDensityEffectDataTable[224].fSternheimerFactor =         2.142;
  fDensityEffectDataTable[224].fParameterC        =        3.4161;
  fDensityEffectDataTable[224].fParameterFitX0    =        0.1648;
  fDensityEffectDataTable[224].fParameterFitX1    =        2.7404;
  fDensityEffectDataTable[224].fParameterFitA     =       0.10606;
  fDensityEffectDataTable[224].fParameterFitM     =        3.4046;
  fDensityEffectDataTable[224].fParameterDelta0   =             0;
  fDensityEffectDataTable[224].fDeltaErrorMax     =         0.073;
  fDensityEffectDataTable[224].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYTRIFLUOROCHLOROETHYLENE --------------------------------------------------
  fDensityEffectDataTable[225].fName              = "NIST_MAT_POLYTRIFLUOROCHLOROETHYLENE";
  fDensityEffectDataTable[225].fPlasmaEnergy      =        28.955 * eV;
  fDensityEffectDataTable[225].fSternheimerFactor =         2.094;
  fDensityEffectDataTable[225].fParameterC        =        3.8551;
  fDensityEffectDataTable[225].fParameterFitX0    =        0.1714;
  fDensityEffectDataTable[225].fParameterFitX1    =        3.0265;
  fDensityEffectDataTable[225].fParameterFitA     =       0.07727;
  fDensityEffectDataTable[225].fParameterFitM     =        3.5085;
  fDensityEffectDataTable[225].fParameterDelta0   =             0;
  fDensityEffectDataTable[225].fDeltaErrorMax     =         0.035;
  fDensityEffectDataTable[225].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYVINYL_ACETATE ------------------------------------------------------------
  fDensityEffectDataTable[226].fName              = "NIST_MAT_POLYVINYL_ACETATE";
  fDensityEffectDataTable[226].fPlasmaEnergy      =        22.978 * eV;
  fDensityEffectDataTable[226].fSternheimerFactor =         2.116;
  fDensityEffectDataTable[226].fParameterC        =        3.3309;
  fDensityEffectDataTable[226].fParameterFitX0    =        0.1769;
  fDensityEffectDataTable[226].fParameterFitX1    =        2.6747;
  fDensityEffectDataTable[226].fParameterFitA     =       0.11442;
  fDensityEffectDataTable[226].fParameterFitM     =        3.3762;
  fDensityEffectDataTable[226].fParameterDelta0   =             0;
  fDensityEffectDataTable[226].fDeltaErrorMax     =         0.055;
  fDensityEffectDataTable[226].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_PLOYVINYL_ALCOHOL ------------------------------------------------------------
  fDensityEffectDataTable[227].fName              = "NIST_MAT_PLOYVINYL_ALCOHOL";
  fDensityEffectDataTable[227].fPlasmaEnergy      =        24.251 * eV;
  fDensityEffectDataTable[227].fSternheimerFactor =         2.071;
  fDensityEffectDataTable[227].fParameterC        =        3.1115;
  fDensityEffectDataTable[227].fParameterFitX0    =        0.1401;
  fDensityEffectDataTable[227].fParameterFitX1    =        2.6315;
  fDensityEffectDataTable[227].fParameterFitA     =       0.11178;
  fDensityEffectDataTable[227].fParameterFitM     =        3.3893;
  fDensityEffectDataTable[227].fParameterDelta0   =             0;
  fDensityEffectDataTable[227].fDeltaErrorMax     =         0.056;
  fDensityEffectDataTable[227].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYVINYL_BUTYRAL ------------------------------------------------------------
  fDensityEffectDataTable[228].fName              = "NIST_MAT_POLYVINYL_BUTYRAL";
  fDensityEffectDataTable[228].fPlasmaEnergy      =        22.521 * eV;
  fDensityEffectDataTable[228].fSternheimerFactor =         2.021;
  fDensityEffectDataTable[228].fParameterC        =        3.1865;
  fDensityEffectDataTable[228].fParameterFitX0    =        0.1555;
  fDensityEffectDataTable[228].fParameterFitX1    =        2.6186;
  fDensityEffectDataTable[228].fParameterFitA     =       0.11544;
  fDensityEffectDataTable[228].fParameterFitM     =        3.3983;
  fDensityEffectDataTable[228].fParameterDelta0   =             0;
  fDensityEffectDataTable[228].fDeltaErrorMax     =         0.054;
  fDensityEffectDataTable[228].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYVINYL_CHLORIDE -----------------------------------------------------------
  fDensityEffectDataTable[229].fName              = "NIST_MAT_POLYVINYL_CHLORIDE";
  fDensityEffectDataTable[229].fPlasmaEnergy      =         23.51 * eV;
  fDensityEffectDataTable[229].fSternheimerFactor =          1.84;
  fDensityEffectDataTable[229].fParameterC        =        4.0532;
  fDensityEffectDataTable[229].fParameterFitX0    =        0.1559;
  fDensityEffectDataTable[229].fParameterFitX1    =        2.9415;
  fDensityEffectDataTable[229].fParameterFitA     =       0.12438;
  fDensityEffectDataTable[229].fParameterFitM     =        3.2104;
  fDensityEffectDataTable[229].fParameterDelta0   =             0;
  fDensityEffectDataTable[229].fDeltaErrorMax     =         0.027;
  fDensityEffectDataTable[229].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYVINYLIDENE_CHLORIDE ------------------------------------------------------
  fDensityEffectDataTable[230].fName              = "NIST_MAT_POLYVINYLIDENE_CHLORIDE";
  fDensityEffectDataTable[230].fPlasmaEnergy      =        26.437 * eV;
  fDensityEffectDataTable[230].fSternheimerFactor =         1.814;
  fDensityEffectDataTable[230].fParameterC        =        4.2506;
  fDensityEffectDataTable[230].fParameterFitX0    =        0.1314;
  fDensityEffectDataTable[230].fParameterFitX1    =        2.9009;
  fDensityEffectDataTable[230].fParameterFitA     =       0.15466;
  fDensityEffectDataTable[230].fParameterFitM     =         3.102;
  fDensityEffectDataTable[230].fParameterDelta0   =             0;
  fDensityEffectDataTable[230].fDeltaErrorMax     =         0.034;
  fDensityEffectDataTable[230].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYVINYLIDENE_FLUORIDE ------------------------------------------------------
  fDensityEffectDataTable[231].fName              = "NIST_MAT_POLYVINYLIDENE_FLUORIDE";
  fDensityEffectDataTable[231].fPlasmaEnergy      =        27.024 * eV;
  fDensityEffectDataTable[231].fSternheimerFactor =          2.16;
  fDensityEffectDataTable[231].fParameterC        =        3.3793;
  fDensityEffectDataTable[231].fParameterFitX0    =        0.1717;
  fDensityEffectDataTable[231].fParameterFitX1    =        2.7375;
  fDensityEffectDataTable[231].fParameterFitA     =       0.10316;
  fDensityEffectDataTable[231].fParameterFitM     =          3.42;
  fDensityEffectDataTable[231].fParameterDelta0   =             0;
  fDensityEffectDataTable[231].fDeltaErrorMax     =         0.067;
  fDensityEffectDataTable[231].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POLYVINYL_PYRROLIDONE --------------------------------------------------------
  fDensityEffectDataTable[232].fName              = "NIST_MAT_POLYVINYL_PYRROLIDONE";
  fDensityEffectDataTable[232].fPlasmaEnergy      =        23.671 * eV;
  fDensityEffectDataTable[232].fSternheimerFactor =         1.989;
  fDensityEffectDataTable[232].fParameterC        =        3.1017;
  fDensityEffectDataTable[232].fParameterFitX0    =        0.1324;
  fDensityEffectDataTable[232].fParameterFitX1    =        2.5867;
  fDensityEffectDataTable[232].fParameterFitA     =       0.12504;
  fDensityEffectDataTable[232].fParameterFitM     =        3.3326;
  fDensityEffectDataTable[232].fParameterDelta0   =             0;
  fDensityEffectDataTable[232].fDeltaErrorMax     =         0.031;
  fDensityEffectDataTable[232].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POTASSIUM_IODIDE -------------------------------------------------------------
  fDensityEffectDataTable[233].fName              = "NIST_MAT_POTASSIUM_IODIDE";
  fDensityEffectDataTable[233].fPlasmaEnergy      =        33.575 * eV;
  fDensityEffectDataTable[233].fSternheimerFactor =         1.784;
  fDensityEffectDataTable[233].fParameterC        =        6.1088;
  fDensityEffectDataTable[233].fParameterFitX0    =        0.1044;
  fDensityEffectDataTable[233].fParameterFitX1    =        3.3442;
  fDensityEffectDataTable[233].fParameterFitA     =       0.22053;
  fDensityEffectDataTable[233].fParameterFitM     =        2.7558;
  fDensityEffectDataTable[233].fParameterDelta0   =             0;
  fDensityEffectDataTable[233].fDeltaErrorMax     =         0.042;
  fDensityEffectDataTable[233].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_POTASSIUM_OXIDE --------------------------------------------------------------
  fDensityEffectDataTable[234].fName              = "NIST_MAT_POTASSIUM_OXIDE";
  fDensityEffectDataTable[234].fPlasmaEnergy      =        30.672 * eV;
  fDensityEffectDataTable[234].fSternheimerFactor =         2.065;
  fDensityEffectDataTable[234].fParameterC        =        4.6463;
  fDensityEffectDataTable[234].fParameterFitX0    =         0.048;
  fDensityEffectDataTable[234].fParameterFitX1    =         3.011;
  fDensityEffectDataTable[234].fParameterFitA     =       0.16789;
  fDensityEffectDataTable[234].fParameterFitM     =        3.0121;
  fDensityEffectDataTable[234].fParameterDelta0   =             0;
  fDensityEffectDataTable[234].fDeltaErrorMax     =         0.027;
  fDensityEffectDataTable[234].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_PROPANE ----------------------------------------------------------------------
  fDensityEffectDataTable[235].fName              = "NIST_MAT_PROPANE";
  fDensityEffectDataTable[235].fPlasmaEnergy      =         0.959 * eV;
  fDensityEffectDataTable[235].fSternheimerFactor =         1.708;
  fDensityEffectDataTable[235].fParameterC        =        8.7878;
  fDensityEffectDataTable[235].fParameterFitX0    =        1.4326;
  fDensityEffectDataTable[235].fParameterFitX1    =        3.7998;
  fDensityEffectDataTable[235].fParameterFitA     =       0.09916;
  fDensityEffectDataTable[235].fParameterFitM     =         3.592;
  fDensityEffectDataTable[235].fParameterDelta0   =             0;
  fDensityEffectDataTable[235].fDeltaErrorMax     =         0.093;
  fDensityEffectDataTable[235].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_lPROPANE ---------------------------------------------------------------------
  fDensityEffectDataTable[236].fName              = "NIST_MAT_lPROPANE";
  fDensityEffectDataTable[236].fPlasmaEnergy      =        14.509 * eV;
  fDensityEffectDataTable[236].fSternheimerFactor =         1.844;
  fDensityEffectDataTable[236].fParameterC        =        3.5529;
  fDensityEffectDataTable[236].fParameterFitX0    =        0.2861;
  fDensityEffectDataTable[236].fParameterFitX1    =        2.6568;
  fDensityEffectDataTable[236].fParameterFitA     =       0.10329;
  fDensityEffectDataTable[236].fParameterFitM     =         3.562;
  fDensityEffectDataTable[236].fParameterDelta0   =             0;
  fDensityEffectDataTable[236].fDeltaErrorMax     =         0.068;
  fDensityEffectDataTable[236].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_ROPYL_ALCOHOL ----------------------------------------------------------------
  fDensityEffectDataTable[237].fName              = "NIST_MAT_ROPYL_ALCOHOL";
  fDensityEffectDataTable[237].fPlasmaEnergy      =        19.429 * eV;
  fDensityEffectDataTable[237].fSternheimerFactor =         1.972;
  fDensityEffectDataTable[237].fParameterC        =        3.2915;
  fDensityEffectDataTable[237].fParameterFitX0    =        0.2046;
  fDensityEffectDataTable[237].fParameterFitX1    =        2.6681;
  fDensityEffectDataTable[237].fParameterFitA     =       0.09644;
  fDensityEffectDataTable[237].fParameterFitM     =        3.5415;
  fDensityEffectDataTable[237].fParameterDelta0   =             0;
  fDensityEffectDataTable[237].fDeltaErrorMax     =          0.07;
  fDensityEffectDataTable[237].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_PYRIDINE ---------------------------------------------------------------------
  fDensityEffectDataTable[238].fName              = "NIST_MAT_PYRIDINE";
  fDensityEffectDataTable[238].fPlasmaEnergy      =        20.807 * eV;
  fDensityEffectDataTable[238].fSternheimerFactor =         1.895;
  fDensityEffectDataTable[238].fParameterC        =        3.3148;
  fDensityEffectDataTable[238].fParameterFitX0    =         0.167;
  fDensityEffectDataTable[238].fParameterFitX1    =        2.5245;
  fDensityEffectDataTable[238].fParameterFitA     =       0.16399;
  fDensityEffectDataTable[238].fParameterFitM     =        3.1977;
  fDensityEffectDataTable[238].fParameterDelta0   =             0;
  fDensityEffectDataTable[238].fDeltaErrorMax     =         0.051;
  fDensityEffectDataTable[238].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_RUBBER_BUTYL -----------------------------------------------------------------
  fDensityEffectDataTable[239].fName              = "NIST_MAT_RUBBER_BUTYL";
  fDensityEffectDataTable[239].fPlasmaEnergy      =        20.873 * eV;
  fDensityEffectDataTable[239].fSternheimerFactor =         1.852;
  fDensityEffectDataTable[239].fParameterC        =        2.9915;
  fDensityEffectDataTable[239].fParameterFitX0    =        0.1347;
  fDensityEffectDataTable[239].fParameterFitX1    =        2.5154;
  fDensityEffectDataTable[239].fParameterFitA     =       0.12108;
  fDensityEffectDataTable[239].fParameterFitM     =        3.4296;
  fDensityEffectDataTable[239].fParameterDelta0   =             0;
  fDensityEffectDataTable[239].fDeltaErrorMax     =         0.051;
  fDensityEffectDataTable[239].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_RUBBER_NATURAL ---------------------------------------------------------------
  fDensityEffectDataTable[240].fName              = "NIST_MAT_RUBBER_NATURAL";
  fDensityEffectDataTable[240].fPlasmaEnergy      =        20.644 * eV;
  fDensityEffectDataTable[240].fSternheimerFactor =         1.889;
  fDensityEffectDataTable[240].fParameterC        =        3.1272;
  fDensityEffectDataTable[240].fParameterFitX0    =        0.1512;
  fDensityEffectDataTable[240].fParameterFitX1    =        2.4815;
  fDensityEffectDataTable[240].fParameterFitA     =       0.15058;
  fDensityEffectDataTable[240].fParameterFitM     =        3.2879;
  fDensityEffectDataTable[240].fParameterDelta0   =             0;
  fDensityEffectDataTable[240].fDeltaErrorMax     =         0.053;
  fDensityEffectDataTable[240].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_RUBBER_NEOPRENE --------------------------------------------------------------
  fDensityEffectDataTable[241].fName              = "NIST_MAT_RUBBER_NEOPRENE";
  fDensityEffectDataTable[241].fPlasmaEnergy      =        23.036 * eV;
  fDensityEffectDataTable[241].fSternheimerFactor =         1.874;
  fDensityEffectDataTable[241].fParameterC        =        3.7911;
  fDensityEffectDataTable[241].fParameterFitX0    =        0.1501;
  fDensityEffectDataTable[241].fParameterFitX1    =        2.9461;
  fDensityEffectDataTable[241].fParameterFitA     =       0.09763;
  fDensityEffectDataTable[241].fParameterFitM     =        3.3632;
  fDensityEffectDataTable[241].fParameterDelta0   =             0;
  fDensityEffectDataTable[241].fDeltaErrorMax     =         0.026;
  fDensityEffectDataTable[241].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SILICON_DIOXIDE --------------------------------------------------------------
  fDensityEffectDataTable[242].fName              = "NIST_MAT_SILICON_DIOXIDE";
  fDensityEffectDataTable[242].fPlasmaEnergy      =        31.014 * eV;
  fDensityEffectDataTable[242].fSternheimerFactor =         2.335;
  fDensityEffectDataTable[242].fParameterC        =        4.0029;
  fDensityEffectDataTable[242].fParameterFitX0    =        0.1385;
  fDensityEffectDataTable[242].fParameterFitX1    =        3.0025;
  fDensityEffectDataTable[242].fParameterFitA     =       0.08408;
  fDensityEffectDataTable[242].fParameterFitM     =        3.5064;
  fDensityEffectDataTable[242].fParameterDelta0   =             0;
  fDensityEffectDataTable[242].fDeltaErrorMax     =         0.018;
  fDensityEffectDataTable[242].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SILVER_BROMIDE ---------------------------------------------------------------
  fDensityEffectDataTable[243].fName              = "NIST_MAT_SILVER_BROMIDE";
  fDensityEffectDataTable[243].fPlasmaEnergy      =        48.448 * eV;
  fDensityEffectDataTable[243].fSternheimerFactor =         2.271;
  fDensityEffectDataTable[243].fParameterC        =        5.6139;
  fDensityEffectDataTable[243].fParameterFitX0    =        0.0352;
  fDensityEffectDataTable[243].fParameterFitX1    =        3.2109;
  fDensityEffectDataTable[243].fParameterFitA     =       0.24582;
  fDensityEffectDataTable[243].fParameterFitM     =         2.682;
  fDensityEffectDataTable[243].fParameterDelta0   =             0;
  fDensityEffectDataTable[243].fDeltaErrorMax     =         0.043;
  fDensityEffectDataTable[243].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SILVER_CHLORIDE --------------------------------------------------------------
  fDensityEffectDataTable[244].fName              = "NIST_MAT_SILVER_CHLORIDE";
  fDensityEffectDataTable[244].fPlasmaEnergy      =        45.405 * eV;
  fDensityEffectDataTable[244].fSternheimerFactor =         2.096;
  fDensityEffectDataTable[244].fParameterC        =        5.3437;
  fDensityEffectDataTable[244].fParameterFitX0    =       -0.0139;
  fDensityEffectDataTable[244].fParameterFitX1    =        3.2022;
  fDensityEffectDataTable[244].fParameterFitA     =       0.22968;
  fDensityEffectDataTable[244].fParameterFitM     =        2.7041;
  fDensityEffectDataTable[244].fParameterDelta0   =             0;
  fDensityEffectDataTable[244].fDeltaErrorMax     =         0.062;
  fDensityEffectDataTable[244].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SILVER_HALIDES ---------------------------------------------------------------
  fDensityEffectDataTable[245].fName              = "NIST_MAT_SILVER_HALIDES";
  fDensityEffectDataTable[245].fPlasmaEnergy      =        48.433 * eV;
  fDensityEffectDataTable[245].fSternheimerFactor =          2.27;
  fDensityEffectDataTable[245].fParameterC        =        5.6166;
  fDensityEffectDataTable[245].fParameterFitX0    =        0.0353;
  fDensityEffectDataTable[245].fParameterFitX1    =        3.2117;
  fDensityEffectDataTable[245].fParameterFitA     =       0.24593;
  fDensityEffectDataTable[245].fParameterFitM     =        2.6814;
  fDensityEffectDataTable[245].fParameterDelta0   =             0;
  fDensityEffectDataTable[245].fDeltaErrorMax     =         0.043;
  fDensityEffectDataTable[245].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SILVER_IODIDE ----------------------------------------------------------------
  fDensityEffectDataTable[246].fName              = "NIST_MAT_SILVER_IODIDE";
  fDensityEffectDataTable[246].fPlasmaEnergy      =        46.105 * eV;
  fDensityEffectDataTable[246].fSternheimerFactor =         1.945;
  fDensityEffectDataTable[246].fParameterC        =        5.9342;
  fDensityEffectDataTable[246].fParameterFitX0    =        0.0148;
  fDensityEffectDataTable[246].fParameterFitX1    =        3.2908;
  fDensityEffectDataTable[246].fParameterFitA     =       0.25059;
  fDensityEffectDataTable[246].fParameterFitM     =        2.6572;
  fDensityEffectDataTable[246].fParameterDelta0   =             0;
  fDensityEffectDataTable[246].fDeltaErrorMax     =         0.071;
  fDensityEffectDataTable[246].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SKIN_ICRP --------------------------------------------------------------------
  fDensityEffectDataTable[247].fName              = "NIST_MAT_SKIN_ICRP";
  fDensityEffectDataTable[247].fPlasmaEnergy      =          22.4 * eV;
  fDensityEffectDataTable[247].fSternheimerFactor =          2.14;
  fDensityEffectDataTable[247].fParameterC        =        3.3546;
  fDensityEffectDataTable[247].fParameterFitX0    =        0.2019;
  fDensityEffectDataTable[247].fParameterFitX1    =        2.7526;
  fDensityEffectDataTable[247].fParameterFitA     =       0.09459;
  fDensityEffectDataTable[247].fParameterFitM     =        3.4643;
  fDensityEffectDataTable[247].fParameterDelta0   =             0;
  fDensityEffectDataTable[247].fDeltaErrorMax     =         0.076;
  fDensityEffectDataTable[247].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SODIUM_CARBONATE -------------------------------------------------------------
  fDensityEffectDataTable[248].fName              = "NIST_MAT_SODIUM_CARBONATE";
  fDensityEffectDataTable[248].fPlasmaEnergy      =        32.117 * eV;
  fDensityEffectDataTable[248].fSternheimerFactor =         2.557;
  fDensityEffectDataTable[248].fParameterC        =        3.7178;
  fDensityEffectDataTable[248].fParameterFitX0    =        0.1287;
  fDensityEffectDataTable[248].fParameterFitX1    =        2.8591;
  fDensityEffectDataTable[248].fParameterFitA     =       0.08715;
  fDensityEffectDataTable[248].fParameterFitM     =        3.5638;
  fDensityEffectDataTable[248].fParameterDelta0   =             0;
  fDensityEffectDataTable[248].fDeltaErrorMax     =         0.074;
  fDensityEffectDataTable[248].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SODIUM_IODIDE ----------------------------------------------------------------
  fDensityEffectDataTable[249].fName              = "NIST_MAT_SODIUM_IODIDE";
  fDensityEffectDataTable[249].fPlasmaEnergy      =        36.057 * eV;
  fDensityEffectDataTable[249].fSternheimerFactor =         1.857;
  fDensityEffectDataTable[249].fParameterC        =        6.0572;
  fDensityEffectDataTable[249].fParameterFitX0    =        0.1203;
  fDensityEffectDataTable[249].fParameterFitX1    =         3.592;
  fDensityEffectDataTable[249].fParameterFitA     =       0.12516;
  fDensityEffectDataTable[249].fParameterFitM     =        3.0398;
  fDensityEffectDataTable[249].fParameterDelta0   =             0;
  fDensityEffectDataTable[249].fDeltaErrorMax     =         0.031;
  fDensityEffectDataTable[249].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SODIUM_MONOXIDE --------------------------------------------------------------
  fDensityEffectDataTable[250].fName              = "NIST_MAT_SODIUM_MONOXIDE";
  fDensityEffectDataTable[250].fPlasmaEnergy      =        30.205 * eV;
  fDensityEffectDataTable[250].fSternheimerFactor =         2.689;
  fDensityEffectDataTable[250].fParameterC        =        4.1892;
  fDensityEffectDataTable[250].fParameterFitX0    =        0.1652;
  fDensityEffectDataTable[250].fParameterFitX1    =        2.9793;
  fDensityEffectDataTable[250].fParameterFitA     =       0.07501;
  fDensityEffectDataTable[250].fParameterFitM     =        3.6943;
  fDensityEffectDataTable[250].fParameterDelta0   =             0;
  fDensityEffectDataTable[250].fDeltaErrorMax     =         0.097;
  fDensityEffectDataTable[250].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SODIUM_NITRATE ---------------------------------------------------------------
  fDensityEffectDataTable[251].fName              = "NIST_MAT_SODIUM_NITRATE";
  fDensityEffectDataTable[251].fPlasmaEnergy      =        30.459 * eV;
  fDensityEffectDataTable[251].fSternheimerFactor =         2.456;
  fDensityEffectDataTable[251].fParameterC        =        3.6502;
  fDensityEffectDataTable[251].fParameterFitX0    =        0.1534;
  fDensityEffectDataTable[251].fParameterFitX1    =        2.8221;
  fDensityEffectDataTable[251].fParameterFitA     =       0.09391;
  fDensityEffectDataTable[251].fParameterFitM     =        3.5097;
  fDensityEffectDataTable[251].fParameterDelta0   =             0;
  fDensityEffectDataTable[251].fDeltaErrorMax     =         0.081;
  fDensityEffectDataTable[251].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_STILBENE ---------------------------------------------------------------------
  fDensityEffectDataTable[252].fName              = "NIST_MAT_STILBENE";
  fDensityEffectDataTable[252].fPlasmaEnergy      =        20.719 * eV;
  fDensityEffectDataTable[252].fSternheimerFactor =         1.963;
  fDensityEffectDataTable[252].fParameterC        =         3.368;
  fDensityEffectDataTable[252].fParameterFitX0    =        0.1734;
  fDensityEffectDataTable[252].fParameterFitX1    =        2.5142;
  fDensityEffectDataTable[252].fParameterFitA     =       0.16659;
  fDensityEffectDataTable[252].fParameterFitM     =        3.2168;
  fDensityEffectDataTable[252].fParameterDelta0   =             0;
  fDensityEffectDataTable[252].fDeltaErrorMax     =         0.052;
  fDensityEffectDataTable[252].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_SUCROSE ----------------------------------------------------------------------
  fDensityEffectDataTable[253].fName              = "NIST_MAT_SUCROSE";
  fDensityEffectDataTable[253].fPlasmaEnergy      =        26.416 * eV;
  fDensityEffectDataTable[253].fSternheimerFactor =         2.167;
  fDensityEffectDataTable[253].fParameterC        =        3.1526;
  fDensityEffectDataTable[253].fParameterFitX0    =        0.1341;
  fDensityEffectDataTable[253].fParameterFitX1    =        2.6558;
  fDensityEffectDataTable[253].fParameterFitA     =       0.11301;
  fDensityEffectDataTable[253].fParameterFitM     =         3.363;
  fDensityEffectDataTable[253].fParameterDelta0   =             0;
  fDensityEffectDataTable[253].fDeltaErrorMax     =         0.057;
  fDensityEffectDataTable[253].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TERPHENYL --------------------------------------------------------------------
  fDensityEffectDataTable[254].fName              = "NIST_MAT_TERPHENYL";
  fDensityEffectDataTable[254].fPlasmaEnergy      =        23.116 * eV;
  fDensityEffectDataTable[254].fSternheimerFactor =         1.976;
  fDensityEffectDataTable[254].fParameterC        =        3.2639;
  fDensityEffectDataTable[254].fParameterFitX0    =        0.1322;
  fDensityEffectDataTable[254].fParameterFitX1    =        2.5429;
  fDensityEffectDataTable[254].fParameterFitA     =       0.14964;
  fDensityEffectDataTable[254].fParameterFitM     =        3.2685;
  fDensityEffectDataTable[254].fParameterDelta0   =             0;
  fDensityEffectDataTable[254].fDeltaErrorMax     =         0.043;
  fDensityEffectDataTable[254].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TESTES_ICRP ------------------------------------------------------------------
  fDensityEffectDataTable[255].fName              = "NIST_MAT_TESTES_ICRP";
  fDensityEffectDataTable[255].fPlasmaEnergy      =        21.815 * eV;
  fDensityEffectDataTable[255].fSternheimerFactor =         2.185;
  fDensityEffectDataTable[255].fParameterC        =        3.4698;
  fDensityEffectDataTable[255].fParameterFitX0    =        0.2274;
  fDensityEffectDataTable[255].fParameterFitX1    =        2.7988;
  fDensityEffectDataTable[255].fParameterFitA     =       0.08533;
  fDensityEffectDataTable[255].fParameterFitM     =        3.5428;
  fDensityEffectDataTable[255].fParameterDelta0   =             0;
  fDensityEffectDataTable[255].fDeltaErrorMax     =         0.091;
  fDensityEffectDataTable[255].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TETRACHLOROETHYLENE ----------------------------------------------------------
  fDensityEffectDataTable[256].fName              = "NIST_MAT_TETRACHLOROETHYLENE";
  fDensityEffectDataTable[256].fPlasmaEnergy      =        25.513 * eV;
  fDensityEffectDataTable[256].fSternheimerFactor =          1.79;
  fDensityEffectDataTable[256].fParameterC        =        4.6619;
  fDensityEffectDataTable[256].fParameterFitX0    =        0.1713;
  fDensityEffectDataTable[256].fParameterFitX1    =        2.9083;
  fDensityEffectDataTable[256].fParameterFitA     =       0.18595;
  fDensityEffectDataTable[256].fParameterFitM     =        3.0156;
  fDensityEffectDataTable[256].fParameterDelta0   =             0;
  fDensityEffectDataTable[256].fDeltaErrorMax     =         0.038;
  fDensityEffectDataTable[256].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_THALIUM_CHLORIDE -------------------------------------------------------------
  fDensityEffectDataTable[257].fName              = "NIST_MAT_THALIUM_CHLORIDE";
  fDensityEffectDataTable[257].fPlasmaEnergy      =        48.749 * eV;
  fDensityEffectDataTable[257].fSternheimerFactor =         1.997;
  fDensityEffectDataTable[257].fParameterC        =        6.3009;
  fDensityEffectDataTable[257].fParameterFitX0    =        0.0705;
  fDensityEffectDataTable[257].fParameterFitX1    =        3.5716;
  fDensityEffectDataTable[257].fParameterFitA     =       0.18599;
  fDensityEffectDataTable[257].fParameterFitM     =         2.769;
  fDensityEffectDataTable[257].fParameterDelta0   =             0;
  fDensityEffectDataTable[257].fDeltaErrorMax     =          0.04;
  fDensityEffectDataTable[257].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TISSUE_SOFT_ICRP -------------------------------------------------------------
  fDensityEffectDataTable[258].fName              = "NIST_MAT_TISSUE_SOFT_ICRP";
  fDensityEffectDataTable[258].fPlasmaEnergy      =        21.394 * eV;
  fDensityEffectDataTable[258].fSternheimerFactor =         2.144;
  fDensityEffectDataTable[258].fParameterC        =        3.4354;
  fDensityEffectDataTable[258].fParameterFitX0    =        0.2211;
  fDensityEffectDataTable[258].fParameterFitX1    =        2.7799;
  fDensityEffectDataTable[258].fParameterFitA     =       0.08926;
  fDensityEffectDataTable[258].fParameterFitM     =         3.511;
  fDensityEffectDataTable[258].fParameterDelta0   =             0;
  fDensityEffectDataTable[258].fDeltaErrorMax     =         0.077;
  fDensityEffectDataTable[258].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TISSUE_SOFT_ICRU-4 -----------------------------------------------------------
  fDensityEffectDataTable[259].fName              = "NIST_MAT_TISSUE_SOFT_ICRU-4";
  fDensityEffectDataTable[259].fPlasmaEnergy      =        21.366 * eV;
  fDensityEffectDataTable[259].fSternheimerFactor =         2.192;
  fDensityEffectDataTable[259].fParameterC        =        3.5087;
  fDensityEffectDataTable[259].fParameterFitX0    =        0.2377;
  fDensityEffectDataTable[259].fParameterFitX1    =        2.7908;
  fDensityEffectDataTable[259].fParameterFitA     =       0.09629;
  fDensityEffectDataTable[259].fParameterFitM     =        3.4371;
  fDensityEffectDataTable[259].fParameterDelta0   =             0;
  fDensityEffectDataTable[259].fDeltaErrorMax     =         0.092;
  fDensityEffectDataTable[259].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TISSUE-METHANE ---------------------------------------------------------------
  fDensityEffectDataTable[260].fName              = "NIST_MAT_TISSUE-METHANE";
  fDensityEffectDataTable[260].fPlasmaEnergy      =         0.697 * eV;
  fDensityEffectDataTable[260].fSternheimerFactor =          1.89;
  fDensityEffectDataTable[260].fParameterC        =          9.95;
  fDensityEffectDataTable[260].fParameterFitX0    =        1.6442;
  fDensityEffectDataTable[260].fParameterFitX1    =        4.1399;
  fDensityEffectDataTable[260].fParameterFitA     =       0.09946;
  fDensityEffectDataTable[260].fParameterFitM     =        3.4708;
  fDensityEffectDataTable[260].fParameterDelta0   =             0;
  fDensityEffectDataTable[260].fDeltaErrorMax     =         0.098;
  fDensityEffectDataTable[260].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TISSUE-PROPANE ---------------------------------------------------------------
  fDensityEffectDataTable[261].fName              = "NIST_MAT_TISSUE-PROPANE";
  fDensityEffectDataTable[261].fPlasmaEnergy      =         0.913 * eV;
  fDensityEffectDataTable[261].fSternheimerFactor =         1.856;
  fDensityEffectDataTable[261].fParameterC        =        9.3529;
  fDensityEffectDataTable[261].fParameterFitX0    =        1.5139;
  fDensityEffectDataTable[261].fParameterFitX1    =        3.9916;
  fDensityEffectDataTable[261].fParameterFitA     =       0.09802;
  fDensityEffectDataTable[261].fParameterFitM     =        3.5159;
  fDensityEffectDataTable[261].fParameterDelta0   =             0;
  fDensityEffectDataTable[261].fDeltaErrorMax     =         0.092;
  fDensityEffectDataTable[261].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TITANIUM_DIOXIDE -------------------------------------------------------------
  fDensityEffectDataTable[262].fName              = "NIST_MAT_TITANIUM_DIOXIDE";
  fDensityEffectDataTable[262].fPlasmaEnergy      =        41.022 * eV;
  fDensityEffectDataTable[262].fSternheimerFactor =         2.307;
  fDensityEffectDataTable[262].fParameterC        =        3.9522;
  fDensityEffectDataTable[262].fParameterFitX0    =       -0.0119;
  fDensityEffectDataTable[262].fParameterFitX1    =        3.1647;
  fDensityEffectDataTable[262].fParameterFitA     =       0.08569;
  fDensityEffectDataTable[262].fParameterFitM     =        3.3267;
  fDensityEffectDataTable[262].fParameterDelta0   =             0;
  fDensityEffectDataTable[262].fDeltaErrorMax     =         0.027;
  fDensityEffectDataTable[262].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TOLUENE ----------------------------------------------------------------------
  fDensityEffectDataTable[263].fName              = "NIST_MAT_TOLUENE";
  fDensityEffectDataTable[263].fPlasmaEnergy      =        19.764 * eV;
  fDensityEffectDataTable[263].fSternheimerFactor =          1.88;
  fDensityEffectDataTable[263].fParameterC        =        3.3026;
  fDensityEffectDataTable[263].fParameterFitX0    =        0.1722;
  fDensityEffectDataTable[263].fParameterFitX1    =        2.5728;
  fDensityEffectDataTable[263].fParameterFitA     =       0.13284;
  fDensityEffectDataTable[263].fParameterFitM     =        3.3558;
  fDensityEffectDataTable[263].fParameterDelta0   =             0;
  fDensityEffectDataTable[263].fDeltaErrorMax     =         0.052;
  fDensityEffectDataTable[263].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TRICHLOROETHYLENE ------------------------------------------------------------
  fDensityEffectDataTable[264].fName              = "NIST_MAT_TRICHLOROETHYLENE";
  fDensityEffectDataTable[264].fPlasmaEnergy      =        24.301 * eV;
  fDensityEffectDataTable[264].fSternheimerFactor =         1.789;
  fDensityEffectDataTable[264].fParameterC        =        4.6148;
  fDensityEffectDataTable[264].fParameterFitX0    =        0.1803;
  fDensityEffectDataTable[264].fParameterFitX1    =         2.914;
  fDensityEffectDataTable[264].fParameterFitA     =       0.18272;
  fDensityEffectDataTable[264].fParameterFitM     =        3.0137;
  fDensityEffectDataTable[264].fParameterDelta0   =             0;
  fDensityEffectDataTable[264].fDeltaErrorMax     =         0.036;
  fDensityEffectDataTable[264].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TRIETHYL_PHOSPHATE -----------------------------------------------------------
  fDensityEffectDataTable[265].fName              = "NIST_MAT_TRIETHYL_PHOSPHATE";
  fDensityEffectDataTable[265].fPlasmaEnergy      =        21.863 * eV;
  fDensityEffectDataTable[265].fSternheimerFactor =           2.1;
  fDensityEffectDataTable[265].fParameterC        =        3.6242;
  fDensityEffectDataTable[265].fParameterFitX0    =        0.2054;
  fDensityEffectDataTable[265].fParameterFitX1    =        2.9428;
  fDensityEffectDataTable[265].fParameterFitA     =       0.06922;
  fDensityEffectDataTable[265].fParameterFitM     =        3.6302;
  fDensityEffectDataTable[265].fParameterDelta0   =             0;
  fDensityEffectDataTable[265].fDeltaErrorMax     =         0.049;
  fDensityEffectDataTable[265].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_TUNGSTEN_HEXAFLUORIDE --------------------------------------------------------
  fDensityEffectDataTable[266].fName              = "NIST_MAT_TUNGSTEN_HEXAFLUORIDE";
  fDensityEffectDataTable[266].fPlasmaEnergy      =        29.265 * eV;
  fDensityEffectDataTable[266].fSternheimerFactor =         2.325;
  fDensityEffectDataTable[266].fParameterC        =        5.9881;
  fDensityEffectDataTable[266].fParameterFitX0    =         0.302;
  fDensityEffectDataTable[266].fParameterFitX1    =        4.2602;
  fDensityEffectDataTable[266].fParameterFitA     =       0.03658;
  fDensityEffectDataTable[266].fParameterFitM     =        3.5134;
  fDensityEffectDataTable[266].fParameterDelta0   =             0;
  fDensityEffectDataTable[266].fDeltaErrorMax     =         0.055;
  fDensityEffectDataTable[266].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_URANIUM_DICARBIDE ------------------------------------------------------------
  fDensityEffectDataTable[267].fName              = "NIST_MAT_URANIUM_DICARBIDE";
  fDensityEffectDataTable[267].fPlasmaEnergy      =        60.969 * eV;
  fDensityEffectDataTable[267].fSternheimerFactor =         1.703;
  fDensityEffectDataTable[267].fParameterC        =        6.0247;
  fDensityEffectDataTable[267].fParameterFitX0    =       -0.2191;
  fDensityEffectDataTable[267].fParameterFitX1    =        3.5208;
  fDensityEffectDataTable[267].fParameterFitA     =        0.2112;
  fDensityEffectDataTable[267].fParameterFitM     =        2.6577;
  fDensityEffectDataTable[267].fParameterDelta0   =             0;
  fDensityEffectDataTable[267].fDeltaErrorMax     =          0.12;
  fDensityEffectDataTable[267].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_URANIUM_MONOCARBIDE ----------------------------------------------------------
  fDensityEffectDataTable[268].fName              = "NIST_MAT_URANIUM_MONOCARBIDE";
  fDensityEffectDataTable[268].fPlasmaEnergy      =        66.602 * eV;
  fDensityEffectDataTable[268].fSternheimerFactor =          1.68;
  fDensityEffectDataTable[268].fParameterC        =         6.121;
  fDensityEffectDataTable[268].fParameterFitX0    =       -0.2524;
  fDensityEffectDataTable[268].fParameterFitX1    =        3.4941;
  fDensityEffectDataTable[268].fParameterFitA     =       0.22972;
  fDensityEffectDataTable[268].fParameterFitM     =        2.6169;
  fDensityEffectDataTable[268].fParameterDelta0   =             0;
  fDensityEffectDataTable[268].fDeltaErrorMax     =         0.132;
  fDensityEffectDataTable[268].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_URANIUM_OXIDE ----------------------------------------------------------------
  fDensityEffectDataTable[269].fName              = "NIST_MAT_URANIUM_OXIDE";
  fDensityEffectDataTable[269].fPlasmaEnergy      =        60.332 * eV;
  fDensityEffectDataTable[269].fSternheimerFactor =          1.76;
  fDensityEffectDataTable[269].fParameterC        =        5.9605;
  fDensityEffectDataTable[269].fParameterFitX0    =       -0.1938;
  fDensityEffectDataTable[269].fParameterFitX1    =        3.5292;
  fDensityEffectDataTable[269].fParameterFitA     =       0.20463;
  fDensityEffectDataTable[269].fParameterFitM     =        2.6711;
  fDensityEffectDataTable[269].fParameterDelta0   =             0;
  fDensityEffectDataTable[269].fDeltaErrorMax     =         0.098;
  fDensityEffectDataTable[269].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_UREA -------------------------------------------------------------------------
  fDensityEffectDataTable[270].fName              = "NIST_MAT_UREA";
  fDensityEffectDataTable[270].fPlasmaEnergy      =        24.194 * eV;
  fDensityEffectDataTable[270].fSternheimerFactor =         2.022;
  fDensityEffectDataTable[270].fParameterC        =        3.2032;
  fDensityEffectDataTable[270].fParameterFitX0    =        0.1603;
  fDensityEffectDataTable[270].fParameterFitX1    =        2.6225;
  fDensityEffectDataTable[270].fParameterFitA     =       0.11609;
  fDensityEffectDataTable[270].fParameterFitM     =        3.3461;
  fDensityEffectDataTable[270].fParameterDelta0   =             0;
  fDensityEffectDataTable[270].fDeltaErrorMax     =          0.06;
  fDensityEffectDataTable[270].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_VALINE -----------------------------------------------------------------------
  fDensityEffectDataTable[271].fName              = "NIST_MAT_VALINE";
  fDensityEffectDataTable[271].fPlasmaEnergy      =        23.622 * eV;
  fDensityEffectDataTable[271].fSternheimerFactor =         2.024;
  fDensityEffectDataTable[271].fParameterC        =        3.1059;
  fDensityEffectDataTable[271].fParameterFitX0    =        0.1441;
  fDensityEffectDataTable[271].fParameterFitX1    =        2.6227;
  fDensityEffectDataTable[271].fParameterFitA     =       0.11386;
  fDensityEffectDataTable[271].fParameterFitM     =        3.3774;
  fDensityEffectDataTable[271].fParameterDelta0   =             0;
  fDensityEffectDataTable[271].fDeltaErrorMax     =         0.056;
  fDensityEffectDataTable[271].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_VITON ------------------------------------------------------------------------
  fDensityEffectDataTable[272].fName              = "NIST_MAT_VITON";
  fDensityEffectDataTable[272].fPlasmaEnergy      =        26.948 * eV;
  fDensityEffectDataTable[272].fSternheimerFactor =         2.227;
  fDensityEffectDataTable[272].fParameterC        =        3.5943;
  fDensityEffectDataTable[272].fParameterFitX0    =        0.2106;
  fDensityEffectDataTable[272].fParameterFitX1    =        2.7874;
  fDensityEffectDataTable[272].fParameterFitA     =       0.09965;
  fDensityEffectDataTable[272].fParameterFitM     =        3.4556;
  fDensityEffectDataTable[272].fParameterDelta0   =             0;
  fDensityEffectDataTable[272].fDeltaErrorMax     =          0.07;
  fDensityEffectDataTable[272].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_WATER ------------------------------------------------------------------------
  fDensityEffectDataTable[273].fName              = "NIST_MAT_WATER";
  fDensityEffectDataTable[273].fPlasmaEnergy      =        21.469 * eV;
  fDensityEffectDataTable[273].fSternheimerFactor =         2.203;
  fDensityEffectDataTable[273].fParameterC        =        3.5017;
  fDensityEffectDataTable[273].fParameterFitX0    =          0.24;
  fDensityEffectDataTable[273].fParameterFitX1    =        2.8004;
  fDensityEffectDataTable[273].fParameterFitA     =       0.09116;
  fDensityEffectDataTable[273].fParameterFitM     =        3.4773;
  fDensityEffectDataTable[273].fParameterDelta0   =             0;
  fDensityEffectDataTable[273].fDeltaErrorMax     =         0.097;
  fDensityEffectDataTable[273].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_WATER_VAPOR ------------------------------------------------------------------
  fDensityEffectDataTable[274].fName              = "NIST_MAT_WATER_VAPOR";
  fDensityEffectDataTable[274].fPlasmaEnergy      =          0.59 * eV;
  fDensityEffectDataTable[274].fSternheimerFactor =         2.175;
  fDensityEffectDataTable[274].fParameterC        =       10.5962;
  fDensityEffectDataTable[274].fParameterFitX0    =        1.7952;
  fDensityEffectDataTable[274].fParameterFitX1    =        4.3437;
  fDensityEffectDataTable[274].fParameterFitA     =       0.08101;
  fDensityEffectDataTable[274].fParameterFitM     =        3.5901;
  fDensityEffectDataTable[274].fParameterDelta0   =             0;
  fDensityEffectDataTable[274].fDeltaErrorMax     =         0.121;
  fDensityEffectDataTable[274].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_XYLENE -----------------------------------------------------------------------
  fDensityEffectDataTable[275].fName              = "NIST_MAT_XYLENE";
  fDensityEffectDataTable[275].fPlasmaEnergy      =        19.866 * eV;
  fDensityEffectDataTable[275].fSternheimerFactor =         1.882;
  fDensityEffectDataTable[275].fParameterC        =        3.2698;
  fDensityEffectDataTable[275].fParameterFitX0    =        0.1695;
  fDensityEffectDataTable[275].fParameterFitX1    =        2.5675;
  fDensityEffectDataTable[275].fParameterFitA     =       0.13216;
  fDensityEffectDataTable[275].fParameterFitM     =        3.3564;
  fDensityEffectDataTable[275].fParameterDelta0   =             0;
  fDensityEffectDataTable[275].fDeltaErrorMax     =         0.051;
  fDensityEffectDataTable[275].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GRAPHITE ---------------------------------------------------------------------
  fDensityEffectDataTable[276].fName              = "NIST_MAT_GRAPHITE";
  fDensityEffectDataTable[276].fPlasmaEnergy      =        30.652 * eV;
  fDensityEffectDataTable[276].fSternheimerFactor =          2.29;
  fDensityEffectDataTable[276].fParameterC        =         2.868;
  fDensityEffectDataTable[276].fParameterFitX0    =       -0.0178;
  fDensityEffectDataTable[276].fParameterFitX1    =        2.3415;
  fDensityEffectDataTable[276].fParameterFitA     =       0.26142;
  fDensityEffectDataTable[276].fParameterFitM     =        2.8697;
  fDensityEffectDataTable[276].fParameterDelta0   =          0.12;
  fDensityEffectDataTable[276].fDeltaErrorMax     =         0.038;
  fDensityEffectDataTable[276].fState             = MaterialState::kStateUndefined;


  // NIST_MAT_GRAPHITE_POROUS --------------------------------------------------------------
  fDensityEffectDataTable[277].fName              = "NIST_MAT_GRAPHITE_POROUS";
  fDensityEffectDataTable[277].fPlasmaEnergy      =        26.555 * eV;
  fDensityEffectDataTable[277].fSternheimerFactor =          2.49;
  fDensityEffectDataTable[277].fParameterC        =         3.155;
  fDensityEffectDataTable[277].fParameterFitX0    =         0.048;
  fDensityEffectDataTable[277].fParameterFitX1    =        2.5387;
  fDensityEffectDataTable[277].fParameterFitA     =       0.20762;
  fDensityEffectDataTable[277].fParameterFitM     =        2.9532;
  fDensityEffectDataTable[277].fParameterDelta0   =          0.14;
  fDensityEffectDataTable[277].fDeltaErrorMax     =         0.038;
  fDensityEffectDataTable[277].fState             = MaterialState::kStateUndefined;
}

}   // namespace geantphysics
