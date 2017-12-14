
#include "NISTMaterialData.h"
#include "PhysicalConstants.h"
#include "MaterialState.h"

namespace geantphysics {


void NISTMaterialData::BuildTable() {
    //
    // Pre-defined Material Properties for 316 Materials:
    //
    // Elemental materials (for Z = 1-98) ==> index [0-97]
    // Compound Materials        ==> index [98-277]
    // HEP and Nuclear Materials ==> index [278-293]
    // Space Science Materials   ==> index [294-296]
    // Biochemical Materials     ==> index [297-315]

    using geant::g;
    using geant::cm3;
    using geant::eV;
    using geant::kelvin;
    using geant::pascal;
    using geant::kUniverseMeanDensity;
    using geant::kNTPTemperature;
    using geant::kSTPPressure;

    // =======================================================================================
    // Elemental materials Z = 1-98 ==> index 0-97
    // =======================================================================================

    // NIST_MAT_H ----------------------------------------------------------------------------
    static const int ElemsZ_0[]        = {  1};
    static const double FractionsZ_0[] = {  1};
    fNISTMaterialDataTable[  0].fName                 = "NIST_MAT_H";
    fNISTMaterialDataTable[  0].fDensity              = 8.3748e-05 * (g/cm3);
    fNISTMaterialDataTable[  0].fMeanExcitationEnergy = 19.2 * eV;
    fNISTMaterialDataTable[  0].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[  0].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[  0].fNumComponents        = 1;
    fNISTMaterialDataTable[  0].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[  0].fElementList          = ElemsZ_0;
    fNISTMaterialDataTable[  0].fElementFraction      = FractionsZ_0;
    fNISTMaterialDataTable[  0].fIsBuiltByAtomCount   = true;


    // NIST_MAT_He ---------------------------------------------------------------------------
    static const int ElemsZ_1[]        = {  2};
    static const double FractionsZ_1[] = {  1};
    fNISTMaterialDataTable[  1].fName                 = "NIST_MAT_He";
    fNISTMaterialDataTable[  1].fDensity              = 0.000166322 * (g/cm3);
    fNISTMaterialDataTable[  1].fMeanExcitationEnergy = 41.8 * eV;
    fNISTMaterialDataTable[  1].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[  1].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[  1].fNumComponents        = 1;
    fNISTMaterialDataTable[  1].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[  1].fElementList          = ElemsZ_1;
    fNISTMaterialDataTable[  1].fElementFraction      = FractionsZ_1;
    fNISTMaterialDataTable[  1].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Li ---------------------------------------------------------------------------
    static const int ElemsZ_2[]        = {  3};
    static const double FractionsZ_2[] = {  1};
    fNISTMaterialDataTable[  2].fName                 = "NIST_MAT_Li";
    fNISTMaterialDataTable[  2].fDensity              = 0.534 * (g/cm3);
    fNISTMaterialDataTable[  2].fMeanExcitationEnergy = 40 * eV;
    fNISTMaterialDataTable[  2].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[  2].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[  2].fNumComponents        = 1;
    fNISTMaterialDataTable[  2].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[  2].fElementList          = ElemsZ_2;
    fNISTMaterialDataTable[  2].fElementFraction      = FractionsZ_2;
    fNISTMaterialDataTable[  2].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Be ---------------------------------------------------------------------------
    static const int ElemsZ_3[]        = {  4};
    static const double FractionsZ_3[] = {  1};
    fNISTMaterialDataTable[  3].fName                 = "NIST_MAT_Be";
    fNISTMaterialDataTable[  3].fDensity              = 1.848 * (g/cm3);
    fNISTMaterialDataTable[  3].fMeanExcitationEnergy = 63.7 * eV;
    fNISTMaterialDataTable[  3].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[  3].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[  3].fNumComponents        = 1;
    fNISTMaterialDataTable[  3].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[  3].fElementList          = ElemsZ_3;
    fNISTMaterialDataTable[  3].fElementFraction      = FractionsZ_3;
    fNISTMaterialDataTable[  3].fIsBuiltByAtomCount   = true;


    // NIST_MAT_B ----------------------------------------------------------------------------
    static const int ElemsZ_4[]        = {  5};
    static const double FractionsZ_4[] = {  1};
    fNISTMaterialDataTable[  4].fName                 = "NIST_MAT_B";
    fNISTMaterialDataTable[  4].fDensity              = 2.37 * (g/cm3);
    fNISTMaterialDataTable[  4].fMeanExcitationEnergy = 76 * eV;
    fNISTMaterialDataTable[  4].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[  4].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[  4].fNumComponents        = 1;
    fNISTMaterialDataTable[  4].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[  4].fElementList          = ElemsZ_4;
    fNISTMaterialDataTable[  4].fElementFraction      = FractionsZ_4;
    fNISTMaterialDataTable[  4].fIsBuiltByAtomCount   = true;


    // NIST_MAT_C ----------------------------------------------------------------------------
    static const int ElemsZ_5[]        = {  6};
    static const double FractionsZ_5[] = {  1};
    fNISTMaterialDataTable[  5].fName                 = "NIST_MAT_C";
    fNISTMaterialDataTable[  5].fDensity              = 2 * (g/cm3);
    fNISTMaterialDataTable[  5].fMeanExcitationEnergy = 81 * eV;
    fNISTMaterialDataTable[  5].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[  5].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[  5].fNumComponents        = 1;
    fNISTMaterialDataTable[  5].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[  5].fElementList          = ElemsZ_5;
    fNISTMaterialDataTable[  5].fElementFraction      = FractionsZ_5;
    fNISTMaterialDataTable[  5].fIsBuiltByAtomCount   = true;


    // NIST_MAT_N ----------------------------------------------------------------------------
    static const int ElemsZ_6[]        = {  7};
    static const double FractionsZ_6[] = {  1};
    fNISTMaterialDataTable[  6].fName                 = "NIST_MAT_N";
    fNISTMaterialDataTable[  6].fDensity              = 0.0011652 * (g/cm3);
    fNISTMaterialDataTable[  6].fMeanExcitationEnergy = 82 * eV;
    fNISTMaterialDataTable[  6].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[  6].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[  6].fNumComponents        = 1;
    fNISTMaterialDataTable[  6].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[  6].fElementList          = ElemsZ_6;
    fNISTMaterialDataTable[  6].fElementFraction      = FractionsZ_6;
    fNISTMaterialDataTable[  6].fIsBuiltByAtomCount   = true;


    // NIST_MAT_O ----------------------------------------------------------------------------
    static const int ElemsZ_7[]        = {  8};
    static const double FractionsZ_7[] = {  1};
    fNISTMaterialDataTable[  7].fName                 = "NIST_MAT_O";
    fNISTMaterialDataTable[  7].fDensity              = 0.00133151 * (g/cm3);
    fNISTMaterialDataTable[  7].fMeanExcitationEnergy = 95 * eV;
    fNISTMaterialDataTable[  7].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[  7].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[  7].fNumComponents        = 1;
    fNISTMaterialDataTable[  7].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[  7].fElementList          = ElemsZ_7;
    fNISTMaterialDataTable[  7].fElementFraction      = FractionsZ_7;
    fNISTMaterialDataTable[  7].fIsBuiltByAtomCount   = true;


    // NIST_MAT_F ----------------------------------------------------------------------------
    static const int ElemsZ_8[]        = {  9};
    static const double FractionsZ_8[] = {  1};
    fNISTMaterialDataTable[  8].fName                 = "NIST_MAT_F";
    fNISTMaterialDataTable[  8].fDensity              = 0.00158029 * (g/cm3);
    fNISTMaterialDataTable[  8].fMeanExcitationEnergy = 115 * eV;
    fNISTMaterialDataTable[  8].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[  8].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[  8].fNumComponents        = 1;
    fNISTMaterialDataTable[  8].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[  8].fElementList          = ElemsZ_8;
    fNISTMaterialDataTable[  8].fElementFraction      = FractionsZ_8;
    fNISTMaterialDataTable[  8].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ne ---------------------------------------------------------------------------
    static const int ElemsZ_9[]        = { 10};
    static const double FractionsZ_9[] = {  1};
    fNISTMaterialDataTable[  9].fName                 = "NIST_MAT_Ne";
    fNISTMaterialDataTable[  9].fDensity              = 0.000838505 * (g/cm3);
    fNISTMaterialDataTable[  9].fMeanExcitationEnergy = 137 * eV;
    fNISTMaterialDataTable[  9].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[  9].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[  9].fNumComponents        = 1;
    fNISTMaterialDataTable[  9].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[  9].fElementList          = ElemsZ_9;
    fNISTMaterialDataTable[  9].fElementFraction      = FractionsZ_9;
    fNISTMaterialDataTable[  9].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Na ---------------------------------------------------------------------------
    static const int ElemsZ_10[]        = { 11};
    static const double FractionsZ_10[] = {  1};
    fNISTMaterialDataTable[ 10].fName                 = "NIST_MAT_Na";
    fNISTMaterialDataTable[ 10].fDensity              = 0.971 * (g/cm3);
    fNISTMaterialDataTable[ 10].fMeanExcitationEnergy = 149 * eV;
    fNISTMaterialDataTable[ 10].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 10].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 10].fNumComponents        = 1;
    fNISTMaterialDataTable[ 10].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 10].fElementList          = ElemsZ_10;
    fNISTMaterialDataTable[ 10].fElementFraction      = FractionsZ_10;
    fNISTMaterialDataTable[ 10].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Mg ---------------------------------------------------------------------------
    static const int ElemsZ_11[]        = { 12};
    static const double FractionsZ_11[] = {  1};
    fNISTMaterialDataTable[ 11].fName                 = "NIST_MAT_Mg";
    fNISTMaterialDataTable[ 11].fDensity              = 1.74 * (g/cm3);
    fNISTMaterialDataTable[ 11].fMeanExcitationEnergy = 156 * eV;
    fNISTMaterialDataTable[ 11].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 11].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 11].fNumComponents        = 1;
    fNISTMaterialDataTable[ 11].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 11].fElementList          = ElemsZ_11;
    fNISTMaterialDataTable[ 11].fElementFraction      = FractionsZ_11;
    fNISTMaterialDataTable[ 11].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Al ---------------------------------------------------------------------------
    static const int ElemsZ_12[]        = { 13};
    static const double FractionsZ_12[] = {  1};
    fNISTMaterialDataTable[ 12].fName                 = "NIST_MAT_Al";
    fNISTMaterialDataTable[ 12].fDensity              = 2.699 * (g/cm3);
    fNISTMaterialDataTable[ 12].fMeanExcitationEnergy = 166 * eV;
    fNISTMaterialDataTable[ 12].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 12].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 12].fNumComponents        = 1;
    fNISTMaterialDataTable[ 12].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 12].fElementList          = ElemsZ_12;
    fNISTMaterialDataTable[ 12].fElementFraction      = FractionsZ_12;
    fNISTMaterialDataTable[ 12].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Si ---------------------------------------------------------------------------
    static const int ElemsZ_13[]        = { 14};
    static const double FractionsZ_13[] = {  1};
    fNISTMaterialDataTable[ 13].fName                 = "NIST_MAT_Si";
    fNISTMaterialDataTable[ 13].fDensity              = 2.33 * (g/cm3);
    fNISTMaterialDataTable[ 13].fMeanExcitationEnergy = 173 * eV;
    fNISTMaterialDataTable[ 13].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 13].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 13].fNumComponents        = 1;
    fNISTMaterialDataTable[ 13].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 13].fElementList          = ElemsZ_13;
    fNISTMaterialDataTable[ 13].fElementFraction      = FractionsZ_13;
    fNISTMaterialDataTable[ 13].fIsBuiltByAtomCount   = true;


    // NIST_MAT_P ----------------------------------------------------------------------------
    static const int ElemsZ_14[]        = { 15};
    static const double FractionsZ_14[] = {  1};
    fNISTMaterialDataTable[ 14].fName                 = "NIST_MAT_P";
    fNISTMaterialDataTable[ 14].fDensity              = 2.2 * (g/cm3);
    fNISTMaterialDataTable[ 14].fMeanExcitationEnergy = 173 * eV;
    fNISTMaterialDataTable[ 14].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 14].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 14].fNumComponents        = 1;
    fNISTMaterialDataTable[ 14].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 14].fElementList          = ElemsZ_14;
    fNISTMaterialDataTable[ 14].fElementFraction      = FractionsZ_14;
    fNISTMaterialDataTable[ 14].fIsBuiltByAtomCount   = true;


    // NIST_MAT_S ----------------------------------------------------------------------------
    static const int ElemsZ_15[]        = { 16};
    static const double FractionsZ_15[] = {  1};
    fNISTMaterialDataTable[ 15].fName                 = "NIST_MAT_S";
    fNISTMaterialDataTable[ 15].fDensity              = 2 * (g/cm3);
    fNISTMaterialDataTable[ 15].fMeanExcitationEnergy = 180 * eV;
    fNISTMaterialDataTable[ 15].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 15].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 15].fNumComponents        = 1;
    fNISTMaterialDataTable[ 15].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 15].fElementList          = ElemsZ_15;
    fNISTMaterialDataTable[ 15].fElementFraction      = FractionsZ_15;
    fNISTMaterialDataTable[ 15].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Cl ---------------------------------------------------------------------------
    static const int ElemsZ_16[]        = { 17};
    static const double FractionsZ_16[] = {  1};
    fNISTMaterialDataTable[ 16].fName                 = "NIST_MAT_Cl";
    fNISTMaterialDataTable[ 16].fDensity              = 0.00299473 * (g/cm3);
    fNISTMaterialDataTable[ 16].fMeanExcitationEnergy = 174 * eV;
    fNISTMaterialDataTable[ 16].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 16].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 16].fNumComponents        = 1;
    fNISTMaterialDataTable[ 16].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[ 16].fElementList          = ElemsZ_16;
    fNISTMaterialDataTable[ 16].fElementFraction      = FractionsZ_16;
    fNISTMaterialDataTable[ 16].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ar ---------------------------------------------------------------------------
    static const int ElemsZ_17[]        = { 18};
    static const double FractionsZ_17[] = {  1};
    fNISTMaterialDataTable[ 17].fName                 = "NIST_MAT_Ar";
    fNISTMaterialDataTable[ 17].fDensity              = 0.00166201 * (g/cm3);
    fNISTMaterialDataTable[ 17].fMeanExcitationEnergy = 188 * eV;
    fNISTMaterialDataTable[ 17].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 17].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 17].fNumComponents        = 1;
    fNISTMaterialDataTable[ 17].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[ 17].fElementList          = ElemsZ_17;
    fNISTMaterialDataTable[ 17].fElementFraction      = FractionsZ_17;
    fNISTMaterialDataTable[ 17].fIsBuiltByAtomCount   = true;


    // NIST_MAT_K ----------------------------------------------------------------------------
    static const int ElemsZ_18[]        = { 19};
    static const double FractionsZ_18[] = {  1};
    fNISTMaterialDataTable[ 18].fName                 = "NIST_MAT_K";
    fNISTMaterialDataTable[ 18].fDensity              = 0.862 * (g/cm3);
    fNISTMaterialDataTable[ 18].fMeanExcitationEnergy = 190 * eV;
    fNISTMaterialDataTable[ 18].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 18].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 18].fNumComponents        = 1;
    fNISTMaterialDataTable[ 18].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 18].fElementList          = ElemsZ_18;
    fNISTMaterialDataTable[ 18].fElementFraction      = FractionsZ_18;
    fNISTMaterialDataTable[ 18].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ca ---------------------------------------------------------------------------
    static const int ElemsZ_19[]        = { 20};
    static const double FractionsZ_19[] = {  1};
    fNISTMaterialDataTable[ 19].fName                 = "NIST_MAT_Ca";
    fNISTMaterialDataTable[ 19].fDensity              = 1.55 * (g/cm3);
    fNISTMaterialDataTable[ 19].fMeanExcitationEnergy = 191 * eV;
    fNISTMaterialDataTable[ 19].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 19].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 19].fNumComponents        = 1;
    fNISTMaterialDataTable[ 19].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 19].fElementList          = ElemsZ_19;
    fNISTMaterialDataTable[ 19].fElementFraction      = FractionsZ_19;
    fNISTMaterialDataTable[ 19].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Sc ---------------------------------------------------------------------------
    static const int ElemsZ_20[]        = { 21};
    static const double FractionsZ_20[] = {  1};
    fNISTMaterialDataTable[ 20].fName                 = "NIST_MAT_Sc";
    fNISTMaterialDataTable[ 20].fDensity              = 2.989 * (g/cm3);
    fNISTMaterialDataTable[ 20].fMeanExcitationEnergy = 216 * eV;
    fNISTMaterialDataTable[ 20].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 20].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 20].fNumComponents        = 1;
    fNISTMaterialDataTable[ 20].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 20].fElementList          = ElemsZ_20;
    fNISTMaterialDataTable[ 20].fElementFraction      = FractionsZ_20;
    fNISTMaterialDataTable[ 20].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ti ---------------------------------------------------------------------------
    static const int ElemsZ_21[]        = { 22};
    static const double FractionsZ_21[] = {  1};
    fNISTMaterialDataTable[ 21].fName                 = "NIST_MAT_Ti";
    fNISTMaterialDataTable[ 21].fDensity              = 4.54 * (g/cm3);
    fNISTMaterialDataTable[ 21].fMeanExcitationEnergy = 233 * eV;
    fNISTMaterialDataTable[ 21].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 21].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 21].fNumComponents        = 1;
    fNISTMaterialDataTable[ 21].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 21].fElementList          = ElemsZ_21;
    fNISTMaterialDataTable[ 21].fElementFraction      = FractionsZ_21;
    fNISTMaterialDataTable[ 21].fIsBuiltByAtomCount   = true;


    // NIST_MAT_V ----------------------------------------------------------------------------
    static const int ElemsZ_22[]        = { 23};
    static const double FractionsZ_22[] = {  1};
    fNISTMaterialDataTable[ 22].fName                 = "NIST_MAT_V";
    fNISTMaterialDataTable[ 22].fDensity              = 6.11 * (g/cm3);
    fNISTMaterialDataTable[ 22].fMeanExcitationEnergy = 245 * eV;
    fNISTMaterialDataTable[ 22].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 22].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 22].fNumComponents        = 1;
    fNISTMaterialDataTable[ 22].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 22].fElementList          = ElemsZ_22;
    fNISTMaterialDataTable[ 22].fElementFraction      = FractionsZ_22;
    fNISTMaterialDataTable[ 22].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Cr ---------------------------------------------------------------------------
    static const int ElemsZ_23[]        = { 24};
    static const double FractionsZ_23[] = {  1};
    fNISTMaterialDataTable[ 23].fName                 = "NIST_MAT_Cr";
    fNISTMaterialDataTable[ 23].fDensity              = 7.18 * (g/cm3);
    fNISTMaterialDataTable[ 23].fMeanExcitationEnergy = 257 * eV;
    fNISTMaterialDataTable[ 23].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 23].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 23].fNumComponents        = 1;
    fNISTMaterialDataTable[ 23].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 23].fElementList          = ElemsZ_23;
    fNISTMaterialDataTable[ 23].fElementFraction      = FractionsZ_23;
    fNISTMaterialDataTable[ 23].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Mn ---------------------------------------------------------------------------
    static const int ElemsZ_24[]        = { 25};
    static const double FractionsZ_24[] = {  1};
    fNISTMaterialDataTable[ 24].fName                 = "NIST_MAT_Mn";
    fNISTMaterialDataTable[ 24].fDensity              = 7.44 * (g/cm3);
    fNISTMaterialDataTable[ 24].fMeanExcitationEnergy = 272 * eV;
    fNISTMaterialDataTable[ 24].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 24].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 24].fNumComponents        = 1;
    fNISTMaterialDataTable[ 24].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 24].fElementList          = ElemsZ_24;
    fNISTMaterialDataTable[ 24].fElementFraction      = FractionsZ_24;
    fNISTMaterialDataTable[ 24].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Fe ---------------------------------------------------------------------------
    static const int ElemsZ_25[]        = { 26};
    static const double FractionsZ_25[] = {  1};
    fNISTMaterialDataTable[ 25].fName                 = "NIST_MAT_Fe";
    fNISTMaterialDataTable[ 25].fDensity              = 7.874 * (g/cm3);
    fNISTMaterialDataTable[ 25].fMeanExcitationEnergy = 286 * eV;
    fNISTMaterialDataTable[ 25].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 25].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 25].fNumComponents        = 1;
    fNISTMaterialDataTable[ 25].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 25].fElementList          = ElemsZ_25;
    fNISTMaterialDataTable[ 25].fElementFraction      = FractionsZ_25;
    fNISTMaterialDataTable[ 25].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Co ---------------------------------------------------------------------------
    static const int ElemsZ_26[]        = { 27};
    static const double FractionsZ_26[] = {  1};
    fNISTMaterialDataTable[ 26].fName                 = "NIST_MAT_Co";
    fNISTMaterialDataTable[ 26].fDensity              = 8.9 * (g/cm3);
    fNISTMaterialDataTable[ 26].fMeanExcitationEnergy = 297 * eV;
    fNISTMaterialDataTable[ 26].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 26].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 26].fNumComponents        = 1;
    fNISTMaterialDataTable[ 26].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 26].fElementList          = ElemsZ_26;
    fNISTMaterialDataTable[ 26].fElementFraction      = FractionsZ_26;
    fNISTMaterialDataTable[ 26].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ni ---------------------------------------------------------------------------
    static const int ElemsZ_27[]        = { 28};
    static const double FractionsZ_27[] = {  1};
    fNISTMaterialDataTable[ 27].fName                 = "NIST_MAT_Ni";
    fNISTMaterialDataTable[ 27].fDensity              = 8.902 * (g/cm3);
    fNISTMaterialDataTable[ 27].fMeanExcitationEnergy = 311 * eV;
    fNISTMaterialDataTable[ 27].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 27].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 27].fNumComponents        = 1;
    fNISTMaterialDataTable[ 27].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 27].fElementList          = ElemsZ_27;
    fNISTMaterialDataTable[ 27].fElementFraction      = FractionsZ_27;
    fNISTMaterialDataTable[ 27].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Cu ---------------------------------------------------------------------------
    static const int ElemsZ_28[]        = { 29};
    static const double FractionsZ_28[] = {  1};
    fNISTMaterialDataTable[ 28].fName                 = "NIST_MAT_Cu";
    fNISTMaterialDataTable[ 28].fDensity              = 8.96 * (g/cm3);
    fNISTMaterialDataTable[ 28].fMeanExcitationEnergy = 322 * eV;
    fNISTMaterialDataTable[ 28].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 28].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 28].fNumComponents        = 1;
    fNISTMaterialDataTable[ 28].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 28].fElementList          = ElemsZ_28;
    fNISTMaterialDataTable[ 28].fElementFraction      = FractionsZ_28;
    fNISTMaterialDataTable[ 28].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Zn ---------------------------------------------------------------------------
    static const int ElemsZ_29[]        = { 30};
    static const double FractionsZ_29[] = {  1};
    fNISTMaterialDataTable[ 29].fName                 = "NIST_MAT_Zn";
    fNISTMaterialDataTable[ 29].fDensity              = 7.133 * (g/cm3);
    fNISTMaterialDataTable[ 29].fMeanExcitationEnergy = 330 * eV;
    fNISTMaterialDataTable[ 29].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 29].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 29].fNumComponents        = 1;
    fNISTMaterialDataTable[ 29].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 29].fElementList          = ElemsZ_29;
    fNISTMaterialDataTable[ 29].fElementFraction      = FractionsZ_29;
    fNISTMaterialDataTable[ 29].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ga ---------------------------------------------------------------------------
    static const int ElemsZ_30[]        = { 31};
    static const double FractionsZ_30[] = {  1};
    fNISTMaterialDataTable[ 30].fName                 = "NIST_MAT_Ga";
    fNISTMaterialDataTable[ 30].fDensity              = 5.904 * (g/cm3);
    fNISTMaterialDataTable[ 30].fMeanExcitationEnergy = 334 * eV;
    fNISTMaterialDataTable[ 30].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 30].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 30].fNumComponents        = 1;
    fNISTMaterialDataTable[ 30].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 30].fElementList          = ElemsZ_30;
    fNISTMaterialDataTable[ 30].fElementFraction      = FractionsZ_30;
    fNISTMaterialDataTable[ 30].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ge ---------------------------------------------------------------------------
    static const int ElemsZ_31[]        = { 32};
    static const double FractionsZ_31[] = {  1};
    fNISTMaterialDataTable[ 31].fName                 = "NIST_MAT_Ge";
    fNISTMaterialDataTable[ 31].fDensity              = 5.323 * (g/cm3);
    fNISTMaterialDataTable[ 31].fMeanExcitationEnergy = 350 * eV;
    fNISTMaterialDataTable[ 31].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 31].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 31].fNumComponents        = 1;
    fNISTMaterialDataTable[ 31].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 31].fElementList          = ElemsZ_31;
    fNISTMaterialDataTable[ 31].fElementFraction      = FractionsZ_31;
    fNISTMaterialDataTable[ 31].fIsBuiltByAtomCount   = true;


    // NIST_MAT_As ---------------------------------------------------------------------------
    static const int ElemsZ_32[]        = { 33};
    static const double FractionsZ_32[] = {  1};
    fNISTMaterialDataTable[ 32].fName                 = "NIST_MAT_As";
    fNISTMaterialDataTable[ 32].fDensity              = 5.73 * (g/cm3);
    fNISTMaterialDataTable[ 32].fMeanExcitationEnergy = 347 * eV;
    fNISTMaterialDataTable[ 32].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 32].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 32].fNumComponents        = 1;
    fNISTMaterialDataTable[ 32].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 32].fElementList          = ElemsZ_32;
    fNISTMaterialDataTable[ 32].fElementFraction      = FractionsZ_32;
    fNISTMaterialDataTable[ 32].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Se ---------------------------------------------------------------------------
    static const int ElemsZ_33[]        = { 34};
    static const double FractionsZ_33[] = {  1};
    fNISTMaterialDataTable[ 33].fName                 = "NIST_MAT_Se";
    fNISTMaterialDataTable[ 33].fDensity              = 4.5 * (g/cm3);
    fNISTMaterialDataTable[ 33].fMeanExcitationEnergy = 348 * eV;
    fNISTMaterialDataTable[ 33].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 33].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 33].fNumComponents        = 1;
    fNISTMaterialDataTable[ 33].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 33].fElementList          = ElemsZ_33;
    fNISTMaterialDataTable[ 33].fElementFraction      = FractionsZ_33;
    fNISTMaterialDataTable[ 33].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Br ---------------------------------------------------------------------------
    static const int ElemsZ_34[]        = { 35};
    static const double FractionsZ_34[] = {  1};
    fNISTMaterialDataTable[ 34].fName                 = "NIST_MAT_Br";
    fNISTMaterialDataTable[ 34].fDensity              = 0.0070721 * (g/cm3);
    fNISTMaterialDataTable[ 34].fMeanExcitationEnergy = 343 * eV;
    fNISTMaterialDataTable[ 34].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 34].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 34].fNumComponents        = 1;
    fNISTMaterialDataTable[ 34].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[ 34].fElementList          = ElemsZ_34;
    fNISTMaterialDataTable[ 34].fElementFraction      = FractionsZ_34;
    fNISTMaterialDataTable[ 34].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Kr ---------------------------------------------------------------------------
    static const int ElemsZ_35[]        = { 36};
    static const double FractionsZ_35[] = {  1};
    fNISTMaterialDataTable[ 35].fName                 = "NIST_MAT_Kr";
    fNISTMaterialDataTable[ 35].fDensity              = 0.00347832 * (g/cm3);
    fNISTMaterialDataTable[ 35].fMeanExcitationEnergy = 352 * eV;
    fNISTMaterialDataTable[ 35].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 35].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 35].fNumComponents        = 1;
    fNISTMaterialDataTable[ 35].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[ 35].fElementList          = ElemsZ_35;
    fNISTMaterialDataTable[ 35].fElementFraction      = FractionsZ_35;
    fNISTMaterialDataTable[ 35].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Rb ---------------------------------------------------------------------------
    static const int ElemsZ_36[]        = { 37};
    static const double FractionsZ_36[] = {  1};
    fNISTMaterialDataTable[ 36].fName                 = "NIST_MAT_Rb";
    fNISTMaterialDataTable[ 36].fDensity              = 1.532 * (g/cm3);
    fNISTMaterialDataTable[ 36].fMeanExcitationEnergy = 363 * eV;
    fNISTMaterialDataTable[ 36].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 36].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 36].fNumComponents        = 1;
    fNISTMaterialDataTable[ 36].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 36].fElementList          = ElemsZ_36;
    fNISTMaterialDataTable[ 36].fElementFraction      = FractionsZ_36;
    fNISTMaterialDataTable[ 36].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Sr ---------------------------------------------------------------------------
    static const int ElemsZ_37[]        = { 38};
    static const double FractionsZ_37[] = {  1};
    fNISTMaterialDataTable[ 37].fName                 = "NIST_MAT_Sr";
    fNISTMaterialDataTable[ 37].fDensity              = 2.54 * (g/cm3);
    fNISTMaterialDataTable[ 37].fMeanExcitationEnergy = 366 * eV;
    fNISTMaterialDataTable[ 37].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 37].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 37].fNumComponents        = 1;
    fNISTMaterialDataTable[ 37].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 37].fElementList          = ElemsZ_37;
    fNISTMaterialDataTable[ 37].fElementFraction      = FractionsZ_37;
    fNISTMaterialDataTable[ 37].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Y ----------------------------------------------------------------------------
    static const int ElemsZ_38[]        = { 39};
    static const double FractionsZ_38[] = {  1};
    fNISTMaterialDataTable[ 38].fName                 = "NIST_MAT_Y";
    fNISTMaterialDataTable[ 38].fDensity              = 4.469 * (g/cm3);
    fNISTMaterialDataTable[ 38].fMeanExcitationEnergy = 379 * eV;
    fNISTMaterialDataTable[ 38].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 38].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 38].fNumComponents        = 1;
    fNISTMaterialDataTable[ 38].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 38].fElementList          = ElemsZ_38;
    fNISTMaterialDataTable[ 38].fElementFraction      = FractionsZ_38;
    fNISTMaterialDataTable[ 38].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Zr ---------------------------------------------------------------------------
    static const int ElemsZ_39[]        = { 40};
    static const double FractionsZ_39[] = {  1};
    fNISTMaterialDataTable[ 39].fName                 = "NIST_MAT_Zr";
    fNISTMaterialDataTable[ 39].fDensity              = 6.506 * (g/cm3);
    fNISTMaterialDataTable[ 39].fMeanExcitationEnergy = 393 * eV;
    fNISTMaterialDataTable[ 39].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 39].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 39].fNumComponents        = 1;
    fNISTMaterialDataTable[ 39].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 39].fElementList          = ElemsZ_39;
    fNISTMaterialDataTable[ 39].fElementFraction      = FractionsZ_39;
    fNISTMaterialDataTable[ 39].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Nb ---------------------------------------------------------------------------
    static const int ElemsZ_40[]        = { 41};
    static const double FractionsZ_40[] = {  1};
    fNISTMaterialDataTable[ 40].fName                 = "NIST_MAT_Nb";
    fNISTMaterialDataTable[ 40].fDensity              = 8.57 * (g/cm3);
    fNISTMaterialDataTable[ 40].fMeanExcitationEnergy = 417 * eV;
    fNISTMaterialDataTable[ 40].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 40].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 40].fNumComponents        = 1;
    fNISTMaterialDataTable[ 40].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 40].fElementList          = ElemsZ_40;
    fNISTMaterialDataTable[ 40].fElementFraction      = FractionsZ_40;
    fNISTMaterialDataTable[ 40].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Mo ---------------------------------------------------------------------------
    static const int ElemsZ_41[]        = { 42};
    static const double FractionsZ_41[] = {  1};
    fNISTMaterialDataTable[ 41].fName                 = "NIST_MAT_Mo";
    fNISTMaterialDataTable[ 41].fDensity              = 10.22 * (g/cm3);
    fNISTMaterialDataTable[ 41].fMeanExcitationEnergy = 424 * eV;
    fNISTMaterialDataTable[ 41].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 41].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 41].fNumComponents        = 1;
    fNISTMaterialDataTable[ 41].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 41].fElementList          = ElemsZ_41;
    fNISTMaterialDataTable[ 41].fElementFraction      = FractionsZ_41;
    fNISTMaterialDataTable[ 41].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Tc ---------------------------------------------------------------------------
    static const int ElemsZ_42[]        = { 43};
    static const double FractionsZ_42[] = {  1};
    fNISTMaterialDataTable[ 42].fName                 = "NIST_MAT_Tc";
    fNISTMaterialDataTable[ 42].fDensity              = 11.5 * (g/cm3);
    fNISTMaterialDataTable[ 42].fMeanExcitationEnergy = 428 * eV;
    fNISTMaterialDataTable[ 42].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 42].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 42].fNumComponents        = 1;
    fNISTMaterialDataTable[ 42].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 42].fElementList          = ElemsZ_42;
    fNISTMaterialDataTable[ 42].fElementFraction      = FractionsZ_42;
    fNISTMaterialDataTable[ 42].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ru ---------------------------------------------------------------------------
    static const int ElemsZ_43[]        = { 44};
    static const double FractionsZ_43[] = {  1};
    fNISTMaterialDataTable[ 43].fName                 = "NIST_MAT_Ru";
    fNISTMaterialDataTable[ 43].fDensity              = 12.41 * (g/cm3);
    fNISTMaterialDataTable[ 43].fMeanExcitationEnergy = 441 * eV;
    fNISTMaterialDataTable[ 43].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 43].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 43].fNumComponents        = 1;
    fNISTMaterialDataTable[ 43].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 43].fElementList          = ElemsZ_43;
    fNISTMaterialDataTable[ 43].fElementFraction      = FractionsZ_43;
    fNISTMaterialDataTable[ 43].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Rh ---------------------------------------------------------------------------
    static const int ElemsZ_44[]        = { 45};
    static const double FractionsZ_44[] = {  1};
    fNISTMaterialDataTable[ 44].fName                 = "NIST_MAT_Rh";
    fNISTMaterialDataTable[ 44].fDensity              = 12.41 * (g/cm3);
    fNISTMaterialDataTable[ 44].fMeanExcitationEnergy = 449 * eV;
    fNISTMaterialDataTable[ 44].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 44].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 44].fNumComponents        = 1;
    fNISTMaterialDataTable[ 44].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 44].fElementList          = ElemsZ_44;
    fNISTMaterialDataTable[ 44].fElementFraction      = FractionsZ_44;
    fNISTMaterialDataTable[ 44].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Pd ---------------------------------------------------------------------------
    static const int ElemsZ_45[]        = { 46};
    static const double FractionsZ_45[] = {  1};
    fNISTMaterialDataTable[ 45].fName                 = "NIST_MAT_Pd";
    fNISTMaterialDataTable[ 45].fDensity              = 12.02 * (g/cm3);
    fNISTMaterialDataTable[ 45].fMeanExcitationEnergy = 470 * eV;
    fNISTMaterialDataTable[ 45].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 45].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 45].fNumComponents        = 1;
    fNISTMaterialDataTable[ 45].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 45].fElementList          = ElemsZ_45;
    fNISTMaterialDataTable[ 45].fElementFraction      = FractionsZ_45;
    fNISTMaterialDataTable[ 45].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ag ---------------------------------------------------------------------------
    static const int ElemsZ_46[]        = { 47};
    static const double FractionsZ_46[] = {  1};
    fNISTMaterialDataTable[ 46].fName                 = "NIST_MAT_Ag";
    fNISTMaterialDataTable[ 46].fDensity              = 10.5 * (g/cm3);
    fNISTMaterialDataTable[ 46].fMeanExcitationEnergy = 470 * eV;
    fNISTMaterialDataTable[ 46].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 46].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 46].fNumComponents        = 1;
    fNISTMaterialDataTable[ 46].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 46].fElementList          = ElemsZ_46;
    fNISTMaterialDataTable[ 46].fElementFraction      = FractionsZ_46;
    fNISTMaterialDataTable[ 46].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Cd ---------------------------------------------------------------------------
    static const int ElemsZ_47[]        = { 48};
    static const double FractionsZ_47[] = {  1};
    fNISTMaterialDataTable[ 47].fName                 = "NIST_MAT_Cd";
    fNISTMaterialDataTable[ 47].fDensity              = 8.65 * (g/cm3);
    fNISTMaterialDataTable[ 47].fMeanExcitationEnergy = 469 * eV;
    fNISTMaterialDataTable[ 47].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 47].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 47].fNumComponents        = 1;
    fNISTMaterialDataTable[ 47].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 47].fElementList          = ElemsZ_47;
    fNISTMaterialDataTable[ 47].fElementFraction      = FractionsZ_47;
    fNISTMaterialDataTable[ 47].fIsBuiltByAtomCount   = true;


    // NIST_MAT_In ---------------------------------------------------------------------------
    static const int ElemsZ_48[]        = { 49};
    static const double FractionsZ_48[] = {  1};
    fNISTMaterialDataTable[ 48].fName                 = "NIST_MAT_In";
    fNISTMaterialDataTable[ 48].fDensity              = 7.31 * (g/cm3);
    fNISTMaterialDataTable[ 48].fMeanExcitationEnergy = 488 * eV;
    fNISTMaterialDataTable[ 48].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 48].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 48].fNumComponents        = 1;
    fNISTMaterialDataTable[ 48].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 48].fElementList          = ElemsZ_48;
    fNISTMaterialDataTable[ 48].fElementFraction      = FractionsZ_48;
    fNISTMaterialDataTable[ 48].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Sn ---------------------------------------------------------------------------
    static const int ElemsZ_49[]        = { 50};
    static const double FractionsZ_49[] = {  1};
    fNISTMaterialDataTable[ 49].fName                 = "NIST_MAT_Sn";
    fNISTMaterialDataTable[ 49].fDensity              = 7.31 * (g/cm3);
    fNISTMaterialDataTable[ 49].fMeanExcitationEnergy = 488 * eV;
    fNISTMaterialDataTable[ 49].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 49].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 49].fNumComponents        = 1;
    fNISTMaterialDataTable[ 49].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 49].fElementList          = ElemsZ_49;
    fNISTMaterialDataTable[ 49].fElementFraction      = FractionsZ_49;
    fNISTMaterialDataTable[ 49].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Sb ---------------------------------------------------------------------------
    static const int ElemsZ_50[]        = { 51};
    static const double FractionsZ_50[] = {  1};
    fNISTMaterialDataTable[ 50].fName                 = "NIST_MAT_Sb";
    fNISTMaterialDataTable[ 50].fDensity              = 6.691 * (g/cm3);
    fNISTMaterialDataTable[ 50].fMeanExcitationEnergy = 487 * eV;
    fNISTMaterialDataTable[ 50].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 50].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 50].fNumComponents        = 1;
    fNISTMaterialDataTable[ 50].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 50].fElementList          = ElemsZ_50;
    fNISTMaterialDataTable[ 50].fElementFraction      = FractionsZ_50;
    fNISTMaterialDataTable[ 50].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Te ---------------------------------------------------------------------------
    static const int ElemsZ_51[]        = { 52};
    static const double FractionsZ_51[] = {  1};
    fNISTMaterialDataTable[ 51].fName                 = "NIST_MAT_Te";
    fNISTMaterialDataTable[ 51].fDensity              = 6.24 * (g/cm3);
    fNISTMaterialDataTable[ 51].fMeanExcitationEnergy = 485 * eV;
    fNISTMaterialDataTable[ 51].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 51].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 51].fNumComponents        = 1;
    fNISTMaterialDataTable[ 51].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 51].fElementList          = ElemsZ_51;
    fNISTMaterialDataTable[ 51].fElementFraction      = FractionsZ_51;
    fNISTMaterialDataTable[ 51].fIsBuiltByAtomCount   = true;


    // NIST_MAT_I ----------------------------------------------------------------------------
    static const int ElemsZ_52[]        = { 53};
    static const double FractionsZ_52[] = {  1};
    fNISTMaterialDataTable[ 52].fName                 = "NIST_MAT_I";
    fNISTMaterialDataTable[ 52].fDensity              = 4.93 * (g/cm3);
    fNISTMaterialDataTable[ 52].fMeanExcitationEnergy = 491 * eV;
    fNISTMaterialDataTable[ 52].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 52].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 52].fNumComponents        = 1;
    fNISTMaterialDataTable[ 52].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 52].fElementList          = ElemsZ_52;
    fNISTMaterialDataTable[ 52].fElementFraction      = FractionsZ_52;
    fNISTMaterialDataTable[ 52].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Xe ---------------------------------------------------------------------------
    static const int ElemsZ_53[]        = { 54};
    static const double FractionsZ_53[] = {  1};
    fNISTMaterialDataTable[ 53].fName                 = "NIST_MAT_Xe";
    fNISTMaterialDataTable[ 53].fDensity              = 0.00548536 * (g/cm3);
    fNISTMaterialDataTable[ 53].fMeanExcitationEnergy = 482 * eV;
    fNISTMaterialDataTable[ 53].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 53].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 53].fNumComponents        = 1;
    fNISTMaterialDataTable[ 53].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[ 53].fElementList          = ElemsZ_53;
    fNISTMaterialDataTable[ 53].fElementFraction      = FractionsZ_53;
    fNISTMaterialDataTable[ 53].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Cs ---------------------------------------------------------------------------
    static const int ElemsZ_54[]        = { 55};
    static const double FractionsZ_54[] = {  1};
    fNISTMaterialDataTable[ 54].fName                 = "NIST_MAT_Cs";
    fNISTMaterialDataTable[ 54].fDensity              = 1.873 * (g/cm3);
    fNISTMaterialDataTable[ 54].fMeanExcitationEnergy = 488 * eV;
    fNISTMaterialDataTable[ 54].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 54].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 54].fNumComponents        = 1;
    fNISTMaterialDataTable[ 54].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 54].fElementList          = ElemsZ_54;
    fNISTMaterialDataTable[ 54].fElementFraction      = FractionsZ_54;
    fNISTMaterialDataTable[ 54].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ba ---------------------------------------------------------------------------
    static const int ElemsZ_55[]        = { 56};
    static const double FractionsZ_55[] = {  1};
    fNISTMaterialDataTable[ 55].fName                 = "NIST_MAT_Ba";
    fNISTMaterialDataTable[ 55].fDensity              = 3.5 * (g/cm3);
    fNISTMaterialDataTable[ 55].fMeanExcitationEnergy = 491 * eV;
    fNISTMaterialDataTable[ 55].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 55].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 55].fNumComponents        = 1;
    fNISTMaterialDataTable[ 55].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 55].fElementList          = ElemsZ_55;
    fNISTMaterialDataTable[ 55].fElementFraction      = FractionsZ_55;
    fNISTMaterialDataTable[ 55].fIsBuiltByAtomCount   = true;


    // NIST_MAT_La ---------------------------------------------------------------------------
    static const int ElemsZ_56[]        = { 57};
    static const double FractionsZ_56[] = {  1};
    fNISTMaterialDataTable[ 56].fName                 = "NIST_MAT_La";
    fNISTMaterialDataTable[ 56].fDensity              = 6.154 * (g/cm3);
    fNISTMaterialDataTable[ 56].fMeanExcitationEnergy = 501 * eV;
    fNISTMaterialDataTable[ 56].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 56].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 56].fNumComponents        = 1;
    fNISTMaterialDataTable[ 56].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 56].fElementList          = ElemsZ_56;
    fNISTMaterialDataTable[ 56].fElementFraction      = FractionsZ_56;
    fNISTMaterialDataTable[ 56].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ce ---------------------------------------------------------------------------
    static const int ElemsZ_57[]        = { 58};
    static const double FractionsZ_57[] = {  1};
    fNISTMaterialDataTable[ 57].fName                 = "NIST_MAT_Ce";
    fNISTMaterialDataTable[ 57].fDensity              = 6.657 * (g/cm3);
    fNISTMaterialDataTable[ 57].fMeanExcitationEnergy = 523 * eV;
    fNISTMaterialDataTable[ 57].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 57].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 57].fNumComponents        = 1;
    fNISTMaterialDataTable[ 57].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 57].fElementList          = ElemsZ_57;
    fNISTMaterialDataTable[ 57].fElementFraction      = FractionsZ_57;
    fNISTMaterialDataTable[ 57].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Pr ---------------------------------------------------------------------------
    static const int ElemsZ_58[]        = { 59};
    static const double FractionsZ_58[] = {  1};
    fNISTMaterialDataTable[ 58].fName                 = "NIST_MAT_Pr";
    fNISTMaterialDataTable[ 58].fDensity              = 6.71 * (g/cm3);
    fNISTMaterialDataTable[ 58].fMeanExcitationEnergy = 535 * eV;
    fNISTMaterialDataTable[ 58].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 58].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 58].fNumComponents        = 1;
    fNISTMaterialDataTable[ 58].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 58].fElementList          = ElemsZ_58;
    fNISTMaterialDataTable[ 58].fElementFraction      = FractionsZ_58;
    fNISTMaterialDataTable[ 58].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Nd ---------------------------------------------------------------------------
    static const int ElemsZ_59[]        = { 60};
    static const double FractionsZ_59[] = {  1};
    fNISTMaterialDataTable[ 59].fName                 = "NIST_MAT_Nd";
    fNISTMaterialDataTable[ 59].fDensity              = 6.9 * (g/cm3);
    fNISTMaterialDataTable[ 59].fMeanExcitationEnergy = 546 * eV;
    fNISTMaterialDataTable[ 59].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 59].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 59].fNumComponents        = 1;
    fNISTMaterialDataTable[ 59].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 59].fElementList          = ElemsZ_59;
    fNISTMaterialDataTable[ 59].fElementFraction      = FractionsZ_59;
    fNISTMaterialDataTable[ 59].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Pm ---------------------------------------------------------------------------
    static const int ElemsZ_60[]        = { 61};
    static const double FractionsZ_60[] = {  1};
    fNISTMaterialDataTable[ 60].fName                 = "NIST_MAT_Pm";
    fNISTMaterialDataTable[ 60].fDensity              = 7.22 * (g/cm3);
    fNISTMaterialDataTable[ 60].fMeanExcitationEnergy = 560 * eV;
    fNISTMaterialDataTable[ 60].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 60].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 60].fNumComponents        = 1;
    fNISTMaterialDataTable[ 60].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 60].fElementList          = ElemsZ_60;
    fNISTMaterialDataTable[ 60].fElementFraction      = FractionsZ_60;
    fNISTMaterialDataTable[ 60].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Sm ---------------------------------------------------------------------------
    static const int ElemsZ_61[]        = { 62};
    static const double FractionsZ_61[] = {  1};
    fNISTMaterialDataTable[ 61].fName                 = "NIST_MAT_Sm";
    fNISTMaterialDataTable[ 61].fDensity              = 7.46 * (g/cm3);
    fNISTMaterialDataTable[ 61].fMeanExcitationEnergy = 574 * eV;
    fNISTMaterialDataTable[ 61].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 61].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 61].fNumComponents        = 1;
    fNISTMaterialDataTable[ 61].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 61].fElementList          = ElemsZ_61;
    fNISTMaterialDataTable[ 61].fElementFraction      = FractionsZ_61;
    fNISTMaterialDataTable[ 61].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Eu ---------------------------------------------------------------------------
    static const int ElemsZ_62[]        = { 63};
    static const double FractionsZ_62[] = {  1};
    fNISTMaterialDataTable[ 62].fName                 = "NIST_MAT_Eu";
    fNISTMaterialDataTable[ 62].fDensity              = 5.243 * (g/cm3);
    fNISTMaterialDataTable[ 62].fMeanExcitationEnergy = 580 * eV;
    fNISTMaterialDataTable[ 62].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 62].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 62].fNumComponents        = 1;
    fNISTMaterialDataTable[ 62].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 62].fElementList          = ElemsZ_62;
    fNISTMaterialDataTable[ 62].fElementFraction      = FractionsZ_62;
    fNISTMaterialDataTable[ 62].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Gd ---------------------------------------------------------------------------
    static const int ElemsZ_63[]        = { 64};
    static const double FractionsZ_63[] = {  1};
    fNISTMaterialDataTable[ 63].fName                 = "NIST_MAT_Gd";
    fNISTMaterialDataTable[ 63].fDensity              = 7.9004 * (g/cm3);
    fNISTMaterialDataTable[ 63].fMeanExcitationEnergy = 591 * eV;
    fNISTMaterialDataTable[ 63].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 63].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 63].fNumComponents        = 1;
    fNISTMaterialDataTable[ 63].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 63].fElementList          = ElemsZ_63;
    fNISTMaterialDataTable[ 63].fElementFraction      = FractionsZ_63;
    fNISTMaterialDataTable[ 63].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Tb ---------------------------------------------------------------------------
    static const int ElemsZ_64[]        = { 65};
    static const double FractionsZ_64[] = {  1};
    fNISTMaterialDataTable[ 64].fName                 = "NIST_MAT_Tb";
    fNISTMaterialDataTable[ 64].fDensity              = 8.229 * (g/cm3);
    fNISTMaterialDataTable[ 64].fMeanExcitationEnergy = 614 * eV;
    fNISTMaterialDataTable[ 64].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 64].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 64].fNumComponents        = 1;
    fNISTMaterialDataTable[ 64].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 64].fElementList          = ElemsZ_64;
    fNISTMaterialDataTable[ 64].fElementFraction      = FractionsZ_64;
    fNISTMaterialDataTable[ 64].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Dy ---------------------------------------------------------------------------
    static const int ElemsZ_65[]        = { 66};
    static const double FractionsZ_65[] = {  1};
    fNISTMaterialDataTable[ 65].fName                 = "NIST_MAT_Dy";
    fNISTMaterialDataTable[ 65].fDensity              = 8.55 * (g/cm3);
    fNISTMaterialDataTable[ 65].fMeanExcitationEnergy = 628 * eV;
    fNISTMaterialDataTable[ 65].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 65].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 65].fNumComponents        = 1;
    fNISTMaterialDataTable[ 65].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 65].fElementList          = ElemsZ_65;
    fNISTMaterialDataTable[ 65].fElementFraction      = FractionsZ_65;
    fNISTMaterialDataTable[ 65].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ho ---------------------------------------------------------------------------
    static const int ElemsZ_66[]        = { 67};
    static const double FractionsZ_66[] = {  1};
    fNISTMaterialDataTable[ 66].fName                 = "NIST_MAT_Ho";
    fNISTMaterialDataTable[ 66].fDensity              = 8.795 * (g/cm3);
    fNISTMaterialDataTable[ 66].fMeanExcitationEnergy = 650 * eV;
    fNISTMaterialDataTable[ 66].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 66].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 66].fNumComponents        = 1;
    fNISTMaterialDataTable[ 66].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 66].fElementList          = ElemsZ_66;
    fNISTMaterialDataTable[ 66].fElementFraction      = FractionsZ_66;
    fNISTMaterialDataTable[ 66].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Er ---------------------------------------------------------------------------
    static const int ElemsZ_67[]        = { 68};
    static const double FractionsZ_67[] = {  1};
    fNISTMaterialDataTable[ 67].fName                 = "NIST_MAT_Er";
    fNISTMaterialDataTable[ 67].fDensity              = 9.066 * (g/cm3);
    fNISTMaterialDataTable[ 67].fMeanExcitationEnergy = 658 * eV;
    fNISTMaterialDataTable[ 67].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 67].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 67].fNumComponents        = 1;
    fNISTMaterialDataTable[ 67].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 67].fElementList          = ElemsZ_67;
    fNISTMaterialDataTable[ 67].fElementFraction      = FractionsZ_67;
    fNISTMaterialDataTable[ 67].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Tm ---------------------------------------------------------------------------
    static const int ElemsZ_68[]        = { 69};
    static const double FractionsZ_68[] = {  1};
    fNISTMaterialDataTable[ 68].fName                 = "NIST_MAT_Tm";
    fNISTMaterialDataTable[ 68].fDensity              = 9.321 * (g/cm3);
    fNISTMaterialDataTable[ 68].fMeanExcitationEnergy = 674 * eV;
    fNISTMaterialDataTable[ 68].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 68].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 68].fNumComponents        = 1;
    fNISTMaterialDataTable[ 68].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 68].fElementList          = ElemsZ_68;
    fNISTMaterialDataTable[ 68].fElementFraction      = FractionsZ_68;
    fNISTMaterialDataTable[ 68].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Yb ---------------------------------------------------------------------------
    static const int ElemsZ_69[]        = { 70};
    static const double FractionsZ_69[] = {  1};
    fNISTMaterialDataTable[ 69].fName                 = "NIST_MAT_Yb";
    fNISTMaterialDataTable[ 69].fDensity              = 6.73 * (g/cm3);
    fNISTMaterialDataTable[ 69].fMeanExcitationEnergy = 684 * eV;
    fNISTMaterialDataTable[ 69].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 69].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 69].fNumComponents        = 1;
    fNISTMaterialDataTable[ 69].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 69].fElementList          = ElemsZ_69;
    fNISTMaterialDataTable[ 69].fElementFraction      = FractionsZ_69;
    fNISTMaterialDataTable[ 69].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Lu ---------------------------------------------------------------------------
    static const int ElemsZ_70[]        = { 71};
    static const double FractionsZ_70[] = {  1};
    fNISTMaterialDataTable[ 70].fName                 = "NIST_MAT_Lu";
    fNISTMaterialDataTable[ 70].fDensity              = 9.84 * (g/cm3);
    fNISTMaterialDataTable[ 70].fMeanExcitationEnergy = 694 * eV;
    fNISTMaterialDataTable[ 70].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 70].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 70].fNumComponents        = 1;
    fNISTMaterialDataTable[ 70].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 70].fElementList          = ElemsZ_70;
    fNISTMaterialDataTable[ 70].fElementFraction      = FractionsZ_70;
    fNISTMaterialDataTable[ 70].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Hf ---------------------------------------------------------------------------
    static const int ElemsZ_71[]        = { 72};
    static const double FractionsZ_71[] = {  1};
    fNISTMaterialDataTable[ 71].fName                 = "NIST_MAT_Hf";
    fNISTMaterialDataTable[ 71].fDensity              = 13.31 * (g/cm3);
    fNISTMaterialDataTable[ 71].fMeanExcitationEnergy = 705 * eV;
    fNISTMaterialDataTable[ 71].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 71].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 71].fNumComponents        = 1;
    fNISTMaterialDataTable[ 71].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 71].fElementList          = ElemsZ_71;
    fNISTMaterialDataTable[ 71].fElementFraction      = FractionsZ_71;
    fNISTMaterialDataTable[ 71].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ta ---------------------------------------------------------------------------
    static const int ElemsZ_72[]        = { 73};
    static const double FractionsZ_72[] = {  1};
    fNISTMaterialDataTable[ 72].fName                 = "NIST_MAT_Ta";
    fNISTMaterialDataTable[ 72].fDensity              = 16.654 * (g/cm3);
    fNISTMaterialDataTable[ 72].fMeanExcitationEnergy = 718 * eV;
    fNISTMaterialDataTable[ 72].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 72].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 72].fNumComponents        = 1;
    fNISTMaterialDataTable[ 72].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 72].fElementList          = ElemsZ_72;
    fNISTMaterialDataTable[ 72].fElementFraction      = FractionsZ_72;
    fNISTMaterialDataTable[ 72].fIsBuiltByAtomCount   = true;


    // NIST_MAT_W ----------------------------------------------------------------------------
    static const int ElemsZ_73[]        = { 74};
    static const double FractionsZ_73[] = {  1};
    fNISTMaterialDataTable[ 73].fName                 = "NIST_MAT_W";
    fNISTMaterialDataTable[ 73].fDensity              = 19.3 * (g/cm3);
    fNISTMaterialDataTable[ 73].fMeanExcitationEnergy = 727 * eV;
    fNISTMaterialDataTable[ 73].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 73].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 73].fNumComponents        = 1;
    fNISTMaterialDataTable[ 73].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 73].fElementList          = ElemsZ_73;
    fNISTMaterialDataTable[ 73].fElementFraction      = FractionsZ_73;
    fNISTMaterialDataTable[ 73].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Re ---------------------------------------------------------------------------
    static const int ElemsZ_74[]        = { 75};
    static const double FractionsZ_74[] = {  1};
    fNISTMaterialDataTable[ 74].fName                 = "NIST_MAT_Re";
    fNISTMaterialDataTable[ 74].fDensity              = 21.02 * (g/cm3);
    fNISTMaterialDataTable[ 74].fMeanExcitationEnergy = 736 * eV;
    fNISTMaterialDataTable[ 74].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 74].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 74].fNumComponents        = 1;
    fNISTMaterialDataTable[ 74].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 74].fElementList          = ElemsZ_74;
    fNISTMaterialDataTable[ 74].fElementFraction      = FractionsZ_74;
    fNISTMaterialDataTable[ 74].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Os ---------------------------------------------------------------------------
    static const int ElemsZ_75[]        = { 76};
    static const double FractionsZ_75[] = {  1};
    fNISTMaterialDataTable[ 75].fName                 = "NIST_MAT_Os";
    fNISTMaterialDataTable[ 75].fDensity              = 22.57 * (g/cm3);
    fNISTMaterialDataTable[ 75].fMeanExcitationEnergy = 746 * eV;
    fNISTMaterialDataTable[ 75].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 75].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 75].fNumComponents        = 1;
    fNISTMaterialDataTable[ 75].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 75].fElementList          = ElemsZ_75;
    fNISTMaterialDataTable[ 75].fElementFraction      = FractionsZ_75;
    fNISTMaterialDataTable[ 75].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ir ---------------------------------------------------------------------------
    static const int ElemsZ_76[]        = { 77};
    static const double FractionsZ_76[] = {  1};
    fNISTMaterialDataTable[ 76].fName                 = "NIST_MAT_Ir";
    fNISTMaterialDataTable[ 76].fDensity              = 22.42 * (g/cm3);
    fNISTMaterialDataTable[ 76].fMeanExcitationEnergy = 757 * eV;
    fNISTMaterialDataTable[ 76].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 76].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 76].fNumComponents        = 1;
    fNISTMaterialDataTable[ 76].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 76].fElementList          = ElemsZ_76;
    fNISTMaterialDataTable[ 76].fElementFraction      = FractionsZ_76;
    fNISTMaterialDataTable[ 76].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Pt ---------------------------------------------------------------------------
    static const int ElemsZ_77[]        = { 78};
    static const double FractionsZ_77[] = {  1};
    fNISTMaterialDataTable[ 77].fName                 = "NIST_MAT_Pt";
    fNISTMaterialDataTable[ 77].fDensity              = 21.45 * (g/cm3);
    fNISTMaterialDataTable[ 77].fMeanExcitationEnergy = 790 * eV;
    fNISTMaterialDataTable[ 77].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 77].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 77].fNumComponents        = 1;
    fNISTMaterialDataTable[ 77].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 77].fElementList          = ElemsZ_77;
    fNISTMaterialDataTable[ 77].fElementFraction      = FractionsZ_77;
    fNISTMaterialDataTable[ 77].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Au ---------------------------------------------------------------------------
    static const int ElemsZ_78[]        = { 79};
    static const double FractionsZ_78[] = {  1};
    fNISTMaterialDataTable[ 78].fName                 = "NIST_MAT_Au";
    fNISTMaterialDataTable[ 78].fDensity              = 19.32 * (g/cm3);
    fNISTMaterialDataTable[ 78].fMeanExcitationEnergy = 790 * eV;
    fNISTMaterialDataTable[ 78].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 78].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 78].fNumComponents        = 1;
    fNISTMaterialDataTable[ 78].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 78].fElementList          = ElemsZ_78;
    fNISTMaterialDataTable[ 78].fElementFraction      = FractionsZ_78;
    fNISTMaterialDataTable[ 78].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Hg ---------------------------------------------------------------------------
    static const int ElemsZ_79[]        = { 80};
    static const double FractionsZ_79[] = {  1};
    fNISTMaterialDataTable[ 79].fName                 = "NIST_MAT_Hg";
    fNISTMaterialDataTable[ 79].fDensity              = 13.546 * (g/cm3);
    fNISTMaterialDataTable[ 79].fMeanExcitationEnergy = 800 * eV;
    fNISTMaterialDataTable[ 79].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 79].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 79].fNumComponents        = 1;
    fNISTMaterialDataTable[ 79].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 79].fElementList          = ElemsZ_79;
    fNISTMaterialDataTable[ 79].fElementFraction      = FractionsZ_79;
    fNISTMaterialDataTable[ 79].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Tl ---------------------------------------------------------------------------
    static const int ElemsZ_80[]        = { 81};
    static const double FractionsZ_80[] = {  1};
    fNISTMaterialDataTable[ 80].fName                 = "NIST_MAT_Tl";
    fNISTMaterialDataTable[ 80].fDensity              = 11.72 * (g/cm3);
    fNISTMaterialDataTable[ 80].fMeanExcitationEnergy = 810 * eV;
    fNISTMaterialDataTable[ 80].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 80].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 80].fNumComponents        = 1;
    fNISTMaterialDataTable[ 80].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 80].fElementList          = ElemsZ_80;
    fNISTMaterialDataTable[ 80].fElementFraction      = FractionsZ_80;
    fNISTMaterialDataTable[ 80].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Pb ---------------------------------------------------------------------------
    static const int ElemsZ_81[]        = { 82};
    static const double FractionsZ_81[] = {  1};
    fNISTMaterialDataTable[ 81].fName                 = "NIST_MAT_Pb";
    fNISTMaterialDataTable[ 81].fDensity              = 11.35 * (g/cm3);
    fNISTMaterialDataTable[ 81].fMeanExcitationEnergy = 823 * eV;
    fNISTMaterialDataTable[ 81].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 81].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 81].fNumComponents        = 1;
    fNISTMaterialDataTable[ 81].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 81].fElementList          = ElemsZ_81;
    fNISTMaterialDataTable[ 81].fElementFraction      = FractionsZ_81;
    fNISTMaterialDataTable[ 81].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Bi ---------------------------------------------------------------------------
    static const int ElemsZ_82[]        = { 83};
    static const double FractionsZ_82[] = {  1};
    fNISTMaterialDataTable[ 82].fName                 = "NIST_MAT_Bi";
    fNISTMaterialDataTable[ 82].fDensity              = 9.747 * (g/cm3);
    fNISTMaterialDataTable[ 82].fMeanExcitationEnergy = 823 * eV;
    fNISTMaterialDataTable[ 82].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 82].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 82].fNumComponents        = 1;
    fNISTMaterialDataTable[ 82].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 82].fElementList          = ElemsZ_82;
    fNISTMaterialDataTable[ 82].fElementFraction      = FractionsZ_82;
    fNISTMaterialDataTable[ 82].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Po ---------------------------------------------------------------------------
    static const int ElemsZ_83[]        = { 84};
    static const double FractionsZ_83[] = {  1};
    fNISTMaterialDataTable[ 83].fName                 = "NIST_MAT_Po";
    fNISTMaterialDataTable[ 83].fDensity              = 9.32 * (g/cm3);
    fNISTMaterialDataTable[ 83].fMeanExcitationEnergy = 830 * eV;
    fNISTMaterialDataTable[ 83].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 83].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 83].fNumComponents        = 1;
    fNISTMaterialDataTable[ 83].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 83].fElementList          = ElemsZ_83;
    fNISTMaterialDataTable[ 83].fElementFraction      = FractionsZ_83;
    fNISTMaterialDataTable[ 83].fIsBuiltByAtomCount   = true;


    // NIST_MAT_At ---------------------------------------------------------------------------
    static const int ElemsZ_84[]        = { 85};
    static const double FractionsZ_84[] = {  1};
    fNISTMaterialDataTable[ 84].fName                 = "NIST_MAT_At";
    fNISTMaterialDataTable[ 84].fDensity              = 9.32 * (g/cm3);
    fNISTMaterialDataTable[ 84].fMeanExcitationEnergy = 825 * eV;
    fNISTMaterialDataTable[ 84].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 84].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 84].fNumComponents        = 1;
    fNISTMaterialDataTable[ 84].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 84].fElementList          = ElemsZ_84;
    fNISTMaterialDataTable[ 84].fElementFraction      = FractionsZ_84;
    fNISTMaterialDataTable[ 84].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Rn ---------------------------------------------------------------------------
    static const int ElemsZ_85[]        = { 86};
    static const double FractionsZ_85[] = {  1};
    fNISTMaterialDataTable[ 85].fName                 = "NIST_MAT_Rn";
    fNISTMaterialDataTable[ 85].fDensity              = 0.00900662 * (g/cm3);
    fNISTMaterialDataTable[ 85].fMeanExcitationEnergy = 794 * eV;
    fNISTMaterialDataTable[ 85].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 85].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 85].fNumComponents        = 1;
    fNISTMaterialDataTable[ 85].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[ 85].fElementList          = ElemsZ_85;
    fNISTMaterialDataTable[ 85].fElementFraction      = FractionsZ_85;
    fNISTMaterialDataTable[ 85].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Fr ---------------------------------------------------------------------------
    static const int ElemsZ_86[]        = { 87};
    static const double FractionsZ_86[] = {  1};
    fNISTMaterialDataTable[ 86].fName                 = "NIST_MAT_Fr";
    fNISTMaterialDataTable[ 86].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[ 86].fMeanExcitationEnergy = 827 * eV;
    fNISTMaterialDataTable[ 86].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 86].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 86].fNumComponents        = 1;
    fNISTMaterialDataTable[ 86].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 86].fElementList          = ElemsZ_86;
    fNISTMaterialDataTable[ 86].fElementFraction      = FractionsZ_86;
    fNISTMaterialDataTable[ 86].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ra ---------------------------------------------------------------------------
    static const int ElemsZ_87[]        = { 88};
    static const double FractionsZ_87[] = {  1};
    fNISTMaterialDataTable[ 87].fName                 = "NIST_MAT_Ra";
    fNISTMaterialDataTable[ 87].fDensity              = 5 * (g/cm3);
    fNISTMaterialDataTable[ 87].fMeanExcitationEnergy = 826 * eV;
    fNISTMaterialDataTable[ 87].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 87].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 87].fNumComponents        = 1;
    fNISTMaterialDataTable[ 87].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 87].fElementList          = ElemsZ_87;
    fNISTMaterialDataTable[ 87].fElementFraction      = FractionsZ_87;
    fNISTMaterialDataTable[ 87].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Ac ---------------------------------------------------------------------------
    static const int ElemsZ_88[]        = { 89};
    static const double FractionsZ_88[] = {  1};
    fNISTMaterialDataTable[ 88].fName                 = "NIST_MAT_Ac";
    fNISTMaterialDataTable[ 88].fDensity              = 10.07 * (g/cm3);
    fNISTMaterialDataTable[ 88].fMeanExcitationEnergy = 841 * eV;
    fNISTMaterialDataTable[ 88].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 88].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 88].fNumComponents        = 1;
    fNISTMaterialDataTable[ 88].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 88].fElementList          = ElemsZ_88;
    fNISTMaterialDataTable[ 88].fElementFraction      = FractionsZ_88;
    fNISTMaterialDataTable[ 88].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Th ---------------------------------------------------------------------------
    static const int ElemsZ_89[]        = { 90};
    static const double FractionsZ_89[] = {  1};
    fNISTMaterialDataTable[ 89].fName                 = "NIST_MAT_Th";
    fNISTMaterialDataTable[ 89].fDensity              = 11.72 * (g/cm3);
    fNISTMaterialDataTable[ 89].fMeanExcitationEnergy = 847 * eV;
    fNISTMaterialDataTable[ 89].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 89].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 89].fNumComponents        = 1;
    fNISTMaterialDataTable[ 89].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 89].fElementList          = ElemsZ_89;
    fNISTMaterialDataTable[ 89].fElementFraction      = FractionsZ_89;
    fNISTMaterialDataTable[ 89].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Pa ---------------------------------------------------------------------------
    static const int ElemsZ_90[]        = { 91};
    static const double FractionsZ_90[] = {  1};
    fNISTMaterialDataTable[ 90].fName                 = "NIST_MAT_Pa";
    fNISTMaterialDataTable[ 90].fDensity              = 15.37 * (g/cm3);
    fNISTMaterialDataTable[ 90].fMeanExcitationEnergy = 878 * eV;
    fNISTMaterialDataTable[ 90].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 90].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 90].fNumComponents        = 1;
    fNISTMaterialDataTable[ 90].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 90].fElementList          = ElemsZ_90;
    fNISTMaterialDataTable[ 90].fElementFraction      = FractionsZ_90;
    fNISTMaterialDataTable[ 90].fIsBuiltByAtomCount   = true;


    // NIST_MAT_U ----------------------------------------------------------------------------
    static const int ElemsZ_91[]        = { 92};
    static const double FractionsZ_91[] = {  1};
    fNISTMaterialDataTable[ 91].fName                 = "NIST_MAT_U";
    fNISTMaterialDataTable[ 91].fDensity              = 18.95 * (g/cm3);
    fNISTMaterialDataTable[ 91].fMeanExcitationEnergy = 890 * eV;
    fNISTMaterialDataTable[ 91].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 91].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 91].fNumComponents        = 1;
    fNISTMaterialDataTable[ 91].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 91].fElementList          = ElemsZ_91;
    fNISTMaterialDataTable[ 91].fElementFraction      = FractionsZ_91;
    fNISTMaterialDataTable[ 91].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Np ---------------------------------------------------------------------------
    static const int ElemsZ_92[]        = { 93};
    static const double FractionsZ_92[] = {  1};
    fNISTMaterialDataTable[ 92].fName                 = "NIST_MAT_Np";
    fNISTMaterialDataTable[ 92].fDensity              = 20.25 * (g/cm3);
    fNISTMaterialDataTable[ 92].fMeanExcitationEnergy = 902 * eV;
    fNISTMaterialDataTable[ 92].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 92].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 92].fNumComponents        = 1;
    fNISTMaterialDataTable[ 92].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 92].fElementList          = ElemsZ_92;
    fNISTMaterialDataTable[ 92].fElementFraction      = FractionsZ_92;
    fNISTMaterialDataTable[ 92].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Pu ---------------------------------------------------------------------------
    static const int ElemsZ_93[]        = { 94};
    static const double FractionsZ_93[] = {  1};
    fNISTMaterialDataTable[ 93].fName                 = "NIST_MAT_Pu";
    fNISTMaterialDataTable[ 93].fDensity              = 19.84 * (g/cm3);
    fNISTMaterialDataTable[ 93].fMeanExcitationEnergy = 921 * eV;
    fNISTMaterialDataTable[ 93].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 93].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 93].fNumComponents        = 1;
    fNISTMaterialDataTable[ 93].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 93].fElementList          = ElemsZ_93;
    fNISTMaterialDataTable[ 93].fElementFraction      = FractionsZ_93;
    fNISTMaterialDataTable[ 93].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Am ---------------------------------------------------------------------------
    static const int ElemsZ_94[]        = { 95};
    static const double FractionsZ_94[] = {  1};
    fNISTMaterialDataTable[ 94].fName                 = "NIST_MAT_Am";
    fNISTMaterialDataTable[ 94].fDensity              = 13.67 * (g/cm3);
    fNISTMaterialDataTable[ 94].fMeanExcitationEnergy = 934 * eV;
    fNISTMaterialDataTable[ 94].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 94].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 94].fNumComponents        = 1;
    fNISTMaterialDataTable[ 94].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 94].fElementList          = ElemsZ_94;
    fNISTMaterialDataTable[ 94].fElementFraction      = FractionsZ_94;
    fNISTMaterialDataTable[ 94].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Cm ---------------------------------------------------------------------------
    static const int ElemsZ_95[]        = { 96};
    static const double FractionsZ_95[] = {  1};
    fNISTMaterialDataTable[ 95].fName                 = "NIST_MAT_Cm";
    fNISTMaterialDataTable[ 95].fDensity              = 13.51 * (g/cm3);
    fNISTMaterialDataTable[ 95].fMeanExcitationEnergy = 939 * eV;
    fNISTMaterialDataTable[ 95].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 95].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 95].fNumComponents        = 1;
    fNISTMaterialDataTable[ 95].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 95].fElementList          = ElemsZ_95;
    fNISTMaterialDataTable[ 95].fElementFraction      = FractionsZ_95;
    fNISTMaterialDataTable[ 95].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Bk ---------------------------------------------------------------------------
    static const int ElemsZ_96[]        = { 97};
    static const double FractionsZ_96[] = {  1};
    fNISTMaterialDataTable[ 96].fName                 = "NIST_MAT_Bk";
    fNISTMaterialDataTable[ 96].fDensity              = 14 * (g/cm3);
    fNISTMaterialDataTable[ 96].fMeanExcitationEnergy = 952 * eV;
    fNISTMaterialDataTable[ 96].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 96].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 96].fNumComponents        = 1;
    fNISTMaterialDataTable[ 96].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 96].fElementList          = ElemsZ_96;
    fNISTMaterialDataTable[ 96].fElementFraction      = FractionsZ_96;
    fNISTMaterialDataTable[ 96].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Cf ---------------------------------------------------------------------------
    static const int ElemsZ_97[]        = { 98};
    static const double FractionsZ_97[] = {  1};
    fNISTMaterialDataTable[ 97].fName                 = "NIST_MAT_Cf";
    fNISTMaterialDataTable[ 97].fDensity              = 10 * (g/cm3);
    fNISTMaterialDataTable[ 97].fMeanExcitationEnergy = 966 * eV;
    fNISTMaterialDataTable[ 97].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 97].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 97].fNumComponents        = 1;
    fNISTMaterialDataTable[ 97].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 97].fElementList          = ElemsZ_97;
    fNISTMaterialDataTable[ 97].fElementFraction      = FractionsZ_97;
    fNISTMaterialDataTable[ 97].fIsBuiltByAtomCount   = true;


    // =======================================================================================
    // Compound Materials :
    // =======================================================================================

    // NIST_MAT_A-150_TISSUE -----------------------------------------------------------------
    static const int ElemsZ_98[]        = {  1,   6,   7,   8,   9,  20};
    static const double FractionsZ_98[] = {  0.101327,   0.775501,   0.035057,
                                              0.052316,   0.017422,   0.018378};
    fNISTMaterialDataTable[ 98].fName                 = "NIST_MAT_A-150_TISSUE";
    fNISTMaterialDataTable[ 98].fDensity              = 1.127 * (g/cm3);
    fNISTMaterialDataTable[ 98].fMeanExcitationEnergy = 65.1 * eV;
    fNISTMaterialDataTable[ 98].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 98].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 98].fNumComponents        = 6;
    fNISTMaterialDataTable[ 98].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 98].fElementList          = ElemsZ_98;
    fNISTMaterialDataTable[ 98].fElementFraction      = FractionsZ_98;
    fNISTMaterialDataTable[ 98].fIsBuiltByAtomCount   = false;


    // NIST_MAT_ACETONE ----------------------------------------------------------------------
    static const int ElemsZ_99[]        = {  6,   1,   8};
    static const double FractionsZ_99[] = {  3,   6,   1};
    fNISTMaterialDataTable[ 99].fName                 = "NIST_MAT_ACETONE";
    fNISTMaterialDataTable[ 99].fDensity              = 0.7899 * (g/cm3);
    fNISTMaterialDataTable[ 99].fMeanExcitationEnergy = 64.2 * eV;
    fNISTMaterialDataTable[ 99].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[ 99].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[ 99].fNumComponents        = 3;
    fNISTMaterialDataTable[ 99].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[ 99].fElementList          = ElemsZ_99;
    fNISTMaterialDataTable[ 99].fElementFraction      = FractionsZ_99;
    fNISTMaterialDataTable[ 99].fIsBuiltByAtomCount   = true;


    // NIST_MAT_ACETYLENE --------------------------------------------------------------------
    static const int ElemsZ_100[]        = {  6,   1};
    static const double FractionsZ_100[] = {  2,   2};
    fNISTMaterialDataTable[100].fName                 = "NIST_MAT_ACETYLENE";
    fNISTMaterialDataTable[100].fDensity              = 0.0010967 * (g/cm3);
    fNISTMaterialDataTable[100].fMeanExcitationEnergy = 58.2 * eV;
    fNISTMaterialDataTable[100].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[100].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[100].fNumComponents        = 2;
    fNISTMaterialDataTable[100].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[100].fElementList          = ElemsZ_100;
    fNISTMaterialDataTable[100].fElementFraction      = FractionsZ_100;
    fNISTMaterialDataTable[100].fIsBuiltByAtomCount   = true;


    // NIST_MAT_ADENINE ----------------------------------------------------------------------
    static const int ElemsZ_101[]        = {  6,   1,   7};
    static const double FractionsZ_101[] = {  5,   5,   5};
    fNISTMaterialDataTable[101].fName                 = "NIST_MAT_ADENINE";
    fNISTMaterialDataTable[101].fDensity              = 1.6 * (g/cm3);
    fNISTMaterialDataTable[101].fMeanExcitationEnergy = 71.4 * eV;
    fNISTMaterialDataTable[101].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[101].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[101].fNumComponents        = 3;
    fNISTMaterialDataTable[101].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[101].fElementList          = ElemsZ_101;
    fNISTMaterialDataTable[101].fElementFraction      = FractionsZ_101;
    fNISTMaterialDataTable[101].fIsBuiltByAtomCount   = true;


    // NIST_MAT_ADIPOSE_TISSUE_ICRP ----------------------------------------------------------
    static const int ElemsZ_102[]        = {  1,   6,   7,   8,  11,  16,
                                             17};
    static const double FractionsZ_102[] = {     0.114,      0.598,      0.007,
                                                 0.278,      0.001,      0.001,
                                                 0.001};
    fNISTMaterialDataTable[102].fName                 = "NIST_MAT_ADIPOSE_TISSUE_ICRP";
    fNISTMaterialDataTable[102].fDensity              = 0.95 * (g/cm3);
    fNISTMaterialDataTable[102].fMeanExcitationEnergy = 63.2 * eV;
    fNISTMaterialDataTable[102].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[102].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[102].fNumComponents        = 7;
    fNISTMaterialDataTable[102].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[102].fElementList          = ElemsZ_102;
    fNISTMaterialDataTable[102].fElementFraction      = FractionsZ_102;
    fNISTMaterialDataTable[102].fIsBuiltByAtomCount   = false;


    // NIST_MAT_AIR --------------------------------------------------------------------------
    static const int ElemsZ_103[]        = {  6,   7,   8,  18};
    static const double FractionsZ_103[] = {  0.000124,   0.755267,   0.231781,
                                              0.012827};
    fNISTMaterialDataTable[103].fName                 = "NIST_MAT_AIR";
    fNISTMaterialDataTable[103].fDensity              = 0.00120479 * (g/cm3);
    fNISTMaterialDataTable[103].fMeanExcitationEnergy = 85.7 * eV;
    fNISTMaterialDataTable[103].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[103].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[103].fNumComponents        = 4;
    fNISTMaterialDataTable[103].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[103].fElementList          = ElemsZ_103;
    fNISTMaterialDataTable[103].fElementFraction      = FractionsZ_103;
    fNISTMaterialDataTable[103].fIsBuiltByAtomCount   = false;


    // NIST_MAT_ALANINE ----------------------------------------------------------------------
    static const int ElemsZ_104[]        = {  6,   1,   7,   8};
    static const double FractionsZ_104[] = {  3,   7,   1,   2};
    fNISTMaterialDataTable[104].fName                 = "NIST_MAT_ALANINE";
    fNISTMaterialDataTable[104].fDensity              = 1.42 * (g/cm3);
    fNISTMaterialDataTable[104].fMeanExcitationEnergy = 71.9 * eV;
    fNISTMaterialDataTable[104].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[104].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[104].fNumComponents        = 4;
    fNISTMaterialDataTable[104].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[104].fElementList          = ElemsZ_104;
    fNISTMaterialDataTable[104].fElementFraction      = FractionsZ_104;
    fNISTMaterialDataTable[104].fIsBuiltByAtomCount   = true;


    // NIST_MAT_ALUMINUM_OXIDE ---------------------------------------------------------------
    static const int ElemsZ_105[]        = { 13,   8};
    static const double FractionsZ_105[] = {  2,   3};
    fNISTMaterialDataTable[105].fName                 = "NIST_MAT_ALUMINUM_OXIDE";
    fNISTMaterialDataTable[105].fDensity              = 3.97 * (g/cm3);
    fNISTMaterialDataTable[105].fMeanExcitationEnergy = 145.2 * eV;
    fNISTMaterialDataTable[105].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[105].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[105].fNumComponents        = 2;
    fNISTMaterialDataTable[105].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[105].fElementList          = ElemsZ_105;
    fNISTMaterialDataTable[105].fElementFraction      = FractionsZ_105;
    fNISTMaterialDataTable[105].fIsBuiltByAtomCount   = true;


    // NIST_MAT_AMBER ------------------------------------------------------------------------
    static const int ElemsZ_106[]        = {  1,   6,   8};
    static const double FractionsZ_106[] = {   0.10593,   0.788973,   0.105096};
    fNISTMaterialDataTable[106].fName                 = "NIST_MAT_AMBER";
    fNISTMaterialDataTable[106].fDensity              = 1.1 * (g/cm3);
    fNISTMaterialDataTable[106].fMeanExcitationEnergy = 63.2 * eV;
    fNISTMaterialDataTable[106].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[106].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[106].fNumComponents        = 3;
    fNISTMaterialDataTable[106].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[106].fElementList          = ElemsZ_106;
    fNISTMaterialDataTable[106].fElementFraction      = FractionsZ_106;
    fNISTMaterialDataTable[106].fIsBuiltByAtomCount   = false;


    // NIST_MAT_AMMONIA ----------------------------------------------------------------------
    static const int ElemsZ_107[]        = {  7,   1};
    static const double FractionsZ_107[] = {  1,   3};
    fNISTMaterialDataTable[107].fName                 = "NIST_MAT_AMMONIA";
    fNISTMaterialDataTable[107].fDensity              = 0.000826019 * (g/cm3);
    fNISTMaterialDataTable[107].fMeanExcitationEnergy = 53.7 * eV;
    fNISTMaterialDataTable[107].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[107].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[107].fNumComponents        = 2;
    fNISTMaterialDataTable[107].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[107].fElementList          = ElemsZ_107;
    fNISTMaterialDataTable[107].fElementFraction      = FractionsZ_107;
    fNISTMaterialDataTable[107].fIsBuiltByAtomCount   = true;


    // NIST_MAT_ANILINE ----------------------------------------------------------------------
    static const int ElemsZ_108[]        = {  6,   1,   7};
    static const double FractionsZ_108[] = {  6,   7,   1};
    fNISTMaterialDataTable[108].fName                 = "NIST_MAT_ANILINE";
    fNISTMaterialDataTable[108].fDensity              = 1.0235 * (g/cm3);
    fNISTMaterialDataTable[108].fMeanExcitationEnergy = 66.2 * eV;
    fNISTMaterialDataTable[108].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[108].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[108].fNumComponents        = 3;
    fNISTMaterialDataTable[108].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[108].fElementList          = ElemsZ_108;
    fNISTMaterialDataTable[108].fElementFraction      = FractionsZ_108;
    fNISTMaterialDataTable[108].fIsBuiltByAtomCount   = true;


    // NIST_MAT_ANTHRACENE -------------------------------------------------------------------
    static const int ElemsZ_109[]        = {  6,   1};
    static const double FractionsZ_109[] = { 14,  10};
    fNISTMaterialDataTable[109].fName                 = "NIST_MAT_ANTHRACENE";
    fNISTMaterialDataTable[109].fDensity              = 1.283 * (g/cm3);
    fNISTMaterialDataTable[109].fMeanExcitationEnergy = 69.5 * eV;
    fNISTMaterialDataTable[109].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[109].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[109].fNumComponents        = 2;
    fNISTMaterialDataTable[109].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[109].fElementList          = ElemsZ_109;
    fNISTMaterialDataTable[109].fElementFraction      = FractionsZ_109;
    fNISTMaterialDataTable[109].fIsBuiltByAtomCount   = true;


    // NIST_MAT_B-100_BONE -------------------------------------------------------------------
    static const int ElemsZ_110[]        = {  1,   6,   7,   8,   9,  20};
    static const double FractionsZ_110[] = {  0.065471,   0.536945,     0.0215,
                                              0.032085,   0.167411,   0.176589};
    fNISTMaterialDataTable[110].fName                 = "NIST_MAT_B-100_BONE";
    fNISTMaterialDataTable[110].fDensity              = 1.45 * (g/cm3);
    fNISTMaterialDataTable[110].fMeanExcitationEnergy = 85.9 * eV;
    fNISTMaterialDataTable[110].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[110].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[110].fNumComponents        = 6;
    fNISTMaterialDataTable[110].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[110].fElementList          = ElemsZ_110;
    fNISTMaterialDataTable[110].fElementFraction      = FractionsZ_110;
    fNISTMaterialDataTable[110].fIsBuiltByAtomCount   = false;


    // NIST_MAT_BAKELITE ---------------------------------------------------------------------
    static const int ElemsZ_111[]        = {  1,   6,   8};
    static const double FractionsZ_111[] = {  0.057441,   0.774591,   0.167968};
    fNISTMaterialDataTable[111].fName                 = "NIST_MAT_BAKELITE";
    fNISTMaterialDataTable[111].fDensity              = 1.25 * (g/cm3);
    fNISTMaterialDataTable[111].fMeanExcitationEnergy = 72.4 * eV;
    fNISTMaterialDataTable[111].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[111].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[111].fNumComponents        = 3;
    fNISTMaterialDataTable[111].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[111].fElementList          = ElemsZ_111;
    fNISTMaterialDataTable[111].fElementFraction      = FractionsZ_111;
    fNISTMaterialDataTable[111].fIsBuiltByAtomCount   = false;


    // NIST_MAT_BARIUM_FLUORIDE --------------------------------------------------------------
    static const int ElemsZ_112[]        = { 56,   9};
    static const double FractionsZ_112[] = {  1,   2};
    fNISTMaterialDataTable[112].fName                 = "NIST_MAT_BARIUM_FLUORIDE";
    fNISTMaterialDataTable[112].fDensity              = 4.89 * (g/cm3);
    fNISTMaterialDataTable[112].fMeanExcitationEnergy = 375.9 * eV;
    fNISTMaterialDataTable[112].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[112].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[112].fNumComponents        = 2;
    fNISTMaterialDataTable[112].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[112].fElementList          = ElemsZ_112;
    fNISTMaterialDataTable[112].fElementFraction      = FractionsZ_112;
    fNISTMaterialDataTable[112].fIsBuiltByAtomCount   = true;


    // NIST_MAT_BARIUM_SULFATE ---------------------------------------------------------------
    static const int ElemsZ_113[]        = { 56,  16,   8};
    static const double FractionsZ_113[] = {  1,   1,   4};
    fNISTMaterialDataTable[113].fName                 = "NIST_MAT_BARIUM_SULFATE";
    fNISTMaterialDataTable[113].fDensity              = 4.5 * (g/cm3);
    fNISTMaterialDataTable[113].fMeanExcitationEnergy = 285.7 * eV;
    fNISTMaterialDataTable[113].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[113].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[113].fNumComponents        = 3;
    fNISTMaterialDataTable[113].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[113].fElementList          = ElemsZ_113;
    fNISTMaterialDataTable[113].fElementFraction      = FractionsZ_113;
    fNISTMaterialDataTable[113].fIsBuiltByAtomCount   = true;


    // NIST_MAT_BENZENE ----------------------------------------------------------------------
    static const int ElemsZ_114[]        = {  6,   1};
    static const double FractionsZ_114[] = {  6,   6};
    fNISTMaterialDataTable[114].fName                 = "NIST_MAT_BENZENE";
    fNISTMaterialDataTable[114].fDensity              = 0.87865 * (g/cm3);
    fNISTMaterialDataTable[114].fMeanExcitationEnergy = 63.4 * eV;
    fNISTMaterialDataTable[114].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[114].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[114].fNumComponents        = 2;
    fNISTMaterialDataTable[114].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[114].fElementList          = ElemsZ_114;
    fNISTMaterialDataTable[114].fElementFraction      = FractionsZ_114;
    fNISTMaterialDataTable[114].fIsBuiltByAtomCount   = true;


    // NIST_MAT_BERYLLIUM_OXIDE --------------------------------------------------------------
    static const int ElemsZ_115[]        = {  4,   8};
    static const double FractionsZ_115[] = {  1,   1};
    fNISTMaterialDataTable[115].fName                 = "NIST_MAT_BERYLLIUM_OXIDE";
    fNISTMaterialDataTable[115].fDensity              = 3.01 * (g/cm3);
    fNISTMaterialDataTable[115].fMeanExcitationEnergy = 93.2 * eV;
    fNISTMaterialDataTable[115].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[115].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[115].fNumComponents        = 2;
    fNISTMaterialDataTable[115].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[115].fElementList          = ElemsZ_115;
    fNISTMaterialDataTable[115].fElementFraction      = FractionsZ_115;
    fNISTMaterialDataTable[115].fIsBuiltByAtomCount   = true;


    // NIST_MAT_BGO --------------------------------------------------------------------------
    static const int ElemsZ_116[]        = { 83,  32,   8};
    static const double FractionsZ_116[] = {  4,   3,  12};
    fNISTMaterialDataTable[116].fName                 = "NIST_MAT_BGO";
    fNISTMaterialDataTable[116].fDensity              = 7.13 * (g/cm3);
    fNISTMaterialDataTable[116].fMeanExcitationEnergy = 534.1 * eV;
    fNISTMaterialDataTable[116].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[116].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[116].fNumComponents        = 3;
    fNISTMaterialDataTable[116].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[116].fElementList          = ElemsZ_116;
    fNISTMaterialDataTable[116].fElementFraction      = FractionsZ_116;
    fNISTMaterialDataTable[116].fIsBuiltByAtomCount   = true;


    // NIST_MAT_BLOOD_ICRP -------------------------------------------------------------------
    static const int ElemsZ_117[]        = {  1,   6,   7,   8,  11,  15,
                                             16,  17,  19,  26};
    static const double FractionsZ_117[] = {     0.102,       0.11,      0.033,
                                                 0.745,      0.001,      0.001,
                                                 0.002,      0.003,      0.002,
                                                 0.001};
    fNISTMaterialDataTable[117].fName                 = "NIST_MAT_BLOOD_ICRP";
    fNISTMaterialDataTable[117].fDensity              = 1.06 * (g/cm3);
    fNISTMaterialDataTable[117].fMeanExcitationEnergy = 75.2 * eV;
    fNISTMaterialDataTable[117].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[117].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[117].fNumComponents        = 10;
    fNISTMaterialDataTable[117].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[117].fElementList          = ElemsZ_117;
    fNISTMaterialDataTable[117].fElementFraction      = FractionsZ_117;
    fNISTMaterialDataTable[117].fIsBuiltByAtomCount   = false;


    // NIST_MAT_BONE_COMPACT_ICRU ------------------------------------------------------------
    static const int ElemsZ_118[]        = {  1,   6,   7,   8,  12,  15,
                                             16,  20};
    static const double FractionsZ_118[] = {     0.064,      0.278,      0.027,
                                                  0.41,      0.002,       0.07,
                                                 0.002,      0.147};
    fNISTMaterialDataTable[118].fName                 = "NIST_MAT_BONE_COMPACT_ICRU";
    fNISTMaterialDataTable[118].fDensity              = 1.85 * (g/cm3);
    fNISTMaterialDataTable[118].fMeanExcitationEnergy = 91.9 * eV;
    fNISTMaterialDataTable[118].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[118].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[118].fNumComponents        = 8;
    fNISTMaterialDataTable[118].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[118].fElementList          = ElemsZ_118;
    fNISTMaterialDataTable[118].fElementFraction      = FractionsZ_118;
    fNISTMaterialDataTable[118].fIsBuiltByAtomCount   = false;


    // NIST_MAT_BONE_CORTICAL_ICRP -----------------------------------------------------------
    // Sceleton Cortical bone for Adult ICRU 46
    static const int ElemsZ_119[]        = {  1,   6,   7,   8,  11,  12,
                                             15,  16,  20};
    static const double FractionsZ_119[] = {     0.034,      0.155,      0.042,
                                                 0.435,      0.001,      0.002,
                                                 0.103,      0.003,      0.225};
    fNISTMaterialDataTable[119].fName                 = "NIST_MAT_BONE_CORTICAL_ICRP";
    fNISTMaterialDataTable[119].fDensity              = 1.92 * (g/cm3);
    fNISTMaterialDataTable[119].fMeanExcitationEnergy = 110 * eV;
    fNISTMaterialDataTable[119].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[119].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[119].fNumComponents        = 9;
    fNISTMaterialDataTable[119].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[119].fElementList          = ElemsZ_119;
    fNISTMaterialDataTable[119].fElementFraction      = FractionsZ_119;
    fNISTMaterialDataTable[119].fIsBuiltByAtomCount   = false;


    // NIST_MAT_BORON_CARBIDE ----------------------------------------------------------------
    static const int ElemsZ_120[]        = {  5,   6};
    static const double FractionsZ_120[] = {  4,   1};
    fNISTMaterialDataTable[120].fName                 = "NIST_MAT_BORON_CARBIDE";
    fNISTMaterialDataTable[120].fDensity              = 2.52 * (g/cm3);
    fNISTMaterialDataTable[120].fMeanExcitationEnergy = 84.7 * eV;
    fNISTMaterialDataTable[120].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[120].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[120].fNumComponents        = 2;
    fNISTMaterialDataTable[120].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[120].fElementList          = ElemsZ_120;
    fNISTMaterialDataTable[120].fElementFraction      = FractionsZ_120;
    fNISTMaterialDataTable[120].fIsBuiltByAtomCount   = true;


    // NIST_MAT_BORON_OXIDE ------------------------------------------------------------------
    static const int ElemsZ_121[]        = {  5,   8};
    static const double FractionsZ_121[] = {  2,   3};
    fNISTMaterialDataTable[121].fName                 = "NIST_MAT_BORON_OXIDE";
    fNISTMaterialDataTable[121].fDensity              = 1.812 * (g/cm3);
    fNISTMaterialDataTable[121].fMeanExcitationEnergy = 99.6 * eV;
    fNISTMaterialDataTable[121].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[121].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[121].fNumComponents        = 2;
    fNISTMaterialDataTable[121].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[121].fElementList          = ElemsZ_121;
    fNISTMaterialDataTable[121].fElementFraction      = FractionsZ_121;
    fNISTMaterialDataTable[121].fIsBuiltByAtomCount   = true;


    // NIST_MAT_BRAIN_ICRP -------------------------------------------------------------------
    static const int ElemsZ_122[]        = {  1,   6,   7,   8,  11,  15,
                                             16,  17,  19};
    static const double FractionsZ_122[] = {     0.107,      0.145,      0.022,
                                                 0.712,      0.002,      0.004,
                                                 0.002,      0.003,      0.003};
    fNISTMaterialDataTable[122].fName                 = "NIST_MAT_BRAIN_ICRP";
    fNISTMaterialDataTable[122].fDensity              = 1.04 * (g/cm3);
    fNISTMaterialDataTable[122].fMeanExcitationEnergy = 73.3 * eV;
    fNISTMaterialDataTable[122].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[122].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[122].fNumComponents        = 9;
    fNISTMaterialDataTable[122].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[122].fElementList          = ElemsZ_122;
    fNISTMaterialDataTable[122].fElementFraction      = FractionsZ_122;
    fNISTMaterialDataTable[122].fIsBuiltByAtomCount   = false;


    // NIST_MAT_BUTANE -----------------------------------------------------------------------
    static const int ElemsZ_123[]        = {  6,   1};
    static const double FractionsZ_123[] = {  4,  10};
    fNISTMaterialDataTable[123].fName                 = "NIST_MAT_BUTANE";
    fNISTMaterialDataTable[123].fDensity              = 0.00249343 * (g/cm3);
    fNISTMaterialDataTable[123].fMeanExcitationEnergy = 48.3 * eV;
    fNISTMaterialDataTable[123].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[123].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[123].fNumComponents        = 2;
    fNISTMaterialDataTable[123].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[123].fElementList          = ElemsZ_123;
    fNISTMaterialDataTable[123].fElementFraction      = FractionsZ_123;
    fNISTMaterialDataTable[123].fIsBuiltByAtomCount   = true;


    // NIST_MAT_N-BUTYL_ALCOHOL --------------------------------------------------------------
    static const int ElemsZ_124[]        = {  6,   1,   8};
    static const double FractionsZ_124[] = {  4,  10,   1};
    fNISTMaterialDataTable[124].fName                 = "NIST_MAT_N-BUTYL_ALCOHOL";
    fNISTMaterialDataTable[124].fDensity              = 0.8098 * (g/cm3);
    fNISTMaterialDataTable[124].fMeanExcitationEnergy = 59.9 * eV;
    fNISTMaterialDataTable[124].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[124].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[124].fNumComponents        = 3;
    fNISTMaterialDataTable[124].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[124].fElementList          = ElemsZ_124;
    fNISTMaterialDataTable[124].fElementFraction      = FractionsZ_124;
    fNISTMaterialDataTable[124].fIsBuiltByAtomCount   = true;


    // NIST_MAT_C-552 ------------------------------------------------------------------------
    static const int ElemsZ_125[]        = {  1,   6,   8,   9,  14};
    static const double FractionsZ_125[] = {   0.02468,    0.50161,   0.004527,
                                              0.465209,   0.003973};
    fNISTMaterialDataTable[125].fName                 = "NIST_MAT_C-552";
    fNISTMaterialDataTable[125].fDensity              = 1.76 * (g/cm3);
    fNISTMaterialDataTable[125].fMeanExcitationEnergy = 86.8 * eV;
    fNISTMaterialDataTable[125].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[125].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[125].fNumComponents        = 5;
    fNISTMaterialDataTable[125].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[125].fElementList          = ElemsZ_125;
    fNISTMaterialDataTable[125].fElementFraction      = FractionsZ_125;
    fNISTMaterialDataTable[125].fIsBuiltByAtomCount   = false;


    // NIST_MAT_CADMIUM_TELLURIDE ------------------------------------------------------------
    static const int ElemsZ_126[]        = { 48,  52};
    static const double FractionsZ_126[] = {  1,   1};
    fNISTMaterialDataTable[126].fName                 = "NIST_MAT_CADMIUM_TELLURIDE";
    fNISTMaterialDataTable[126].fDensity              = 6.2 * (g/cm3);
    fNISTMaterialDataTable[126].fMeanExcitationEnergy = 539.3 * eV;
    fNISTMaterialDataTable[126].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[126].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[126].fNumComponents        = 2;
    fNISTMaterialDataTable[126].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[126].fElementList          = ElemsZ_126;
    fNISTMaterialDataTable[126].fElementFraction      = FractionsZ_126;
    fNISTMaterialDataTable[126].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CADMIUM_TUNGSTATE ------------------------------------------------------------
    static const int ElemsZ_127[]        = { 48,  74,   8};
    static const double FractionsZ_127[] = {  1,   1,   4};
    fNISTMaterialDataTable[127].fName                 = "NIST_MAT_CADMIUM_TUNGSTATE";
    fNISTMaterialDataTable[127].fDensity              = 7.9 * (g/cm3);
    fNISTMaterialDataTable[127].fMeanExcitationEnergy = 468.3 * eV;
    fNISTMaterialDataTable[127].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[127].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[127].fNumComponents        = 3;
    fNISTMaterialDataTable[127].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[127].fElementList          = ElemsZ_127;
    fNISTMaterialDataTable[127].fElementFraction      = FractionsZ_127;
    fNISTMaterialDataTable[127].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CALCIUM_CARBONATE ------------------------------------------------------------
    static const int ElemsZ_128[]        = { 20,   6,   8};
    static const double FractionsZ_128[] = {  1,   1,   3};
    fNISTMaterialDataTable[128].fName                 = "NIST_MAT_CALCIUM_CARBONATE";
    fNISTMaterialDataTable[128].fDensity              = 2.8 * (g/cm3);
    fNISTMaterialDataTable[128].fMeanExcitationEnergy = 136.4 * eV;
    fNISTMaterialDataTable[128].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[128].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[128].fNumComponents        = 3;
    fNISTMaterialDataTable[128].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[128].fElementList          = ElemsZ_128;
    fNISTMaterialDataTable[128].fElementFraction      = FractionsZ_128;
    fNISTMaterialDataTable[128].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CALCIUM_FLUORIDE -------------------------------------------------------------
    static const int ElemsZ_129[]        = { 20,   9};
    static const double FractionsZ_129[] = {  1,   2};
    fNISTMaterialDataTable[129].fName                 = "NIST_MAT_CALCIUM_FLUORIDE";
    fNISTMaterialDataTable[129].fDensity              = 3.18 * (g/cm3);
    fNISTMaterialDataTable[129].fMeanExcitationEnergy = 166 * eV;
    fNISTMaterialDataTable[129].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[129].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[129].fNumComponents        = 2;
    fNISTMaterialDataTable[129].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[129].fElementList          = ElemsZ_129;
    fNISTMaterialDataTable[129].fElementFraction      = FractionsZ_129;
    fNISTMaterialDataTable[129].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CALCIUM_OXIDE ----------------------------------------------------------------
    static const int ElemsZ_130[]        = { 20,   8};
    static const double FractionsZ_130[] = {  1,   1};
    fNISTMaterialDataTable[130].fName                 = "NIST_MAT_CALCIUM_OXIDE";
    fNISTMaterialDataTable[130].fDensity              = 3.3 * (g/cm3);
    fNISTMaterialDataTable[130].fMeanExcitationEnergy = 176.1 * eV;
    fNISTMaterialDataTable[130].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[130].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[130].fNumComponents        = 2;
    fNISTMaterialDataTable[130].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[130].fElementList          = ElemsZ_130;
    fNISTMaterialDataTable[130].fElementFraction      = FractionsZ_130;
    fNISTMaterialDataTable[130].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CALCIUM_SULFATE --------------------------------------------------------------
    static const int ElemsZ_131[]        = { 20,  16,   8};
    static const double FractionsZ_131[] = {  1,   1,   4};
    fNISTMaterialDataTable[131].fName                 = "NIST_MAT_CALCIUM_SULFATE";
    fNISTMaterialDataTable[131].fDensity              = 2.96 * (g/cm3);
    fNISTMaterialDataTable[131].fMeanExcitationEnergy = 152.3 * eV;
    fNISTMaterialDataTable[131].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[131].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[131].fNumComponents        = 3;
    fNISTMaterialDataTable[131].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[131].fElementList          = ElemsZ_131;
    fNISTMaterialDataTable[131].fElementFraction      = FractionsZ_131;
    fNISTMaterialDataTable[131].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CALCIUM_TUNGSTATE ------------------------------------------------------------
    static const int ElemsZ_132[]        = { 20,  74,   8};
    static const double FractionsZ_132[] = {  1,   1,   4};
    fNISTMaterialDataTable[132].fName                 = "NIST_MAT_CALCIUM_TUNGSTATE";
    fNISTMaterialDataTable[132].fDensity              = 6.062 * (g/cm3);
    fNISTMaterialDataTable[132].fMeanExcitationEnergy = 395 * eV;
    fNISTMaterialDataTable[132].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[132].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[132].fNumComponents        = 3;
    fNISTMaterialDataTable[132].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[132].fElementList          = ElemsZ_132;
    fNISTMaterialDataTable[132].fElementFraction      = FractionsZ_132;
    fNISTMaterialDataTable[132].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CARBON_DIOXIDE ---------------------------------------------------------------
    static const int ElemsZ_133[]        = {  6,   8};
    static const double FractionsZ_133[] = {  1,   2};
    fNISTMaterialDataTable[133].fName                 = "NIST_MAT_CARBON_DIOXIDE";
    fNISTMaterialDataTable[133].fDensity              = 0.00184212 * (g/cm3);
    fNISTMaterialDataTable[133].fMeanExcitationEnergy = 85 * eV;
    fNISTMaterialDataTable[133].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[133].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[133].fNumComponents        = 2;
    fNISTMaterialDataTable[133].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[133].fElementList          = ElemsZ_133;
    fNISTMaterialDataTable[133].fElementFraction      = FractionsZ_133;
    fNISTMaterialDataTable[133].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CARBON_TETRACHLORIDE ---------------------------------------------------------
    static const int ElemsZ_134[]        = {  6,  17};
    static const double FractionsZ_134[] = {  1,   4};
    fNISTMaterialDataTable[134].fName                 = "NIST_MAT_CARBON_TETRACHLORIDE";
    fNISTMaterialDataTable[134].fDensity              = 1.594 * (g/cm3);
    fNISTMaterialDataTable[134].fMeanExcitationEnergy = 166.3 * eV;
    fNISTMaterialDataTable[134].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[134].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[134].fNumComponents        = 2;
    fNISTMaterialDataTable[134].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[134].fElementList          = ElemsZ_134;
    fNISTMaterialDataTable[134].fElementFraction      = FractionsZ_134;
    fNISTMaterialDataTable[134].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CELLULOSE_CELLOPHANE ---------------------------------------------------------
    static const int ElemsZ_135[]        = {  6,   1,   8};
    static const double FractionsZ_135[] = {  6,  10,   5};
    fNISTMaterialDataTable[135].fName                 = "NIST_MAT_CELLULOSE_CELLOPHANE";
    fNISTMaterialDataTable[135].fDensity              = 1.42 * (g/cm3);
    fNISTMaterialDataTable[135].fMeanExcitationEnergy = 77.6 * eV;
    fNISTMaterialDataTable[135].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[135].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[135].fNumComponents        = 3;
    fNISTMaterialDataTable[135].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[135].fElementList          = ElemsZ_135;
    fNISTMaterialDataTable[135].fElementFraction      = FractionsZ_135;
    fNISTMaterialDataTable[135].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CELLULOSE_BUTYRATE -----------------------------------------------------------
    static const int ElemsZ_136[]        = {  1,   6,   8};
    static const double FractionsZ_136[] = {  0.067125,   0.545403,   0.387472};
    fNISTMaterialDataTable[136].fName                 = "NIST_MAT_CELLULOSE_BUTYRATE";
    fNISTMaterialDataTable[136].fDensity              = 1.2 * (g/cm3);
    fNISTMaterialDataTable[136].fMeanExcitationEnergy = 74.6 * eV;
    fNISTMaterialDataTable[136].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[136].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[136].fNumComponents        = 3;
    fNISTMaterialDataTable[136].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[136].fElementList          = ElemsZ_136;
    fNISTMaterialDataTable[136].fElementFraction      = FractionsZ_136;
    fNISTMaterialDataTable[136].fIsBuiltByAtomCount   = false;


    // NIST_MAT_CELLULOSE_NITRATE ------------------------------------------------------------
    static const int ElemsZ_137[]        = {  1,   6,   7,   8};
    static const double FractionsZ_137[] = {  0.029216,   0.271296,   0.121276,
                                              0.578212};
    fNISTMaterialDataTable[137].fName                 = "NIST_MAT_CELLULOSE_NITRATE";
    fNISTMaterialDataTable[137].fDensity              = 1.49 * (g/cm3);
    fNISTMaterialDataTable[137].fMeanExcitationEnergy = 87 * eV;
    fNISTMaterialDataTable[137].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[137].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[137].fNumComponents        = 4;
    fNISTMaterialDataTable[137].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[137].fElementList          = ElemsZ_137;
    fNISTMaterialDataTable[137].fElementFraction      = FractionsZ_137;
    fNISTMaterialDataTable[137].fIsBuiltByAtomCount   = false;


    // NIST_MAT_CERIC_SULFATE ----------------------------------------------------------------
    static const int ElemsZ_138[]        = {  1,   7,   8,  16,  58};
    static const double FractionsZ_138[] = {  0.107596,     0.0008,   0.874976,
                                              0.014627,   0.002001};
    fNISTMaterialDataTable[138].fName                 = "NIST_MAT_CERIC_SULFATE";
    fNISTMaterialDataTable[138].fDensity              = 1.03 * (g/cm3);
    fNISTMaterialDataTable[138].fMeanExcitationEnergy = 76.7 * eV;
    fNISTMaterialDataTable[138].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[138].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[138].fNumComponents        = 5;
    fNISTMaterialDataTable[138].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[138].fElementList          = ElemsZ_138;
    fNISTMaterialDataTable[138].fElementFraction      = FractionsZ_138;
    fNISTMaterialDataTable[138].fIsBuiltByAtomCount   = false;


    // NIST_MAT_CESIUM_FLUORIDE --------------------------------------------------------------
    static const int ElemsZ_139[]        = { 55,   9};
    static const double FractionsZ_139[] = {  1,   1};
    fNISTMaterialDataTable[139].fName                 = "NIST_MAT_CESIUM_FLUORIDE";
    fNISTMaterialDataTable[139].fDensity              = 4.115 * (g/cm3);
    fNISTMaterialDataTable[139].fMeanExcitationEnergy = 440.7 * eV;
    fNISTMaterialDataTable[139].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[139].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[139].fNumComponents        = 2;
    fNISTMaterialDataTable[139].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[139].fElementList          = ElemsZ_139;
    fNISTMaterialDataTable[139].fElementFraction      = FractionsZ_139;
    fNISTMaterialDataTable[139].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CESIUM_IODIDE ----------------------------------------------------------------
    static const int ElemsZ_140[]        = { 55,  53};
    static const double FractionsZ_140[] = {  1,   1};
    fNISTMaterialDataTable[140].fName                 = "NIST_MAT_CESIUM_IODIDE";
    fNISTMaterialDataTable[140].fDensity              = 4.51 * (g/cm3);
    fNISTMaterialDataTable[140].fMeanExcitationEnergy = 553.1 * eV;
    fNISTMaterialDataTable[140].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[140].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[140].fNumComponents        = 2;
    fNISTMaterialDataTable[140].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[140].fElementList          = ElemsZ_140;
    fNISTMaterialDataTable[140].fElementFraction      = FractionsZ_140;
    fNISTMaterialDataTable[140].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CHLOROBENZENE ----------------------------------------------------------------
    static const int ElemsZ_141[]        = {  6,   1,  17};
    static const double FractionsZ_141[] = {  6,   5,   1};
    fNISTMaterialDataTable[141].fName                 = "NIST_MAT_CHLOROBENZENE";
    fNISTMaterialDataTable[141].fDensity              = 1.1058 * (g/cm3);
    fNISTMaterialDataTable[141].fMeanExcitationEnergy = 89.1 * eV;
    fNISTMaterialDataTable[141].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[141].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[141].fNumComponents        = 3;
    fNISTMaterialDataTable[141].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[141].fElementList          = ElemsZ_141;
    fNISTMaterialDataTable[141].fElementFraction      = FractionsZ_141;
    fNISTMaterialDataTable[141].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CHLOROFORM -------------------------------------------------------------------
    static const int ElemsZ_142[]        = {  6,   1,  17};
    static const double FractionsZ_142[] = {  1,   1,   3};
    fNISTMaterialDataTable[142].fName                 = "NIST_MAT_CHLOROFORM";
    fNISTMaterialDataTable[142].fDensity              = 1.4832 * (g/cm3);
    fNISTMaterialDataTable[142].fMeanExcitationEnergy = 156 * eV;
    fNISTMaterialDataTable[142].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[142].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[142].fNumComponents        = 3;
    fNISTMaterialDataTable[142].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[142].fElementList          = ElemsZ_142;
    fNISTMaterialDataTable[142].fElementFraction      = FractionsZ_142;
    fNISTMaterialDataTable[142].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CONCRETE ---------------------------------------------------------------------
    static const int ElemsZ_143[]        = {  1,   6,   8,  11,  12,  13,
                                             14,  19,  20,  26};
    static const double FractionsZ_143[] = {      0.01,      0.001,   0.529107,
                                                 0.016,      0.002,   0.033872,
                                              0.337021,      0.013,      0.044,
                                                 0.014};
    fNISTMaterialDataTable[143].fName                 = "NIST_MAT_CONCRETE";
    fNISTMaterialDataTable[143].fDensity              = 2.3 * (g/cm3);
    fNISTMaterialDataTable[143].fMeanExcitationEnergy = 135.2 * eV;
    fNISTMaterialDataTable[143].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[143].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[143].fNumComponents        = 10;
    fNISTMaterialDataTable[143].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[143].fElementList          = ElemsZ_143;
    fNISTMaterialDataTable[143].fElementFraction      = FractionsZ_143;
    fNISTMaterialDataTable[143].fIsBuiltByAtomCount   = false;


    // NIST_MAT_CYCLOHEXANE ------------------------------------------------------------------
    static const int ElemsZ_144[]        = {  6,   1};
    static const double FractionsZ_144[] = {  6,  12};
    fNISTMaterialDataTable[144].fName                 = "NIST_MAT_CYCLOHEXANE";
    fNISTMaterialDataTable[144].fDensity              = 0.779 * (g/cm3);
    fNISTMaterialDataTable[144].fMeanExcitationEnergy = 56.4 * eV;
    fNISTMaterialDataTable[144].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[144].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[144].fNumComponents        = 2;
    fNISTMaterialDataTable[144].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[144].fElementList          = ElemsZ_144;
    fNISTMaterialDataTable[144].fElementFraction      = FractionsZ_144;
    fNISTMaterialDataTable[144].fIsBuiltByAtomCount   = true;


    // NIST_MAT_2-DICHLOROBENZENE ------------------------------------------------------------
    static const int ElemsZ_145[]        = {  6,   1,  17};
    static const double FractionsZ_145[] = {  6,   4,   2};
    fNISTMaterialDataTable[145].fName                 = "NIST_MAT_2-DICHLOROBENZENE";
    fNISTMaterialDataTable[145].fDensity              = 1.3048 * (g/cm3);
    fNISTMaterialDataTable[145].fMeanExcitationEnergy = 106.5 * eV;
    fNISTMaterialDataTable[145].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[145].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[145].fNumComponents        = 3;
    fNISTMaterialDataTable[145].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[145].fElementList          = ElemsZ_145;
    fNISTMaterialDataTable[145].fElementFraction      = FractionsZ_145;
    fNISTMaterialDataTable[145].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DICHLORODIETHYL_ETHER --------------------------------------------------------
    static const int ElemsZ_146[]        = {  6,   1,   8,  17};
    static const double FractionsZ_146[] = {  4,   8,   1,   2};
    fNISTMaterialDataTable[146].fName                 = "NIST_MAT_DICHLORODIETHYL_ETHER";
    fNISTMaterialDataTable[146].fDensity              = 1.2199 * (g/cm3);
    fNISTMaterialDataTable[146].fMeanExcitationEnergy = 103.3 * eV;
    fNISTMaterialDataTable[146].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[146].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[146].fNumComponents        = 4;
    fNISTMaterialDataTable[146].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[146].fElementList          = ElemsZ_146;
    fNISTMaterialDataTable[146].fElementFraction      = FractionsZ_146;
    fNISTMaterialDataTable[146].fIsBuiltByAtomCount   = true;


    // NIST_MAT_2-DICHLOROETHANE -------------------------------------------------------------
    static const int ElemsZ_147[]        = {  6,   1,  17};
    static const double FractionsZ_147[] = {  2,   4,   2};
    fNISTMaterialDataTable[147].fName                 = "NIST_MAT_2-DICHLOROETHANE";
    fNISTMaterialDataTable[147].fDensity              = 1.2351 * (g/cm3);
    fNISTMaterialDataTable[147].fMeanExcitationEnergy = 111.9 * eV;
    fNISTMaterialDataTable[147].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[147].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[147].fNumComponents        = 3;
    fNISTMaterialDataTable[147].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[147].fElementList          = ElemsZ_147;
    fNISTMaterialDataTable[147].fElementFraction      = FractionsZ_147;
    fNISTMaterialDataTable[147].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DIETHYL_ETHER ----------------------------------------------------------------
    static const int ElemsZ_148[]        = {  6,   1,   8};
    static const double FractionsZ_148[] = {  4,  10,   1};
    fNISTMaterialDataTable[148].fName                 = "NIST_MAT_DIETHYL_ETHER";
    fNISTMaterialDataTable[148].fDensity              = 0.71378 * (g/cm3);
    fNISTMaterialDataTable[148].fMeanExcitationEnergy = 60 * eV;
    fNISTMaterialDataTable[148].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[148].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[148].fNumComponents        = 3;
    fNISTMaterialDataTable[148].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[148].fElementList          = ElemsZ_148;
    fNISTMaterialDataTable[148].fElementFraction      = FractionsZ_148;
    fNISTMaterialDataTable[148].fIsBuiltByAtomCount   = true;


    // NIST_MAT_N-DIMETHYL_FORMAMIDE ---------------------------------------------------------
    static const int ElemsZ_149[]        = {  6,   1,   7,   8};
    static const double FractionsZ_149[] = {  3,   7,   1,   1};
    fNISTMaterialDataTable[149].fName                 = "NIST_MAT_N-DIMETHYL_FORMAMIDE";
    fNISTMaterialDataTable[149].fDensity              = 0.9487 * (g/cm3);
    fNISTMaterialDataTable[149].fMeanExcitationEnergy = 66.6 * eV;
    fNISTMaterialDataTable[149].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[149].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[149].fNumComponents        = 4;
    fNISTMaterialDataTable[149].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[149].fElementList          = ElemsZ_149;
    fNISTMaterialDataTable[149].fElementFraction      = FractionsZ_149;
    fNISTMaterialDataTable[149].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DIMETHYL_SULFOXIDE -----------------------------------------------------------
    static const int ElemsZ_150[]        = {  6,   1,   8,  16};
    static const double FractionsZ_150[] = {  2,   6,   1,   1};
    fNISTMaterialDataTable[150].fName                 = "NIST_MAT_DIMETHYL_SULFOXIDE";
    fNISTMaterialDataTable[150].fDensity              = 1.1014 * (g/cm3);
    fNISTMaterialDataTable[150].fMeanExcitationEnergy = 98.6 * eV;
    fNISTMaterialDataTable[150].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[150].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[150].fNumComponents        = 4;
    fNISTMaterialDataTable[150].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[150].fElementList          = ElemsZ_150;
    fNISTMaterialDataTable[150].fElementFraction      = FractionsZ_150;
    fNISTMaterialDataTable[150].fIsBuiltByAtomCount   = true;


    // NIST_MAT_ETHANE -----------------------------------------------------------------------
    static const int ElemsZ_151[]        = {  6,   1};
    static const double FractionsZ_151[] = {  2,   6};
    fNISTMaterialDataTable[151].fName                 = "NIST_MAT_ETHANE";
    fNISTMaterialDataTable[151].fDensity              = 0.00125324 * (g/cm3);
    fNISTMaterialDataTable[151].fMeanExcitationEnergy = 45.4 * eV;
    fNISTMaterialDataTable[151].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[151].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[151].fNumComponents        = 2;
    fNISTMaterialDataTable[151].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[151].fElementList          = ElemsZ_151;
    fNISTMaterialDataTable[151].fElementFraction      = FractionsZ_151;
    fNISTMaterialDataTable[151].fIsBuiltByAtomCount   = true;


    // NIST_MAT_ETHYL_ALCOHOL ----------------------------------------------------------------
    static const int ElemsZ_152[]        = {  6,   1,   8};
    static const double FractionsZ_152[] = {  2,   6,   1};
    fNISTMaterialDataTable[152].fName                 = "NIST_MAT_ETHYL_ALCOHOL";
    fNISTMaterialDataTable[152].fDensity              = 0.7893 * (g/cm3);
    fNISTMaterialDataTable[152].fMeanExcitationEnergy = 62.9 * eV;
    fNISTMaterialDataTable[152].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[152].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[152].fNumComponents        = 3;
    fNISTMaterialDataTable[152].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[152].fElementList          = ElemsZ_152;
    fNISTMaterialDataTable[152].fElementFraction      = FractionsZ_152;
    fNISTMaterialDataTable[152].fIsBuiltByAtomCount   = true;


    // NIST_MAT_ETHYL_CELLULOSE --------------------------------------------------------------
    static const int ElemsZ_153[]        = {  1,   6,   8};
    static const double FractionsZ_153[] = {  0.090027,   0.585182,   0.324791};
    fNISTMaterialDataTable[153].fName                 = "NIST_MAT_ETHYL_CELLULOSE";
    fNISTMaterialDataTable[153].fDensity              = 1.13 * (g/cm3);
    fNISTMaterialDataTable[153].fMeanExcitationEnergy = 69.3 * eV;
    fNISTMaterialDataTable[153].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[153].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[153].fNumComponents        = 3;
    fNISTMaterialDataTable[153].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[153].fElementList          = ElemsZ_153;
    fNISTMaterialDataTable[153].fElementFraction      = FractionsZ_153;
    fNISTMaterialDataTable[153].fIsBuiltByAtomCount   = false;


    // NIST_MAT_ETHYLENE ---------------------------------------------------------------------
    static const int ElemsZ_154[]        = {  6,   1};
    static const double FractionsZ_154[] = {  2,   4};
    fNISTMaterialDataTable[154].fName                 = "NIST_MAT_ETHYLENE";
    fNISTMaterialDataTable[154].fDensity              = 0.00117497 * (g/cm3);
    fNISTMaterialDataTable[154].fMeanExcitationEnergy = 50.7 * eV;
    fNISTMaterialDataTable[154].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[154].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[154].fNumComponents        = 2;
    fNISTMaterialDataTable[154].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[154].fElementList          = ElemsZ_154;
    fNISTMaterialDataTable[154].fElementFraction      = FractionsZ_154;
    fNISTMaterialDataTable[154].fIsBuiltByAtomCount   = true;


    // NIST_MAT_EYE_LENS_ICRP ----------------------------------------------------------------
    static const int ElemsZ_155[]        = {  1,   6,   7,   8,  11,  15,
                                             16,  17};
    static const double FractionsZ_155[] = {     0.096,      0.195,      0.057,
                                                 0.646,      0.001,      0.001,
                                                 0.003,      0.001};
    fNISTMaterialDataTable[155].fName                 = "NIST_MAT_EYE_LENS_ICRP";
    fNISTMaterialDataTable[155].fDensity              = 1.07 * (g/cm3);
    fNISTMaterialDataTable[155].fMeanExcitationEnergy = 73.3 * eV;
    fNISTMaterialDataTable[155].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[155].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[155].fNumComponents        = 8;
    fNISTMaterialDataTable[155].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[155].fElementList          = ElemsZ_155;
    fNISTMaterialDataTable[155].fElementFraction      = FractionsZ_155;
    fNISTMaterialDataTable[155].fIsBuiltByAtomCount   = false;


    // NIST_MAT_FERRIC_OXIDE -----------------------------------------------------------------
    static const int ElemsZ_156[]        = { 26,   8};
    static const double FractionsZ_156[] = {  2,   3};
    fNISTMaterialDataTable[156].fName                 = "NIST_MAT_FERRIC_OXIDE";
    fNISTMaterialDataTable[156].fDensity              = 5.2 * (g/cm3);
    fNISTMaterialDataTable[156].fMeanExcitationEnergy = 227.3 * eV;
    fNISTMaterialDataTable[156].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[156].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[156].fNumComponents        = 2;
    fNISTMaterialDataTable[156].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[156].fElementList          = ElemsZ_156;
    fNISTMaterialDataTable[156].fElementFraction      = FractionsZ_156;
    fNISTMaterialDataTable[156].fIsBuiltByAtomCount   = true;


    // NIST_MAT_FERROBORIDE ------------------------------------------------------------------
    static const int ElemsZ_157[]        = { 26,   5};
    static const double FractionsZ_157[] = {  1,   1};
    fNISTMaterialDataTable[157].fName                 = "NIST_MAT_FERROBORIDE";
    fNISTMaterialDataTable[157].fDensity              = 7.15 * (g/cm3);
    fNISTMaterialDataTable[157].fMeanExcitationEnergy = 261 * eV;
    fNISTMaterialDataTable[157].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[157].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[157].fNumComponents        = 2;
    fNISTMaterialDataTable[157].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[157].fElementList          = ElemsZ_157;
    fNISTMaterialDataTable[157].fElementFraction      = FractionsZ_157;
    fNISTMaterialDataTable[157].fIsBuiltByAtomCount   = true;


    // NIST_MAT_FERROUS_OXIDE ----------------------------------------------------------------
    static const int ElemsZ_158[]        = { 26,   8};
    static const double FractionsZ_158[] = {  1,   1};
    fNISTMaterialDataTable[158].fName                 = "NIST_MAT_FERROUS_OXIDE";
    fNISTMaterialDataTable[158].fDensity              = 5.7 * (g/cm3);
    fNISTMaterialDataTable[158].fMeanExcitationEnergy = 248.6 * eV;
    fNISTMaterialDataTable[158].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[158].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[158].fNumComponents        = 2;
    fNISTMaterialDataTable[158].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[158].fElementList          = ElemsZ_158;
    fNISTMaterialDataTable[158].fElementFraction      = FractionsZ_158;
    fNISTMaterialDataTable[158].fIsBuiltByAtomCount   = true;


    // NIST_MAT_FERROUS_SULFATE --------------------------------------------------------------
    static const int ElemsZ_159[]        = {  1,   7,   8,  11,  16,  17,
                                             26};
    static const double FractionsZ_159[] = {  0.108259,    2.7e-05,   0.878636,
                                               2.2e-05,   0.012968,    3.4e-05,
                                               5.4e-05};
    fNISTMaterialDataTable[159].fName                 = "NIST_MAT_FERROUS_SULFATE";
    fNISTMaterialDataTable[159].fDensity              = 1.024 * (g/cm3);
    fNISTMaterialDataTable[159].fMeanExcitationEnergy = 76.4 * eV;
    fNISTMaterialDataTable[159].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[159].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[159].fNumComponents        = 7;
    fNISTMaterialDataTable[159].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[159].fElementList          = ElemsZ_159;
    fNISTMaterialDataTable[159].fElementFraction      = FractionsZ_159;
    fNISTMaterialDataTable[159].fIsBuiltByAtomCount   = false;


    // NIST_MAT_FREON-12 ---------------------------------------------------------------------
    static const int ElemsZ_160[]        = {  6,   9,  17};
    static const double FractionsZ_160[] = {  0.099335,   0.314247,   0.586418};
    fNISTMaterialDataTable[160].fName                 = "NIST_MAT_FREON-12";
    fNISTMaterialDataTable[160].fDensity              = 1.12 * (g/cm3);
    fNISTMaterialDataTable[160].fMeanExcitationEnergy = 143 * eV;
    fNISTMaterialDataTable[160].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[160].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[160].fNumComponents        = 3;
    fNISTMaterialDataTable[160].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[160].fElementList          = ElemsZ_160;
    fNISTMaterialDataTable[160].fElementFraction      = FractionsZ_160;
    fNISTMaterialDataTable[160].fIsBuiltByAtomCount   = false;


    // NIST_MAT_FREON-12B2 -------------------------------------------------------------------
    static const int ElemsZ_161[]        = {  6,   9,  35};
    static const double FractionsZ_161[] = {  0.057245,   0.181096,   0.761659};
    fNISTMaterialDataTable[161].fName                 = "NIST_MAT_FREON-12B2";
    fNISTMaterialDataTable[161].fDensity              = 1.8 * (g/cm3);
    fNISTMaterialDataTable[161].fMeanExcitationEnergy = 284.9 * eV;
    fNISTMaterialDataTable[161].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[161].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[161].fNumComponents        = 3;
    fNISTMaterialDataTable[161].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[161].fElementList          = ElemsZ_161;
    fNISTMaterialDataTable[161].fElementFraction      = FractionsZ_161;
    fNISTMaterialDataTable[161].fIsBuiltByAtomCount   = false;


    // NIST_MAT_FREON-13 ---------------------------------------------------------------------
    static const int ElemsZ_162[]        = {  6,   9,  17};
    static const double FractionsZ_162[] = {  0.114983,   0.545622,   0.339396};
    fNISTMaterialDataTable[162].fName                 = "NIST_MAT_FREON-13";
    fNISTMaterialDataTable[162].fDensity              = 0.95 * (g/cm3);
    fNISTMaterialDataTable[162].fMeanExcitationEnergy = 126.6 * eV;
    fNISTMaterialDataTable[162].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[162].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[162].fNumComponents        = 3;
    fNISTMaterialDataTable[162].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[162].fElementList          = ElemsZ_162;
    fNISTMaterialDataTable[162].fElementFraction      = FractionsZ_162;
    fNISTMaterialDataTable[162].fIsBuiltByAtomCount   = false;


    // NIST_MAT_FREON-13B1 -------------------------------------------------------------------
    static const int ElemsZ_163[]        = {  6,   9,  35};
    static const double FractionsZ_163[] = {  1,   3,   1};
    fNISTMaterialDataTable[163].fName                 = "NIST_MAT_FREON-13B1";
    fNISTMaterialDataTable[163].fDensity              = 1.5 * (g/cm3);
    fNISTMaterialDataTable[163].fMeanExcitationEnergy = 210.5 * eV;
    fNISTMaterialDataTable[163].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[163].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[163].fNumComponents        = 3;
    fNISTMaterialDataTable[163].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[163].fElementList          = ElemsZ_163;
    fNISTMaterialDataTable[163].fElementFraction      = FractionsZ_163;
    fNISTMaterialDataTable[163].fIsBuiltByAtomCount   = true;


    // NIST_MAT_FREON-13I1 -------------------------------------------------------------------
    static const int ElemsZ_164[]        = {  6,   9,  53};
    static const double FractionsZ_164[] = {  0.061309,   0.290924,   0.647767};
    fNISTMaterialDataTable[164].fName                 = "NIST_MAT_FREON-13I1";
    fNISTMaterialDataTable[164].fDensity              = 1.8 * (g/cm3);
    fNISTMaterialDataTable[164].fMeanExcitationEnergy = 293.5 * eV;
    fNISTMaterialDataTable[164].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[164].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[164].fNumComponents        = 3;
    fNISTMaterialDataTable[164].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[164].fElementList          = ElemsZ_164;
    fNISTMaterialDataTable[164].fElementFraction      = FractionsZ_164;
    fNISTMaterialDataTable[164].fIsBuiltByAtomCount   = false;


    // NIST_MAT_GADOLINIUM_OXYSULFIDE --------------------------------------------------------
    static const int ElemsZ_165[]        = { 64,   8,  16};
    static const double FractionsZ_165[] = {  2,   2,   1};
    fNISTMaterialDataTable[165].fName                 = "NIST_MAT_GADOLINIUM_OXYSULFIDE";
    fNISTMaterialDataTable[165].fDensity              = 7.44 * (g/cm3);
    fNISTMaterialDataTable[165].fMeanExcitationEnergy = 493.3 * eV;
    fNISTMaterialDataTable[165].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[165].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[165].fNumComponents        = 3;
    fNISTMaterialDataTable[165].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[165].fElementList          = ElemsZ_165;
    fNISTMaterialDataTable[165].fElementFraction      = FractionsZ_165;
    fNISTMaterialDataTable[165].fIsBuiltByAtomCount   = true;


    // NIST_MAT_GALLIUM_ARSENIDE -------------------------------------------------------------
    static const int ElemsZ_166[]        = { 31,  33};
    static const double FractionsZ_166[] = {  1,   1};
    fNISTMaterialDataTable[166].fName                 = "NIST_MAT_GALLIUM_ARSENIDE";
    fNISTMaterialDataTable[166].fDensity              = 5.31 * (g/cm3);
    fNISTMaterialDataTable[166].fMeanExcitationEnergy = 384.9 * eV;
    fNISTMaterialDataTable[166].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[166].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[166].fNumComponents        = 2;
    fNISTMaterialDataTable[166].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[166].fElementList          = ElemsZ_166;
    fNISTMaterialDataTable[166].fElementFraction      = FractionsZ_166;
    fNISTMaterialDataTable[166].fIsBuiltByAtomCount   = true;


    // NIST_MAT_GEL_PHOTO_EMULSION -----------------------------------------------------------
    static const int ElemsZ_167[]        = {  1,   6,   7,   8,  16};
    static const double FractionsZ_167[] = {   0.08118,    0.41606,    0.11124,
                                               0.38064,    0.01088};
    fNISTMaterialDataTable[167].fName                 = "NIST_MAT_GEL_PHOTO_EMULSION";
    fNISTMaterialDataTable[167].fDensity              = 1.2914 * (g/cm3);
    fNISTMaterialDataTable[167].fMeanExcitationEnergy = 74.8 * eV;
    fNISTMaterialDataTable[167].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[167].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[167].fNumComponents        = 5;
    fNISTMaterialDataTable[167].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[167].fElementList          = ElemsZ_167;
    fNISTMaterialDataTable[167].fElementFraction      = FractionsZ_167;
    fNISTMaterialDataTable[167].fIsBuiltByAtomCount   = false;


    // NIST_MAT_Pyrex_Glass ------------------------------------------------------------------
    static const int ElemsZ_168[]        = {  5,   8,  11,  13,  14,  19};
    static const double FractionsZ_168[] = {  0.040064,   0.539562,   0.028191,
                                              0.011644,    0.37722,   0.003321};
    fNISTMaterialDataTable[168].fName                 = "NIST_MAT_Pyrex_Glass";
    fNISTMaterialDataTable[168].fDensity              = 2.23 * (g/cm3);
    fNISTMaterialDataTable[168].fMeanExcitationEnergy = 134 * eV;
    fNISTMaterialDataTable[168].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[168].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[168].fNumComponents        = 6;
    fNISTMaterialDataTable[168].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[168].fElementList          = ElemsZ_168;
    fNISTMaterialDataTable[168].fElementFraction      = FractionsZ_168;
    fNISTMaterialDataTable[168].fIsBuiltByAtomCount   = false;


    // NIST_MAT_GLASS_LEAD -------------------------------------------------------------------
    static const int ElemsZ_169[]        = {  8,  14,  22,  33,  82};
    static const double FractionsZ_169[] = {  0.156453,   0.080866,   0.008092,
                                              0.002651,   0.751938};
    fNISTMaterialDataTable[169].fName                 = "NIST_MAT_GLASS_LEAD";
    fNISTMaterialDataTable[169].fDensity              = 6.22 * (g/cm3);
    fNISTMaterialDataTable[169].fMeanExcitationEnergy = 526.4 * eV;
    fNISTMaterialDataTable[169].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[169].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[169].fNumComponents        = 5;
    fNISTMaterialDataTable[169].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[169].fElementList          = ElemsZ_169;
    fNISTMaterialDataTable[169].fElementFraction      = FractionsZ_169;
    fNISTMaterialDataTable[169].fIsBuiltByAtomCount   = false;


    // NIST_MAT_GLASS_PLATE ------------------------------------------------------------------
    static const int ElemsZ_170[]        = {  8,  11,  14,  20};
    static const double FractionsZ_170[] = {    0.4598,   0.096441,   0.336553,
                                              0.107205};
    fNISTMaterialDataTable[170].fName                 = "NIST_MAT_GLASS_PLATE";
    fNISTMaterialDataTable[170].fDensity              = 2.4 * (g/cm3);
    fNISTMaterialDataTable[170].fMeanExcitationEnergy = 145.4 * eV;
    fNISTMaterialDataTable[170].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[170].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[170].fNumComponents        = 4;
    fNISTMaterialDataTable[170].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[170].fElementList          = ElemsZ_170;
    fNISTMaterialDataTable[170].fElementFraction      = FractionsZ_170;
    fNISTMaterialDataTable[170].fIsBuiltByAtomCount   = false;


    // NIST_MAT_GLUTAMINE --------------------------------------------------------------------
    static const int ElemsZ_171[]        = {  6,   1,   7,   8};
    static const double FractionsZ_171[] = {  5,  10,   2,   3};
    fNISTMaterialDataTable[171].fName                 = "NIST_MAT_GLUTAMINE";
    fNISTMaterialDataTable[171].fDensity              = 1.46 * (g/cm3);
    fNISTMaterialDataTable[171].fMeanExcitationEnergy = 73.3 * eV;
    fNISTMaterialDataTable[171].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[171].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[171].fNumComponents        = 4;
    fNISTMaterialDataTable[171].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[171].fElementList          = ElemsZ_171;
    fNISTMaterialDataTable[171].fElementFraction      = FractionsZ_171;
    fNISTMaterialDataTable[171].fIsBuiltByAtomCount   = true;


    // NIST_MAT_GLYCEROL ---------------------------------------------------------------------
    static const int ElemsZ_172[]        = {  6,   1,   8};
    static const double FractionsZ_172[] = {  3,   8,   3};
    fNISTMaterialDataTable[172].fName                 = "NIST_MAT_GLYCEROL";
    fNISTMaterialDataTable[172].fDensity              = 1.2613 * (g/cm3);
    fNISTMaterialDataTable[172].fMeanExcitationEnergy = 72.6 * eV;
    fNISTMaterialDataTable[172].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[172].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[172].fNumComponents        = 3;
    fNISTMaterialDataTable[172].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[172].fElementList          = ElemsZ_172;
    fNISTMaterialDataTable[172].fElementFraction      = FractionsZ_172;
    fNISTMaterialDataTable[172].fIsBuiltByAtomCount   = true;


    // NIST_MAT_GUANINE ----------------------------------------------------------------------
    static const int ElemsZ_173[]        = {  6,   1,   7,   8};
    static const double FractionsZ_173[] = {  5,   5,   5,   1};
    fNISTMaterialDataTable[173].fName                 = "NIST_MAT_GUANINE";
    fNISTMaterialDataTable[173].fDensity              = 2.2 * (g/cm3);
    fNISTMaterialDataTable[173].fMeanExcitationEnergy = 75 * eV;
    fNISTMaterialDataTable[173].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[173].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[173].fNumComponents        = 4;
    fNISTMaterialDataTable[173].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[173].fElementList          = ElemsZ_173;
    fNISTMaterialDataTable[173].fElementFraction      = FractionsZ_173;
    fNISTMaterialDataTable[173].fIsBuiltByAtomCount   = true;


    // NIST_MAT_GYPSUM -----------------------------------------------------------------------
    static const int ElemsZ_174[]        = { 20,  16,   8,   1};
    static const double FractionsZ_174[] = {  1,   1,   6,   4};
    fNISTMaterialDataTable[174].fName                 = "NIST_MAT_GYPSUM";
    fNISTMaterialDataTable[174].fDensity              = 2.32 * (g/cm3);
    fNISTMaterialDataTable[174].fMeanExcitationEnergy = 129.7 * eV;
    fNISTMaterialDataTable[174].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[174].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[174].fNumComponents        = 4;
    fNISTMaterialDataTable[174].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[174].fElementList          = ElemsZ_174;
    fNISTMaterialDataTable[174].fElementFraction      = FractionsZ_174;
    fNISTMaterialDataTable[174].fIsBuiltByAtomCount   = true;


    // NIST_MAT_N-HEPTANE --------------------------------------------------------------------
    static const int ElemsZ_175[]        = {  6,   1};
    static const double FractionsZ_175[] = {  7,  16};
    fNISTMaterialDataTable[175].fName                 = "NIST_MAT_N-HEPTANE";
    fNISTMaterialDataTable[175].fDensity              = 0.68376 * (g/cm3);
    fNISTMaterialDataTable[175].fMeanExcitationEnergy = 54.4 * eV;
    fNISTMaterialDataTable[175].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[175].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[175].fNumComponents        = 2;
    fNISTMaterialDataTable[175].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[175].fElementList          = ElemsZ_175;
    fNISTMaterialDataTable[175].fElementFraction      = FractionsZ_175;
    fNISTMaterialDataTable[175].fIsBuiltByAtomCount   = true;


    // NIST_MAT_N-HEXANE ---------------------------------------------------------------------
    static const int ElemsZ_176[]        = {  6,   1};
    static const double FractionsZ_176[] = {  6,  14};
    fNISTMaterialDataTable[176].fName                 = "NIST_MAT_N-HEXANE";
    fNISTMaterialDataTable[176].fDensity              = 0.6603 * (g/cm3);
    fNISTMaterialDataTable[176].fMeanExcitationEnergy = 54 * eV;
    fNISTMaterialDataTable[176].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[176].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[176].fNumComponents        = 2;
    fNISTMaterialDataTable[176].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[176].fElementList          = ElemsZ_176;
    fNISTMaterialDataTable[176].fElementFraction      = FractionsZ_176;
    fNISTMaterialDataTable[176].fIsBuiltByAtomCount   = true;


    // NIST_MAT_KAPTON -----------------------------------------------------------------------
    static const int ElemsZ_177[]        = {  6,   1,   7,   8};
    static const double FractionsZ_177[] = { 22,  10,   2,   5};
    fNISTMaterialDataTable[177].fName                 = "NIST_MAT_KAPTON";
    fNISTMaterialDataTable[177].fDensity              = 1.42 * (g/cm3);
    fNISTMaterialDataTable[177].fMeanExcitationEnergy = 79.6 * eV;
    fNISTMaterialDataTable[177].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[177].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[177].fNumComponents        = 4;
    fNISTMaterialDataTable[177].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[177].fElementList          = ElemsZ_177;
    fNISTMaterialDataTable[177].fElementFraction      = FractionsZ_177;
    fNISTMaterialDataTable[177].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LANTHANUM_OXYBROMIDE ---------------------------------------------------------
    static const int ElemsZ_178[]        = { 57,  35,   8};
    static const double FractionsZ_178[] = {  1,   1,   1};
    fNISTMaterialDataTable[178].fName                 = "NIST_MAT_LANTHANUM_OXYBROMIDE";
    fNISTMaterialDataTable[178].fDensity              = 6.28 * (g/cm3);
    fNISTMaterialDataTable[178].fMeanExcitationEnergy = 439.7 * eV;
    fNISTMaterialDataTable[178].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[178].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[178].fNumComponents        = 3;
    fNISTMaterialDataTable[178].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[178].fElementList          = ElemsZ_178;
    fNISTMaterialDataTable[178].fElementFraction      = FractionsZ_178;
    fNISTMaterialDataTable[178].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LANTHANUM_OXYSULFIDE ---------------------------------------------------------
    static const int ElemsZ_179[]        = { 57,   8,  16};
    static const double FractionsZ_179[] = {  2,   2,   1};
    fNISTMaterialDataTable[179].fName                 = "NIST_MAT_LANTHANUM_OXYSULFIDE";
    fNISTMaterialDataTable[179].fDensity              = 5.86 * (g/cm3);
    fNISTMaterialDataTable[179].fMeanExcitationEnergy = 421.2 * eV;
    fNISTMaterialDataTable[179].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[179].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[179].fNumComponents        = 3;
    fNISTMaterialDataTable[179].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[179].fElementList          = ElemsZ_179;
    fNISTMaterialDataTable[179].fElementFraction      = FractionsZ_179;
    fNISTMaterialDataTable[179].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LEAD_OXIDE -------------------------------------------------------------------
    static const int ElemsZ_180[]        = {  8,  82};
    static const double FractionsZ_180[] = {  0.071682,   0.928318};
    fNISTMaterialDataTable[180].fName                 = "NIST_MAT_LEAD_OXIDE";
    fNISTMaterialDataTable[180].fDensity              = 9.53 * (g/cm3);
    fNISTMaterialDataTable[180].fMeanExcitationEnergy = 766.7 * eV;
    fNISTMaterialDataTable[180].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[180].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[180].fNumComponents        = 2;
    fNISTMaterialDataTable[180].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[180].fElementList          = ElemsZ_180;
    fNISTMaterialDataTable[180].fElementFraction      = FractionsZ_180;
    fNISTMaterialDataTable[180].fIsBuiltByAtomCount   = false;


    // NIST_MAT_LITHIUM_AMIDE ----------------------------------------------------------------
    static const int ElemsZ_181[]        = {  3,   7,   1};
    static const double FractionsZ_181[] = {  1,   1,   2};
    fNISTMaterialDataTable[181].fName                 = "NIST_MAT_LITHIUM_AMIDE";
    fNISTMaterialDataTable[181].fDensity              = 1.178 * (g/cm3);
    fNISTMaterialDataTable[181].fMeanExcitationEnergy = 55.5 * eV;
    fNISTMaterialDataTable[181].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[181].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[181].fNumComponents        = 3;
    fNISTMaterialDataTable[181].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[181].fElementList          = ElemsZ_181;
    fNISTMaterialDataTable[181].fElementFraction      = FractionsZ_181;
    fNISTMaterialDataTable[181].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LITHIUM_CARBONATE ------------------------------------------------------------
    static const int ElemsZ_182[]        = {  3,   6,   8};
    static const double FractionsZ_182[] = {  2,   1,   3};
    fNISTMaterialDataTable[182].fName                 = "NIST_MAT_LITHIUM_CARBONATE";
    fNISTMaterialDataTable[182].fDensity              = 2.11 * (g/cm3);
    fNISTMaterialDataTable[182].fMeanExcitationEnergy = 87.9 * eV;
    fNISTMaterialDataTable[182].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[182].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[182].fNumComponents        = 3;
    fNISTMaterialDataTable[182].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[182].fElementList          = ElemsZ_182;
    fNISTMaterialDataTable[182].fElementFraction      = FractionsZ_182;
    fNISTMaterialDataTable[182].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LITHIUM_FLUORIDE -------------------------------------------------------------
    static const int ElemsZ_183[]        = {  3,   9};
    static const double FractionsZ_183[] = {  1,   1};
    fNISTMaterialDataTable[183].fName                 = "NIST_MAT_LITHIUM_FLUORIDE";
    fNISTMaterialDataTable[183].fDensity              = 2.635 * (g/cm3);
    fNISTMaterialDataTable[183].fMeanExcitationEnergy = 94 * eV;
    fNISTMaterialDataTable[183].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[183].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[183].fNumComponents        = 2;
    fNISTMaterialDataTable[183].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[183].fElementList          = ElemsZ_183;
    fNISTMaterialDataTable[183].fElementFraction      = FractionsZ_183;
    fNISTMaterialDataTable[183].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LITHIUM_HYDRIDE --------------------------------------------------------------
    static const int ElemsZ_184[]        = {  3,   1};
    static const double FractionsZ_184[] = {  1,   1};
    fNISTMaterialDataTable[184].fName                 = "NIST_MAT_LITHIUM_HYDRIDE";
    fNISTMaterialDataTable[184].fDensity              = 0.82 * (g/cm3);
    fNISTMaterialDataTable[184].fMeanExcitationEnergy = 36.5 * eV;
    fNISTMaterialDataTable[184].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[184].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[184].fNumComponents        = 2;
    fNISTMaterialDataTable[184].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[184].fElementList          = ElemsZ_184;
    fNISTMaterialDataTable[184].fElementFraction      = FractionsZ_184;
    fNISTMaterialDataTable[184].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LITHIUM_IODIDE ---------------------------------------------------------------
    static const int ElemsZ_185[]        = {  3,  53};
    static const double FractionsZ_185[] = {  1,   1};
    fNISTMaterialDataTable[185].fName                 = "NIST_MAT_LITHIUM_IODIDE";
    fNISTMaterialDataTable[185].fDensity              = 3.494 * (g/cm3);
    fNISTMaterialDataTable[185].fMeanExcitationEnergy = 485.1 * eV;
    fNISTMaterialDataTable[185].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[185].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[185].fNumComponents        = 2;
    fNISTMaterialDataTable[185].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[185].fElementList          = ElemsZ_185;
    fNISTMaterialDataTable[185].fElementFraction      = FractionsZ_185;
    fNISTMaterialDataTable[185].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LITHIUM_OXIDE ----------------------------------------------------------------
    static const int ElemsZ_186[]        = {  3,   8};
    static const double FractionsZ_186[] = {  2,   1};
    fNISTMaterialDataTable[186].fName                 = "NIST_MAT_LITHIUM_OXIDE";
    fNISTMaterialDataTable[186].fDensity              = 2.013 * (g/cm3);
    fNISTMaterialDataTable[186].fMeanExcitationEnergy = 73.6 * eV;
    fNISTMaterialDataTable[186].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[186].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[186].fNumComponents        = 2;
    fNISTMaterialDataTable[186].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[186].fElementList          = ElemsZ_186;
    fNISTMaterialDataTable[186].fElementFraction      = FractionsZ_186;
    fNISTMaterialDataTable[186].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LITHIUM_TETRABORATE ----------------------------------------------------------
    static const int ElemsZ_187[]        = {  3,   5,   8};
    static const double FractionsZ_187[] = {  2,   4,   7};
    fNISTMaterialDataTable[187].fName                 = "NIST_MAT_LITHIUM_TETRABORATE";
    fNISTMaterialDataTable[187].fDensity              = 2.44 * (g/cm3);
    fNISTMaterialDataTable[187].fMeanExcitationEnergy = 94.6 * eV;
    fNISTMaterialDataTable[187].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[187].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[187].fNumComponents        = 3;
    fNISTMaterialDataTable[187].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[187].fElementList          = ElemsZ_187;
    fNISTMaterialDataTable[187].fElementFraction      = FractionsZ_187;
    fNISTMaterialDataTable[187].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LUNG_ICRP --------------------------------------------------------------------
    //Adult Lung congested
    static const int ElemsZ_188[]        = {  1,   6,   7,   8,  11,  15,
                                             16,  17,  19};
    static const double FractionsZ_188[] = {     0.105,      0.083,      0.023,
                                                 0.779,      0.002,      0.001,
                                                 0.002,      0.003,      0.002};
    fNISTMaterialDataTable[188].fName                 = "NIST_MAT_LUNG_ICRP";
    fNISTMaterialDataTable[188].fDensity              = 1.04 * (g/cm3);
    fNISTMaterialDataTable[188].fMeanExcitationEnergy = 75.3 * eV;
    fNISTMaterialDataTable[188].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[188].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[188].fNumComponents        = 9;
    fNISTMaterialDataTable[188].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[188].fElementList          = ElemsZ_188;
    fNISTMaterialDataTable[188].fElementFraction      = FractionsZ_188;
    fNISTMaterialDataTable[188].fIsBuiltByAtomCount   = false;


    // NIST_MAT_M3_WAX -----------------------------------------------------------------------
    static const int ElemsZ_189[]        = {  1,   6,   8,  12,  20};
    static const double FractionsZ_189[] = {  0.114318,   0.655823,   0.092183,
                                              0.134792,   0.002883};
    fNISTMaterialDataTable[189].fName                 = "NIST_MAT_M3_WAX";
    fNISTMaterialDataTable[189].fDensity              = 1.05 * (g/cm3);
    fNISTMaterialDataTable[189].fMeanExcitationEnergy = 67.9 * eV;
    fNISTMaterialDataTable[189].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[189].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[189].fNumComponents        = 5;
    fNISTMaterialDataTable[189].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[189].fElementList          = ElemsZ_189;
    fNISTMaterialDataTable[189].fElementFraction      = FractionsZ_189;
    fNISTMaterialDataTable[189].fIsBuiltByAtomCount   = false;


    // NIST_MAT_MAGNESIUM_CARBONATE ----------------------------------------------------------
    static const int ElemsZ_190[]        = { 12,   6,   8};
    static const double FractionsZ_190[] = {  1,   1,   3};
    fNISTMaterialDataTable[190].fName                 = "NIST_MAT_MAGNESIUM_CARBONATE";
    fNISTMaterialDataTable[190].fDensity              = 2.958 * (g/cm3);
    fNISTMaterialDataTable[190].fMeanExcitationEnergy = 118 * eV;
    fNISTMaterialDataTable[190].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[190].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[190].fNumComponents        = 3;
    fNISTMaterialDataTable[190].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[190].fElementList          = ElemsZ_190;
    fNISTMaterialDataTable[190].fElementFraction      = FractionsZ_190;
    fNISTMaterialDataTable[190].fIsBuiltByAtomCount   = true;


    // NIST_MAT_MAGNESIUM_FLUORIDE -----------------------------------------------------------
    static const int ElemsZ_191[]        = { 12,   9};
    static const double FractionsZ_191[] = {  1,   2};
    fNISTMaterialDataTable[191].fName                 = "NIST_MAT_MAGNESIUM_FLUORIDE";
    fNISTMaterialDataTable[191].fDensity              = 3 * (g/cm3);
    fNISTMaterialDataTable[191].fMeanExcitationEnergy = 134.3 * eV;
    fNISTMaterialDataTable[191].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[191].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[191].fNumComponents        = 2;
    fNISTMaterialDataTable[191].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[191].fElementList          = ElemsZ_191;
    fNISTMaterialDataTable[191].fElementFraction      = FractionsZ_191;
    fNISTMaterialDataTable[191].fIsBuiltByAtomCount   = true;


    // NIST_MAT_MAGNESIUM_OXIDE --------------------------------------------------------------
    static const int ElemsZ_192[]        = { 12,   8};
    static const double FractionsZ_192[] = {  1,   1};
    fNISTMaterialDataTable[192].fName                 = "NIST_MAT_MAGNESIUM_OXIDE";
    fNISTMaterialDataTable[192].fDensity              = 3.58 * (g/cm3);
    fNISTMaterialDataTable[192].fMeanExcitationEnergy = 143.8 * eV;
    fNISTMaterialDataTable[192].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[192].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[192].fNumComponents        = 2;
    fNISTMaterialDataTable[192].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[192].fElementList          = ElemsZ_192;
    fNISTMaterialDataTable[192].fElementFraction      = FractionsZ_192;
    fNISTMaterialDataTable[192].fIsBuiltByAtomCount   = true;


    // NIST_MAT_MAGNESIUM_TETRABORATE --------------------------------------------------------
    static const int ElemsZ_193[]        = { 12,   5,   8};
    static const double FractionsZ_193[] = {  1,   4,   7};
    fNISTMaterialDataTable[193].fName                 = "NIST_MAT_MAGNESIUM_TETRABORATE";
    fNISTMaterialDataTable[193].fDensity              = 2.53 * (g/cm3);
    fNISTMaterialDataTable[193].fMeanExcitationEnergy = 108.3 * eV;
    fNISTMaterialDataTable[193].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[193].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[193].fNumComponents        = 3;
    fNISTMaterialDataTable[193].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[193].fElementList          = ElemsZ_193;
    fNISTMaterialDataTable[193].fElementFraction      = FractionsZ_193;
    fNISTMaterialDataTable[193].fIsBuiltByAtomCount   = true;


    // NIST_MAT_MERCURIC_IODIDE --------------------------------------------------------------
    static const int ElemsZ_194[]        = { 80,  53};
    static const double FractionsZ_194[] = {  1,   2};
    fNISTMaterialDataTable[194].fName                 = "NIST_MAT_MERCURIC_IODIDE";
    fNISTMaterialDataTable[194].fDensity              = 6.36 * (g/cm3);
    fNISTMaterialDataTable[194].fMeanExcitationEnergy = 684.5 * eV;
    fNISTMaterialDataTable[194].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[194].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[194].fNumComponents        = 2;
    fNISTMaterialDataTable[194].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[194].fElementList          = ElemsZ_194;
    fNISTMaterialDataTable[194].fElementFraction      = FractionsZ_194;
    fNISTMaterialDataTable[194].fIsBuiltByAtomCount   = true;


    // NIST_MAT_METHANE ----------------------------------------------------------------------
    static const int ElemsZ_195[]        = {  6,   1};
    static const double FractionsZ_195[] = {  1,   4};
    fNISTMaterialDataTable[195].fName                 = "NIST_MAT_METHANE";
    fNISTMaterialDataTable[195].fDensity              = 0.000667151 * (g/cm3);
    fNISTMaterialDataTable[195].fMeanExcitationEnergy = 41.7 * eV;
    fNISTMaterialDataTable[195].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[195].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[195].fNumComponents        = 2;
    fNISTMaterialDataTable[195].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[195].fElementList          = ElemsZ_195;
    fNISTMaterialDataTable[195].fElementFraction      = FractionsZ_195;
    fNISTMaterialDataTable[195].fIsBuiltByAtomCount   = true;


    // NIST_MAT_METHANOL ---------------------------------------------------------------------
    static const int ElemsZ_196[]        = {  6,   1,   8};
    static const double FractionsZ_196[] = {  1,   4,   1};
    fNISTMaterialDataTable[196].fName                 = "NIST_MAT_METHANOL";
    fNISTMaterialDataTable[196].fDensity              = 0.7914 * (g/cm3);
    fNISTMaterialDataTable[196].fMeanExcitationEnergy = 67.6 * eV;
    fNISTMaterialDataTable[196].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[196].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[196].fNumComponents        = 3;
    fNISTMaterialDataTable[196].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[196].fElementList          = ElemsZ_196;
    fNISTMaterialDataTable[196].fElementFraction      = FractionsZ_196;
    fNISTMaterialDataTable[196].fIsBuiltByAtomCount   = true;


    // NIST_MAT_MIX_D_WAX --------------------------------------------------------------------
    static const int ElemsZ_197[]        = {  1,   6,   8,  12,  22};
    static const double FractionsZ_197[] = {   0.13404,    0.77796,    0.03502,
                                              0.038594,   0.014386};
    fNISTMaterialDataTable[197].fName                 = "NIST_MAT_MIX_D_WAX";
    fNISTMaterialDataTable[197].fDensity              = 0.99 * (g/cm3);
    fNISTMaterialDataTable[197].fMeanExcitationEnergy = 60.9 * eV;
    fNISTMaterialDataTable[197].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[197].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[197].fNumComponents        = 5;
    fNISTMaterialDataTable[197].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[197].fElementList          = ElemsZ_197;
    fNISTMaterialDataTable[197].fElementFraction      = FractionsZ_197;
    fNISTMaterialDataTable[197].fIsBuiltByAtomCount   = false;


    // NIST_MAT_MS20_TISSUE ------------------------------------------------------------------
    static const int ElemsZ_198[]        = {  1,   6,   7,   8,  12,  17};
    static const double FractionsZ_198[] = {  0.081192,   0.583442,   0.017798,
                                              0.186381,   0.130287,     0.0009};
    fNISTMaterialDataTable[198].fName                 = "NIST_MAT_MS20_TISSUE";
    fNISTMaterialDataTable[198].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[198].fMeanExcitationEnergy = 75.1 * eV;
    fNISTMaterialDataTable[198].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[198].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[198].fNumComponents        = 6;
    fNISTMaterialDataTable[198].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[198].fElementList          = ElemsZ_198;
    fNISTMaterialDataTable[198].fElementFraction      = FractionsZ_198;
    fNISTMaterialDataTable[198].fIsBuiltByAtomCount   = false;


    // NIST_MAT_MUSCLE_SKELETAL_ICRP ---------------------------------------------------------
    static const int ElemsZ_199[]        = {  1,   6,   7,   8,  11,  15,
                                             16,  17,  19};
    static const double FractionsZ_199[] = {     0.102,      0.143,      0.034,
                                                  0.71,      0.001,      0.002,
                                                 0.003,      0.001,      0.004};
    fNISTMaterialDataTable[199].fName                 = "NIST_MAT_MUSCLE_SKELETAL_ICRP";
    fNISTMaterialDataTable[199].fDensity              = 1.05 * (g/cm3);
    fNISTMaterialDataTable[199].fMeanExcitationEnergy = 75.3 * eV;
    fNISTMaterialDataTable[199].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[199].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[199].fNumComponents        = 9;
    fNISTMaterialDataTable[199].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[199].fElementList          = ElemsZ_199;
    fNISTMaterialDataTable[199].fElementFraction      = FractionsZ_199;
    fNISTMaterialDataTable[199].fIsBuiltByAtomCount   = false;


    // NIST_MAT_MUSCLE_STRIATED_ICRU ---------------------------------------------------------
    // from old ICRU report
    static const int ElemsZ_200[]        = {  1,   6,   7,   8,  11,  15,
                                             16,  19};
    static const double FractionsZ_200[] = {     0.102,      0.123,      0.035,
                                                 0.729,      0.001,      0.002,
                                                 0.004,      0.003};
    fNISTMaterialDataTable[200].fName                 = "NIST_MAT_MUSCLE_STRIATED_ICRU";
    fNISTMaterialDataTable[200].fDensity              = 1.04 * (g/cm3);
    fNISTMaterialDataTable[200].fMeanExcitationEnergy = 74.7 * eV;
    fNISTMaterialDataTable[200].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[200].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[200].fNumComponents        = 8;
    fNISTMaterialDataTable[200].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[200].fElementList          = ElemsZ_200;
    fNISTMaterialDataTable[200].fElementFraction      = FractionsZ_200;
    fNISTMaterialDataTable[200].fIsBuiltByAtomCount   = false;


    // NIST_MAT_MUSCLE_WITH_SUCROSE ----------------------------------------------------------
    static const int ElemsZ_201[]        = {  1,   6,   7,   8};
    static const double FractionsZ_201[] = {  0.098234,   0.156214,   0.035451,
                                                0.7101};
    fNISTMaterialDataTable[201].fName                 = "NIST_MAT_MUSCLE_WITH_SUCROSE";
    fNISTMaterialDataTable[201].fDensity              = 1.11 * (g/cm3);
    fNISTMaterialDataTable[201].fMeanExcitationEnergy = 74.3 * eV;
    fNISTMaterialDataTable[201].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[201].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[201].fNumComponents        = 4;
    fNISTMaterialDataTable[201].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[201].fElementList          = ElemsZ_201;
    fNISTMaterialDataTable[201].fElementFraction      = FractionsZ_201;
    fNISTMaterialDataTable[201].fIsBuiltByAtomCount   = false;


    // NIST_MAT_MUSCLE_WITHOUT_SUCROSE -------------------------------------------------------
    static const int ElemsZ_202[]        = {  1,   6,   7,   8};
    static const double FractionsZ_202[] = {  0.101969,   0.120058,   0.035451,
                                              0.742522};
    fNISTMaterialDataTable[202].fName                 = "NIST_MAT_MUSCLE_WITHOUT_SUCROSE";
    fNISTMaterialDataTable[202].fDensity              = 1.07 * (g/cm3);
    fNISTMaterialDataTable[202].fMeanExcitationEnergy = 74.2 * eV;
    fNISTMaterialDataTable[202].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[202].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[202].fNumComponents        = 4;
    fNISTMaterialDataTable[202].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[202].fElementList          = ElemsZ_202;
    fNISTMaterialDataTable[202].fElementFraction      = FractionsZ_202;
    fNISTMaterialDataTable[202].fIsBuiltByAtomCount   = false;


    // NIST_MAT_NAPHTHALENE ------------------------------------------------------------------
    static const int ElemsZ_203[]        = {  6,   1};
    static const double FractionsZ_203[] = { 10,   8};
    fNISTMaterialDataTable[203].fName                 = "NIST_MAT_NAPHTHALENE";
    fNISTMaterialDataTable[203].fDensity              = 1.145 * (g/cm3);
    fNISTMaterialDataTable[203].fMeanExcitationEnergy = 68.4 * eV;
    fNISTMaterialDataTable[203].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[203].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[203].fNumComponents        = 2;
    fNISTMaterialDataTable[203].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[203].fElementList          = ElemsZ_203;
    fNISTMaterialDataTable[203].fElementFraction      = FractionsZ_203;
    fNISTMaterialDataTable[203].fIsBuiltByAtomCount   = true;


    // NIST_MAT_NITROBENZENE -----------------------------------------------------------------
    static const int ElemsZ_204[]        = {  6,   1,   7,   8};
    static const double FractionsZ_204[] = {  6,   5,   1,   2};
    fNISTMaterialDataTable[204].fName                 = "NIST_MAT_NITROBENZENE";
    fNISTMaterialDataTable[204].fDensity              = 1.19867 * (g/cm3);
    fNISTMaterialDataTable[204].fMeanExcitationEnergy = 75.8 * eV;
    fNISTMaterialDataTable[204].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[204].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[204].fNumComponents        = 4;
    fNISTMaterialDataTable[204].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[204].fElementList          = ElemsZ_204;
    fNISTMaterialDataTable[204].fElementFraction      = FractionsZ_204;
    fNISTMaterialDataTable[204].fIsBuiltByAtomCount   = true;


    // NIST_MAT_NITROUS_OXIDE ----------------------------------------------------------------
    static const int ElemsZ_205[]        = {  7,   8};
    static const double FractionsZ_205[] = {  2,   1};
    fNISTMaterialDataTable[205].fName                 = "NIST_MAT_NITROUS_OXIDE";
    fNISTMaterialDataTable[205].fDensity              = 0.00183094 * (g/cm3);
    fNISTMaterialDataTable[205].fMeanExcitationEnergy = 84.9 * eV;
    fNISTMaterialDataTable[205].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[205].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[205].fNumComponents        = 2;
    fNISTMaterialDataTable[205].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[205].fElementList          = ElemsZ_205;
    fNISTMaterialDataTable[205].fElementFraction      = FractionsZ_205;
    fNISTMaterialDataTable[205].fIsBuiltByAtomCount   = true;


    // NIST_MAT_NYLON-8062 -------------------------------------------------------------------
    static const int ElemsZ_206[]        = {  1,   6,   7,   8};
    static const double FractionsZ_206[] = {  0.103509,   0.648415,   0.099536,
                                              0.148539};
    fNISTMaterialDataTable[206].fName                 = "NIST_MAT_NYLON-8062";
    fNISTMaterialDataTable[206].fDensity              = 1.08 * (g/cm3);
    fNISTMaterialDataTable[206].fMeanExcitationEnergy = 64.3 * eV;
    fNISTMaterialDataTable[206].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[206].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[206].fNumComponents        = 4;
    fNISTMaterialDataTable[206].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[206].fElementList          = ElemsZ_206;
    fNISTMaterialDataTable[206].fElementFraction      = FractionsZ_206;
    fNISTMaterialDataTable[206].fIsBuiltByAtomCount   = false;


    // NIST_MAT_NYLON-6-6 --------------------------------------------------------------------
    static const int ElemsZ_207[]        = {  6,   1,   7,   8};
    static const double FractionsZ_207[] = {  6,  11,   1,   1};
    fNISTMaterialDataTable[207].fName                 = "NIST_MAT_NYLON-6-6";
    fNISTMaterialDataTable[207].fDensity              = 1.14 * (g/cm3);
    fNISTMaterialDataTable[207].fMeanExcitationEnergy = 63.9 * eV;
    fNISTMaterialDataTable[207].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[207].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[207].fNumComponents        = 4;
    fNISTMaterialDataTable[207].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[207].fElementList          = ElemsZ_207;
    fNISTMaterialDataTable[207].fElementFraction      = FractionsZ_207;
    fNISTMaterialDataTable[207].fIsBuiltByAtomCount   = true;


    // NIST_MAT_NYLON-6-10 -------------------------------------------------------------------
    static const int ElemsZ_208[]        = {  1,   6,   7,   8};
    static const double FractionsZ_208[] = {  0.107062,   0.680449,   0.099189,
                                                0.1133};
    fNISTMaterialDataTable[208].fName                 = "NIST_MAT_NYLON-6-10";
    fNISTMaterialDataTable[208].fDensity              = 1.14 * (g/cm3);
    fNISTMaterialDataTable[208].fMeanExcitationEnergy = 63.2 * eV;
    fNISTMaterialDataTable[208].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[208].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[208].fNumComponents        = 4;
    fNISTMaterialDataTable[208].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[208].fElementList          = ElemsZ_208;
    fNISTMaterialDataTable[208].fElementFraction      = FractionsZ_208;
    fNISTMaterialDataTable[208].fIsBuiltByAtomCount   = false;


    // NIST_MAT_NYLON-11_RILSAN --------------------------------------------------------------
    static const int ElemsZ_209[]        = {  1,   6,   7,   8};
    static const double FractionsZ_209[] = {  0.115476,   0.720819,   0.076417,
                                              0.087289};
    fNISTMaterialDataTable[209].fName                 = "NIST_MAT_NYLON-11_RILSAN";
    fNISTMaterialDataTable[209].fDensity              = 1.425 * (g/cm3);
    fNISTMaterialDataTable[209].fMeanExcitationEnergy = 61.6 * eV;
    fNISTMaterialDataTable[209].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[209].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[209].fNumComponents        = 4;
    fNISTMaterialDataTable[209].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[209].fElementList          = ElemsZ_209;
    fNISTMaterialDataTable[209].fElementFraction      = FractionsZ_209;
    fNISTMaterialDataTable[209].fIsBuiltByAtomCount   = false;


    // NIST_MAT_OCTANE -----------------------------------------------------------------------
    static const int ElemsZ_210[]        = {  6,   1};
    static const double FractionsZ_210[] = {  8,  18};
    fNISTMaterialDataTable[210].fName                 = "NIST_MAT_OCTANE";
    fNISTMaterialDataTable[210].fDensity              = 0.7026 * (g/cm3);
    fNISTMaterialDataTable[210].fMeanExcitationEnergy = 54.7 * eV;
    fNISTMaterialDataTable[210].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[210].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[210].fNumComponents        = 2;
    fNISTMaterialDataTable[210].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[210].fElementList          = ElemsZ_210;
    fNISTMaterialDataTable[210].fElementFraction      = FractionsZ_210;
    fNISTMaterialDataTable[210].fIsBuiltByAtomCount   = true;


    // NIST_MAT_PARAFFIN ---------------------------------------------------------------------
    static const int ElemsZ_211[]        = {  6,   1};
    static const double FractionsZ_211[] = { 25,  52};
    fNISTMaterialDataTable[211].fName                 = "NIST_MAT_PARAFFIN";
    fNISTMaterialDataTable[211].fDensity              = 0.93 * (g/cm3);
    fNISTMaterialDataTable[211].fMeanExcitationEnergy = 55.9 * eV;
    fNISTMaterialDataTable[211].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[211].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[211].fNumComponents        = 2;
    fNISTMaterialDataTable[211].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[211].fElementList          = ElemsZ_211;
    fNISTMaterialDataTable[211].fElementFraction      = FractionsZ_211;
    fNISTMaterialDataTable[211].fIsBuiltByAtomCount   = true;


    // NIST_MAT_N-PENTANE --------------------------------------------------------------------
    static const int ElemsZ_212[]        = {  6,   1};
    static const double FractionsZ_212[] = {  5,  12};
    fNISTMaterialDataTable[212].fName                 = "NIST_MAT_N-PENTANE";
    fNISTMaterialDataTable[212].fDensity              = 0.6262 * (g/cm3);
    fNISTMaterialDataTable[212].fMeanExcitationEnergy = 53.6 * eV;
    fNISTMaterialDataTable[212].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[212].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[212].fNumComponents        = 2;
    fNISTMaterialDataTable[212].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[212].fElementList          = ElemsZ_212;
    fNISTMaterialDataTable[212].fElementFraction      = FractionsZ_212;
    fNISTMaterialDataTable[212].fIsBuiltByAtomCount   = true;


    // NIST_MAT_PHOTO_EMULSION ---------------------------------------------------------------
    static const int ElemsZ_213[]        = {  1,   6,   7,   8,  16,  35,
                                             47,  53};
    static const double FractionsZ_213[] = {    0.0141,   0.072261,    0.01932,
                                              0.066101,    0.00189,   0.349103,
                                              0.474105,    0.00312};
    fNISTMaterialDataTable[213].fName                 = "NIST_MAT_PHOTO_EMULSION";
    fNISTMaterialDataTable[213].fDensity              = 3.815 * (g/cm3);
    fNISTMaterialDataTable[213].fMeanExcitationEnergy = 331 * eV;
    fNISTMaterialDataTable[213].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[213].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[213].fNumComponents        = 8;
    fNISTMaterialDataTable[213].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[213].fElementList          = ElemsZ_213;
    fNISTMaterialDataTable[213].fElementFraction      = FractionsZ_213;
    fNISTMaterialDataTable[213].fIsBuiltByAtomCount   = false;


    // NIST_MAT_PLASTIC_SC_VINYLTOLUENE ------------------------------------------------------
    static const int ElemsZ_214[]        = {  6,   1};
    static const double FractionsZ_214[] = {  9,  10};
    fNISTMaterialDataTable[214].fName                 = "NIST_MAT_PLASTIC_SC_VINYLTOLUENE";
    fNISTMaterialDataTable[214].fDensity              = 1.032 * (g/cm3);
    fNISTMaterialDataTable[214].fMeanExcitationEnergy = 64.7 * eV;
    fNISTMaterialDataTable[214].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[214].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[214].fNumComponents        = 2;
    fNISTMaterialDataTable[214].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[214].fElementList          = ElemsZ_214;
    fNISTMaterialDataTable[214].fElementFraction      = FractionsZ_214;
    fNISTMaterialDataTable[214].fIsBuiltByAtomCount   = true;


    // NIST_MAT_PLUTONIUM_DIOXIDE ------------------------------------------------------------
    static const int ElemsZ_215[]        = { 94,   8};
    static const double FractionsZ_215[] = {  1,   2};
    fNISTMaterialDataTable[215].fName                 = "NIST_MAT_PLUTONIUM_DIOXIDE";
    fNISTMaterialDataTable[215].fDensity              = 11.46 * (g/cm3);
    fNISTMaterialDataTable[215].fMeanExcitationEnergy = 746.5 * eV;
    fNISTMaterialDataTable[215].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[215].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[215].fNumComponents        = 2;
    fNISTMaterialDataTable[215].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[215].fElementList          = ElemsZ_215;
    fNISTMaterialDataTable[215].fElementFraction      = FractionsZ_215;
    fNISTMaterialDataTable[215].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYACRYLONITRILE ------------------------------------------------------------
    static const int ElemsZ_216[]        = {  6,   1,   7};
    static const double FractionsZ_216[] = {  3,   3,   1};
    fNISTMaterialDataTable[216].fName                 = "NIST_MAT_POLYACRYLONITRILE";
    fNISTMaterialDataTable[216].fDensity              = 1.17 * (g/cm3);
    fNISTMaterialDataTable[216].fMeanExcitationEnergy = 69.6 * eV;
    fNISTMaterialDataTable[216].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[216].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[216].fNumComponents        = 3;
    fNISTMaterialDataTable[216].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[216].fElementList          = ElemsZ_216;
    fNISTMaterialDataTable[216].fElementFraction      = FractionsZ_216;
    fNISTMaterialDataTable[216].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYCARBONATE ----------------------------------------------------------------
    static const int ElemsZ_217[]        = {  6,   1,   8};
    static const double FractionsZ_217[] = { 16,  14,   3};
    fNISTMaterialDataTable[217].fName                 = "NIST_MAT_POLYCARBONATE";
    fNISTMaterialDataTable[217].fDensity              = 1.2 * (g/cm3);
    fNISTMaterialDataTable[217].fMeanExcitationEnergy = 73.1 * eV;
    fNISTMaterialDataTable[217].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[217].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[217].fNumComponents        = 3;
    fNISTMaterialDataTable[217].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[217].fElementList          = ElemsZ_217;
    fNISTMaterialDataTable[217].fElementFraction      = FractionsZ_217;
    fNISTMaterialDataTable[217].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYCHLOROSTYRENE ------------------------------------------------------------
    static const int ElemsZ_218[]        = {  6,   1,  17};
    static const double FractionsZ_218[] = {  8,   7,   1};
    fNISTMaterialDataTable[218].fName                 = "NIST_MAT_POLYCHLOROSTYRENE";
    fNISTMaterialDataTable[218].fDensity              = 1.3 * (g/cm3);
    fNISTMaterialDataTable[218].fMeanExcitationEnergy = 81.7 * eV;
    fNISTMaterialDataTable[218].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[218].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[218].fNumComponents        = 3;
    fNISTMaterialDataTable[218].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[218].fElementList          = ElemsZ_218;
    fNISTMaterialDataTable[218].fElementFraction      = FractionsZ_218;
    fNISTMaterialDataTable[218].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYETHYLENE -----------------------------------------------------------------
    static const int ElemsZ_219[]        = {  6,   1};
    static const double FractionsZ_219[] = {  1,   2};
    fNISTMaterialDataTable[219].fName                 = "NIST_MAT_POLYETHYLENE";
    fNISTMaterialDataTable[219].fDensity              = 0.94 * (g/cm3);
    fNISTMaterialDataTable[219].fMeanExcitationEnergy = 57.4 * eV;
    fNISTMaterialDataTable[219].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[219].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[219].fNumComponents        = 2;
    fNISTMaterialDataTable[219].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[219].fElementList          = ElemsZ_219;
    fNISTMaterialDataTable[219].fElementFraction      = FractionsZ_219;
    fNISTMaterialDataTable[219].fIsBuiltByAtomCount   = true;


    // NIST_MAT_MYLAR ------------------------------------------------------------------------
    static const int ElemsZ_220[]        = {  6,   1,   8};
    static const double FractionsZ_220[] = { 10,   8,   4};
    fNISTMaterialDataTable[220].fName                 = "NIST_MAT_MYLAR";
    fNISTMaterialDataTable[220].fDensity              = 1.4 * (g/cm3);
    fNISTMaterialDataTable[220].fMeanExcitationEnergy = 78.7 * eV;
    fNISTMaterialDataTable[220].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[220].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[220].fNumComponents        = 3;
    fNISTMaterialDataTable[220].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[220].fElementList          = ElemsZ_220;
    fNISTMaterialDataTable[220].fElementFraction      = FractionsZ_220;
    fNISTMaterialDataTable[220].fIsBuiltByAtomCount   = true;


    // NIST_MAT_PLEXIGLASS -------------------------------------------------------------------
    static const int ElemsZ_221[]        = {  6,   1,   8};
    static const double FractionsZ_221[] = {  5,   8,   2};
    fNISTMaterialDataTable[221].fName                 = "NIST_MAT_PLEXIGLASS";
    fNISTMaterialDataTable[221].fDensity              = 1.19 * (g/cm3);
    fNISTMaterialDataTable[221].fMeanExcitationEnergy = 74 * eV;
    fNISTMaterialDataTable[221].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[221].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[221].fNumComponents        = 3;
    fNISTMaterialDataTable[221].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[221].fElementList          = ElemsZ_221;
    fNISTMaterialDataTable[221].fElementFraction      = FractionsZ_221;
    fNISTMaterialDataTable[221].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYOXYMETHYLENE -------------------------------------------------------------
    static const int ElemsZ_222[]        = {  6,   1,   8};
    static const double FractionsZ_222[] = {  1,   2,   1};
    fNISTMaterialDataTable[222].fName                 = "NIST_MAT_POLYOXYMETHYLENE";
    fNISTMaterialDataTable[222].fDensity              = 1.425 * (g/cm3);
    fNISTMaterialDataTable[222].fMeanExcitationEnergy = 77.4 * eV;
    fNISTMaterialDataTable[222].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[222].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[222].fNumComponents        = 3;
    fNISTMaterialDataTable[222].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[222].fElementList          = ElemsZ_222;
    fNISTMaterialDataTable[222].fElementFraction      = FractionsZ_222;
    fNISTMaterialDataTable[222].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYPROPYLENE ----------------------------------------------------------------
    static const int ElemsZ_223[]        = {  6,   1};
    static const double FractionsZ_223[] = {  2,   4};
    fNISTMaterialDataTable[223].fName                 = "NIST_MAT_POLYPROPYLENE";
    fNISTMaterialDataTable[223].fDensity              = 0.9 * (g/cm3);
    fNISTMaterialDataTable[223].fMeanExcitationEnergy = 56.5 * eV;
    fNISTMaterialDataTable[223].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[223].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[223].fNumComponents        = 2;
    fNISTMaterialDataTable[223].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[223].fElementList          = ElemsZ_223;
    fNISTMaterialDataTable[223].fElementFraction      = FractionsZ_223;
    fNISTMaterialDataTable[223].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYSTYRENE ------------------------------------------------------------------
    static const int ElemsZ_224[]        = {  6,   1};
    static const double FractionsZ_224[] = {  8,   8};
    fNISTMaterialDataTable[224].fName                 = "NIST_MAT_POLYSTYRENE";
    fNISTMaterialDataTable[224].fDensity              = 1.06 * (g/cm3);
    fNISTMaterialDataTable[224].fMeanExcitationEnergy = 68.7 * eV;
    fNISTMaterialDataTable[224].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[224].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[224].fNumComponents        = 2;
    fNISTMaterialDataTable[224].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[224].fElementList          = ElemsZ_224;
    fNISTMaterialDataTable[224].fElementFraction      = FractionsZ_224;
    fNISTMaterialDataTable[224].fIsBuiltByAtomCount   = true;


    // NIST_MAT_TEFLON -----------------------------------------------------------------------
    static const int ElemsZ_225[]        = {  6,   9};
    static const double FractionsZ_225[] = {  2,   4};
    fNISTMaterialDataTable[225].fName                 = "NIST_MAT_TEFLON";
    fNISTMaterialDataTable[225].fDensity              = 2.2 * (g/cm3);
    fNISTMaterialDataTable[225].fMeanExcitationEnergy = 99.1 * eV;
    fNISTMaterialDataTable[225].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[225].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[225].fNumComponents        = 2;
    fNISTMaterialDataTable[225].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[225].fElementList          = ElemsZ_225;
    fNISTMaterialDataTable[225].fElementFraction      = FractionsZ_225;
    fNISTMaterialDataTable[225].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYTRIFLUOROCHLOROETHYLENE --------------------------------------------------
    static const int ElemsZ_226[]        = {  6,   9,  17};
    static const double FractionsZ_226[] = {  2,   3,   1};
    fNISTMaterialDataTable[226].fName                 = "NIST_MAT_POLYTRIFLUOROCHLOROETHYLENE";
    fNISTMaterialDataTable[226].fDensity              = 2.1 * (g/cm3);
    fNISTMaterialDataTable[226].fMeanExcitationEnergy = 120.7 * eV;
    fNISTMaterialDataTable[226].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[226].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[226].fNumComponents        = 3;
    fNISTMaterialDataTable[226].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[226].fElementList          = ElemsZ_226;
    fNISTMaterialDataTable[226].fElementFraction      = FractionsZ_226;
    fNISTMaterialDataTable[226].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYVINYL_ACETATE ------------------------------------------------------------
    static const int ElemsZ_227[]        = {  6,   1,   8};
    static const double FractionsZ_227[] = {  4,   6,   2};
    fNISTMaterialDataTable[227].fName                 = "NIST_MAT_POLYVINYL_ACETATE";
    fNISTMaterialDataTable[227].fDensity              = 1.19 * (g/cm3);
    fNISTMaterialDataTable[227].fMeanExcitationEnergy = 73.7 * eV;
    fNISTMaterialDataTable[227].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[227].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[227].fNumComponents        = 3;
    fNISTMaterialDataTable[227].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[227].fElementList          = ElemsZ_227;
    fNISTMaterialDataTable[227].fElementFraction      = FractionsZ_227;
    fNISTMaterialDataTable[227].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYVINYL_ALCOHOL ------------------------------------------------------------
    static const int ElemsZ_228[]        = {  6,   1,   8};
    static const double FractionsZ_228[] = {  2,   4,   1};
    fNISTMaterialDataTable[228].fName                 = "NIST_MAT_POLYVINYL_ALCOHOL";
    fNISTMaterialDataTable[228].fDensity              = 1.3 * (g/cm3);
    fNISTMaterialDataTable[228].fMeanExcitationEnergy = 69.7 * eV;
    fNISTMaterialDataTable[228].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[228].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[228].fNumComponents        = 3;
    fNISTMaterialDataTable[228].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[228].fElementList          = ElemsZ_228;
    fNISTMaterialDataTable[228].fElementFraction      = FractionsZ_228;
    fNISTMaterialDataTable[228].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYVINYL_BUTYRAL ------------------------------------------------------------
    static const int ElemsZ_229[]        = {  6,   1,   8};
    static const double FractionsZ_229[] = {  8,  14,   2};
    fNISTMaterialDataTable[229].fName                 = "NIST_MAT_POLYVINYL_BUTYRAL";
    fNISTMaterialDataTable[229].fDensity              = 1.12 * (g/cm3);
    fNISTMaterialDataTable[229].fMeanExcitationEnergy = 67.2 * eV;
    fNISTMaterialDataTable[229].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[229].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[229].fNumComponents        = 3;
    fNISTMaterialDataTable[229].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[229].fElementList          = ElemsZ_229;
    fNISTMaterialDataTable[229].fElementFraction      = FractionsZ_229;
    fNISTMaterialDataTable[229].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYVINYL_CHLORIDE -----------------------------------------------------------
    static const int ElemsZ_230[]        = {  6,   1,  17};
    static const double FractionsZ_230[] = {  2,   3,   1};
    fNISTMaterialDataTable[230].fName                 = "NIST_MAT_POLYVINYL_CHLORIDE";
    fNISTMaterialDataTable[230].fDensity              = 1.3 * (g/cm3);
    fNISTMaterialDataTable[230].fMeanExcitationEnergy = 108.2 * eV;
    fNISTMaterialDataTable[230].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[230].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[230].fNumComponents        = 3;
    fNISTMaterialDataTable[230].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[230].fElementList          = ElemsZ_230;
    fNISTMaterialDataTable[230].fElementFraction      = FractionsZ_230;
    fNISTMaterialDataTable[230].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYVINYLIDENE_CHLORIDE ------------------------------------------------------
    static const int ElemsZ_231[]        = {  6,   1,  17};
    static const double FractionsZ_231[] = {  2,   2,   2};
    fNISTMaterialDataTable[231].fName                 = "NIST_MAT_POLYVINYLIDENE_CHLORIDE";
    fNISTMaterialDataTable[231].fDensity              = 1.7 * (g/cm3);
    fNISTMaterialDataTable[231].fMeanExcitationEnergy = 134.3 * eV;
    fNISTMaterialDataTable[231].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[231].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[231].fNumComponents        = 3;
    fNISTMaterialDataTable[231].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[231].fElementList          = ElemsZ_231;
    fNISTMaterialDataTable[231].fElementFraction      = FractionsZ_231;
    fNISTMaterialDataTable[231].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYVINYLIDENE_FLUORIDE ------------------------------------------------------
    static const int ElemsZ_232[]        = {  6,   1,   9};
    static const double FractionsZ_232[] = {  2,   2,   2};
    fNISTMaterialDataTable[232].fName                 = "NIST_MAT_POLYVINYLIDENE_FLUORIDE";
    fNISTMaterialDataTable[232].fDensity              = 1.76 * (g/cm3);
    fNISTMaterialDataTable[232].fMeanExcitationEnergy = 88.8 * eV;
    fNISTMaterialDataTable[232].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[232].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[232].fNumComponents        = 3;
    fNISTMaterialDataTable[232].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[232].fElementList          = ElemsZ_232;
    fNISTMaterialDataTable[232].fElementFraction      = FractionsZ_232;
    fNISTMaterialDataTable[232].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POLYVINYL_PYRROLIDONE --------------------------------------------------------
    static const int ElemsZ_233[]        = {  6,   1,   7,   8};
    static const double FractionsZ_233[] = {  6,   9,   1,   1};
    fNISTMaterialDataTable[233].fName                 = "NIST_MAT_POLYVINYL_PYRROLIDONE";
    fNISTMaterialDataTable[233].fDensity              = 1.25 * (g/cm3);
    fNISTMaterialDataTable[233].fMeanExcitationEnergy = 67.7 * eV;
    fNISTMaterialDataTable[233].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[233].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[233].fNumComponents        = 4;
    fNISTMaterialDataTable[233].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[233].fElementList          = ElemsZ_233;
    fNISTMaterialDataTable[233].fElementFraction      = FractionsZ_233;
    fNISTMaterialDataTable[233].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POTASSIUM_IODIDE -------------------------------------------------------------
    static const int ElemsZ_234[]        = { 19,  53};
    static const double FractionsZ_234[] = {  1,   1};
    fNISTMaterialDataTable[234].fName                 = "NIST_MAT_POTASSIUM_IODIDE";
    fNISTMaterialDataTable[234].fDensity              = 3.13 * (g/cm3);
    fNISTMaterialDataTable[234].fMeanExcitationEnergy = 431.9 * eV;
    fNISTMaterialDataTable[234].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[234].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[234].fNumComponents        = 2;
    fNISTMaterialDataTable[234].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[234].fElementList          = ElemsZ_234;
    fNISTMaterialDataTable[234].fElementFraction      = FractionsZ_234;
    fNISTMaterialDataTable[234].fIsBuiltByAtomCount   = true;


    // NIST_MAT_POTASSIUM_OXIDE --------------------------------------------------------------
    static const int ElemsZ_235[]        = { 19,   8};
    static const double FractionsZ_235[] = {  2,   1};
    fNISTMaterialDataTable[235].fName                 = "NIST_MAT_POTASSIUM_OXIDE";
    fNISTMaterialDataTable[235].fDensity              = 2.32 * (g/cm3);
    fNISTMaterialDataTable[235].fMeanExcitationEnergy = 189.9 * eV;
    fNISTMaterialDataTable[235].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[235].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[235].fNumComponents        = 2;
    fNISTMaterialDataTable[235].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[235].fElementList          = ElemsZ_235;
    fNISTMaterialDataTable[235].fElementFraction      = FractionsZ_235;
    fNISTMaterialDataTable[235].fIsBuiltByAtomCount   = true;


    // NIST_MAT_PROPANE ----------------------------------------------------------------------
    static const int ElemsZ_236[]        = {  6,   1};
    static const double FractionsZ_236[] = {  3,   8};
    fNISTMaterialDataTable[236].fName                 = "NIST_MAT_PROPANE";
    fNISTMaterialDataTable[236].fDensity              = 0.00187939 * (g/cm3);
    fNISTMaterialDataTable[236].fMeanExcitationEnergy = 47.1 * eV;
    fNISTMaterialDataTable[236].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[236].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[236].fNumComponents        = 2;
    fNISTMaterialDataTable[236].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[236].fElementList          = ElemsZ_236;
    fNISTMaterialDataTable[236].fElementFraction      = FractionsZ_236;
    fNISTMaterialDataTable[236].fIsBuiltByAtomCount   = true;


    // NIST_MAT_lPROPANE ---------------------------------------------------------------------
    static const int ElemsZ_237[]        = {  6,   1};
    static const double FractionsZ_237[] = {  3,   8};
    fNISTMaterialDataTable[237].fName                 = "NIST_MAT_lPROPANE";
    fNISTMaterialDataTable[237].fDensity              = 0.43 * (g/cm3);
    fNISTMaterialDataTable[237].fMeanExcitationEnergy = 52 * eV;
    fNISTMaterialDataTable[237].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[237].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[237].fNumComponents        = 2;
    fNISTMaterialDataTable[237].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[237].fElementList          = ElemsZ_237;
    fNISTMaterialDataTable[237].fElementFraction      = FractionsZ_237;
    fNISTMaterialDataTable[237].fIsBuiltByAtomCount   = true;


    // NIST_MAT_N-PROPYL_ALCOHOL -------------------------------------------------------------
    static const int ElemsZ_238[]        = {  6,   1,   8};
    static const double FractionsZ_238[] = {  3,   8,   1};
    fNISTMaterialDataTable[238].fName                 = "NIST_MAT_N-PROPYL_ALCOHOL";
    fNISTMaterialDataTable[238].fDensity              = 0.8035 * (g/cm3);
    fNISTMaterialDataTable[238].fMeanExcitationEnergy = 61.1 * eV;
    fNISTMaterialDataTable[238].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[238].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[238].fNumComponents        = 3;
    fNISTMaterialDataTable[238].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[238].fElementList          = ElemsZ_238;
    fNISTMaterialDataTable[238].fElementFraction      = FractionsZ_238;
    fNISTMaterialDataTable[238].fIsBuiltByAtomCount   = true;


    // NIST_MAT_PYRIDINE ---------------------------------------------------------------------
    static const int ElemsZ_239[]        = {  6,   1,   7};
    static const double FractionsZ_239[] = {  5,   5,   1};
    fNISTMaterialDataTable[239].fName                 = "NIST_MAT_PYRIDINE";
    fNISTMaterialDataTable[239].fDensity              = 0.9819 * (g/cm3);
    fNISTMaterialDataTable[239].fMeanExcitationEnergy = 66.2 * eV;
    fNISTMaterialDataTable[239].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[239].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[239].fNumComponents        = 3;
    fNISTMaterialDataTable[239].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[239].fElementList          = ElemsZ_239;
    fNISTMaterialDataTable[239].fElementFraction      = FractionsZ_239;
    fNISTMaterialDataTable[239].fIsBuiltByAtomCount   = true;


    // NIST_MAT_RUBBER_BUTYL -----------------------------------------------------------------
    static const int ElemsZ_240[]        = {  1,   6};
    static const double FractionsZ_240[] = {  0.143711,   0.856289};
    fNISTMaterialDataTable[240].fName                 = "NIST_MAT_RUBBER_BUTYL";
    fNISTMaterialDataTable[240].fDensity              = 0.92 * (g/cm3);
    fNISTMaterialDataTable[240].fMeanExcitationEnergy = 56.5 * eV;
    fNISTMaterialDataTable[240].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[240].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[240].fNumComponents        = 2;
    fNISTMaterialDataTable[240].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[240].fElementList          = ElemsZ_240;
    fNISTMaterialDataTable[240].fElementFraction      = FractionsZ_240;
    fNISTMaterialDataTable[240].fIsBuiltByAtomCount   = false;


    // NIST_MAT_RUBBER_NATURAL ---------------------------------------------------------------
    static const int ElemsZ_241[]        = {  1,   6};
    static const double FractionsZ_241[] = {  0.118371,   0.881629};
    fNISTMaterialDataTable[241].fName                 = "NIST_MAT_RUBBER_NATURAL";
    fNISTMaterialDataTable[241].fDensity              = 0.92 * (g/cm3);
    fNISTMaterialDataTable[241].fMeanExcitationEnergy = 59.8 * eV;
    fNISTMaterialDataTable[241].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[241].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[241].fNumComponents        = 2;
    fNISTMaterialDataTable[241].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[241].fElementList          = ElemsZ_241;
    fNISTMaterialDataTable[241].fElementFraction      = FractionsZ_241;
    fNISTMaterialDataTable[241].fIsBuiltByAtomCount   = false;


    // NIST_MAT_RUBBER_NEOPRENE --------------------------------------------------------------
    static const int ElemsZ_242[]        = {  1,   6,  17};
    static const double FractionsZ_242[] = {   0.05692,   0.542646,   0.400434};
    fNISTMaterialDataTable[242].fName                 = "NIST_MAT_RUBBER_NEOPRENE";
    fNISTMaterialDataTable[242].fDensity              = 1.23 * (g/cm3);
    fNISTMaterialDataTable[242].fMeanExcitationEnergy = 93 * eV;
    fNISTMaterialDataTable[242].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[242].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[242].fNumComponents        = 3;
    fNISTMaterialDataTable[242].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[242].fElementList          = ElemsZ_242;
    fNISTMaterialDataTable[242].fElementFraction      = FractionsZ_242;
    fNISTMaterialDataTable[242].fIsBuiltByAtomCount   = false;


    // NIST_MAT_SILICON_DIOXIDE --------------------------------------------------------------
    static const int ElemsZ_243[]        = { 14,   8};
    static const double FractionsZ_243[] = {  1,   2};
    fNISTMaterialDataTable[243].fName                 = "NIST_MAT_SILICON_DIOXIDE";
    fNISTMaterialDataTable[243].fDensity              = 2.32 * (g/cm3);
    fNISTMaterialDataTable[243].fMeanExcitationEnergy = 139.2 * eV;
    fNISTMaterialDataTable[243].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[243].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[243].fNumComponents        = 2;
    fNISTMaterialDataTable[243].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[243].fElementList          = ElemsZ_243;
    fNISTMaterialDataTable[243].fElementFraction      = FractionsZ_243;
    fNISTMaterialDataTable[243].fIsBuiltByAtomCount   = true;


    // NIST_MAT_SILVER_BROMIDE ---------------------------------------------------------------
    static const int ElemsZ_244[]        = { 47,  35};
    static const double FractionsZ_244[] = {  1,   1};
    fNISTMaterialDataTable[244].fName                 = "NIST_MAT_SILVER_BROMIDE";
    fNISTMaterialDataTable[244].fDensity              = 6.473 * (g/cm3);
    fNISTMaterialDataTable[244].fMeanExcitationEnergy = 486.6 * eV;
    fNISTMaterialDataTable[244].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[244].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[244].fNumComponents        = 2;
    fNISTMaterialDataTable[244].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[244].fElementList          = ElemsZ_244;
    fNISTMaterialDataTable[244].fElementFraction      = FractionsZ_244;
    fNISTMaterialDataTable[244].fIsBuiltByAtomCount   = true;


    // NIST_MAT_SILVER_CHLORIDE --------------------------------------------------------------
    static const int ElemsZ_245[]        = { 47,  17};
    static const double FractionsZ_245[] = {  1,   1};
    fNISTMaterialDataTable[245].fName                 = "NIST_MAT_SILVER_CHLORIDE";
    fNISTMaterialDataTable[245].fDensity              = 5.56 * (g/cm3);
    fNISTMaterialDataTable[245].fMeanExcitationEnergy = 398.4 * eV;
    fNISTMaterialDataTable[245].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[245].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[245].fNumComponents        = 2;
    fNISTMaterialDataTable[245].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[245].fElementList          = ElemsZ_245;
    fNISTMaterialDataTable[245].fElementFraction      = FractionsZ_245;
    fNISTMaterialDataTable[245].fIsBuiltByAtomCount   = true;


    // NIST_MAT_SILVER_HALIDES ---------------------------------------------------------------
    static const int ElemsZ_246[]        = { 35,  47,  53};
    static const double FractionsZ_246[] = {  0.422895,   0.573748,   0.003357};
    fNISTMaterialDataTable[246].fName                 = "NIST_MAT_SILVER_HALIDES";
    fNISTMaterialDataTable[246].fDensity              = 6.47 * (g/cm3);
    fNISTMaterialDataTable[246].fMeanExcitationEnergy = 487.1 * eV;
    fNISTMaterialDataTable[246].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[246].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[246].fNumComponents        = 3;
    fNISTMaterialDataTable[246].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[246].fElementList          = ElemsZ_246;
    fNISTMaterialDataTable[246].fElementFraction      = FractionsZ_246;
    fNISTMaterialDataTable[246].fIsBuiltByAtomCount   = false;


    // NIST_MAT_SILVER_IODIDE ----------------------------------------------------------------
    static const int ElemsZ_247[]        = { 47,  53};
    static const double FractionsZ_247[] = {  1,   1};
    fNISTMaterialDataTable[247].fName                 = "NIST_MAT_SILVER_IODIDE";
    fNISTMaterialDataTable[247].fDensity              = 6.01 * (g/cm3);
    fNISTMaterialDataTable[247].fMeanExcitationEnergy = 543.5 * eV;
    fNISTMaterialDataTable[247].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[247].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[247].fNumComponents        = 2;
    fNISTMaterialDataTable[247].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[247].fElementList          = ElemsZ_247;
    fNISTMaterialDataTable[247].fElementFraction      = FractionsZ_247;
    fNISTMaterialDataTable[247].fIsBuiltByAtomCount   = true;


    // NIST_MAT_SKIN_ICRP --------------------------------------------------------------------
    static const int ElemsZ_248[]        = {  1,   6,   7,   8,  11,  15,
                                             16,  17,  19};
    static const double FractionsZ_248[] = {       0.1,      0.204,      0.042,
                                                 0.645,      0.002,      0.001,
                                                 0.002,      0.003,      0.001};
    fNISTMaterialDataTable[248].fName                 = "NIST_MAT_SKIN_ICRP";
    fNISTMaterialDataTable[248].fDensity              = 1.09 * (g/cm3);
    fNISTMaterialDataTable[248].fMeanExcitationEnergy = 72.7 * eV;
    fNISTMaterialDataTable[248].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[248].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[248].fNumComponents        = 9;
    fNISTMaterialDataTable[248].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[248].fElementList          = ElemsZ_248;
    fNISTMaterialDataTable[248].fElementFraction      = FractionsZ_248;
    fNISTMaterialDataTable[248].fIsBuiltByAtomCount   = false;


    // NIST_MAT_SODIUM_CARBONATE -------------------------------------------------------------
    static const int ElemsZ_249[]        = { 11,   6,   8};
    static const double FractionsZ_249[] = {  2,   1,   3};
    fNISTMaterialDataTable[249].fName                 = "NIST_MAT_SODIUM_CARBONATE";
    fNISTMaterialDataTable[249].fDensity              = 2.532 * (g/cm3);
    fNISTMaterialDataTable[249].fMeanExcitationEnergy = 125 * eV;
    fNISTMaterialDataTable[249].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[249].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[249].fNumComponents        = 3;
    fNISTMaterialDataTable[249].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[249].fElementList          = ElemsZ_249;
    fNISTMaterialDataTable[249].fElementFraction      = FractionsZ_249;
    fNISTMaterialDataTable[249].fIsBuiltByAtomCount   = true;


    // NIST_MAT_SODIUM_IODIDE ----------------------------------------------------------------
    static const int ElemsZ_250[]        = { 11,  53};
    static const double FractionsZ_250[] = {  1,   1};
    fNISTMaterialDataTable[250].fName                 = "NIST_MAT_SODIUM_IODIDE";
    fNISTMaterialDataTable[250].fDensity              = 3.667 * (g/cm3);
    fNISTMaterialDataTable[250].fMeanExcitationEnergy = 452 * eV;
    fNISTMaterialDataTable[250].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[250].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[250].fNumComponents        = 2;
    fNISTMaterialDataTable[250].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[250].fElementList          = ElemsZ_250;
    fNISTMaterialDataTable[250].fElementFraction      = FractionsZ_250;
    fNISTMaterialDataTable[250].fIsBuiltByAtomCount   = true;


    // NIST_MAT_SODIUM_MONOXIDE --------------------------------------------------------------
    static const int ElemsZ_251[]        = { 11,   8};
    static const double FractionsZ_251[] = {  2,   1};
    fNISTMaterialDataTable[251].fName                 = "NIST_MAT_SODIUM_MONOXIDE";
    fNISTMaterialDataTable[251].fDensity              = 2.27 * (g/cm3);
    fNISTMaterialDataTable[251].fMeanExcitationEnergy = 148.8 * eV;
    fNISTMaterialDataTable[251].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[251].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[251].fNumComponents        = 2;
    fNISTMaterialDataTable[251].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[251].fElementList          = ElemsZ_251;
    fNISTMaterialDataTable[251].fElementFraction      = FractionsZ_251;
    fNISTMaterialDataTable[251].fIsBuiltByAtomCount   = true;


    // NIST_MAT_SODIUM_NITRATE ---------------------------------------------------------------
    static const int ElemsZ_252[]        = { 11,   7,   8};
    static const double FractionsZ_252[] = {  1,   1,   3};
    fNISTMaterialDataTable[252].fName                 = "NIST_MAT_SODIUM_NITRATE";
    fNISTMaterialDataTable[252].fDensity              = 2.261 * (g/cm3);
    fNISTMaterialDataTable[252].fMeanExcitationEnergy = 114.6 * eV;
    fNISTMaterialDataTable[252].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[252].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[252].fNumComponents        = 3;
    fNISTMaterialDataTable[252].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[252].fElementList          = ElemsZ_252;
    fNISTMaterialDataTable[252].fElementFraction      = FractionsZ_252;
    fNISTMaterialDataTable[252].fIsBuiltByAtomCount   = true;


    // NIST_MAT_STILBENE ---------------------------------------------------------------------
    static const int ElemsZ_253[]        = {  6,   1};
    static const double FractionsZ_253[] = { 14,  12};
    fNISTMaterialDataTable[253].fName                 = "NIST_MAT_STILBENE";
    fNISTMaterialDataTable[253].fDensity              = 0.9707 * (g/cm3);
    fNISTMaterialDataTable[253].fMeanExcitationEnergy = 67.7 * eV;
    fNISTMaterialDataTable[253].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[253].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[253].fNumComponents        = 2;
    fNISTMaterialDataTable[253].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[253].fElementList          = ElemsZ_253;
    fNISTMaterialDataTable[253].fElementFraction      = FractionsZ_253;
    fNISTMaterialDataTable[253].fIsBuiltByAtomCount   = true;


    // NIST_MAT_SUCROSE ----------------------------------------------------------------------
    static const int ElemsZ_254[]        = {  6,   1,   8};
    static const double FractionsZ_254[] = { 12,  22,  11};
    fNISTMaterialDataTable[254].fName                 = "NIST_MAT_SUCROSE";
    fNISTMaterialDataTable[254].fDensity              = 1.5805 * (g/cm3);
    fNISTMaterialDataTable[254].fMeanExcitationEnergy = 77.5 * eV;
    fNISTMaterialDataTable[254].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[254].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[254].fNumComponents        = 3;
    fNISTMaterialDataTable[254].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[254].fElementList          = ElemsZ_254;
    fNISTMaterialDataTable[254].fElementFraction      = FractionsZ_254;
    fNISTMaterialDataTable[254].fIsBuiltByAtomCount   = true;


    // NIST_MAT_TERPHENYL --------------------------------------------------------------------
    static const int ElemsZ_255[]        = {  6,   1};
    static const double FractionsZ_255[] = { 18,  14};
    fNISTMaterialDataTable[255].fName                 = "NIST_MAT_TERPHENYL";
    fNISTMaterialDataTable[255].fDensity              = 1.24 * (g/cm3);
    fNISTMaterialDataTable[255].fMeanExcitationEnergy = 71.7 * eV;
    fNISTMaterialDataTable[255].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[255].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[255].fNumComponents        = 2;
    fNISTMaterialDataTable[255].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[255].fElementList          = ElemsZ_255;
    fNISTMaterialDataTable[255].fElementFraction      = FractionsZ_255;
    fNISTMaterialDataTable[255].fIsBuiltByAtomCount   = true;


    // NIST_MAT_TESTIS_ICRP ------------------------------------------------------------------
    static const int ElemsZ_256[]        = {  1,   6,   7,   8,  11,  15,
                                             16,  17,  19};
    static const double FractionsZ_256[] = {     0.106,      0.099,       0.02,
                                                 0.766,      0.002,      0.001,
                                                 0.002,      0.002,      0.002};
    fNISTMaterialDataTable[256].fName                 = "NIST_MAT_TESTIS_ICRP";
    fNISTMaterialDataTable[256].fDensity              = 1.04 * (g/cm3);
    fNISTMaterialDataTable[256].fMeanExcitationEnergy = 75 * eV;
    fNISTMaterialDataTable[256].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[256].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[256].fNumComponents        = 9;
    fNISTMaterialDataTable[256].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[256].fElementList          = ElemsZ_256;
    fNISTMaterialDataTable[256].fElementFraction      = FractionsZ_256;
    fNISTMaterialDataTable[256].fIsBuiltByAtomCount   = false;


    // NIST_MAT_TETRACHLOROETHYLENE ----------------------------------------------------------
    static const int ElemsZ_257[]        = {  6,  17};
    static const double FractionsZ_257[] = {  2,   4};
    fNISTMaterialDataTable[257].fName                 = "NIST_MAT_TETRACHLOROETHYLENE";
    fNISTMaterialDataTable[257].fDensity              = 1.625 * (g/cm3);
    fNISTMaterialDataTable[257].fMeanExcitationEnergy = 159.2 * eV;
    fNISTMaterialDataTable[257].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[257].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[257].fNumComponents        = 2;
    fNISTMaterialDataTable[257].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[257].fElementList          = ElemsZ_257;
    fNISTMaterialDataTable[257].fElementFraction      = FractionsZ_257;
    fNISTMaterialDataTable[257].fIsBuiltByAtomCount   = true;


    // NIST_MAT_THALLIUM_CHLORIDE ------------------------------------------------------------
    static const int ElemsZ_258[]        = { 81,  17};
    static const double FractionsZ_258[] = {  1,   1};
    fNISTMaterialDataTable[258].fName                 = "NIST_MAT_THALLIUM_CHLORIDE";
    fNISTMaterialDataTable[258].fDensity              = 7.004 * (g/cm3);
    fNISTMaterialDataTable[258].fMeanExcitationEnergy = 690.3 * eV;
    fNISTMaterialDataTable[258].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[258].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[258].fNumComponents        = 2;
    fNISTMaterialDataTable[258].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[258].fElementList          = ElemsZ_258;
    fNISTMaterialDataTable[258].fElementFraction      = FractionsZ_258;
    fNISTMaterialDataTable[258].fIsBuiltByAtomCount   = true;


    // NIST_MAT_TISSUE_SOFT_ICRP -------------------------------------------------------------
    // TISSUE_SOFT_MALE ICRU-44/46 (1989)
    static const int ElemsZ_259[]        = {  1,   6,   7,   8,  11,  15,
                                             16,  17,  19};
    static const double FractionsZ_259[] = {     0.105,      0.256,      0.027,
                                                 0.602,      0.001,      0.002,
                                                 0.003,      0.002,      0.002};
    fNISTMaterialDataTable[259].fName                 = "NIST_MAT_TISSUE_SOFT_ICRP";
    fNISTMaterialDataTable[259].fDensity              = 1.03 * (g/cm3);
    fNISTMaterialDataTable[259].fMeanExcitationEnergy = 72.3 * eV;
    fNISTMaterialDataTable[259].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[259].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[259].fNumComponents        = 9;
    fNISTMaterialDataTable[259].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[259].fElementList          = ElemsZ_259;
    fNISTMaterialDataTable[259].fElementFraction      = FractionsZ_259;
    fNISTMaterialDataTable[259].fIsBuiltByAtomCount   = false;


    // NIST_MAT_TISSUE_SOFT_ICRU-4 -----------------------------------------------------------
    // Tissue soft adult ICRU-33 (1980)
    static const int ElemsZ_260[]        = {  1,   6,   7,   8};
    static const double FractionsZ_260[] = {     0.101,      0.111,      0.026,
                                                 0.762};
    fNISTMaterialDataTable[260].fName                 = "NIST_MAT_TISSUE_SOFT_ICRU-4";
    fNISTMaterialDataTable[260].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[260].fMeanExcitationEnergy = 74.9 * eV;
    fNISTMaterialDataTable[260].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[260].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[260].fNumComponents        = 4;
    fNISTMaterialDataTable[260].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[260].fElementList          = ElemsZ_260;
    fNISTMaterialDataTable[260].fElementFraction      = FractionsZ_260;
    fNISTMaterialDataTable[260].fIsBuiltByAtomCount   = false;


    // NIST_MAT_TISSUE-METHANE ---------------------------------------------------------------
    static const int ElemsZ_261[]        = {  1,   6,   7,   8};
    static const double FractionsZ_261[] = {  0.101869,   0.456179,   0.035172,
                                               0.40678};
    fNISTMaterialDataTable[261].fName                 = "NIST_MAT_TISSUE-METHANE";
    fNISTMaterialDataTable[261].fDensity              = 0.00106409 * (g/cm3);
    fNISTMaterialDataTable[261].fMeanExcitationEnergy = 61.2 * eV;
    fNISTMaterialDataTable[261].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[261].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[261].fNumComponents        = 4;
    fNISTMaterialDataTable[261].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[261].fElementList          = ElemsZ_261;
    fNISTMaterialDataTable[261].fElementFraction      = FractionsZ_261;
    fNISTMaterialDataTable[261].fIsBuiltByAtomCount   = false;


    // NIST_MAT_TISSUE-PROPANE ---------------------------------------------------------------
    static const int ElemsZ_262[]        = {  1,   6,   7,   8};
    static const double FractionsZ_262[] = {  0.102672,    0.56894,   0.035022,
                                              0.293366};
    fNISTMaterialDataTable[262].fName                 = "NIST_MAT_TISSUE-PROPANE";
    fNISTMaterialDataTable[262].fDensity              = 0.00182628 * (g/cm3);
    fNISTMaterialDataTable[262].fMeanExcitationEnergy = 59.5 * eV;
    fNISTMaterialDataTable[262].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[262].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[262].fNumComponents        = 4;
    fNISTMaterialDataTable[262].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[262].fElementList          = ElemsZ_262;
    fNISTMaterialDataTable[262].fElementFraction      = FractionsZ_262;
    fNISTMaterialDataTable[262].fIsBuiltByAtomCount   = false;


    // NIST_MAT_TITANIUM_DIOXIDE -------------------------------------------------------------
    static const int ElemsZ_263[]        = { 22,   8};
    static const double FractionsZ_263[] = {  1,   2};
    fNISTMaterialDataTable[263].fName                 = "NIST_MAT_TITANIUM_DIOXIDE";
    fNISTMaterialDataTable[263].fDensity              = 4.26 * (g/cm3);
    fNISTMaterialDataTable[263].fMeanExcitationEnergy = 179.5 * eV;
    fNISTMaterialDataTable[263].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[263].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[263].fNumComponents        = 2;
    fNISTMaterialDataTable[263].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[263].fElementList          = ElemsZ_263;
    fNISTMaterialDataTable[263].fElementFraction      = FractionsZ_263;
    fNISTMaterialDataTable[263].fIsBuiltByAtomCount   = true;


    // NIST_MAT_TOLUENE ----------------------------------------------------------------------
    static const int ElemsZ_264[]        = {  6,   1};
    static const double FractionsZ_264[] = {  7,   8};
    fNISTMaterialDataTable[264].fName                 = "NIST_MAT_TOLUENE";
    fNISTMaterialDataTable[264].fDensity              = 0.8669 * (g/cm3);
    fNISTMaterialDataTable[264].fMeanExcitationEnergy = 62.5 * eV;
    fNISTMaterialDataTable[264].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[264].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[264].fNumComponents        = 2;
    fNISTMaterialDataTable[264].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[264].fElementList          = ElemsZ_264;
    fNISTMaterialDataTable[264].fElementFraction      = FractionsZ_264;
    fNISTMaterialDataTable[264].fIsBuiltByAtomCount   = true;


    // NIST_MAT_TRICHLOROETHYLENE ------------------------------------------------------------
    static const int ElemsZ_265[]        = {  6,   1,  17};
    static const double FractionsZ_265[] = {  2,   1,   3};
    fNISTMaterialDataTable[265].fName                 = "NIST_MAT_TRICHLOROETHYLENE";
    fNISTMaterialDataTable[265].fDensity              = 1.46 * (g/cm3);
    fNISTMaterialDataTable[265].fMeanExcitationEnergy = 148.1 * eV;
    fNISTMaterialDataTable[265].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[265].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[265].fNumComponents        = 3;
    fNISTMaterialDataTable[265].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[265].fElementList          = ElemsZ_265;
    fNISTMaterialDataTable[265].fElementFraction      = FractionsZ_265;
    fNISTMaterialDataTable[265].fIsBuiltByAtomCount   = true;


    // NIST_MAT_TRIETHYL_PHOSPHATE -----------------------------------------------------------
    static const int ElemsZ_266[]        = {  6,   1,   8,  15};
    static const double FractionsZ_266[] = {  6,  15,   4,   1};
    fNISTMaterialDataTable[266].fName                 = "NIST_MAT_TRIETHYL_PHOSPHATE";
    fNISTMaterialDataTable[266].fDensity              = 1.07 * (g/cm3);
    fNISTMaterialDataTable[266].fMeanExcitationEnergy = 81.2 * eV;
    fNISTMaterialDataTable[266].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[266].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[266].fNumComponents        = 4;
    fNISTMaterialDataTable[266].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[266].fElementList          = ElemsZ_266;
    fNISTMaterialDataTable[266].fElementFraction      = FractionsZ_266;
    fNISTMaterialDataTable[266].fIsBuiltByAtomCount   = true;


    // NIST_MAT_TUNGSTEN_HEXAFLUORIDE --------------------------------------------------------
    static const int ElemsZ_267[]        = { 74,   9};
    static const double FractionsZ_267[] = {  1,   6};
    fNISTMaterialDataTable[267].fName                 = "NIST_MAT_TUNGSTEN_HEXAFLUORIDE";
    fNISTMaterialDataTable[267].fDensity              = 2.4 * (g/cm3);
    fNISTMaterialDataTable[267].fMeanExcitationEnergy = 354.4 * eV;
    fNISTMaterialDataTable[267].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[267].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[267].fNumComponents        = 2;
    fNISTMaterialDataTable[267].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[267].fElementList          = ElemsZ_267;
    fNISTMaterialDataTable[267].fElementFraction      = FractionsZ_267;
    fNISTMaterialDataTable[267].fIsBuiltByAtomCount   = true;


    // NIST_MAT_URANIUM_DICARBIDE ------------------------------------------------------------
    static const int ElemsZ_268[]        = { 92,   6};
    static const double FractionsZ_268[] = {  1,   2};
    fNISTMaterialDataTable[268].fName                 = "NIST_MAT_URANIUM_DICARBIDE";
    fNISTMaterialDataTable[268].fDensity              = 11.28 * (g/cm3);
    fNISTMaterialDataTable[268].fMeanExcitationEnergy = 752 * eV;
    fNISTMaterialDataTable[268].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[268].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[268].fNumComponents        = 2;
    fNISTMaterialDataTable[268].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[268].fElementList          = ElemsZ_268;
    fNISTMaterialDataTable[268].fElementFraction      = FractionsZ_268;
    fNISTMaterialDataTable[268].fIsBuiltByAtomCount   = true;


    // NIST_MAT_URANIUM_MONOCARBIDE ----------------------------------------------------------
    static const int ElemsZ_269[]        = { 92,   6};
    static const double FractionsZ_269[] = {  1,   1};
    fNISTMaterialDataTable[269].fName                 = "NIST_MAT_URANIUM_MONOCARBIDE";
    fNISTMaterialDataTable[269].fDensity              = 13.63 * (g/cm3);
    fNISTMaterialDataTable[269].fMeanExcitationEnergy = 862 * eV;
    fNISTMaterialDataTable[269].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[269].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[269].fNumComponents        = 2;
    fNISTMaterialDataTable[269].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[269].fElementList          = ElemsZ_269;
    fNISTMaterialDataTable[269].fElementFraction      = FractionsZ_269;
    fNISTMaterialDataTable[269].fIsBuiltByAtomCount   = true;


    // NIST_MAT_URANIUM_OXIDE ----------------------------------------------------------------
    static const int ElemsZ_270[]        = { 92,   8};
    static const double FractionsZ_270[] = {  1,   2};
    fNISTMaterialDataTable[270].fName                 = "NIST_MAT_URANIUM_OXIDE";
    fNISTMaterialDataTable[270].fDensity              = 10.96 * (g/cm3);
    fNISTMaterialDataTable[270].fMeanExcitationEnergy = 720.6 * eV;
    fNISTMaterialDataTable[270].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[270].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[270].fNumComponents        = 2;
    fNISTMaterialDataTable[270].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[270].fElementList          = ElemsZ_270;
    fNISTMaterialDataTable[270].fElementFraction      = FractionsZ_270;
    fNISTMaterialDataTable[270].fIsBuiltByAtomCount   = true;


    // NIST_MAT_UREA -------------------------------------------------------------------------
    static const int ElemsZ_271[]        = {  6,   1,   7,   8};
    static const double FractionsZ_271[] = {  1,   4,   2,   1};
    fNISTMaterialDataTable[271].fName                 = "NIST_MAT_UREA";
    fNISTMaterialDataTable[271].fDensity              = 1.323 * (g/cm3);
    fNISTMaterialDataTable[271].fMeanExcitationEnergy = 72.8 * eV;
    fNISTMaterialDataTable[271].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[271].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[271].fNumComponents        = 4;
    fNISTMaterialDataTable[271].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[271].fElementList          = ElemsZ_271;
    fNISTMaterialDataTable[271].fElementFraction      = FractionsZ_271;
    fNISTMaterialDataTable[271].fIsBuiltByAtomCount   = true;


    // NIST_MAT_VALINE -----------------------------------------------------------------------
    static const int ElemsZ_272[]        = {  6,   1,   7,   8};
    static const double FractionsZ_272[] = {  5,  11,   1,   2};
    fNISTMaterialDataTable[272].fName                 = "NIST_MAT_VALINE";
    fNISTMaterialDataTable[272].fDensity              = 1.23 * (g/cm3);
    fNISTMaterialDataTable[272].fMeanExcitationEnergy = 67.7 * eV;
    fNISTMaterialDataTable[272].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[272].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[272].fNumComponents        = 4;
    fNISTMaterialDataTable[272].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[272].fElementList          = ElemsZ_272;
    fNISTMaterialDataTable[272].fElementFraction      = FractionsZ_272;
    fNISTMaterialDataTable[272].fIsBuiltByAtomCount   = true;


    // NIST_MAT_VITON ------------------------------------------------------------------------
    static const int ElemsZ_273[]        = {  1,   6,   9};
    static const double FractionsZ_273[] = {  0.009417,   0.280555,   0.710028};
    fNISTMaterialDataTable[273].fName                 = "NIST_MAT_VITON";
    fNISTMaterialDataTable[273].fDensity              = 1.8 * (g/cm3);
    fNISTMaterialDataTable[273].fMeanExcitationEnergy = 98.6 * eV;
    fNISTMaterialDataTable[273].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[273].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[273].fNumComponents        = 3;
    fNISTMaterialDataTable[273].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[273].fElementList          = ElemsZ_273;
    fNISTMaterialDataTable[273].fElementFraction      = FractionsZ_273;
    fNISTMaterialDataTable[273].fIsBuiltByAtomCount   = false;


    // NIST_MAT_WATER ------------------------------------------------------------------------
    static const int ElemsZ_274[]        = {  1,   8};
    static const double FractionsZ_274[] = {  2,   1};
    fNISTMaterialDataTable[274].fName                 = "NIST_MAT_WATER";
    fNISTMaterialDataTable[274].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[274].fMeanExcitationEnergy = 78 * eV;
    fNISTMaterialDataTable[274].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[274].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[274].fNumComponents        = 2;
    fNISTMaterialDataTable[274].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[274].fElementList          = ElemsZ_274;
    fNISTMaterialDataTable[274].fElementFraction      = FractionsZ_274;
    fNISTMaterialDataTable[274].fIsBuiltByAtomCount   = true;


    // NIST_MAT_WATER_VAPOR ------------------------------------------------------------------
    static const int ElemsZ_275[]        = {  1,   8};
    static const double FractionsZ_275[] = {  2,   1};
    fNISTMaterialDataTable[275].fName                 = "NIST_MAT_WATER_VAPOR";
    fNISTMaterialDataTable[275].fDensity              = 0.000756182 * (g/cm3);
    fNISTMaterialDataTable[275].fMeanExcitationEnergy = 71.6 * eV;
    fNISTMaterialDataTable[275].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[275].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[275].fNumComponents        = 2;
    fNISTMaterialDataTable[275].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[275].fElementList          = ElemsZ_275;
    fNISTMaterialDataTable[275].fElementFraction      = FractionsZ_275;
    fNISTMaterialDataTable[275].fIsBuiltByAtomCount   = true;


    // NIST_MAT_XYLENE -----------------------------------------------------------------------
    static const int ElemsZ_276[]        = {  6,   1};
    static const double FractionsZ_276[] = {  8,  10};
    fNISTMaterialDataTable[276].fName                 = "NIST_MAT_XYLENE";
    fNISTMaterialDataTable[276].fDensity              = 0.87 * (g/cm3);
    fNISTMaterialDataTable[276].fMeanExcitationEnergy = 61.8 * eV;
    fNISTMaterialDataTable[276].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[276].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[276].fNumComponents        = 2;
    fNISTMaterialDataTable[276].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[276].fElementList          = ElemsZ_276;
    fNISTMaterialDataTable[276].fElementFraction      = FractionsZ_276;
    fNISTMaterialDataTable[276].fIsBuiltByAtomCount   = true;


    // NIST_MAT_GRAPHITE ---------------------------------------------------------------------
    static const int ElemsZ_277[]        = {  6};
    static const double FractionsZ_277[] = {  1};
    fNISTMaterialDataTable[277].fName                 = "NIST_MAT_GRAPHITE";
    fNISTMaterialDataTable[277].fDensity              = 2.21 * (g/cm3);
    fNISTMaterialDataTable[277].fMeanExcitationEnergy = 78 * eV;
    fNISTMaterialDataTable[277].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[277].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[277].fNumComponents        = 1;
    fNISTMaterialDataTable[277].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[277].fElementList          = ElemsZ_277;
    fNISTMaterialDataTable[277].fElementFraction      = FractionsZ_277;
    fNISTMaterialDataTable[277].fIsBuiltByAtomCount   = true;



    // =======================================================================================
    // HEP and Nuclear Materials :
    // =======================================================================================

    // NIST_MAT_lH2 --------------------------------------------------------------------------
    static const int ElemsZ_278[]        = {  1};
    static const double FractionsZ_278[] = {  1};
    fNISTMaterialDataTable[278].fName                 = "NIST_MAT_lH2";
    fNISTMaterialDataTable[278].fDensity              = 0.0708 * (g/cm3);
    fNISTMaterialDataTable[278].fMeanExcitationEnergy = 21.8 * eV;
    fNISTMaterialDataTable[278].fTemperature          = 0.0 * kelvin;
    fNISTMaterialDataTable[278].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[278].fNumComponents        = 1;
    fNISTMaterialDataTable[278].fState                = MaterialState::kStateLiquid;
    fNISTMaterialDataTable[278].fElementList          = ElemsZ_278;
    fNISTMaterialDataTable[278].fElementFraction      = FractionsZ_278;
    fNISTMaterialDataTable[278].fIsBuiltByAtomCount   = true;


    // NIST_MAT_lN2 --------------------------------------------------------------------------
    static const int ElemsZ_279[]        = {  7};
    static const double FractionsZ_279[] = {  1};
    fNISTMaterialDataTable[279].fName                 = "NIST_MAT_lN2";
    fNISTMaterialDataTable[279].fDensity              = 0.807 * (g/cm3);
    fNISTMaterialDataTable[279].fMeanExcitationEnergy = 82 * eV;
    fNISTMaterialDataTable[279].fTemperature          = 0.0 * kelvin;
    fNISTMaterialDataTable[279].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[279].fNumComponents        = 1;
    fNISTMaterialDataTable[279].fState                = MaterialState::kStateLiquid;
    fNISTMaterialDataTable[279].fElementList          = ElemsZ_279;
    fNISTMaterialDataTable[279].fElementFraction      = FractionsZ_279;
    fNISTMaterialDataTable[279].fIsBuiltByAtomCount   = true;


    // NIST_MAT_lO2 --------------------------------------------------------------------------
    static const int ElemsZ_280[]        = {  8};
    static const double FractionsZ_280[] = {  1};
    fNISTMaterialDataTable[280].fName                 = "NIST_MAT_lO2";
    fNISTMaterialDataTable[280].fDensity              = 1.141 * (g/cm3);
    fNISTMaterialDataTable[280].fMeanExcitationEnergy = 95 * eV;
    fNISTMaterialDataTable[280].fTemperature          = 0.0 * kelvin;
    fNISTMaterialDataTable[280].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[280].fNumComponents        = 1;
    fNISTMaterialDataTable[280].fState                = MaterialState::kStateLiquid;
    fNISTMaterialDataTable[280].fElementList          = ElemsZ_280;
    fNISTMaterialDataTable[280].fElementFraction      = FractionsZ_280;
    fNISTMaterialDataTable[280].fIsBuiltByAtomCount   = true;


    // NIST_MAT_lAr --------------------------------------------------------------------------
    static const int ElemsZ_281[]        = { 18};
    static const double FractionsZ_281[] = {  1};
    fNISTMaterialDataTable[281].fName                 = "NIST_MAT_lAr";
    fNISTMaterialDataTable[281].fDensity              = 1.396 * (g/cm3);
    fNISTMaterialDataTable[281].fMeanExcitationEnergy = 188 * eV;
    fNISTMaterialDataTable[281].fTemperature          = 0.0 * kelvin;
    fNISTMaterialDataTable[281].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[281].fNumComponents        = 1;
    fNISTMaterialDataTable[281].fState                = MaterialState::kStateLiquid;
    fNISTMaterialDataTable[281].fElementList          = ElemsZ_281;
    fNISTMaterialDataTable[281].fElementFraction      = FractionsZ_281;
    fNISTMaterialDataTable[281].fIsBuiltByAtomCount   = true;


    // NIST_MAT_lBr --------------------------------------------------------------------------
    static const int ElemsZ_282[]        = { 35};
    static const double FractionsZ_282[] = {  1};
    fNISTMaterialDataTable[282].fName                 = "NIST_MAT_lBr";
    fNISTMaterialDataTable[282].fDensity              = 3.1028 * (g/cm3);
    fNISTMaterialDataTable[282].fMeanExcitationEnergy = 343 * eV;
    fNISTMaterialDataTable[282].fTemperature          = 0.0 * kelvin;
    fNISTMaterialDataTable[282].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[282].fNumComponents        = 1;
    fNISTMaterialDataTable[282].fState                = MaterialState::kStateLiquid;
    fNISTMaterialDataTable[282].fElementList          = ElemsZ_282;
    fNISTMaterialDataTable[282].fElementFraction      = FractionsZ_282;
    fNISTMaterialDataTable[282].fIsBuiltByAtomCount   = true;


    // NIST_MAT_lKr --------------------------------------------------------------------------
    static const int ElemsZ_283[]        = { 36};
    static const double FractionsZ_283[] = {  1};
    fNISTMaterialDataTable[283].fName                 = "NIST_MAT_lKr";
    fNISTMaterialDataTable[283].fDensity              = 2.418 * (g/cm3);
    fNISTMaterialDataTable[283].fMeanExcitationEnergy = 352 * eV;
    fNISTMaterialDataTable[283].fTemperature          = 0.0 * kelvin;
    fNISTMaterialDataTable[283].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[283].fNumComponents        = 1;
    fNISTMaterialDataTable[283].fState                = MaterialState::kStateLiquid;
    fNISTMaterialDataTable[283].fElementList          = ElemsZ_283;
    fNISTMaterialDataTable[283].fElementFraction      = FractionsZ_283;
    fNISTMaterialDataTable[283].fIsBuiltByAtomCount   = true;


    // NIST_MAT_lXe --------------------------------------------------------------------------
    static const int ElemsZ_284[]        = { 54};
    static const double FractionsZ_284[] = {  1};
    fNISTMaterialDataTable[284].fName                 = "NIST_MAT_lXe";
    fNISTMaterialDataTable[284].fDensity              = 2.953 * (g/cm3);
    fNISTMaterialDataTable[284].fMeanExcitationEnergy = 482 * eV;
    fNISTMaterialDataTable[284].fTemperature          = 0.0 * kelvin;
    fNISTMaterialDataTable[284].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[284].fNumComponents        = 1;
    fNISTMaterialDataTable[284].fState                = MaterialState::kStateLiquid;
    fNISTMaterialDataTable[284].fElementList          = ElemsZ_284;
    fNISTMaterialDataTable[284].fElementFraction      = FractionsZ_284;
    fNISTMaterialDataTable[284].fIsBuiltByAtomCount   = true;


    // NIST_MAT_PbWO4 ------------------------------------------------------------------------
    static const int ElemsZ_285[]        = {  8,  82,  74};
    static const double FractionsZ_285[] = {  4,   1,   1};
    fNISTMaterialDataTable[285].fName                 = "NIST_MAT_PbWO4";
    fNISTMaterialDataTable[285].fDensity              = 8.28 * (g/cm3);
    fNISTMaterialDataTable[285].fMeanExcitationEnergy = 0 * eV;
    fNISTMaterialDataTable[285].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[285].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[285].fNumComponents        = 3;
    fNISTMaterialDataTable[285].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[285].fElementList          = ElemsZ_285;
    fNISTMaterialDataTable[285].fElementFraction      = FractionsZ_285;
    fNISTMaterialDataTable[285].fIsBuiltByAtomCount   = true;


    // NIST_MAT_Galactic ---------------------------------------------------------------------
    static const int ElemsZ_286[]        = {  1};
    static const double FractionsZ_286[] = {  1};
    fNISTMaterialDataTable[286].fName                 = "NIST_MAT_Galactic";
    fNISTMaterialDataTable[286].fDensity              = kUniverseMeanDensity;
    fNISTMaterialDataTable[286].fMeanExcitationEnergy = 21.8 * eV;
    fNISTMaterialDataTable[286].fTemperature          = 2.73 * kelvin;
    fNISTMaterialDataTable[286].fPressure             = 3.e-18 * pascal;
    fNISTMaterialDataTable[286].fNumComponents        = 1;
    fNISTMaterialDataTable[286].fState                = MaterialState::kStateGas;
    fNISTMaterialDataTable[286].fElementList          = ElemsZ_286;
    fNISTMaterialDataTable[286].fElementFraction      = FractionsZ_286;
    fNISTMaterialDataTable[286].fIsBuiltByAtomCount   = true;


    // NIST_MAT_GRAPHITE_POROUS --------------------------------------------------------------
    static const int ElemsZ_287[]        = {  6};
    static const double FractionsZ_287[] = {  1};
    fNISTMaterialDataTable[287].fName                 = "NIST_MAT_GRAPHITE_POROUS";
    fNISTMaterialDataTable[287].fDensity              = 1.7 * (g/cm3);
    fNISTMaterialDataTable[287].fMeanExcitationEnergy = 78 * eV;
    fNISTMaterialDataTable[287].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[287].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[287].fNumComponents        = 1;
    fNISTMaterialDataTable[287].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[287].fElementList          = ElemsZ_287;
    fNISTMaterialDataTable[287].fElementFraction      = FractionsZ_287;
    fNISTMaterialDataTable[287].fIsBuiltByAtomCount   = true;


    // NIST_MAT_LUCITE -----------------------------------------------------------------------
    // LUCITE is equal to plustiglass
    static const int ElemsZ_288[]        = {  1,   6,   8};
    static const double FractionsZ_288[] = {  0.080538,   0.599848,   0.319614};
    fNISTMaterialDataTable[288].fName                 = "NIST_MAT_LUCITE";
    fNISTMaterialDataTable[288].fDensity              = 1.19 * (g/cm3);
    fNISTMaterialDataTable[288].fMeanExcitationEnergy = 74 * eV;
    fNISTMaterialDataTable[288].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[288].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[288].fNumComponents        = 3;
    fNISTMaterialDataTable[288].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[288].fElementList          = ElemsZ_288;
    fNISTMaterialDataTable[288].fElementFraction      = FractionsZ_288;
    fNISTMaterialDataTable[288].fIsBuiltByAtomCount   = false;


    // NIST_MAT_BRASS ------------------------------------------------------------------------
    // SRIM-2008 materials
    static const int ElemsZ_289[]        = { 29,  30,  82};
    static const double FractionsZ_289[] = { 62,  35,   3};
    fNISTMaterialDataTable[289].fName                 = "NIST_MAT_BRASS";
    fNISTMaterialDataTable[289].fDensity              = 8.52 * (g/cm3);
    fNISTMaterialDataTable[289].fMeanExcitationEnergy = 0 * eV;
    fNISTMaterialDataTable[289].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[289].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[289].fNumComponents        = 3;
    fNISTMaterialDataTable[289].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[289].fElementList          = ElemsZ_289;
    fNISTMaterialDataTable[289].fElementFraction      = FractionsZ_289;
    fNISTMaterialDataTable[289].fIsBuiltByAtomCount   = true;


    // NIST_MAT_BRONZE -----------------------------------------------------------------------
    static const int ElemsZ_290[]        = { 29,  30,  82};
    static const double FractionsZ_290[] = { 89,   9,   2};
    fNISTMaterialDataTable[290].fName                 = "NIST_MAT_BRONZE";
    fNISTMaterialDataTable[290].fDensity              = 8.82 * (g/cm3);
    fNISTMaterialDataTable[290].fMeanExcitationEnergy = 0 * eV;
    fNISTMaterialDataTable[290].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[290].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[290].fNumComponents        = 3;
    fNISTMaterialDataTable[290].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[290].fElementList          = ElemsZ_290;
    fNISTMaterialDataTable[290].fElementFraction      = FractionsZ_290;
    fNISTMaterialDataTable[290].fIsBuiltByAtomCount   = true;


    // NIST_MAT_STAINLESS-STEEL --------------------------------------------------------------
    // parameters are taken from
    // http://www.azom.com/article.aspx?ArticleID=965
    static const int ElemsZ_291[]        = { 26,  24,  28};
    static const double FractionsZ_291[] = { 74,  18,   8};
    fNISTMaterialDataTable[291].fName                 = "NIST_MAT_STAINLESS-STEEL";
    fNISTMaterialDataTable[291].fDensity              = 8 * (g/cm3);
    fNISTMaterialDataTable[291].fMeanExcitationEnergy = 0 * eV;
    fNISTMaterialDataTable[291].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[291].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[291].fNumComponents        = 3;
    fNISTMaterialDataTable[291].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[291].fElementList          = ElemsZ_291;
    fNISTMaterialDataTable[291].fElementFraction      = FractionsZ_291;
    fNISTMaterialDataTable[291].fIsBuiltByAtomCount   = true;


    // NIST_MAT_CR39 -------------------------------------------------------------------------
    static const int ElemsZ_292[]        = {  1,   6,   8};
    static const double FractionsZ_292[] = { 18,  12,   7};
    fNISTMaterialDataTable[292].fName                 = "NIST_MAT_CR39";
    fNISTMaterialDataTable[292].fDensity              = 1.32 * (g/cm3);
    fNISTMaterialDataTable[292].fMeanExcitationEnergy = 0 * eV;
    fNISTMaterialDataTable[292].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[292].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[292].fNumComponents        = 3;
    fNISTMaterialDataTable[292].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[292].fElementList          = ElemsZ_292;
    fNISTMaterialDataTable[292].fElementFraction      = FractionsZ_292;
    fNISTMaterialDataTable[292].fIsBuiltByAtomCount   = true;


    // NIST_MAT_OCTADECANOL ------------------------------------------------------------------
    static const int ElemsZ_293[]        = {  1,   6,   8};
    static const double FractionsZ_293[] = { 38,  18,   1};
    fNISTMaterialDataTable[293].fName                 = "NIST_MAT_OCTADECANOL";
    fNISTMaterialDataTable[293].fDensity              = 0.812 * (g/cm3);
    fNISTMaterialDataTable[293].fMeanExcitationEnergy = 0 * eV;
    fNISTMaterialDataTable[293].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[293].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[293].fNumComponents        = 3;
    fNISTMaterialDataTable[293].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[293].fElementList          = ElemsZ_293;
    fNISTMaterialDataTable[293].fElementFraction      = FractionsZ_293;
    fNISTMaterialDataTable[293].fIsBuiltByAtomCount   = true;



    // =======================================================================================
    // Space Science Materials :
    // =======================================================================================

    // NIST_MAT_KEVLAR -----------------------------------------------------------------------
    static const int ElemsZ_294[]        = {  6,   1,   8,   7};
    static const double FractionsZ_294[] = { 14,  10,   2,   2};
    fNISTMaterialDataTable[294].fName                 = "NIST_MAT_KEVLAR";
    fNISTMaterialDataTable[294].fDensity              = 1.44 * (g/cm3);
    fNISTMaterialDataTable[294].fMeanExcitationEnergy = 0 * eV;
    fNISTMaterialDataTable[294].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[294].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[294].fNumComponents        = 4;
    fNISTMaterialDataTable[294].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[294].fElementList          = ElemsZ_294;
    fNISTMaterialDataTable[294].fElementFraction      = FractionsZ_294;
    fNISTMaterialDataTable[294].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DACRON -----------------------------------------------------------------------
    // POLYETHYLENE_TEREPHTALATE
    static const int ElemsZ_295[]        = {  6,   1,   8};
    static const double FractionsZ_295[] = { 10,   8,   4};
    fNISTMaterialDataTable[295].fName                 = "NIST_MAT_DACRON";
    fNISTMaterialDataTable[295].fDensity              = 1.4 * (g/cm3);
    fNISTMaterialDataTable[295].fMeanExcitationEnergy = 0 * eV;
    fNISTMaterialDataTable[295].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[295].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[295].fNumComponents        = 3;
    fNISTMaterialDataTable[295].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[295].fElementList          = ElemsZ_295;
    fNISTMaterialDataTable[295].fElementFraction      = FractionsZ_295;
    fNISTMaterialDataTable[295].fIsBuiltByAtomCount   = true;


    // NIST_MAT_NEOPRENE ---------------------------------------------------------------------
    // POLYCLOROPRENE
    static const int ElemsZ_296[]        = {  6,   1,  17};
    static const double FractionsZ_296[] = {  4,   5,   1};
    fNISTMaterialDataTable[296].fName                 = "NIST_MAT_NEOPRENE";
    fNISTMaterialDataTable[296].fDensity              = 1.23 * (g/cm3);
    fNISTMaterialDataTable[296].fMeanExcitationEnergy = 0 * eV;
    fNISTMaterialDataTable[296].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[296].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[296].fNumComponents        = 3;
    fNISTMaterialDataTable[296].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[296].fElementList          = ElemsZ_296;
    fNISTMaterialDataTable[296].fElementFraction      = FractionsZ_296;
    fNISTMaterialDataTable[296].fIsBuiltByAtomCount   = true;



    // =======================================================================================
    // Biochemical Materials :
    // =======================================================================================

    // NIST_MAT_CYTOSINE ---------------------------------------------------------------------
    static const int ElemsZ_297[]        = {  1,   6,   7,   8};
    static const double FractionsZ_297[] = {  5,   4,   3,   1};
    fNISTMaterialDataTable[297].fName                 = "NIST_MAT_CYTOSINE";
    fNISTMaterialDataTable[297].fDensity              = 1.55 * (g/cm3);
    fNISTMaterialDataTable[297].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[297].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[297].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[297].fNumComponents        = 4;
    fNISTMaterialDataTable[297].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[297].fElementList          = ElemsZ_297;
    fNISTMaterialDataTable[297].fElementFraction      = FractionsZ_297;
    fNISTMaterialDataTable[297].fIsBuiltByAtomCount   = true;


    // NIST_MAT_THYMINE ----------------------------------------------------------------------
    static const int ElemsZ_298[]        = {  1,   6,   7,   8};
    static const double FractionsZ_298[] = {  6,   5,   2,   2};
    fNISTMaterialDataTable[298].fName                 = "NIST_MAT_THYMINE";
    fNISTMaterialDataTable[298].fDensity              = 1.23 * (g/cm3);
    fNISTMaterialDataTable[298].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[298].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[298].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[298].fNumComponents        = 4;
    fNISTMaterialDataTable[298].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[298].fElementList          = ElemsZ_298;
    fNISTMaterialDataTable[298].fElementFraction      = FractionsZ_298;
    fNISTMaterialDataTable[298].fIsBuiltByAtomCount   = true;


    // NIST_MAT_URACIL -----------------------------------------------------------------------
    static const int ElemsZ_299[]        = {  1,   6,   7,   8};
    static const double FractionsZ_299[] = {  4,   4,   2,   2};
    fNISTMaterialDataTable[299].fName                 = "NIST_MAT_URACIL";
    fNISTMaterialDataTable[299].fDensity              = 1.32 * (g/cm3);
    fNISTMaterialDataTable[299].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[299].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[299].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[299].fNumComponents        = 4;
    fNISTMaterialDataTable[299].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[299].fElementList          = ElemsZ_299;
    fNISTMaterialDataTable[299].fElementFraction      = FractionsZ_299;
    fNISTMaterialDataTable[299].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_ADENINE ------------------------------------------------------------------
    // DNA_Nucleobase (Nucleobase-1H)
    static const int ElemsZ_300[]        = {  1,   6,   7};
    static const double FractionsZ_300[] = {  4,   5,   5};
    fNISTMaterialDataTable[300].fName                 = "NIST_MAT_DNA_ADENINE";
    fNISTMaterialDataTable[300].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[300].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[300].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[300].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[300].fNumComponents        = 3;
    fNISTMaterialDataTable[300].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[300].fElementList          = ElemsZ_300;
    fNISTMaterialDataTable[300].fElementFraction      = FractionsZ_300;
    fNISTMaterialDataTable[300].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_GUANINE ------------------------------------------------------------------
    static const int ElemsZ_301[]        = {  1,   6,   7,   8};
    static const double FractionsZ_301[] = {  4,   5,   5,   1};
    fNISTMaterialDataTable[301].fName                 = "NIST_MAT_DNA_GUANINE";
    fNISTMaterialDataTable[301].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[301].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[301].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[301].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[301].fNumComponents        = 4;
    fNISTMaterialDataTable[301].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[301].fElementList          = ElemsZ_301;
    fNISTMaterialDataTable[301].fElementFraction      = FractionsZ_301;
    fNISTMaterialDataTable[301].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_CYTOSINE -----------------------------------------------------------------
    static const int ElemsZ_302[]        = {  1,   6,   7,   8};
    static const double FractionsZ_302[] = {  4,   4,   3,   1};
    fNISTMaterialDataTable[302].fName                 = "NIST_MAT_DNA_CYTOSINE";
    fNISTMaterialDataTable[302].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[302].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[302].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[302].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[302].fNumComponents        = 4;
    fNISTMaterialDataTable[302].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[302].fElementList          = ElemsZ_302;
    fNISTMaterialDataTable[302].fElementFraction      = FractionsZ_302;
    fNISTMaterialDataTable[302].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_THYMINE ------------------------------------------------------------------
    static const int ElemsZ_303[]        = {  1,   6,   7,   8};
    static const double FractionsZ_303[] = {  5,   5,   2,   2};
    fNISTMaterialDataTable[303].fName                 = "NIST_MAT_DNA_THYMINE";
    fNISTMaterialDataTable[303].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[303].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[303].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[303].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[303].fNumComponents        = 4;
    fNISTMaterialDataTable[303].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[303].fElementList          = ElemsZ_303;
    fNISTMaterialDataTable[303].fElementFraction      = FractionsZ_303;
    fNISTMaterialDataTable[303].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_URACIL -------------------------------------------------------------------
    static const int ElemsZ_304[]        = {  1,   6,   7,   8};
    static const double FractionsZ_304[] = {  3,   4,   2,   2};
    fNISTMaterialDataTable[304].fName                 = "NIST_MAT_DNA_URACIL";
    fNISTMaterialDataTable[304].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[304].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[304].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[304].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[304].fNumComponents        = 4;
    fNISTMaterialDataTable[304].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[304].fElementList          = ElemsZ_304;
    fNISTMaterialDataTable[304].fElementFraction      = FractionsZ_304;
    fNISTMaterialDataTable[304].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_ADENOSINE ----------------------------------------------------------------
    // DNA_Nucleoside (Nucleoside-3H)
    static const int ElemsZ_305[]        = {  1,   6,   7,   8};
    static const double FractionsZ_305[] = { 10,  10,   5,   4};
    fNISTMaterialDataTable[305].fName                 = "NIST_MAT_DNA_ADENOSINE";
    fNISTMaterialDataTable[305].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[305].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[305].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[305].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[305].fNumComponents        = 4;
    fNISTMaterialDataTable[305].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[305].fElementList          = ElemsZ_305;
    fNISTMaterialDataTable[305].fElementFraction      = FractionsZ_305;
    fNISTMaterialDataTable[305].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_GUANOSINE ----------------------------------------------------------------
    static const int ElemsZ_306[]        = {  1,   6,   7,   8};
    static const double FractionsZ_306[] = { 10,  10,   5,   5};
    fNISTMaterialDataTable[306].fName                 = "NIST_MAT_DNA_GUANOSINE";
    fNISTMaterialDataTable[306].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[306].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[306].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[306].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[306].fNumComponents        = 4;
    fNISTMaterialDataTable[306].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[306].fElementList          = ElemsZ_306;
    fNISTMaterialDataTable[306].fElementFraction      = FractionsZ_306;
    fNISTMaterialDataTable[306].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_CYTIDINE -----------------------------------------------------------------
    static const int ElemsZ_307[]        = {  1,   6,   7,   8};
    static const double FractionsZ_307[] = { 10,   9,   3,   5};
    fNISTMaterialDataTable[307].fName                 = "NIST_MAT_DNA_CYTIDINE";
    fNISTMaterialDataTable[307].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[307].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[307].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[307].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[307].fNumComponents        = 4;
    fNISTMaterialDataTable[307].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[307].fElementList          = ElemsZ_307;
    fNISTMaterialDataTable[307].fElementFraction      = FractionsZ_307;
    fNISTMaterialDataTable[307].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_URIDINE ------------------------------------------------------------------
    static const int ElemsZ_308[]        = {  1,   6,   7,   8};
    static const double FractionsZ_308[] = {  9,   9,   2,   6};
    fNISTMaterialDataTable[308].fName                 = "NIST_MAT_DNA_URIDINE";
    fNISTMaterialDataTable[308].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[308].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[308].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[308].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[308].fNumComponents        = 4;
    fNISTMaterialDataTable[308].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[308].fElementList          = ElemsZ_308;
    fNISTMaterialDataTable[308].fElementFraction      = FractionsZ_308;
    fNISTMaterialDataTable[308].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_METHYLURIDINE ------------------------------------------------------------
    static const int ElemsZ_309[]        = {  1,   6,   7,   8};
    static const double FractionsZ_309[] = { 11,  10,   2,   6};
    fNISTMaterialDataTable[309].fName                 = "NIST_MAT_DNA_METHYLURIDINE";
    fNISTMaterialDataTable[309].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[309].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[309].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[309].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[309].fNumComponents        = 4;
    fNISTMaterialDataTable[309].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[309].fElementList          = ElemsZ_309;
    fNISTMaterialDataTable[309].fElementFraction      = FractionsZ_309;
    fNISTMaterialDataTable[309].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_MONOPHOSPHATE ------------------------------------------------------------
    static const int ElemsZ_310[]        = { 15,   8};
    static const double FractionsZ_310[] = {  1,   3};
    fNISTMaterialDataTable[310].fName                 = "NIST_MAT_DNA_MONOPHOSPHATE";
    fNISTMaterialDataTable[310].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[310].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[310].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[310].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[310].fNumComponents        = 2;
    fNISTMaterialDataTable[310].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[310].fElementList          = ElemsZ_310;
    fNISTMaterialDataTable[310].fElementFraction      = FractionsZ_310;
    fNISTMaterialDataTable[310].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_A ------------------------------------------------------------------------
    // Adenine base
    static const int ElemsZ_311[]        = {  1,   6,   7,   8,  15};
    static const double FractionsZ_311[] = { 10,  10,   5,   7,   1};
    fNISTMaterialDataTable[311].fName                 = "NIST_MAT_DNA_A";
    fNISTMaterialDataTable[311].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[311].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[311].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[311].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[311].fNumComponents        = 5;
    fNISTMaterialDataTable[311].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[311].fElementList          = ElemsZ_311;
    fNISTMaterialDataTable[311].fElementFraction      = FractionsZ_311;
    fNISTMaterialDataTable[311].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_G ------------------------------------------------------------------------
    // Guanine base
    static const int ElemsZ_312[]        = {  1,   6,   7,   8,  15};
    static const double FractionsZ_312[] = { 10,  10,   5,   8,   1};
    fNISTMaterialDataTable[312].fName                 = "NIST_MAT_DNA_G";
    fNISTMaterialDataTable[312].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[312].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[312].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[312].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[312].fNumComponents        = 5;
    fNISTMaterialDataTable[312].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[312].fElementList          = ElemsZ_312;
    fNISTMaterialDataTable[312].fElementFraction      = FractionsZ_312;
    fNISTMaterialDataTable[312].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_C ------------------------------------------------------------------------
    // Cytosine base
    static const int ElemsZ_313[]        = {  1,   6,   7,   8,  15};
    static const double FractionsZ_313[] = { 10,   9,   3,   8,   1};
    fNISTMaterialDataTable[313].fName                 = "NIST_MAT_DNA_C";
    fNISTMaterialDataTable[313].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[313].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[313].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[313].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[313].fNumComponents        = 5;
    fNISTMaterialDataTable[313].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[313].fElementList          = ElemsZ_313;
    fNISTMaterialDataTable[313].fElementFraction      = FractionsZ_313;
    fNISTMaterialDataTable[313].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_U ------------------------------------------------------------------------
    // Uracil base
    static const int ElemsZ_314[]        = {  1,   6,   7,   8,  15};
    static const double FractionsZ_314[] = {  9,   9,   2,   9,   1};
    fNISTMaterialDataTable[314].fName                 = "NIST_MAT_DNA_U";
    fNISTMaterialDataTable[314].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[314].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[314].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[314].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[314].fNumComponents        = 5;
    fNISTMaterialDataTable[314].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[314].fElementList          = ElemsZ_314;
    fNISTMaterialDataTable[314].fElementFraction      = FractionsZ_314;
    fNISTMaterialDataTable[314].fIsBuiltByAtomCount   = true;


    // NIST_MAT_DNA_MU -----------------------------------------------------------------------
    // MethaUracil base
    static const int ElemsZ_315[]        = {  1,   6,   7,   8,  15};
    static const double FractionsZ_315[] = { 11,  10,   2,   9,   1};
    fNISTMaterialDataTable[315].fName                 = "NIST_MAT_DNA_MU";
    fNISTMaterialDataTable[315].fDensity              = 1 * (g/cm3);
    fNISTMaterialDataTable[315].fMeanExcitationEnergy = 72 * eV;
    fNISTMaterialDataTable[315].fTemperature          = kNTPTemperature;
    fNISTMaterialDataTable[315].fPressure             = kSTPPressure;
    fNISTMaterialDataTable[315].fNumComponents        = 5;
    fNISTMaterialDataTable[315].fState                = MaterialState::kStateSolid;
    fNISTMaterialDataTable[315].fElementList          = ElemsZ_315;
    fNISTMaterialDataTable[315].fElementFraction      = FractionsZ_315;
    fNISTMaterialDataTable[315].fIsBuiltByAtomCount   = true;
}

} // namespace geantphysics
