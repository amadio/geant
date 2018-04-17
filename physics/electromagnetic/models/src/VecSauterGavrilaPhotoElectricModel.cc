#include <Geant/PhysicalConstants.h>
#include <TString.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include "Geant/VecSauterGavrilaPhotoElectricModel.h"
#include "Geant/AliasTable.h"

namespace geantphysics {

void VecSauterGavrilaPhotoElectricModel::Initialize()
{
  SauterGavrilaPhotoElectricModel::Initialize();
  //some other operations if needed
}
    
    

void VecSauterGavrilaPhotoElectricModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
  
}
PhysDI VecSauterGavrilaPhotoElectricModel::SampleShellAliasVec(PhysDV egamma, PhysDI zed, PhysDV r1, PhysDV r2)
{
    //this will be directly stored in the track
    PhysDV lGammaEnergy_v = vecCore::math::Log(egamma);
    // LOWER bin in the BASE vector
    PhysDI tableIndexBase_v = (PhysDI) ((lGammaEnergy_v - fShellPrimEnLMin) * PhysDV(fShellPrimEnILDelta));
    PhysDI sampledShells;
    //These are static informations that can be passed as an argument in Real_v form - TO DO
    PhysDV kBindingEn_v;
    for (int k=0; k<kPhysDVWidth; k++){
        vecCore::Set(kBindingEn_v, k, fBindingEn[(int)zed[k]][0]);
    }
    PhysDM lowEn (egamma<kBindingEn_v);
    PhysDI tableIndexBinding_v;
    PhysDV baseEn_v(999),bindEn_v(999);
    PhysDI indexBaseEn_v, indexBindingEn_v(-1);
    
    for (int k=0; k<kPhysDVWidth; k++){
        vecCore::Set(indexBaseEn_v, k, fIndexBaseEn[(int)zed[k]][tableIndexBase_v[k]+1]-1);
    }
    PhysDI tableIndex_v(indexBaseEn_v);
    
    //Only the values of tableIndex_v that need to be changed
    if(lowEn.isNotEmpty())
    {
        for (int k=0; k<kPhysDVWidth; k++){
            
            if(lowEn[k]){
                vecCore::Set(baseEn_v, k, fShellSamplingPrimEnergies[tableIndexBase_v[k]+1]); //UPPER VALUE (it could be the last meaningful value in the vector (for those that have the mask set to FALSE)
                //LOWER bin in the BINDING vector
                int tableIndexBinding=std::lower_bound(fSortedDoubledBindingEn[(int)zed[k]].begin(), fSortedDoubledBindingEn[(int)zed[k]].end(), egamma[k]) - fSortedDoubledBindingEn[(int)zed[k]].begin()-1;
                vecCore::Set(bindEn_v, k, fSortedDoubledBindingEn[(int)zed[k]][tableIndexBinding+1]);//UPPER VALUE
                vecCore::Set(indexBindingEn_v, k, fIndexSortedDoubledBindingEn[(int)zed[k]][tableIndexBinding+1]-1);
            }
        }
        PhysDM checkMinVal(baseEn_v>bindEn_v); //If TRUE, take bindingEnergies, otherwise keep the base energies
        vecCore::MaskedAssign (tableIndex_v, checkMinVal&&lowEn, indexBindingEn_v);
    
    }
    PhysDI lastSSAliasIndex_v;
    for (int k=0; k<kPhysDVWidth; k++){
        vecCore::Set(lastSSAliasIndex_v, k, fLastSSAliasIndex[(int)zed[k]-1]);
        
    }
    
    PhysDV val        = (lGammaEnergy_v - fShellPrimEnLMin) * fShellPrimEnILDelta; //To correct - inverse of delta of the log of real en
    // LOWER energy bin index
    PhysDI indxEgamma = (PhysDI)val;
    PhysDV pIndxHigh  = val - indxEgamma;
    PhysDM check(r1<=pIndxHigh);
    vecCore::MaskedAssign (tableIndex_v, check, tableIndex_v+1);
    PhysDI indxTable_v = lastSSAliasIndex_v+tableIndex_v;

    //NB: the SCALAR and the VECTORIZED are almost equivalent
    //SCALAR
    for (int i=0; i<kPhysDVWidth; i++){
        int xsampl = fShellAliasSampler->SampleDiscrete(fShellAliasData[indxTable_v[i]]->fAliasW, fShellAliasData[indxTable_v[i]]->fAliasIndx, fShellAliasData[indxTable_v[i]]->fNumdata, r2[i]);
            vecCore::Set(sampledShells, i, xsampl);
        
    }
    //END SCALAR
        
//    //VECTORIZED
//    Real_v aliasW_v, aliasIndx_v, aliasNumdata_v;
//    //first I need the numData
//    for (size_t kk=0; kk<kRealSize; kk++)
//        vecCore::Set(aliasNumdata_v, kk, fShellAliasData[indxTable_v[kk]]->fNumdata);
//
//    Real_v rest_v  = r2*aliasNumdata_v;
//    Real_v indxBin = vecCore::math::Floor(rest_v);
//    RIndex indxBin_v(indxBin);
//
//    for (size_t kk=0; kk<kRealSize; kk++){
//        vecCore::Set(aliasW_v, kk, fShellAliasData[indxTable_v[kk]]->fAliasW[indxBin_v[kk]]);
//        vecCore::Set(aliasIndx_v, kk, fShellAliasData[indxTable_v[kk]]->fAliasIndx[indxBin_v[kk]]);
//
//    }
//    RMask check3(aliasW_v<rest_v-indxBin_v);
//    vecCore::MaskedAssign(indxBin,check3,aliasIndx_v);
//    //std::cout<<indxBin_v<<std::endl;
//    RIndex temp(indxBin);
//    sampledShells=temp;
//    //END VECTORIZED
    return sampledShells;
}

PhysDV VecSauterGavrilaPhotoElectricModel::SamplePhotoElectronDirectionAliasVec(PhysDV egamma, PhysDV r1, PhysDV r2, PhysDV r3)
{
    PhysDV ecosT(1.);
    PhysDM checkEnergy = egamma > (PhysDV)100 * geant::units::MeV;
    //if some energy is below 100 MeV
    if(!checkEnergy.isFull())
    {
        PhysDV legamma = vecCore::math::Log(egamma);
//        //With this implementation is assumed that the model is never built for energies outside the range where the alias tables where built
//        PhysDV val        = (legamma - fPrimEnLMin) * fPrimEnILDelta;
//        PhysDI gammaEnergyIndx = (PhysDI)val; // lower electron energy bin index
//        PhysDV pIndxHigh  = val - gammaEnergyIndx;
//        PhysDM checkIndex = r1 < pIndxHigh;
//        if (!checkIndex.isEmpty()) {
//            vecCore::MaskedAssign(gammaEnergyIndx, checkIndex, gammaEnergyIndx + 1);
//        }
        
        PhysDI gammaEnergyIndx = (PhysDI)((legamma - fPrimEnLMin) * fPrimEnILDelta);
        PhysDM check1 = gammaEnergyIndx >= fNumSamplingPrimEnergies - 1;
        vecCore::MaskedAssign(gammaEnergyIndx, check1, PhysDI(fNumSamplingPrimEnergies - 2));
        PhysDV fLSamplingPrimEnergies_v;
        for (int i=0; i<kPhysDVWidth; i++)
            vecCore::Set(fLSamplingPrimEnergies_v, i,fLSamplingPrimEnergies[gammaEnergyIndx[i] + 1]);
        PhysDV pLowerGammaEner = (fLSamplingPrimEnergies_v - legamma) * fPrimEnILDelta;
        PhysDM check2 = r1 > pLowerGammaEner;
        vecCore::MaskedAssign(gammaEnergyIndx, check2, gammaEnergyIndx+1);

        //going scalar here
        for (int i=0; i < kPhysDVWidth ; ++i)
        {
            if(egamma[i] < 100 * geant::units::MeV){
                double ecosTheta = fAliasSampler->SampleLinear(fAliasData[gammaEnergyIndx[i]]->fXdata, fAliasData[gammaEnergyIndx[i]]->fYdata, fAliasData[gammaEnergyIndx[i]]->fAliasW,fAliasData[gammaEnergyIndx[i]]->fAliasIndx, fAliasData[gammaEnergyIndx[i]]->fNumdata, r2[i], r3[i]);
                vecCore::Set(ecosT, i, ecosTheta);
            }
        }
    }
    return ecosT;
}
 
void VecSauterGavrilaPhotoElectricModel::SamplePhotoElectronDirectionRejVec(const double *egamma, double *cosTheta, int N, const geant::TaskData *td)
{
    
    int currN        = 0;
    PhysDM lanesDone = PhysDM::Zero(); //no lanes done
    PhysDI idx;
    for (int l = 0; l < kPhysDVWidth; ++l) {
        idx[l] = currN++; //indexes initialization
    }
    while (currN < N || !lanesDone.isFull()) {
     
        // 1) initialize energy-dependent variables
        // Variable naming according to Eq. (2.24) of Penelope Manual
        // (pag. 44)
        //std::cout<<"Gathering: "<<idx<<"\n";
        //std::cout<<"lanesDone: "<<lanesDone<<"\n";
        PhysDV gamma  = 1.0 + vecCore::Gather<PhysDV>(egamma, idx) / geant::units::kElectronMassC2;
        //std::cout<<"gamma: "<<gamma<<"\n";
        PhysDV gamma2 = gamma * gamma;
        PhysDV beta   = std::sqrt((gamma2 - 1.0) / gamma2);
    
        // ac corresponds to "A" of Eq. (2.31)
        //
        PhysDV ac = (1.0 / beta) - 1.0;
        PhysDV a1 = 0.5 * beta * gamma * (gamma - 1.0) * (gamma - 2.0);
        PhysDV a2 = ac + 2.0;
        // gtmax = maximum of the rejection function according to Eq. (2.28), obtained for tsam=0
        PhysDV gtmax = 2.0 * (a1 + 1.0 / ac);
        //std::cout<<"gtmax: "<<gtmax<<"\n";
        PhysDV tsam = 0;
        PhysDV gtr  = 0;
        
        // 2) sampling. Eq. (2.31) of Penelope Manual
        // tsam = 1-std::cos(theta)
        // gtr = rejection function according to Eq. (2.28)
        PhysDV rnd1 = td->fRndm->uniformV(); //here, wasting random numbers  - to be handled for reproducibility issues
        PhysDV rnd2 = td->fRndm->uniformV();
        
        tsam = 2.0 * ac * (2.0 * rnd1 + a2 * std::sqrt(rnd1)) / (a2 * a2 - 4.0 * rnd1);
        //std::cout<<"tsam: "<<tsam<<"\n";
        gtr  = (2.0 - tsam) * (a1 + 1.0 / (ac + tsam));
        PhysDM cond1 = rnd2 * gtmax > gtr;
        PhysDV cosTheta_v = 1.0 - tsam;
        //std::cout<<"cosTheta_v: "<<cosTheta_v<<"\n";
        
        //Scatter anyway all the values, but if cond1 is false the lane will be processed again
        if (cond1.isNotEmpty()) {
            //std::cout<<"Scattering: "<<cosTheta_v<<" on indx: "<<idx<<"\n";
            vecCore::Scatter(cosTheta_v, cosTheta, idx);
        }
        lanesDone = lanesDone || cond1;
        for (int l = 0; l < kPhysDVWidth; ++l) {

            auto laneDone = cond1[l];
            if (laneDone) {
                //std::cout<<"currN: "<<currN<<" and N: "<<N<<std::endl;
                if (currN < N) {
                    idx[l]       = currN++;
                    lanesDone[l] = false;
                    //std::cout<<"idx["<<l<<"]: "<<idx[l]<<" and lanesDone["<<l<<"]: "<<lanesDone[l]<<std::endl;
                } else {
                    idx[l] = N;
                    //std::cout<<"ELSE: idx["<<l<<"]: "<<idx[l]<<" and lanesDone["<<l<<"]: "<<lanesDone[l]<<std::endl;
                }
            }
            //std::cout<<"out1\n";
        }
        //std::cout<<"out2\n";
    }
    //std::cout<<"out3\n";
}

bool VecSauterGavrilaPhotoElectricModel::IsModelUsable(const MaterialCuts *, double ekin)
{
  return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit();
}
}
