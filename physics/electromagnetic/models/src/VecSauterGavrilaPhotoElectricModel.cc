#include <Geant/PhysicalConstants.h>
#include <TString.h>
#include <Geant/TaskData.h>
#include <Geant/PhysicsData.h>
#include "Geant/VecSauterGavrilaPhotoElectricModel.h"
#include "Geant/AliasTable.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/Material.h"
#include "Geant/Element.h"
#include "Geant/MaterialProperties.h"

#include "Geant/MaterialCuts.h"
#include "Geant/AliasTable.h"
#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"
#include "Geant/XSectionsVector.h"
#include "Geant/Spline.h"


namespace geantphysics {

void VecSauterGavrilaPhotoElectricModel::Initialize()
{
  SauterGavrilaPhotoElectricModel::Initialize();
  //some other operations if needed
}
    
    

void VecSauterGavrilaPhotoElectricModel::SampleSecondariesVector(LightTrack_v &tracks, geant::TaskData *td)
{
    int N  = tracks.GetNtracks();
    double *kin = tracks.GetKinEArr();
    int zed[N];
    int nshells[N];
    size_t targetElemIndx[N];
    
    int sampledShells[N];
    double cosTheta[N];
    
    for (int i = 0; i < N; ++i) {
        const MaterialCuts *matCut   =  MaterialCuts::GetMaterialCut(tracks.GetMaterialCutCoupleIndex(i));
        const Vector_t<Element *> &theElements = matCut->GetMaterial()->GetElementVector();
        if (theElements.size() > 1) {
            targetElemIndx[i] = SampleTargetElementIndex(matCut, kin[i], td);
            zed[i]=(int)theElements[targetElemIndx[i]]->GetZ();
        }
        else{
            zed[i]=(int)theElements[0]->GetZ();
        }
        nshells[i]=fNShellsUsed[zed[i]];
    }
    if(GetUseSamplingTables()){
        for (int i = 0; i < N; i += kVecLenD) {
            Double_v gammaekin;// = tracks.GetKinEVec(i);
            //MaskD_v chechEnValidity(gammaekin<GetLowEnergyUsageLimit() || gammaekin>GetHighEnergyUsageLimit());
            //if(chechEnValidity.isFull() &&vecCore::EarlyReturnAllowed()) return 0 secondaries -> do nothing
        
            //IndexD_v targetElemIndx_v;
            Double_v zed_v;
            IndexD_v nshells_v;
            //std::cout<<"Load\n";
            vecCore::Load(gammaekin, kin+i);
            //vecCore::Load(zed_v, zed+i);
            vecCore::Load(nshells_v, nshells+i);
        
            IndexD_v z;//(IndexD_v)zed_v;
            vecCore::Load(z, zed+i);
//        MaskD_v elementInitialized() --> to be seen
//        // if element was not initialised, gamma should be absorbed
//        if (!fCrossSectionLE[Z] && !fCrossSection[Z]) {
//            track.SetEnergyDeposit(gammaekin0);
//            // std::cout<<"Model not initialized, Exiting!\n";
//            return 0;
//        }
            //SAMPLING OF THE SHELL WITH ALIAS
            IndexD_v shellIdx(0), sampledShells_v(0), tmp;
            MaskDI_v activateSamplingShells = nshells_v > 1;
            Double_v r1  = td->fRndm->uniformV();
            Double_v r2  = td->fRndm->uniformV();
            (IndexD_v) tmp = SampleShellAliasVec(gammaekin, z, r1, r2);
            vecCore::MaskedAssign(sampledShells_v, activateSamplingShells, tmp);
            vecCore::Store(sampledShells_v, sampledShells+i);
            //SAMPLING OF THE ANGLE WITH ALIAS
            Double_v cosTheta_v;
            MaskDI_v activateSamplingAngle(gammaekin <= 100 * geant::units::MeV);
            if(activateSamplingAngle.isNotEmpty()){
                Double_v r1  = td->fRndm->uniformV();
                Double_v r2  = td->fRndm->uniformV();
                Double_v r3  = td->fRndm->uniformV();
                cosTheta_v = SamplePhotoElectronDirectionAliasVec(gammaekin,r1, r2, r3);
                vecCore::Store(cosTheta_v, cosTheta+i);
            }
        }
       
    }
    else{
        //std::cout<<"Sampling with rejection\n";
        double rands[N];
        //NB: generating random numbers here just for reproducibility issues
        for (int i=0; i<N; i++)
            rands[i]= td->fRndm->uniform();
       // std::cout<<"1\n";
        SampleShellVec(kin, zed, sampledShells, N, td, rands);
       // std::cout<<"2\n";
        SamplePhotoElectronDirectionRejVec(kin, cosTheta, N, td);
        //std::cout<<"3\n";
    
    }
    for (int i = 0; i < N; i += kVecLenD){
        //std::cout<<"4\n";
        Double_v gammaekin_v, cosTheta_v;
        IndexD_v zed_v;
        vecCore::Load(gammaekin_v, kin+i);
        vecCore::Load(zed_v, zed+i);
        vecCore::Load(cosTheta_v, cosTheta+i);
       
        // Retrieving ionized shell bindingEnergy
        Double_v bindingEnergy_v;
        for (int k = 0; k < kVecLenD; ++k) {
            vecCore::Set(bindingEnergy_v, k, (*(fParamHigh[zed_v[k]]))[sampledShells[k+i] * 7 + 1]);
        }

        // Create the secondary particle e-
        Double_v eekin = gammaekin_v - bindingEnergy_v;
        MaskDI_v activateSamplingAngle(gammaekin_v <= 100 * geant::units::MeV);
        
        Double_v eDirX1;
        Double_v eDirY1;
        Double_v eDirZ1;
        if(activateSamplingAngle.isNotEmpty())
        {

            Double_v sinTheta = vecCore::math::Sqrt((1.0 - cosTheta_v) * (1.0 + cosTheta_v));
            Double_v phi  = geant::units::kTwoPi * td->fRndm->uniformV();
        
            // new photoelectron direction in the scattering frame
             eDirX1 = sinTheta * vecCore::math::Cos(phi);
             eDirY1 = sinTheta * vecCore::math::Sin(phi);
             eDirZ1 = cosTheta_v;
        
            Math::RotateToLabFrame(eDirX1, eDirY1, eDirZ1, tracks.GetDirXVec(i), tracks.GetDirYVec(i), tracks.GetDirZVec(i));
        }
        vecCore::MaskedAssign(eDirX1, !activateSamplingAngle, tracks.GetDirXVec(i));
        vecCore::MaskedAssign(eDirY1, !activateSamplingAngle, tracks.GetDirYVec(i));
        vecCore::MaskedAssign(eDirZ1, !activateSamplingAngle, tracks.GetDirZVec(i));
        LightTrack_v &secondaries = td->fPhysicsData->GetSecondarySOA();
        for (int l = 0; l < kVecLenD; ++l) {
            int idx = secondaries.InsertTrack();
            secondaries.SetKinE(eekin[l], idx);
            secondaries.SetDirX(eDirX1[l], idx);
            secondaries.SetDirY(eDirY1[l], idx);
            secondaries.SetDirZ(eDirZ1[l], idx);
            secondaries.SetGVcode(fSecondaryInternalCode, idx); // electron GV code
            secondaries.SetMass(geant::units::kElectronMassC2, idx);
            secondaries.SetTrackIndex(tracks.GetTrackIndex(i + l), idx); // parent Track index
        }
        // update primary track - always kill primary photon
        for (int l = 0; l < kVecLenD; ++l) {
            tracks.SetTrackStatus(LTrackStatus::kKill, i+l);
            tracks.SetKinE(0.0, i + l);
            if (bindingEnergy_v[l] > 0.0) {
                tracks.SetEnergyDeposit(bindingEnergy_v[l], i+l);
            }
        }
        //std::cout<<"stranizza..\n";
    }
    
}

IndexD_v VecSauterGavrilaPhotoElectricModel::SampleShellAliasVec(Double_v egamma, IndexD_v zed, Double_v r1, Double_v r2)
{
    IndexD_v sampledShells(0);
    MaskDI_v enableSamplingShells = (zed != 1) && (zed != 2);
    if(enableSamplingShells.isNotEmpty()){
        //std::cout<<zed<<std::endl;
        //this will be directly stored in the track
        Double_v lGammaEnergy_v = vecCore::math::Log(egamma);
    
        // LOWER bin in the BASE vector
        IndexD_v tableIndexBase_v = (IndexD_v) ((lGammaEnergy_v - fShellPrimEnLMin) * Double_v(fShellPrimEnILDelta));
    
        //These are static informations that can be passed as an argument in Real_v form - TO DO
        Double_v kBindingEn_v;
        for (int k=0; k<kVecLenD; k++){
            vecCore::Set(kBindingEn_v, k, fBindingEn[(int)zed[k]][0]);
        }
        MaskDI_v lowEn (egamma<kBindingEn_v);
        IndexD_v tableIndexBinding_v;
        Double_v baseEn_v(999),bindEn_v(999);
        IndexD_v indexBaseEn_v, indexBindingEn_v(-1);
    
        for (int k=0; k<kVecLenD; k++){
            if(enableSamplingShells[k])
                vecCore::Set(indexBaseEn_v, k, fIndexBaseEn[(int)zed[k]][tableIndexBase_v[k]+1]-1);
            
        }
   
        IndexD_v tableIndex_v(indexBaseEn_v);
    
        //Only the values of tableIndex_v that need to be changed
        if(lowEn.isNotEmpty())
        {
            for (int k=0; k<kVecLenD; k++){
            
                if(lowEn[k]&&enableSamplingShells[k]){
                    vecCore::Set(baseEn_v, k, fShellSamplingPrimEnergies[tableIndexBase_v[k]+1]); //UPPER VALUE (it could be the last meaningful value in the vector (for those that have the mask set to FALSE)
                    
                    //LOWER bin in the BINDING vector
                    int tableIndexBinding=std::lower_bound(fSortedDoubledBindingEn[(int)zed[k]].begin(), fSortedDoubledBindingEn[(int)zed[k]].end(), egamma[k]) - fSortedDoubledBindingEn[(int)zed[k]].begin()-1;
                    vecCore::Set(bindEn_v, k, fSortedDoubledBindingEn[(int)zed[k]][tableIndexBinding+1]);//UPPER VALUE
                    vecCore::Set(indexBindingEn_v, k, fIndexSortedDoubledBindingEn[(int)zed[k]][tableIndexBinding+1]-1);
                    
                }
                
            }
            MaskDI_v checkMinVal(baseEn_v>bindEn_v); //If TRUE, take bindingEnergies, otherwise keep the base energies
            vecCore::MaskedAssign (tableIndex_v, checkMinVal&&lowEn, indexBindingEn_v);
        }
        IndexD_v lastSSAliasIndex_v;
        for (int k=0; k<kVecLenD; k++){
            if(enableSamplingShells[k])
                vecCore::Set(lastSSAliasIndex_v, k, fLastSSAliasIndex[(int)zed[k]-1]);
        }
    
        Double_v val        = (lGammaEnergy_v - fShellPrimEnLMin) * fShellPrimEnILDelta; //To correct - inverse of delta of the log of real en
        // LOWER energy bin index
        IndexD_v indxEgamma = (IndexD_v)val;
        Double_v pIndxHigh  = val - indxEgamma;
        MaskDI_v check(r1<=pIndxHigh);
        vecCore::MaskedAssign (tableIndex_v, check, tableIndex_v+1);
        IndexD_v indxTable_v = lastSSAliasIndex_v+tableIndex_v;

        //NB: the SCALAR and the VECTORIZED are almost equivalent
        //SCALAR
        for (int i=0; i<kVecLenD; i++){
            if(enableSamplingShells[i]){
                int xsampl = fShellAliasSampler->SampleDiscrete(fShellAliasData[indxTable_v[i]]->fAliasW, fShellAliasData[indxTable_v[i]]->fAliasIndx, fShellAliasData[indxTable_v[i]]->fNumdata, r2[i]);
                vecCore::Set(sampledShells, i, xsampl);
            }
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
        
    }
    return sampledShells;
}

Double_v VecSauterGavrilaPhotoElectricModel::SamplePhotoElectronDirectionAliasVec(Double_v egamma, Double_v r1, Double_v r2, Double_v r3)
{
    Double_v ecosT(1.);
    MaskDI_v checkEnergy = egamma > (Double_v)100 * geant::units::MeV;
    //if some energy is below 100 MeV
    if(!checkEnergy.isFull())
    {
        Double_v legamma = vecCore::math::Log(egamma);
//        //With this implementation is assumed that the model is never built for energies outside the range where the alias tables where built
//        Double_v val        = (legamma - fPrimEnLMin) * fPrimEnILDelta;
//        IndexD_v gammaEnergyIndx = (IndexD_v)val; // lower electron energy bin index
//        Double_v pIndxHigh  = val - gammaEnergyIndx;
//        MaskD_v checkIndex = r1 < pIndxHigh;
//        if (!checkIndex.isEmpty()) {
//            vecCore::MaskedAssign(gammaEnergyIndx, checkIndex, gammaEnergyIndx + 1);
//        }
        
        IndexD_v gammaEnergyIndx = (IndexD_v)((legamma - fPrimEnLMin) * fPrimEnILDelta);
        MaskDI_v check1 = (gammaEnergyIndx >= (fNumSamplingPrimEnergies - (Double_v)1));
        vecCore::MaskedAssign(gammaEnergyIndx, check1, IndexD_v(fNumSamplingPrimEnergies - 2));
        Double_v fLSamplingPrimEnergies_v;
        for (int i=0; i<kVecLenD; i++)
            vecCore::Set(fLSamplingPrimEnergies_v, i,fLSamplingPrimEnergies[gammaEnergyIndx[i] + 1]);
        Double_v pLowerGammaEner = (fLSamplingPrimEnergies_v - legamma) * fPrimEnILDelta;
        MaskDI_v check2 = r1 > pLowerGammaEner;
        vecCore::MaskedAssign(gammaEnergyIndx, check2, gammaEnergyIndx+1);

        //going scalar here
        for (int i=0; i < kVecLenD ; ++i)
        {
            if(egamma[i] < 100 * geant::units::MeV){
                double ecosTheta = fAliasSampler->SampleLinear(fAliasData[gammaEnergyIndx[i]]->fXdata, fAliasData[gammaEnergyIndx[i]]->fYdata, fAliasData[gammaEnergyIndx[i]]->fAliasW,fAliasData[gammaEnergyIndx[i]]->fAliasIndx, fAliasData[gammaEnergyIndx[i]]->fNumdata, r2[i], r3[i]);
                vecCore::Set(ecosT, i, ecosTheta);
            }
        }
    }
    return ecosT;
}
    
void VecSauterGavrilaPhotoElectricModel::SampleShellVec(double *egamma, int * zed, int* ss, int N, const geant::TaskData */*td*/, double* randoms)
{

    IndexD_v nshells_v;
    IndexD_v sampledShells_v(0);
    
    std::vector<double> ehep;
    std::vector<double> elep;
    std::vector<double> etab;
    
    std::vector<int> zhep;
    std::vector<int> zlep;
    std::vector<int> ztab;
    
    std::vector<int> indexhep;
    std::vector<int> indexlep;
    std::vector<int> indextab;

    
    std::vector<int> nshellshep;
    std::vector<int> nshellslep;
    std::vector<int> nshellstab;
    
    std::vector<double> randhep;
    std::vector<double> randlep;
    std::vector<double> randtab;
    
    //PREFILTERING
    for (int i=0; i<N; ++i){
       
        if(egamma[i] >= (*(fParamHigh[zed[i]]))[0]){
            ehep.push_back(egamma[i]);
            zhep.push_back(zed[i]);
            nshellshep.push_back(fNShells[zed[i]]);
            indexhep.push_back(i);
            randhep.push_back(randoms[i]);
        }else if(egamma[i] >= (*(fParamLow[zed[i]]))[0]){
            elep.push_back(egamma[i]);
            zlep.push_back(zed[i]);
            nshellslep.push_back(fNShells[zed[i]]);
            indexlep.push_back(i);
            randlep.push_back(randoms[i]);
        }
        else {
            etab.push_back(egamma[i]);
            ztab.push_back(zed[i]);
            indextab.push_back(i);
            randtab.push_back(randoms[i]);
        }
    }
    
    //**** PROCESS THE LEP
    size_t currlep       = 0;
    MaskDI_v lanesDonelep = MaskDI_v::Zero(); //no lanes done
    IndexD_v idxlep;
    
    for (int l = 0; l < kVecLenD; ++l) {
        idxlep[l] = currlep++; //indexes initialization
    }
    IndexD_v idxForLoop(0);
    int sampledShellslep[elep.size()];
    
    while (elep.size()>3 && (currlep < elep.size() || !lanesDonelep.isFull())) {
        IndexD_v zeds     = vecCore::Gather<IndexD_v>(zlep.data(), idxlep);
        Double_v egamma_v = vecCore::Gather<Double_v>(elep.data(), idxlep);
        Double_v iegamma  = (geant::units::MeV)/egamma_v;
        Double_v iegamma2 = iegamma*iegamma;
        Double_v iegamma3 = iegamma2*iegamma;
        Double_v iegamma4 = iegamma2*iegamma2;
        Double_v iegamma5 = iegamma2*iegamma3;
        
        IndexD_v totShells =vecCore::Gather<IndexD_v>(nshellslep.data(), idxlep);
        IndexD_v idx_v   = totShells * 7 - 5;
        
        Double_v pm1, p0, p1, p2, p3, p4, p5;
        
        for(int k=0; k<kVecLenD && !lanesDonelep[k]; k++)
        {
            if(!lanesDonelep[k]){
                //std::cout<<"uea\n";
                vecCore::Set(p0, k, (*(fParamLow[zeds[k]]))[idx_v[k]]);
                vecCore::Set(p1, k, (*(fParamLow[zeds[k]]))[idx_v[k]+1]);
                vecCore::Set(p2, k, (*(fParamLow[zeds[k]]))[idx_v[k]+2]);
                vecCore::Set(p3, k, (*(fParamLow[zeds[k]]))[idx_v[k]+3]);
                vecCore::Set(p4, k, (*(fParamLow[zeds[k]]))[idx_v[k]+4]);
                vecCore::Set(p5, k, (*(fParamLow[zeds[k]]))[idx_v[k]+5]);
                
            }
        }
        
        Double_v rand_v=vecCore::Gather<Double_v>(randlep.data(), idxlep);
        //td->fRndm->uniformV();
        Double_v cs0 = rand_v * (p0+iegamma*p1+iegamma2*p2+iegamma3*p3+iegamma4*p4+iegamma5*p5);
        //std::cout<<"Calculated cs0: \t"<<cs0<<std::endl;
        Double_v idxShells = idxForLoop*7 + 2;
        for(int k=0; k<kVecLenD ; k++)
        {
            if(!lanesDonelep[k]){
                vecCore::Set(pm1, k, (*(fParamLow[zeds[k]]))[idxShells[k]]-1);
                vecCore::Set(p0, k, (*(fParamLow[zeds[k]]))[idxShells[k]]);
                vecCore::Set(p1, k, (*(fParamLow[zeds[k]]))[idxShells[k]+1]);
                vecCore::Set(p2, k, (*(fParamLow[zeds[k]]))[idxShells[k]+2]);
                vecCore::Set(p3, k, (*(fParamLow[zeds[k]]))[idxShells[k]+3]);
                vecCore::Set(p4, k, (*(fParamLow[zeds[k]]))[idxShells[k]+4]);
                vecCore::Set(p5, k, (*(fParamLow[zeds[k]]))[idxShells[k]+5]);
                
            }
        }
        
        MaskDI_v checkKinE(egamma_v > pm1);
        Double_v cs = p0+iegamma*p1+iegamma2*p2+iegamma3*p3+iegamma4*p4+iegamma5*p5;
        //std::cout<<"Calculated cs: \t"<<cs<<std::endl;
        MaskDI_v accepted(cs>=cs0);
        //std::cout<<"accepted: \t"<<accepted<<std::endl;
        MaskDI_v lastShell(idxForLoop==totShells-1);
        MaskDI_v checkOut(accepted||lastShell);
        vecCore::MaskedAssign(idxForLoop, !checkOut, idxForLoop+1);
        //I could scatter directly to the original indexes
        vecCore::Scatter(idxForLoop, sampledShellslep, idxlep);
        vecCore::MaskedAssign(idxForLoop, checkOut, (IndexD_v)0);
        lanesDonelep = lanesDonelep || checkOut;
        for (int l = 0; l < kVecLenD; ++l) {
            auto laneDone = checkOut[l];
            if (laneDone) {
                //std::cout<<"Lane "<<l<<"was accepted"<<std::endl;
                //std::cout<<"Sampled shell is "<<idxForLoop[l]<<std::endl;
                if (currlep < elep.size()) {
                    idxlep[l]       = currlep++;
                    //std::cout<<"idxlep[l] "<<idxlep[l]<<", e si ricomincia."<<std::endl;
                    lanesDonelep[l] = false;
                } else {
                    idxlep[l] = elep.size();
                }
            }
        }
    }
    //std::cout<<"Sampled shells for lep: "<<elep.size()<<std::endl;
    for(size_t i=0; i<elep.size(); i++){
        //std::cout<<i<<"\t"<<sampledShellslep[i]<<std::endl;
        ss[indexlep[i]]=sampledShellslep[i];
    }
    
    
    //**** PROCESS THE HEP
    size_t currhep       = 0;
    MaskDI_v lanesDonehep = MaskDI_v::Zero(); //no lanes done
    IndexD_v idxhep;
    
    for (int l = 0; l < kVecLenD; ++l) {
        idxhep[l] = currhep++; //indexes initialization
    }
    idxForLoop=0;
    int sampledShellshep[ehep.size()];
    
    while (ehep.size()>3 && (currhep < ehep.size() || !lanesDonehep.isFull())) {
        
        //std::cout<<"currhep: "<<currhep<<" and ehep.size(): "<<ehep.size()<<std::endl;
        //std::cout<<"Current indexes: \t"<<idxhep<<std::endl;
        IndexD_v zeds     = vecCore::Gather<IndexD_v>(zhep.data(), idxhep);
        Double_v egamma_v = vecCore::Gather<Double_v>(ehep.data(), idxhep);
        Double_v rand_v   = vecCore::Gather<Double_v>(randhep.data(), idxhep);

        Double_v iegamma  = (geant::units::MeV)/egamma_v;
        Double_v iegamma2 = iegamma*iegamma;
        Double_v iegamma3 = iegamma2*iegamma;
        Double_v iegamma4 = iegamma2*iegamma2;
        Double_v iegamma5 = iegamma2*iegamma3;
        IndexD_v totShells =vecCore::Gather<IndexD_v>(nshellshep.data(), idxhep);
        IndexD_v idx_v   = totShells * 7 - 5;

//        std::cout<<"zeds: "<<zeds<<"\n";
//        std::cout<<"totShells: "<<totShells<<"\n";
//        std::cout<<"rand_v: "<<rand_v<<"\n";
        Double_v pm1, p0, p1, p2, p3, p4, p5;
        
        for(int k=0; k<kVecLenD; k++)
        {
            if(!lanesDonehep[k]){
                vecCore::Set(p0, k, (*(fParamHigh[zeds[k]]))[idx_v[k]]);
                vecCore::Set(p1, k, (*(fParamHigh[zeds[k]]))[idx_v[k]+1]);
                vecCore::Set(p2, k, (*(fParamHigh[zeds[k]]))[idx_v[k]+2]);
                vecCore::Set(p3, k, (*(fParamHigh[zeds[k]]))[idx_v[k]+3]);
                vecCore::Set(p4, k, (*(fParamHigh[zeds[k]]))[idx_v[k]+4]);
                vecCore::Set(p5, k, (*(fParamHigh[zeds[k]]))[idx_v[k]+5]);
                
            }
        }
        //td->fRndm->uniformV();
        
        Double_v cs0 = rand_v * (p0+iegamma*p1+iegamma2*p2+iegamma3*p3+iegamma4*p4+iegamma5*p5);
        //std::cout<<"Calculated cs0: \t"<<cs0<<std::endl;
        
        Double_v idxShells = idxForLoop*7 + 2;
        for(int k=0; k<kVecLenD; k++)
        {
            if(!lanesDonehep[k]){
                vecCore::Set(pm1, k, (*(fParamHigh[zeds[k]]))[idxShells[k]-1]);
                vecCore::Set(p0, k, (*(fParamHigh[zeds[k]]))[idxShells[k]]);
                vecCore::Set(p1, k, (*(fParamHigh[zeds[k]]))[idxShells[k]+1]);
                vecCore::Set(p2, k, (*(fParamHigh[zeds[k]]))[idxShells[k]+2]);
                vecCore::Set(p3, k, (*(fParamHigh[zeds[k]]))[idxShells[k]+3]);
                vecCore::Set(p4, k, (*(fParamHigh[zeds[k]]))[idxShells[k]+4]);
                vecCore::Set(p5, k, (*(fParamHigh[zeds[k]]))[idxShells[k]+5]);
            }
        }
        
        MaskDI_v checkKinE(egamma_v > pm1);
        Double_v cs = p0+iegamma*p1+iegamma2*p2+iegamma3*p3+iegamma4*p4+iegamma5*p5;
        //std::cout<<"Calculated cs: \t"<<cs<<std::endl;
        MaskDI_v accepted(cs>=cs0);
        //std::cout<<"accepted: \t"<<accepted<<std::endl;
        MaskDI_v lastShell(idxForLoop==totShells-1);
        MaskDI_v checkOut(accepted||lastShell);
        vecCore::MaskedAssign(idxForLoop, !checkOut, idxForLoop+1);
        //I could scatter directly to the original indexes
        vecCore::Scatter(idxForLoop, sampledShellshep, idxhep);
        vecCore::MaskedAssign(idxForLoop, checkOut, (IndexD_v)0);
        lanesDonehep = lanesDonehep || checkOut;
        for (int l = 0; l < kVecLenD; ++l) {
            auto laneDone = checkOut[l];
            if (laneDone) {
                //std::cout<<"Lane "<<l<<" was accepted"<<std::endl;
                //std::cout<<"Sampled shell is "<<idxForLoop[l]<<std::endl;
                if (currhep < ehep.size()) {
                    idxhep[l]       = currhep++;
                    //std::cout<<"idxhep[l] "<<idxhep[l]<<", e si ricomincia."<<std::endl;
                    lanesDonehep[l] = false;
                } else {
                    idxhep[l] = ehep.size();
                }
            }
        }
    }
    //std::cout<<"Sampled shells for hep: "<<ehep.size()<<std::endl;
    for(size_t i=0; i<ehep.size(); i++){
        //std::cout<<i<<"\t"<<sampledShellshep[i]<<std::endl;
        ss[indexhep[i]]=sampledShellshep[i];
    }

    //**** PROCESS THE LET && HET in old SCALAR MODE
    for(size_t i=0; i<etab.size(); i++){
        size_t nshells = fNShells[ztab[i]];
        size_t shellIdx=0;
        if(nshells > 1){
            double cs = randtab[i];
            size_t idx= 0;
            // (***) Tabulated values above k-shell ionization energy
            if(etab[i] >= (*(fParamHigh[ztab[i]]))[1]) {
                idx=fCSVector[ztab[i]]->FindCSBinLocation(etab[i], idx);
                cs*=fCSVector[ztab[i]]->fSplineInt->GetValueAt(etab[i], idx);
            }
            //(****) Tabulated values below k-shell ionization energy
            else{
                cs*=fLECSVector[ztab[i]]->GetValue(etab[i], idx);
                
            }
            for(size_t j=0; j<nshells; ++j){
                shellIdx=(size_t)fShellVector[ztab[i]][j]->fCompID;
                if(etab[i] > (*(fParamLow[ztab[i]]))[7*shellIdx+1]) {
                    size_t idx=0;
                    cs-=fShellVector[ztab[i]][j]->GetValue(etab[i], idx);
                    
                }
                if(cs <= 0.0 || j+1 == nshells)
                break;
            }
            ss[indextab[i]]=shellIdx;
        }
        
    }
//    std::cout<<"*******"<<std::endl;
//    std::cout<<"N: "<<N<<std::endl;
//    std::cout<<"lep.size(): "<<elep.size()<<std::endl;
//    std::cout<<"hep.size(): "<<ehep.size()<<std::endl;
//    std::cout<<"tab.size(): "<<etab.size()<<std::endl;

}
 
void VecSauterGavrilaPhotoElectricModel::SamplePhotoElectronDirectionRejVec(const double *egamma, double *cosTheta, int N, const geant::TaskData *td)
{
    
    int currN        = 0;
    MaskDI_v lanesDone = MaskDI_v::Zero(); //no lanes done
    IndexD_v idx;
    for (int l = 0; l < kVecLenD; ++l) {
        idx[l] = currN++; //indexes initialization
    }
    while (currN < N || !lanesDone.isFull()) {
     
        // 1) initialize energy-dependent variables
        // Variable naming according to Eq. (2.24) of Penelope Manual
        // (pag. 44)
        Double_v gamma  = 1.0 + vecCore::Gather<Double_v>(egamma, idx) / geant::units::kElectronMassC2;
        Double_v gamma2 = gamma * gamma;
        Double_v beta   = std::sqrt((gamma2 - 1.0) / gamma2);
    
        // ac corresponds to "A" of Eq. (2.31)
        //
        Double_v ac = (1.0 / beta) - 1.0;
        Double_v a1 = 0.5 * beta * gamma * (gamma - 1.0) * (gamma - 2.0);
        Double_v a2 = ac + 2.0;
        // gtmax = maximum of the rejection function according to Eq. (2.28), obtained for tsam=0
        Double_v gtmax = 2.0 * (a1 + 1.0 / ac);
        Double_v tsam = 0;
        Double_v gtr  = 0;
        
        // 2) sampling. Eq. (2.31) of Penelope Manual
        // tsam = 1-std::cos(theta)
        // gtr = rejection function according to Eq. (2.28)
        Double_v rnd1 = td->fRndm->uniformV(); //here, wasting random numbers  - to be handled for reproducibility issues
        Double_v rnd2 = td->fRndm->uniformV();
        
        tsam = 2.0 * ac * (2.0 * rnd1 + a2 * std::sqrt(rnd1)) / (a2 * a2 - 4.0 * rnd1);
        gtr  = (2.0 - tsam) * (a1 + 1.0 / (ac + tsam));
        MaskDI_v cond1 = rnd2 * gtmax > gtr;
        Double_v cosTheta_v = 1.0 - tsam;
        
        //Scatter anyway all the values, but if cond1 is false the lane will be processed again
        if (cond1.isNotEmpty()) {
            vecCore::Scatter(cosTheta_v, cosTheta, idx);
        }
        lanesDone = lanesDone || cond1;
        for (int l = 0; l < kVecLenD; ++l) {

            auto laneDone = cond1[l];
            if (laneDone) {
                if (currN < N) {
                    idx[l]       = currN++;
                    lanesDone[l] = false;
                } else {
                    idx[l] = N;
                }
            }
        }
    }
}

bool VecSauterGavrilaPhotoElectricModel::IsModelUsable(const MaterialCuts *, double ekin)
{
  return ekin < GetHighEnergyUsageLimit() && ekin > GetLowEnergyUsageLimit();
}
}
