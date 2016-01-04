#ifndef _READVECTORDATA_H_
#define _READVECTORDATA_H_


#include "iostream"
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include "base/Global.h"
#include <Vc/Vc>
#include "backend/vc/Backend.h"
#include "backend/Backend.h"

//namespace vecgeom{

class ReadVectorData
{
public: 
    std::vector<double> fRadius, fPhi, fZ, fBr, fBz, fBphi;
    ReadVectorData(std::string inputMap);

    template<class Backend>
    void GatherAlias(typename Backend::Index_t index, typename Backend::Double_t &probNA);
    
    ~ReadVectorData();

private:
    std::string dataFile;
    void PleaseReadData();
};

ReadVectorData::ReadVectorData(std::string inputMap){
        dataFile = inputMap;
        PleaseReadData();
};

ReadVectorData::~ReadVectorData(){

};


void ReadVectorData::PleaseReadData(){
    std::string line;
    std::string s1,s2,s3,s4,s5,s0;
    double d1,d2,d3,d4,d5,d0;
    ifstream pFile(dataFile);
    if (pFile.is_open()){
        while(getline(pFile,line)){
            std::stringstream ss(line);
            ss>> s0>> d1>> s1>> d0>> s2>> d2>> s3>> d3>> s4>> d4>> s5>> d5;
            fRadius.push_back(d1);
            fPhi.push_back(d0);
            fZ.push_back(d2);
            fBz.push_back(d3);
            fBr.push_back(d4);
            fBphi.push_back(d5);
        }
        pFile.close();
    }
    else{
        cout<<"Unable to open file"<<endl;
    }
}


// Scalar method - to be used ONLY for 'scalar-type' backends
//                  i.e. currently: scalar & CUDA

template<class Backend>
void ReadVectorData::GatherAlias(typename Backend::Index_t    index, typename Backend::Double_t  &probNA){
    int     intIndex= (int) index;
    probNA = fBr[intIndex];


}




// Specialisation for all vector-type backends - Vc for now
template<>
void ReadVectorData::GatherAlias<vecgeom::kVc>(typename vecgeom::kVc::Index_t    index, typename vecgeom::kVc::Double_t  &probNA){
  //gather for alias table lookups - (backend type has no ptr arithmetic)
  //vecgeom::kVectorSize == vecphys::kVc::kSize
  //ReadVectorData dataMap("/home/ananya/Work/MagFieldRoutine/cms2015.txt");
  for(int i = 0; i < vecgeom::kVectorSize ; ++i) 
  {
    int ind = index[i];
    //std::cout<<ind<<std::endl;
    if(ind < 0) {
      // printf("Warning: negative index! - in GUPhotoElectronSauterGavrila\n");
      ind = 0;
    }

    //assert( z > 0  && z <= fMaxZelement );
    //    assert( ind >= 0 && ind < fAliasTable[z]->SizeOfGrid() );
    
    //Say index has 5 indices a,b,c,d,e i.e. 5 particles
    //ind=a.  Now for first particle, need to fill probNA with Br[a.. .. e]
    //so say
    probNA[i] = fBz[ind];
    //... but access to fBr needed which we will have obviously
    


    // int tableIndex = fAliasTableManager->GetTableIndex(z);
    // probNA[i]=   (fAliasTableManager->GetAliasTable(tableIndex))->fProbQ[ ind ];
    // aliasInd[i]= (fAliasTableManager->GetAliasTable(tableIndex))->fAlias[ ind ];
  }
};


//}




#endif