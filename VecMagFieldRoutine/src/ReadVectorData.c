#include "ReadVectorData.h"
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>


using namespace std;

ReadVectorData::ReadVectorData(std::string inputMap){
        dataFile = inputMap;
        PleaseReadData();
};

ReadVectorData::~ReadVectorData(){

};

void ReadVectorData::PleaseReadData(){
    string line;
    string s1,s2,s3,s4,s5,s0;
    double d1,d2,d3,d4,d5,d0;
    ifstream pFile(dataFile);
    if (pFile.is_open()){
        while(getline(pFile,line)){
            stringstream ss(line);
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