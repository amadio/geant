#ifndef _STRUCTS_H_
#define _STRUCTS_H_


typedef float dataType;

struct MagVector{
public:
    dataType Br   =0.;
    dataType Bz   =0.;
    dataType Bphi =0.;
public:
    void SetBr(dataType a){
        Br = a;
    }

    void SetBz(dataType a){
        Bz = a;
    }

    void SetBphi(dataType a){
        Bphi=a;
    }

    dataType GetBr(){
        return Br;
    }

    dataType GetBz(){
        return Bz;
    }

    dataType GetBphi(){
        return Bphi;
    }
};


struct MagCellStructs{
public:
    MagVector m0,m1, m2,m3;
};


struct MagCellArrays{
    dataType sBr[4];
    dataType sBphi[4];
    dataType sBz[4];
};

#endif