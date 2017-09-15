

#ifndef TVecUniMF_H
#define TVecUniMF_H

#include <iostream>
#include <base/Vector3D.h>
#include <base/AlignedBase.h>
#include <Geant/VectorTypes.h>

class TVecUniMF : public vecgeom::AlignedBase {  
    public:  

        TVecUniMF()
        {
          std::cout<<"\n---- entered TVecUniMF constructor ---"<<std::endl;
        }

        ~TVecUniMF() {}

    private:
        Geant::Double_v fTestMember; 
};

#endif
