

#ifndef TVecUniMF_H
#define TVecUniMF_H

#include <iostream>
#include "base/Vector3D.h"
#include "AlignedBase.h"

class TVecUniMF : public AlignedBase {  
    public:  

        TVecUniMF()
        {
          std::cout<<"\n---- entered TVecUniMF constructor ---"<<std::endl;
        }

        ~TVecUniMF() {}

    private:
        // vecgeom::Vector3D<typename Vc::Vector<float>> fFieldComponents;
        typename Vc::Vector<double> fTestMember; 
        //typename vecgeom::kVc::precision_v fTestMember2;
};

#endif
