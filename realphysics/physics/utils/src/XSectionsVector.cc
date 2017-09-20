#include <XSectionsVector.h>

using namespace std;
namespace geantphysics {

    XSectionsVector::XSectionsVector()
    {
        fBinVector = std::vector<double>();
        fDataVector = std::vector<double>();
        
        fBinVector.clear();
        fDataVector.clear();
        
    }
    
    XSectionsVector::~XSectionsVector(){}
    
    //Given and energy, retrieve the bin index corresponding to that energy
    size_t XSectionsVector::FindCSBinLocation(double energy, size_t idx)const{
        
        size_t id = idx;
        if(energy < fBinVector[1]) {
            id=0;
            //return std::min(bin, numberOfNodes-2);
        } else if(energy >= fBinVector[numberOfNodes-2]) {
            id = numberOfNodes - 2;
        } //else
        if(idx >= numberOfNodes || energy < fBinVector[idx]
           || energy > fBinVector[idx+1])
        {
            // Bin location proposed by K.Genser (FNAL) from G4
            id = std::lower_bound(fBinVector.begin(), fBinVector.end(), energy) - fBinVector.begin() - 1;
            id = std::min(id, numberOfNodes-2);
        }
        return id;
    }
    
    
    
    //_____________________________
    //Given an energy, first retrieve the binIndex corresponding to that energy and then calculate the interpolated value (Linear Interpolation) corresponding to the data stored at that bin index
    double XSectionsVector::GetValue(double energy, size_t& shellIdx){
        
        
        if(energy <= edgeMin)
        { shellIdx = 0; return fDataVector[0];}
        if(energy >= edgeMax) {
            shellIdx= numberOfNodes-1;
            return fDataVector[shellIdx];
        }
        shellIdx=FindCSBinLocation(energy, shellIdx);
        return LinearInterpolation(energy, shellIdx);

    }
    
    
    // Linear interpolation is used to get the interpolated value for lowEnergy cross sections (below K-shell binding energy).
    //Before this method is called it is ensured that the energy is inside the bin
    double XSectionsVector::LinearInterpolation(double energy, size_t idx)
    {
        return fDataVector[idx] +( fDataVector[idx + 1]-fDataVector[idx] ) * (energy - fBinVector[idx]) /( fBinVector[idx + 1]-fBinVector[idx] );
    
    }

}
