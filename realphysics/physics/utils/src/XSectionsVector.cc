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
        } else if(energy >= fBinVector[fNumberOfNodes-2]) {
            id = fNumberOfNodes - 2;
        } //else
        if(idx >= fNumberOfNodes || energy < fBinVector[idx]
           || energy > fBinVector[idx+1])
        {
            // Bin location proposed by K.Genser (FNAL) from G4
            id = std::lower_bound(fBinVector.begin(), fBinVector.end(), energy) - fBinVector.begin() - 1;
            id = std::min(id, fNumberOfNodes-2);
        }
        return id;
    }
    
    //Binary search algorithm to find the binIndex corresponding to 'energy'and retreiving the corresponding linear interpolated value.
    double XSectionsVector::GetValueAt(double energy) const {
        
        size_t idx;
        int    first,last,middle;
        int    upperm2 = fNumberOfNodes-2;
        // check if 'energy' is above/below the highes/lowest value
        if (energy>=fBinVector[upperm2]) { //optimal
            idx = upperm2;
        } else if (energy<=fBinVector[1]) { //optimal
            idx = 0;
        } else {
            // Perform a binary search to find the binIndex corresponding to 'energy'
            first = 1;
            last = upperm2;
            while (std::abs(last-first)>1) {
                middle = (first+last)/2.;
                if (energy<fBinVector[middle])
                    last = middle;
                else
                    first = middle;
            }
            idx = last-1;
        }
        // check if the 2 grid point is 0,0
        if ((fDataVector[0]+idx)+(fDataVector[0]+idx+1)==0.0) {
            return 0.0;
        }
        //std::cout<<"GetValueAt: "<<idx<<"\t"<<fBinVector[idx]<<"\t"<<fDataVector[idx]<<std::endl;
        return fDataVector[idx] +( fDataVector[idx + 1]-fDataVector[idx] ) * (energy - fBinVector[idx]) /( fBinVector[idx + 1]-fBinVector[idx] );

    }
    
    //_____________________________
    //Given an energy, first retrieve the binIndex corresponding to that energy and then calculate the interpolated value (Linear Interpolation) corresponding to the data stored at that bin index
    double XSectionsVector::GetValue(double energy, size_t& shellIdx) const{
        
        if(energy <= fEdgeMin)
        { shellIdx = 0;
            //std::cout<<"GetValue:   "<<shellIdx<<"\t"<<fBinVector[shellIdx]<<"\t"<<fDataVector[shellIdx]<<std::endl;
            return fDataVector[0];}
        if(energy >= fEdgeMax) {
            shellIdx= fNumberOfNodes-2;
            //std::cout<<"GetValue:   "<<shellIdx<<"\t"<<fBinVector[shellIdx]<<"\t"<<fDataVector[shellIdx]<<std::endl;
            return fDataVector[shellIdx];
        }
        shellIdx=FindCSBinLocation(energy, shellIdx);
        //std::cout<<"GetValue:    "<<shellIdx<<"\t"<<fBinVector[shellIdx]<<"\t"<<fDataVector[shellIdx]<<std::endl;
        return LinearInterpolation(energy, shellIdx);

    }
    
    
    // Linear interpolation is used to get the interpolated value for lowEnergy cross sections (below K-shell binding energy).
    //Before this method is called it is ensured that the energy is inside the bin
    double XSectionsVector::LinearInterpolation(double energy, size_t idx) const
    {
        return fDataVector[idx] +( fDataVector[idx + 1]-fDataVector[idx] ) * (energy - fBinVector[idx]) /( fBinVector[idx + 1]-fBinVector[idx] );
    
    }

}//end namespace geantphysics
