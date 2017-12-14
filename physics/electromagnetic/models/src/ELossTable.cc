
#include "ELossTable.h"
// from material
#include "Types.h"

#include "ELossTableRegister.h"
#include "PhysicsParameters.h"
#include "EMPhysicsProcess.h"
#include "Material.h"
#include "MaterialCuts.h"
#include "Particle.h"

#include "Spline.h"
#include "GLIntegral.h"

#include <cmath>
#include <iostream>
#include <iomanip>

namespace geantphysics {

ELossTable::ELossTable(PhysicsParameters *physpar) : fPhysicsParameters(physpar), fEnergyGrid(nullptr) {
  fIsComputeCSDARange    = false;
  fNGL                   =  8;
  fNumLossTableBins      = -1;
  fMinLossTableEnergy    = 0.0;
  fMaxLossTableEnergy    = 0.0;
  fLogMinLossTableEnergy = 1.0;
  fEnergyILDelta         = 0.0;
  fGL                    = nullptr;
}


ELossTable::~ELossTable() {
  // clear all allocated memory and reset vector sizes
  Clear();
  if (fGL) {
    delete fGL;
  }
}


double ELossTable::GetRestrictedDEDX(int matcutindx, int partindx, double kinenergy) {
  double dedx = 0.0;
  // find the ELossData for the Particle and MaterialCuts by their provided indices
  ELossData *lossData = nullptr;
  // should work properly without this check; we need to put it under verbose build
  if (int(fELossDataPerMaterialCutsPerParticle[matcutindx].size())>partindx && (lossData = fELossDataPerMaterialCutsPerParticle[matcutindx][partindx])) {
    // compute lower index of the kinetic energy bin
    if (kinenergy<=fMinLossTableEnergy) {
      dedx = lossData->fRestrictedDEDXData[0];
      // low energy extrapolation
      dedx *= std::sqrt(kinenergy/fMinLossTableEnergy);
    } else if (kinenergy>=fMaxLossTableEnergy) {
      dedx = lossData->fRestrictedDEDXData[fNumLossTableBins-1];
    } else {
      double logE     = std::log(kinenergy);
      int    lowEIndx = (int) ((logE-fLogMinLossTableEnergy)*fEnergyILDelta);
      // we might put it under verbose build since
      // protection against very small numerical uncertainties
//      if (lowEIndx>0 && kinenergy<fEnergyGrid[lowEIndx]) {
//        --lowEIndx;
//      } else if (kinenergy>fEnergyGrid[lowEIndx+1]) {
//        ++lowEIndx;
//      }
      dedx = lossData->fSplineRestrictedDEDX->GetValueAt(kinenergy,lowEIndx);
    }
  } else {
    std::cerr << " ======= ELossTable::GetRestrictedDEDX() no restricted dedx for Particle = "
              << Particle::GetParticleByInternalCode(partindx)->GetName()
              << "  in MaterialCuts: \n";
    std::cerr << MaterialCuts::GetTheMaterialCutsTable()[matcutindx] <<std::endl;
    exit(-1);
  }
  return dedx;
}


double ELossTable::GetRestrictedRange(int matcutindx, int partindx, double kinenergy) {
  double range = 1.0e+20;
  // find the ELossData for the Particle and MaterialCuts by their provided indices
  ELossData *lossData = nullptr;
  // should work properly without this check; we need to put it under verbose build
  if (int(fELossDataPerMaterialCutsPerParticle[matcutindx].size())>partindx && (lossData = fELossDataPerMaterialCutsPerParticle[matcutindx][partindx])) {
    // compute lower index of the kinetic energy bin
    if (kinenergy<=fMinLossTableEnergy) {
      range = lossData->fRestrictedRangeData[0];
      // loww energy extrapolation
      range *= std::sqrt(kinenergy/fMinLossTableEnergy);
    } else if (kinenergy>=fMaxLossTableEnergy) {
      range = lossData->fRestrictedRangeData[fNumLossTableBins-1];
    } else {
      double logE     = std::log(kinenergy);
      int    lowEIndx = (int) ((logE-fLogMinLossTableEnergy)*fEnergyILDelta);
      // we might put it under verbose build since
      // protection against very small numerical uncertainties
//      if (lowEIndx>0 && kinenergy<fEnergyGrid[lowEIndx]) {
//        --lowEIndx;
//      } else if (kinenergy>fEnergyGrid[lowEIndx+1]) {
//        ++lowEIndx;
//      }
      range = lossData->fSplineRestrictedRange->GetValueAt(kinenergy,lowEIndx);
    }
  } else {
    std::cerr << " ======= ELossTable::GetRestrictedRange() no restricted range table for Particle = "
              << Particle::GetParticleByInternalCode(partindx)->GetName()
              << "  in MaterialCuts: \n";
    std::cerr << MaterialCuts::GetTheMaterialCutsTable()[matcutindx] <<std::endl;
    exit(-1); // fatal
  }
  return range;
}


double ELossTable::GetEnergyForRestrictedRange(int matcutindx, int partindx, double range) {
  double energy = 0.0;
  // find the ELossData for the Particle and MaterialCuts by their provided indices
  ELossData *lossData = nullptr;
  // should work properly without this check; we need to put it under verbose build
  if (int(fELossDataPerMaterialCutsPerParticle[matcutindx].size())>partindx && (lossData = fELossDataPerMaterialCutsPerParticle[matcutindx][partindx])) {
    double minRange = lossData->fRestrictedRangeData[0];
    double maxRange = lossData->fRestrictedRangeData[lossData->fNumData-1];
    if (range>=minRange) {
      // if range is higher than the maximum range value in the table set energy to the corresponding max energy
      if (range>=maxRange) {
        energy = fMaxLossTableEnergy;
      } else { // get the corresponding energy by interpolation; NOTE: a binary search will be included if indx = -1
        energy = lossData->fSplineRestrictedInvRange->GetValueAt(range);
      }
    } else if (range>0.0) { // low energy extrapolation
      double dummy = range/minRange;
      energy = fMinLossTableEnergy*dummy*dummy;
    }
  } else {
    std::cerr << " ======= ELossTable::GetEnergyForRestrictedRange() no inverse range table for Particle = "
              << Particle::GetParticleByInternalCode(partindx)->GetName()
              << "  in MaterialCuts: \n";
    std::cerr << MaterialCuts::GetTheMaterialCutsTable()[matcutindx] <<std::endl;
    exit(-1);
  }
  return energy;
}


double ELossTable::GetRange(int matindx, int partindx, double kinenergy) {
  double range = 1.e+20;
  // find the ELossData for the Particle and MaterialCuts by their provided indices
  ELossData *lossData = nullptr;
  // should work properly without this check; we need to put it under verbose build
  if (int(fELossDataPerMaterialPerParticle[matindx].size())>partindx
      && (lossData = fELossDataPerMaterialPerParticle[matindx][partindx])
      && (lossData->fRangeData)) {
    // compute lower index of the kinetic energy bin
    if (kinenergy<=fMinLossTableEnergy) {
      range = lossData->fRangeData[0];
    } else if (kinenergy>=fMaxLossTableEnergy) {
      range = lossData->fRangeData[fNumLossTableBins-1];
    } else {
      double logE     = std::log(kinenergy);
      int    lowEIndx = (int) ((logE-fLogMinLossTableEnergy)*fEnergyILDelta);
      // we might put it under verbose build since
      // protection against very small numerical uncertainties
  //      if (lowEIndx>0 && kinenergy<fEnergyGrid[lowEIndx]) {
  //        --lowEIndx;
  //      } else if (kinenergy>fEnergyGrid[lowEIndx+1]) {
  //        ++lowEIndx;
  //      }
      range = lossData->fSplineRange->GetValueAt(kinenergy,lowEIndx);
    }
  } else {
    std::cerr << " ======= ELossTable::GetRange() no Range data for Particle = "
              << Particle::GetParticleByInternalCode(partindx)->GetName()
              << "  in Material: \n";
    std::cerr << Material::GetTheMaterialTable()[matindx] <<std::endl;
    exit(-1);
  }
  return range;
}


void ELossTable::BuildELossTable(std::vector<ELossTable*> &elosstablespermatcut) {
  // initialize the energy grid
  InitializeEnergyGrid();
  // get some other parameters from the PhysicsParameters object that this table belongs to
  fIsComputeCSDARange = fPhysicsParameters->GetIsComputeCSDARange();
  // get lists of registered EnegyLoss EMPhysicsProcess-es per particle
  const std::vector<std::vector<EMPhysicsProcess*> > eLossProcs = ELossTableRegister::Instance().GetListEnergyLossProcesses();
  // get list of MaterialCuts and Material-s
  const std::vector<MaterialCuts*> theMatCutTable   = MaterialCuts::GetTheMaterialCutsTable();
  const Vector_t<Material*>        theMaterialTable = Material::GetTheMaterialTable();
  const std::vector<bool>          isActiveList     = fPhysicsParameters->GetListActiveRegions();
  int numMutCuts    = theMatCutTable.size();   // number of MaterialCuts in the MaterialCuts table
  int numMaterials  = theMaterialTable.size(); //
  int numParticles  = eLossProcs.size();       // maximum internal code of particle that has EnergyLoss
                                               // EMPhysicsProcess registered in the ELossTableRegister
  // clear all ELossData; normaly everything is clean but make sure;
  Clear();
  // initialise the size of the fELossDataPerMaterialCutsPerParticle and fELossDataPerMaterialPerParticle
  fELossDataPerMaterialCutsPerParticle.resize(numMutCuts);
  if (fIsComputeCSDARange) {
    fELossDataPerMaterialPerParticle.resize(numMaterials);
  }
  //loop over MaterialCuts and take each of those one by one that belongs to the region where the current table is active
  for (int imatcut=0; imatcut<numMutCuts; ++imatcut) {
    const MaterialCuts *matCut = theMatCutTable[imatcut];
    // check if the current MatrialCuts belongs to any of the regions where this ELossTable is active:
    if (!fPhysicsParameters->IsActiveRegion(matCut->GetRegionIndex())) {
      continue;
    }
    // set the elosstablespermatcut pointer for the current MaterialCuts to point this ELossTable since
    // the corresponding MaterialCuts is handled by this ELossTable i.e. the MatrialCuts belongs to a region
    // where this ELossTable is active
    elosstablespermatcut[matCut->GetIndex()] = this;
    // loop over the partciles and for each:
    // 1. check if there were any EnergyLoss EMPhysicsProcess registered for the particle in the ELossTableRegister
    // 2. if yes, collect those of the registered EnergyLoss EMPhysicsProcess-es that are active in the regions where
    //    the current ELossTable is active
    // 3. if there is at least one active EnergyLoss process in this collection then set up the ELossData structure
    //    for this MaterialCuts, partcile and registered EnergyLoss processes
    // 4. check if the total eloss data (i.e. data per material per partcile) has already built for the corresponding
    //    Material and partcile; request to be built if not yet
    // 1.
    for (int ipart=0; ipart<numParticles; ++ipart) {
      // get the current particle by internal partcile code; EnergyLoss EMPhysicsProcess-es were registered by this
      // index in the ELossTableRegister
      const Particle *particle = Particle::GetParticleByInternalCode(ipart);
      // create a collection for possible EnergyLoss processes for this partcile with initial size of zero
      std::vector<EMPhysicsProcess*>  theELossProcessList(0);
      //1.
      if (eLossProcs[ipart].size()>0) {  // if the particle has EnergyLoss process registered
        // 2.
        // loop over the EnergyLoss EMPhysicsProcess-es registered for thsi particle and check activations
        for (unsigned long iproc=0; iproc<eLossProcs[ipart].size(); ++iproc) {
          // check if the current process is active in the same regions where the current ELossTable
          bool isOk   = true;
          for (unsigned long ir=0; ir<isActiveList.size(); ++ir) {
            if (isActiveList[ir] != eLossProcs[ipart][iproc]->IsActiveRegion(ir)) {
              isOk = false;
              break;
            }
          }
          // register the EnergyLoss process for this particle if it is active where the current ELossTable
          if (isOk) {
            theELossProcessList.push_back(eLossProcs[ipart][iproc]);
          }
        }
      }
      // 3.
      if (theELossProcessList.size()>0) {
        ELossData    *lossData                = new ELossData();
        lossData->fNumData                    = fNumLossTableBins;
        lossData->fMaterialCuts               = matCut;
        lossData->fParticle                   = particle;
        lossData->fEnergyGridData             = fEnergyGrid;
        lossData->fRestrictedDEDXData         = nullptr;
        lossData->fRestrictedRangeData        = nullptr;
        lossData->fRangeData                  = nullptr;
        lossData->fSplineRestrictedDEDX       = nullptr;
        lossData->fSplineRestrictedRange      = nullptr;
        lossData->fSplineRestrictedInvRange   = nullptr;
        lossData->fSplineRange                = nullptr;
        // TODO: we do not need later these process pointers so they should not be member of the ELossData
        //       because it is enough just to pass them to the dedx computation
        for (unsigned long ip=0; ip<theELossProcessList.size(); ++ip) {
          (lossData->fLossProcesses).push_back(theELossProcessList[ip]);
        }
        // 4.
        bool isComputeTotalData = false;  // for this material and particle;
        if (fIsComputeCSDARange && fELossDataPerMaterialPerParticle[matCut->GetMaterial()->GetIndex()].size()==0) {
          // this is the first partcile for this material so set the size of the corresponding per particle ELossData*
          // vector to numParticles and init with nullptr-s
          fELossDataPerMaterialPerParticle[matCut->GetMaterial()->GetIndex()].resize(numParticles,nullptr);
          // for sure we need to compute total data as well
          isComputeTotalData = true;
        } else if (fIsComputeCSDARange && !(fELossDataPerMaterialPerParticle[matCut->GetMaterial()->GetIndex()][ipart])) {
          // the ELossData* element, that correspond to the current material and current particle is still nullptr
          // so we need to build the corresponding total data as well
          isComputeTotalData =true;
        }
        //
        // Compute loss data
        //
        BuildOneELossData(lossData,isComputeTotalData);  // will compute dedx, range, inverse range,...
        // set this ELossData structure for the current MatrialCuts in fELossDataPerMaterialCutsPerParticle vector;
        // index will be the Particle internal code that this ELossData belongs to
        // if this is the first ELossData structure to store for this MaterialCuts then set the size to #partciles
        if (fELossDataPerMaterialCutsPerParticle[matCut->GetIndex()].size()==0) {
          fELossDataPerMaterialCutsPerParticle[matCut->GetIndex()].resize(numParticles,nullptr);
          //elosstablespermatcut[matCut->GetIndex()] = this;
        }
        fELossDataPerMaterialCutsPerParticle[matCut->GetIndex()][ipart] = lossData;
        if (isComputeTotalData) {
          fELossDataPerMaterialPerParticle[matCut->GetMaterial()->GetIndex()][ipart] = lossData;
        }
      }
    }  // end loop over the Particle-s
  } // end loop over the MatrialCuts-s
}


//
// private method implementations
//
void ELossTable::InitializeEnergyGrid() {
  fNumLossTableBins          = fPhysicsParameters->GetNumLossTableBins()+1;
  fMinLossTableEnergy        = fPhysicsParameters->GetMinLossTableEnergy();
  fMaxLossTableEnergy        = fPhysicsParameters->GetMaxLossTableEnergy();
  if (!fEnergyGrid) {
    delete [] fEnergyGrid;
    fEnergyGrid = nullptr;
  }
  fEnergyGrid                      = new double[fNumLossTableBins]();
  fLogMinLossTableEnergy           = std::log(fMinLossTableEnergy);
  double delta                     = std::log(fMaxLossTableEnergy/fMinLossTableEnergy)/(fNumLossTableBins-1.0);
  fEnergyILDelta                   = 1.0/delta;
  fEnergyGrid[0]                   = fMinLossTableEnergy;
  fEnergyGrid[fNumLossTableBins-1] = fMaxLossTableEnergy;
  for (int i=1; i<fNumLossTableBins-1; ++i) {
    fEnergyGrid[i] = std::exp(fLogMinLossTableEnergy+i*delta);
  }
  if (!fGL) {
    fGL = new GLIntegral(fNGL,0.0,1.0);
  }
}


void ELossTable::BuildOneELossData(ELossData *lossdata, bool iscomputetotaldata) {
/*
  std::cerr<<"  ---------------------------------------------------------------\n";
  std::cerr<<"  ====  Building ELossData for: \n"
           <<"      particle     = " << (lossdata->fParticle)->GetName() <<"\n"
           <<"      processes    = ";
           for (unsigned long ip=0; ip<(lossdata->fLossProcesses).size(); ++ip)
             std::cerr<< lossdata->fLossProcesses[ip]->GetName() <<"  ";
  std::cerr<<std::endl;
  std::cerr<<"      materialcuts = " << lossdata->fMaterialCuts;
  std::cerr<<std::endl;
  std::cerr<<"   --------- building restricted dxdx .............";
*/
  BuildRestrictedDEDXTable(lossdata);
//  std::cerr<<"  is done! " << std::endl;
//  std::cerr<<"   --------- building restricted range ............";
  BuildRestrictedRangeTable(lossdata);
//  std::cerr<<"  is done! " << std::endl;
  if (iscomputetotaldata) {
//    std::cerr<<"   --------- building total(CSDA) range ........... ";
    BuildTotalRangeTable(lossdata);
//    std::cerr<<" is done! " << std::endl;
  }
}


void ELossTable::BuildRestrictedDEDXTable(ELossData *lossdata) {
  // compute dedx by summing up contributions from each EnergyLoss EMProcesses
  // allocate the space for the dedx
  lossdata->fRestrictedDEDXData = new double[fNumLossTableBins]();
  for (int i=0; i<lossdata->fNumData; ++i) {
    double energy = lossdata->fEnergyGridData[i];
    double dedx   = 0.0;
    for (unsigned long ip=0; ip<(lossdata->fLossProcesses).size(); ++ip) {
      dedx += lossdata->fLossProcesses[ip]->ComputeDEDX(lossdata->fMaterialCuts, energy, lossdata->fParticle);
    }
    if (dedx>0.0) {
      lossdata->fRestrictedDEDXData[i] = dedx;
    }
  }
  // create cubic spline for dedx
  lossdata->fSplineRestrictedDEDX = new Spline(lossdata->fEnergyGridData, lossdata->fRestrictedDEDXData, lossdata->fNumData);
}


void ELossTable::BuildRestrictedRangeTable(ELossData *lossdata) {
  // TODO: - should have some protection against all zero dedx!!!
  //       - if dedx = [0,0,0,x,...] then up to element x range should be set to ??
  //       - if dedx = [0,0,0,x,0,0,...] then range is what for the second set of zeros?
  //
  // the low energy tail i.e. between E=[0,E_0] by assuming that dedx proportional to beta
  // now assume that the i=0 dedx is already not zero; anyway it should be a problem having zero dedx at any low energies
  int    ifirst    = 0;
  double dedx      = lossdata->fRestrictedDEDXData[ifirst];// if it is non zero
  double rangeTail = 2.0*lossdata->fEnergyGridData[ifirst]/dedx;
  // allocate the space for the range
  lossdata->fRestrictedRangeData = new double[fNumLossTableBins]();
  //set first
  lossdata->fRestrictedRangeData[ifirst] = rangeTail;
  // create a 16 point GL integral on [0,1]; integral between [E_i,E_i+1] will be transformed to [0,1]
  const std::vector<double> &glX = fGL->GetAbscissas();
  const std::vector<double> &glW = fGL->GetWeights();
  for (int i=ifirst; i<lossdata->fNumData-1; ++i) {
    // for each E_i, E_i+1 interval apply the 16 point GL by substitution
    double emin  = lossdata->fEnergyGridData[i];
    double emax  = lossdata->fEnergyGridData[i+1];
    double delta = (emax-emin);
    double res   = 0.0;
    for (int j=0; j<fNGL; ++j) {
      double xi = delta*glX[j]+emin;
      dedx = lossdata->fSplineRestrictedDEDX->GetValueAt(xi,i); // i is the low Energy bin index
      if (dedx>0.0) {
        res += glW[j]/dedx;
      }
    }
    res *= delta;
    lossdata->fRestrictedRangeData[i+1] = res+lossdata->fRestrictedRangeData[i];
  }
  // create spline for range interpolation
  lossdata->fSplineRestrictedRange    = new Spline(lossdata->fEnergyGridData, lossdata->fRestrictedRangeData, lossdata->fNumData);
  // create spline for inverse range interpolation;
  lossdata->fSplineRestrictedInvRange = new Spline(lossdata->fRestrictedRangeData, lossdata->fEnergyGridData, lossdata->fNumData);
}


void ELossTable::BuildTotalRangeTable(ELossData *lossdata) {
  // compute dedx by summing up contributions from each EnergyLoss EMProcesses
  // allocate the space for the total dedx
  std::vector<double>  theTotalDEDX(lossdata->fNumData,0.0);
  for (int i=0; i<lossdata->fNumData; ++i) {
    double energy = lossdata->fEnergyGridData[i];
    double dedx   = 0.0;
    for (unsigned long ip=0; ip<(lossdata->fLossProcesses).size(); ++ip) {
      dedx += lossdata->fLossProcesses[ip]->ComputeDEDX(lossdata->fMaterialCuts, energy, lossdata->fParticle, true);
    }
    if (dedx>0.0) {
      theTotalDEDX[i] = dedx;
    }
  }
  // create a spline for the integration
  Spline *sp = new Spline(lossdata->fEnergyGridData, &(theTotalDEDX[0]), lossdata->fNumData);

  // integrate to get the range
  // TODO: - should have some protection against all zero dedx!!!
  //       - if dedx = [0,0,0,x,...] then up to element x range should be set to ??
  //       - if dedx = [0,0,0,x,0,0,...] then range is what for the second set of zeros?
  //
  // the low energy tail i.e. between E=[0,E_0] by assuming that dedx proportional to beta
  // now assume that the i=0 dedx is already not zero; anyway it should be a problem having zero dedx at any low energies
  int    ifirst    = 0;
  double dedx      = theTotalDEDX[ifirst];// if it is non zero
  double rangeTail = 2.0*lossdata->fEnergyGridData[ifirst]/dedx;
  // allocate the space for the range
  lossdata->fRangeData = new double[fNumLossTableBins]();
  //set first
  lossdata->fRangeData[ifirst] = rangeTail;
  // create a 16 point GL integral on [0,1]; integral between [E_i,E_i+1] will be transformed to [0,1]
  const std::vector<double> &glX = fGL->GetAbscissas();
  const std::vector<double> &glW = fGL->GetWeights();
  for (int i=ifirst; i<lossdata->fNumData-1; ++i) {
    // for each E_i, E_i+1 interval apply the 16 point GL by substitution
    double emin  = lossdata->fEnergyGridData[i];
    double emax  = lossdata->fEnergyGridData[i+1];
    double delta = (emax-emin);
    double res   = 0.0;
    for (int j=0; j<fNGL; ++j) {
      double xi = delta*glX[j]+emin;
      dedx = sp->GetValueAt(xi,i); // i is the low Energy bin index
      if (dedx>0.0) {
        res += glW[j]/dedx;
      }
    }
    res *= delta;
    lossdata->fRangeData[i+1] = res+lossdata->fRangeData[i];
  }
  // delete the spline that was set up on the total dedx for the integration
  delete sp;
  // create spline for range interpolation
  lossdata->fSplineRange = new Spline(lossdata->fEnergyGridData, lossdata->fRangeData, lossdata->fNumData);
}


void ELossTable::Clear() {
  for (unsigned long i=0; i<fELossDataPerMaterialCutsPerParticle.size(); ++i) {
    for (unsigned long j=0; j<fELossDataPerMaterialCutsPerParticle[i].size(); ++j) {
      if (fELossDataPerMaterialCutsPerParticle[i][j]) {
        ELossData *lossData = fELossDataPerMaterialCutsPerParticle[i][j];
//        std::cerr<<"  ++++  Deleting ELossData data for \n"
//                 <<"        particle     =  " <<lossData->fParticle->GetName() <<"\n"
//                 <<"        materialcut  = " <<lossData->fMaterialCuts<<std::endl;
        delete [] lossData->fRestrictedDEDXData;
        delete [] lossData->fRestrictedRangeData;
        if (lossData->fRangeData) {    // if total data was also computed
          delete [] lossData->fRangeData;
        }
        if (lossData->fSplineRestrictedDEDX) {
          delete lossData->fSplineRestrictedDEDX;
        }
        if (lossData->fSplineRestrictedRange) {
          delete lossData->fSplineRestrictedRange;
        }
        if (lossData->fSplineRestrictedInvRange) {
          delete lossData->fSplineRestrictedInvRange;
        }
        if (lossData->fSplineRange) {
          delete lossData->fSplineRange;
        }
        delete lossData;
        fELossDataPerMaterialCutsPerParticle[i][j] = nullptr;
      }
    }
    fELossDataPerMaterialCutsPerParticle[i].clear();
  }
  fELossDataPerMaterialCutsPerParticle.clear();

  // clear per material per particle vectors; ELossData objects were already deleted above since pointers are shared
  for (unsigned long i=0; i<fELossDataPerMaterialPerParticle.size(); ++i) {
    fELossDataPerMaterialPerParticle[i].clear();
  }
  fELossDataPerMaterialPerParticle.clear();

  if (!fEnergyGrid) {
    delete [] fEnergyGrid;
    fEnergyGrid = nullptr;
  }
}


} // namespace geantphysics
