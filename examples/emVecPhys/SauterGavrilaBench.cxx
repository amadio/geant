#include <benchmark/benchmark.h>

#include <VecCore/VecCore>

#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <math.h>

#include <getopt.h>
#include <err.h>

#include "Geant/Material.h"
#include "Geant/Element.h"
#include "Geant/MaterialCuts.h"

// vecgeom includes
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"

#include "Geant/Region.h"
#include "Geant/PhysicsParameters.h"

// just to clear them
#include "Geant/ELossTableManager.h"
#include "Geant/ELossTableRegister.h"

#include "Geant/Particle.h"
#include "Geant/Electron.h"
#include "Geant/Positron.h"
#include "Geant/Gamma.h"

#include "Geant/EMModel.h"
#include "Geant/SauterGavrilaPhotoElectricModel.h"

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

// from geantV
#include "Geant/Typedefs.h"
#include "Geant/TaskData.h"

using namespace geantphysics;

const int kMinBasket = 16;
const int kMaxBasket = 256;

geant::TaskData *PrepareTaskData()
{
    auto PhysData    = new geantphysics::PhysicsData();
    auto Td          = new geant::TaskData(1, kMaxBasket);
    Td->fPhysicsData = PhysData;
    return Td;
}

class EMModelTest : public geantphysics::EMModel{
public:
    
    using Real_v = vecCore::backend::VcVector::Real_v;
    using Real = double;
    //    using Real_v = vecCore::backend::VcVector::Float_v;
    //    using Real = float;
    EMModel *emModel;
    const MaterialCuts *matCut;
    Particle * particle;
    int kBasketSize = 4*10;
    Real kPi = 3.1415;
    
    
    EMModelTest() : geantphysics::EMModel("test") {
        //std::cout<<"Costruisco!\n";
        this->setUpSimulation();
        this->setUpModel();
    }
    ~EMModelTest(){
        delete emModel;
        //std::cout<<"tutto distrutto e deallocato!\n";
        //delete matCut;
    }
    void RotateToLabFrameScalar(std::vector<Real> &u, std::vector<Real> &v, std::vector<Real> &w,
                                std::vector<Real> u1, std::vector<Real> u2, std::vector<Real> u3) {
        for(size_t i = 0; i < u.size(); i++){
            this->RotateToLabFrame(u[i],v[i],w[i],u1[i],u2[i],u3[i]);
        }
        
    }
    
    void setUpSimulation(){
        
        // default values of the input parameters
        static std::string   particleName("gamma");                  // primary particle is gamma
        static std::string   materialName("NIST_MAT_Pb");         // material is lead
        static std::string   photoElectricModelName("SauterGavrilaPhotoElectric");      // name of the sauterGavrilaPhotoElectric model to test
        //static int           numHistBins       = 100;             // number of histogram bins between min/max values
        //static int           numSamples        = 1.e+06;           // number of required final state samples - 1.e+07;
        //static Real        primaryEnergy     = 0.01;             // primary particle energy in [GeV] - 0.01=10MeV
        static Real        prodCutValue      = 0.1;             // by default in length and internal units i.e. [cm]
        static bool          isProdCutInLength = true;            // is the production cut value given in length ?
        //static bool          isAlias           = false;            // is the Alias sampling active ?
        
        
        //============================= Set user defined input data =================================//
        // Create target material: which is supposed to be a NIST Material
        Material *matDetector = Material::NISTMaterial(materialName);
        //
        // Set particle kinetic energy
        //Real kineticEnergy    = primaryEnergy;
        //
        // Set production cuts if needed
        bool iscutinlength    = isProdCutInLength;
        Real gcut             = prodCutValue;
        Real emcut            = prodCutValue;
        Real epcut            = prodCutValue;
        //===========================================================================================//
        
        // Create primary particle
        //Particle *
        particle = Gamma::Definition();
        std::string  pname;
        
        
        //============= Initialization i.e. building up and init the physics ========================//
        // Create a dummy vecgeom::geometry:
        //  - with only one volume i.e. world
        //  - create a region and set production cuts
        //
        // create a vecgeom::LogicalVolume
        vecgeom::UnplacedBox worldParams = vecgeom::UnplacedBox(1.,1.,1.);
        vecgeom::LogicalVolume  worldl(&worldParams);
        // create one region and assigne to the logical volume
        vecgeom::Region *aRegion = new vecgeom::Region("ARegion",iscutinlength, gcut, emcut, epcut);
        worldl.SetRegion(aRegion);
        // set the material pointer in the world logical volume
        worldl.SetMaterialPtr((void*)matDetector);
        vecgeom::GeoManager::Instance().SetWorld(worldl.Place());
        //vecgeom::GeoManager::Instance().CloseGeometry();
        // Create all(we have only one) MaterialCuts
        MaterialCuts::CreateAll();
        //===========================================================================================//
        
        // if primary particle energy < gamma production cut => there is no secondary gamma production
        // So get the MaterialCuts of the target: we have only one
        //const MaterialCuts *
        matCut = MaterialCuts::GetMaterialCut(aRegion->GetIndex(),matDetector->GetIndex());
        //===========================================================================================//
        
    }
    
    
    void setUpModel(){
        //std::cout<<"Creating the model SauterGavrilaPhotoElectricModel\n";
        emModel = new SauterGavrilaPhotoElectricModel("SauterGavrilaPhotoElectricModel", true);
        // - Set low/high energy usage limits to their min/max possible values
        emModel->SetLowEnergyUsageLimit ( 0.01*geant::units::keV);
        emModel->SetHighEnergyUsageLimit(100.0*geant::units::GeV);
        
        //=========== Set the active regions of the model and one physics-parameter object ==========//
        // - Set the model to be active in region index 0
        (emModel->GetListActiveRegions()).resize(1); // one region
        (emModel->GetListActiveRegions())[0] = true; // make it active there
        // - Create one PhysicsParameters object (with defult values)
        PhysicsParameters *physPars = new PhysicsParameters();
        // - Set it to be active in region index 0
        (physPars->GetListActiveRegions()).resize(1);
        (physPars->GetListActiveRegions())[0] = true;
        // - Initialisation of the model
        emModel->Initialize();
        //===========================================================================================//
        
    }
    
    int GetVectorSize(){
        return vecCore::VectorSize<Real_v>();
    }
    
    
    
    
    //Starts from vectors of Reals
    void Alias_vector(std::vector<Real> &kinE, std::vector<Real> &eDirX, std::vector<Real> &eDirY, std::vector<Real> &eDirZ){
        
        // Set up a dummy geant::TaskData and its geantphysics::PhysicsData member: they are needed in the final state
        // sampling
        geant::TaskData *td = new geant::TaskData(1,1);
        PhysicsData *phd = new PhysicsData();
        td->fPhysicsData = phd;
        double *rndArray = td->fDblArray;
        size_t vecsize = (kinE.size() / vecCore::VectorSize<Real_v>())*vecCore::VectorSize<Real_v>() ;
        size_t i=0;
        Real_v rand_v, rand2_v, rand3_v, rand4_v, kinE_v;
        
        for(; i < vecsize; i+=vecCore::VectorSize<Real_v>() ){
            
            td->fRndm->uniform_array(vecCore::VectorSize<Real_v>(), rndArray);
            Real *real_rndArray= (Real*)rndArray;
            vecCore::Load(rand_v, real_rndArray);
            //vecCore::Load(rand_v, rndArray);
            td->fRndm->uniform_array(vecCore::VectorSize<Real_v>(), rndArray);
            real_rndArray= (Real*)rndArray;
            vecCore::Load(rand2_v, real_rndArray);
            //vecCore::Load(rand2_v, rndArray);
            td->fRndm->uniform_array(vecCore::VectorSize<Real_v>(), rndArray);
            real_rndArray= (Real*)rndArray;
            vecCore::Load(rand3_v, real_rndArray);
            //vecCore::Load(rand3_v, rndArray);
            td->fRndm->uniform_array(vecCore::VectorSize<Real_v>(), rndArray);
            real_rndArray= (Real*)rndArray;
            vecCore::Load(rand4_v, real_rndArray);
            //vecCore::Load(rand4_v, rndArray);
            vecCore::Load(kinE_v, &kinE[i]);
            
            Real_v eDirX1_v,  eDirY1_v , eDirZ1_v;
            //this->emModel->SamplePhotoElectronDirection_Alias_v (kinE_v, rand_v, rand2_v, rand3_v, rand4_v, eDirX1_v, eDirY1_v, eDirZ1_v);
            vecCore::Store(eDirX1_v, &eDirX[i]);
            vecCore::Store(eDirY1_v, &eDirY[i]);
            vecCore::Store(eDirZ1_v, &eDirZ[i]);
        }
        
        //leftOver part
        for(;i<kinE.size();i++){
            td->fRndm->uniform_array(4, rndArray);
            
            Real cosTheta ;//= this->emModel->SamplePhotoElectronDirection_Alias(kinE[i], rndArray[0], rndArray[1], rndArray[2]);
            Real senTheta = std::sqrt((1 - cosTheta)*(1 + cosTheta));
            Real phi = (geant::units::kTwoPi * rndArray[3]);
            eDirX.push_back(senTheta * std::cos(phi));
            eDirY.push_back(senTheta * std::sin(phi));
            eDirZ.push_back(cosTheta);
        }
    }
    
    
    
    //Starts from vectors of Reals
    void Alias_vector_noRandom(std::vector<Real> &kinE, std::vector<Real> &eDirX, std::vector<Real> &eDirY, std::vector<Real> &eDirZ, std::vector<Real> &randoms){
        
        size_t vecsize = (kinE.size() / vecCore::VectorSize<Real_v>())*vecCore::VectorSize<Real_v>() ;
        size_t i=0, k=0;
        Real_v rand_v, rand2_v, rand3_v, rand4_v, kinE_v;
        //std::cout<<"OCCHIO QUI : "<<vecCore::VectorSize<Real_v>()<<std::endl;
        for(; i < vecsize; i+=vecCore::VectorSize<Real_v>() ){
            
            vecCore::Load(kinE_v, &kinE[i]);
            vecCore::Load(rand_v,  randoms.data()+k);
            vecCore::Load(rand2_v, randoms.data()+(k+  vecCore::VectorSize<Real_v>()) );
            vecCore::Load(rand3_v, randoms.data()+(k+2*vecCore::VectorSize<Real_v>()) );
            vecCore::Load(rand4_v, randoms.data()+(k+3*vecCore::VectorSize<Real_v>()) );
            k+=4*vecCore::VectorSize<Real_v>();
            
            Real_v eDirX1_v,  eDirY1_v , eDirZ1_v;
            //this->emModel->SamplePhotoElectronDirection_Alias_v (kinE_v, rand_v, rand2_v, rand3_v, rand4_v, eDirX1_v, eDirY1_v, eDirZ1_v);
            vecCore::Store(eDirX1_v, &eDirX[i]);
            vecCore::Store(eDirY1_v, &eDirY[i]);
            vecCore::Store(eDirZ1_v, &eDirZ[i]);
        }
        
        //leftOver part
        for(;i<kinE.size();i++){
            Real cosTheta;// = this->emModel->SamplePhotoElectronDirection_Alias(kinE[i], randoms[k], randoms[k+1], randoms[k+2]);
            Real senTheta = std::sqrt((1 - cosTheta)*(1 + cosTheta));
            Real phi = (geant::units::kTwoPi * randoms[k+3]);
            eDirX.push_back(senTheta * std::cos(phi));
            eDirY.push_back(senTheta * std::sin(phi));
            eDirZ.push_back(cosTheta);
            k+=4;
        }
    }
    
    void Alias_scalar(std::vector<Real> &kinE, std::vector<Real> &eDirX, std::vector<Real> &eDirY, std::vector<Real> &eDirZ){
        
        // Set up a dummy geant::TaskData and its geantphysics::PhysicsData member: they are needed in the final state
        // sampling
        geant::TaskData *td = new geant::TaskData(1,1);
        PhysicsData *phd = new PhysicsData();
        td->fPhysicsData = phd;
        double *rndArray = td->fDblArray;
        Real cosTheta,senTheta,phi;
        
        for(size_t i = 0; i<kinE.size(); i++){
            
            td->fRndm->uniform_array(4, rndArray);
            cosTheta ;//= this->emModel->SamplePhotoElectronDirection_Alias(kinE[i], rndArray[0], rndArray[1], rndArray[2]);
            senTheta = std::sqrt((1 - cosTheta)*(1 + cosTheta));
            phi = (geant::units::kTwoPi * rndArray[3]);
            eDirX.push_back(senTheta * std::cos(phi));
            eDirY.push_back(senTheta * std::sin(phi));
            eDirZ.push_back(cosTheta);
        }
    }
    
    void Alias_scalar_noRandom(std::vector<Real> &kinE, std::vector<Real> &eDirX, std::vector<Real> &eDirY, std::vector<Real> &eDirZ, std::vector<Real> &randoms){
        
        Real cosTheta,senTheta,phi;
        
        for(size_t i = 0; i<kinE.size(); i++){
            
            size_t k = i*4;
            cosTheta;// = this->emModel->SamplePhotoElectronDirection_Alias(kinE[i], randoms[k], randoms[k+1], randoms[k+2]);
            senTheta = std::sqrt((1 - cosTheta)*(1 + cosTheta));
            phi = (geant::units::kTwoPi * randoms[k+3]);
            eDirX.push_back(senTheta * std::cos(phi));
            eDirY.push_back(senTheta * std::sin(phi));
            eDirZ.push_back(cosTheta);
        }
    }
    void SampleShell_vector(std::vector<Real>&kinE, std::vector<int>&zed, std::vector<int>&shells, std::vector<Real>&randoms, std::vector<int>&sampledShells){
        
        size_t vecsize = (kinE.size() / vecCore::VectorSize<Real_v>())*vecCore::VectorSize<Real_v>() ;
        size_t i=0;
        vecCore::Index_v<Real_v> shells_v, zed_v, sampledShells_v;
        Real_v rand_v, kinE_v;
        
        for(; i < vecsize; i+=vecCore::VectorSize<Real_v>() ){
            
            vecCore::Load(kinE_v, &kinE[i]);
            vecCore::Load(zed_v, &zed[i]);
            vecCore::Load(shells_v, &shells[i]);
            vecCore::Load(rand_v, &randoms[i]);
            //this->emModel->SampleShell_v(kinE_v, zed_v, shells_v, rand_v, sampledShells_v);
            vecCore::Store(sampledShells_v, &sampledShells[i]);
        }
        //exit(-1);
        
        //leftOver part
        for(;i<kinE.size();i++){
            Real sampledShell;
            size_t zeta=zed[i];
            size_t shell=shells[i];
            size_t sampledShella;
            //this->emModel->SampleShell(kinE[i], zeta, shell, randoms[i], sampledShella);
            sampledShells.push_back(sampledShella);
            
        }
    }
    
    void SampleShell_scalar(std::vector<Real>& kinE, std::vector<size_t>&zed, std::vector<size_t>&shells, std::vector<Real>&randoms, std::vector<Real>&sampledShells){
        
        for(size_t i=0;i<kinE.size();i++){
            size_t sampledShell=10; //for the moment
            //this->emModel->SampleShell(kinE[i], zed[i], shells[i], randoms[i], sampledShell);
            //std::cout<<"---> DOING NOTHING: "<<sampledShell<<std::endl;
            sampledShells.push_back(sampledShell);
            
        }
    }
    
    
    void SampleSec_vector(std::vector<LightTrack> lt){
        
        //Set up the model
        // Set up a dummy geant::TaskData and its geantphysics::PhysicsData member: they are needed in the final state
        // sampling
        //        geant::TaskData *td = new geant::TaskData(1,1);
        //        PhysicsData *phd = new PhysicsData();
        //        td->fPhysicsData = phd;
        auto td = PrepareTaskData();
        //this->emModel->SampleSecondariesVectorized (lt,td);
        //leftOver part
        size_t i = ((size_t)lt.size()/vecCore::VectorSize<Real_v>())*vecCore::VectorSize<Real_v>();
        for(;i<lt.size();i++){
            this->emModel->SampleSecondaries(lt[i], td);
        }
        
    }
    void SampleSec_scalar(std::vector<LightTrack> lt){
        // Set up a dummy geant::TaskData and its geantphysics::PhysicsData member: they are needed in the final state
        // sampling
        auto td = PrepareTaskData();
        //        geant::TaskData *td = new geant::TaskData(1,1);
        //        PhysicsData *phd = new PhysicsData();
        //        td->fPhysicsData = phd;
        for(size_t i=0; i<lt.size(); i++)
            this->SampleSecondaries(lt[i], td);
    }
    
    
    void RotateToLabFrameVector(std::vector<Real> &u, std::vector<Real> &v, std::vector<Real> &w,
                                std::vector<Real> u1, std::vector<Real> u2, std::vector<Real> u3) {
        size_t vecsize = (u.size() / vecCore::VectorSize<Real_v>())*vecCore::VectorSize<Real_v>() ;
        size_t i;
        for(i = 0; i < vecsize; i+=vecCore::VectorSize<Real_v>() ){
            Real_v u_v,v_v,w_v,u1_v,u2_v,u3_v;
            vecCore::Load(u_v,&u[i]);
            vecCore::Load(v_v,&v[i]);
            vecCore::Load(w_v,&w[i]);
            vecCore::Load(u1_v,&u1[i]);
            vecCore::Load(u2_v,&u2[i]);
            vecCore::Load(u3_v,&u3[i]);
            this->RotateToLabFrame_v(u_v,v_v,w_v,u1_v,u2_v,u3_v);
            vecCore::Store(u_v,&u[i]);
            vecCore::Store(v_v,&v[i]);
            vecCore::Store(w_v,&w[i]);
        }
        for(;i<u.size();i++){
            this->RotateToLabFrame(u[i],v[i],w[i],u1[i],u2[i],u3[i]);
        }
    }
    
    void RotateToLabFrameVector2(std::vector<Real> &u, std::vector<Real> &v, std::vector<Real> &w,
                                 std::vector<Real> u1, std::vector<Real> u2, std::vector<Real> u3) {
        size_t vecsize = (u.size() / vecCore::VectorSize<Real_v>())*vecCore::VectorSize<Real_v>() ;
        size_t i;
        for(i = 0; i < vecsize; i+=vecCore::VectorSize<Real_v>() ){
            Real_v u_v,v_v,w_v,u1_v,u2_v,u3_v;
            vecCore::Load(u_v,&u[i]);
            vecCore::Load(v_v,&v[i]);
            vecCore::Load(w_v,&w[i]);
            vecCore::Load(u1_v,&u1[i]);
            vecCore::Load(u2_v,&u2[i]);
            vecCore::Load(u3_v,&u3[i]);
            //this->RotateToLabFrame_v2(u_v,v_v,w_v,u1_v,u2_v,u3_v);
            vecCore::Store(u_v,&u[i]);
            vecCore::Store(v_v,&v[i]);
            vecCore::Store(w_v,&w[i]);
        }
        for(;i<u.size();i++){
            this->RotateToLabFrame(u[i],v[i],w[i],u1[i],u2[i],u3[i]);
        }
    }
    
//    void GenerateDirectionVectors(std::vector<Real>& x,std::vector<Real>& y,std::vector<Real>& z,int N){
//        std::default_random_engine generator;
//        std::uniform_real_distribution<Real> distribution(0.0,1.0);
//        for(int i = 0; i < N; i++){
//            Real phi1 = 2*kPi*distribution(generator);
//            Real theta1 = kPi*distribution(generator);
//            x.push_back(std::sin(theta1)*std::cos(phi1));
//            y.push_back(std::sin(theta1)*std::sin(phi1));
//            z.push_back(std::cos(theta1));
//        }
//    }
//
    void GenerateLightTracks(std::vector<geantphysics::LightTrack>& lt, Real primaryEnergy, int N){
        
        Real ekin       = primaryEnergy;
        std::vector<Real> dirx; //     = 0.0;   // direction
        std::vector<Real>  diry;//      = 0.0;
        std::vector<Real>  dirz;//       = 1.0;
        GenerateDirectionVectors(dirx,diry,dirz,N);
        int    gvcode     = particle->GetInternalCode();// internal code of the primary particle i.e. e-
        
        for (long int i=0; i<N; ++i){
            geantphysics::LightTrack primaryLT;
            primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
            primaryLT.SetKinE(ekin);
            primaryLT.SetGVcode(gvcode);
            primaryLT.SetDirX(dirx[i]);
            primaryLT.SetDirY(diry[i]);
            primaryLT.SetDirZ(dirz[i]);
            lt.push_back(primaryLT);
        }
    }
    
    
}; //end of class EMModelTest





/// :::::: START DEFINITION OF BENCHMARKS ::::::: /////

////::: ROTATE
//static void BM_RotateToLabFrameScalar(benchmark::State& state) {
//    for (auto _ : state) {
//        state.PauseTiming();
//        EMModelTest emModelTest;
//        std::vector<Real> dirX,dirY,dirZ,trX,trY,trZ;
//        emModelTest.GenerateDirectionVectors(dirX,dirY,dirZ,state.range(0));
//        emModelTest.GenerateDirectionVectors(trX,trY,trZ,state.range(0));
//        state.ResumeTiming();
//
//        emModelTest.RotateToLabFrameScalar(dirX,dirY,dirZ,trX,trY,trZ);
//    }
//}
//
//BENCHMARK(BM_RotateToLabFrameScalar)->RangeMultiplier(2)->Range(1<<1,1<<10);
//
//static void BM_RotateToLabFrameVector(benchmark::State& state) {
//    for (auto _ : state) {
//        state.PauseTiming();
//        EMModelTest emModelTest;
//        std::vector<Real> dirX,dirY,dirZ,trX,trY,trZ;
//        emModelTest.GenerateDirectionVectors(dirX,dirY,dirZ,state.range(0));
//        emModelTest.GenerateDirectionVectors(trX,trY,trZ,state.range(0));
//        state.ResumeTiming();
//
//        emModelTest.RotateToLabFrameVector(dirX,dirY,dirZ,trX,trY,trZ);
//    }
//}
////
//BENCHMARK(BM_RotateToLabFrameVector)->RangeMultiplier(2)->Range(1<<1,1<<10);

//static void BM_RotateToLabFrameVector2(benchmark::State& state) {
//    for (auto _ : state) {
//        state.PauseTiming();
//        EMModelTest emModelTest;
//        std::vector<Real> dirX,dirY,dirZ,trX,trY,trZ;
//        emModelTest.GenerateDirectionVectors(dirX,dirY,dirZ,state.range(0));
//        emModelTest.GenerateDirectionVectors(trX,trY,trZ,state.range(0));
//        state.ResumeTiming();
//
//        emModelTest.RotateToLabFrameVector2(dirX,dirY,dirZ,trX,trY,trZ);
//    }
//}
//
//BENCHMARK(BM_RotateToLabFrameVector2)->RangeMultiplier(2)->Range(1<<1,1<<3);
//BENCHMARK(BM_RotateToLabFrameVector2)->Repetitions(2);

//// ::: SAMPLESEC
//static void BM_SampleSec_vector(benchmark::State& state) {
//
//    // Set particle kinetic energy
//    static Real        primaryEnergy     = 0.01;             // primary particle energy in [GeV] - 0.01=10MeV
//    EMModelTest emModelTest;
//    for (auto _ : state) {
//        state.PauseTiming();
//        std::vector<LightTrack> lt;
//        emModelTest.GenerateLightTracks(lt, primaryEnergy, state.range(0));
//        state.ResumeTiming();
//        emModelTest.SampleSec_vector(lt);
//        //benchmark::Repetitions(10);
//    }
//}
//
//BENCHMARK(BM_SampleSec_vector)->RangeMultiplier(2)->Range(1<<1,1<<10);
//
//
//static void BM_SampleSec_scalar(benchmark::State& state) {
//
//    // Set particle kinetic energy
//    static Real        primaryEnergy     = 0.01;             // primary particle energy in [GeV] - 0.01=10MeV
//    EMModelTest emModelTest;
//    for (auto _ : state) {
//        state.PauseTiming();
//        std::vector<LightTrack> lt;
//        emModelTest.GenerateLightTracks(lt, primaryEnergy, state.range(0));
//        state.ResumeTiming();
//        emModelTest.SampleSec_scalar(lt);
//    }
//}
////1<<1 equivale a 1 spostato di un bit---> 2
////1<<10 equivale a 1 spostato di 10 bit --> 1*2^10 = 1024
//BENCHMARK(BM_SampleSec_scalar)->RangeMultiplier(2)->Range(1<<1,1<<10);
////BENCHMARK(BM_SampleSec_scalar)->Repetitions(2);



// ::: ALIAS SAMPLING


//static void BM_Alias_scalar(benchmark::State& state) {
//
//    // Set particle kinetic energy
//    static Real        primaryEnergy     = 0.01;             // primary particle energy in [GeV] - 0.01=10MeV
//    EMModelTest emModelTest;
//    for (auto _ : state) {
//        state.PauseTiming();
//        std::vector<Real> kinE, eDirX, eDirY, eDirZ ;
//        for (int i=0; i<state.range(0) ; i++)
//            kinE.push_back(primaryEnergy);
//        state.ResumeTiming();
//        emModelTest.Alias_scalar(kinE, eDirX, eDirY, eDirZ);
//    }
//}
//
//BENCHMARK(BM_Alias_scalar)->RangeMultiplier(2)->Range(1<<1,1<<10);
//
//static void BM_Alias_vector(benchmark::State& state) {
//
//    // Set particle kinetic energy
//    static Real        primaryEnergy     = 0.01;             // primary particle energy in [GeV] - 0.01=10MeV
//    EMModelTest emModelTest;
//    for (auto _ : state) {
//        state.PauseTiming();
//        std::vector<Real> kinE, eDirX, eDirY, eDirZ ;
//        eDirX.resize(state.range(0));
//        eDirY.resize(state.range(0));
//        eDirZ.resize(state.range(0));
//        for (int i=0; i<state.range(0) ; i++)
//            kinE.push_back(primaryEnergy);
//        state.ResumeTiming();
//        emModelTest.Alias_vector(kinE, eDirX, eDirY, eDirZ);
//
//    }
//
//}
//
//BENCHMARK(BM_Alias_vector)->RangeMultiplier(2)->Range(1<<1,1<<10);

//static void BM_Alias_scalar_noRandom(benchmark::State& state) {
//
//    using  Real = EMModelTest::Real;
//    // Set particle kinetic energy
//    static Real        primaryEnergy     = 0.01;             // primary particle energy in [GeV] - 0.01=10MeV
//    EMModelTest emModelTest;
//    for (auto _ : state) {
//        state.PauseTiming();
//
//        //////generate energies vector and random numbers vector
//        geant::TaskData *td = new geant::TaskData(1,1);
//        PhysicsData *phd = new PhysicsData();
//        td->fPhysicsData = phd;
//        double *rndArray = td->fDblArray;
//        std::vector<Real> kinE, eDirX, eDirY, eDirZ, randoms;
//        eDirX.resize(state.range(0));
//        eDirY.resize(state.range(0));
//        eDirZ.resize(state.range(0));
//        for (int i=0; i<state.range(0) ; i++){
//            td->fRndm->uniform_array(1, rndArray);
//            kinE.push_back(0.6+rndArray[0]*primaryEnergy);
//            td->fRndm->uniform_array(4, rndArray);
//            for(int k=0; k<4 ;k++)
//                randoms.push_back(rndArray[k]);
//        }
//        state.ResumeTiming();
//        emModelTest.Alias_scalar_noRandom(kinE, eDirX, eDirY, eDirZ, randoms);
//    }
//}
//
//BENCHMARK(BM_Alias_scalar_noRandom)->Repetitions(5)->RangeMultiplier(2)->Range(1<<1,1<<10);
//
//
//static void BM_Alias_vector_noRandom(benchmark::State& state) {
//
//    using  Real = EMModelTest::Real;
//    // Set particle kinetic energy
//    static Real        primaryEnergy     = 0.01;             // primary particle energy in [GeV] - 0.01=10MeV
//    EMModelTest emModelTest;
//    for (auto _ : state) {
//        state.PauseTiming();
//        geant::TaskData *td = new geant::TaskData(1,1);
//        PhysicsData *phd = new PhysicsData();
//        td->fPhysicsData = phd;
//        double *rndArray = td->fDblArray;
//        std::vector<Real> kinE, eDirX, eDirY, eDirZ, randoms;
//        eDirX.resize(state.range(0));
//        eDirY.resize(state.range(0));
//        eDirZ.resize(state.range(0));
//        for (int i=0; i<state.range(0) ; i++){
//            td->fRndm->uniform_array(1, rndArray);
//            kinE.push_back(0.6+rndArray[0]*primaryEnergy);
//            td->fRndm->uniform_array(4, rndArray);
//            for(int k=0; k<4 ;k++)
//                randoms.push_back(rndArray[k]);
//        }
//        state.ResumeTiming();
//        emModelTest.Alias_vector_noRandom(kinE, eDirX, eDirY, eDirZ, randoms);
//
//    }
//
//}
//
//BENCHMARK(BM_Alias_vector_noRandom)->Repetitions(5)->RangeMultiplier(2)->Range(2048,1<<12);

// **** END: SETTING UP BENCHMARKS FOR ALIAS SAMPLING

static void BM_SampleShell_scalar(benchmark::State& state) {
    
    //std::cout<<"Ok, executing this\n";
    // Set particle kinetic energy
    using  Real = EMModelTest::Real;
    static Real  primaryEnergy= 0.01;             // primary particle energy in [GeV] - 0.01=10MeV
    static int     z            = 82;
    static int     shells      = 24;
    
    EMModelTest emModelTest;
    for (auto _ : state) {
        state.PauseTiming();
        auto td = PrepareTaskData();
        //        geant::TaskData *td = new geant::TaskData(1,1);
        //        //PhysicsData *phd = new PhysicsData();
        //        td->fPhysicsData = phd;
        double *rndArray = td->fDblArray;
        std::vector<Real> kinE, randoms, sampledShells;
        std::vector<size_t> zed, nshells;
        for (int i=0; i<state.range(0) ; i++){
            td->fRndm->uniform_array(2, rndArray);
            kinE.push_back(0.6+rndArray[0]*primaryEnergy);
            zed.push_back(z);
            nshells.push_back(shells);
            randoms.push_back(rndArray[1]);
        }
        sampledShells.resize(state.range(0));
        state.ResumeTiming();
        emModelTest.SampleShell_scalar(kinE, zed, nshells, randoms, sampledShells);
    }
}

BENCHMARK(BM_SampleShell_scalar)->RangeMultiplier(2)->Range(1<<1,1<<5);


//static void BM_SampleShell_vector(benchmark::State& state) {
//
//    using  Real = EMModelTest::Real;
//    // Set particle kinetic energy
//    static Real  primaryEnergy= 0.01;             // primary particle energy in [GeV] - 0.01=10MeV
//    static int     z            = 82;
//    static int     shells      = 24;
//    //std::cout<<"Starting test\n";
//    EMModelTest emModelTest;
//    for (auto _ : state) {
//        state.PauseTiming();
//        auto td = PrepareTaskData();
//        //geant::TaskData *td = new geant::TaskData(1,1);
//        //PhysicsData *phd = new PhysicsData();
//        //td->fPhysicsData = phd;
//        double *rndArray = td->fDblArray;
//        std::vector<Real> kinE, randoms;
//        std::vector<int>  zed, nshells, sampledShells ;
//        for (int i=0; i<state.range(0) ; i++){
//            td->fRndm->uniform_array(2, rndArray);
//            kinE.push_back(0.6+rndArray[0]*primaryEnergy);
//            zed.push_back(z);
//            nshells.push_back(shells);
//            randoms.push_back(rndArray[1]);
//        }
//        sampledShells.resize(state.range(0));
//        state.ResumeTiming();
//        emModelTest.SampleShell_vector(kinE, zed, nshells, randoms, sampledShells);
//    }
//}
//
//BENCHMARK(BM_SampleShell_vector)->RangeMultiplier(2)->Range(10, 10);


void TestCorrectness(int numTracks){
    using  Real = EMModelTest::Real;
    EMModelTest emModelTest;
    std::cerr << "Testing correctness. Vector Size: " << emModelTest.GetVectorSize() << std::endl;
    
    //emModelTest.emModel->funzione(); //in this way I call the SAUTERGAVRILA IMPLEMENTATION
    //emModelTest.funzione(); //in this way I call the emModel implementation
    
    static Real        primaryEnergy     = 0.01;
    std::vector<Real> kinE, eDirX, eDirY, eDirZ ;
    for (int i=0; i<numTracks ; i++)
        kinE.push_back(primaryEnergy);
    emModelTest.Alias_scalar(kinE, eDirX, eDirY, eDirZ);
}

int main(int argc, char** argv) {
    
        TestCorrectness(2);
        std::cout<<"\n Continua da qui, sistemare i benchmark, ma compila\n";
        TestCorrectness(8);
    
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
    ::benchmark::RunSpecifiedBenchmarks();
    
    return 0;
}

