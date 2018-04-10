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
#include "Geant/SauterGavrilaPhotoElectricModelVec.h"

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

// from geantV
#include "Geant/Typedefs.h"
#include "Geant/TaskData.h"

using namespace geantphysics;

//const int kMinBasket = 16;
const int kMaxBasket = 256;

const double kShellMinPrimEnergy    = 100 * geant::units::eV;
const double kShellMaxPrimEnergy    = 50 *   geant::units::keV;

//const double kShellMinPrimEnergy    = 50 * geant::units::keV;
//const double kShellMaxPrimEnergy    = 10 *   geant::units::MeV;

geant::TaskData *PrepareTaskData()
{
    auto PhysData    = new geantphysics::PhysicsData();
    auto Td          = new geant::TaskData(1, kMaxBasket);
    Td->fPhysicsData = PhysData;
    return Td;
}

class EMModelTest : public geantphysics::EMModel{
public:
    
    using Real_v = vecCore::backend::VcVector::Double_v;
    using Real = double;
    EMModel *emModel;
    EMModel *emModelVec;
    const MaterialCuts *matCut;
    Particle * particle;
    int kBasketSize = 4*10;
    Real kPi = 3.1415;
    
    
    EMModelTest() : geantphysics::EMModel("test") {
        //std::cout<<"Costruttore\n";
        //this->setUpSimulation();
        this->setUpModel();
    }
    ~EMModelTest(){
        delete emModel;
        //delete emModelVec;
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
        //std::cout<<"Initializing the model SauterGavrilaPhotoElectricModel\n";
        emModel->Initialize();
        //===========================================================================================//
    }
    
    int GetVectorSize(){
        return vecCore::VectorSize<Real_v>();
    }
    
    void SampleShellAlias_vector(std::vector<Real>&kinE, std::vector<double>&zed, std::vector<Real>&r1, std::vector<Real>&r2, std::vector<int>&sampledShells){
        this->emModel->SampleShellAliasVec(kinE.data(), zed.data(), r1.data(), r2.data(), sampledShells.data(), kinE.size());
        //for (size_t i=0; i<kinE.size(); i++)
        //    std::cout<<"sampledShells["<<i<<"]: "<<sampledShells[i]<<std::endl;

    }
    
    void SampleShell_scalar(std::vector<Real>& kinE, std::vector<size_t>&zed, std::vector<Real>&randoms, std::vector<Real>&sampledShells){
        
        for(size_t i=0;i<kinE.size();i++){
            //std::cout<<"SampleShell_scalar: "<<i<<std::endl;
            size_t sampledShell;
            this->emModel->SampleShell(kinE[i], zed[i], randoms[i], sampledShell);
            sampledShells.push_back(sampledShell);
        }
    }
    
    void SampleShellAlias_scalar(std::vector<Real>& kinE, std::vector<size_t>&zed, std::vector<Real>&r1,std::vector<Real>&r2, std::vector<Real>&sampledShells){
        
        for(size_t i=0;i<kinE.size();i++){
            //std::cout<<"SampleShell_scalar: "<<i<<std::endl;
            size_t sampledShell;
            this->emModel->SampleShellAlias(kinE[i], zed[i], r1[i],r2[i], sampledShell);
            //std::cout<<"TORNATO: SampleShell_scalar: "<<sampledShell<<std::endl;
            sampledShells.push_back(sampledShell);
        }
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
    
    void GenerateDirectionVectors(std::vector<Real>& x,std::vector<Real>& y,std::vector<Real>& z,int N){
        std::default_random_engine generator;
        std::uniform_real_distribution<Real> distribution(0.0,1.0);
        for(int i = 0; i < N; i++){
            Real phi1 = 2*kPi*distribution(generator);
            Real theta1 = kPi*distribution(generator);
            x.push_back(std::sin(theta1)*std::cos(phi1));
            y.push_back(std::sin(theta1)*std::sin(phi1));
            z.push_back(std::cos(theta1));
        }
    }
    
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

static void BM_Base_SampleShell_scalar(benchmark::State& state) {
    
    using  Real = EMModelTest::Real;
    EMModelTest emModelTest;
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(3,98); // guaranteed unbiased
    
    for (auto _ : state) {
        //state.PauseTiming();
        auto td = PrepareTaskData();
        double *rndArray = td->fDblArray;
        std::vector<Real> kinE, randoms, sampledShells;
        std::vector<size_t> zed;//, nshells;
        for (int i=0; i<state.range(0) ; i++){
            td->fRndm->uniform_array(2, rndArray);
            double energy = kShellMinPrimEnergy+rndArray[0]*(kShellMaxPrimEnergy-kShellMinPrimEnergy);
            kinE.push_back(energy);
            auto random_integer = uni(rng);
            zed.push_back(random_integer);
            randoms.push_back(rndArray[1]);
        }
        //state.ResumeTiming();
        //emModelTest.SampleShell_scalar(kinE, zed, randoms, sampledShells);
        
    }
}

//BENCHMARK(BM_Base_SampleShell_scalar)->RangeMultiplier(2)->Range(16,1<<8);

static void BM_SampleShell_scalar(benchmark::State& state) {

    
    using  Real = EMModelTest::Real;
    EMModelTest emModelTest;
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(3,98); // guaranteed unbiased
    
    for (auto _ : state) {
        state.PauseTiming();
        auto td = PrepareTaskData();
        double *rndArray = td->fDblArray;
        std::vector<Real> kinE, randoms, sampledShells;
        std::vector<size_t> zed;//, nshells;
        for (int i=0; i<state.range(0) ; i++){
            td->fRndm->uniform_array(2, rndArray);
            double energy = kShellMinPrimEnergy+rndArray[0]*(kShellMaxPrimEnergy-kShellMinPrimEnergy);
            kinE.push_back(energy);
            auto random_integer = uni(rng);
            zed.push_back(random_integer);
            randoms.push_back(rndArray[1]);
        }
        state.ResumeTiming();
        emModelTest.SampleShell_scalar(kinE, zed, randoms, sampledShells);

    }
}

//BENCHMARK(BM_SampleShell_scalar)->RangeMultiplier(2)->Range(16,1<<8);


//static void BM_Base_SampleShellAlias(benchmark::State& state) {
//
//    using  Real = EMModelTest::Real;
//    //std::cout<<"BM_SampleShellAlias_scalar\n";
//    EMModelTest emModelTest;
//    //std::cout<<"Model Created\n";
//
//    std::random_device rd;     // only used once to initialise (seed) engine
//    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
//    std::uniform_int_distribution<int> uni(3,98); // guaranteed unbiased
//
//    for (auto _ : state) {
//        //state.PauseTiming();
//        auto td = PrepareTaskData();
//        double *rndArray = td->fDblArray;
//        std::vector<Real> kinE, r1, r2, sampledShells;
//        std::vector<size_t> zed;//, nshells;
//        for (int i=0; i<state.range(0) ; i++){
//            //std::cout<< " ---- "<<i<<std::endl;
//            td->fRndm->uniform_array(3, rndArray);
//            double energy = kShellMinPrimEnergy+rndArray[0]*(kShellMaxPrimEnergy-kShellMinPrimEnergy);
//            kinE.push_back(energy);
//            auto random_integer = uni(rng);
//            zed.push_back(random_integer);
//            r1.push_back(rndArray[1]);
//            r2.push_back(rndArray[2]);
//        }
//        //state.ResumeTiming();
//        //emModelTest.SampleShellAlias_scalar(kinE, zed, r1, r2, sampledShells);
//    }
//}
//
////BENCHMARK(BM_Base_SampleShellAlias)->RangeMultiplier(2)->Range(16,1<<8);

static void BM_SampleShellAlias_scalar(benchmark::State& state) {
    
    
    using  Real = EMModelTest::Real;
    //std::cout<<"BM_SampleShellAlias_scalar\n";
    EMModelTest emModelTest;
    //std::cout<<"Model Created\n";

    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(3,98); // guaranteed unbiased

    for (auto _ : state) {
        state.PauseTiming();
        auto td = PrepareTaskData();
        double *rndArray = td->fDblArray;
        std::vector<Real> kinE, r1, r2, sampledShells;
        std::vector<size_t> zed;//, nshells;
        for (int i=0; i<state.range(0) ; i++){
            //std::cout<< " ---- "<<i<<std::endl;
            td->fRndm->uniform_array(3, rndArray);
            double energy = kShellMinPrimEnergy+rndArray[0]*(kShellMaxPrimEnergy-kShellMinPrimEnergy);
            kinE.push_back(energy);
            auto random_integer = uni(rng);
            zed.push_back(random_integer);
            r1.push_back(rndArray[1]);
            r2.push_back(rndArray[2]);
        }
        state.ResumeTiming();
        emModelTest.SampleShellAlias_scalar(kinE, zed, r1, r2, sampledShells);
    }
}

//BENCHMARK(BM_SampleShellAlias_scalar)->RangeMultiplier(2)->Range(16,1<<8);




//
static void BM_SampleShellAlias_vector(benchmark::State& state) {

    using  Real = EMModelTest::Real;
    EMModelTest emModelTest;
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(3,98); // guaranteed unbiased
    
    for (auto _ : state) {
        state.PauseTiming();
        auto td = PrepareTaskData();
        double *rndArray = td->fDblArray;
        std::vector<Real> zed, kinE, r1, r2;
        std::vector<int>  nshells, sampledShells ;
        for (int i=0; i<state.range(0) ; i++){
            td->fRndm->uniform_array(3, rndArray);
            double energy = kShellMinPrimEnergy+rndArray[0]*(kShellMaxPrimEnergy-kShellMinPrimEnergy);
            kinE.push_back(energy);
            auto random_integer = uni(rng);
            zed.push_back(random_integer);
            r1.push_back(rndArray[1]);
            r2.push_back(rndArray[2]);
        }
        sampledShells.resize(state.range(0));
        state.ResumeTiming();
        emModelTest.SampleShellAlias_vector(kinE, zed, r1, r2, sampledShells);
    }
}

BENCHMARK(BM_SampleShellAlias_vector)->RangeMultiplier(2)->Range(16,1<<8);



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
    //emModelTest.Alias_scalar(kinE, eDirX, eDirY, eDirZ);
}

int main(int argc, char** argv) {
    
        //TestCorrectness(2);
        //TestCorrectness(8);
    
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
    ::benchmark::RunSpecifiedBenchmarks();
    
    return 0;
}

