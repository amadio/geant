#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <dirent.h>
#include "Geant/TNudyENDF.h"
#include "Geant/TNudyEndfSigma.h"
using namespace Nudy;
int main(int, char *argv[]) {
  std::clock_t startTime = clock();
  std::ofstream out, outtotal;
  std::string   str,str2,tmp,tmp2,tmp3,tmp4,tmp5,tmp6, name;
  const char* fENDF = argv[1];
  str = argv[1];
  std::size_t found = str.find_last_of(".");
  name = str.substr(0,found) + ".root";
  const char* rENDF = name.c_str();
  Nudy::TNudyENDF *proc = new Nudy::TNudyENDF(fENDF, rENDF, "recreate");
  proc->SetPreProcess (0) ;
  proc->Process();
  tmp6 = (str.substr(1,found+5));
  tmp5 = "nfy"+tmp6 ;
  std::string fENDFSUB = tmp5;
  proc->SetEndfSub(fENDFSUB);
  proc->Process();
  
  NudyPhysics::TNudyEndfSigma *sigma = new NudyPhysics::TNudyEndfSigma();
  sigma->outstring = str.substr(0,found) +"_300.txt";
  tmp5 = str.substr(0,found)+"_total.txt";
  sigma->outstringTotal.assign(tmp5);
  if(tmp2=="m1")sigma->outstringTotal = tmp5+"m"+"_total.txt";
  sigma->out.open(sigma->outstring.c_str(),std::ios::out);
  if(sigma->out.is_open()){std::cout<<" yes reco open "<< std::endl;}
  sigma->outtotal.open(sigma->outstringTotal.c_str(),std::ios::out);
  if(sigma->outtotal.is_open()){std::cout<<" yes total open "<< std::endl;}
  sigma->SetPreProcess (0) ;
  sigma->SetInitTempDop(0);
  sigma->SetOutTempDop(293.6);
  sigma ->GetData(rENDF, 1E-3);  
  std::clock_t endTime = clock();
  std::clock_t clockTicksTaken = endTime - startTime;
  double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
  std::cout<<"Time in seconds = "<< timeInSeconds <<std::endl;
  sigma->out<<"Time in seconds = "<< timeInSeconds <<std::endl;
  sigma->out.close();
  sigma->outtotal.close();
}
