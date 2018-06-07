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
using namespace Nudy;
int main(int, char *argv[]) {
  std::clock_t startTime = clock();
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
  std::clock_t endTime = clock();
  std::clock_t clockTicksTaken = endTime - startTime;
  double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
  std::cout<<"Time in seconds = "<< timeInSeconds <<std::endl;
}
