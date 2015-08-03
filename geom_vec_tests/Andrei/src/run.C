{
  TString arch = gSystem->GetBuildArch();
  if (arch.Contains("icc")) {
//    gSystem->SetFlagsOpt("-O0");
    gSystem->SetFlagsOpt("-O3 -vec-report2");
  } else {
    gSystem->SetFlagsOpt("-O3 -march=amdfam10 -msse3 -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=2");
  }  
  gROOT->LoadMacro("vectors.C++O");
  vectors(1000000);
}
  
