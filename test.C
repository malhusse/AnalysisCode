#include "interface/NvtxReWeighting.h"
void test(){

  TString _data = "resources/data/allData2018.root";
  TString _mc = "resources/corrections/2018/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8__RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1__MINIAODSIM.root";
  // TString _mc = "GluGluHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8__RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1__MINIAODSIM.root";

  NvtxReWeighting *w = new NvtxReWeighting(_mc.Data(),_data.Data(),"nvtx","nvtx");
  TCanvas *c1 = new TCanvas("c1");
  w->weights_->Draw("");
  c1->Draw();
}
