//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 18 13:08:17 2019 by ROOT version 6.10/09
// from TChain ntuplemaker_H2DiMuonMaker/Events/
//////////////////////////////////////////////////////////

#ifndef hmumuSelector_h
#define hmumuSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <random>
// #include <TH1.h>
#include <TLorentzVector.h>
#include <vector>
#include "TParameter.h"
// #include <TNtuple.h>
// #include <TProofOutputFile.h>
#include <TMVA/Reader.h>
// #include <TF1.h>

#include "interface/LumiReweightingStandAlone.h"
#include "interface/ZptReweight.h"
// #include "interface/NvtxReWeighting.h"
#include "interface/MetCorrection.h"
#include "interface/jetIDs.h"
#include "interface/Muon.h"
#include "interface/Jet.h"
#include "interface/Vertex.h"
#include "interface/Event.h"
#include "interface/MET.h"
#include "interface/MetaHiggs.h"
#include "interface/Electron.h"
#include "interface/histogramCollection.h"
#include "interface/jetUncertainties.h"
// #include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
// #include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"


class hmumuSelector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   Int_t valueSumEvents = 0;
   Int_t valueSumEventsWeighted = 0;
   Int_t year = 0;
   Int_t mcLabel = -99;
  
   // const TF1 *nvtxFunc = 0;

   reweight::LumiReWeighting *weighter;
   reweight::LumiReWeighting *weighter_up;
   reweight::LumiReWeighting *weighter_down;


   zptutils *zptweighter;

   TString _dataCorrFile;
   TString _mcCorrFile;

   TString _outputRoot;
   TString _outputName;
   
   TString _dataPUfile;
   TString _dataPUfile_up;
   TString _dataPUfile_down;

   TString _mcPUfile;
   TString _bdtxml;
   
   // TString _jetUncFile;
   
   TMVA::Reader* reader_bdt = 0;
   // TMVA::Reader* reader_2jet = 0;

   // make these pointers and initialize them in
   // SlaveBegin() then use them to feed the BDT
   // reader variables..
   float *hmmpt, *hmmrap, *hmmthetacs, *hmmphics, *j1pt, *j1eta, *j2pt, *detajj, *dphijj;
   float *mjj, *zepen, *njets, *dphimmj, *detammj, *m1ptOverMass, *m2ptOverMass, *m1eta, *m2eta; 
   // spectators
   float *event, *hmass, *hmerr, *weight, *bdtucsd_2jet_nonvbf;

   vector<float> bound;

   std::mt19937 m_random_generator;
   // std::vector<std::string> jetUncertainties_;
   std::map<std::string, JetCorrectionUncertainty*> jetUncertainties_;
   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<analysis::core::Muon> Muons = {fReader, "Muons"};
   TTreeReaderArray<analysis::core::Jet> Jets = {fReader, "Jets"};
   TTreeReaderArray<analysis::core::Vertex> Vertices = {fReader, "Vertices"};
   TTreeReaderArray<analysis::core::Electron> Electrons = {fReader,"Electrons"};

   TTreeReaderValue<Int_t> _nPU = {fReader, "_nPU"};
   TTreeReaderValue<Int_t> _genWeight = {fReader, "_genWeight"};
   TTreeReaderValue<vector<bool>> _hasHLTFired = {fReader, "_hasHLTFired"};

   // This is the MET pt, phi, px and py variables
   TTreeReaderValue<Float_t> _pt = {fReader, "_pt"};
   TTreeReaderValue<Float_t> _phi = {fReader, "_phi"};
   // TTreeReaderValue<Float_t> _px = {fReader, "_px"};
   // TTreeReaderValue<Float_t> _py = {fReader, "_py"};

  
   TTreeReaderValue<Int_t> _run = {fReader, "_run"};
   TTreeReaderValue<Int_t> _lumi = {fReader, "_lumi"};
   TTreeReaderValue<Long64_t> _event = {fReader, "_event"};
   TTreeReaderValue<Int_t> _nvtx = {fReader, "_nvtx"};

   TTreeReaderValue<Double_t> _prefiringweight = {fReader, "_prefiringweight"};
   TTreeReaderValue<Double_t> _prefiringweightup = {fReader, "_prefiringweightup"};
   TTreeReaderValue<Double_t> _prefiringweightdown = {fReader, "_prefiringweightdown"};

   TTreeReaderValue<Float_t> _trigEffSF = {fReader, "_trigEffSF"};
   TTreeReaderValue<Float_t> _trigEffSF_up = {fReader, "_trigEffSF_up"};
   TTreeReaderValue<Float_t> _trigEffSF_down = {fReader, "_trigEffSF_down"};

   TTreeReaderValue<Float_t> _idSF = {fReader, "_idSF"};
   TTreeReaderValue<Float_t> _idSF_up = {fReader, "_idSF_up"};
   TTreeReaderValue<Float_t> _idSF_down = {fReader, "_idSF_down"};

   TTreeReaderValue<Float_t> _isoSF = {fReader, "_isoSF"};
   TTreeReaderValue<Float_t> _isoSF_up = {fReader, "_isoSF_up"};
   TTreeReaderValue<Float_t> _isoSF_down = {fReader, "_isoSF_down"};

   // TTreeReaderValue<Float_t> _btagSF = {fReader, "_btagSF"};

   // This is now fixed, so we can just use it instead of the individual filters..
   TTreeReaderValue<Bool_t> _passedMetFilters = {fReader, "_passedMetFilters"};
   // TTreeReaderArray<pair<string,int>> _metFilterBits = {fReader, "_metFilterBits"};


   histogramCollection* outputHistograms = 0;

   hmumuSelector(TTree * /*tree*/ =0) { }
   virtual ~hmumuSelector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   bool passVertex(std::vector<analysis::core::Vertex>& vertexCol);
   bool passMuon(analysis::core::Muon const &m, bool useMiniIso);
   bool passElectron(analysis::core::Electron const &e);
   bool passMuonHLT(analysis::core::Muon const &m, int year);
   bool passMuons(analysis::core::Muon const &mu1, analysis::core::Muon const &mu2, int year);
   float jetMuondR(float jeta, float jphi, float meta, float mphi);
   // bool passTightJetID(analysis::core::Jet j);
   // bool passLoosePUID(analysis::core::Jet j);
   // bool passNoiseJet(analysis::core::Jet j);
   bool passFSR(TLorentzVector& fsrP4);
   // should already be correct using _passedMetFilters
   // bool passMetFilters(std::vector<std::pair<string,int>> filterBits);
   float deltaPhi(float phiOne, float phiTwo);
   float getCsTheta(TLorentzVector& v1, TLorentzVector& v2);
   float getCsPhi(TLorentzVector& v1, TLorentzVector& v2);
   double CosThetaStar(TLorentzVector& v1, TLorentzVector& v2);
   double CosThetaCSPos(analysis::core::Muon& a, analysis::core::Muon& b);
   double PhiCSPos(analysis::core::Muon& a, analysis::core::Muon& b);
   std::vector<analysis::core::Jet> shiftJets(std::vector<analysis::core::Jet> inputJets, std::string uncSource, bool shiftUp);
   std::vector<analysis::core::Jet> smearJets(std::vector<analysis::core::Jet> inputJets);

   eventEntry analyze(analysis::core::Muon* leadMuon, analysis::core::Muon* subMuon, std::vector<analysis::core::Muon> goodMuons, std::vector<analysis::core::Electron> goodElectrons, std::vector<analysis::core::Jet> preselectedJets);
   ClassDef(hmumuSelector,0);

};

#endif

#ifdef hmumuSelector_cxx
void hmumuSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t hmumuSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef hmumuSelector_cxx
