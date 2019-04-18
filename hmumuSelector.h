//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 24 14:13:01 2018 by ROOT version 6.10/09
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
#include <TH1.h>
#include <TLorentzVector.h>
#include <vector>
#include "TParameter.h"
#include <TNtuple.h>

#include "interface/LumiReweightingStandAlone.h"
#include "interface/Muon.h"
#include "interface/Jet.h"
#include "interface/Vertex.h"
#include "interface/Event.h"
#include "interface/MET.h"
//#include "interface/MetaHiggs.h"
#include "interface/Electron.h"

class hmumuSelector : public TSelector
{
 private:
   TNtuple *ntuple = 0;
 public:
   TTreeReader fReader; //!the tree reader
   TTree *fChain = 0;   //!pointer to the analyzed TTree or TChain

   TH1 *h_eweight = 0;
   TH1 *h_muon_pt = 0;
   TH1 *h_muon_corrpt = 0;

   TH1 *h_leadMuon_pt = 0;
   TH1 *h_leadMuon_phi = 0;
   TH1 *h_leadMuon_eta = 0;

   TH1 *h_subMuon_pt = 0;
   TH1 *h_subMuon_phi = 0;
   TH1 *h_subMuon_eta = 0;

   TH1 *h_dimuon_mass = 0;
   TH1 *h_dimuon_pt = 0;
   TH1 *h_dimuon_eta = 0;
   TH1 *h_dimuon_phi = 0;
   TH1 *h_dimuon_deta = 0;
   TH1 *h_dimuon_dphi = 0;

   TH1 *h_num_jets = 0;
   TH1 *h_num_bjets = 0;

   TH1 *h_leadjet_pt = 0;
   TH1 *h_leadjet_eta = 0;
   TH1 *h_leadjet_phi = 0;

   TH1 *h_subjet_pt = 0;
   TH1 *h_subjet_eta = 0;
   TH1 *h_subjet_phi = 0;

   TH1 *h_dijet_pt = 0;
   TH1 *h_dijet_mass = 0;
   TH1 *h_dijet_eta = 0;
   TH1 *h_dijet_phi = 0;
   TH1 *h_dijet_deta = 0;
   TH1 *h_dijet_dphi = 0;

   TH1 *h_met_pt = 0;

   TH1 *h_num_vertices = 0;

   TH1 *h_numEventsWeighted = 0;
   TH1 *h_numEvents = 0;

   Int_t valueSumEvents = 0;
   Int_t valueSumEventsWeighted = 0;
   Int_t ntupleModeM = -1;
   Int_t ntupleModeS = -1;
   Int_t mcLabel = -99;
   Double_t xsec = 0;


   /* std::vector<TH1F*> vec_dimuon_mass_jets; */
   /* std::vector<TH1F*> vec_dimuon_mass_jets_r; */

   reweight::LumiReWeighting *weighter;
   TString _outputRoot;
   TString _outputNameFinal;
   TString _dataPUfile;
   TString _mcPUfile;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<analysis::core::Muon> Muons = {fReader, "Muons"};
   TTreeReaderArray<analysis::core::Jet> Jets = {fReader, "Jets"};
   TTreeReaderArray<analysis::core::Vertex> Vertices = {fReader, "Vertices"};
   TTreeReaderValue<Int_t> _nPU = {fReader, "_nPU"};
   TTreeReaderValue<Int_t> _genWeight = {fReader, "_genWeight"};
   TTreeReaderValue<Bool_t> _passedMetFilters = {fReader, "_passedMetFilters"};
   TTreeReaderValue<vector<bool>> _hasHLTFired = {fReader, "_hasHLTFired"};

   // This is the MET pt. 
   TTreeReaderValue<Float_t> _pt = {fReader, "_pt"};
  
   // TTreeReaderArray<analysis::core::Electron> Electrons = {fReader,"Electrons"};
   TTreeReaderValue<Int_t> _run = {fReader, "_run"};
   TTreeReaderValue<Int_t> _lumi = {fReader, "_lumi"};
   TTreeReaderValue<Long64_t> _event = {fReader, "_event"};
   // TTreeReaderValue<Int_t> _bx = {fReader, "_bx"};
   // TTreeReaderValue<Int_t> _orbit = {fReader, "_orbit"};

   hmumuSelector(TTree * /*tree*/ = 0) {}
   virtual ~hmumuSelector() {}
   virtual Int_t Version() const { return 2; }
   virtual void Begin(TTree *tree);
   virtual void SlaveBegin(TTree *tree);
   virtual void Init(TTree *tree);
   virtual Bool_t Notify();
   virtual Bool_t Process(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void SetOption(const char *option) { fOption = option; }
   virtual void SetObject(TObject *obj) { fObject = obj; }
   virtual void SetInputList(TList *input) { fInput = input; }
   virtual TList *GetOutputList() const { return fOutput; }
   virtual void SlaveTerminate();
   virtual void Terminate();
   bool passVertex(std::vector<analysis::core::Vertex> vertexCol);
   bool passMuon(analysis::core::Muon const &m);
   bool passMuonHLT(analysis::core::Muon const &m);
   bool passMuons(analysis::core::Muon const &mu1, analysis::core::Muon const &mu2);
   float jetMuondR(float jeta, float jphi, float meta, float mphi);
   bool passTightJetID(analysis::core::Jet j);
   bool passLoosePUID(analysis::core::Jet j);
   bool passNoiseJet(analysis::core::Jet j);

   ClassDef(hmumuSelector, 0);
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
