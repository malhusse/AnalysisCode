//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr 28 20:03:07 2019 by ROOT version 6.16/00
// from TChain ntupledData/
//////////////////////////////////////////////////////////

#ifndef plotSelec_h
#define plotSelec_h

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

// Headers needed by this particular selector

class plotSelec : public TSelector
{
public:
   TTreeReader fReader; //!the tree reader
   TTree *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Z window plots
   // 70 - 110 GeV
   TH1 *h_leadMuon_pt_onZ = 0;
   TH1 *h_leadMuon_eta_onZ = 0;
   TH1 *h_subMuon_pt_onZ = 0;
   TH1 *h_subMuon_eta_onZ = 0;
   TH1 *h_dimuon_mass_onZ = 0;
   TH1 *h_dimuon_pt_onZ = 0;
   TH1 *h_dimuon_eta_onZ = 0;
   TH1 *h_dimuon_phi_onZ = 0;
   TH1 *h_dimuon_deta_onZ = 0;
   TH1 *h_dimuon_dphi_onZ = 0;
   TH1 *h_num_jets_onZ = 0;
   TH1 *h_num_bjets_onZ = 0;
   TH1 *h_leadjet_pt_onZ = 0;
   TH1 *h_leadjet_eta_onZ = 0;
   TH1 *h_subjet_pt_onZ = 0;
   TH1 *h_subjet_eta_onZ = 0;
   TH1 *h_dijet_mass_onZ = 0;
   TH1 *h_dijet_deta_onZ = 0;
   TH1 *h_met_pt_onZ = 0;
   TH1 *h_mindrmj_onZ = 0;
   TH1 *h_zeppen_onZ = 0;
   TH1 *h_csTheta_onZ = 0;
   TH1 *h_csPhi_onZ = 0;

   // H window plots
   // 110-150 GeV
   TH1 *h_leadMuon_pt_onH = 0;
   TH1 *h_leadMuon_eta_onH = 0;
   TH1 *h_subMuon_pt_onH = 0;
   TH1 *h_subMuon_eta_onH = 0;
   TH1 *h_dimuon_mass_onH = 0;
   TH1 *h_dimuon_pt_onH = 0;
   TH1 *h_dimuon_eta_onH = 0;
   TH1 *h_dimuon_phi_onH = 0;
   TH1 *h_dimuon_deta_onH = 0;
   TH1 *h_dimuon_dphi_onH = 0;
   TH1 *h_num_jets_onH = 0;
   TH1 *h_num_bjets_onH = 0;
   TH1 *h_leadjet_pt_onH = 0;
   TH1 *h_leadjet_eta_onH = 0;
   TH1 *h_subjet_pt_onH = 0;
   TH1 *h_subjet_eta_onH = 0;
   TH1 *h_dijet_mass_onH = 0;
   TH1 *h_dijet_deta_onH = 0;
   TH1 *h_met_pt_onH = 0;
   TH1 *h_mindrmj_onH = 0;
   TH1 *h_zeppen_onH = 0;
   TH1 *h_csTheta_onH = 0;
   TH1 *h_csPhi_onH = 0;

   // H window only
   TH1 *h_bdtScore01jet = 0;
   TH1 *h_bdtScore2jet = 0;
   TH1 *h_dimuon_0 = 0;
   TH1 *h_dimuon_1 = 0;
   TH1 *h_dimuon_2 = 0;
   TH1 *h_dimuon_3 = 0;
   TH1 *h_dimuon_4 = 0;
   TH1 *h_dimuon_5 = 0;
   TH1 *h_dimuon_6 = 0;
   TH1 *h_dimuon_7 = 0;
   TH1 *h_dimuon_8 = 0;
   TH1 *h_dimuon_9 = 0;
   TH1 *h_dimuon_10 = 0;
   TH1 *h_dimuon_11 = 0;
   TH1 *h_dimuon_12 = 0;
   TH1 *h_dimuon_13 = 0;
   TH1 *h_dimuon_14 = 0;
   // TH1 *h_dimuon_01jet = 0;
   // TH1 *h_dimuon_2jet = 0;
   
   TString _outputName;
   Int_t collectionYear;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Float_t> mclabel = {fReader, "mclabel"};
   // TTreeReaderValue<Float_t> totalSF = {fReader, "totalSF"};
   // TTreeReaderValue<Float_t> eWeight = {fReader, "eWeight"};
   // TTreeReaderValue<Float_t> totalWeight = {fReader, "totalWeight"};
   TTreeReaderValue<Float_t>  puW= {fReader, "puW"};
   TTreeReaderValue<Float_t>  genXsOverN= {fReader, "genXsOverN"};
   TTreeReaderValue<Float_t>  zptW= {fReader, "zptW"};
   TTreeReaderValue<Float_t>  nvtxW= {fReader, "nvtxW"};
   TTreeReaderValue<Float_t>  l1preW = {fReader, "l1preW"};
   TTreeReaderValue<Float_t>  btagSF = {fReader, "btagSF"};
   TTreeReaderValue<Float_t>  idSF = {fReader, "idSF"};
   TTreeReaderValue<Float_t>  isoSF = {fReader, "isoSF"};
   TTreeReaderValue<Float_t>  trigSF = {fReader, "trigSF"};

   TTreeReaderValue<Float_t> muPtC_1 = {fReader, "muPtC_1"};
   TTreeReaderValue<Float_t> muEtaC_1 = {fReader, "muEtaC_1"};
   TTreeReaderValue<Float_t> muPtC_2 = {fReader, "muPtC_2"};
   TTreeReaderValue<Float_t> muEtaC_2 = {fReader, "muEtaC_2"};
   TTreeReaderValue<Float_t> h_mass = {fReader, "h_mass"};
   TTreeReaderValue<Float_t> h_pt = {fReader, "h_pt"};
   TTreeReaderValue<Float_t> h_eta = {fReader, "h_eta"};
   TTreeReaderValue<Float_t> h_phi = {fReader, "h_phi"};
   TTreeReaderValue<Float_t> h_deta = {fReader, "h_deta"};
   TTreeReaderValue<Float_t> h_dphi = {fReader, "h_dphi"};
   TTreeReaderValue<Float_t> njets = {fReader, "njets"};
   TTreeReaderValue<Float_t> nbtagJets = {fReader, "nbtagJets"};
   TTreeReaderValue<Float_t> jetpt_1 = {fReader, "jetpt_1"};
   TTreeReaderValue<Float_t> jeteta_1 = {fReader, "jeteta_1"};
   TTreeReaderValue<Float_t> jetpt_2 = {fReader, "jetpt_2"};
   TTreeReaderValue<Float_t> jeteta_2 = {fReader, "jeteta_2"};
   TTreeReaderValue<Float_t> mjj_1 = {fReader, "mjj_1"};
   TTreeReaderValue<Float_t> detajj_1 = {fReader, "detajj_1"};
   TTreeReaderValue<Float_t> metpt = {fReader, "metpt"};
   TTreeReaderValue<Float_t> mindrmj = {fReader, "mindrmj"};
   TTreeReaderValue<Float_t> zeppen = {fReader, "zeppen"};
   TTreeReaderValue<Float_t> csTheta = {fReader, "csTheta"};
   TTreeReaderValue<Float_t> csPhi = {fReader, "csPhi"};
   TTreeReaderValue<Float_t> category = {fReader, "category"};
   TTreeReaderValue<Float_t> bdtScore = {fReader, "bdtScore"};

   plotSelec(TTree * /*tree*/ = 0) {}
   virtual ~plotSelec() {}
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

   ClassDef(plotSelec, 0);
};

#endif

#ifdef plotSelec_cxx
void plotSelec::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t plotSelec::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef plotSelec_cxx
