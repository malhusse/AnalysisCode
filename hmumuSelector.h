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
// #include "interface/Electron.h"


class hmumuSelector : public TSelector {
private:
   TNtuple *ntuple = 0;
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   Int_t valueSumEvents = 0;
   Int_t valueSumEventsWeighted = 0;
   Int_t yearM = 0;
   Int_t yearS = 0;
   Int_t mcLabel = -99;
   Double_t xsec = 0;

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
  
   TTreeReaderValue<Int_t> _run = {fReader, "_run"};
   TTreeReaderValue<Int_t> _lumi = {fReader, "_lumi"};
   TTreeReaderValue<Long64_t> _event = {fReader, "_event"};
   TTreeReaderValue<Int_t> _nvtx = {fReader, "_nvtx"};
   TTreeReaderValue<Double_t> _prefiringweight = {fReader, "_prefiringweight"};
   TTreeReaderValue<Float_t> _trigEffSF = {fReader, "_trigEffSF"};
   TTreeReaderValue<Float_t> _idSF = {fReader, "_idSF"};
   TTreeReaderValue<Float_t> _isoSF = {fReader, "_isoSF"};
   TTreeReaderArray<pair<string,int>> _metFilterBits = {fReader, "_metFilterBits"};

   // TTreeReaderValue<Float_t> _bareMCWeight = {fReader, "_bareMCWeight"};
   // TTreeReaderValue<Double_t> _prefiringweightup = {fReader, "_prefiringweightup"};
   // TTreeReaderValue<Double_t> _prefiringweightdown = {fReader, "_prefiringweightdown"};
   // TTreeReaderValue<Float_t> _trigEffSF_up = {fReader, "_trigEffSF_up"};
   // TTreeReaderValue<Float_t> _trigEffSF_down = {fReader, "_trigEffSF_down"};
   // TTreeReaderValue<Float_t> _idSF_up = {fReader, "_idSF_up"};
   // TTreeReaderValue<Float_t> _idSF_down = {fReader, "_idSF_down"};
   // TTreeReaderValue<Float_t> _isoSF_up = {fReader, "_isoSF_up"};
   // TTreeReaderValue<Float_t> _isoSF_down = {fReader, "_isoSF_down"};
   // This is wrong since 0s are not ignored
   // TTreeReaderValue<Float_t> _btagSF = {fReader, "_btagSF"};

   // Should avoid using the Flag_BadChargedCandidateFilter, but is used in 2017
   // so recalculate this...
   // TTreeReaderValue<Bool_t> _passedMetFilters = {fReader, "_passedMetFilters"};

   // Unused Branches -- For Now
   // TTreeReaderArray<analysis::core::Electron> Electrons = {fReader,"Electrons"};
   // TTreeReaderValue<Int_t> _bx = {fReader, "_bx"};
   // TTreeReaderValue<Int_t> _orbit = {fReader, "_orbit"};
   // TTreeReaderArray<Int_t> Muons__charge = {fReader, "Muons._charge"};
   // TTreeReaderArray<Float_t> Muons__pt = {fReader, "Muons._pt"};
   // TTreeReaderArray<Float_t> Muons__pterr = {fReader, "Muons._pterr"};
   // TTreeReaderArray<Float_t> Muons__eta = {fReader, "Muons._eta"};
   // TTreeReaderArray<Float_t> Muons__phi = {fReader, "Muons._phi"};
   // TTreeReaderArray<Bool_t> Muons__isTracker = {fReader, "Muons._isTracker"};
   // TTreeReaderArray<Int_t> Muons__track__charge = {fReader, "Muons._track._charge"};
   // TTreeReaderArray<Float_t> Muons__track__pt = {fReader, "Muons._track._pt"};
   // TTreeReaderArray<Float_t> Muons__track__pterr = {fReader, "Muons._track._pterr"};
   // TTreeReaderArray<Float_t> Muons__track__eta = {fReader, "Muons._track._eta"};
   // TTreeReaderArray<Float_t> Muons__track__phi = {fReader, "Muons._track._phi"};
   // TTreeReaderArray<Float_t> Jets__px = {fReader, "Jets._px"};
   // TTreeReaderArray<Float_t> Jets__py = {fReader, "Jets._py"};
   // TTreeReaderArray<Float_t> Jets__pz = {fReader, "Jets._pz"};
   // TTreeReaderArray<Float_t> Jets__pt = {fReader, "Jets._pt"};
   // TTreeReaderArray<Float_t> Jets__eta = {fReader, "Jets._eta"};
   // TTreeReaderArray<Float_t> Jets__phi = {fReader, "Jets._phi"};
   // TTreeReaderArray<Float_t> Jets__mass = {fReader, "Jets._mass"};
   // TTreeReaderArray<Float_t> Jets__charge = {fReader, "Jets._charge"};
   // TTreeReaderArray<Float_t> Jets__partonFlavour = {fReader, "Jets._partonFlavour"};
   // TTreeReaderArray<Float_t> Jets__chf = {fReader, "Jets._chf"};
   // TTreeReaderArray<Float_t> Jets__nhf = {fReader, "Jets._nhf"};
   // TTreeReaderArray<Float_t> Jets__cef = {fReader, "Jets._cef"};
   // TTreeReaderArray<Float_t> Jets__nef = {fReader, "Jets._nef"};
   // TTreeReaderArray<Float_t> Jets__muf = {fReader, "Jets._muf"};
   // TTreeReaderArray<Float_t> Jets__hfhf = {fReader, "Jets._hfhf"};
   // TTreeReaderArray<Float_t> Jets__hfef = {fReader, "Jets._hfef"};
   // TTreeReaderArray<Float_t> Jets__cm = {fReader, "Jets._cm"};
   // TTreeReaderArray<Float_t> Jets__nm = {fReader, "Jets._nm"};
   // TTreeReaderArray<Float_t> Jets__chm = {fReader, "Jets._chm"};
   // TTreeReaderArray<Float_t> Jets__nhm = {fReader, "Jets._nhm"};
   // TTreeReaderArray<Float_t> Jets__cem = {fReader, "Jets._cem"};
   // TTreeReaderArray<Float_t> Jets__nem = {fReader, "Jets._nem"};
   // TTreeReaderArray<Float_t> Jets__mum = {fReader, "Jets._mum"};
   // TTreeReaderArray<Float_t> Jets__hfhm = {fReader, "Jets._hfhm"};
   // TTreeReaderArray<Float_t> Jets__hfem = {fReader, "Jets._hfem"};
   // TTreeReaderArray<Float_t> Jets__jecf = {fReader, "Jets._jecf"};
   // TTreeReaderArray<Float_t> Jets__jecu = {fReader, "Jets._jecu"};
   // TTreeReaderArray<vector<float>> Jets__btag = {fReader, "Jets._btag"};
   // TTreeReaderArray<Double_t> Jets__btag_sf = {fReader, "Jets._btag_sf"};
   // TTreeReaderArray<Double_t> Jets__btag_sf_up = {fReader, "Jets._btag_sf_up"};
   // TTreeReaderArray<Double_t> Jets__btag_sf_down = {fReader, "Jets._btag_sf_down"};
   // TTreeReaderArray<Float_t> Jets__puid = {fReader, "Jets._puid"};
   // TTreeReaderArray<Int_t> Jets__fullid = {fReader, "Jets._fullid"};
   // TTreeReaderArray<Double_t> Jets__uncAK4 = {fReader, "Jets._uncAK4"};
   // TTreeReaderArray<Double_t> Jets__pt_upAK4 = {fReader, "Jets._pt_upAK4"};
   // TTreeReaderArray<Double_t> Jets__pt_downAK4 = {fReader, "Jets._pt_downAK4"};
   // TTreeReaderArray<Float_t> Jets__jer = {fReader, "Jets._jer"};
   // TTreeReaderArray<Float_t> Jets__jerSF = {fReader, "Jets._jerSF"};
   // TTreeReaderArray<Float_t> Jets__jerSF_up = {fReader, "Jets._jerSF_up"};
   // TTreeReaderArray<Float_t> Jets__jerSF_down = {fReader, "Jets._jerSF_down"};
   // TTreeReaderArray<Float_t> Jets__genjet__px = {fReader, "Jets._genjet._px"};
   // TTreeReaderArray<Float_t> Jets__genjet__py = {fReader, "Jets._genjet._py"};
   // TTreeReaderArray<Float_t> Jets__genjet__pz = {fReader, "Jets._genjet._pz"};
   // TTreeReaderArray<Float_t> Jets__genjet__pt = {fReader, "Jets._genjet._pt"};
   // TTreeReaderArray<Float_t> Jets__genjet__eta = {fReader, "Jets._genjet._eta"};
   // TTreeReaderArray<Float_t> Jets__genjet__phi = {fReader, "Jets._genjet._phi"};
   // TTreeReaderArray<Float_t> Jets__genjet__mass = {fReader, "Jets._genjet._mass"};
   // TTreeReaderArray<Float_t> Jets__genjet__charge = {fReader, "Jets._genjet._charge"};
   // TTreeReaderArray<Int_t> Vertices__isValid = {fReader, "Vertices._isValid"};
   // TTreeReaderArray<Float_t> Vertices__x = {fReader, "Vertices._x"};
   // TTreeReaderArray<Float_t> Vertices__y = {fReader, "Vertices._y"};
   // TTreeReaderArray<Float_t> Vertices__z = {fReader, "Vertices._z"};
   // TTreeReaderArray<Float_t> Vertices__xerr = {fReader, "Vertices._xerr"};
   // TTreeReaderArray<Float_t> Vertices__yerr = {fReader, "Vertices._yerr"};
   // TTreeReaderArray<Float_t> Vertices__zerr = {fReader, "Vertices._zerr"};
   // TTreeReaderArray<Float_t> Vertices__chi2 = {fReader, "Vertices._chi2"};
   // TTreeReaderArray<Float_t> Vertices__ndf = {fReader, "Vertices._ndf"};
   // TTreeReaderArray<Float_t> Vertices__normChi2 = {fReader, "Vertices._normChi2"};
   // TTreeReaderValue<Int_t> _run = {fReader, "_run"};
   // TTreeReaderValue<Int_t> _lumi = {fReader, "_lumi"};
   // TTreeReaderValue<Long64_t> _event = {fReader, "_event"};
   // TTreeReaderValue<Int_t> _bx = {fReader, "_bx"};
   // TTreeReaderValue<Int_t> _orbit = {fReader, "_orbit"};
   


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
   bool passVertex(std::vector<analysis::core::Vertex> vertexCol);
   bool passMuon(analysis::core::Muon const &m);
   bool passMuonHLT(analysis::core::Muon const &m);
   bool passMuons(analysis::core::Muon const &mu1, analysis::core::Muon const &mu2);
   float jetMuondR(float jeta, float jphi, float meta, float mphi);
   bool passTightJetID(analysis::core::Jet j);
   bool passLoosePUID(analysis::core::Jet j);
   bool passNoiseJet(analysis::core::Jet j);
   bool passMetFilters(std::vector<std::pair<string,int>> filterBits);
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
