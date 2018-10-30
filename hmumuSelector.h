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

#include "interface/Muon.h"
#include "interface/Jet.h"
#include "interface/Vertex.h"
#include "interface/Event.h"
#include "interface/MET.h"
#include "interface/MetaHiggs.h"
#include "interface/Electron.h"

class hmumuSelector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   TH1   *h_muon_pt = 0;
   TH1   *h_muon_corrpt = 0;

   TH1   *h_leadMuon_pt = 0;
   TH1   *h_leadMuon_phi = 0;
   TH1   *h_leadMuon_eta = 0;

   TH1   *h_subMuon_pt = 0;
   TH1   *h_subMuon_phi = 0;
   TH1   *h_subMuon_eta = 0;

   TH1   *h_dimuon_mass = 0;
   TH1   *h_dimuon_pt = 0;
   TH1   *h_dimuon_eta = 0;
   TH1   *h_dimuon_phi = 0;
   TH1   *h_dimuon_deta = 0;
   TH1   *h_dimuon_dphi = 0;

   TH1   *h_num_jets = 0;
   TH1   *h_num_bjets = 0;

   TH1   *h_leadjet_pt = 0;
   TH1   *h_leadjet_eta = 0;
   TH1   *h_leadjet_phi = 0;

   TH1   *h_subjet_pt = 0;
   TH1   *h_subjet_eta = 0;
   TH1   *h_subjet_phi = 0;

   TH1   *h_dijet_pt = 0;
   TH1   *h_dijet_mass = 0;
   TH1   *h_dijet_eta = 0;
   TH1   *h_dijet_phi = 0;
   TH1   *h_dijet_deta = 0;
   TH1   *h_dijet_dphi = 0;

   TH1   *h_met_pt = 0;
   TH1   *h_met_eta = 0;
   TH1   *h_met_phi = 0;

   TH1   *h_num_vertices = 0;

   // Readers to access the data (delete the ones you do not need).
   //   TTreeReaderArray<Muon> 
   TTreeReaderArray<analysis::core::Muon> Muons = {fReader,"Muons"};
   TTreeReaderArray<analysis::core::Jet> Jets = {fReader,"Jets"};
   TTreeReaderArray<analysis::core::Vertex> Vertices = {fReader,"Vertices"};
//    TTreeReaderArray<analysis::core::Event> Event = {fReader,"Event"};
//    TTreeReaderArray<analysis::core::EventAuxiliary> EventAux = {fReader,"EventAuxiliary"}
//    TTreeReaderArray<analysis::core::MET> MET = {fReader,"MET"};
//    TTreeReaderArray<analysis::core::MetaHiggs> Meta = {fReader,"Meta"};
   TTreeReaderArray<analysis::core::Electron> Electrons = {fReader,"Electrons"};
//    TTreeReaderArray<analysis::core::
   //   TTreeReaderArray<Int_t> Muons__charge = {fReader, "Muons._charge"};
   //TTreeReaderArray<Float_t> Muons__pt = {fReader, "Muons._pt"};
   // TTreeReaderArray<Float_t> Muons__corrPT = {fReader, "Muons._corrPT"};
   // TTreeReaderArray<Float_t> Muons__pterr = {fReader, "Muons._pterr"};
   //TTreeReaderArray<Float_t> Muons__eta = {fReader, "Muons._eta"};
   //TTreeReaderArray<Float_t> Muons__phi = {fReader, "Muons._phi"};
   //TTreeReaderArray<Bool_t> Muons__isTracker = {fReader, "Muons._isTracker"};
   //TTreeReaderArray<Bool_t> Muons__isGlobal = {fReader, "Muons._isGlobal"};
   // TTreeReaderArray<vector<bool>>
   //TTreeReaderArray<Int_t> Muons__track__charge = {fReader, "Muons._track._charge"};
   //TTreeReaderArray<Float_t> Muons__track__pt = {fReader, "Muons._track._pt"};
   //TTreeReaderArray<Float_t> Muons__track__pterr = {fReader, "Muons._track._pterr"};
   //TTreeReaderArray<Float_t> Muons__track__eta = {fReader, "Muons._track._eta"};
   //TTreeReaderArray<Float_t> Muons__track__phi = {fReader, "Muons._track._phi"};
//    TTreeReaderArray<Float_t> Jets__px = {fReader, "Jets._px"};
//    TTreeReaderArray<Float_t> Jets__py = {fReader, "Jets._py"};
//    TTreeReaderArray<Float_t> Jets__pz = {fReader, "Jets._pz"};
//    TTreeReaderArray<Float_t> Jets__pt = {fReader, "Jets._pt"};
//    TTreeReaderArray<Float_t> Jets__eta = {fReader, "Jets._eta"};
//    TTreeReaderArray<Float_t> Jets__phi = {fReader, "Jets._phi"};
//    TTreeReaderArray<Float_t> Jets__mass = {fReader, "Jets._mass"};
//    TTreeReaderArray<Float_t> Jets__charge = {fReader, "Jets._charge"};
//    TTreeReaderArray<Float_t> Jets__partonFlavour = {fReader, "Jets._partonFlavour"};
//    TTreeReaderArray<Float_t> Jets__chf = {fReader, "Jets._chf"};
//    TTreeReaderArray<Float_t> Jets__nhf = {fReader, "Jets._nhf"};
//    TTreeReaderArray<Float_t> Jets__cef = {fReader, "Jets._cef"};
//    TTreeReaderArray<Float_t> Jets__nef = {fReader, "Jets._nef"};
//    TTreeReaderArray<Float_t> Jets__muf = {fReader, "Jets._muf"};
//    TTreeReaderArray<Float_t> Jets__hfhf = {fReader, "Jets._hfhf"};
//    TTreeReaderArray<Float_t> Jets__hfef = {fReader, "Jets._hfef"};
//    TTreeReaderArray<Float_t> Jets__cm = {fReader, "Jets._cm"};
//    TTreeReaderArray<Float_t> Jets__nm = {fReader, "Jets._nm"};
//    TTreeReaderArray<Float_t> Jets__chm = {fReader, "Jets._chm"};
//    TTreeReaderArray<Float_t> Jets__nhm = {fReader, "Jets._nhm"};
//    TTreeReaderArray<Float_t> Jets__cem = {fReader, "Jets._cem"};
//    TTreeReaderArray<Float_t> Jets__nem = {fReader, "Jets._nem"};
//    TTreeReaderArray<Float_t> Jets__mum = {fReader, "Jets._mum"};
//    TTreeReaderArray<Float_t> Jets__hfhm = {fReader, "Jets._hfhm"};
//    TTreeReaderArray<Float_t> Jets__hfem = {fReader, "Jets._hfem"};
//    TTreeReaderArray<Float_t> Jets__jecf = {fReader, "Jets._jecf"};
//    TTreeReaderArray<Float_t> Jets__jecu = {fReader, "Jets._jecu"};
//    TTreeReaderArray<vector<float>> Jets__btag = {fReader, "Jets._btag"};
//    TTreeReaderArray<Float_t> Jets__puid = {fReader, "Jets._puid"};
//    TTreeReaderArray<Int_t> Jets__fullid = {fReader, "Jets._fullid"};
//    TTreeReaderArray<Double_t> Jets__uncAK5 = {fReader, "Jets._uncAK5"};
//    TTreeReaderArray<Double_t> Jets__uncAK4 = {fReader, "Jets._uncAK4"};
//    TTreeReaderArray<Double_t> Jets__pt_upAK5 = {fReader, "Jets._pt_upAK5"};
//    TTreeReaderArray<Double_t> Jets__pt_upAK4 = {fReader, "Jets._pt_upAK4"};
//    TTreeReaderArray<Double_t> Jets__pt_downAK5 = {fReader, "Jets._pt_downAK5"};
//    TTreeReaderArray<Double_t> Jets__pt_downAK4 = {fReader, "Jets._pt_downAK4"};
//    TTreeReaderArray<Float_t> Jets__genjet__px = {fReader, "Jets._genjet._px"};
//    TTreeReaderArray<Float_t> Jets__genjet__py = {fReader, "Jets._genjet._py"};
//    TTreeReaderArray<Float_t> Jets__genjet__pz = {fReader, "Jets._genjet._pz"};
//    TTreeReaderArray<Float_t> Jets__genjet__pt = {fReader, "Jets._genjet._pt"};
//    TTreeReaderArray<Float_t> Jets__genjet__eta = {fReader, "Jets._genjet._eta"};
//    TTreeReaderArray<Float_t> Jets__genjet__phi = {fReader, "Jets._genjet._phi"};
//    TTreeReaderArray<Float_t> Jets__genjet__mass = {fReader, "Jets._genjet._mass"};
//    TTreeReaderArray<Float_t> Jets__genjet__charge = {fReader, "Jets._genjet._charge"};
//    TTreeReaderArray<Int_t> Vertices__isValid = {fReader, "Vertices._isValid"};
//    TTreeReaderArray<Float_t> Vertices__x = {fReader, "Vertices._x"};
//    TTreeReaderArray<Float_t> Vertices__y = {fReader, "Vertices._y"};
//    TTreeReaderArray<Float_t> Vertices__z = {fReader, "Vertices._z"};
//    TTreeReaderArray<Float_t> Vertices__xerr = {fReader, "Vertices._xerr"};
//    TTreeReaderArray<Float_t> Vertices__yerr = {fReader, "Vertices._yerr"};
//    TTreeReaderArray<Float_t> Vertices__zerr = {fReader, "Vertices._zerr"};
//    TTreeReaderArray<Float_t> Vertices__chi2 = {fReader, "Vertices._chi2"};
//    TTreeReaderArray<Float_t> Vertices__ndf = {fReader, "Vertices._ndf"};
//    TTreeReaderArray<Float_t> Vertices__normChi2 = {fReader, "Vertices._normChi2"};
   TTreeReaderValue<Int_t> _run = {fReader, "_run"};
   TTreeReaderValue<Int_t> _lumi = {fReader, "_lumi"};
   TTreeReaderValue<Long64_t> _event = {fReader, "_event"};
   TTreeReaderValue<Int_t> _bx = {fReader, "_bx"};
   TTreeReaderValue<Int_t> _orbit = {fReader, "_orbit"};
   TTreeReaderValue<Int_t> _nPU = {fReader, "_nPU"};
   TTreeReaderValue<Int_t> _genWeight = {fReader, "_genWeight"};
   TTreeReaderValue<Bool_t> _passedMetFilters = {fReader, "_passedMetFilters"};
   TTreeReaderValue<vector<bool>> _hasHLTFired = {fReader, "_hasHLTFired"};
//    TTreeReaderArray<pair<string,int>> _metFilterBits = {fReader, "_metFilterBits"};
   TTreeReaderValue<Float_t> _px = {fReader, "_px"};
   TTreeReaderValue<Float_t> _py = {fReader, "_py"};
   TTreeReaderValue<Float_t> _pt = {fReader, "_pt"};
   TTreeReaderValue<Float_t> _phi = {fReader, "_phi"};
   TTreeReaderValue<Float_t> _sumEt = {fReader, "_sumEt"};
   TTreeReaderValue<Int_t> _nEventsProcessed = {fReader, "_nEventsProcessed"};
   TTreeReaderValue<Int_t> _sumEventWeights = {fReader, "_sumEventWeights"};
   TTreeReaderValue<Bool_t> _isMC = {fReader, "_isMC"};
   TTreeReaderArray<string> _triggerNames = {fReader, "_triggerNames"};
   TTreeReaderArray<string> _btagNames = {fReader, "_btagNames"};
   TTreeReaderArray<string> _tauIDNames = {fReader, "_tauIDNames"};
   TTreeReaderArray<string> _metFilterNames = {fReader, "_metFilterNames"};
   TTreeReaderValue<unsigned int> _nMuons = {fReader, "_nMuons"};
   TTreeReaderValue<Bool_t> _checkTrigger = {fReader, "_checkTrigger"};
   TTreeReaderValue<Bool_t> _isGlobalMuon = {fReader, "_isGlobalMuon"};
   TTreeReaderValue<Bool_t> _isTrackerMuon = {fReader, "_isTrackerMuon"};
   TTreeReaderValue<Bool_t> _isStandAloneMuon = {fReader, "_isStandAloneMuon"};
   TTreeReaderValue<Float_t> _minPt = {fReader, "_minPt"};
   TTreeReaderValue<Float_t> _maxeta = {fReader, "_maxeta"};
   TTreeReaderValue<Float_t> _maxNormChi2 = {fReader, "_maxNormChi2"};
   TTreeReaderValue<Int_t> _minMuonHits = {fReader, "_minMuonHits"};
   TTreeReaderValue<Int_t> _minPixelHits = {fReader, "_minPixelHits"};
   TTreeReaderValue<Int_t> _minStripHits = {fReader, "_minStripHits"};
   TTreeReaderValue<Int_t> _minTrackerHits = {fReader, "_minTrackerHits"};
   TTreeReaderValue<Int_t> _minSegmentMatches = {fReader, "_minSegmentMatches"};
   TTreeReaderValue<Int_t> _minMatchedStations = {fReader, "_minMatchedStations"};
   TTreeReaderValue<Int_t> _minPixelLayers = {fReader, "_minPixelLayers"};
   TTreeReaderValue<Int_t> _minTrackerLayers = {fReader, "_minTrackerLayers"};
   TTreeReaderValue<Int_t> _minStripLayers = {fReader, "_minStripLayers"};
   TTreeReaderValue<Float_t> _minValidFractionTracker = {fReader, "_minValidFractionTracker"};
   TTreeReaderValue<Float_t> _maxd0 = {fReader, "_maxd0"};
   TTreeReaderValue<Float_t> _maxTrackIsoSumPt = {fReader, "_maxTrackIsoSumPt"};
   TTreeReaderValue<Float_t> _maxRelCombIso = {fReader, "_maxRelCombIso"};
   TTreeReaderValue<Int_t> m_numHT = {fReader, "m_numHT"};
   TTreeReaderValue<Double_t> m_ptSum = {fReader, "m_ptSum"};
//    TTreeReaderArray<Int_t> Electrons__charge = {fReader, "Electrons._charge"};
//    TTreeReaderArray<Float_t> Electrons__pt = {fReader, "Electrons._pt"};
//    TTreeReaderArray<Float_t> Electrons__pterr = {fReader, "Electrons._pterr"};
//    TTreeReaderArray<Float_t> Electrons__eta = {fReader, "Electrons._eta"};
//    TTreeReaderArray<Float_t> Electrons__phi = {fReader, "Electrons._phi"};
//    TTreeReaderArray<vector<bool>> Electrons__ids = {fReader, "Electrons._ids"};
//    TTreeReaderArray<Float_t> Electrons__sumChargedHadronPt = {fReader, "Electrons._sumChargedHadronPt"};
//    TTreeReaderArray<Float_t> Electrons__sumNeutralHadronEt = {fReader, "Electrons._sumNeutralHadronEt"};
//    TTreeReaderArray<Float_t> Electrons__sumPhotonEt = {fReader, "Electrons._sumPhotonEt"};
//    TTreeReaderArray<Float_t> Electrons__sumPUPt = {fReader, "Electrons._sumPUPt"};
//    TTreeReaderArray<Float_t> Electrons__sumChargedParticlePt = {fReader, "Electrons._sumChargedParticlePt"};
//    TTreeReaderArray<Double_t> Electrons__dz = {fReader, "Electrons._dz"};
//    TTreeReaderArray<Bool_t> Electrons__isPF = {fReader, "Electrons._isPF"};
//    TTreeReaderArray<Bool_t> Electrons__convVeto = {fReader, "Electrons._convVeto"};

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
   bool passVertex(TTreeReaderArray<analysis::core::Vertex> vertexCol);
   // bool passMuon();
   // bool passMuonHLT();
   bool passMuons(TTreeReaderArray<analysis::core::Muon> muonCol);
   // float jetMuondR();
   // bool passElectronVeto();
   // bool passBTaggedJetVeto();
   // bool passTightJetID();
   // bool passLoosePUID();

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
