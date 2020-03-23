#define hmumuSelector_cxx
// The class definition in hmumuSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("hmumuSelector.C")
// root> T->Process("hmumuSelector.C","some options")
// root> T->Process("hmumuSelector.C+")
//

#include "hmumuSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TList.h>
#include <TMath.h>
#include <TString.h>
#include <TNamed.h>

double const PDG_MASS_Mu = 0.1056583745;
double const _muonMatchedEta = 2.4;
double const _muonPt = 20.;
double const _muonEta = 2.4;
double const _muonIso = 0.25;
double const _muonMiniIso = .4;

double const _electronPt = 20.;
double const _electronMiniIso = .4;
double const _electronEta = 2.4;

double const _dimuonMinMass = 50.;
double const _dimuonMaxMass = 200.;
double const _JetPt = 25.;
double const _JetEta = 4.7;

void hmumuSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   TString option = GetOption();
   
   TNamed *name = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));
   _outputRoot = name->GetTitle();
   
   TParameter<Int_t> *pYear = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getYear"));
   year = pYear->GetVal();

   _outputName = "histoFiles/";
   _outputName += std::to_string(year);
   _outputName += "/";
   _outputName += _outputRoot;

}

void hmumuSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   
   TNamed *name = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));
   _outputRoot = name->GetTitle();
   
   TParameter<Int_t> *pYear = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getYear"));
   year = pYear->GetVal();

   _outputName = "histoFiles/";
   _outputName += std::to_string(year);
   _outputName += "/";
   _outputName += _outputRoot;

   // fProofFile = new TProofOutputFile(_outputName, TProofOutputFile::kMerge, TProofOutputFile::kLocal);

   // fFile = fProofFile->OpenFile("RECREATE");
   // if (fFile and fFile->IsZombie())
   //    SafeDelete(fFile);
   // // Cannot continue
   // if (!fFile)
   // {
   //    Info("SlaveBegin", "could not create '%s': instance is invalid!", fProofFile->GetName());
   //    return;
   // }

   TParameter<Int_t> *pSumEvents = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getSumEvents"));
   TParameter<Int_t> *pSumEventsWeighted = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getSumEventsWeighted"));

   valueSumEvents = pSumEvents->GetVal();
   valueSumEventsWeighted = pSumEventsWeighted->GetVal();

   TParameter<Int_t> *pmcLabel = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getmcLabel"));
   mcLabel = pmcLabel->GetVal();


   // mcLabel should be 0 for data, so this is true only for mc
   if (mcLabel)
   {
      _dataPUfile = "/uscms_data/d1/malhusse/analysis/AnalysisCode/resources/data/pileup/pu_data_";
      _dataPUfile += std::to_string(year);
      _dataPUfile += ".root";
      // _dataPUfile = "/uscms_data/d1/malhusse/build/AnalysisCode/pileup/pu_data_2017.root";
      _mcPUfile = "/uscms_data/d1/malhusse/analysis/AnalysisCode/resources/corrections/";
      _mcPUfile += std::to_string(year);
      _mcPUfile += "/";
      _mcPUfile += _outputRoot;

      weighter = new reweight::LumiReWeighting(_mcPUfile.Data(), _dataPUfile.Data(), "pileup", "pileup");
      zptweighter = new zptutils;
   }

   if (year == 2016)
      bound = {-1.00001, -0.125, 0.185, 0.309, .471, 1.00001};
   else if (year == 2017)
      bound = {-1.00001, -0.113, 0.191, 0.315, .498, 1.00001};
   else if (year == 2018)
      bound = {-1.00001, -0.115, 0.191, 0.312, .496, 1.00001};

   hmmpt = new float(0.0);
   hmmrap = new float(0.0);
   hmmthetacs = new float(0.0);
   hmmphics = new float(0.0);
   j1pt = new float(0.0);
   j1eta = new float(0.0);
   j2pt = new float(0.0);
   detajj = new float(0.0);
   dphijj = new float(-1.0);
   mjj = new float(0.0);
   zepen = new float(0.0);
   njets = new float(0.0);
   dphimmj = new float(0.0);
   detammj = new float(0.0);
   m1ptOverMass = new float(0.0);
   m2ptOverMass = new float(0.0);
   m1eta = new float(0.0);
   m2eta = new float(0.0);
   event = new float(1.0);
   hmass = new float(125);
   hmerr = new float(1.0);
   weight = new float(1.0);
   bdtucsd_2jet_nonvbf = new float(0.0);

   // construct readers
   reader_bdt = new TMVA::Reader();

   // build string to load ucsd BDT training files
   // _01jetxml = "/uscms_data/d1/malhusse/analysis/AnalysisCode/resources/Hmm_BDT_xml/";
   // _01jetxml += std::to_string(year);
   // _01jetxml += "/TMVAClassification_BDTG.weights.01jet.xml";

   _bdtxml = "/uscms_data/d1/malhusse/analysis/AnalysisCode/resources/Hmm_BDT_xml/";
   _bdtxml += std::to_string(year);
   _bdtxml += "/TMVAClassification_BDTG.weights.nonvbf.xml";

   _jetUncFile = "/uscms/home/malhusse/nobackup/analysis/AnalysisCode/resources/jetunc/Regrouped_";
   _jetUncFile += std::to_string(year);
   _jetUncFile += "_AK4PFchs.txt";

   // // reader variables..
   // float hmmpt, hmmrap, hmmthetacs, hmmphics, j1pt, j1eta, j2pt, detajj, dphijj;
   // float mjj, met, zepen, njets, drmj, m1ptOverMass, m2ptOverMass, m1eta, m2eta;
   // // spectators
   // float hmerr, weight, hmass, nbjets, bdtucsd_inclusive, bdtucsd_01jet, bdtucsd_2jet;

   // Add Variables to readers
   // 2 Jet Reader
   reader_bdt->AddVariable("hmmpt", hmmpt);
   reader_bdt->AddVariable("hmmrap", hmmrap);
   reader_bdt->AddVariable("hmmthetacs", hmmthetacs);
   reader_bdt->AddVariable("hmmphics", hmmphics);
   reader_bdt->AddVariable("j1pt", j1pt);
   reader_bdt->AddVariable("j1eta", j1eta);
   reader_bdt->AddVariable("j2pt", j2pt);
   reader_bdt->AddVariable("detajj", detajj);
   reader_bdt->AddVariable("dphijj", dphijj);
   reader_bdt->AddVariable("mjj", mjj);
   reader_bdt->AddVariable("zepen", zepen);
   reader_bdt->AddVariable("njets", njets);
   reader_bdt->AddVariable("dphimmj", dphimmj);
   reader_bdt->AddVariable("detammj", detammj);
   reader_bdt->AddVariable("m1ptOverMass", m1ptOverMass);
   reader_bdt->AddVariable("m2ptOverMass", m2ptOverMass);
   reader_bdt->AddVariable("m1eta", m1eta);
   reader_bdt->AddVariable("m2eta", m2eta);

   reader_bdt->AddSpectator("event", event);
   reader_bdt->AddSpectator("hmass", hmass);
   reader_bdt->AddSpectator("hmerr", hmerr);
   reader_bdt->AddSpectator("weight", weight);
   reader_bdt->AddSpectator("bdtucsd_2jet_nonvbf", bdtucsd_2jet_nonvbf);

   // Book Reader
   reader_bdt->BookMVA("BDTreader", _bdtxml);


   // string vars = "year:run:lumi:event:mclabel:eWeight:puWeight:zptWeight:"
   //               "prefireSF:idSF:isoSF:isoUp:isoDown:trigSF:trigUp:trigDown:btagSF:"
   //               "muEtaC_1:muEtaC_2:"
   //               "h_mass:h_pt:h_eta:"
   //               "njets:nbtagJets:"
   //               "jetpt_1:jeteta_1:jetpt_2:"
   //               "mjj:detajj:dphijj:"
   //               "metpt:metphi:zeppen:csTheta:csPhi:category:bdtScore";

   // ntuple = new TNtuple("ntupledData", "Data TNtuple", vars.c_str());

   // ntuple->SetDirectory(fFile);
   // ntuple->AutoSave();
   outputHistograms = new histogramCollection();
   // GetOutputList()->Add(outputHistograms);
}

Bool_t hmumuSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);

   // std::cout << "event number " << *_event << std::endl;
   Long64_t eventDebug = -1;

   // create vector of vertices
   // can change this to pass an TTreeReaderArray,
   // dont make vector anymore
   std::vector<analysis::core::Vertex> vecVertices;
   for (int iVer = 0, nVer = Vertices.GetSize(); iVer < nVer; ++iVer)
   {
      vecVertices.push_back(Vertices[iVer]);
   }

   // make sure to reset the bdt vars so last event info isn't carried over
   *hmmpt = 0.0;
   *hmmrap = 0.0;
   *hmmthetacs = 0.0;
   *hmmphics = 0.0;
   *j1pt = 0.0;
   *j1eta = 0.0;
   *j2pt = 0.0;
   *detajj = 0.0;
   *dphijj = -1.0;
   *mjj = 0.0;
   *zepen = 0.0;
   *njets = 0.0;
   *dphimmj = 0.0;
   *detammj = 0.0;
   *m1ptOverMass = 0.0;
   *m2ptOverMass = 0.0;
   *m1eta = 0.0;
   *m2eta = 0.0;

   // use _passedMetFilters flag..
   // for (int iBit = 0, nBits = _metFilterBits.GetSize(); iBit < nBits; ++iBit)
   // {
   //    metFilterBits.push_back(_metFilterBits[iBit]);
   // }

   // Check for a primary vertex , Trigger and Met Filters
   if (!passVertex(vecVertices))
      return kFALSE;
   if (!(std::any_of(_hasHLTFired->begin(), _hasHLTFired->end(), [](bool v) { return v; })))
      return kFALSE;
   if (!*_passedMetFilters)
      return kFALSE;

   // The following vectors will hold physics objects that
   // pass the required cuts
   std::vector<analysis::core::Muon> goodMuons;
   std::vector<analysis::core::Electron> goodElectrons;
   std::vector<analysis::core::Jet> goodJets;

   // start with muons
   // correct Muons using geofit after rochester + fsr?
   for (analysis::core::Muon iMu : Muons)
   {
      // Try to FSR correct the muon, set pt to roch Cor one or roch Cor + FSR
      if (iMu._corrPT > 20 and passFSR(iMu.fsrP4))
      {
         TLorentzVector p4mu; 
         p4mu.SetPtEtaPhiM(iMu._corrPT, iMu._eta, iMu._phi, PDG_MASS_Mu);
         if ( (iMu.fsrP4.Et() / p4mu.Et()) < 0.4)
         {            
            p4mu = p4mu + iMu.fsrP4;
            iMu._pt = p4mu.Pt();
            iMu._eta = p4mu.Eta();
            iMu._phi = p4mu.Phi();
         }         
      }
      else
      {
            iMu._pt = iMu._geoPT;
      }
      
      if (passMuon(iMu, true)) // corrected should be more than 20?
      {
         goodMuons.push_back(iMu);
      }
   }
   // skip event if didn't find two good muons
   if (goodMuons.size() < 2)
      return kFALSE;

   std::sort(goodMuons.begin(), goodMuons.end(), [](analysis::core::Muon a, analysis::core::Muon b) -> bool { return a._pt > b._pt; });

   int nMuons = goodMuons.size();

   // std::cout << "number of good muons: " << goodMuons.size() << std::endl;
   // pt sort muons?

   // now do electrons?
   for (analysis::core::Electron iEle : Electrons)
   {
      if (passElectron(iEle))
         goodElectrons.push_back(iEle);
   }
   // pt sort Electrons?
   std::sort(goodElectrons.begin(), goodElectrons.end(), [](analysis::core::Electron a, analysis::core::Electron b) -> bool { return a._pt > b._pt; });

   int nElectrons = goodElectrons.size();
   // int nLeptons = nMuons + nElectrons;

   // do jets - cleaned against both muons and electrons!
   int _btagJetsL = 0;
   int _btagJetsM = 0;

   int _ncentJets = 0;
   int _nfwdJets = 0;
   // int _numJets = 0;

   float btagSF = 1.0;
   // float maxBDisc = 0.0;

   for (analysis::core::Jet iJet : Jets)
   {
      if (iJet._pt > _JetPt and TMath::Abs(iJet._eta) < _JetEta and passJetID(iJet, year) and passPUID(iJet, year))
      {
         // clean against leptons..
         bool cleanJet = true;
         for (analysis::core::Muon iM : goodMuons)
         {
            if (jetMuondR(iJet._eta, iJet._phi, iM._eta, iM._phi) < 0.4)
               cleanJet = false;
         }
         if (!cleanJet)
            continue;

         for (analysis::core::Electron iE : goodElectrons)
         {
            if (jetMuondR(iJet._eta, iJet._phi, iE._eta, iE._phi) < 0.4)
               cleanJet = false;
         }
         if (!cleanJet)
            continue;

         // _numJets++;
         
         btagSF *= iJet._btag_sf ? iJet._btag_sf : 1.0;

         // now the jet is clean and can be used further

         // maxBDisc = std::max(maxBDisc, iJet._btag);
         if (iJet._dscvLoose)
            _btagJetsL++;
         if (iJet._dcsvMedium)
            _btagJetsM++;
         if (TMath::Abs(iJet._eta) <= 2.4)
            _ncentJets++;
         else
            _nfwdJets++;
         goodJets.push_back(iJet);
      }
   }
   std::sort(goodJets.begin(), goodJets.end(), [](analysis::core::Jet a, analysis::core::Jet b) -> bool { return a._pt > b._pt; });

   // now we have jets muons and electrons..
   // start looking into muon pair

   // We can start this as it's own function that takes
   // goodJets, goodElectrons, goodMuons


   float pileupWeight = 1.0;
   float zptWeight = 1.0;
   float eWeight = 1.0;

   if (((string)_outputRoot.Data()).find("DY") != string::npos) //is Drell-Yan Sample
   {
      zptWeight = zptweighter->getZPtReweight(*hmmpt, year, std::min(2,  (_ncentJets + _nfwdJets )));
   }

   if (mcLabel)
   {
      // nvtxWeight = (*_nvtx <= 60) ? nvtxFunc->Eval(*_nvtx) : 1.0;
      pileupWeight = weighter->weight(*_nPU);
      eWeight =  (float)*_genWeight / valueSumEventsWeighted;
   }

   float totalWeight = (eWeight) * (pileupWeight) * (zptWeight) * (*_prefiringweight) * (*_idSF) * (*_isoSF) * (*_trigEffSF) * (btagSF);
   float totalWeightIsoUp = (eWeight) * (pileupWeight) * (zptWeight) * (*_prefiringweight) * (*_idSF) * (*_isoSF_up) * (*_trigEffSF) * (btagSF);
   float totalWeightIsoDown = (eWeight) * (pileupWeight) * (zptWeight) * (*_prefiringweight) * (*_idSF) * (*_isoSF_down) * (*_trigEffSF) * (btagSF);

   eventEntry baselineResult = analyze(goodMuons, goodElectrons, goodJets, 
                               nMuons, nElectrons, _ncentJets, _nfwdJets, _btagJetsL, _btagJetsM);

   
   if (!baselineResult.isValid) return kFALSE;

   outputHistograms->fillHistograms(baselineResult, "noSyst", totalWeight);
   outputHistograms->fillHistograms(baselineResult, "isoUp", totalWeightIsoUp);
   outputHistograms->fillHistograms(baselineResult, "isoDown", totalWeightIsoDown);


   return kTRUE;
}

void hmumuSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   // if (outputHistograms)
   // {
      std::cout << "crashing here?" << std::endl;
      std::vector<TH1F*> outputHistos = outputHistograms->getHistograms();

      for (auto h: outputHistos)
      {
         GetOutputList()->Add(h);
      }
   // }
   // if (fFile)
   // {
   //    if (!ntuple)
   //    {
   //       Error("SlaveTerminate", "'ntuple' is undefined!");
   //       return;
   //    }
   //    Bool_t cleanup = kFALSE;
   //    TDirectory *savedir = gDirectory;
   //    if (ntuple->GetEntries() > 0)
   //    {
   //       fFile->cd();
   //       ntuple->Write(0, TObject::kOverwrite);
   //       fProofFile->Print();
   //       fOutput->Add(fProofFile);
   //    }
   //    else
   //    {
   //       cleanup = kTRUE;
   //    }
   //    ntuple->SetDirectory(0);
   //    gDirectory = savedir;
   //    fFile->Close();
   //    // Cleanup, if needed
   //    if (cleanup)
   //    {
   //       TUrl uf(*(fFile->GetEndpointUrl()));
   //       SafeDelete(fFile);
   //       // gSystem->Unlink(uf.GetFile());
   //       SafeDelete(fProofFile);
   //    }
   // }
}

void hmumuSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TFile fout(_outputName, "recreate");
   TList *output_list = (TList *)GetOutputList();
   for (const auto &&obj : *output_list)
   {
      if (obj->IsA()->InheritsFrom("TH1"))
         obj->Write();
   }
   fout.Close();
   // ntuple = dynamic_cast<TNtuple *>(fOutput->FindObject("ntupledData"));
   // if (ntuple)
   // {
   //    ntuple->Write();
   // }

   // fout.Close();
}

// Pass Vertices
bool hmumuSelector::passVertex(std::vector<analysis::core::Vertex>& vertexCol)
{
   if (vertexCol.size() == 0)
      return false;

   for (const analysis::core::Vertex &iV : vertexCol)
   {
      if (TMath::Abs(iV._z) < 24 and iV._ndf > 4)
         return true;
   }
   return false;
}

bool hmumuSelector::passElectron(analysis::core::Electron const &e)
{
   if (e._isLoose and e._pt > _electronPt and e._eta < _electronEta and
       e._miniIso < _electronMiniIso)
      return true;
   return false;
}

bool hmumuSelector::passMuon(analysis::core::Muon const &m, bool useMiniIso)
{
   // This should now exist as m._pfIso
   // double muonIsolation = (m._sumChargedHadronPtR04 + std::max(0., m._sumNeutralHadronEtR04 + m._sumPhotonEtR04 - 0.5 * m._sumPUPtR04)) / m._pt;
   // miniIso is found as m._miniIso .. _muonMiniIso
   bool isIsoMuon = useMiniIso ? m._miniIso < _muonMiniIso : m._pfIso < _muonIso;

   if (m._isGlobal and m._isTracker and
       m._pt > _muonPt and TMath::Abs(m._eta) < _muonEta and
       m._isMedium and isIsoMuon)
      return true;
   return false;
}

bool hmumuSelector::passMuonHLT(analysis::core::Muon const &m, int year)
{
   std::map<int, int> _muonMatchedPt = { {2016, 26}, {2017, 29}, {2018, 26}};
   if ((m._isHLTMatched[1] || m._isHLTMatched[0]) and m._pt > _muonMatchedPt[year] and TMath::Abs(m._eta) < _muonMatchedEta)
      return true;

   return false;
  // bypass for the moment since the trigger paths are not correct?
  return true;
}

bool hmumuSelector::passMuons(analysis::core::Muon const &m1, analysis::core::Muon const &m2, int year)
{
   if (m1._charge != m2._charge) // Opposite Sign
   {
      if (passMuonHLT(m1, year) || passMuonHLT(m2, year)) // at least one hlt triggered muon..
         return true;
   }
   return false;
}

float hmumuSelector::jetMuondR(float jeta, float jphi, float meta, float mphi)
{
   TLorentzVector p4j, p4m;
   p4j.SetPtEtaPhiM(10, jeta, jphi, 0);
   p4m.SetPtEtaPhiM(10, meta, mphi, 0);
   return p4j.DeltaR(p4m);
}

// bool hmumuSelector::passMetFilters(std::vector<std::pair<string, int>> filterBits)
// {
//    bool pass = true;
//    for (const std::pair<string, int> &bit : filterBits)
//    {
//       if (bit.first.compare("Flag_BadChargedCandidateFilter")) // 0 if equal, so gets skipped!
//       {
//          pass = pass and bit.second;
//       }
//    }
//    return pass;
// }

float hmumuSelector::getCsTheta(TLorentzVector& v1, TLorentzVector& v2)
{
   TLorentzVector b1, b2, v;
   b1.SetPx(0);
   b1.SetPy(0);
   b2.SetPx(0);
   b2.SetPy(0);
   double beamE = 6500;    // 1/2 sqrtS in GeV
   double pMass = .938272; // proton mass in GeV
   b1.SetPz(std::hypot(beamE, pMass));
   b1.SetE(beamE);
   b2.SetPz(-std::hypot(beamE, pMass));
   b2.SetE(beamE);

   v = v1 + v2;
   TVector3 boostToVFrame = -v.BoostVector();
   // Boost to higgs frame
   TLorentzVector refV_v1 = v1;
   refV_v1.Boost(boostToVFrame);
   TLorentzVector refV_b1 = b1;
   refV_b1.Boost(boostToVFrame);
   TLorentzVector refV_b2 = b2;
   refV_b2.Boost(boostToVFrame);
   // Getting beam 3-vector from 4-vectors
   TVector3 refV_vb1_direction = refV_b1.Vect().Unit();
   TVector3 refV_vb2_direction = refV_b2.Vect().Unit();
   // Definition of zz directions
   TVector3 direction_cs = (refV_vb1_direction - refV_vb2_direction).Unit(); // CS direction

   return TMath::Cos(direction_cs.Angle(refV_v1.Vect()));
}

float hmumuSelector::getCsPhi(TLorentzVector& v1, TLorentzVector& v2)
{
   TLorentzVector b1, b2, v;
   b1.SetPx(0);
   b1.SetPy(0);
   b2.SetPx(0);
   b2.SetPy(0);
   double beamE = 6500;     // 1/2 sqrtS in GeV
   double pMass = 0.938272; // proton mass in GeV
   b1.SetPz(std::hypot(beamE, pMass));
   b1.SetE(beamE);
   b2.SetPz(-std::hypot(beamE, pMass));
   b2.SetE(beamE);

   v = v1 + v2;
   TVector3 boostToVFrame = -v.BoostVector();
   // Boost to higgs frame
   TLorentzVector refV_v1 = v1;
   refV_v1.Boost(boostToVFrame);
   TLorentzVector refV_b1 = b1;
   refV_b1.Boost(boostToVFrame);
   TLorentzVector refV_b2 = b2;
   refV_b2.Boost(boostToVFrame);
   // Getting beam 3-vector from 4-vectors
   TVector3 refV_vb1_direction = refV_b1.Vect().Unit();
   TVector3 refV_vb2_direction = refV_b2.Vect().Unit();
   // Definition of zz directions
   TVector3 direction_cs = (refV_vb1_direction - refV_vb2_direction).Unit();
   TVector3 xAxis_cs = (refV_vb1_direction + refV_vb2_direction).Unit();
   TVector3 yAxis_cs = direction_cs.Cross(xAxis_cs);

   return std::atan2(refV_v1.Vect() * yAxis_cs, refV_v1.Vect() * xAxis_cs);
}

double hmumuSelector::CosThetaStar(TLorentzVector& v1, TLorentzVector& v2)
{
    TLorentzVector v= v1 + v2;
    const double deta= fabs(v1.Eta() - v2.Eta() );
    const double M=v.M();
    const double Pt = v.Pt();
    return fabs(std::sinh(deta)) / std::sqrt( 1+ std::pow(Pt/M,2)) * 2 * v1.Pt() * v2.Pt() / std::pow(M,2);
}

double hmumuSelector::CosThetaCSPos(analysis::core::Muon& a, analysis::core::Muon& b){ 
   TLorentzVector av, bv;
   av.SetPtEtaPhiM(a._pt, a._eta, a._phi, PDG_MASS_Mu);
   bv.SetPtEtaPhiM(b._pt, b._eta, b._phi, PDG_MASS_Mu);

   if (a._charge < 0) return getCsTheta(av ,bv) ;
   else return getCsTheta(bv, av);
}
double hmumuSelector::PhiCSPos(analysis::core::Muon& a, analysis::core::Muon& b){ 
   TLorentzVector av, bv;
   av.SetPtEtaPhiM(a._pt, a._eta, a._phi, PDG_MASS_Mu);
   bv.SetPtEtaPhiM(b._pt, b._eta, b._phi, PDG_MASS_Mu);

   if (a._charge < 0) return getCsPhi(av ,bv) ; 
   else return getCsPhi(bv, av);
}

bool hmumuSelector::passFSR(TLorentzVector& fsrP4)
{
   if (fsrP4.Pt() > 2)
   {
      if (TMath::Abs(fsrP4.Eta()) < 2.4)
      {
         if (TMath::Abs(fsrP4.Eta()) > 1.4 and TMath::Abs(fsrP4.Eta()) < 1.6)            
            return false;
         return true;
      }
   }
   
   return false;
}


eventEntry hmumuSelector::analyze(std::vector<analysis::core::Muon> goodMuons, std::vector<analysis::core::Electron> goodElectrons, std::vector<analysis::core::Jet> goodJets, 
                                  int numMuons, int numElectrons, int numCentJets, int numFwdJets, int numBtagL, int numBtagM) {

   int _numJets = numCentJets + numFwdJets;
   int nLeptons = numMuons + numElectrons;
   
   // Start of categorization
   analysis::core::Muon *leadMuon = nullptr;
   analysis::core::Muon *subMuon = nullptr;
   for (analysis::core::Muon &im : goodMuons)
   {
      if (leadMuon == nullptr)
         leadMuon = &im;
      if (subMuon == nullptr and (leadMuon->_charge * im._charge == -1))
         subMuon = &im;
   }

   if (leadMuon == nullptr or subMuon == nullptr)
      return eventEntry {false};

   TLorentzVector higgsCandidate;
   TLorentzVector p4m1, p4m2;
   p4m1.SetPtEtaPhiM(leadMuon->_pt, leadMuon->_eta, leadMuon->_phi, PDG_MASS_Mu);
   p4m2.SetPtEtaPhiM(subMuon->_pt, subMuon->_eta, subMuon->_phi, PDG_MASS_Mu);
   higgsCandidate = p4m1 + p4m2;

   TLorentzVector diJet, leadJetV, subJetV;
   analysis::core::Jet leadJet, subJet;

   if (_numJets > 0)
      leadJet = goodJets.at(0);
 
   if (_numJets > 1)
   {
     subJet = goodJets.at(1);

      TLorentzVector p4lead, p4sub, p4dijet;
      
      for (unsigned int i = 0; i < goodJets.size(); ++i)
      {
         for (unsigned int j = i + 1; j < goodJets.size(); ++j)
         {
            analysis::core::Jet j1 = goodJets.at(i);
            analysis::core::Jet j2 = goodJets.at(j);

            p4lead.SetPtEtaPhiM(j1._pt, j1._eta, j1._phi, j1._mass);
            p4sub.SetPtEtaPhiM(j2._pt, j2._eta, j2._phi, j2._mass);
            p4dijet = p4lead + p4sub;

            if ((p4dijet.M() > 400) and (TMath::Abs(j1._eta - j2._eta) > 2.5) and (p4dijet.M() > diJet.M())) // highest mass jet pair > 400 GeV means VBF category.. otherwise use highest pT jets
            {
               leadJet = j1;
               subJet = j2;
               diJet = p4dijet;
            }
         }
      }
   }

   if (_numJets > 0)
      leadJetV.SetPtEtaPhiM(leadJet._pt, leadJet._eta, leadJet._phi, leadJet._mass);
   if (_numJets > 1)
   {      
      subJetV.SetPtEtaPhiM(subJet._pt, subJet._eta, subJet._phi, subJet._mass);
      diJet = leadJetV + subJetV;
   }  
   
   int category = -99;


   if ((category < 0) and (numBtagL > 1 or numBtagM > 0) and nLeptons > 2)
      category = 6; // ttH Leptonic
   
   if ((category < 0) and (numBtagL > 1 or numBtagM > 0) and _numJets > 4)
      category = 5; // ttH Hadronic

   if ((category < 0) and numElectrons > 1) // ZH, Z to electrons
   {
      analysis::core::Electron *e1 = nullptr;
      analysis::core::Electron *e2 = nullptr;

      for (analysis::core::Electron &ie : goodElectrons)
      {
         if (e1 == nullptr)
            e1 = &ie;
         if (e2 == nullptr and (e1->_charge * ie._charge == -1))
            e2 = &ie;
      }

      if (e1 != nullptr and e2 != nullptr)
      {
         TLorentzVector e1LV, e2LV;
         e1LV.SetPtEtaPhiM(e1->_pt, e1->_eta, e1->_phi, 0.511);
         e2LV.SetPtEtaPhiM(e2->_pt, e2->_eta, e2->_phi, 0.511);
         if ((e1LV + e2LV).M() > 71 and (e1LV + e2LV).M() < 111) category = 7;
      }
   }

   

   if ((category < 0) and numMuons > 3) // ZH Z to mumu
   {

      analysis::core::Muon *muZ1 = nullptr;
      analysis::core::Muon *muZ2 = nullptr;

      float zpt = -1;
      for (unsigned i = 0; i < goodMuons.size(); ++i)
      {
         analysis::core::Muon *m1 = &goodMuons[i];
         for (unsigned j = 0; j < goodMuons.size(); ++j)
         {
            analysis::core::Muon *m2 = &goodMuons[j];
            if (m1->_charge * m2->_charge != -1)
               continue;
            TLorentzVector m1LV, m2LV;
            m1LV.SetPtEtaPhiM(m1->_pt, m1->_eta, m1->_phi, PDG_MASS_Mu);
            m2LV.SetPtEtaPhiM(m2->_pt, m2->_eta, m2->_phi, PDG_MASS_Mu);

            if ( (m1LV + m2LV).M() > 71 and (m1LV + m2LV).M() < 111 and (zpt < 0 or (m1LV + m2LV).Pt() < zpt ))
            {
               muZ1 = m1;
               muZ2 = m2;
               zpt = (m1LV + m2LV).Pt();
            }
         }
      }
   
      


      if (zpt > 0)
      {
         analysis::core::Muon *mt1 = nullptr;
         analysis::core::Muon *mt2 = nullptr;

         for (analysis::core::Muon &m : goodMuons)
         {
            
            if (&m == muZ1)
               {
               
                  continue;}
            if (&m == muZ2)
               {
            
                  continue;}
            if (mt1 == nullptr)
               {
                  
                  mt1 = &m;}
            if (mt2 == nullptr and mt1->_charge * m._charge == -1)
               {
                
                  mt2 = &m;}
         }

         if (mt1 != nullptr and mt2 != nullptr)
         {
            
            leadMuon = mt1;
            subMuon = mt2;
            category = 7;
         }
      }
   }

   

   if ((category < 0) and nLeptons >= 3) // WH leptonic
   {
      category = 7;
   }

   

   if ((category < 0) and _numJets >= 2 and numFwdJets == 0 and numBtagM == 0 and
       diJet.M() >= 64 and diJet.M() < 106 and
       TMath::Abs(higgsCandidate.DeltaR(diJet) - TMath::Pi()) < 0.4 and
       std::min(leadJet._qgLikelihood, subJet._qgLikelihood) > .25 and
       std::max(leadJet._qgLikelihood, subJet._qgLikelihood) > .65 and
       diJet.Pt() > 75)
   {
      category = 8; // VH Hadronic
   }


   

   if (category and category != 8) // change the Higgs muons and p4 if exclusive
   {
      higgsCandidate.SetPtEtaPhiM(0, 0, 0, 0);
      TLorentzVector p4mu1, p4mu2;
      p4mu1.SetPtEtaPhiM(leadMuon->_pt, leadMuon->_eta, leadMuon->_phi, PDG_MASS_Mu);
      p4mu2.SetPtEtaPhiM(subMuon->_pt, subMuon->_eta, subMuon->_phi, PDG_MASS_Mu);
      higgsCandidate = p4mu1 + p4mu2;
   }

   
   
   *hmass = higgsCandidate.M();
   *hmmpt = higgsCandidate.Pt();
   *hmmrap = higgsCandidate.Rapidity();
   // float h_phi = higgsCandidate.Phi();
   // float h_deta = TMath::Abs(mu1LV.Eta() - mu2LV.Eta());
   // float h_dphi = TMath::Abs(mu1LV.Eta() - mu2LV.Eta());

   // float m1pt = leadMuon->_pt;
   *m1eta = leadMuon->_eta;
   // float m1phi = leadMuon->_phi;
   // float mupt_2 = subMuon->_pt;
   *m2eta = subMuon->_eta;
   // float muphi_2 = subMuon->_phi;

   if (_numJets > 0)
   {
      *j1pt = leadJet._pt;
      *j1eta = leadJet._eta;
      *dphimmj = TMath::Abs(higgsCandidate.Phi() - leadJet._phi);
      *detammj = TMath::Abs(higgsCandidate.Eta() - leadJet._eta);

      if( _numJets > 1)
      {
         *j2pt = subJet._pt;
         *mjj = diJet.M();
         *detajj = TMath::Abs(leadJet._eta - subJet._eta);
         *dphijj = TMath::Abs(leadJet._phi - subJet._phi);
         *dphimmj = std::min( TMath::Abs(higgsCandidate.Phi() - leadJet._phi) , TMath::Abs(higgsCandidate.Phi() - subJet._phi) );
	      *detammj = std::min( TMath::Abs(higgsCandidate.Eta() - leadJet._eta) , TMath::Abs(higgsCandidate.Eta() - subJet._eta) );
      }

   }

   if (_numJets > 1)
      // *zepen = (higgsCandidate.Rapidity() - .5 * (leadJet._eta + subJet._eta)) / *detajj;
      *zepen = ( higgsCandidate.Rapidity() -  ( leadJetV.Rapidity() + subJetV.Rapidity() ) / 2 ) / ( leadJetV.Rapidity() - subJetV.Rapidity() );


   // 

   // 
   // 



   // TLorentzVector leadMuonP4, subMuonP4;
   // leadMuonP4.SetPtEtaPhiM(hMu1->_pt, hMu1->_eta, hMu1->_phi, PDG_MASS_Mu);
   // subMuonP4.SetPtEtaPhiM(hMu2->_pt, hMu2->_eta, hMu2->_phi, PDG_MASS_Mu);


   *hmmthetacs = CosThetaCSPos(*leadMuon, *subMuon);

   *hmmphics = PhiCSPos(*leadMuon, *subMuon);

   *njets = _numJets;

   *m1ptOverMass = (leadMuon->_pt / higgsCandidate.M());
   *m2ptOverMass = (subMuon->_pt / higgsCandidate.M());

   float bdtScore = -99;
   // nbjets = _btagJetsM;

   std::pair<double,double> met_pt_phi = METXYCorr_Met_MetPhi( *_pt, *_phi, *_run, year, bool(mcLabel), *_nvtx);
   float met_pt = met_pt_phi.first;
   float met_phi = met_pt_phi.second;
   
   if ((category < 0) and *mjj > 400 and *detajj > 2.5 and *j1pt > 35)
   {
      category = 9; // VBF
   }


   if ((category < 0))
   {
      bdtScore = reader_bdt->EvaluateMVA("BDTreader");
      auto ib = std::lower_bound(bound.begin(), bound.end(), bdtScore);
      category = (ib - bound.begin()) - 1;
   }
 
   eventEntry toReturn = {true, leadMuon->_pt, *m1eta, leadMuon->_phi, subMuon->_pt, *m2eta, subMuon->_phi, 
                          *hmass, *hmmpt, static_cast<float>(higgsCandidate.Eta()), static_cast<float>(higgsCandidate.Phi()),
                          static_cast<float>(higgsCandidate.Rapidity()), TMath::Abs(*m1eta - *m2eta),
                          TMath::Abs(leadMuon->_phi - subMuon->_phi), *hmmthetacs,
                          *hmmphics, _numJets, numCentJets, numFwdJets, numBtagL, numBtagM, *j1pt, *j1eta, leadJet._phi,
                          *j2pt, subJet._eta, subJet._phi, *mjj, *detajj, *dphijj, *zepen, *dphimmj, *detammj, met_pt, met_phi,
                          bdtScore, category};
   return toReturn;
}