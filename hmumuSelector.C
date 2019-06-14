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
double _muonMatchedPt = 30.;
double _muonMatchedEta = 2.4;
double _muonPt = 20.;
double _muonEta = 2.4;
double _muonIso = 0.25;
double _muonMiniIso = .4;

double _electronPt = 20.;
double _electronMiniIso = .4;
double _electronEta = 2.4;

double _dimuonMinMass = 50.;
double _dimuonMaxMass = 200.;
double _JetPt = 30.;
double _JetEta = 4.7;

void hmumuSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   TString option = GetOption();
   std::cout << "TEST" << std::endl;
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

   // _outputNameFinal = "ntupleFiles/";
   // _outputNameFinal += std::to_string(year);
   // _outputNameFinal += "/";
   _outputNameFinal = _outputRoot;

   fProofFile = new TProofOutputFile(_outputNameFinal, TProofOutputFile::kMerge, TProofOutputFile::kLocal);

   fFile = fProofFile->OpenFile("RECREATE");
   if (fFile && fFile->IsZombie())
      SafeDelete(fFile);
   // Cannot continue
   if (!fFile)
   {
      Info("SlaveBegin", "could not create '%s': instance is invalid!", fProofFile->GetName());
      return;
   }

   TParameter<Int_t> *pSumEvents = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getSumEvents"));
   TParameter<Int_t> *pSumEventsWeighted = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getSumEventsWeighted"));

   valueSumEvents = pSumEvents->GetVal();
   valueSumEventsWeighted = pSumEventsWeighted->GetVal();

   TParameter<Double_t> *pxsec = dynamic_cast<TParameter<Double_t> *>(fInput->FindObject("getxsec"));
   TParameter<Int_t> *pmcLabel = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getmcLabel"));
   mcLabel = pmcLabel->GetVal();
   xsec = static_cast<Float_t>(pxsec->GetVal());

   // mcLabel should be 0 for data, so this is true only for mc
   if (mcLabel)
   {
      _dataPUfile = "/uscms_data/d1/malhusse/analysis/AnalysisCode/data/pileup/pu_data_";
      _dataPUfile += std::to_string(year);
      _dataPUfile += ".root";
      // _dataPUfile = "/uscms_data/d1/malhusse/build/AnalysisCode/pileup/pu_data_2017.root";
      _mcPUfile = "/uscms_data/d1/malhusse/analysis/AnalysisCode/data/pileup/";
      _mcPUfile += std::to_string(year);
      _mcPUfile += "/";
      _mcPUfile += _outputRoot;

      weighter = new reweight::LumiReWeighting(_mcPUfile.Data(), _dataPUfile.Data(), "pileup", "pileup");
   }
   if ( ((string) _outputRoot.Data()).find("DY") != string::npos ) //is Drell-Yan Sample
   {
      _dataZptFile = "/uscms_data/d1/malhusse/analysis/AnalysisCode/data/zpt/" ;
      _dataZptFile += std::to_string(year);
      _dataZptFile += "/";
      _dataZptFile += "allData";
      _dataZptFile += std::to_string(year);
      _dataZptFile += ".root";

      _mcZptFile = "/uscms_data/d1/malhusse/analysis/AnalysisCode/data/zpt/";
      _mcZptFile += std::to_string(year);
      _mcZptFile += "/";
      _mcZptFile += _outputRoot;

      zptweighter = new ZptReWeighting(_mcZptFile.Data(), _dataZptFile.Data(), "zpt", "zpt");
   }
   // construct readers
   reader_01jet = new TMVA::Reader();
   reader_2jet = new TMVA::Reader();

   // build string to load ucsd BDT training files
   _01jetxml = "/uscms_data/d1/malhusse/analysis/AnalysisCode/data/Hmm_BDT_xml/";
   _01jetxml += std::to_string(year);
   _01jetxml += "/TMVAClassification_BDTG.weights.01jet.xml";

   _2jetxml = "/uscms_data/d1/malhusse/analysis/AnalysisCode/data/Hmm_BDT_xml/";
   _2jetxml += std::to_string(year);
   _2jetxml += "/TMVAClassification_BDTG.weights.2jet_bveto.xml";

   // // reader variables..
   // float hmmpt, hmmrap, hmmthetacs, hmmphics, j1pt, j1eta, j2pt, detajj, dphijj;
   // float mjj, met, zepen, njets, drmj, m1ptOverMass, m2ptOverMass, m1eta, m2eta;
   // // spectators
   // float hmerr, weight, hmass, nbjets, bdtucsd_inclusive, bdtucsd_01jet, bdtucsd_2jet;

   // Add Variables to readers
   // 2 Jet Reader
   reader_2jet->AddVariable("hmmpt", &hmmpt);
   reader_2jet->AddVariable("hmmrap", &hmmrap);
   reader_2jet->AddVariable("hmmthetacs", &hmmthetacs);
   reader_2jet->AddVariable("hmmphics", &hmmphics);
   reader_2jet->AddVariable("j1pt", &j1pt);
   reader_2jet->AddVariable("j1eta", &j1eta);
   reader_2jet->AddVariable("j2pt", &j2pt);
   reader_2jet->AddVariable("detajj", &detajj);
   reader_2jet->AddVariable("dphijj", &dphijj);
   reader_2jet->AddVariable("mjj", &mjj);
   reader_2jet->AddVariable("met", &met);
   reader_2jet->AddVariable("zepen", &zepen);
   reader_2jet->AddVariable("njets", &njets);
   reader_2jet->AddVariable("drmj", &drmj);
   reader_2jet->AddVariable("m1ptOverMass", &m1ptOverMass);
   reader_2jet->AddVariable("m2ptOverMass", &m2ptOverMass);
   reader_2jet->AddVariable("m1eta", &m1eta);
   reader_2jet->AddVariable("m2eta", &m2eta);
   reader_2jet->AddSpectator("hmerr", &hmerr);
   reader_2jet->AddSpectator("weight", &weight);
   reader_2jet->AddSpectator("hmass", &hmass);
   reader_2jet->AddSpectator("nbjets", &nbjets);
   reader_2jet->AddSpectator("bdtucsd_inclusive", &bdtucsd_inclusive);
   reader_2jet->AddSpectator("bdtucsd_01jet", &bdtucsd_01jet);
   reader_2jet->AddSpectator("bdtucsd_2jet", &bdtucsd_2jet);

   // 01 Jet Reader
   reader_01jet->AddVariable("hmmpt", &hmmpt);
   reader_01jet->AddVariable("hmmrap", &hmmrap);
   reader_01jet->AddVariable("hmmthetacs", &hmmthetacs);
   reader_01jet->AddVariable("hmmphics", &hmmphics);
   reader_01jet->AddVariable("j1pt", &j1pt);
   reader_01jet->AddVariable("j1eta", &j1eta);
   reader_01jet->AddVariable("met", &met);
   reader_01jet->AddVariable("nbjets", &nbjets);
   reader_01jet->AddVariable("drmj", &drmj);
   reader_01jet->AddVariable("njets", &njets);
   reader_01jet->AddVariable("m1ptOverMass", &m1ptOverMass);
   reader_01jet->AddVariable("m2ptOverMass", &m2ptOverMass);
   reader_01jet->AddVariable("m1eta", &m1eta);
   reader_01jet->AddVariable("m2eta", &m2eta);
   reader_01jet->AddSpectator("hmerr", &hmerr);
   reader_01jet->AddSpectator("weight", &weight);
   reader_01jet->AddSpectator("hmass", &hmass);
   reader_01jet->AddSpectator("bdtucsd_inclusive", &bdtucsd_inclusive);
   reader_01jet->AddSpectator("bdtucsd_01jet", &bdtucsd_01jet);
   reader_01jet->AddSpectator("bdtucsd_2jet", &bdtucsd_2jet);

   // Book Reader
   reader_01jet->BookMVA("BDT01jets", _01jetxml);
   reader_2jet->BookMVA("BDT2jets", _2jetxml);

   //   "muRoch_1:muRoch_2:mu_kinfit_pt_1:mu_kinfit_pt_2"
   //   "fsrPt_1:fsrEta_1:fsrPhi_1:fsrPt_2:fsrEta_2:fsrPhi_2"
   //   "muPt_1:muEta_1:muPhi_1:muPt_2:muEta_2:muPhi_2:"

   // string vars = "year:run:lumi:event:mclabel:genWeight:numGen:xsec:puWeight:"
   //               "idSF:isoSF:btagSF:trigSF:prefireSF:totalSF:eWeight:totalWeight:"
   //               "muPtC_1:muEtaC_1:muPhiC_1:"
   //               "muPtC_2:muEtaC_2:muPhiC_2:"
   //               "h_mass:h_pt:h_eta:h_phi:h_deta:h_dphi:"
   //               "njets:ncentJets:nfwdJets:nbtagJets:maxbdisc:"
   //               "jetpt_1:jetmass_1:jeteta_1:jetpt_2:jetmass_2:jeteta_2:"
   //               "mjj_1:detajj_1:mjj_2:detajj_2:"
   //               "metpt:mindrmj:zeppen:csTheta:csPhi:bdtScore";
   //
   // define output TNtuples..
   //
   string vars = "mclabel:totalSF:eWeight:totalWeight:zptWeight:"
                 "muPtC_1:muEtaC_1:"
                 "muPtC_2:muEtaC_2:"
                 "h_mass:h_pt:h_eta:h_phi:h_deta:h_dphi:"
                 "njets:nbtagJets:"
                 "jetpt_1:jeteta_1:jetpt_2:jeteta_2:"
                 "mjj_1:detajj_1:"
                 "metpt:mindrmj:zeppen:csTheta:csPhi:category:bdtScore";

   ntuple = new TNtuple("ntupledData", "Data TNtuple", vars.c_str());

   ntuple->SetDirectory(fFile);
   ntuple->AutoSave();

   GetOutputList()->Add(ntuple);
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

   // create vector of vertices
   // can change this to pass an TTreeReaderArray,
   // dont make vector anymore
   std::vector<analysis::core::Vertex> vecVertices;
   for (int iVer = 0, nVer = Vertices.GetSize(); iVer < nVer; ++iVer)
   {
      vecVertices.push_back(Vertices[iVer]);
   }

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
   // correct Muons using kinfit pt * rochester factor + fsr photon
   // The correct method is to apply rochester w.r.t kinfit values - update this
   for (analysis::core::Muon iMu : Muons)
   {
      if (passMuon(iMu, true))
      {
         TLorentzVector p4mu, p4pho;
         p4mu.SetPtEtaPhiM((iMu._pt_kinfit * iMu._roccCor), iMu._eta, iMu._phi, PDG_MASS_Mu);
         p4pho = iMu.fsrP4;
         p4mu += p4pho;
         iMu._pt = p4mu.Pt();
         iMu._eta = p4mu.Eta();
         iMu._phi = p4mu.Phi();
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
   int nLeptons = nMuons + nElectrons;

   // do jets - cleaned against both muons and electrons!
   int _btagJetsL = 0;
   int _btagJetsM = 0;

   int _ncentJets = 0;
   int _nfwdJets = 0;
   int _numJets = 0;

   float rebtagSF = 1.0;
   float maxBDisc = 0.0;
   float mindrmj = 99;

   for (analysis::core::Jet iJet : Jets)
   {
      if (iJet._pt > _JetPt && TMath::Abs(iJet._eta) < _JetEta && passTightJetID(iJet) && passLoosePUID(iJet) && passNoiseJet(iJet))
      {
         // clean against leptons..
         bool cleanJet = true;
         for (analysis::core::Electron iE : goodElectrons)
         {
            if (jetMuondR(iJet._eta, iJet._phi, iE._eta, iE._phi) < 0.4)
               cleanJet = false;
         }
         if (!cleanJet)
            continue;

         for (analysis::core::Muon iM : goodMuons)
         {
            float drmj = jetMuondR(iJet._eta, iJet._phi, iM._eta, iM._phi);
            if (drmj < 0.4)
               cleanJet = false;
            else
               mindrmj = std::min(mindrmj, drmj);
         }
         if (!cleanJet)
            continue;

         // now the jet is clean and can be used further

         maxBDisc = std::max(maxBDisc, iJet._btag);
         if (iJet._dscvLoose)
            _btagJetsL++;
         if (iJet._dcsvMedium)
            _btagJetsM++;
         _numJets++;
         if (TMath::Abs(iJet._eta) <= 2.4)
            _ncentJets++;
         else
            _nfwdJets++;
         rebtagSF *= iJet._btag_sf ? iJet._btag_sf : 1.0;
         goodJets.push_back(iJet);
      }
   }
   std::sort(goodJets.begin(), goodJets.end(), [](analysis::core::Jet a, analysis::core::Jet b) -> bool { return a._pt > b._pt; });

   // now we have jets muons and electrons..
   // start lookin into muon pairs

   bool goodMuonPair = false;
   float highestPtSum = 0;
   TLorentzVector higgsCandidate;
   analysis::core::Muon leadMuon, subMuon;
   TLorentzVector leadMuonP4, subMuonP4;

   for (unsigned int im = 0; im < goodMuons.size(); ++im)
   {
      for (unsigned int jm = im + 1; jm < goodMuons.size(); ++jm)
      {
         analysis::core::Muon mu1 = goodMuons.at(im);
         analysis::core::Muon mu2 = goodMuons.at(jm);

         if (passMuons(mu1, mu2))
         {
            // muon pair passed..
            TLorentzVector p4m1, p4m2;
            p4m1.SetPtEtaPhiM(mu1._pt, mu1._eta, mu1._phi, PDG_MASS_Mu);
            p4m2.SetPtEtaPhiM(mu2._pt, mu2._eta, mu2._phi, PDG_MASS_Mu);
            TLorentzVector p4dimuon = p4m1 + p4m2;
            if (p4dimuon.Pt() > highestPtSum)
            {
               highestPtSum = p4dimuon.Pt();
               higgsCandidate = p4dimuon;
               // already sorted so (mu1._pt > mu2._pt) should be true
               leadMuon = mu1;
               leadMuonP4 = p4m1;
               subMuon = mu2;
               subMuonP4 = p4m2;
            }
         }
      }
   }

   // remove the leading muons from goodMuons
   goodMuons.erase(std::remove(goodMuons.begin(), goodMuons.end(), leadMuon));
   goodMuons.erase(std::remove(goodMuons.begin(), goodMuons.end(), subMuon));

   // std::cout << "after removal size " << goodMuons.size() << std::endl;

   // sort remaining muons and electrons (jets?)
   //  sort(muonsSelected.begin(), muonsSelected.end(), [](pat::Muon it, pat::Muon jt) -> bool { return it.pt() > jt.pt(); });

   // std::vector<analysis::core::Jet> goodJets;

   // // set to 0 if no good jets in event (thus mindrmj should still be 99)
   if (mindrmj == 99)
      mindrmj = 0;

   TLorentzVector diJet;
   analysis::core::Jet leadJet, subJet;

   if (_numJets == 1)
      leadJet = goodJets.at(0);

   else if (_numJets > 1)
   {
      for (unsigned int i = 0; i < goodJets.size(); ++i)
      {
         for (unsigned int j = i + 1; j < goodJets.size(); ++j)
         {
            analysis::core::Jet j1 = goodJets.at(i);
            analysis::core::Jet j2 = goodJets.at(j);

            TLorentzVector p4lead, p4sub;
            p4lead.SetPtEtaPhiM(j1._pt, j1._eta, j1._phi, j1._mass);
            p4sub.SetPtEtaPhiM(j2._pt, j2._eta, j2._phi, j2._mass);
            TLorentzVector p4dijet = p4lead + p4sub;

            if (p4dijet.M() > diJet.M())
            {
               leadJet = j1;
               subJet = j2;
               diJet = p4dijet;
            }
         }
      }
      // remove the leading jets in case we want 3rd, 4th etc
      goodJets.erase(std::remove(goodJets.begin(), goodJets.end(), leadJet));
      goodJets.erase(std::remove(goodJets.begin(), goodJets.end(), subJet));
   }

   // no need for second jet pair anymore
   // TLorentzVector leadJet2, subJet2, diJet2;
   // // look for next highest mass dijet pair, save as leadJet2, subJet2 and diJet2
   // // First jet pair should be deleted by now!
   // if (p4jets.size() >= 2)
   // {
   //    for (unsigned int i = 0; i < p4jets.size(); ++i)
   //    {
   //       for (unsigned int j = i + 1; j < p4jets.size(); ++j)
   //       {
   //          TLorentzVector p4lead2 = p4jets[i];
   //          TLorentzVector p4sub2 = p4jets[j];
   //          TLorentzVector p4dijet2 = p4lead2 + p4sub2;
   //          if (p4dijet2.M() > diJet2.M())
   //          {
   //             leadJet2 = p4lead2;
   //             subJet2 = p4sub2;
   //             diJet2 = p4dijet2;
   //          }
   //       }
   //    }
   // }
   float h_mass = higgsCandidate.M();
   float h_pt = higgsCandidate.Pt();
   float h_eta = higgsCandidate.Eta();
   float h_phi = higgsCandidate.Phi();
   float h_deta = TMath::Abs(leadMuon._eta - subMuon._eta);
   float h_dphi = TMath::Abs(leadMuon._phi - subMuon._phi);

   // //event weight stuff goes here?
   float pileupWeight = 1.0;
   float totalSF = 1.0;
   float eWeight = 1.0;
   float totalWeight = 1.0;
   float zptWeight = 1.0;

   if (((string) _outputRoot.Data()).find("DY") != string::npos ) //is Drell-Yan Sample
   {
      zptWeight = zptweighter->weight(h_pt);
   }

   if (mcLabel)
   {
      pileupWeight = weighter->weight(*_nPU);
      totalSF = (*_idSF) * (*_isoSF) * (*_trigEffSF) * (*_prefiringweight) * rebtagSF * zptWeight;
      eWeight = pileupWeight * (*_genWeight) * xsec / valueSumEventsWeighted;
      totalWeight = totalSF * eWeight;
   }

   float mupt_1 = leadMuon._pt;
   float mueta_1 = leadMuon._eta;
   // float muphi_1 = leadMuon._phi;
   float mupt_2 = subMuon._pt;
   float mueta_2 = subMuon._eta;
   // float muphi_2 = subMuon._phi;

   float jetpt_1 = leadJet._pt;
   float jetmass_1 = leadJet._mass;
   float jeteta_1 = leadJet._pt ? leadJet._eta : -5;

   float jetpt_2 = subJet._pt;
   float jetmass_2 = subJet._mass;
   float jeteta_2 = subJet._pt ? subJet._eta : -5;

   float mjj_1 = diJet.M() ? diJet.M() : 0;
   float detajj_1 = diJet.M() ? TMath::Abs(leadJet._eta - subJet._eta) : -1;
   float dphijj_1 = diJet.M() ? TMath::Abs(leadJet._phi - subJet._phi) : 0;

   // float mjj_2 = diJet2.M() ? diJet2.M() : 0;
   // float detajj_2 = diJet2.M() ? TMath::Abs(leadJet2.Eta() - subJet2.Eta()) : -1;
   // float dphijj_2 = diJet2.M() ? TMath::Abs(leadJet2.Phi() - subJet2.Phi()) : 0;

   float zeppen = 0;

   if (_numJets > 1 && detajj_1 > 0)
      zeppen = h_eta - (.5 * (jeteta_1 + jeteta_2) / detajj_1);

   float csTheta = getCsTheta(leadMuonP4, subMuonP4);
   float csPhi = getCsPhi(leadMuonP4, subMuonP4);

   float zero = 0.0;
   hmmpt = h_pt;
   hmmrap = h_eta;
   hmmthetacs = csTheta;
   hmmphics = csPhi;
   j1pt = jetpt_1;
   j1eta = std::max(zero, jeteta_1);
   j2pt = jetpt_2;
   detajj = std::max(zero, detajj_1);
   dphijj = dphijj_1;
   mjj = mjj_1;
   met = *_pt;
   zepen = zeppen;
   njets = _numJets;
   nbjets = _btagJetsM;
   drmj = mindrmj;
   m1ptOverMass = (mupt_1 / h_mass);
   m2ptOverMass = (mupt_2 / h_mass);
   m1eta = mueta_1;
   m2eta = mueta_2;

   float bdtScore = -99;
   float category = 0;

   if ((_btagJetsL > 1 || _btagJetsM > 0) && nLeptons > 2)
   {
      category = 2; // ttH Leptonic
   }

   if (!category && _btagJetsM > 1 && _numJets > 5)
   {
      category = 4; // ttH Hadronic Tight
   }
   
   if (!category && _btagJetsM > 0 && _numJets >= 2)
   {
      category = 6; // ttH Hadronic Loose
   }
   
   if (!category && nElectrons > 0)
   {
      category = 8; // VH leptonic electron
   }

   if (!category && nMuons > 2)
   {
      category = 10; // VH leptonic muon
   }

   if (!category && njets == 2 && _nfwdJets == 0 &&
       mjj_1 >= 64 && mjj_1 < 106 &&
       TMath::Abs(higgsCandidate.DeltaR(diJet) - TMath::Pi()) < 0.4 &&
       std::min(leadJet._qgLikelihood, subJet._qgLikelihood) > .25 &&
       std::max(leadJet._qgLikelihood, subJet._qgLikelihood) > .65 &&
       diJet.Pt() > 75)
   {
      category = 12; // VH Hadronic
   }

   if (!category && mjj_1 > 400)
   {
      category = 14; // VBF
   }

   if (!category && njets < 2)
   {
      bdtScore = reader_01jet->EvaluateMVA("BDT01jets");
      category = 16; // 01 jet inclusive
   }

   if (!category && njets >= 2 && nbjets == 0)
   {
      bdtScore = reader_2jet->EvaluateMVA("BDT2jets");
      category = 18; // 02 jet inclusive
   }
   if (!category)
   {
      category = -99; // These are unicorns..
   }

   float toFill[] = {
       static_cast<float>(mcLabel),
       totalSF,
       eWeight,
       totalWeight,
       zptWeight,
       mupt_1,
       mueta_1,
       mupt_2,
       mueta_2,
       h_mass,
       h_pt,
       h_eta,
       h_phi,
       h_deta,
       h_dphi,
       static_cast<float>(_numJets),
       static_cast<float>(_btagJetsM),
       jetpt_1,
       jeteta_1,
       jetpt_2,
       jeteta_2,
       mjj_1,
       detajj_1,
       *_pt,
       mindrmj,
       zeppen,
       csTheta,
       csPhi,
       category,
       bdtScore};

   // float toFill[] = { static_cast<float>(year), static_cast<float>(*_run), static_cast<float>(*_lumi), static_cast<float>(*_event),
   //  static_cast<float>(mcLabel), static_cast<float>(*_genWeight), static_cast<float>(valueSumEventsWeighted), static_cast<float>(xsec),
   //  pileupWeight, *_idSF, *_isoSF, rebtagSF, *_trigEffSF, static_cast<float>(*_prefiringweight), totalSF, eWeight, totalWeight, mupt_1,
   //  mueta_1, muphi_1, mupt_2, mueta_2, muphi_2, h_mass, h_pt, h_eta, h_phi, h_deta, h_dphi, static_cast<float>(_numJets),
   //  static_cast<float>(_ncentJets), static_cast<float>(_nfwdJets), static_cast<float>(_btagJets), maxBDisc, jetpt_1, jetmass_1, jeteta_1,
   //  jetpt_2, jetmass_2, jeteta_2, mjj_1, detajj_1, mjj_2, detajj_2, *_pt, drmj, zepen, hmmthetacs, hmmphics, bdtScore };

   ntuple->Fill(toFill);

   return kTRUE;
}

void hmumuSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   if (fFile)
   {
      if (!ntuple)
      {
         Error("SlaveTerminate", "'ntuple' is undefined!");
         return;
      }
      Bool_t cleanup = kFALSE;
      TDirectory *savedir = gDirectory;
      if (ntuple->GetEntries() > 0)
      {
         fFile->cd();
         ntuple->Write(0, TObject::kOverwrite);
         fProofFile->Print();
         fOutput->Add(fProofFile);
      }
      else
      {
         cleanup = kTRUE;
      }
      ntuple->SetDirectory(0);
      gDirectory = savedir;
      fFile->Close();
      // Cleanup, if needed
      if (cleanup)
      {
         TUrl uf(*(fFile->GetEndpointUrl()));
         SafeDelete(fFile);
         // gSystem->Unlink(uf.GetFile());
         SafeDelete(fProofFile);
      }
   }
}

void hmumuSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   // TFile fout(_outputNameFinal, "recreate");
   // ntuple = dynamic_cast<TNtuple *>(fOutput->FindObject("ntupledData"));
   // if (ntuple)
   // {
   //    ntuple->Write();
   // }

   // fout.Close();
}

// Pass Vertices
bool hmumuSelector::passVertex(std::vector<analysis::core::Vertex> vertexCol)
{
   if (vertexCol.size() == 0)
      return false;

   for (const analysis::core::Vertex &iV : vertexCol)
   {
      if (TMath::Abs(iV._z) < 24 && iV._ndf > 4)
         return true;
   }
   return false;
}

bool hmumuSelector::passElectron(analysis::core::Electron const &e)
{
   if (e._isLoose && e._pt > _electronPt && e._eta < _electronEta &&
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

   if (m._isGlobal && m._isTracker &&
       m._pt > _muonPt && TMath::Abs(m._eta) < _muonEta &&
       m._isMedium && isIsoMuon)
      return true;
   return false;
}

bool hmumuSelector::passMuonHLT(analysis::core::Muon const &m)
{
   if ((m._isHLTMatched[1] || m._isHLTMatched[0]) && m._pt > _muonMatchedPt && TMath::Abs(m._eta) < _muonMatchedEta)
      return true;

   return false;
}

bool hmumuSelector::passMuons(analysis::core::Muon const &m1, analysis::core::Muon const &m2)
{
   if (m1._charge != m2._charge) // Opposite Sign
   {
      if (passMuonHLT(m1) || passMuonHLT(m2)) // at least one hlt triggered muon..
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

bool hmumuSelector::passTightJetID(analysis::core::Jet j)
{
   bool tightID = false;
   double jeta = TMath::Abs(j._eta);
   int numConst = j._cm + j._nm;

   if (jeta <= 2.7)
   {
      tightID = (j._nhf < 0.90 && j._nef < 0.90 && numConst > 1);

      if (jeta < 2.4)
      {
         tightID &= (j._chf > 0 && j._cm > 0);
      }
   }
   else if (jeta <= 3.0)
   {
      tightID = (j._nef > 0.02 && j._nef < 0.99 && j._nm > 2);
   }
   else
   {
      tightID = (j._nef < 0.90 && j._nhf > 0.02 && j._nm > 10);
   }

   return tightID;
}

bool hmumuSelector::passLoosePUID(analysis::core::Jet j)
{
   float jeta = TMath::Abs(j._eta);
   float jpt = j._pt;
   float jpuid = j._puid;

   if (jeta < 2.5)
   {
      if (jpt >= 30 and jpt < 50 and jpuid < -0.89)
         return false;
      if (jpt >= 10 and jpt < 30 and jpuid < -0.97)
         return false;
   }
   else if (jeta < 2.75)
   {
      if (jpt >= 30 and jpt < 50 and jpuid < -0.52)
         return false;
      if (jpt >= 10 and jpt < 30 and jpuid < -0.68)
         return false;
   }
   else if (jeta < 3.0)
   {
      if (jpt >= 30 and jpt < 50 and jpuid < -0.38)
         return false;
      if (jpt >= 10 and jpt < 30 and jpuid < -0.53)
         return false;
   }
   else if (jeta < 5)
   {
      if (jpt >= 30 and jpt < 50 and jpuid < -0.30)
         return false;
      if (jpt >= 10 and jpt < 30 and jpuid < -0.47)
         return false;
   }
   return true;
}

bool hmumuSelector::passNoiseJet(analysis::core::Jet j)
{
   if (!mcLabel && year == 2017)
   {
      float jeta = TMath::Abs(j._eta);
      float jpt = j._pt;
      if (jeta >= 2.65 and jeta <= 3.139 and jpt < 50)
         return false;
   }
   return true;
}

// bool hmumuSelector::passMetFilters(std::vector<std::pair<string, int>> filterBits)
// {
//    bool pass = true;
//    for (const std::pair<string, int> &bit : filterBits)
//    {
//       if (bit.first.compare("Flag_BadChargedCandidateFilter")) // 0 if equal, so gets skipped!
//       {
//          pass = pass && bit.second;
//       }
//    }
//    return pass;
// }

float hmumuSelector::getCsTheta(TLorentzVector v1, TLorentzVector v2)
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

float hmumuSelector::getCsPhi(TLorentzVector v1, TLorentzVector v2)
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
