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
double _dimuonMinMass = 100.;
double _dimuonMaxMass = 200.;
double _JetPt = 30.;
double _JetEta = 4.7;

void hmumuSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   TNamed *name = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));

   TParameter<Int_t> *pNtupleMode = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getNtupleMode"));
   ntupleModeM = pNtupleMode->GetVal();

   TParameter<Int_t> *pSumEvents = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getSumEvents"));
   TParameter<Int_t> *pSumEventsWeighted = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getSumEventsWeighted"));
   valueSumEvents = pSumEvents->GetVal();
   valueSumEventsWeighted = pSumEventsWeighted->GetVal();

   if (ntupleModeM)
   {
      _outputNameFinal = "ntupleFiles/";
   }
   else
   {
      _outputNameFinal = "histoFiles/";
      h_numEventsWeighted = new TH1I("numEventsWeighted", "Weighted numEvents Proccessed", 2, 0, 2);
      h_numEvents = new TH1I("numEvents", "numEvents Processed", 2, 0, 2);
      GetOutputList()->Add(h_numEventsWeighted);
      GetOutputList()->Add(h_numEvents);
      h_numEventsWeighted->Fill(1, valueSumEventsWeighted);
      h_numEvents->Fill(1, valueSumEvents);
   }
   _outputNameFinal += name ? name->GetTitle() : "outputName.root";
}

void hmumuSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   TNamed *name2 = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));
   _outputRoot = name2->GetTitle();

   TParameter<Int_t> *pSumEventsWeightedAgain = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getSumEventsWeighted"));
   valueSumEventsWeighted = pSumEventsWeightedAgain->GetVal();

   TParameter<Int_t> *pNtupleModeAgain = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getNtupleMode"));
   TParameter<Double_t> *pxsec = dynamic_cast<TParameter<Double_t> *>(fInput->FindObject("getxsec"));
   TParameter<Int_t> *pmcLabel = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getmcLabel"));

   // ntupleMode = 0 or 1
   // let 0 be histogram creation
   // and 1 be TNtuple creation
   ntupleModeS = pNtupleModeAgain->GetVal();
   mcLabel = pmcLabel->GetVal();
   xsec = static_cast<Float_t>(pxsec->GetVal());

   // _isMC = true;
   // if (mcLabel)
   // _isMC = false;

   // mcLabel should be 0 for data, so this is true only for mc!
   if (mcLabel)
   {
      _dataPUfile = "/uscms_data/d1/malhusse/build/AnalysisCode/pileup/pu_data_2017.root";
      _mcPUfile = "/uscms_data/d1/malhusse/build/AnalysisCode/pileup/";
      _mcPUfile += _outputRoot;

      weighter = new reweight::LumiReWeighting(_mcPUfile.Data(), _dataPUfile.Data(), "pileup", "pileup");
   }

   TString option = GetOption();

   if (ntupleModeS)
   {
      // TODO: split this over multiple lines (cleaner code)
      string vars = "run:lumi:event:mclabel:eweight:h_mass:h_pt:h_eta:h_phi:h_deta:h_dphi:mupt_1:mueta_1:muphi_1:mupt_2:mueta_2:muphi_2:njets:ncentJets:nfwdJets:nbtagJets:jetpt_1:jetmass_1:jeteta_1:jetpt_2:jetmass_2:jeteta_2:mjj_1:detajj_1:mjj_2:detajj_2:metpt";
      ntuple = new TNtuple("ntupledData", "Data TNtuple", vars.c_str());
      GetOutputList()->Add(ntuple);
   }
   else
   {
      h_muon_pt = new TH1F("muon_pt", "Muon p_{T};p_{T}  (GeV);Events ", 200, 0, 400);
      h_muon_corrpt = new TH1F("muon_corrpt", "Muon corrected p_{T};p_{T}  (GeV);Events ", 200, 0, 400);
      h_leadMuon_pt = new TH1F("lead_muon_pt", "Leading Muon p_{T};p_{T}  (GeV);Events ", 100, 0, 200);
      h_leadMuon_eta = new TH1F("lead_muon_eta", "Leading Muon \\eta;\\eta;Events ", 50, -2.5, 2.5);
      h_leadMuon_phi = new TH1F("lead_muon_phi", "Leading Muon \\phi;\\phi;Events ", 36, -3.6, 3.6);
      h_subMuon_pt = new TH1F("sub_muon_pt", "Subleading Muon p_{T};p_{T}  (GeV);Events ", 100, 0, 200);
      h_subMuon_eta = new TH1F("sub_muon_eta", "Subleading Muon eta;\\eta;Events ", 50, -2.5, 2.5);
      h_subMuon_phi = new TH1F("sub_muon_phi", "Subleading Muon \\phi;\\phi;Events ", 36, -3.6, 3.6);
      h_dimuon_mass = new TH1F("dimuon_mass", "Dimuon Mass;M_{\\mu \\mu}  (Gev);Events ", 100, 100, 200);
      h_dimuon_pt = new TH1F("dimuon_pt", "Dimuon p_{T};p_{T}  (GeV);Events ", 200, 0, 400);
      h_dimuon_eta = new TH1F("dimuon_eta", "Dimuon \\eta;\\eta;Events ", 100, -5.0, 5.0);
      h_dimuon_phi = new TH1F("dimuon_phi", "Dimuon \\phi;\\phi;Events ", 36, -3.6, 3.6);
      h_dimuon_deta = new TH1F("dimuon_deta", "Dimuon deta;deta;Events ", 100, -5.0, 5.0);
      h_dimuon_dphi = new TH1F("dimuon_dphi", "Dimuon dphi;dphi;Events ", 36, -3.6, 3.6);
      h_num_jets = new TH1F("num_jets", "Number of Jets;nJets;Events ", 10, 0, 10);
      h_num_bjets = new TH1F("num_bjets", "Number of B Jets;nBJets;Events ", 10, 0, 10);
      h_leadjet_pt = new TH1F("leadjet_pt", "Leading Jet p_{T};p_{T}  (GeV);Events ", 500, 0, 500);
      h_leadjet_eta = new TH1F("leadjet_eta", "Leading Jet \\eta;\\eta;Events ", 94, -4.7, 4.7);
      h_leadjet_phi = new TH1F("leadjet_phi", "Leading Jet \\phi;\\phi;Events ", 36, -3.6, 3.6);
      h_subjet_pt = new TH1F("subjet_pt", "Subleading Jet p_{T};p_{T}  (GeV);Events ", 500, 0, 500);
      h_subjet_eta = new TH1F("subjet_eta", "Subleading Jet \\eta;\\eta;Events ", 94, -4.7, 4.7);
      h_subjet_phi = new TH1F("subjet_phi", "Subleading Jet \\phi;\\phi;Events ", 36, -3.6, 3.6);
      h_dijet_pt = new TH1F("dijet_pt", "DiJet p_{T};p_{T}  (GeV),Events ", 1000, 0, 1000);
      h_dijet_mass = new TH1F("dijet_mass", "DiJet Mass;M_{jj}  (GeV);Events ", 1000, 0, 1000);
      h_dijet_eta = new TH1F("dijet_eta", "DiJet \\eta;\\eta;Events ", 188, -9.4, 9.4);
      h_dijet_phi = new TH1F("dijet_phi", "DiJet \\phi;\\phi;Events ", 36, -3.6, 3.6);
      h_dijet_dphi = new TH1F("dijet_dphi", "DiJet dphi;dphi;Events ", 36, -3.6, 3.6);
      h_dijet_deta = new TH1F("dijet_deta", "DiJet deta;deta;Events ", 188, -9.4, 9.4);
      h_met_pt = new TH1F("met_pt", "MET p_{T};p_{T}  (GeV);Events ", 250, 0, 500);
      h_num_vertices = new TH1F("num_vertices", "Number of Vertices;NPV;Events ", 50, 0, 50);
      h_eweight = new TH1F("pu_weight", "Pileup Weight", 99, 0, 99);
      h_muon_pt->Sumw2();
      h_muon_corrpt->Sumw2();
      h_leadMuon_pt->Sumw2();
      h_leadMuon_eta->Sumw2();
      h_leadMuon_phi->Sumw2();
      h_subMuon_pt->Sumw2();
      h_subMuon_eta->Sumw2();
      h_subMuon_phi->Sumw2();
      h_dimuon_mass->Sumw2();
      h_dimuon_pt->Sumw2();
      h_dimuon_eta->Sumw2();
      h_dimuon_phi->Sumw2();
      h_dimuon_deta->Sumw2();
      h_dimuon_dphi->Sumw2();
      h_num_jets->Sumw2();
      h_num_bjets->Sumw2();
      h_leadjet_pt->Sumw2();
      h_leadjet_eta->Sumw2();
      h_leadjet_phi->Sumw2();
      h_subjet_pt->Sumw2();
      h_subjet_eta->Sumw2();
      h_subjet_phi->Sumw2();
      h_dijet_pt->Sumw2();
      h_dijet_mass->Sumw2();
      h_dijet_eta->Sumw2();
      h_dijet_phi->Sumw2();
      h_dijet_dphi->Sumw2();
      h_dijet_deta->Sumw2();
      h_met_pt->Sumw2();
      h_num_vertices->Sumw2();

      GetOutputList()->Add(h_muon_pt);
      GetOutputList()->Add(h_muon_corrpt);
      GetOutputList()->Add(h_leadMuon_pt);
      GetOutputList()->Add(h_leadMuon_eta);
      GetOutputList()->Add(h_leadMuon_phi);
      GetOutputList()->Add(h_subMuon_pt);
      GetOutputList()->Add(h_subMuon_eta);
      GetOutputList()->Add(h_subMuon_phi);
      GetOutputList()->Add(h_dimuon_mass);
      GetOutputList()->Add(h_dimuon_pt);
      GetOutputList()->Add(h_dimuon_eta);
      GetOutputList()->Add(h_dimuon_phi);
      GetOutputList()->Add(h_dimuon_deta);
      GetOutputList()->Add(h_dimuon_dphi);
      GetOutputList()->Add(h_num_jets);
      GetOutputList()->Add(h_num_bjets);
      GetOutputList()->Add(h_leadjet_pt);
      GetOutputList()->Add(h_leadjet_eta);
      GetOutputList()->Add(h_leadjet_phi);
      GetOutputList()->Add(h_subjet_pt);
      GetOutputList()->Add(h_subjet_eta);
      GetOutputList()->Add(h_subjet_phi);
      GetOutputList()->Add(h_dijet_pt);
      GetOutputList()->Add(h_dijet_mass);
      GetOutputList()->Add(h_dijet_eta);
      GetOutputList()->Add(h_dijet_phi);
      GetOutputList()->Add(h_dijet_deta);
      GetOutputList()->Add(h_dijet_dphi);
      GetOutputList()->Add(h_met_pt);
      GetOutputList()->Add(h_num_vertices);
      GetOutputList()->Add(h_eweight);

      for (Int_t i = 0; i <= 1000; i += 10)
      {
         TString histname = Form("dimuon_mass_jet_%d", i);
         TString histTitle = Form("Dimuon Mass, Jet Mass > %s ;M_{\\mu \\mu}  (Gev);Events / bin", to_string(i).c_str());
         vec_dimuon_mass_jets.push_back(new TH1F(histname, histTitle, 100, 100, 150));
	 vec_dimuon_mass_jets.at(i / 10)->Sumw2();
         GetOutputList()->Add(vec_dimuon_mass_jets.at(i / 10));

         TString histname_r = Form("dimuon_mass_jet_r_%d", i);
         TString histTitle_r = Form("Dimuon Mass, Jet Mass < %s ;M_{\\mu \\mu}  (Gev);Events / bin", to_string(i).c_str());
         vec_dimuon_mass_jets_r.push_back(new TH1F(histname_r, histTitle_r, 100, 100, 150));
	 vec_dimuon_mass_jets_r.at(i / 10)->Sumw2();
         GetOutputList()->Add(vec_dimuon_mass_jets_r.at(i / 10));
      }
   }
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

   std::vector<analysis::core::Vertex> vecVertices;
   for (int iVer = 0, nVer = Vertices.GetSize(); iVer < nVer; ++iVer)
   {
      vecVertices.push_back(Vertices[iVer]);
   }
   if (!passVertex(vecVertices))
      return kFALSE;
   if (!(std::any_of(_hasHLTFired->begin(), _hasHLTFired->end(), [](bool v) { return v; })))
      return kFALSE;
   if (!*_passedMetFilters)
      return kFALSE;

   std::vector<std::pair<analysis::core::Muon, analysis::core::Muon>> muonPairs;
   for (int im = 0, nMuons = Muons.GetSize(); im < nMuons; ++im)
   {
      for (int jm = im + 1; jm < nMuons; ++jm)
      {
         if (passMuons(Muons[im], Muons[jm]))
            muonPairs.push_back(std::make_pair(Muons[im], Muons[jm]));
      }
   }

   if (muonPairs.size() == 0)
      return kFALSE;

   float highestPtSum = 0;
   std::pair<analysis::core::Muon, analysis::core::Muon> highestPtMuonPair;
   TLorentzVector highestPtMuonsP4;

   for (const std::pair<analysis::core::Muon, analysis::core::Muon> &twoMuons : muonPairs)
   {
      TLorentzVector p4m1, p4m2;
      p4m1.SetPtEtaPhiM(twoMuons.first._corrPT, twoMuons.first._eta, twoMuons.first._phi, PDG_MASS_Mu);
      p4m2.SetPtEtaPhiM(twoMuons.second._corrPT, twoMuons.second._eta, twoMuons.second._phi, PDG_MASS_Mu);
      TLorentzVector p4dimuon = p4m1 + p4m2;

      if (p4dimuon.Pt() > highestPtSum)
      {
         highestPtSum = p4dimuon.Pt();
         highestPtMuonPair = twoMuons;
         highestPtMuonsP4 = p4dimuon;
      }
   }

   if (highestPtMuonsP4.M() < _dimuonMinMass || highestPtMuonsP4.M() > _dimuonMaxMass)
      return kFALSE;

   float eweight = 1;

   if (mcLabel)
   {
      eweight = ((weighter->weight(*_nPU)) * (*_genWeight) * xsec / valueSumEventsWeighted);
   }

   if (!ntupleModeS)
   {
      h_eweight->Fill(eweight);
      h_num_vertices->Fill(Vertices.GetSize(), eweight);
      h_muon_pt->Fill(highestPtMuonPair.first._pt, eweight);
      h_muon_pt->Fill(highestPtMuonPair.second._pt, eweight);
      h_muon_corrpt->Fill(highestPtMuonPair.first._corrPT, eweight);
      h_muon_corrpt->Fill(highestPtMuonPair.second._corrPT, eweight);

      h_leadMuon_pt->Fill(highestPtMuonPair.first._corrPT, eweight);
      h_leadMuon_phi->Fill(highestPtMuonPair.first._phi, eweight);
      h_leadMuon_eta->Fill(highestPtMuonPair.first._eta, eweight);

      h_subMuon_pt->Fill(highestPtMuonPair.second._corrPT, eweight);
      h_subMuon_phi->Fill(highestPtMuonPair.second._phi, eweight);
      h_subMuon_eta->Fill(highestPtMuonPair.second._eta, eweight);

      h_dimuon_mass->Fill(highestPtMuonsP4.M(), eweight);
      h_dimuon_pt->Fill(highestPtMuonsP4.Pt(), eweight);
      h_dimuon_eta->Fill(highestPtMuonsP4.Eta(), eweight);
      h_dimuon_phi->Fill(highestPtMuonsP4.Phi(), eweight);
      h_dimuon_deta->Fill(highestPtMuonPair.first._eta - highestPtMuonPair.second._eta, eweight);
      h_dimuon_dphi->Fill(highestPtMuonPair.first._phi - highestPtMuonPair.second._phi, eweight);

      h_met_pt->Fill(*_pt, eweight);
   }

   // Jet Selection

   std::vector<TLorentzVector> p4jets;

   Int_t _btagJets = 0;
   Int_t _ncentJets = 0;
   Int_t _nfwdJets = 0;
   Int_t _numJets = 0;

   for (analysis::core::Jet iJet : Jets)
   {
      if (iJet._pt > _JetPt && TMath::Abs(iJet._eta) < _JetEta && passTightJetID(iJet) && passLoosePUID(iJet._fullid))
      {
         if ((jetMuondR(iJet._eta, iJet._phi, highestPtMuonPair.first._eta, highestPtMuonPair.first._phi) > 0.4) && (jetMuondR(iJet._eta, iJet._phi, highestPtMuonPair.second._eta, highestPtMuonPair.second._phi) > 0.4))
         {
            if (iJet._btag[0] > 0.4941)
               _btagJets++;
            if (TMath::Abs(iJet._eta) <= 2.4)
               _ncentJets++;
            else
               _nfwdJets++;
            TLorentzVector p4;
            p4.SetPtEtaPhiM(iJet._pt, iJet._eta, iJet._phi, iJet._mass);
            p4jets.push_back(p4);
         }
      }
   }
   _numJets = p4jets.size();

   if (!ntupleModeS)
   {
      h_num_jets->Fill(_numJets, eweight);
      h_num_bjets->Fill(_btagJets, eweight);
   }
   TLorentzVector leadJet, subJet, diJet;

   if (_numJets == 1)
      leadJet = p4jets[0];

   else if (_numJets >= 2)
   {
      for (unsigned int i = 0; i < p4jets.size(); ++i)
      {
         for (unsigned int j = i + 1; j < p4jets.size(); ++j)
         {
            TLorentzVector p4lead = p4jets[i];
            TLorentzVector p4sub = p4jets[j];
            TLorentzVector p4dijet = p4lead + p4sub;
            if (p4dijet.M() > diJet.M())
            {
               leadJet = p4lead;
               subJet = p4sub;
               diJet = p4dijet;
            }
         }
      }
      // remove the jets, to later find second jet pair with high mass
      p4jets.erase(std::remove(p4jets.begin(), p4jets.end(), leadJet));
      p4jets.erase(std::remove(p4jets.begin(), p4jets.end(), subJet));
   }

   TLorentzVector leadJet2, subJet2, diJet2;
   // look for next highest mass dijet pair, save as leadJet2, subJet2 and diJet2
   // First jet should be deleted by now!
   // TO DO: Create Histograms for the second jet pair
   if (p4jets.size() >= 2)
   {
      for (unsigned int i = 0; i < p4jets.size(); ++i)
      {
         for (unsigned int j = i + 1; j < p4jets.size(); ++j)
         {
            TLorentzVector p4lead2 = p4jets[i];
            TLorentzVector p4sub2 = p4jets[j];
            TLorentzVector p4dijet2 = p4lead2 + p4sub2;
            if (p4dijet2.M() > diJet2.M())
            {
               leadJet2 = p4lead2;
               subJet2 = p4sub2;
               diJet2 = p4dijet2;
            }
         }
      }
   }

   if (ntupleModeS)
   {
      float h_mass = highestPtMuonsP4.M();
      float h_pt = highestPtMuonsP4.Pt();
      float h_eta = highestPtMuonsP4.Eta();
      float h_phi = highestPtMuonsP4.Phi();
      float h_deta = TMath::Abs(highestPtMuonPair.first._eta - highestPtMuonPair.second._eta);
      float h_dphi = TMath::Abs(highestPtMuonPair.first._phi - highestPtMuonPair.second._phi);
      float mupt_1 = highestPtMuonPair.first._corrPT;
      float mueta_1 = highestPtMuonPair.first._eta;
      float muphi_1 = highestPtMuonPair.first._phi;
      float mupt_2 = highestPtMuonPair.second._corrPT;
      float mueta_2 = highestPtMuonPair.second._eta;
      float muphi_2 = highestPtMuonPair.second._phi;
      // Int_t njets = p4jets.size();

      float jetpt_1 = leadJet.Pt();
      float jetmass_1 = leadJet.M();
      float jeteta_1 = leadJet.Pt() ? leadJet.Eta() : -5;

      float jetpt_2 = subJet.Pt();
      float jetmass_2 = subJet.M();
      float jeteta_2 = subJet.Pt() ? subJet.Eta() : -5;

      float mjj_1 = diJet.M() ? diJet.M() : 0;
      float detajj_1 = diJet.M() ? TMath::Abs(leadJet.Eta() - subJet.Eta()) : -1;

      float mjj_2 = diJet2.M() ? diJet2.M() : 0;
      float detajj_2 = diJet2.M() ? TMath::Abs(leadJet2.Eta() - subJet2.Eta()) : -1;

      float toFill[] = {
          static_cast<float>(*_run),
          static_cast<float>(*_lumi),
          static_cast<float>(*_event),
          static_cast<float>(mcLabel),
          eweight,
          h_mass,
          h_pt,
          h_eta,
          h_phi,
          h_deta,
          h_dphi,
          mupt_1,
          mueta_1,
          muphi_1,
          mupt_2,
          mueta_2,
          muphi_2,
          static_cast<float>(_numJets),
          static_cast<float>(_ncentJets),
          static_cast<float>(_nfwdJets),
          static_cast<float>(_btagJets),
          jetpt_1,
          jetmass_1,
          jeteta_1,
          jetpt_2,
          jetmass_2,
          jeteta_2,
          mjj_1,
          detajj_1,
          mjj_2,
          detajj_2,
          *_pt};

      ntuple->Fill(toFill);
   }

   if (!ntupleModeS)
   {
      // make sure the jets are not 0 (so they exist)
      // can just do jet.Pt() to know it exists as well
      if (leadJet.Pt())
      {
         h_leadjet_pt->Fill(leadJet.Pt(), eweight);
         h_leadjet_eta->Fill(leadJet.Eta(), eweight);
         h_leadjet_phi->Fill(leadJet.Phi(), eweight);
         if (subJet.Pt())
         {
            h_subjet_pt->Fill(subJet.Pt(), eweight);
            h_subjet_eta->Fill(subJet.Eta(), eweight);
            h_subjet_phi->Fill(subJet.Phi(), eweight);

            if (diJet.M() > 0)
            {
               h_dijet_pt->Fill(diJet.Pt(), eweight);
               h_dijet_mass->Fill(diJet.M(), eweight);
               h_dijet_eta->Fill(diJet.Eta(), eweight);
               h_dijet_phi->Fill(diJet.Phi(), eweight);
               h_dijet_deta->Fill(leadJet.Eta() - subJet.Eta(), eweight);
               h_dijet_dphi->Fill(leadJet.Phi() - subJet.Phi(), eweight);
            }
         }
      }

      for (Int_t i = 0; i <= 1000; i += 10)
      {
         if (diJet.M() > i)
         {
            vec_dimuon_mass_jets.at(i / 10)->Fill(highestPtMuonsP4.M(), eweight);
         }
         if (diJet.M() > 0 && diJet.M() < i)
         {
            vec_dimuon_mass_jets_r.at(i / 10)->Fill(highestPtMuonsP4.M(), eweight);
         }
      }
   }

   return kTRUE;
}

void hmumuSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void hmumuSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   TFile fout(_outputNameFinal, "recreate");

   if (ntupleModeM)
   {
      ntuple = dynamic_cast<TNtuple *>(fOutput->FindObject("ntupledData"));
      if (ntuple)
      {
         ntuple->Write();
      }
   }
   else
   {
      TList *output_list = (TList *)GetOutputList();
      for (const auto &&obj : *output_list)
      {
         if (obj->IsA()->InheritsFrom("TH1"))
            obj->Write();
      }
   }

   fout.Close();
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

bool hmumuSelector::passMuon(analysis::core::Muon const &m)
{
   double muonIsolation = (m._sumChargedHadronPtR04 + std::max(0., m._sumNeutralHadronEtR04 + m._sumPhotonEtR04 - 0.5 * m._sumPUPtR04)) / m._corrPT;

   if (m._isGlobal && m._isTracker &&
       m._corrPT > _muonPt && TMath::Abs(m._eta) < _muonEta &&
       m._isMedium && muonIsolation < _muonIso)
      return true;
   return false;
}

bool hmumuSelector::passMuonHLT(analysis::core::Muon const &m)
{
   if ((m._isHLTMatched[1] || m._isHLTMatched[0]) && m._corrPT > _muonMatchedPt && TMath::Abs(m._eta) < _muonMatchedEta)
      return true;

   return false;
}

bool hmumuSelector::passMuons(analysis::core::Muon const &m1, analysis::core::Muon const &m2)
{
   if ((m1._charge != m2._charge) && passMuon(m1) && passMuon(m2))
   {
      if (passMuonHLT(m1) || passMuonHLT(m2))
      {
         TLorentzVector p4m1, p4m2;
         p4m1.SetPtEtaPhiM(m1._pt, m1._eta, m1._phi, PDG_MASS_Mu);
         p4m2.SetPtEtaPhiM(m2._pt, m2._eta, m2._phi, PDG_MASS_Mu);
         TLorentzVector p4dimuon = p4m1 + p4m2;

         if (p4dimuon.M() > _dimuonMinMass && p4dimuon.M() < _dimuonMaxMass)
            return true;
      }
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

bool hmumuSelector::passLoosePUID(int jetfullID)
{
   return bool(jetfullID & (1 << 2));
}
