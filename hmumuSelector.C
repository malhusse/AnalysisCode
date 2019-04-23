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
   TNamed *name = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));

   TParameter<Int_t> *pYearM = dynamic_cast<TParameter<Int_t>*>(fInput->FindObject("getYear"));
   yearM = pYearM->GetVal();

   _outputNameFinal = "ntupleFiles/";
   _outputNameFinal += std::to_string(yearM);
   _outputNameFinal += "/";
   _outputNameFinal += name ? name->GetTitle() : "outputName.root";
}

void hmumuSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   TNamed *name2 = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));
   _outputRoot = name2->GetTitle();

   TParameter<Int_t> *pSumEvents = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getSumEvents"));
   TParameter<Int_t> *pSumEventsWeighted = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getSumEventsWeighted"));

   valueSumEvents = pSumEvents->GetVal();
   valueSumEventsWeighted = pSumEventsWeighted->GetVal();

   TParameter<Double_t> *pxsec = dynamic_cast<TParameter<Double_t> *>(fInput->FindObject("getxsec"));
   TParameter<Int_t> *pmcLabel = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getmcLabel"));
   TParameter<Int_t> *pYearS = dynamic_cast<TParameter<Int_t>*>(fInput->FindObject("getYear"));

   mcLabel = pmcLabel->GetVal();
   xsec = static_cast<Float_t>(pxsec->GetVal());
   yearS = pYearS->GetVal();

      // mcLabel should be 0 for data, so this is true only for mc!
   if (mcLabel)
   {
      _dataPUfile = "/uscms_data/d1/malhusse/analysis/AnalysisCode/data/pileup/pu_data_";
      _dataPUfile += std::to_string(yearS);
      _dataPUfile += ".root";
      // _dataPUfile = "/uscms_data/d1/malhusse/build/AnalysisCode/pileup/pu_data_2017.root";
      _mcPUfile = "/uscms_data/d1/malhusse/analysis/AnalysisCode/data/pileup/";
      _mcPUfile += std::to_string(yearS);
      _mcPUfile += _outputRoot;

      weighter = new reweight::LumiReWeighting(_mcPUfile.Data(), _dataPUfile.Data(), "pileup", "pileup");
   }
               //   "muRoch_1:muRoch_2:mu_kinfit_pt_1:mu_kinfit_pt_2"
               //   "fsrPt_1:fsrEta_1:fsrPhi_1:fsrPt_2:fsrEta_2:fsrPhi_2"
               //   "muPt_1:muEta_1:muPhi_1:muPt_2:muEta_2:muPhi_2:"

   string vars = "year:run:lumi:event:mclabel:genWeight:numGen:xsec:puWeight"
                 "idSF:isoSF:btagSF:trigSF:prefireSF:totalSF:eWeight:totalWeight"
                 "muPtC_1:muEtaC_1:muPhiC_1"
                 "muPtC_2:muEtaC_2:muPhiC_2"
                 "h_mass:h_pt:h_eta:h_phi:h_deta:h_dphi:"
                 "njets:ncentJets:nfwdJets:nbtagJets:maxbdisc"
                 "jetpt_1:jetmass_1:jeteta_1:jetpt_2:jetmass_2:jeteta_2:"
                 "mjj_1:detajj_1:mjj_2:detajj_2:"
                 "metpt";
   
   ntuple = new TNtuple("ntupledData", "Data TNtuple", vars.c_str());
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

   std::vector<analysis::core::Vertex> vecVertices;
   std::vector<std::pair<string,int>> metFilterBits;

   for (int iVer = 0, nVer = Vertices.GetSize(); iVer < nVer; ++iVer)
   {
      vecVertices.push_back(Vertices[iVer]);
   }
   
   for (int iBit = 0, nBits = _metFilterBits.GetSize(); iBit < nBits; ++iBit)
   {
      metFilterBits.push_back(_metFilterBits[iBit]);
   }

   if (!passVertex(vecVertices))
      return kFALSE;
   if (!(std::any_of(_hasHLTFired->begin(), _hasHLTFired->end(), [](bool v) { return v; })))
      return kFALSE;
   if (!passMetFilters(metFilterBits))
      return kFALSE;

   // correct Muons using kinfit pt * rochester factor + fsr photon
   for (int iMu = 0, nMu = Muons.GetSize(); iMu < nMu; ++iMu)
   {
      TLorentzVector p4mu, p4pho;
      p4mu.SetPtEtaPhiM((Muons[iMu]._pt_kinfit * Muons[iMu]._roccCor), Muons[iMu]._eta,Muons[iMu]._phi,PDG_MASS_Mu);
      p4pho = Muons[iMu].fsrP4;
      p4mu += p4pho;
      Muons[iMu]._pt = p4mu.Pt();
      Muons[iMu]._eta = p4mu.Eta();
      Muons[iMu]._phi = p4mu.Phi();
   }

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
      p4m1.SetPtEtaPhiM(twoMuons.first._pt, twoMuons.first._eta, twoMuons.first._phi, PDG_MASS_Mu);
      p4m2.SetPtEtaPhiM(twoMuons.second._pt, twoMuons.second._eta, twoMuons.second._phi, PDG_MASS_Mu);
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
   
   // Jets

   std::vector<TLorentzVector> p4jets;

   int _btagJets = 0;
   int _ncentJets = 0;
   int _nfwdJets = 0;
   int _numJets = 0;

   float maxBDisc = 0.0;
   float btagSF = 1.0;

   for (analysis::core::Jet iJet : Jets)
   {
   
     if (iJet._btag_sf != 0) btagSF*= iJet._btag_sf;

     if (iJet._pt > _JetPt && TMath::Abs(iJet._eta) < _JetEta && passTightJetID(iJet) && passLoosePUID(iJet) && passNoiseJet(iJet))
       // if (iJet._pt > _JetPt && TMath::Abs(iJet._eta) < _JetEta && passTightJetID(iJet))
      {
         if ((jetMuondR(iJet._eta, iJet._phi, highestPtMuonPair.first._eta, highestPtMuonPair.first._phi) > 0.4) && (jetMuondR(iJet._eta, iJet._phi, highestPtMuonPair.second._eta, highestPtMuonPair.second._phi) > 0.4))
         {
            float ibtagD = iJet._btag[0];
            if (ibtagD > 0.4941)
               _btagJets++;
               if (ibtagD > maxBDisc)
                  maxBDisc = ibtagD;
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
   // First jet pair should be deleted by now!
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

   //event weight stuff goes here?
   float pileupWeight = 1.0;
   float totalSF = 1.0;
   float eWeight = 1.0;
   float totalWeight = 1.0;
   if (mcLabel)
   {
      pileupWeight = weighter->weight(*_nPU);
      totalSF = (*_idSF) * (*_isoSF) * btagSF * (*_trigEffSF) * (*_prefiringweight);
      eWeight = pileupWeight * (*_genWeight) * xsec / valueSumEventsWeighted;
      totalWeight = totalSF * eWeight;
   }

   float h_mass = highestPtMuonsP4.M();
   float h_pt = highestPtMuonsP4.Pt();
   float h_eta = highestPtMuonsP4.Eta();
   float h_phi = highestPtMuonsP4.Phi();
   float h_deta = TMath::Abs(highestPtMuonPair.first._eta - highestPtMuonPair.second._eta);
   float h_dphi = TMath::Abs(highestPtMuonPair.first._phi - highestPtMuonPair.second._phi);
   
   float mupt_1 = highestPtMuonPair.first._pt;
   float mueta_1 = highestPtMuonPair.first._eta;
   float muphi_1 = highestPtMuonPair.first._phi;
   float mupt_2 = highestPtMuonPair.second._pt;
   float mueta_2 = highestPtMuonPair.second._eta;
   float muphi_2 = highestPtMuonPair.second._phi;

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
         static_cast<float>(yearS),
         static_cast<float>(*_run),
         static_cast<float>(*_lumi),
         static_cast<float>(*_event),
         static_cast<float>(mcLabel),
         static_cast<float>(*_genWeight),
         static_cast<float>(valueSumEventsWeighted),
         static_cast<float>(xsec),
         pileupWeight,
         *_idSF,
         *_isoSF,
         btagSF,
         *_trigEffSF,
         static_cast<float>(*_prefiringweight),
         totalSF,
         eWeight,
         totalWeight,
         mupt_1,
         mueta_1,
         muphi_1,
         mupt_2,
         mueta_2,
         muphi_2,
         h_mass,
         h_pt,
         h_eta,
         h_phi,
         h_deta,
         h_dphi,
         static_cast<float>(_numJets),
         static_cast<float>(_ncentJets),
         static_cast<float>(_nfwdJets),
         static_cast<float>(_btagJets),
         maxBDisc,
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
   ntuple = dynamic_cast<TNtuple *>(fOutput->FindObject("ntupledData"));
   if (ntuple)
   {
      ntuple->Write();
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
   double muonIsolation = (m._sumChargedHadronPtR04 + std::max(0., m._sumNeutralHadronEtR04 + m._sumPhotonEtR04 - 0.5 * m._sumPUPtR04)) / m._pt;

   if (m._isGlobal && m._isTracker &&
       m._pt > _muonPt && TMath::Abs(m._eta) < _muonEta &&
       m._isMedium && muonIsolation < _muonIso)
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

bool hmumuSelector::passLoosePUID(analysis::core::Jet j)
{
  float jeta = TMath::Abs(j._eta);
  float jpt = j._pt;
  float jpuid = j._puid;

  if(jeta < 2.5)
    {
      if (jpt >=30 and jpt < 50 and jpuid <-0.89)  return false;
      if (jpt >=10 and jpt < 30 and jpuid <-0.97)  return false;
    }
  else if(jeta < 2.75)
    {
      if (jpt >=30 and jpt < 50 and jpuid <-0.52)  return false;
      if (jpt >=10 and jpt < 30 and jpuid <-0.68)  return false;
    }
  else if(jeta < 3.0)
    {
      if (jpt >=30 and jpt < 50 and jpuid <-0.38)  return false;
      if (jpt >=10 and jpt < 30 and jpuid <-0.53)  return false;
    }
  else if(jeta < 5)
    {
      if (jpt >=30 and jpt < 50 and jpuid <-0.30)  return false;
      if (jpt >=10 and jpt < 30 and jpuid <-0.47)  return false;
    }
    return true;
}

bool hmumuSelector::passNoiseJet(analysis::core::Jet j)
{
  float jeta = TMath::Abs(j._eta);
  float jpt = j._pt;

  if (jeta >= 2.65 and jeta <= 3.139 and jpt < 50) return false;
  return true;
}

bool hmumuSelector::passMetFilters(std::vector<std::pair<string,int>> filterBits)
{
   bool pass = true;
   for (const std::pair<string,int> &bit : filterBits)
   {
      if (bit.first.compare("Flag_BadChargedCandidateFilter")) // 0 if equal, so gets skipped!
      {
         pass = pass && bit.second;
      }
   }
   return pass;
}