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
}

void hmumuSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   h_muon_pt = new TH1F("muon_pt", "Muon p_{T};p_{T}  (GeV);Events / bin", 200,0,400);
   h_muon_corrpt = new TH1F("muon_corrpt", "Muon corrected p_{T};p_{T}  (GeV);Events / bin", 200,0,400);
   h_leadMuon_pt = new TH1F("lead_muon_pt", "Leading Muon p_{T};p_{T}  (GeV);Events / bin", 100,0,200);
   h_leadMuon_eta = new TH1F("lead_muon_eta", "Leading Muon \\eta;\\eta;Events / bin", 50, -2.5, 2.5);
   h_leadMuon_phi = new TH1F("lead_muon_phi", "Leading Muon \\phi;\\phi;Events / bin",36,-3.6,3.6);
   h_subMuon_pt = new TH1F("sub_muon_pt", "Subleading Muon p_{T};p_{T}  (GeV);Events / bin", 100, 0,200);
   h_subMuon_eta = new TH1F("sub_muon_eta", "Subleading Muon eta;\\eta;Events / bin", 50, -2.5, 2.5);
   h_subMuon_phi = new TH1F("sub_muon_phi", "Subleading Muon \\phi;\\phi;Events / bin",36,-3.6,3.6);
   h_dimuon_mass = new TH1F("dimuon_mass","Dimuon Mass;M_{\\mu \\mu}  (Gev);Events / bin",100,100,200);
   h_dimuon_pt = new TH1F("dimuon_pt","Dimuon p_{T};p_{T}  (GeV);Events / bin",200,0,400);
   h_dimuon_eta = new TH1F("dimuon_eta","Dimuon \\eta;\\eta;Events / bin",50,-2.5,2.5);
   h_dimuon_phi = new TH1F("dimuon_phi","Dimuon \\phi;\\phi;Events / bin",36,-3.6,3.6);
   h_dimuon_deta = new TH1F("dimuon_deta","Dimuon deta;deta;Events / bin",25,0,2.5);
   h_dimuon_dphi = new TH1F("dimuon_dphi","Dimuon dphi;dphi;Events / bin",18, 0, 3.6);
   h_num_jets = new TH1F("num_jets","Number of Jets;nJets;Events / bin",10,0,10);
   h_num_bjets = new TH1F("num_bjets","Number of B Jets;nBJets;Events / bin",10,0,10);
   h_leadjet_pt = new TH1F("leadjet_pt","Leading Jet p_{T};p_{T}  (GeV);Events / bin",250,0,500);
   h_leadjet_eta = new TH1F("leadjet_eta","Leading Jet \\eta;\\eta;Events / bin",94,-4.7,4.7);
   h_leadjet_phi = new TH1F("leadjet_phi", "Leading Jet \\phi;\\phi;Events / bin",36,-3.6,3.6);
   h_subjet_pt = new TH1F("subjet_pt","Subleading Jet p_{T};p_{T}  (GeV);Events / bin",250,0,500);
   h_subjet_eta = new TH1F("subjet_eta","Subleading Jet \\eta;\\eta;Events / bin",94,-4.7,4.7);
   h_subjet_phi = new TH1F("subjet_eta", "Subleading Jet \\phi;\\phi;Events / bin",36,-3.6,3.6);
   h_dijet_pt = new TH1F("dijet_pt","DiJet p_{T};p_{T}  (GeV),Events / bin",500,0,1000);
   h_dijet_mass = new TH1F("dijet_mass","DiJet Mass;M_{jj}  (GeV);Events / bin",500,0,1000);
   h_dijet_eta = new TH1F("dijet_eta", "DiJet \\eta;\\eta;Events / bin",94,-4.7,4.7);
   h_dijet_phi = new TH1F("dijet_phi","DiJet \\phi;\\phi;Events / bin",36,-3.6,3.6);
   h_dijet_dphi = new TH1F("dijet_dphi", "DiJet dphi;dphi;Events / bin",18,0,3.6);
   h_dijet_deta = new TH1F("dijet_deta","DiJet deta;deta;Events / bin",47,0,4.7);
   h_met_pt = new TH1F("met_pt","MET p_{T};p_{T}  (GeV);Events / bin",250,0,500);
   h_num_vertices = new TH1F("num_vertices","Number of Vertices;NPV;Events / bin",50,0,50);

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

   // int nVert = Vertices__z.GetSize();
   if (!passVertex(Vertices))
      return kFALSE;
   if(!(std::any_of(_hasHLTFired->begin(), _hasHLTFired->end(), [](bool v) {return v;})))
      return kFALSE;
   if(!*_passedMetFilters)
      return kFALSE;


   std::vector<std::pair<analysis::core::Muon, analysis::core::Muon>> muonPairs;
   for(int im = 0, nMuons = Muons.GetSize(); im < nMuons; ++im){
     for(int jm = im+1; jm < nMuons; ++jm){
       if(passMuons(Muons[im],Muons[jm]))
	 muonPairs.push_back(std::make_pair(Muons[im],Muons[jm]));
     }
   }
     // for(analysis::core::Muon it = (analysis::core::Muon) Muons.begin(); it != Muons.end(); ++it)
     //for(analysis::core::Muon jt = (analysis::core::Muon) (it + 1); jt != Muons.end(); ++jt)
     //if (passMuons(it,jt))
     // if (it->_corrPT > jt->_corrPT)
       //	   muonPairs.push_back(std::make_pair(it, jt));
       //else
       //  muonPairs.push_back(std::make_pair(jt, it));
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

   h_num_vertices->Fill(Vertices.GetSize());
   h_muon_pt->Fill(highestPtMuonPair.first._pt);
   h_muon_pt->Fill(highestPtMuonPair.second._pt);
   h_muon_corrpt->Fill(highestPtMuonPair.first._corrPT);
   h_muon_corrpt->Fill(highestPtMuonPair.second._corrPT);

   int lead_muon_id = 0;
   int sub_muon_id = 1;

   h_leadMuon_pt->Fill(highestPtMuonPair.first._corrPT);
   h_leadMuon_phi->Fill(highestPtMuonPair.first._phi);
   h_leadMuon_eta->Fill(highestPtMuonPair.first._eta);

   h_subMuon_pt->Fill(highestPtMuonPair.second._corrPT);
   h_subMuon_phi->Fill(highestPtMuonPair.second._phi);
   h_subMuon_eta->Fill(highestPtMuonPair.second._eta);

   h_dimuon_mass->Fill(highestPtMuonsP4.M());
   h_dimuon_pt->Fill(highestPtMuonsP4.Pt());
   h_dimuon_eta->Fill(highestPtMuonsP4.Eta());
   h_dimuon_phi->Fill(highestPtMuonsP4.Phi());
   h_dimuon_deta->Fill(TMath::Abs(highestPtMuonPair.first._eta - highestPtMuonPair.second._eta));
   h_dimuon_dphi->Fill(TMath::Abs(highestPtMuonPair.first._phi - highestPtMuonPair.second._phi));

   h_met_pt->Fill(*_pt);
   // Jet Selection

   std::vector<TLorentzVector> p4jets;

   int _btagJets = 0;
   for (analysis::core::Jet iJet: Jets){
     if(iJet._pt > _JetPt && TMath::Abs(iJet._eta) < _JetEta && passTightJetID(iJet) && passLoosePUID(iJet._fullid)){
       if((jetMuondR(iJet._eta,iJet._phi,highestPtMuonPair.first._eta,highestPtMuonPair.first._phi) > 0.4) && (jetMuondR(iJet._eta,iJet._phi,highestPtMuonPair.second._eta,highestPtMuonPair.second._phi) > 0.4))
	 {
	   if (iJet._btag[0] > 0.4941)
	     _btagJets++;
	   TLorentzVector p4;
	   p4.SetPtEtaPhiM(iJet._pt, iJet._eta, iJet._phi, iJet._mass);
	   p4jets.push_back(p4);
	 }
     } 
   }

   h_num_jets->Fill(p4jets.size());
   h_num_bjets->Fill(_btagJets);
   
   TLorentzVector leadJet, subJet, diJet;

   if (p4jets.size() == 1)
     leadJet = p4jets[0];
   
   else if (p4jets.size() >= 2){
     for (unsigned int i = 0; i < p4jets.size(); ++i){
       for (unsigned int j = i + 1; j < p4jets.size(); ++j){
	 TLorentzVector p4lead = p4jets[i];
	 TLorentzVector p4sub = p4jets[j];
	 TLorentzVector p4dijet = p4lead + p4sub;
	 if (p4dijet.M() > diJet.M()){
	   leadJet = p4lead;
	   subJet = p4sub;
	   diJet = p4dijet;
	 }
       }
     }
   }
   
   
   if (leadJet.Pt() > 30){
     h_leadjet_pt->Fill(leadJet.Pt());
     h_leadjet_eta->Fill(leadJet.Eta());
     h_leadjet_phi->Fill(leadJet.Phi());
     if (subJet.Pt() > 30){
       h_subjet_pt->Fill(subJet.Pt());
       h_subjet_eta->Fill(subJet.Eta());
       h_subjet_phi->Fill(subJet.Phi());
       
       if (diJet.M() > 1){
	 h_dijet_pt->Fill(diJet.Pt());
	 h_dijet_mass->Fill(diJet.M());
	 h_dijet_eta->Fill(diJet.Eta());
	 h_dijet_phi->Fill(diJet.Phi());
	 h_dijet_deta->Fill(TMath::Abs(leadJet.Eta()-subJet.Eta()));
	 h_dijet_dphi->Fill(TMath::Abs(leadJet.Phi()-subJet.Phi()));
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
  //hmuon_pt = dynamic_cast<TH1D*>(fOutput->FindObject(Form("muon_pt")));
  //TFile fout("test.root","recreate");
  //hmuon_pt->Write();
  //fout.Close();
  // TFile output("processed_ntuples.root","recreate");
   //output.Write();
   TList *output_list = (TList*)GetOutputList();
   TFile fout("processed.root","recreate");
   //TIter iter(output_list);
   //std::for_each(iter.Begin(), TIter::End(), writeObj());
   for(const auto&& obj: *output_list){
      if(obj->IsA()->InheritsFrom("TH1"))
      obj->Write();
   }
   fout.Close();
}


// implement the following::


// pass vertices
bool hmumuSelector::passVertex(TTreeReaderArray<analysis::core::Vertex> vertexCol)
{
	if (vertexCol.GetSize() == 0)
		return false;

	for (unsigned int iVert = 0; iVert < vertexCol.GetSize(); ++iVert)
	{
		if (TMath::Abs(vertexCol[iVert]._z) < 24 &&
			vertexCol[iVert]._ndf > 4)
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
  return (jetfullID & (1 << 2));
}
