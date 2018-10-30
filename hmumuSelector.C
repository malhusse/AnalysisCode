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

double _muonMatchedPt = 30.;
double _muonMatchedEta = 2.4;
double _muonPt = 20.;
double _muonEta = 2.4;
double _muonIso = 0.25;
// double _dimuonMinMass = 80.;
// double _dimuonMaxMass = 85.;
// double _JetPt = 30.;
// double _JetEta = 4.7;

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
   h_dimuon_deta = new TH1F("dimuon_deta","Dimuon deta;deta;Events / bin",50,-2.5,2.5);
   h_dimuon_dphi = new TH1F("dimuon_dphi","Dimuon dphi;dphi;Events / bin",36, -3.6, 3.6);
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
   h_dijet_dphi = new TH1F("dijet_dphi", "DiJet dphi;dphi;Events / bin",36,-3.6,3.6);
   h_dijet_deta = new TH1F("dijet_deta","DiJet deta;deta;Events / bin",94,-4.7,4.7);
   h_met_pt = new TH1F("met_pt","MET p_{T};p_{T}  (GeV);Events / bin",250,0,500);
   h_met_eta = new TH1F("met_eta","MET \\eta;\\eta;Events / bin",94,-4.7,4.7);
   h_met_phi = new TH1F("met_phi","MET \\phi;\\phi;Events / bin",36,-3.6,3.6);
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
   GetOutputList()->Add(h_met_eta);
   GetOutputList()->Add(h_met_phi);
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


   if (!passMuons(Muons))
      return kFALSE;

   // only two muons...
   // for (int iMuon = 0, nMuons = Muons__charge.GetSize(); iMuon < nMuons; ++iMuon){

   h_num_vertices->Fill(Vertices.GetSize());
   h_muon_pt->Fill(Muons[0]._pt);
   h_muon_pt->Fill(Muons[1]._pt);
   h_muon_corrpt->Fill(Muons[0]._corrPT);
   h_muon_corrpt->Fill(Muons[1]._corrPT);
   // }


   int lead_muon_id = 0;
   int sub_muon_id = 1;

   if (Muons[1]._corrPT > Muons[0]._corrPT){
      lead_muon_id = 1;
      sub_muon_id = 0;
   }

   h_leadMuon_pt->Fill(Muons[lead_muon_id]._corrPT);
   h_subMuon_pt->Fill(Muons[sub_muon_id]._corrPT);

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

   for (int iVert = 0; iVert < vertexCol.GetSize(); ++iVert)
	{
		if (TMath::Abs(vertexCol[iVert]._z) < 24 &&
			vertexCol[iVert]._ndf > 4)
			return true;
	}

	return false;
}


bool passMuon(analysis::core::Muon m)
{
	double muonIsolation = (m._sumChargedHadronPtR04 + std::max(0.,
																m._sumNeutralHadronEtR04 + m._sumPhotonEtR04 - 0.5 * m._sumPUPtR04)) /
						   m._corrPT;

	if (m._isGlobal && m._isTracker &&
		m._corrPT > _muonPt && TMath::Abs(m._eta) < _muonEta &&
		m._isMedium && muonIsolation < _muonIso)
		return true;
	return false;
}

bool passMuonHLT(analysis::core::Muon m)
{
	if ((m._isHLTMatched[1] || m._isHLTMatched[0]) && m._corrPT > _muonMatchedPt && TMath::Abs(m._eta) < _muonMatchedEta)
		return true;

	return false;
}

bool passMuons(TTreeReaderArray<analysis::core::Muon> muonCol);
{
	if ((muonCol[0]._charge != muonCol[1]._charge) && passMuon(muonCol[0]) && passMuon(muonCol[1]))
	{
		if (passMuonHLT(muonCol[0]) || passMuonHLT(muonCol[0]))
		{
			// TLorentzVector p4m1, p4m2;
			// p4m1.SetPtEtaPhiM(m1._pt, m1._eta, m1._phi, PDG_MASS_Mu);
			// p4m2.SetPtEtaPhiM(m2._pt, m2._eta, m2._phi, PDG_MASS_Mu);
			// TLorentzVector p4dimuon = p4m1 + p4m2;

			// if (p4dimuon.M() > _dimuonMinMass && p4dimuon.M() < _dimuonMaxMass)
			return true;
		}
	}

	return false;
}

// float jetMuondR(float jeta, float jphi, float meta, float mphi)
// {
// 	TLorentzVector p4j, p4m;
// 	p4j.SetPtEtaPhiM(10, jeta, jphi, 0);
// 	p4m.SetPtEtaPhiM(10, meta, mphi, 0);
// 	return p4j.DeltaR(p4m);
// }

// bool passElectronVeto(Electrons *electrons, Muon m1, Muon m2)
// {
// 	for (Electrons::const_iterator it = electrons->begin();
// 		 it != electrons->end(); ++it)
// 	{
// 		if (it->_ids[2] && it->_pt > 10. &&
// 			(TMath::Abs(it->_eta) < 1.4442 || (1.566 < TMath::Abs(it->_eta) && TMath::Abs(it->_eta) < 2.5)) && jetMuondR(it->_eta, it->_phi, m1._eta, m1._phi) > 0.4 && jetMuondR(it->_eta, it->_phi, m2._eta, m2._phi) > 0.4)
// 			return false;
// 	}

// 	return true;
// }

// bool passBTaggedJetVeto(Jets *jets)
// {
// 	for (Jets::const_iterator it = jets->begin(); it != jets->end(); ++it)
// 		if (it->_pt > 30. & TMath::Abs(it->_eta) < 2.4 && it->_btag[0] > 0.4941)
// 			return false;
// 	return true;
// }
// bool passTightJetID(Jet const &j)
// {
// 	bool tightID = false;
// 	double jeta = TMath::Abs(j._eta);
// 	int numConst = j._cm + j._nm;

// 	if (jeta <= 2.7)
// 	{
// 		tightID = (j._nhf < 0.90 && j._nef < 0.90 && numConst > 1);

// 		if (jeta < 2.4)
// 		{
// 			tightID &= (j._chf > 0 && j._cm > 0);
// 		}
// 	}
// 	else if (jeta <= 3.0)
// 	{
// 		tightID = (j._nef > 0.02 && j._nef < 0.99 && j._nm > 2);
// 	}
// 	else
// 	{
// 		tightID = (j._nef < 0.90 && j._nhf > 0.02 && j._nm > 10);
// 	}

// 	return tightID;
// }

// bool passLoosePUID(Jet const &j)
// {
// 	return (j._fullid & (1 << 2));
// }
