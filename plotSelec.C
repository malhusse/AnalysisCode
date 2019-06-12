#define plotSelec_cxx
// The class definition in plotSelec.h has been generated automatically
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
// root> T->Process("plotSelec.C")
// root> T->Process("plotSelec.C","some options")
// root> T->Process("plotSelec.C+")
//

#include "plotSelec.h"
#include <TH2.h>
#include <TStyle.h>
#include <TNamed.h>

void plotSelec::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   TNamed *name = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));

   TParameter<Int_t> *pYear = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getYear"));
   collectionYear = pYear->GetVal();

   _outputName = "histoFiles/";
   _outputName += std::to_string(collectionYear);
   _outputName += "/";
   _outputName += name->GetTitle();
}

void plotSelec::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   h_leadMuon_pt = new TH1F("lead_muon_pt", "Leading Muon p_{T};p_{T}  (GeV);Events ", 100, 0, 200);
   h_leadMuon_eta = new TH1F("lead_muon_eta", "Leading Muon \\eta;\\eta;Events ", 50, -2.5, 2.5);

   h_subMuon_pt = new TH1F("sub_muon_pt", "Subleading Muon p_{T};p_{T}  (GeV);Events ", 100, 0, 200);
   h_subMuon_eta = new TH1F("sub_muon_eta", "Subleading Muon eta;\\eta;Events ", 50, -2.5, 2.5);

   h_dimuon_mass = new TH1F("dimuon_mass", "Dimuon Mass;M_{\\mu \\mu}  (Gev);Events ", 100, 100, 150);
   h_dimuon_pt = new TH1F("dimuon_pt", "Dimuon p_{T};p_{T}  (GeV);Events ", 200, 0, 400);
   h_dimuon_eta = new TH1F("dimuon_eta", "Dimuon \\eta;\\eta;Events ", 100, -5.0, 5.0);
   h_dimuon_phi = new TH1F("dimuon_phi", "Dimuon \\phi;\\phi;Events ", 36, -3.6, 3.6);
   h_dimuon_deta = new TH1F("dimuon_deta", "Dimuon deta;deta;Events ", 100, -5.0, 5.0);
   h_dimuon_dphi = new TH1F("dimuon_dphi", "Dimuon dphi;dphi;Events ", 36, -3.6, 3.6);

   h_num_jets = new TH1F("num_jets", "Number of Jets;nJets;Events ", 10, 0, 10);
   h_num_bjets = new TH1F("num_bjets", "Number of B Jets;nBJets;Events ", 10, 0, 10);

   h_leadjet_pt = new TH1F("leadjet_pt", "Leading Jet p_{T};p_{T}  (GeV);Events ", 500, 0, 500);
   h_leadjet_eta = new TH1F("leadjet_eta", "Leading Jet \\eta;\\eta;Events ", 94, -4.7, 4.7);

   h_subjet_pt = new TH1F("subjet_pt", "Subleading Jet p_{T};p_{T}  (GeV);Events ", 500, 0, 500);
   h_subjet_eta = new TH1F("subjet_eta", "Subleading Jet \\eta;\\eta;Events ", 94, -4.7, 4.7);

   h_dijet_mass = new TH1F("dijet_mass1", "DiJet Mass;M_{jj}  (GeV);Events ", 1000, 0, 1000);
   h_dijet_deta = new TH1F("dijet_deta1", "DiJet deta;deta;Events ", 188, -9.4, 9.4);

   h_met_pt = new TH1F("met_pt", "MET p_{T};p_{T}  (GeV) ", 250, 0, 500);
   h_mindrmj = new TH1F("mindrmj", "min(dR(m,j));dR(m,j) ", 50, -5, 5);
   h_zeppen = new TH1F("zeppen", "; zeppen; ", 100, -10, 10);
   h_csTheta = new TH1F("csTheta", " ; csTheta; ", 10, -1, 1);
   h_csPhi = new TH1F("csPhi", " ; csPhi ; ", 100, -10, 10);
   h_bdtScore01jet = new TH1F("h_bdtScore01jet", "; ; ", 20, -1, 1);
   h_bdtScore2jet = new TH1F("h_bdtScore2jet", ";;", 20, -1, 1);

   h_leadMuon_pt->Sumw2();
   h_leadMuon_eta->Sumw2();
   h_subMuon_pt->Sumw2();
   h_subMuon_eta->Sumw2();
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
   h_subjet_pt->Sumw2();
   h_subjet_eta->Sumw2();
   h_dijet_mass->Sumw2();
   h_dijet_deta->Sumw2();
   h_met_pt->Sumw2();
   h_mindrmj->Sumw2();
   h_zeppen->Sumw2();
   h_csTheta->Sumw2();
   h_csPhi->Sumw2();
   h_bdtScore01jet->Sumw2();
   h_bdtScore2jet->Sumw2();

   GetOutputList()->Add(h_leadMuon_pt);
   GetOutputList()->Add(h_leadMuon_eta);
   GetOutputList()->Add(h_subMuon_pt);
   GetOutputList()->Add(h_subMuon_eta);
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
   GetOutputList()->Add(h_subjet_pt);
   GetOutputList()->Add(h_subjet_eta);
   GetOutputList()->Add(h_dijet_mass);
   GetOutputList()->Add(h_dijet_deta);
   GetOutputList()->Add(h_met_pt);
   GetOutputList()->Add(h_mindrmj);
   GetOutputList()->Add(h_zeppen);
   GetOutputList()->Add(h_csTheta);
   GetOutputList()->Add(h_csPhi);
   GetOutputList()->Add(h_bdtScore01jet);
   GetOutputList()->Add(h_bdtScore2jet);
}

Bool_t plotSelec::Process(Long64_t entry)
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
   if (*h_mass < 100 || *h_mass > 150)
      return kTRUE;

   float weight = *totalWeight;
   h_leadMuon_pt->Fill(*muPtC_1, weight);
   h_leadMuon_eta->Fill(*muEtaC_1, weight);

   h_subMuon_pt->Fill(*muPtC_2, weight);
   h_subMuon_eta->Fill(*muEtaC_2, weight);

   h_dimuon_mass->Fill(*h_mass, weight);
   h_dimuon_pt->Fill(*h_pt, weight);
   h_dimuon_eta->Fill(*h_eta, weight);
   h_dimuon_phi->Fill(*h_phi, weight);
   h_dimuon_deta->Fill(*h_deta, weight);
   h_dimuon_dphi->Fill(*h_dphi, weight);

   h_num_jets->Fill(*njets, weight);
   h_num_bjets->Fill(*nbtagJets, weight);

   h_leadjet_pt->Fill(*jetpt_1, weight);
   h_leadjet_eta->Fill(*jeteta_1, weight);

   h_subjet_pt->Fill(*jetpt_2, weight);
   h_subjet_eta->Fill(*jeteta_2, weight);

   h_dijet_mass->Fill(*mjj_1, weight);
   h_dijet_deta->Fill(*detajj_1, weight);

   h_met_pt->Fill(*metpt, weight);
   h_mindrmj->Fill(*mindrmj,weight);
   h_zeppen->Fill(*zeppen,weight);
   h_csTheta->Fill(*csTheta,weight);
   h_csPhi->Fill(*csPhi,weight);

   // insert if here for different event category!

   if (*category == 3)
      h_bdtScore01jet->Fill(*bdtScore,weight);
   if (*category == 13)
      h_bdtScore2jet->Fill(*bdtScore,weight);
      
   return kTRUE;
}

void plotSelec::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void plotSelec::Terminate()
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
}
