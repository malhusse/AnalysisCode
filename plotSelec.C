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

   //   _outputName = "histoFilesH/";
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
   // h_leadMuon_pt_onZ = new TH1F("lead_muon_pt_onZ", "Leading Muon p_{T};p_{T}  (GeV);Events ", 100, 0, 200);
   h_leadMuon_eta_onZ = new TH1F("lead_muon_eta_onZ", "Leading Muon \\eta;\\eta;Events ", 50, -2.5, 2.5);
   // h_subMuon_pt_onZ = new TH1F("sub_muon_pt_onZ", "Subleading Muon p_{T};p_{T}  (GeV);Events ", 100, 0, 200);
   h_subMuon_eta_onZ = new TH1F("sub_muon_eta_onZ", "Subleading Muon eta;\\eta;Events ", 50, -2.5, 2.5);
   h_dimuon_mass_onZ = new TH1F("dimuon_mass_onZ", "Dimuon Mass;M_{\\mu \\mu}  (Gev);Events ", 100, 70, 110);
   h_dimuon_pt_onZ = new TH1F("dimuon_pt_onZ", "Dimuon p_{T};p_{T}  (GeV);Events ", 200, 0, 400);
   h_dimuon_eta_onZ = new TH1F("dimuon_eta_onZ", "Dimuon \\eta;\\eta;Events ", 100, -5.0, 5.0);
   // h_dimuon_phi_onZ = new TH1F("dimuon_phi_onZ", "Dimuon \\phi;\\phi;Events ", 36, -3.6, 3.6);
   // h_dimuon_deta_onZ = new TH1F("dimuon_deta_onZ", "Dimuon deta;deta;Events ", 50, 0, 5.0);
   // h_dimuon_dphi_onZ = new TH1F("dimuon_dphi_onZ", "Dimuon dphi;dphi;Events ", 18, 0, 3.6);
   h_num_jets_onZ = new TH1F("num_jets_onZ", "Number of Jets;nJets;Events ", 8, 0, 8);
   h_num_bjets_onZ = new TH1F("num_bjets_onZ", "Number of B Jets;nBJets;Events ", 6, 0, 6);
   h_leadjet_pt_onZ = new TH1F("leadjet_pt_onZ", "Leading Jet p_{T};p_{T}  (GeV);Events ", 250, 0, 500);
   h_leadjet_eta_onZ = new TH1F("leadjet_eta_onZ", "Leading Jet \\eta;\\eta;Events ", 94, -4.7, 4.7);
   h_subjet_pt_onZ = new TH1F("subjet_pt_onZ", "Subleading Jet p_{T};p_{T}  (GeV);Events ", 250, 0, 500);
   // h_subjet_eta_onZ = new TH1F("subjet_eta_onZ", "Subleading Jet \\eta;\\eta;Events ", 94, -4.7, 4.7);
   h_dijet_mass_onZ = new TH1F("dijet_mass1_onZ", "DiJet Mass;M_{jj}  (GeV);Events ", 60, 0, 600);
   h_dijet_deta_onZ = new TH1F("dijet_deta1_onZ", "DiJet deta;deta;Events ", 94, 0, 9.4);
   h_met_pt_onZ = new TH1F("met_pt_onZ", "MET p_{T};p_{T}  (GeV) ", 100, 0, 200);
   h_met_phi_onZ = new TH1F("met_phi_onZ", "MET phi", 36, -3.6, 3.6);
   // h_zeppen_onZ = new TH1F("zeppen_onZ", "; zeppen; ", 100, -10, 10);
   // h_csTheta_onZ = new TH1F("csTheta_onZ", " ; csTheta; ", 10, -1, 1);
   // h_csPhi_onZ = new TH1F("csPhi_onZ", " ; csPhi ; ", 100, -10, 10);

   // h_leadMuon_pt_onH = new TH1F("lead_muon_pt_onH", "Leading Muon p_{T};p_{T}  (GeV);Events ", 100, 0, 200);
   h_leadMuon_eta_onH = new TH1F("lead_muon_eta_onH", "Leading Muon \\eta;\\eta;Events ", 50, -2.5, 2.5);
   // h_subMuon_pt_onH = new TH1F("sub_muon_pt_onH", "Subleading Muon p_{T};p_{T}  (GeV);Events ", 100, 0, 200);
   h_subMuon_eta_onH = new TH1F("sub_muon_eta_onH", "Subleading Muon eta;\\eta;Events ", 50, -2.5, 2.5);
   h_dimuon_mass_onH = new TH1F("dimuon_mass_onH", "Dimuon Mass;M_{\\mu \\mu}  (Gev);Events ", 80, 110, 150);
   h_dimuon_pt_onH = new TH1F("dimuon_pt_onH", "Dimuon p_{T};p_{T}  (GeV);Events ", 200, 0, 400);
   h_dimuon_eta_onH = new TH1F("dimuon_eta_onH", "Dimuon \\eta;\\eta;Events ", 100, -5.0, 5.0);
   // h_dimuon_phi_onH = new TH1F("dimuon_phi_onH", "Dimuon \\phi;\\phi;Events ", 36, -3.6, 3.6);
   // h_dimuon_deta_onH = new TH1F("dimuon_deta_onH", "Dimuon deta;deta;Events ", 50, 0, 5.0);
   // h_dimuon_dphi_onH = new TH1F("dimuon_dphi_onH", "Dimuon dphi;dphi;Events ", 18, 0, 3.6);
   h_num_jets_onH = new TH1F("num_jets_onH", "Number of Jets;nJets;Events ", 8, 0, 8);
   h_num_bjets_onH = new TH1F("num_bjets_onH", "Number of B Jets;nBJets;Events ", 6, 0, 6);
   h_leadjet_pt_onH = new TH1F("leadjet_pt_onH", "Leading Jet p_{T};p_{T}  (GeV);Events ", 250, 0, 500);
   h_leadjet_eta_onH = new TH1F("leadjet_eta_onH", "Leading Jet \\eta;\\eta;Events ", 94, -4.7, 4.7);
   h_subjet_pt_onH = new TH1F("subjet_pt_onH", "Subleading Jet p_{T};p_{T}  (GeV);Events ", 250, 0, 500);
   // h_subjet_eta_onH = new TH1F("subjet_eta_onH", "Subleading Jet \\eta;\\eta;Events ", 94, -4.7, 4.7);
   h_dijet_mass_onH = new TH1F("dijet_mass1_onH", "DiJet Mass;M_{jj}  (GeV);Events ", 60, 0, 600);
   h_dijet_deta_onH = new TH1F("dijet_deta1_onH", "DiJet deta;deta;Events ", 94, 0, 9.4);
   h_met_pt_onH = new TH1F("met_pt_onH", "MET p_{T};p_{T}  (GeV) ", 100, 0, 200);
   h_met_phi_onH = new TH1F("met_phi_onH", "MET phi", 36, -3.6, 3.6);

   // h_zeppen_onH = new TH1F("zeppen_onH", "; zeppen; ", 100, -10, 10);
   // h_csTheta_onH = new TH1F("csTheta_onH", " ; csTheta; ", 10, -1, 1);
   // h_csPhi_onH = new TH1F("csPhi_onH", " ; csPhi ; ", 100, -10, 10);

   h_bdtScore = new TH1F("h_bdtScore01jet", "; ; ", 20, -1, 1);
   h_dimuon_0 = new TH1F("dimuon_mass_cat0", " cat 0 ; dimuon mass (GeV); Events", 80, 110, 150);
   h_dimuon_1 = new TH1F("dimuon_mass_cat1", " cat 1 ; dimuon mass (GeV); Events", 80, 110, 150);
   h_dimuon_2 = new TH1F("dimuon_mass_cat2", " cat 2 ; dimuon mass (GeV); Events", 80, 110, 150);
   h_dimuon_3 = new TH1F("dimuon_mass_cat3", " cat 3 ; dimuon mass (GeV); Events", 80, 110, 150);
   h_dimuon_4 = new TH1F("dimuon_mass_cat4", " cat 4 ; dimuon mass (GeV); Events", 80, 110, 150);
   h_dimuon_5 = new TH1F("dimuon_mass_cat5", " cat 5 ; dimuon mass (GeV); Events", 80, 110, 150);
   h_dimuon_6 = new TH1F("dimuon_mass_cat6", " cat 6 ; dimuon mass (GeV); Events", 80, 110, 150);
   h_dimuon_7 = new TH1F("dimuon_mass_cat7", " cat 7 ; dimuon mass (GeV); Events", 80, 110, 150);
   h_dimuon_8 = new TH1F("dimuon_mass_cat8", " cat 8 ; dimuon mass (GeV); Events", 80, 110, 150);
   h_dimuon_9 = new TH1F("dimuon_mass_cat9", " cat 9 ; dimuon mass (GeV); Events", 80, 110, 150);

   // h_leadMuon_pt_onZ->Sumw2();
   h_leadMuon_eta_onZ->Sumw2();
   // h_subMuon_pt_onZ->Sumw2();
   h_subMuon_eta_onZ->Sumw2();
   h_dimuon_mass_onZ->Sumw2();
   h_dimuon_pt_onZ->Sumw2();
   h_dimuon_eta_onZ->Sumw2();
   // h_dimuon_phi_onZ->Sumw2();
   // h_dimuon_deta_onZ->Sumw2();
   // h_dimuon_dphi_onZ->Sumw2();
   h_num_jets_onZ->Sumw2();
   h_num_bjets_onZ->Sumw2();
   h_leadjet_pt_onZ->Sumw2();
   h_leadjet_eta_onZ->Sumw2();
   h_subjet_pt_onZ->Sumw2();
   // h_subjet_eta_onZ->Sumw2();
   h_dijet_mass_onZ->Sumw2();
   h_dijet_deta_onZ->Sumw2();
   h_met_pt_onZ->Sumw2();
   h_met_phi_onZ->Sumw2();
   // h_zeppen_onZ->Sumw2();
   // h_csTheta_onZ->Sumw2();
   // h_csPhi_onZ->Sumw2();

   // h_leadMuon_pt_onH->Sumw2();
   h_leadMuon_eta_onH->Sumw2();
   // h_subMuon_pt_onH->Sumw2();
   h_subMuon_eta_onH->Sumw2();
   h_dimuon_mass_onH->Sumw2();
   h_dimuon_pt_onH->Sumw2();
   h_dimuon_eta_onH->Sumw2();
   // h_dimuon_phi_onH->Sumw2();
   // h_dimuon_deta_onH->Sumw2();
   // h_dimuon_dphi_onH->Sumw2();
   h_num_jets_onH->Sumw2();
   h_num_bjets_onH->Sumw2();
   h_leadjet_pt_onH->Sumw2();
   h_leadjet_eta_onH->Sumw2();
   h_subjet_pt_onH->Sumw2();
   // h_subjet_eta_onH->Sumw2();
   h_dijet_mass_onH->Sumw2();
   h_dijet_deta_onH->Sumw2();
   h_met_pt_onH->Sumw2();
   h_met_phi_onH->Sumw2();
   // h_zeppen_onH->Sumw2();
   // h_csTheta_onH->Sumw2();
   // h_csPhi_onH->Sumw2();

   h_bdtScore->Sumw2();

   h_dimuon_0->Sumw2();
   h_dimuon_1->Sumw2();
   h_dimuon_2->Sumw2();
   h_dimuon_3->Sumw2();
   h_dimuon_4->Sumw2();
   h_dimuon_5->Sumw2();
   h_dimuon_6->Sumw2();
   h_dimuon_7->Sumw2();
   h_dimuon_8->Sumw2();
   h_dimuon_9->Sumw2();


   // GetOutputList()->Add(h_leadMuon_pt_onZ);
   GetOutputList()->Add(h_leadMuon_eta_onZ);
   // GetOutputList()->Add(h_subMuon_pt_onZ);
   GetOutputList()->Add(h_subMuon_eta_onZ);
   GetOutputList()->Add(h_dimuon_mass_onZ);
   GetOutputList()->Add(h_dimuon_pt_onZ);
   GetOutputList()->Add(h_dimuon_eta_onZ);
   // GetOutputList()->Add(h_dimuon_phi_onZ);
   // GetOutputList()->Add(h_dimuon_deta_onZ);
   // GetOutputList()->Add(h_dimuon_dphi_onZ);
   GetOutputList()->Add(h_num_jets_onZ);
   GetOutputList()->Add(h_num_bjets_onZ);
   GetOutputList()->Add(h_leadjet_pt_onZ);
   GetOutputList()->Add(h_leadjet_eta_onZ);
   GetOutputList()->Add(h_subjet_pt_onZ);
   // GetOutputList()->Add(h_subjet_eta_onZ);
   GetOutputList()->Add(h_dijet_mass_onZ);
   GetOutputList()->Add(h_dijet_deta_onZ);
   GetOutputList()->Add(h_met_pt_onZ);
   GetOutputList()->Add(h_met_phi_onZ);

   // GetOutputList()->Add(h_zeppen_onZ);
   // GetOutputList()->Add(h_csTheta_onZ);
   // GetOutputList()->Add(h_csPhi_onZ);

   // GetOutputList()->Add(h_leadMuon_pt_onH);
   GetOutputList()->Add(h_leadMuon_eta_onH);
   // GetOutputList()->Add(h_subMuon_pt_onH);
   GetOutputList()->Add(h_subMuon_eta_onH);
   GetOutputList()->Add(h_dimuon_mass_onH);
   GetOutputList()->Add(h_dimuon_pt_onH);
   GetOutputList()->Add(h_dimuon_eta_onH);
   // GetOutputList()->Add(h_dimuon_phi_onH);
   // GetOutputList()->Add(h_dimuon_deta_onH);
   // GetOutputList()->Add(h_dimuon_dphi_onH);
   GetOutputList()->Add(h_num_jets_onH);
   GetOutputList()->Add(h_num_bjets_onH);
   GetOutputList()->Add(h_leadjet_pt_onH);
   GetOutputList()->Add(h_leadjet_eta_onH);
   GetOutputList()->Add(h_subjet_pt_onH);
   // GetOutputList()->Add(h_subjet_eta_onH);
   GetOutputList()->Add(h_dijet_mass_onH);
   GetOutputList()->Add(h_dijet_deta_onH);
   GetOutputList()->Add(h_met_pt_onH);
   GetOutputList()->Add(h_met_phi_onH);
   // GetOutputList()->Add(h_zeppen_onH);
   // GetOutputList()->Add(h_csTheta_onH);
   // GetOutputList()->Add(h_csPhi_onH);
   GetOutputList()->Add(h_bdtScore);

   GetOutputList()->Add(h_dimuon_0);
   GetOutputList()->Add(h_dimuon_1);
   GetOutputList()->Add(h_dimuon_2);
   GetOutputList()->Add(h_dimuon_3);
   GetOutputList()->Add(h_dimuon_4);
   GetOutputList()->Add(h_dimuon_5);
   GetOutputList()->Add(h_dimuon_6);
   GetOutputList()->Add(h_dimuon_7);
   GetOutputList()->Add(h_dimuon_8);
   GetOutputList()->Add(h_dimuon_9);

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
   //   if (*h_mass < 110 || *h_mass > 150)


   float weight = 1.0;

   if (*h_mass > 70 && *h_mass < 110)
   {
      if (*mclabel)
      {
         weight = (*eWeight) * (*puWeight) * (*prefireSF) * (*idSF) * (*isoSF) * (*trigSF) * (*btagSF);
      }

      // h_leadMuon_pt_onZ->Fill(*muPtC_1, weight);
      h_leadMuon_eta_onZ->Fill(*muEtaC_1, weight);

      // h_subMuon_pt_onZ->Fill(*muPtC_2, weight);
      h_subMuon_eta_onZ->Fill(*muEtaC_2, weight);

      h_dimuon_mass_onZ->Fill(*h_mass, weight);
      h_dimuon_pt_onZ->Fill(*h_pt, weight);
      h_dimuon_eta_onZ->Fill(*h_eta, weight);
      // h_dimuon_phi_onZ->Fill(*h_phi, weight);
      // h_dimuon_deta_onZ->Fill(*h_deta, weight);
      // h_dimuon_dphi_onZ->Fill(*h_dphi, weight);

      h_num_jets_onZ->Fill(*njets, weight);
      h_num_bjets_onZ->Fill(*nbtagJets, weight);

      if (*njets > 0)
      {
         h_leadjet_pt_onZ->Fill(*jetpt_1, weight);
         h_leadjet_eta_onZ->Fill(*jeteta_1, weight);

         if (*njets > 1)
         {
            h_subjet_pt_onZ->Fill(*jetpt_2, weight);
            // h_subjet_eta_onZ->Fill(*jeteta_2, weight);

            h_dijet_mass_onZ->Fill(*mjj, weight);
            h_dijet_deta_onZ->Fill(*detajj, weight);
         }
      }
      h_met_pt_onZ->Fill(*metpt, weight);
      h_met_phi_onZ->Fill(*metphi, weight);

      // h_zeppen_onZ->Fill(*zeppen, weight);
      // h_csTheta_onZ->Fill(*csTheta, weight);
      // h_csPhi_onZ->Fill(*csPhi, weight);
   }

   if (*h_mass > 110 && *h_mass < 150)
   {
      if (*mclabel)
      {
         weight = (*eWeight) * (*zptWeight) * (*puWeight) * (*prefireSF) * (*idSF) * (*isoSF) * (*trigSF) * (*btagSF);
      }

      // h_leadMuon_pt_onH->Fill(*muPtC_1, weight);
      h_leadMuon_eta_onH->Fill(*muEtaC_1, weight);

      // h_subMuon_pt_onH->Fill(*muPtC_2, weight);
      h_subMuon_eta_onH->Fill(*muEtaC_2, weight);

      h_dimuon_mass_onH->Fill(*h_mass, weight);
      h_dimuon_pt_onH->Fill(*h_pt, weight);
      h_dimuon_eta_onH->Fill(*h_eta, weight);
      // h_dimuon_phi_onH->Fill(*h_phi, weight);
      // h_dimuon_deta_onH->Fill(*h_deta, weight);
      // h_dimuon_dphi_onH->Fill(*h_dphi, weight);

      h_num_jets_onH->Fill(*njets, weight);
      h_num_bjets_onH->Fill(*nbtagJets, weight);

      if (*njets > 0)
      {
         h_leadjet_pt_onH->Fill(*jetpt_1, weight);
         h_leadjet_eta_onH->Fill(*jeteta_1, weight);

         if (*njets > 1)
         {
            h_subjet_pt_onH->Fill(*jetpt_2, weight);
            // h_subjet_eta_onH->Fill(*jeteta_2, weight);

            h_dijet_mass_onH->Fill(*mjj, weight);
            h_dijet_deta_onH->Fill(*detajj, weight);
         }
      }
      h_met_pt_onH->Fill(*metpt, weight);
      h_met_phi_onH->Fill(*metphi, weight);

      // h_zeppen_onH->Fill(*zeppen, weight);
      // h_csTheta_onH->Fill(*csTheta, weight);
      // h_csPhi_onH->Fill(*csPhi, weight);

      if (*category == 0)
      {
         h_bdtScore->Fill(*bdtScore, weight);
         h_dimuon_0->Fill(*h_mass, weight);
      }
      if (*category == 1)
      {
         h_bdtScore->Fill(*bdtScore, weight);
         h_dimuon_1->Fill(*h_mass, weight);
      }
      if (*category == 2)
      {
         h_bdtScore->Fill(*bdtScore, weight);
         h_dimuon_2->Fill(*h_mass, weight);
      }
      if (*category == 3)
      {
         h_bdtScore->Fill(*bdtScore, weight);
         h_dimuon_3->Fill(*h_mass, weight);
      }
      if (*category == 4)
      {
         h_bdtScore->Fill(*bdtScore, weight);
         h_dimuon_4->Fill(*h_mass, weight);
      }
      if (*category == 5)
      {
         h_dimuon_5->Fill(*h_mass, weight);
      }
      if (*category == 6)
      {
         h_dimuon_6->Fill(*h_mass, weight);
      }
      if (*category == 7)
      {
         h_dimuon_7->Fill(*h_mass, weight);
      }
      if (*category == 8)
      {
         h_dimuon_8->Fill(*h_mass, weight);
      }
      if (*category == 9)
      {
         h_dimuon_9->Fill(*h_mass, weight);
      }
      // if (*category == 10)
      // {
      //    h_dimuon_10->Fill(*h_mass, weight);
      // }
      // if (*category == 11)
      // {
      //    h_dimuon_11->Fill(*h_mass, weight);
      // }
      // if (*category == 12)
      // {
      //    h_dimuon_12->Fill(*h_mass, weight);
      // }
      // if (*category == 13)
      // {
      //    h_dimuon_13->Fill(*h_mass, weight);
      // }
      // if (*category == 14)
      // {
      //    h_dimuon_14->Fill(*h_mass, weight);
      // }
   }
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
