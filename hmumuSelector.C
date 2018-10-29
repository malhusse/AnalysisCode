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
   //   outputTree = new TTree("processedNtuples","processedNtuplesTree");
    h_muon_pt = new TH1F("muon_pt", "Muon pT", 500,0,500);
    h_muon_corrPT = new TH1F("muon_corrPT", "Muon corrected pT", 500,0,500);
    h_leadMuon_pt = new TH1F("lead_muon_pt", "Leading Muon pT", 500,0,500);
    h_subMuon_pt = new TH1F("sub_muon_pt", "Subleading Muon pT", 500, 0,500);

    GetOutputList()->Add(h_muon_pt);
    GetOutputList()->Add(h_muon_corrPT);
    GetOutputList()->Add(h_leadMuon_pt);
    GetOutputList()->Add(h_subMuon_pt);

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
    for (int iMuon = 0, nMuons = Muons__charge.GetSize(); iMuon < nMuons; ++iMuon){
       h_muon_pt->Fill(Muons__pt[iMuon]);
       h_muon_corrPT->Fill(Muons__corrPT[iMuon]);
    }


    int lead_muon_id = 0;
    int sub_muon_id = 1;

    if (Muons__corrPT[1] > Muons__corrPT[0]){
        lead_muon_id = 1;
        sub_muon_id = 0;
    }

    h_leadMuon_pt->Fill(Muons__corrPT[lead_muon_id]);
    h_subMuon_pt->Fill(Muons__corrPT[sub_muon_id]);

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

