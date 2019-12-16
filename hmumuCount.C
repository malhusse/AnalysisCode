#define hmumuCount_cxx
// The class definition in hmumuCount.h has been generated automatically
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
// root> T->Process("hmumuCount.C")
// root> T->Process("hmumuCount.C","some options")
// root> T->Process("hmumuCount.C+")
//


#include "hmumuCount.h"
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TList.h>

TString _outputName;

void hmumuCount::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   TNamed *name = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));
   TParameter<Int_t> *pYear = dynamic_cast<TParameter<Int_t>*>(fInput->FindObject("getYear"));
   year = pYear->GetVal();

   _outputName = "resources/count/";
   _outputName += std::to_string(year);
   _outputName += "/";
   _outputName += name ? name->GetTitle() : "outputName.root";  
}

void hmumuCount::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   hc_numEventsWeighted = new TH1I("cnumEventsWeighted","Weighted numEvents Proccessed", 2, 0, 2);
   hc_numEvents = new TH1I("cnumEvents","numEvents Processed",2,0,2);
   GetOutputList()->Add(hc_numEventsWeighted);
   GetOutputList()->Add(hc_numEvents);
}

Bool_t hmumuCount::Process(Long64_t entry)
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
   hc_numEvents->Fill(1,*_nEventsProcessed);
   hc_numEventsWeighted->Fill(1,*_sumEventWeights);
   return kTRUE;
}

void hmumuCount::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void hmumuCount::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  TList *output_list = (TList*)GetOutputList();
  TFile fout(_outputName,"recreate");
  for(const auto&& obj: *output_list){
    if(obj->IsA()->InheritsFrom("TH1"))
      obj->Write();
  }
  fout.Close();
}
