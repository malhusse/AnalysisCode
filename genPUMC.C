#define genPUMC_cxx
// The class definition in genPUMC.h has been generated automatically
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
// root> T->Process("genPUMC.C")
// root> T->Process("genPUMC.C","some options")
// root> T->Process("genPUMC.C+")
//


#include "genPUMC.h"
#include <TH2.h>
#include <TStyle.h>

TString _outputName;

void genPUMC::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   TNamed *name = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));
   TParameter<Int_t> *pYear = dynamic_cast<TParameter<Int_t>*>(fInput->FindObject("getYear"));
   year = pYear->GetVal();

   _outputName = "data/pileup/";
   _outputName += std::to_string(year);
   _outputName += "/";
   _outputName += name ? name->GetTitle() : "outputName.root";  
}

void genPUMC::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   h_pileup = new TH1D("pileup","pileup",99,0,99);
   GetOutputList()->Add(h_pileup);

}

Bool_t genPUMC::Process(Long64_t entry)
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
   h_pileup->Fill(*_nPU,*_genWeight);

   return kTRUE;
}

void genPUMC::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void genPUMC::Terminate()
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
