#define genZPT_cxx
// The class definition in genZPT.h has been generated automatically
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
// root> T->Process("genZPT.C")
// root> T->Process("genZPT.C","some options")
// root> T->Process("genZPT.C+")
//

#include "genZPT.h"
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TList.h>

double const PDG_MASS_Mu = 0.1056583745;

TString _outputName;

void genZPT::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   TNamed *name = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));
   TParameter<Int_t> *pYear = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getYear"));
   year = pYear->GetVal();

   _outputName = "data/zpt/";
   _outputName += std::to_string(year);
   _outputName += "/";
   _outputName += name ? name->GetTitle() : "outputName.root";
}

void genZPT::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   h_zpt = new TH1D("zpt", "zpt", 100, 0, 1000);
   GetOutputList()->Add(h_zpt);
}

Bool_t genZPT::Process(Long64_t entry)
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
   if (Muons.GetSize() < 2)
      return kFALSE;
   if (Muons[0]._charge * Muons[1]._charge != -1)
      return kFALSE;
   if ((!Muons[0]._isMedium || Muons[0]._pfIso > 0.25) || (!Muons[1]._isMedium || Muons[1]._pfIso > 0.25))
      return kFALSE;
   TLorentzVector p4m1, p4m2;
   p4m1.SetPtEtaPhiM(Muons[0]._pt, Muons[0]._eta, Muons[0]._phi, PDG_MASS_Mu);
   p4m2.SetPtEtaPhiM(Muons[1]._pt, Muons[1]._eta, Muons[1]._phi, PDG_MASS_Mu);
   if ((p4m1 + p4m2).M() < 70 || (p4m1 + p4m2).M() > 110)
      return kFALSE;

   h_zpt->Fill((p4m1 + p4m2).Pt(), *_genWeight);

   return kTRUE;
}

void genZPT::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void genZPT::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   TList *output_list = (TList *)GetOutputList();
   TFile fout(_outputName, "recreate");
   for (const auto &&obj : *output_list)
   {
      if (obj->IsA()->InheritsFrom("TH1"))
         obj->Write();
   }
   fout.Close();
}
