#define genCorr_cxx
// The class definition in genCorr.h has been generated automatically
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
// root> T->Process("genCorr.C")
// root> T->Process("genCorr.C","some options")
// root> T->Process("genCorr.C+")
//

#include "genCorr.h"
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TList.h>

double const PDG_MASS_Mu = 0.1056583745;
TString _outputName;

void genCorr::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   TNamed *name = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));

   TParameter<Int_t> *pYear = dynamic_cast<TParameter<Int_t> *>(fInput->FindObject("getYear"));
   year = pYear->GetVal();

   _outputName = "resources/corrections/";
   _outputName += std::to_string(year);
   _outputName += "/";
   _outputName += name->GetTitle();

   std::cout << "Correction Histograms" << std::endl;
}

void genCorr::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   TNamed *nameSlave = dynamic_cast<TNamed *>(fInput->FindObject("outputName"));
   TString _sname = nameSlave->GetTitle();

   if (((string)_sname.Data()).find("SingleMuon") != string::npos) //is Data
      isData = true;

   // if (((string)_sname.Data()).find("DY") != string::npos) //is Drell-Yan Sample
      // isDY = true;

   // if (((string)_sname.Data()).find("GluGlu") != string::npos)
      // isGGH = true; // for nnlops

   if (!isData)
   {
      h_pileup = new TH1D("pileup", "pileup", 99, 0, 99);
      h_pileup->Sumw2();
      GetOutputList()->Add(h_pileup);
   }
   // if (isDY || isData)
   // {
   //    h_zpt_0j = new TH1D("zpt_0j", "zpt 0j", 100, 0, 1000);
   //    h_zpt_1j = new TH1D("zpt_1j", "zpt 1j", 100, 0, 1000);
   //    h_zpt_2j = new TH1D("zpt_2j", "zpt 2j", 100, 0, 1000);
   //    h_zpt_0j->Sumw2();
   //    h_zpt_1j->Sumw2();
   //    h_zpt_2j->Sumw2();
   //    GetOutputList()->Add(h_zpt_0j);
   //    GetOutputList()->Add(h_zpt_1j);
   //    GetOutputList()->Add(h_zpt_2j);
   // }

}

Bool_t genCorr::Process(Long64_t entry)
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


   if (!isData)
   {
      h_pileup->Fill(*_nPU, *_genWeight);
   }

   // double dimupt = -10;

   // for (unsigned int i = 1; i < Muons.GetSize(); ++i)
   // {
   //    if (Muons[i]._charge != Muons[0]._charge )
   //    {
   //       TLorentzVector p4m1, p4m2;
   //       p4m1.SetPtEtaPhiM(Muons[0]._pt, Muons[0]._eta, Muons[0]._phi, PDG_MASS_Mu);
   //       p4m2.SetPtEtaPhiM(Muons[i]._pt, Muons[i]._eta, Muons[i]._phi, PDG_MASS_Mu);
         
   //       if ((p4m1 + p4m2).M() > 70 && (p4m1 + p4m2).M() < 110)
   //       {
   //          dimupt = (p4m1 + p4m2).Pt();
   //       }

   //    }
   // }


   // if (dimupt > 0 && (isDY || isData))
   // {
   //    if (Jets.GetSize() == 0)
   //       h_zpt_0j->Fill(dimupt, *_genWeight);
   //    else if (Jets.GetSize() == 1)
   //       h_zpt_1j->Fill(dimupt, *_genWeight);
   //    else if (Jets.GetSize() >= 2)
   //       h_zpt_2j->Fill(dimupt, *_genWeight);
   // }

   return kTRUE;
}

void genCorr::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void genCorr::Terminate()
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
