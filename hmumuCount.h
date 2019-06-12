//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Nov 11 02:13:09 2018 by ROOT version 6.10/09
// from TChain ntuplemaker_H2DiMuonMaker/Meta/
//////////////////////////////////////////////////////////

#ifndef hmumuCount_h
#define hmumuCount_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "interface/MetaHiggs.h"



class hmumuCount : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain
   
   Int_t year = 0;
   TH1 *hc_numEventsWeighted = 0;
   TH1 *hc_numEvents = 0;

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> _nEventsProcessed = {fReader, "_nEventsProcessed"};
   TTreeReaderValue<Int_t> _sumEventWeights = {fReader, "_sumEventWeights"};

   hmumuCount(TTree * /*tree*/ =0) { }
   virtual ~hmumuCount() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(hmumuCount,0);

};

#endif

#ifdef hmumuCount_cxx
void hmumuCount::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t hmumuCount::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef hmumuCount_cxx
