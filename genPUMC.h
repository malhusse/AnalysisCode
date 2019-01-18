//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 15 11:48:37 2019 by ROOT version 6.10/09
// from TChain ntuplemaker_H2DiMuonMaker/Events/
//////////////////////////////////////////////////////////

#ifndef genPUMC_h
#define genPUMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <vector>

#include "interface/Muon.h"
#include "interface/Jet.h"
#include "interface/Vertex.h"
#include "interface/Event.h"
#include "interface/MET.h"
//#include "interface/MetaHiggs.h"
#include "interface/Electron.h"


class genPUMC : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain


   TH1 *h_pileup = 0;
   
   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> _nPU = {fReader, "_nPU"};
   // TTreeReaderValue<Int_t> _nvtx = {fReader, "_nvtx"};
   TTreeReaderValue<Int_t> _genWeight = {fReader, "_genWeight"};
   // TTreeReaderValue<Bool_t> _passedMetFilters = {fReader, "_passedMetFilters"};
   // TTreeReaderValue<vector<bool>> _hasHLTFired = {fReader, "_hasHLTFired"};
   // TTreeReaderArray<pair<string,int>> _metFilterBits = {fReader, "_metFilterBits"};


   genPUMC(TTree * /*tree*/ =0) { }
   virtual ~genPUMC() { }
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

   ClassDef(genPUMC,0);

};

#endif

#ifdef genPUMC_cxx
void genPUMC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t genPUMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef genPUMC_cxx
