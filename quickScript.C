void quickScript(string inputList, string outputFile)
{
   cout << "input filelist: " << inputList << " , output filename: " << outputFile << endl;

   TProof::Open("workers=4");
   gProof->AddInput(new TNamed("outputName", outputFile));
   gProof->Exec("gSystem->Load(\"/uscms/home/malhusse/nobackup/build/libAnalysisCore.so\")");
   gProof->Exec("gSystem->Load(\"/uscms/home/malhusse/nobackup/build/libAnalysisAuxTools.so\")");
   TString countFile = "count/";
   countFile += outputFile;

   TFile *f = TFile::Open(countFile);
   TH1I *hsumEvents = (TH1I *)f->Get("cnumEvents");
   TH1I *hsumEventsWeighted = (TH1I *)f->Get("cnumEventsWeighted");
   Int_t sumEvents = hsumEvents->GetMaximum();
   Int_t sumEventsWeighted = hsumEventsWeighted->GetMaximum();

   cout << sumEvents << endl;
   cout << sumEventsWeighted << endl;

   gProof->SetParameter("getSumEvents", sumEvents);
   gProof->SetParameter("getSumEventsWeighted", sumEventsWeighted);

   TChain *eventsChain = new TChain("ntuplemaker_H2DiMuonMaker/Events");
   std::ifstream inputAgain(inputList);
   for (std::string line; std::getline(inputAgain, line);)
   {
      eventsChain->Add(TString(line.c_str()));
   }

   eventsChain->SetProof();
   eventsChain->Process("hmumuSelector.C+");
}
