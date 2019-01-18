void quickPUMC(string inputList, string outputFile)
{
  cout << "input filelist: "<< inputList << " , output filename: " << outputFile << endl;
  cout << "loading count chain" << endl;

  TChain *eventsChain = new TChain("ntuplemaker_H2DiMuonMaker/Events");
  std::ifstream inputAgain(inputList);
  for (std::string line; std::getline(inputAgain,line); ){
    eventsChain->Add(TString(line.c_str()));
  }

  TProof::Open("workers=8");
  gProof->AddInput(new TNamed("outputName",outputFile));
  gProof->Exec("gSystem->Load(\"/uscms/home/malhusse/nobackup/build/libAnalysisCore.so\")");

  eventsChain->SetProof();
  eventsChain->Process("genPUMC.C+");
}
