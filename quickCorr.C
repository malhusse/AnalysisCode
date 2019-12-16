void quickCorr(string inputList, string outputFile, Int_t year)
{
  cout << "input filelist: "<< inputList << " , output filename: " << outputFile << endl;
  cout << "loading count chain" << endl;

  gSystem->Load("/uscms_data/d1/malhusse/analysis/libAnalysisCore.so");

  std::string fileList = "filelists/";
  fileList += std::to_string(year);
  fileList += "/";
  fileList += inputList;
  fileList += ".files";

  TChain *eventsChain = new TChain("ntuplemaker_H2DiMuonMaker/Events");
  std::ifstream inputAgain(fileList);
  for (std::string line; std::getline(inputAgain,line); ){
    eventsChain->Add(TString(line.c_str()));
  }

  TProof::Open("workers=4");
  gProof->AddInput(new TNamed("outputName",outputFile));
  gProof->Exec("gSystem->Load(\"/uscms_data/d1/malhusse/analysis/libAnalysisCore.so\")");
  gProof->SetParameter("getYear", year);

  eventsChain->SetProof();
  eventsChain->Process("genCorr.C+");
}
