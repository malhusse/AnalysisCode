void quickCount(string inputList, string outputFile, Int_t year)
{

  gSystem->Load("/uscms_data/d1/malhusse/analysis/libAnalysisCore.so");
  cout << "input filelist: "<< inputList << " , output filename: " << outputFile << endl;
  cout << "loading count chain" << endl;

  std::string fileList = "filelists/";
  fileList += std::to_string(year);
  fileList += "/";
  fileList += inputList;
  fileList += ".files";

  std::cout << fileList << std::endl;

  TChain *metaChain = new TChain("ntuplemaker_H2DiMuonMaker/Meta");
  std::ifstream input(fileList);
  for (std::string line; std::getline(input,line); ){
    metaChain->Add(TString(line.c_str()));
  }

  std::cout << metaChain->GetEntries() << std::endl;

  TProof::Open("workers=4");
  gProof->AddInput(new TNamed("outputName",outputFile));
  gProof->Exec("gSystem->Load(\"/uscms_data/d1/malhusse/analysis/libAnalysisCore.so\")");
  gProof->SetParameter("getYear", year);

  metaChain->SetProof();
  metaChain->Process("hmumuCount.C+");

}
