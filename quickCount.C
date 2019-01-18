void quickCount(string inputList, string outputFile)
{

  cout << "input filelist: "<< inputList << " , output filename: " << outputFile << endl;
  cout << "loading count chain" << endl;

  TChain *metaChain = new TChain("ntuplemaker_H2DiMuonMaker/Meta");
  std::ifstream input(inputList);
  for (std::string line; std::getline(input,line); ){
    metaChain->Add(TString(line.c_str()));
  }


  TProof::Open("workers=8");
  gProof->AddInput(new TNamed("outputName",outputFile));
  gProof->Exec("gSystem->Load(\"/uscms/home/malhusse/nobackup/build/libAnalysisCore.so\")");

  metaChain->SetProof();
  metaChain->Process("hmumuCount.C+");

}
