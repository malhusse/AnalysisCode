// ntupleMode = 0 to create histograms
// ntupldeMode = 1 to create TNtuple (to be used with TMVA)
void quickScript(string inputList, string outputFile, Int_t ntupleMode)
{
   std::map<std::string, Float_t> xsec =
    {{"Run2017B-31Mar2018-v1.root", 1.0},
     {"Run2017C-31Mar2018-v1.root", 1.0},
     {"Run2017D-31Mar2018-v1.root", 1.0},
     {"Run2017E-31Mar2018-v1.root", 1.0},
     {"Run2017F-31Mar2018-v1.root", 1.0},
     {"GluGluHToMuMu_M125_13TeV_amcatnloFXFX_pythia8.root", 0.01057},
     {"VBFHToMuMu_M125_13TeV_amcatnlo_pythia8.root", .0008230},
     {"WminusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8.root", 0.0001160},
     {"WplusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8.root", 0.0001852},
     {"ZH_HToMuMu_ZToAll_M125_13TeV_powheg_pythia8.root", 0.0001923},
     {"ttHToMuMu_M125_TuneCP5_13TeV-powheg-pythia8.root", 0.00011034},
     {"DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root", 5765.0},
     {"ZZ_TuneCP5_13TeV-pythia8.root", 16.523},
     {"WWTo2L2Nu_NNPDF31_TuneCP5Up_13TeV-powheg-pythia8.root", 12.46},
     {"WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8.root", 45.99},
     {"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root", 11.61},
     {"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root", 5.595},
     {"WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8.root", 4.42965},
     {"WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8.root", 0.2086},
     {"WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8.root", 0.1651},
     {"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root", 85.656},
     {"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root", 687.1}};
 
 std::map<std::string, Int_t> mc =
    {{"Run2017B-31Mar2018-v1.root", 0},
     {"Run2017C-31Mar2018-v1.root", 0},
     {"Run2017D-31Mar2018-v1.root", 0},
     {"Run2017E-31Mar2018-v1.root", 0},
     {"Run2017F-31Mar2018-v1.root", 0},
     {"GluGluHToMuMu_M125_13TeV_amcatnloFXFX_pythia8.root", -10},
     {"VBFHToMuMu_M125_13TeV_amcatnlo_pythia8.root", -20},
     {"WminusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8.root", -40},
     {"WplusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8.root", -50},
     {"ZH_HToMuMu_ZToAll_M125_13TeV_powheg_pythia8.root", -30},
     {"ttHToMuMu_M125_TuneCP5_13TeV-powheg-pythia8.root", -60},
     {"DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root", 10},
     {"ZZ_TuneCP5_13TeV-pythia8.root", 32},
     {"WWTo2L2Nu_NNPDF31_TuneCP5Up_13TeV-powheg-pythia8.root", 30},
     {"WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8.root", 30},
     {"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root", 31},
     {"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root", 31},
     {"WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8.root", 31},
     {"WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8.root",33 },
     {"WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8.root", 34},
     {"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root", 20},
     {"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root", 20}};


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
   cout << "ntuple mode? " <<  ntupleMode << endl;
   cout << "xsec: " << xsec[outputFile] << endl;
   cout << "mcLabel: " << mc[outputFile] << endl;

   gProof->SetParameter("getNtupleMode", ntupleMode);
   gProof->SetParameter("getSumEvents", sumEvents);
   gProof->SetParameter("getSumEventsWeighted", sumEventsWeighted);
   gProof->SetParameter("getxsec", xsec[outputFile]);
   gProof->SetParameter("getmcLabel", mc[outputFile]);

   TChain *eventsChain = new TChain("ntuplemaker_H2DiMuonMaker/Events");
   std::ifstream inputAgain(inputList);
   for (std::string line; std::getline(inputAgain, line);)
   {
      eventsChain->Add(TString(line.c_str()));
   }

   eventsChain->SetProof();
   eventsChain->Process("hmumuSelector.C+");
}
