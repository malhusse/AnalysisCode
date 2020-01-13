// ntupleMode = 0 to create histograms
// ntupldeMode = 1 to create TNtuple (to be used with TMVA)
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

void quickScript(string inputList, string outputFile, Int_t year)
{
   std::ifstream jsonFile("resources/cross_sections.json");
   std::ifstream labelFile("resources/mc_labels.json");
   // boost::property_tree::ptree _xsecJson;
   boost::property_tree::ptree _labelJson;

   // boost::property_tree::json_parser::read_json(jsonFile, _xsecJson);
   boost::property_tree::json_parser::read_json(labelFile, _labelJson);

   typedef boost::property_tree::ptree::path_type path;

   std::string inputDataset = std::to_string(year);
   inputDataset += ".";
   inputDataset += inputList;

   std::string fileList = "filelists/";
   fileList += std::to_string(year);
   fileList += "/";
   fileList += inputList;
   fileList += ".files";

   std::cout << "input filelist: " << fileList << " , output filename: " << outputFile << endl;
   // std::cout << _xsecJson.get<double>(inputDataset) << std::endl;
   std::cout << _labelJson.get<int>(inputDataset) << std::endl;
   // double xsec = _xsecJson.get<double>(inputDataset);
   int mcLabel = _labelJson.get<int>(inputDataset);

   gSystem->Load("/uscms_data/d1/malhusse/analysis/libAnalysisCore.so");
   TProof::Open("workers=4");
   gProof->AddInput(new TNamed("outputName", outputFile));
   gProof->Exec("gSystem->Load(\"/uscms_data/d1/malhusse/analysis/libAnalysisCore.so\")");
   //   gProof->Exec("gSystem->Load(\"/uscms_data/d1/malhusse/analysis/libAnalysisAuxTools.so\")");

   Int_t sumEvents = 1;
   Int_t sumEventsWeighted = 1;
   if (mcLabel)
   {
      TString countFile = "resources/count/";
      countFile += std::to_string(year);
      countFile += "/";
      countFile += outputFile;

      TFile *f = TFile::Open(countFile);
      TH1I *hsumEvents = (TH1I *)f->Get("cnumEvents");
      TH1I *hsumEventsWeighted = (TH1I *)f->Get("cnumEventsWeighted");
      sumEvents = hsumEvents->GetMaximum();
      sumEventsWeighted = hsumEventsWeighted->GetMaximum();
   }

   cout << sumEvents << endl;
   cout << sumEventsWeighted << endl;
   // exit(1);

   gProof->SetParameter("getSumEvents", sumEvents);
   gProof->SetParameter("getSumEventsWeighted", sumEventsWeighted);
   // gProof->SetParameter("getxsec", xsec);
   gProof->SetParameter("getmcLabel", mcLabel);
   gProof->SetParameter("getYear", year);

   TChain *eventsChain = new TChain("ntuplemaker_H2DiMuonMaker/Events");
   std::ifstream inputAgain(fileList);
   for (std::string line; std::getline(inputAgain, line);)
   {
      eventsChain->Add(TString(line.c_str()));
   }

   eventsChain->SetProof();
   eventsChain->Process("hmumuSelector.C+");
}
