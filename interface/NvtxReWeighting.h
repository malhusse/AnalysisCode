#include <boost/shared_ptr.hpp>
#include <cmath>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1.h"

class NvtxReWeighting
{

public:
  NvtxReWeighting(std::string generatedFile,
                  std::string dataFile,
                  std::string GenHistName,
                  std::string DataHistName);

  NvtxReWeighting(){};

  double weight(int nvertices);
  boost::shared_ptr<TH1> weights_;
protected:
  std::string generatedFileName_;
  std::string dataFileName_;
  std::string GenHistName_;
  std::string DataHistName_;
  boost::shared_ptr<TFile> generatedFile_;
  boost::shared_ptr<TFile> dataFile_;
  //  boost::shared_ptr<TH1> weights_;

  // keep copies of normalized distributions:
  boost::shared_ptr<TH1> MC_distr_;
  boost::shared_ptr<TH1> Data_distr_;
};

NvtxReWeighting::NvtxReWeighting(std::string generatedFile, std::string dataFile, std::string GenHistName, std::string DataHistName)
    : generatedFileName_(generatedFile), dataFileName_(dataFile), GenHistName_(GenHistName), DataHistName_(DataHistName)
{
  generatedFile_ = boost::shared_ptr<TFile>(new TFile(generatedFileName_.c_str()));
  dataFile_ = boost::shared_ptr<TFile>(new TFile(dataFileName_.c_str()));

  Data_distr_ = boost::shared_ptr<TH1>((static_cast<TH1 *>(dataFile_->Get(DataHistName_.c_str())->Clone())));
  MC_distr_ = boost::shared_ptr<TH1>((static_cast<TH1 *>(generatedFile_->Get(GenHistName_.c_str())->Clone())));

  // normalize both histograms first //
  Data_distr_->Scale(1.0 / Data_distr_->Integral());
  MC_distr_->Scale(1.0 / MC_distr_->Integral());

  weights_ = boost::shared_ptr<TH1>(static_cast<TH1 *>(Data_distr_->Clone()));
  weights_->SetName("NvtxWeights");

  TH1 *den = dynamic_cast<TH1 *>(MC_distr_->Clone());
  weights_->Divide(den);
  int NBins = weights_->GetNbinsX();

  // compute avarage Nvtx re-weight //
  double SumMC = 0;
  double SumNvtxReWeightTimesMC = 0;
  for (int ibin = 0; ibin < NBins; ++ibin)
  {
    double MC = den->GetBinContent(ibin);
    double NvtxReWeight = weights_->GetBinContent(ibin);
    SumMC += MC;
    SumNvtxReWeightTimesMC += NvtxReWeight * MC;
  }
  double SumWeight = SumNvtxReWeightTimesMC / SumMC;
  weights_->Scale(1 / SumWeight);
  ///////////////////////////////

  double GetWeight = 0;
  for (int ibin = 1; ibin < NBins + 1; ++ibin)
    GetWeight += weights_->GetBinContent(ibin) * den->GetBinContent(ibin);
}

double NvtxReWeighting::weight(int nvertices)
{
  int bin = weights_->GetXaxis()->FindBin(nvertices);
  return weights_->GetBinContent(bin);
}