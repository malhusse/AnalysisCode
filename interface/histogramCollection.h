#include "histosMap.h"

class histogramCollection {
    
    private:
    std::map <std::string, histosMap*> histosCollection_;

    public:

    histogramCollection(int mcLabel) {
        histosCollection_["noSyst"] = new histosMap("noSyst");
        if (mcLabel) {
            histosCollection_["isoUp"] = new histosMap("isoUp");
            histosCollection_["isoDown"] = new histosMap("isoDown");
            histosCollection_["idUp"] = new histosMap("idUp");
            histosCollection_["idDown"] = new histosMap("idDown");
            histosCollection_["trigUp"] = new histosMap("trigUp");
            histosCollection_["trigDown"] = new histosMap("trigDown");
            histosCollection_["pileupUp"] = new histosMap("pileupUp");
            histosCollection_["pileupDown"] = new histosMap("pileupDown");
            histosCollection_["prefireUp"] = new histosMap("prefireUp");
            histosCollection_["prefireDown"] = new histosMap("prefireDown");
            histosCollection_["jetTotal_Up"] = new histosMap("jetTotal_Up");
            histosCollection_["jetAbsolute_Up"] = new histosMap("jetAbsolute_Up");
            histosCollection_["jetBBEC1_Up"] = new histosMap("jetBBEC1_Up");
            histosCollection_["jetEC2_Up"] = new histosMap("jetEC2_Up");
            histosCollection_["jetFlavorQCD_Up"] = new histosMap("jetFlavorQCD_Up");
            histosCollection_["jetHF_Up"] = new histosMap("jetHF_Up");
            histosCollection_["jetRelativeBal_Up"] = new histosMap("jetRelativeBal_Up");
            histosCollection_["jetAbsolute_2016_Up"] = new histosMap("jetAbsolute_2016_Up");
            histosCollection_["jetAbsolute_2017_Up"] = new histosMap("jetAbsolute_2017_Up");
            histosCollection_["jetAbsolute_2018_Up"] = new histosMap("jetAbsolute_2018_Up");
            histosCollection_["jetBBEC1_2016_Up"] = new histosMap("jetBBEC1_2016_Up");
            histosCollection_["jetBBEC1_2017_Up"] = new histosMap("jetBBEC1_2017_Up");
            histosCollection_["jetBBEC1_2018_Up"] = new histosMap("jetBBEC1_2018_Up");
            histosCollection_["jetEC2_2016_Up"] = new histosMap("jetEC2_2016_Up");
            histosCollection_["jetEC2_2017_Up"] = new histosMap("jetEC2_2017_Up");
            histosCollection_["jetEC2_2018_Up"] = new histosMap("jetEC2_2018_Up");
            histosCollection_["jetHF_2016_Up"] = new histosMap("jetHF_2016_Up");
            histosCollection_["jetHF_2017_Up"] = new histosMap("jetHF_2017_Up");
            histosCollection_["jetHF_2018_Up"] = new histosMap("jetHF_2018_Up");
            histosCollection_["jetRelativeSample_2016_Up"] = new histosMap("jetRelativeSample_2016_Up");
            histosCollection_["jetRelativeSample_2017_Up"] = new histosMap("jetRelativeSample_2017_Up");
            histosCollection_["jetRelativeSample_2018_Up"] = new histosMap("jetRelativeSample_2018_Up");
            histosCollection_["jetTotal_Down"] = new histosMap("jetTotal_Down");
            histosCollection_["jetAbsolute_Down"] = new histosMap("jetAbsolute_Down");
            histosCollection_["jetBBEC1_Down"] = new histosMap("jetBBEC1_Down");
            histosCollection_["jetEC2_Down"] = new histosMap("jetEC2_Down");
            histosCollection_["jetFlavorQCD_Down"] = new histosMap("jetFlavorQCD_Down");
            histosCollection_["jetHF_Down"] = new histosMap("jetHF_Down");
            histosCollection_["jetRelativeBal_Down"] = new histosMap("jetRelativeBal_Down");
            histosCollection_["jetAbsolute_2016_Down"] = new histosMap("jetAbsolute_2016_Down");
            histosCollection_["jetAbsolute_2017_Down"] = new histosMap("jetAbsolute_2017_Down");
            histosCollection_["jetAbsolute_2018_Down"] = new histosMap("jetAbsolute_2018_Down");
            histosCollection_["jetBBEC1_2016_Down"] = new histosMap("jetBBEC1_2016_Down");
            histosCollection_["jetBBEC1_2017_Down"] = new histosMap("jetBBEC1_2017_Down");
            histosCollection_["jetBBEC1_2018_Down"] = new histosMap("jetBBEC1_2018_Down");
            histosCollection_["jetEC2_2016_Down"] = new histosMap("jetEC2_2016_Down");
            histosCollection_["jetEC2_2017_Down"] = new histosMap("jetEC2_2017_Down");
            histosCollection_["jetEC2_2018_Down"] = new histosMap("jetEC2_2018_Down");
            histosCollection_["jetHF_2016_Down"] = new histosMap("jetHF_2016_Down");
            histosCollection_["jetHF_2017_Down"] = new histosMap("jetHF_2017_Down");
            histosCollection_["jetHF_2018_Down"] = new histosMap("jetHF_2018_Down");
            histosCollection_["jetRelativeSample_2016_Down"] = new histosMap("jetRelativeSample_2016_Down");
            histosCollection_["jetRelativeSample_2017_Down"] = new histosMap("jetRelativeSample_2017_Down");
            histosCollection_["jetRelativeSample_2018_Down"] = new histosMap("jetRelativeSample_2018_Down");
        }
    }

    void fillHistograms(eventEntry dataToFill, std::string systName, float totalWeight) {
        if (systName == "noSyst") histosCollection_[systName]->fillEntry(dataToFill, totalWeight);
        else if (systName.find("jet") != std::string::npos) histosCollection_[systName]->fillJetSystEntry(dataToFill, totalWeight);
        else histosCollection_[systName]->fillEventSystEntry(dataToFill, totalWeight);; 

        
    }

    std::vector<TH1F*> getAllHistograms() {
        std::vector<TH1F*> histosVec;
        
        for ( std::map< std::string, histosMap*>::iterator m = histosCollection_.begin() ; m != histosCollection_.end(); ++m)
        {
            std::vector<TH1F*> temp = m->second->getHistograms();
            histosVec.insert( histosVec.end(),  temp.begin(), temp.end());
        }
        return histosVec;
    }

    TH1F* getHistogram(std::string systName, std::string histoName) {
        return histosCollection_[systName]->getHistogram(histoName);
    }

    std::vector<std::string> getSystematics() {
        std::vector<std::string> systVec;

        for ( std::map< std::string, histosMap*>::iterator m = histosCollection_.begin() ; m != histosCollection_.end(); ++m)
        {
            systVec.push_back(m->first);
        }

        return systVec;
    }

};