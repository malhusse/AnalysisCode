#include "histosMap.h"

class histogramCollection {
    
    private:
    std::map <std::string, histosMap*> histosCollection_;

    public:

    histogramCollection() {
        histosCollection_["noSyst"] = new histosMap("noSyst");
        histosCollection_["isoUp"] = new histosMap("isoUp");
        histosCollection_["isoDown"] = new histosMap("isoDown");

    }

    void fillHistograms(eventEntry dataToFill, std::string systName, float totalWeight) {
        if (! histosCollection_.count(systName)) {
            // For some reason creating new maps dynamically is not working at the moment..
            // histosCollection_[systName] = new histosMap(systName);
            return;
        }

        histosCollection_[systName]->fillEntry(dataToFill, totalWeight);
    }

    std::vector<TH1F*> getHistograms() {
        std::vector<TH1F*> histosVec;
        
        for ( std::map< std::string, histosMap*>::iterator m = histosCollection_.begin() ; m != histosCollection_.end(); ++m)
        {
            std::vector<TH1F*> temp = m->second->getHistograms();
            histosVec.insert( histosVec.end(),  temp.begin(), temp.end());
        }
        return histosVec;
    }


};