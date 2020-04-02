#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

std::map<std::string, JetCorrectionUncertainty *> getUncertainties(int year)
{
    const int nsrc = 12;
    std::map<std::string, JetCorrectionUncertainty *> vsrc;
    if (year == 2016)
    {
        const char *srcnames[nsrc] = {"Absolute", "Absolute_2016", "BBEC1", "BBEC1_2016", "EC2", "EC2_2016", "FlavorQCD", "HF", "HF_2016", "RelativeBal", "RelativeSample_2016", "Total"};
        for (int isrc = 0; isrc < nsrc; isrc++)
        {

            const char *name = srcnames[isrc];
            JetCorrectorParameters *p = new JetCorrectorParameters("/uscms/home/malhusse/nobackup/analysis/AnalysisCode/resources/jetunc/Regrouped_2016_AK4PFchs.txt", name);
            JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
            vsrc[name] = unc;
        }
    }

    else if (year == 2017)
    {
        const char *srcnames[nsrc] = {"Absolute", "Absolute_2017", "BBEC1", "BBEC1_2017", "EC2", "EC2_2017", "FlavorQCD", "HF", "HF_2017", "RelativeBal", "RelativeSample_2017", "Total"};
        for (int isrc = 0; isrc < nsrc; isrc++)
        {
            const char *name = srcnames[isrc];
            JetCorrectorParameters *p = new JetCorrectorParameters("/uscms/home/malhusse/nobackup/analysis/AnalysisCode/resources/jetunc/Regrouped_2017_AK4PFchs.txt", name);
            JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
            vsrc[name] = unc;
        }
    }

    else if (year == 2018)
    {
        const char *srcnames[nsrc] = {"Absolute", "Absolute_2018", "BBEC1", "BBEC1_2018", "EC2", "EC2_2018", "FlavorQCD", "HF", "HF_2018", "RelativeBal", "RelativeSample_2018", "Total"};
        for (int isrc = 0; isrc < nsrc; isrc++)
        {

            const char *name = srcnames[isrc];
            JetCorrectorParameters *p = new JetCorrectorParameters("/uscms/home/malhusse/nobackup/analysis/AnalysisCode/resources/jetunc/Regrouped_2018_AK4PFchs.txt", name);
            JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
            vsrc[name] = unc;
        }
    }

    return vsrc;
}