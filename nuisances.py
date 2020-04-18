import ROOT as R
import os
import sys
import json
import pandas as pd
from math import trunc

os.environ["ANALYSISHOME"] = "/Users/mo/hep/Analysis/HMuMu"
sys.path.append(os.path.join(
    os.environ["ANALYSISHOME"], "Configuration", "higgs"))
import Samples as S

def truncate(f, n):
    return trunc(f * 10 ** n) / 10 ** n

def get_systematics(year):
    systs = ["iso", "id", "trig", "pileup", "prefire", "jetAbsolute",
             "jetBBEC1", "jetEC2", "jetFlavorQCD", "jetHF", "jetRelativeBal"]
    systs += ["jetAbsolute_{}".format(year), "jetBBEC1_{}".format(year), "jetEC2_{}".format(
        year), "jetHF_{}".format(year), "jetRelativeSample_{}".format(year)]
    return systs


def calculate_nuisances(cat, era):

    # I dont need a stack, all I need is for each year to get the sample at 120, 125, 130 histogram for each cateogry
    # cat is catX
    # histo is of form dimuon_mass_cat0_noSyst
    # this is the data histogram for a specific year and
    # specific variable (dimuon_mass_catX_noSyst)

    rootFile.cd()

    nominal_histo_name = "dimuon_mass_{}_noSyst".format(cat)
    hdata = dataTFile.Get(nominal_histo_name).Clone()
    hdata.SetName("data_{}_{}".format(cat, era))
    hdata.Write()

    for v in root_dic:  # This loops over the samples...
        # print(v)
        central_mass = False
        if v[1]["isSignal"]:
            nominal_histo = v[1]["TFile"].Get(nominal_histo_name).Clone()
            nominal_histo.Scale(lumi[v[1]["year"]] * v[1]["xsec"])

            key = ""
            if "ttH" in v[0]:
                key += "ttH_"
            elif "WminusH" in v[0]:
                key += "WHminus_"
            elif "WplusH" in v[0]:
                key += "WHplus_"
            elif "ZH" in v[0]:
                key += "ZH_"
            elif "GluGluH" in v[0]:
                key += "ggH_"
            elif "VBFH" in v[0]:
                key += "qqH_"

            if "M120" in v[0]:
                key += "120_"
            elif "M125" in v[0]:
                key += "125_"
                central_mass = True
            elif "M130" in v[0]:
                key += "130_"

            key += "{}_{}".format(cat, era)

            nominal_histo.SetName(key)
            nominal_histo.Write()

            if (central_mass):
                for sys in systematics:
                    if "jet" in sys:
                        sys_name = sys + "_"
                    else:
                        sys_name = sys
                    histo_up_name = "dimuon_mass_{}_{}Up".format(cat, sys_name)
                    histo_up = v[1]["TFile"].Get(histo_up_name).Clone()
                    histo_up.Scale(lumi[v[1]["year"]] * v[1]["xsec"])

                    histo_down_name = "dimuon_mass_{}_{}Down".format(cat, sys_name)
                    histo_down = v[1]["TFile"].Get(histo_down_name).Clone()
                    histo_down.Scale(lumi[v[1]["year"]] * v[1]["xsec"])

                    kappa_down = histo_down.Integral() / nominal_histo.Integral()
                    kappa_up = histo_up.Integral() / nominal_histo.Integral()

                    process = key.split("_")[0]
                    nuisance_dataframe[process][sys] = "{}/{}".format(truncate(kappa_down,3),truncate(kappa_up,3))
                # this is a M125 sample, use it to calculate the up/down fluctuations....


if __name__ == "__main__":
    R.gROOT.SetBatch(R.kTRUE)
    years = ["2016", "2017", "2018"]

    if len(sys.argv) > 1:
        modifier = sys.argv[1]
    else:
        print("need modifier!")
        exit(1)

    # if not os.path.isdir("plots_{}/{}".format(modifier,year)):
        # os.makedirs("plots_{}/{}".format(modifier,year))
    # print(datafile)

    mc2016 = {}
    # mc2016.update(S.mc_background_2016)
    mc2016.update(S.mc_signal_2016)
    mc2016.update(S.mc_signal_2016_extra)
    # mc2016.update(S.mc_background_2016_extra)
    # mc2016.update(S.mc_background_2016_extra_2)

    mc2017 = {}
    # mc2017.update(S.mc_background_2017)
    mc2017.update(S.mc_signal_2017)
    mc2017.update(S.mc_signal_2017_extra)
    # mc2017.update(S.mc_background_2017_extra)
    # mc2017.update(S.mc_background_2017_extra_2)

    mc2018 = {}
    # mc2018.update(S.mc_background_2018)
    mc2018.update(S.mc_signal_2018)
    mc2018.update(S.mc_signal_2018_extra)
    # mc2018.update(S.mc_background_2018_extra)
    # mc2018.update(S.mc_background_2018_extra_2)

    mc_datasets_dic = {
        "2016": mc2016,
        "2017": mc2017,
        "2018": mc2018
    }

    lumi = {
        "2016": 35922,  # 35.922
        "2017": 41529,  # 41.529
        "2018": 59740  # 59.74
    }

    with open("resources/run2xsecs.json") as json_file:
        xsecDic = json.load(json_file)

    for year in years:

        datafile = "/Users/mo/hep/Analysis/AnalysisCode/histoFiles_{}/{}/allData{}.root".format(
            modifier, year, year)

        mc_samples = mc_datasets_dic[year]
        root_dic = []
        # lumi = lumi[year]

        # bcColors = [40, 30, 41, 42, 43, 35, 46, 47, 38, 28, 29]
        # sigColors = [2, 3, 4, 6, 8, 9]
        for file in os.listdir("/Users/mo/hep/Analysis/AnalysisCode/histoFiles_{}/{}/".format(modifier, year)):
            # if "105To160" in file:
            #     continue
            if file.endswith(".root") and not file.startswith("all"):
                nickname = file.replace(".root", "")
                xsec = xsecDic[nickname]
                nickname = nickname.replace("__", "/")
                dic = {}

                for k, v in mc_samples.items():
                    if nickname in k:
                        # print(nickname)
                        dic["isSignal"] = v.isSignal
                        dic["year"] = str(v.year)

                        f = R.TFile.Open(
                            "/Users/mo/hep/Analysis/AnalysisCode/histoFiles_{}/{}/".format(modifier, year)+file)
                        dic["TFile"] = f
                        dic["xsec"] = xsec
                        root_dic.append((nickname, dic))

        systematics = get_systematics(year)
        channels = ["ggH", "qqH", "ttH", "WHminus", "WHplus", "ZH"]

        
        
        # root_dic_sorted = sorted(root_dic, key=lambda el: el[1]['wEvents'])

        dataTFile = R.TFile.Open(datafile)
        categories = ["cat0", "cat1", "cat2", "cat3", "cat4", "cat5", "cat6", "cat7", "cat8", "cat9"]

        # if "2019" in year:
        # outfile = R.TFile.Open("FitHistosAll.root", "RECREATE")
        rootFile = R.TFile.Open("FitHistos_{}.root".format(year), "RECREATE")

        for cat in categories:
            # print(variable)
            # calculate_unc(variable)
            nuisance_dataframe = pd.DataFrame(columns=channels, index=systematics)
            calculate_nuisances(cat, year)
            nuisance_dataframe.to_csv("nuisanceCSV/nuisances_{}_{}.csv".format(cat, year))
            

        rootFile.Close()

    # plot_muon_corr()
