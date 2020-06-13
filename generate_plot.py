import os, sys
import json

os.environ["ANALYSISHOME"] = "/Users/mo/hep/Analysis/HMuMu"
sys.path.append(os.path.join(
    os.environ["ANALYSISHOME"], "Configuration", "higgs"))

import Samples as S
import ROOT as R
from array import array
from math import sqrt
import numpy as np
from variable_dics import *


runBlind = False

doUncertainties = True
doStat = True
doJES = True
doJER = True
doTheory = True
doMuon = True

def calculate_unc(variable, nominalRatio, rebinFactor = 0):

    hdata = dataTFile.Get(variable).Clone()

    if rebinFactor:
        hdata.Rebin(rebinFactor)

    combined_systematics_dic = {}

    # If I split the uncertainties into theory / stat / expertimen then I would neeed 3 graphs

    x = array('d')
    y = array('d')
    exl = array('d')
    exh = array('d')
    nbins = hdata.GetNbinsX()


    for i in range(1, nbins + 1):
        x.append( hdata.GetBinLowEdge(i)+.5*hdata.GetBinWidth(i) )
        y.append(1.0)
        exl.append(.5*hdata.GetBinWidth(i))
        exh.append(.5*hdata.GetBinWidth(i))

    if (doTheory):
        bkgStack_up = R.THStack("bkgStackUp", "")
        bkgStack_down = R.THStack("bkgStackDown","")

        for v in root_dic:
            if "DY" in v[0] or "EWK" in v[0]:
                shift_up = 1.05
                shift_down = .95
            else:
                shift_up = 1.07
                shift_down = 0.93

            histoUp = v[1]["TFile"].Get(variable).Clone()
            histoDown = v[1]["TFile"].Get(variable).Clone()
            histoUp.Scale(lumi[v[1]["year"]] * v[1]["xsec"] * shift_up)
            histoDown.Scale(lumi[v[1]["year"]] * v[1]["xsec"] * shift_down)
            if rebinFactor:
                histoUp.Rebin(rebinFactor)
                histoDown.Rebin(rebinFactor)
            bkgStack_up.Add(histoUp)
            bkgStack_down.Add(histoDown)

        stackUp = hdata.Clone()
        stackUp.Divide(bkgStack_up.GetStack().Last())
        stackDown = hdata.Clone()
        stackDown.Divide(bkgStack_down.GetStack().Last())

        combined_systematics_up = []
        combined_systematics_down = []

        for i in range(1,nbins+1):
            if stackUp.GetBinContent(i):
                if "dimuon_mass" in variable and stackUp.GetBinCenter(i) > 120 and stackUp.GetBinCenter(i) < 130 and runBlind:
                    error_up = 0
                    error_down = 0

                else:
                    error_up = abs(stackUp.GetBinContent(i) - nominalRatio.GetBinContent(i))
                    error_down = abs(stackDown.GetBinContent(i) - nominalRatio.GetBinContent(i))
            else:
                error_up = 0
                error_down = 0

            combined_systematics_up.append(error_up)
            combined_systematics_down.append(error_down)

        combined_systematics_dic["theory"] = (np.asarray(combined_systematics_up), np.asarray(combined_systematics_down))
    
    if (doMuon):
        bkgStack_up = R.THStack("bkgStackUp", "")
        bkgStack_down = R.THStack("bkgStackDown","")

        for v in root_dic:
            shift_up = 1.01
            shift_down = 0.99

            histoUp = v[1]["TFile"].Get(variable).Clone()
            histoDown = v[1]["TFile"].Get(variable).Clone()
            histoUp.Scale(lumi[v[1]["year"]] * v[1]["xsec"] * shift_up)
            histoDown.Scale(lumi[v[1]["year"]] * v[1]["xsec"] * shift_down)
            if rebinFactor:
                histoUp.Rebin(rebinFactor)
                histoDown.Rebin(rebinFactor)
            bkgStack_up.Add(histoUp)
            bkgStack_down.Add(histoDown)

        stackUp = hdata.Clone()
        stackUp.Divide(bkgStack_up.GetStack().Last())
        stackDown = hdata.Clone()
        stackDown.Divide(bkgStack_down.GetStack().Last())

        combined_systematics_up = []
        combined_systematics_down = []

        for i in range(1,nbins+1):
            if stackUp.GetBinContent(i):
                if "dimuon_mass" in variable and stackUp.GetBinCenter(i) > 120 and stackUp.GetBinCenter(i) < 130 and runBlind:
                    error_up = 0
                    error_down = 0

                else:
                    error_up = abs(stackUp.GetBinContent(i) - nominalRatio.GetBinContent(i))
                    error_down = abs(stackDown.GetBinContent(i) - nominalRatio.GetBinContent(i))
            else:
                error_up = 0
                error_down = 0

            combined_systematics_up.append(error_up)
            combined_systematics_down.append(error_down)

        combined_systematics_dic["muon"] = (np.asarray(combined_systematics_up), np.asarray(combined_systematics_down))

    if (doStat):

        combined_systematics_up = []
        combined_systematics_down = []
        
        for i in range(1,nbins+1):
            if nominalRatio.GetBinContent(i):
                if "dimuon_mass" in variable and nominalRatio.GetBinCenter(i) > 120 and nominalRatio.GetBinCenter(i) < 130 and runBlind:
                    error = 0
                else:
                    error = abs(nominalRatio.GetBinError(i))
            else:
                error = 0

            combined_systematics_up.append(error)
            combined_systematics_down.append(error)

        combined_systematics_dic["stat"] = (np.asarray(combined_systematics_up), np.asarray(combined_systematics_down))
        
    if (doJES):
        bkgStack_up = R.THStack("bkgStackUp", "")
        bkgStack_down = R.THStack("bkgStackDown","")

        varUp = variable.split("noSyst")[0] + "jetTotal_Up"
        varDown = variable.split("noSyst")[0] + "jetTotal_Down"
    
        # print(variable, varUp, varDown)
        for v in root_dic:

            histoUp = v[1]["TFile"].Get(varUp).Clone()
            histoDown = v[1]["TFile"].Get(varDown).Clone()
            histoUp.Scale(lumi[v[1]["year"]] * v[1]["xsec"])
            histoDown.Scale(lumi[v[1]["year"]] * v[1]["xsec"])
            if rebinFactor:
                histoUp.Rebin(rebinFactor)
                histoDown.Rebin(rebinFactor)
            bkgStack_up.Add(histoUp)
            bkgStack_down.Add(histoDown)

        stackUp = hdata.Clone()
        stackUp.Divide(bkgStack_up.GetStack().Last())
        stackDown = hdata.Clone()
        stackDown.Divide(bkgStack_down.GetStack().Last())

        combined_systematics_up = []
        combined_systematics_down = []

        for i in range(1,nbins+1):
            if stackUp.GetBinContent(i):
                if "dimuon_mass" in variable and stackUp.GetBinCenter(i) > 120 and stackUp.GetBinCenter(i) < 130 and runBlind:
                    error_up = 0
                    error_down = 0

                else:
                    error_up = abs(stackUp.GetBinContent(i) - nominalRatio.GetBinContent(i))
                    error_down = abs(stackDown.GetBinContent(i) - nominalRatio.GetBinContent(i))
            
            else:
                error_up = 0
                error_down = 0

            combined_systematics_up.append(error_up)
            combined_systematics_down.append(error_down)

        combined_systematics_dic["jetTotal"] = (np.asarray(combined_systematics_up), np.asarray(combined_systematics_down))
    
    if (doJER):
        bkgStackJER = R.THStack("bkgStackJER", "")

        varJER = variable.split("noSyst")[0] + "jetResolution"
    
        # print(variable, varUp, varDown)
        for v in root_dic:

            histoJER = v[1]["TFile"].Get(varJER).Clone()
            histoJER.Scale(lumi[v[1]["year"]] * v[1]["xsec"])
            if rebinFactor:
                histoJER.Rebin(rebinFactor)

            bkgStackJER.Add(histoJER)

        stackJER = hdata.Clone()
        stackJER.Divide(bkgStackJER.GetStack().Last())

        combined_systematics_up = []
        combined_systematics_down = []

        for i in range(1,nbins+1):
            if stackJER.GetBinContent(i):
                if "dimuon_mass" in variable and stackJER.GetBinCenter(i) > 120 and stackJER.GetBinCenter(i) < 130 and runBlind:
                    error_up = 0
                    error_down = 0

                else:
                    if stackJER.GetBinContent(i) < 1:
                        error_up = abs(stackJER.GetBinContent(i) - nominalRatio.GetBinContent(i))
                        error_down = 0
                    else:
                        error_up = 0
                        error_down = abs(stackJER.GetBinContent(i) - nominalRatio.GetBinContent(i))

            else:
                error_up = 0
                error_down = 0

            combined_systematics_up.append(error_up)
            combined_systematics_down.append(error_down)

        combined_systematics_dic["jetResolution"] = (np.asarray(combined_systematics_up), np.asarray(combined_systematics_down))

    total_up = np.zeros(nbins)
    total_down = np.zeros(nbins)
    for _, combined_systematics in combined_systematics_dic.items():
        total_up += combined_systematics[0]**2
        total_down += combined_systematics[1]**2

    g = R.TGraphAsymmErrors(nbins, x, y, exl, exh, np.sqrt(total_up), np.sqrt(total_down))
    g.SetLineColor(0)
    g.SetMarkerColor(0)
    g.SetMarkerSize(0)
    g.SetFillColor(R.kOrange)

    return g

CMS_Text = R.TPaveText(0.13,0.91,0.18,.96,"NB NDC")
CMS_Text.SetTextSize(0.04)
CMS_Text.SetFillColor(0)
CMS_Text.AddText("#font[62]{CMS}") #font[12]{work in progress}")

def plot_variable(variable, modifier, logy=True):
    canvas = R.TCanvas("c1", "c1", 400,700)
    canvas.SetCanvasSize(500, 600)
    canvas.cd()
    pad1 = R.TPad("pad1", "pad1", 0, 0.24, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.SetTicks(1,1)
    pad1.Draw()
    canvas.cd()
    pad2 = R.TPad("pad2", "pad2", 0, 0, 1, 0.24)
    pad2.Draw()
    pad2.cd()
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.3)
    pad2.SetGridy()
    pad1.cd()

    bkgStack = R.THStack("bkgStack", "")
    ttStack = R.THStack("ttStack", "")
    dyStack = R.THStack("dyStack", "")
    # ewkStack = R.THStack("ewkStack", "")
    diBosonStack = R.THStack("diBosonStack", "")
    otherStack = R.THStack("otherStack","")

    sigStack_125 = R.THStack("sigStack_125", "")
    gghStack_125 = R.THStack("gghStack_125","ggH")
    vbfStack_125 = R.THStack("vbfStack_125","VBF")
    vAndttHStack_125 = R.THStack("vAndttHStack_125","VH and ttH")
    # tthStack_125 = R.THStack("tthStack_125","ttH")
    # wmHStack_125 = R.THStack("wmHStack_125","wmH")
    # wpHStack_125 = R.THStack("wpHStack_125","wpH")
    # zHStack_125 = R.THStack("zHStack_125","ZH")

    hdata = dataTFile.Get(variable).Clone()
    hdata.SetBinErrorOption(R.TH1.kPoisson)

    leg = R.TLegend(0.45, 0.70, .88, .88)
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    lumis = float(lumi[year])/1000
    lumi_Text = R.TPaveText(0.65,0.91,0.85,.96,"NB NDC")
    lumi_Text.SetTextSize(0.04)
    lumi_Text.SetFillColor(0)
    lumi_Text.AddText("#font[62]{{ {:.1f} fb^{{-1}} (13 TeV)}}".format(lumis))

    if "2019" in year:
        leg.AddEntry(hdata, "Run II Data")
    else:
        leg.AddEntry(hdata, "Data")

    for v in root_dic:
        if v[1]["isSignal"]:
            histo = v[1]["TFile"].Get(variable).Clone()
            histo.Scale(lumi[v[1]["year"]] * v[1]["xsec"])

            if "M120" in v[0] or "M130" in v[0]:
                continue

            elif "M125" in v[0]:
                if "ttH" in v[0]:
                    vAndttHStack_125.Add(histo)
                elif "WminusH" in v[0]:
                    vAndttHStack_125.Add(histo)
                elif "WplusH" in v[0]:
                    vAndttHStack_125.Add(histo)
                elif "ZH" in v[0]:
                    vAndttHStack_125.Add(histo)
                elif "GluGluH" in v[0]:
                    gghStack_125.Add(histo)
                elif "VBFH" in v[0]:
                    vbfStack_125.Add(histo)
                sigStack_125.Add(histo)
                # allStack.Add(histo)

        if not v[1]["isSignal"]:
            histo = v[1]["TFile"].Get(variable).Clone()
            histo.Scale(lumi[v[1]["year"]] * v[1]["xsec"])
            if "ST" in v[0]:
                ttStack.Add(histo)
            elif "TTW" in v[0] or "TTZ" in v[0]:
                ttStack.Add(histo)
            elif "TT" in v[0]:
                ttStack.Add(histo)
            elif "WWW" in v[0] or "WWZ" in v[0] or "ZZZ" in v[0]:
                otherStack.Add(histo)
            elif "ZZ" in v[0] or "WZ" in v[0] or "WW" in v[0]:
                diBosonStack.Add(histo)
            elif "DY" in v[0]:
                if "onZ" in variable and "M-105" in v[0]:
                    continue
                elif "M-50" in v[0] and ("onH" in variable or "onGGH" in variable or "cat" in variable):
                    continue
                dyStack.Add(histo)

            elif "EWK" in v[0]:
                otherStack.Add(histo)
            elif "tZq" in v[0]:
                otherStack.Add(histo)
            elif "GluGluToContin" in v[0]:
                otherStack.Add(histo)

            else:
                print(v[0])
                print("This sample is not in any stack??")

    # make this None, then get the .Last() using a
    # try: get() except: pass (keep moving)
    # then later on do if is not None

    hist_ttStack = None
    hist_dyStack = None
    # hist_ewkStack = None
    hist_diBosonStack = None
    hist_otherStack = None

    # hist_tthStack_125 = tthStack_125.GetStack().Last()
    # hist_tthStack_125.SetName("ttH_125_{}".format(variable))
    hist_vbfStack_125 = vbfStack_125.GetStack().Last()
    hist_vbfStack_125.SetName("vbf_125_{}".format(variable))
    # hist_wmHStack_125 = wmHStack_125.GetStack().Last()
    # hist_wmHStack_125.SetName("wmH_125_{}".format(variable))
    # hist_wpHStack_125 = wpHStack_125.GetStack().Last()
    # hist_wpHStack_125.SetName("wpH_125_{}".format(variable))
    hist_gghStack_125 = gghStack_125.GetStack().Last()
    hist_gghStack_125.SetName("ggH_125_{}".format(variable))
    # hist_zHStack_125 = zHStack_125.GetStack().Last()
    # hist_zHStack_125.SetName("zH_125_{}".format(variable))
    hist_vAndttHStack_125 = vAndttHStack_125.GetStack().Last()
    hist_vAndttHStack_125.SetName("vAndttH_125_{}".format(variable))

        
    try:
        hist_ttStack = ttStack.GetStack().Last()
        leg.AddEntry(hist_ttStack, "Top", "F")
    except:
        pass
  
    try:
        hist_dyStack = dyStack.GetStack().Last()
        leg.AddEntry(hist_dyStack, "DY", "F")

    except:
        pass
    
    # try:
    #     hist_ewkStack = ewkStack.GetStack().Last()
    #     leg.AddEntry(hist_ewkStack, "EWK", "F")

    # except:
    #     pass

    try:
        hist_diBosonStack = diBosonStack.GetStack().Last()
        leg.AddEntry(hist_diBosonStack, "Diboson", "F")

    except:
        pass
  
    try:
        hist_otherStack = otherStack.GetStack().Last()
        leg.AddEntry(hist_otherStack, "Other Bkg.", "F")
    except:
        pass

    stackList = [ hist_otherStack, hist_diBosonStack, hist_ttStack, hist_dyStack]

    for i in range(0, len(stackList)):
        if stackList[i] is not None:
            stackList[i].SetFillColor(bcColors[(i % len(bcColors))])
            stackList[i].SetLineColor(R.kBlack)
            if variable in variable_rebin:
                stackList[i].Rebin(variable_rebin[variable])
            bkgStack.Add(stackList[i])

    # blind data in 120-130 GeV bins on mass plots
    if "dimuon_mass" in variable and runBlind:
        for i in range(hdata.GetNbinsX()):
            if hdata.GetBinCenter(i + 1) > 120 and hdata.GetBinCenter(i + 1) < 130:
                hdata.SetBinContent(i + 1, 0)
                hdata.SetBinError(i + 1, 0)

    if variable in variable_rebin:
        hdata.Rebin(variable_rebin[variable])

    hdata.SetMarkerStyle(20)
    hdata.SetMarkerSize(.7)

    if logy:
        pad1.SetLogy()
        bkgStack.SetMinimum(.1)

    bkgStack.SetMaximum(hdata.GetMaximum() * 1000)
    hdata.SetStats(R.kFALSE)
    hdata.GetXaxis().SetLabelSize(0)
    bkgStack.Draw("hist")

    j = 0
    sigStackList = [hist_vAndttHStack_125, hist_vbfStack_125, hist_gghStack_125]
    for h in sigStackList:
            h.SetLineColor(sigColors[(j % len(sigColors))])
            h.SetLineWidth(3)

            if variable in variable_rebin:
                h.Rebin(variable_rebin[variable])
            h.Draw("hist same")
            j+=1

    leg.AddEntry(hist_vbfStack_125,"VBF","l")
    leg.AddEntry(hist_vAndttHStack_125, "VH and ttH", "l")
    leg.AddEntry(hist_gghStack_125,"ggH","l")

    hdata.Draw("SAME P")
    bkgStack.GetYaxis().SetTitle("Events")
    leg.Draw()
    CMS_Text.Draw()
    lumi_Text.Draw()

    R.gPad.Modified()

    pad2.cd()
    hratio = hdata.Clone()
    hratio.SetTitle("")
    hratio.GetYaxis().SetTitle("Data / Pred")
    hratio.GetXaxis().SetTitle(variable_label[variable])
    hratio.GetYaxis().SetRangeUser(0.55,1.45)
    hratio.GetYaxis().SetNdivisions(505)
    hratio.GetYaxis().SetTitleSize(15)
    hratio.GetYaxis().SetTitleFont(43)
    hratio.GetYaxis().SetTitleOffset(1.5)
    hratio.GetYaxis().SetLabelFont(43)
    hratio.GetYaxis().SetLabelSize(12)
    hratio.GetXaxis().SetTitleSize(15)
    hratio.GetXaxis().SetTitleFont(43)
    hratio.GetXaxis().SetTitleOffset(4.0)
    hratio.GetXaxis().SetLabelFont(43)
    hratio.GetXaxis().SetLabelOffset(.007)
    hratio.GetXaxis().SetLabelSize(13)
    hratio.Divide(bkgStack.GetStack().Last())
    hratio.SetStats(R.kFALSE)
    # hratio.SetMaximum(1.45)
    # hratio.SetMinimum(0.55)
    hratio.SetMarkerStyle(20)
    hratio.SetMarkerSize(0.7)
    # hratio.Draw("E P")
    
    if doUncertainties:
        rebinFactor = 0 if variable not in variable_rebin else variable_rebin[variable]
        g = calculate_unc(variable, hratio.Clone(), rebinFactor)
        g.SetTitle("")
        g.GetYaxis().SetTitle("Data / Pred")
        g.GetXaxis().SetTitle(variable_label[variable])
        g.GetYaxis().SetNdivisions(4, R.kFALSE)
        g.GetYaxis().SetTitleSize(15)
        g.GetYaxis().SetTitleFont(43)
        g.GetYaxis().SetTitleOffset(1.5)
        g.GetYaxis().SetLabelFont(43)
        g.GetYaxis().SetLabelSize(12)
        g.GetXaxis().SetTitleSize(15)
        g.GetXaxis().SetTitleFont(43)
        g.GetXaxis().SetTitleOffset(4.0)
        g.GetXaxis().SetLabelFont(43)
        g.GetXaxis().SetLabelOffset(.007)
        g.GetXaxis().SetLabelSize(13)
        g.GetYaxis().SetRangeUser(0.55,1.45)
        g.GetYaxis().SetNdivisions(505)

        # g.SetMaximum(1.45)
        # g.SetMinimum(0.55)
        g.GetXaxis().SetLimits(hratio.GetBinLowEdge(1), hratio.GetBinLowEdge(hratio.GetNbinsX())+hratio.GetBinWidth(1))
        g.Draw("a2")
        # g.Draw("p")

    unhist = hratio.Clone("unhist")
    for i in range(0, unhist.GetNbinsX()+1):
        unhist.SetBinContent(i,1)
        unhist.SetBinError(i,0)

    unhist.SetMarkerSize(0)
    unhist.SetLineWidth(2)
    unhist.SetLineColor(R.kRed)
    unhist.SetLineStyle(7)
    unhist.SetFillColor(0)

    unhist.Draw("hist ][ same")
    hratio.Draw("P E0 same")
    
    R.gPad.Modified()
    canvas.Draw()
    canvas.SaveAs("plots_{}/{}/".format(modifier, year) +
                  variable + "_{}.pdf".format(year))
    
    del canvas, pad1, pad2, hratio, bkgStack, hdata


if __name__ == "__main__":
    R.gROOT.SetBatch(R.kTRUE)
    years = ["2016", "2017", "2018"]
    years = ["2019"]

    if len(sys.argv) > 1:
        modifier = sys.argv[1]

    else:
        print("need modifier!")
        exit(1)

    for year in years:
        if not os.path.isdir("plots_{}/{}".format(modifier,year)):
            os.makedirs("plots_{}/{}".format(modifier,year))
        datafile = "/Users/mo/hep/Analysis/AnalysisCode/histoFiles_{}/{}/allData{}.root".format(
            modifier, year, year)

        mc2016 = S.mc_background_2016
        mc2016.update(S.mc_signal_2016)
        mc2016.update(S.mc_signal_2016_extra)
        mc2016.update(S.mc_background_2016_extra)
        mc2016.update(S.mc_background_2016_extra_2)


        mc2017 = S.mc_background_2017
        mc2017.update(S.mc_signal_2017)
        mc2017.update(S.mc_signal_2017_extra)
        mc2017.update(S.mc_background_2017_extra)
        mc2017.update(S.mc_background_2017_extra_2)


        mc2018 = S.mc_background_2018
        mc2018.update(S.mc_signal_2018)
        mc2018.update(S.mc_signal_2018_extra)
        mc2018.update(S.mc_background_2018_extra)
        mc2018.update(S.mc_background_2018_extra_2)

        mc2019 = {}
        mc2019.update(mc2016)
        mc2019.update(mc2017)
        mc2019.update(mc2018)

        mc_datasets_dic = {
            "2016": mc2016,
            "2017": mc2017,
            "2018": mc2018,
            "2019": mc2019
        }

        lumi = {
            "2016": 35922, # 35.922
            "2017": 41529, # 41.529
            "2018": 59740, # 59.74
            "2019": 137191
        }

        mc_samples = mc_datasets_dic[year]
        root_dic = []

        
        bcColors = [R.kGray, R.kGreen+1, R.kBlue-7, R.kOrange+1 ]
        sigColors = [2, 4, 1]

        with open("resources/run2xsecs.json") as json_file:
            xsecDic = json.load(json_file)

        for file in os.listdir("/Users/mo/hep/Analysis/AnalysisCode/histoFiles_{}/{}/".format(modifier,year)):
            if file.endswith(".root") and not file.startswith("all"):
                nickname = file.replace(".root", "")
                xsec = xsecDic[nickname]
                nickname = nickname.replace("__", "/")
                dic = {}
                dic["fullname"] = file

                for k, v in mc_samples.items():
                    if nickname in k:
                        dic["isSignal"] = v.isSignal
                        dic["year"] = str(v.year)
                f = R.TFile.Open(
                    "/Users/mo/hep/Analysis/AnalysisCode/histoFiles_{}/{}/".format(modifier,year)+file)
                dic["TFile"] = f
                dic["xsec"] = xsec
                root_dic.append((nickname, dic))

        dataTFile = R.TFile.Open(datafile)
        variables = [
                    "bdtScore_onGGH_noSyst",
                    "csPhi_onGGH_noSyst",
                    "csTheta_onGGH_noSyst",
                    "detammj_onGGH_noSyst",
                    "dijet_deta_onGGH_noSyst",
                    "dijet_dphi_onGGH_noSyst",
                    "dijet_mass_onGGH_noSyst",
                    "dimuon_mass_cat0_noSyst",
                    "dimuon_mass_cat1_noSyst",
                    "dimuon_mass_cat2_noSyst",
                    "dimuon_mass_cat3_noSyst",
                    "dimuon_mass_cat4_noSyst",
                    "dimuon_mass_cat5_noSyst",
                    "dimuon_mass_cat6_noSyst",
                    "dimuon_mass_cat7_noSyst",
                    "dimuon_mass_cat8_noSyst",
                    "dimuon_deta_onGGH_noSyst",
                    "dimuon_dphi_onGGH_noSyst",
                    "dimuon_eta_onGGH_noSyst",
                    "dimuon_mass_onGGH_noSyst",
                    "dimuon_phi_onGGH_noSyst",
                    "dimuon_pt_onGGH_noSyst",
                    "dimuon_rap_onGGH_noSyst",
                    "dphimmj_onGGH_noSyst",
                    "lead_muon_eta_onGGH_noSyst",
                    "lead_muon_phi_onGGH_noSyst",
                    "leadMuon_ptOverM_onGGH_noSyst",
                    "lead_muon_pt_onGGH_noSyst",
                    "leadjet_eta_onGGH_noSyst",
                    "leadjet_pt_onGGH_noSyst",
                    "met_phi_onGGH_noSyst",
                    "met_pt_onGGH_noSyst",
                    "num_bjets_onGGH_noSyst",
                    "num_jets_onGGH_noSyst",
                    "sub_muon_eta_onGGH_noSyst",
                    "sub_muon_phi_onGGH_noSyst",
                    "subMuon_ptOverM_onGGH_noSyst",
                    "sub_muon_pt_onGGH_noSyst",
                    "subjet_eta_onGGH_noSyst",
                    "subjet_pt_onGGH_noSyst",
                    "zeppen_onGGH_noSyst",
                    ]
            
        for variable in variables:
            plot_variable(variable, modifier)