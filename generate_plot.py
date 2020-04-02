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

variable_label = {
    "lead_muon_eta_onH" : " \eta (\mu_{1}) ",
    "sub_muon_eta_onH" : "\eta (\mu_{2})",
    "dimuon_mass_onH" : "M_{\mu \mu}",
    "dimuon_pt_onH" : "p_{T} (\mu,\mu)",
    "dimuon_eta_onH" : "\eta (\mu,\mu)",
    "num_jets_onH" : "# jets",
    "num_bjets_onH" : "# bjets)",
    "leadjet_pt_onH" : "p_{T} (j_{1})",
    "leadjet_eta_onH" : "\eta (j_1)",
    "subjet_pt_onH" : "p_{T} (j_{2})",
    "dijet_mass1_onH" : "M_{j j}",
    "dijet_deta1_onH" : "\Delta \eta_{j j}",
    "met_pt_onH" : "p_{MET}",
    "met_phi_onH" : "\phi_{MET}",
    "zeppen_onH" : "zepen",
    "csTheta_onH" : "cos(\theta_{cs})",
    "csPhi_onH" : "\phi_{CS}",
    "h_bdtScore01jet" : "BDT Output",
    "dimuon_mass_cat0": "M_{\mu \mu}", 
    "dimuon_mass_cat1": "M_{\mu \mu}",
    "dimuon_mass_cat2": "M_{\mu \mu}", 
    "dimuon_mass_cat3": "M_{\mu \mu}",
    "dimuon_mass_cat4": "M_{\mu \mu}", 
    "dimuon_mass_cat5": "M_{\mu \mu}",
    "dimuon_mass_cat6": "M_{\mu \mu}", 
    "dimuon_mass_cat7": "M_{\mu \mu}",
    "dimuon_mass_cat8": "M_{\mu \mu}", 
    "dimuon_mass_cat9": "M_{\mu \mu}"
    }

# systematic_TGraphErrors = {}
variable_systematics = {
    "bdtScore_onZ_": 1,
    "bdtScore_onH_": 1,
    "bdtScore_onGGH_": 1,
    "dimuon_mass_cat0_": 1,
    "dimuon_mass_cat1_": 1,
    "dimuon_mass_cat2_": 1,
    "dimuon_mass_cat3_": 1,
    "dimuon_mass_cat4_": 1,
    "dimuon_mass_cat5_": 1,
    "dimuon_mass_cat6_": 1,
    "dimuon_mass_cat7_": 1,
    "dimuon_mass_cat8_": 1,
    "dimuon_mass_cat9_": 1,
    "leadjet_pt_onZ_": 2,
    "subjet_pt_onZ_": 2,
    "dijet_mass1_onZ_": 2,
    "leadjet_pt_onH_": 2,
    "subjet_pt_onH_": 2,
    "dijet_mass1_onH_": 2,
    "leadjet_pt_onGGH_": 2,
    "subjet_pt_onGGH_": 2,
    "dijet_mass1_onGGH_" : 2,
    "dimuon_mass_onZ_": 3,
    "dimuon_mass_onH_": 3,
    "dimuon_mass_onGGH_": 3
}

systematics_map = {
    1: [
        "id",
        "iso",
        "jetTotal_",
        #  "pileup",
        "prefire",
        "trig"
        ],

    2: [
        "jetTotal_"
        ],
    
    3: [
        "id",
        "iso",
        # "jetTotal_",
        #  "pileup",
        "prefire",
        "trig"
        ]

}

def calculate_unc(variable, nominalRatio):
    systematics = systematics_map[variable_systematics[variable.split("noSyst")[0]]]

    # try:
    hdata = dataTFile.Get(variable).Clone()
    # except:
        # print(variable)
        # exit()

    combined_systematics_dic = {}

    # bkgStack = R.THStack("bkgStack", "")
    # for v in root_dic:
    #     histo = v[1]["TFile"].Get(variable).Clone()
    #     # print(v[0], histo.Integral())
    #     histo.Scale(lumi[v[1]["year"]] * v[1]["xsec"])
    #     bkgStack.Add(histo)

    # stack = bkgStack.GetStack().Last()


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

    for sys in systematics:
        bkgStack_up = R.THStack("bkgStackUp", "")
        bkgStack_down = R.THStack("bkgStackDown","")

        varUp = variable.split("noSyst")[0] + sys + "Up"
        varDown = variable.split("noSyst")[0] + sys + "Down"
    
        # print(variable, varUp, varDown)
        for v in root_dic:
            # if not v[1]["isSignal"]:
            # print(v[0])
            histoUp = v[1]["TFile"].Get(varUp).Clone()
            histoDown = v[1]["TFile"].Get(varDown).Clone()
            histoUp.Scale(lumi[v[1]["year"]] * v[1]["xsec"])
            histoDown.Scale(lumi[v[1]["year"]] * v[1]["xsec"])
            bkgStack_up.Add(histoUp)
            bkgStack_down.Add(histoDown)

        stackUp = hdata.Clone()
        stackUp.Divide(bkgStack_up.GetStack().Last())
        stackDown = hdata.Clone()
        stackDown.Divide(bkgStack_down.GetStack().Last())

        # print( stackUp.Integral() )
        # print( stack.Integral() )
        # print( stackDown.Integral() )

    
            # now the (bgkStack_up - bgkStack_down) / 2 should estimate the uncertainty?
        combined_systematics_up = []
        combined_systematics_down = []
        # print(len(combined_systematics_up))
        # print(nbins)
        # total_syst_error = 0
        for i in range(1,nbins+1):
            if (stackUp.GetBinContent(i)):
                combined_systematics_up.append(stackUp.GetBinError(i) / stackUp.GetBinContent(i))
                combined_systematics_down.append(stackDown.GetBinError(i) / stackDown.GetBinContent(i))
            else:
                combined_systematics_up.append(0)
                combined_systematics_down.append(0)

        combined_systematics_dic['stat'] = (np.asarray(combined_systematics_up), np.asarray(combined_systematics_down))
        # if (hnominal->GetBinContent(i) > 0) g_ratio_err->SetBinError(i, hnominal->GetBinError(i)/hnominal->GetBinContent(i));
        
        combined_systematics_up = []
        combined_systematics_down = []

        for i in range(1,nbins+1):
            if stackUp.GetBinContent(i):
                # print(sys, variable, stackDown.GetBinContent(i), stack.GetBinContent(i), stackUp.GetBinContent(i))
                # print(sys, variable, stackUp.GetBinContent(i)/stack.GetBinContent(i), abs(stackDown.GetBinContent(i) / stack.GetBinContent(i) ))
                # print( abs(stackUp.GetBinContent(i) - stackDown.GetBinContent(i)) / 2)
                error_up = abs(stackUp.GetBinContent(i) - nominalRatio.GetBinContent(i))
                error_down = abs(stackDown.GetBinContent(i) - nominalRatio.GetBinContent(i))

                # syst_error = max(error_up, error_down)
                # total_syst_error += syst_error
                # eyl.append( stackUp.GetBinError(i) / stackUp.GetBinContent(i) )
                # eyh.append( stackDown.GetBinError(i) / stackDown.GetBinContent(i) )
                # eyl.append( (stackDown.GetBinContent(i) / stack.GetBinContent(i)) )
                # eyh.append( (stackUp.GetBinContent(i) / stack.GetBinContent(i)) )    
                combined_systematics_up.append(error_up)
                combined_systematics_down.append(error_down)
                # eyl.append(0)
                # eyh.append(0)
            else:
                combined_systematics_up.append(0)
                combined_systematics_down.append(0)

        # print(len(combined_systematics))
        combined_systematics_dic[sys] = (np.asarray(combined_systematics_up), np.asarray(combined_systematics_down))
        # combined_systematics_dic[sys+"Down"] = combined_systematics_down


    total_up = np.zeros(nbins)
    total_down = np.zeros(nbins)
    # print(combined_systematics_dic)
    for _, combined_systematics in combined_systematics_dic.items():
        total_up += combined_systematics[0]**2
        total_down += combined_systematics[1]**2

        # print(combined_systematics[0])
        # print(combined_systematics[1])
        
        
        # for i in range(len(combined_systematics)):
            # total_error = eyl[i]
            # eyl[i] = ( sqrt(total_error * total_error + combined_systematics[i] * combined_systematics[i]) )
            # eyh[i] = ( sqrt(total_error * total_error + combined_systematics[i] * combined_systematics[i]) )

    # print( nbins, len(x), len(y) , len(exl), len(exh), len(eyl), len(eyh) )

    print(np.sqrt(total_up))
    print(np.sqrt(total_down))
    # print(eyl)
    # print(eyh)
    # exit()
    g = R.TGraphAsymmErrors(nbins, x, y, exl, exh, np.sqrt(total_up), np.sqrt(total_down))
    g.SetFillColor(2)
    g.SetFillStyle(3144)

    return g

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
    stStack = R.THStack("stStack", "")
    ttStack = R.THStack("ttStack", "")
    ttVStack = R.THStack("ttVStack", "")
    dyStack = R.THStack("dyStack", "")
    ewkStack = R.THStack("ewkStack", "")
    diBosonStack = R.THStack("diBosonStack", "")
    triBosonStack = R.THStack("triBosonStack", "")
    ggVVStack = R.THStack("ggVVStack","")
    tzqStack = R.THStack("tzqStack","")

    sigStack_120 = R.THStack("sigStack_120", "")
    tthStack_120 = R.THStack("tthStack_120","ttH")
    vbfStack_120 = R.THStack("vbfStack_120","VBF")
    wmHStack_120 = R.THStack("wmHStack_120","wmH")
    wpHStack_120 = R.THStack("wpHStack_120","wpH")
    gghStack_120 = R.THStack("gghStack_120","ggH")
    zHStack_120 = R.THStack("zHStack_120","ZH")

    sigStack_125 = R.THStack("sigStack_125", "")
    tthStack_125 = R.THStack("tthStack_125","ttH")
    vbfStack_125 = R.THStack("vbfStack_125","VBF")
    wmHStack_125 = R.THStack("wmHStack_125","wmH")
    wpHStack_125 = R.THStack("wpHStack_125","wpH")
    gghStack_125 = R.THStack("gghStack_125","ggH")
    zHStack_125 = R.THStack("zHStack_125","ZH")

    sigStack_130 = R.THStack("sigStack_130", "")
    tthStack_130 = R.THStack("tthStack_130","ttH")
    vbfStack_130 = R.THStack("vbfStack_130","VBF")
    wmHStack_130 = R.THStack("wmHStack_130","wmH")
    wpHStack_130 = R.THStack("wpHStack_130","wpH")
    gghStack_130 = R.THStack("gghStack_130","ggH")
    zHStack_130 = R.THStack("zHStack_130","ZH")


    allStack = R.THStack("allStack", "")
    print(variable)
    hdata = dataTFile.Get(variable).Clone()
    # Nbins = hdata.GetNbinsX()
    # lowbin = hdata.GetBinLowEdge(0)+hdata.GetBinWidth(0)
    # highbin = hdata.GetBinLowEdge(Nbins)+hdata.GetBinWidth(Nbins)
    # tthStack = R.TH1F("tthStack","",Nbins, lowbin, highbin)
    # vbfStack  = R.TH1F("vbfStack","",Nbins, lowbin, highbin)
    # wmHStack  = R.TH1F("wmHStack","",Nbins, lowbin, highbin)
    # wpHStack  = R.TH1F("wpHStack","",Nbins, lowbin, highbin)
    # gghStack  = R.TH1F("gghStack","",Nbins, lowbin, highbin)
    # zHStack  = R.TH1F("zHStack","",Nbins, lowbin, highbin)

    leg = R.TLegend(0.45, 0.70, .88, .88)
    # leg.SetHeader("Samples")
    leg.SetNColumns(3)
    # leg.SetTextSize(12)
    leg.SetBorderSize(0)
    fb = "fb^{-1}"
    lumis = float(lumi[year])/1000
    if "2019" in year:
        # print(year)
        leg.AddEntry(hdata, "Run II Data")
    else:
        leg.AddEntry(hdata, "{} Data".format(year))

    # sigHistos = []
    for v in root_dic:
        # print(v)
        if v[1]["isSignal"]:
            histo = v[1]["TFile"].Get(variable).Clone()
            histo.Scale(lumi[v[1]["year"]] * v[1]["xsec"])

            if "M120" in v[0]:
                if "ttH" in v[0]:
                    tthStack_120.Add(histo)
                elif "WminusH" in v[0]:
                    wmHStack_120.Add(histo)
                elif "WplusH" in v[0]:
                    wpHStack_120.Add(histo)
                elif "ZH" in v[0]:
                    zHStack_120.Add(histo)
                elif "GluGluH" in v[0]:
                    gghStack_120.Add(histo)
                elif "VBFH" in v[0]:
                    vbfStack_120.Add(histo)
                sigStack_120.Add(histo)

            elif "M125" in v[0]:
                if "ttH" in v[0]:
                    tthStack_125.Add(histo)
                elif "WminusH" in v[0]:
                    wmHStack_125.Add(histo)
                elif "WplusH" in v[0]:
                    wpHStack_125.Add(histo)
                elif "ZH" in v[0]:
                    zHStack_125.Add(histo)
                elif "GluGluH" in v[0]:
                    gghStack_125.Add(histo)
                elif "VBFH" in v[0]:
                    vbfStack_125.Add(histo)
                sigStack_125.Add(histo)
                allStack.Add(histo)

            elif "M130" in v[0]:
                if "ttH" in v[0]:
                    tthStack_130.Add(histo)
                elif "WminusH" in v[0]:
                    wmHStack_130.Add(histo)
                elif "WplusH" in v[0]:
                    wpHStack_130.Add(histo)
                elif "ZH" in v[0]:
                    zHStack_130.Add(histo)
                elif "GluGluH" in v[0]:
                    gghStack_130.Add(histo)
                elif "VBFH" in v[0]:
                    vbfStack_130.Add(histo)
                sigStack_130.Add(histo)

            # histo.SetLineColor(sigColors[(j % len(sigColors))])
            # leg.AddEntry(histo, v[0].split("_")[0], "l")
            # sigHistos.append(histo)

        if not v[1]["isSignal"]:
            histo = v[1]["TFile"].Get(variable).Clone()
            histo.Scale(lumi[v[1]["year"]] * v[1]["xsec"])
            if "ST" in v[0]:
                # print(v[0])
                # print("ST")
                stStack.Add(histo)
            elif "TTW" in v[0] or "TTZ" in v[0]:
                # print(v[0])
                # print("TTV")
                ttVStack.Add(histo)
            elif "TT" in v[0]:
                # print(v[0])
                # print("TT")
                ttStack.Add(histo)
            elif "WWW" in v[0] or "WWZ" in v[0] or "ZZZ" in v[0]:
                # print(v[0])
                # print("VVV")
                triBosonStack.Add(histo)
            elif "ZZ" in v[0] or "WZ" in v[0] or "WW" in v[0]:
                # print(v[0])
                # print("VV")
                diBosonStack.Add(histo)
            elif "DY" in v[0]:
#                # if "onZ" in variable and "M-50" in v[0]:
                dyStack.Add(histo)
                # elif "onH" in variable and "M-150" in v[0]:
                    # dyStack.Add(histo)

            elif "EWK" in v[0]:
                ewkStack.Add(histo)
            elif "tZq" in v[0]:
                tzqStack.Add(histo)
            elif "GluGluToContin" in v[0]:
                ggVVStack.Add(histo)

            else:
                print(v[0])
                print("This sample is not in any stack??")
            # histo.SetLineWidth(0)
            # histo.SetLineColor(R.kBlack)
            # histo.SetFillColor(bcColors[(i % len(bcColors))])
            allStack.Add(histo)

    # make this None, then get the .Last() using a
    # try: get() except: pass (keep moving)
    # then later on do if is not None
    hist_stStack = None
    hist_ttStack = None
    hist_ttVStack = None
    hist_dyStack = None
    hist_ewkStack = None
    hist_diBosonStack = None
    hist_triBosonStack = None
    hist_tzqStack = None
    hist_ggVVStack = None

    # hist_tthStack_120 = tthStack_120.GetStack().Last()
    # hist_tthStack_120.SetName("ttH_120_{}".format(variable))
    # hist_vbfStack_120 = vbfStack_120.GetStack().Last()
    # hist_vbfStack_120.SetName("vbf_120_{}".format(variable))
    # hist_wmHStack_120 = wmHStack_120.GetStack().Last()
    # hist_wmHStack_120.SetName("wmH_120_{}".format(variable))
    # hist_wpHStack_120 = wpHStack_120.GetStack().Last()
    # hist_wpHStack_120.SetName("wpH_120_{}".format(variable))
    # hist_gghStack_120 = gghStack_120.GetStack().Last()
    # hist_gghStack_120.SetName("ggH_120_{}".format(variable))
    # hist_zHStack_120 = zHStack_120.GetStack().Last()
    # hist_zHStack_120.SetName("zH_120_{}".format(variable))

    hist_tthStack_125 = tthStack_125.GetStack().Last()
    hist_tthStack_125.SetName("ttH_125_{}".format(variable))
    hist_vbfStack_125 = vbfStack_125.GetStack().Last()
    hist_vbfStack_125.SetName("vbf_125_{}".format(variable))
    hist_wmHStack_125 = wmHStack_125.GetStack().Last()
    hist_wmHStack_125.SetName("wmH_125_{}".format(variable))
    hist_wpHStack_125 = wpHStack_125.GetStack().Last()
    hist_wpHStack_125.SetName("wpH_125_{}".format(variable))
    hist_gghStack_125 = gghStack_125.GetStack().Last()
    hist_gghStack_125.SetName("ggH_125_{}".format(variable))
    hist_zHStack_125 = zHStack_125.GetStack().Last()
    hist_zHStack_125.SetName("zH_125_{}".format(variable))

    # hist_tthStack_130 = tthStack_130.GetStack().Last()
    # hist_tthStack_130.SetName("ttH_130_{}".format(variable))
    # hist_vbfStack_130 = vbfStack_130.GetStack().Last()
    # hist_vbfStack_130.SetName("vbf_130_{}".format(variable))
    # hist_wmHStack_130 = wmHStack_130.GetStack().Last()
    # hist_wmHStack_130.SetName("wmH_130_{}".format(variable))
    # hist_wpHStack_130 = wpHStack_130.GetStack().Last()
    # hist_wpHStack_130.SetName("wpH_130_{}".format(variable))
    # hist_gghStack_130 = gghStack_130.GetStack().Last()
    # hist_gghStack_130.SetName("ggH_130_{}".format(variable))
    # hist_zHStack_130 = zHStack_130.GetStack().Last()
    # hist_zHStack_130.SetName("zH_130_{}".format(variable))

    try:
        hist_stStack = stStack.GetStack().Last()
        leg.AddEntry(hist_stStack, "st", "F")
    except:
        pass
        
    try:
        hist_ttStack = ttStack.GetStack().Last()
        leg.AddEntry(hist_ttStack, "tt", "F")
    except:
        pass

    try:
        hist_ttVStack = ttVStack.GetStack().Last()
        leg.AddEntry(hist_ttVStack, "ttV", "F")

    except:
        pass
  
    try:
        hist_dyStack = dyStack.GetStack().Last()
        leg.AddEntry(hist_dyStack, "DY", "F")

    except:
        pass
    
    try:
        hist_ewkStack = ewkStack.GetStack().Last()
        leg.AddEntry(hist_ewkStack, "EWK", "F")

    except:
        pass

    try:
        hist_diBosonStack = diBosonStack.GetStack().Last()
        leg.AddEntry(hist_diBosonStack, "VV", "F")

    except:
        pass
  
    try:
        hist_triBosonStack = triBosonStack.GetStack().Last()
        leg.AddEntry(hist_triBosonStack, "VVV", "F")
    except:
        pass
    
    try:
        hist_tzqStack = tzqStack.GetStack().Last()
        leg.AddEntry(hist_tzqStack, "tZq", "F")

    except:
        pass
   
    try:
        hist_ggVVStack = ggVVStack.GetStack().Last()
        leg.AddEntry(hist_ggVVStack, "ggVV", "F")

    except:
        pass
    stackList = [hist_tzqStack, hist_ggVVStack, hist_diBosonStack, hist_triBosonStack, hist_ewkStack,
                 hist_ttVStack, hist_stStack, hist_ttStack, hist_dyStack]
    # stackList = [hist_ttStack]
    for i in range(0, len(stackList)):
        if stackList[i] is not None:
            stackList[i].SetFillColor(bcColors[(i % len(bcColors))])
            bkgStack.Add(stackList[i])

    # blind data in 120-130 GeV bins on mass plots
    if "dimuon_mass" in variable:
        for i in range(hdata.GetNbinsX()):
            if hdata.GetBinCenter(i + 1) > 120 and hdata.GetBinCenter(i + 1) < 130:
                hdata.SetBinContent(i + 1, 0)
                hdata.SetBinError(i + 1, 0)

    hdata.SetMarkerStyle(20)
    hdata.SetMarkerSize(.5)

    if logy:
        pad1.SetLogy()
        bkgStack.SetMinimum(.1)

    bkgStack.SetMaximum(hdata.GetMaximum() * 100)
    # bkgStack.SetMaximum(1000000)
    hdata.SetStats(R.kFALSE)
    hdata.GetXaxis().SetLabelSize(0)
    bkgStack.Draw("hist")

    j = 0
    sigStackList = [hist_tthStack_125, hist_vbfStack_125, hist_wmHStack_125, hist_wpHStack_125, hist_gghStack_125, hist_zHStack_125]
    for h in sigStackList:
            h.SetLineColor(sigColors[(j % len(sigColors))])
            h.Draw("hist same")
            j+=1

    # sigStackList += [hist_tthStack_120, hist_vbfStack_120, hist_wmHStack_120, hist_wpHStack_120, hist_gghStack_120, hist_zHStack_120]
    # sigStackList += [hist_tthStack_130, hist_vbfStack_130, hist_wmHStack_130, hist_wpHStack_130, hist_gghStack_130, hist_zHStack_130]    
    leg.AddEntry(hist_tthStack_125,"ttH","l")
    leg.AddEntry(hist_vbfStack_125,"VBF","l")
    leg.AddEntry(hist_wmHStack_125,"W^{-}H","l")
    leg.AddEntry(hist_wpHStack_125,"W^{+}H","l")
    leg.AddEntry(hist_gghStack_125,"ggH","l")
    leg.AddEntry(hist_zHStack_125,"ZH","l")

    hdata.Draw("SAME P")
    # bkgStack.GetXaxis().SetTitle(variable)
    bkgStack.GetYaxis().SetTitle("# Events")
    leg.Draw()
    # pad1.BuildLegend()
    R.gPad.Modified()

    pad2.cd()
    hratio = hdata.Clone()
    hratio.SetTitle("")
    hratio.GetYaxis().SetTitle("Data / MC")
    hratio.GetXaxis().SetTitle(variable)
    hratio.GetYaxis().SetNdivisions(4, R.kFALSE)
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
    hratio.Divide(allStack.GetStack().Last())
    hratio.SetStats(R.kFALSE)
    hratio.SetMaximum(1.4)
    hratio.SetMinimum(0.6)
    hratio.SetMarkerStyle(20)
    hratio.SetMarkerSize(0.5)
    # hratio.Draw("E P")

    if variable.split("noSyst")[0] in variable_systematics:
        print("systematics apply..?")
        g = calculate_unc(variable, hratio.Clone())
        g.SetTitle("")
        g.GetYaxis().SetTitle("Data / MC")
        g.GetXaxis().SetTitle(variable)
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
        g.SetMaximum(1.4)
        g.SetMinimum(0.6)
        g.GetXaxis().SetLimits(hratio.GetBinLowEdge(1), hratio.GetBinLowEdge(hratio.GetNbinsX())+hratio.GetBinWidth(1))
        g.Draw("a2")
        # g.Draw("p")


    hratio.Draw("E P same")

    R.gPad.Modified()
    canvas.Draw()
    canvas.SaveAs("plots_{}/{}/".format(modifier, year) +
                  variable + "_{}.pdf".format(year))
    
    if "2019" in year and "dimuon_mass" in variable:
        outfile.cd()
        hdata.Write()
        for h in sigStackList:
            h.Write()

            
    del canvas, pad1, pad2, hratio, bkgStack, hdata


if __name__ == "__main__":
    R.gROOT.SetBatch(R.kTRUE)
    years = ["2016", "2017", "2018"]
    # years = ["2019"] 
    years = ["2018"]
    # 5405018017057433 04/21 356
    if len(sys.argv) > 1:
        # years = [sys.argv[1]]
        modifier = sys.argv[1]
    else:
        print("need modifier!")
        exit(1)

    for year in years:
        if not os.path.isdir("plots_{}/{}".format(modifier,year)):
            os.makedirs("plots_{}/{}".format(modifier,year))
        datafile = "/Users/mo/hep/Analysis/AnalysisCode/histoFiles_{}/{}/allData{}.root".format(
            modifier, year, year)
        # print(datafile)

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
        # lumi = lumi[year]

        bcColors = [40, 30, 41, 42, 43, 35, 46, 47, 38, 28, 29]
        sigColors = [2, 3, 4, 6, 8, 9]

        with open("resources/run2xsecs.json") as json_file:
            xsecDic = json.load(json_file)

        for file in os.listdir("/Users/mo/hep/Analysis/AnalysisCode/histoFiles_{}/{}/".format(modifier,year)):
            # if "105To160" in file:
            #     continue
            if file.endswith(".root") and not file.startswith("all"):
                nickname = file.replace(".root", "")
                xsec = xsecDic[nickname]
                nickname = nickname.replace("__", "/")
                dic = {}
                dic["fullname"] = file

                for k, v in mc_samples.items():
                    if nickname in k:
                        # print(nickname)
                        dic["isSignal"] = v.isSignal
                        dic["year"] = str(v.year)
                f = R.TFile.Open(
                    "/Users/mo/hep/Analysis/AnalysisCode/histoFiles_{}/{}/".format(modifier,year)+file)
                dic["TFile"] = f
                dic["xsec"] = xsec
                root_dic.append((nickname, dic))
        # exit()
        # root_dic_sorted = sorted(root_dic, key=lambda el: el[1]['wEvents'])
        # root_dic_sorted = root_dic
        dataTFile = R.TFile.Open(datafile)
        variables = [
                    "bdtScore_onGGH_noSyst",
                    "bdtScore_onH_noSyst",
#                    "bdtScore_onZ_noSyst",
                    "csPhi_onGGH_noSyst",
                    "csPhi_onH_noSyst",
#                    "csPhi_onZ_noSyst",
                    "csTheta_onGGH_noSyst",
                    "csTheta_onH_noSyst",
#                    "csTheta_onZ_noSyst",
                    "dijet_deta1_onGGH_noSyst",
                    "dijet_deta1_onH_noSyst",
#                    "dijet_deta1_onZ_noSyst",
                    "dijet_mass1_onGGH_noSyst",
                    "dijet_mass1_onH_noSyst",
#                    "dijet_mass1_onZ_noSyst",
                    "dimuon_mass_cat0_noSyst",
                    "dimuon_mass_cat1_noSyst",
                    "dimuon_mass_cat2_noSyst",
                    "dimuon_mass_cat3_noSyst",
                    "dimuon_mass_cat4_noSyst",
                    "dimuon_mass_cat5_noSyst",
                    "dimuon_mass_cat6_noSyst",
                    "dimuon_mass_cat7_noSyst",
                    "dimuon_mass_cat8_noSyst",
                    "dimuon_mass_cat9_noSyst",
                    "dimuon_deta_onGGH_noSyst",
                    "dimuon_deta_onH_noSyst",
#                    "dimuon_deta_onZ_noSyst",
                    "dimuon_dphi_onGGH_noSyst",
                    "dimuon_dphi_onH_noSyst",
#                    "dimuon_dphi_onZ_noSyst",
                    "dimuon_eta_onGGH_noSyst",
                    "dimuon_eta_onH_noSyst",
#                    "dimuon_eta_onZ_noSyst",
                    "dimuon_mass_onGGH_noSyst",
                    "dimuon_mass_onH_noSyst",
                   "dimuon_mass_onZ_noSyst",
                    "dimuon_phi_onGGH_noSyst",
                    "dimuon_phi_onH_noSyst",
#                    "dimuon_phi_onZ_noSyst",
                    "dimuon_pt_onGGH_noSyst",
                    "dimuon_pt_onH_noSyst",
#                    "dimuon_pt_onZ_noSyst",
                    "lead_muon_eta_onGGH_noSyst",
                    "lead_muon_eta_onH_noSyst",
#                    "lead_muon_eta_onZ_noSyst",
                    "lead_muon_phi_onGGH_noSyst",
                    "lead_muon_phi_onH_noSyst",
                   "lead_muon_phi_onZ_noSyst",
                    "lead_muon_pt_onGGH_noSyst",
                    "lead_muon_pt_onH_noSyst",
#                    "lead_muon_pt_onZ_noSyst",
                    "leadjet_eta_onGGH_noSyst",
                    "leadjet_eta_onH_noSyst",
#                    "leadjet_eta_onZ_noSyst",
                    "leadjet_pt_onGGH_noSyst",
                    "leadjet_pt_onH_noSyst",
#                    "leadjet_pt_onZ_noSyst",
                    "met_phi_onGGH_noSyst",
                    "met_phi_onH_noSyst",
#                    "met_phi_onZ_noSyst",
                    "met_pt_onGGH_noSyst",
                    "met_pt_onH_noSyst",
#                    "met_pt_onZ_noSyst",
                    "num_bjets_onGGH_noSyst",
                    "num_bjets_onH_noSyst",
#                    "num_bjets_onZ_noSyst",
                    "num_jets_onGGH_noSyst",
                    "num_jets_onH_noSyst",
#                    "num_jets_onZ_noSyst",
                    "sub_muon_eta_onGGH_noSyst",
                    "sub_muon_eta_onH_noSyst",
#                    "sub_muon_eta_onZ_noSyst",
                    "sub_muon_phi_onGGH_noSyst",
                    "sub_muon_phi_onH_noSyst",
#                    "sub_muon_phi_onZ_noSyst",
                    "sub_muon_pt_onGGH_noSyst",
                    "sub_muon_pt_onH_noSyst",
#                    "sub_muon_pt_onZ_noSyst",
                    "subjet_eta_onGGH_noSyst",
                    "subjet_eta_onH_noSyst",
#                    "subjet_eta_onZ_noSyst",
                    "subjet_pt_onGGH_noSyst",
                    "subjet_pt_onH_noSyst",
#                    "subjet_pt_onZ_noSyst",
                    "zeppen_onGGH_noSyst",
                    "zeppen_onH_noSyst"
#                   "zeppen_onZ_noSyst"]
]
            

        # variables += ["dimuon_mass_cat0", "dimuon_mass_cat1",
                    #   "dimuon_mass_cat2", "dimuon_mass_cat3",
                    #   "dimuon_mass_cat4", "dimuon_mass_cat5",
                    #   "dimuon_mass_cat6", "dimuon_mass_cat7",
                    #   "dimuon_mass_cat8", "dimuon_mass_cat9"]
        # for jetMass in range(0,510,10):
        #     histo_name = 'dimuon_mass_jet_{}'.format(jetMass)
        #     variables.append(histo_name)

        if "2019" in year:
            outfile = R.TFile.Open("FitHistosAll.root", "RECREATE")
        

        for variable in variables:
            # print(variable)
            # calculate_unc(variable)
            plot_variable(variable, modifier)
        if "2019" in year:
            outfile.Close()

    # plot_muon_corr()