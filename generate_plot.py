import os, sys
import json

os.environ["ANALYSISHOME"] = "/Users/mo/hep/Analysis/HMuMu"
sys.path.append(os.path.join(
    os.environ["ANALYSISHOME"], "Configuration", "higgs"))

import Samples as S
import ROOT as R

def plot_variable(variable, modifier, logy=True):
    canvas = R.TCanvas("c1", "c1")
    canvas.cd()
    pad1 = R.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.Draw()
    canvas.cd()
    pad2 = R.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.Draw()
    pad2.cd()
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.2)
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



    sigStack = R.THStack("sigStack", "")
    tthStack = R.THStack("tthStack","")
    vbfStack = R.THStack("vbfStack","")
    wmHStack = R.THStack("wmHStack","")
    wpHStack = R.THStack("wpHStack","")
    gghStack = R.THStack("gghStack","")
    zHStack  = R.THStack("zHStack","")

    allStack = R.THStack("allStack", "")
    # print(variable)
    hdata = dataTFile.Get(variable)
    # Nbins = hdata.GetNbinsX()
    # lowbin = hdata.GetBinLowEdge(0)+hdata.GetBinWidth(0)
    # highbin = hdata.GetBinLowEdge(Nbins)+hdata.GetBinWidth(Nbins)
    # tthStack = R.TH1F("tthStack","",Nbins, lowbin, highbin)
    # vbfStack  = R.TH1F("vbfStack","",Nbins, lowbin, highbin)
    # wmHStack  = R.TH1F("wmHStack","",Nbins, lowbin, highbin)
    # wpHStack  = R.TH1F("wpHStack","",Nbins, lowbin, highbin)
    # gghStack  = R.TH1F("gghStack","",Nbins, lowbin, highbin)
    # zHStack  = R.TH1F("zHStack","",Nbins, lowbin, highbin)

    leg = R.TLegend(0.6, 0.7, .89, .89)
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
        if v[1]["isSignal"] and "M125" in v[0]:
            histo = v[1]["TFile"].Get(variable)
            histo.Scale(lumi[v[1]["year"]] * v[1]["xsec"])
            if "ttH" in v[0]:
                tthStack.Add(histo)
            elif "WminusH" in v[0]:
                wmHStack.Add(histo)
            elif "WplusH" in v[0]:
                wpHStack.Add(histo)
            elif "ZH" in v[0]:
                zHStack.Add(histo)
            elif "GluGluH" in v[0]:
                gghStack.Add(histo)
            elif "VBFH" in v[0]:
                vbfStack.Add(histo)
            # histo.SetLineColor(sigColors[(j % len(sigColors))])
            # leg.AddEntry(histo, v[0].split("_")[0], "l")
            # sigHistos.append(histo)
            sigStack.Add(histo)
            allStack.Add(histo)

    for v in root_dic:
        if not v[1]["isSignal"]:
            histo = v[1]["TFile"].Get(variable)
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
                dyStack.Add(histo)
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

    hist_tthStack = tthStack.GetStack().Last()
    hist_vbfStack = vbfStack.GetStack().Last()
    hist_wmHStack = wmHStack.GetStack().Last()
    hist_wpHStack = wpHStack.GetStack().Last()
    hist_gghStack = gghStack.GetStack().Last()
    hist_zHStack  = zHStack.GetStack().Last()

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
    if "mass" in variable:
        for i in range(hdata.GetNbinsX()):
            if hdata.GetBinCenter(i + 1) > 120 and hdata.GetBinCenter(i + 1) < 130:
                hdata.SetBinContent(i + 1, 0)
                hdata.SetBinError(i + 1, 0)

    hdata.SetMarkerStyle(20)
    hdata.SetMarkerSize(.5)

    if logy:
        pad1.SetLogy()
        bkgStack.SetMinimum(.001)

    bkgStack.SetMaximum(hdata.GetMaximum() * 100)
    # bkgStack.SetMaximum(1000000)
    hdata.SetStats(R.kFALSE)
    hdata.GetXaxis().SetLabelSize(0)
    bkgStack.Draw("hist")

    j = 0
    sigStackList = [hist_tthStack, hist_vbfStack, hist_wmHStack, hist_wpHStack, hist_gghStack, hist_zHStack ]
    for h in sigStackList:
            h.SetLineColor(sigColors[(j % len(sigColors))])
            h.Draw("hist same")
            j+=1

    leg.AddEntry(hist_tthStack,"ttH","l")
    leg.AddEntry(hist_vbfStack,"VBF","l")
    leg.AddEntry(hist_wmHStack,"W^{-}H","l")
    leg.AddEntry(hist_wpHStack,"W^{+}H","l")
    leg.AddEntry(hist_gghStack,"ggH","l")
    leg.AddEntry(hist_zHStack,"ZH","l")

    hdata.Draw("SAME P")
    bkgStack.GetXaxis().SetTitle(variable)
    bkgStack.GetYaxis().SetTitle("# Events")
    leg.Draw()
    # pad1.BuildLegend()
    R.gPad.Modified()

    pad2.cd()
    hratio = hdata.Clone()
    hratio.SetTitle("")
    hratio.GetYaxis().SetTitle("Data / MC")
    hratio.GetXaxis().SetTitle(variable)
    hratio.GetYaxis().SetNdivisions(6, R.kFALSE)
    hratio.GetYaxis().SetTitleSize(10)
    hratio.GetYaxis().SetTitleFont(43)
    hratio.GetYaxis().SetTitleOffset(1.55)
    hratio.GetYaxis().SetLabelFont(43)
    hratio.GetYaxis().SetLabelSize(15)
    hratio.GetXaxis().SetTitleSize(10)
    hratio.GetXaxis().SetTitleFont(43)
    hratio.GetXaxis().SetTitleOffset(4)
    hratio.GetXaxis().SetLabelFont(43)
    hratio.GetXaxis().SetLabelSize(15)
    hratio.Divide(bkgStack.GetStack().Last())
    hratio.SetStats(R.kFALSE)
    hratio.Draw("E P")
    hratio.SetMaximum(1.6)
    hratio.SetMinimum(0.4)
    hratio.SetMarkerStyle(20)
    hratio.SetMarkerSize(0.5)
    R.gPad.Modified()
    canvas.Draw()
    canvas.SaveAs("plots_{}/{}/".format(modifier, year) +
                  variable + "_{}.pdf".format(year))
    # outfile.cd()
#    canvas.Write()
    del canvas, pad1, pad2, hratio, bkgStack, sigStack, hdata


if __name__ == "__main__":
    years = ["2016", "2017", "2018"]
    # years = ["2019"] 
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
        mc2016.update(S.mc_background_2016_extra)
        mc2016.update(S.mc_background_2016_extra_2)


        mc2017 = S.mc_background_2017
        mc2017.update(S.mc_signal_2017)
        mc2017.update(S.mc_background_2017_extra)
        mc2017.update(S.mc_background_2017_extra_2)


        mc2018 = S.mc_background_2018
        mc2018.update(S.mc_signal_2018)
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
            "2016": 35867,
            "2017": 41860,
            "2018": 58860,
            "2019": 136587
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
            # 'lead_muon_pt_onZ', 'lead_muon_eta_onZ', 'sub_muon_pt_onZ', 'sub_muon_eta_onZ',
            # 'dimuon_mass_onZ', 'dimuon_pt_onZ', 'dimuon_eta_onZ', 'dimuon_phi_onZ',
            # 'dimuon_deta_onZ', 'dimuon_dphi_onZ', 'num_jets_onZ', 'num_bjets_onZ',
            # 'leadjet_pt_onZ', 'leadjet_eta_onZ', 'subjet_pt_onZ', 'subjet_eta_onZ',
            # 'dijet_mass1_onZ', 'dijet_deta1_onZ', 'met_pt_onZ', 'met_phi_onZ',
            # 'zeppen_onZ', 'csTheta_onZ', 'csPhi_onZ',
            # 'lead_muon_pt_onH', 'lead_muon_eta_onH', 'sub_muon_pt_onH', 'sub_muon_eta_onH',
            # 'dimuon_mass_onH', 'dimuon_pt_onH', 'dimuon_eta_onH', 'dimuon_phi_onH',
            # 'dimuon_deta_onH', 'dimuon_dphi_onH', 'num_jets_onH', 'num_bjets_onH',
            # 'leadjet_pt_onH', 'leadjet_eta_onH', 'subjet_pt_onH', 'subjet_eta_onH',
            # 'dijet_mass1_onH', 'dijet_deta1_onH', 'met_pt_onH', 'met_phi_onH',
            'lead_muon_eta_onH', 'sub_muon_eta_onH',
            'dimuon_mass_onH', 'dimuon_pt_onH', 'dimuon_eta_onH', 
            'num_jets_onH', 'num_bjets_onH',
            'leadjet_pt_onH', 'leadjet_eta_onH', 'subjet_pt_onH',
            'dijet_mass1_onH', 'dijet_deta1_onH', 'met_pt_onH', 'met_phi_onH',
            # 'zeppen_onH', 'csTheta_onH', 'csPhi_onH',
            'h_bdtScore01jet']

        variables += ["dimuon_mass_cat0", "dimuon_mass_cat1",
                      "dimuon_mass_cat2", "dimuon_mass_cat3",
                      "dimuon_mass_cat4", "dimuon_mass_cat5",
                      "dimuon_mass_cat6", "dimuon_mass_cat7",
                      "dimuon_mass_cat8", "dimuon_mass_cat9"]
        # for jetMass in range(0,510,10):
        #     histo_name = 'dimuon_mass_jet_{}'.format(jetMass)
        #     variables.append(histo_name)

        # outfile = R.TFile.Open("plots/{}/plots_{}.root".format(year,year), "RECREATE")
        for variable in variables:
            # print(variable)
            plot_variable(variable, modifier)
        # outfile.Close()

    # plot_muon_corr()
