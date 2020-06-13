#include <TH1.h>
#include "eventEntry.h"

class histosMap
{
private:
    std::map<std::string, TH1F *> histosMap_;

public:
    histosMap(std::string systName)
    {
        addMap(systName);
    }

    void addMap(std::string systName)
    {
        TH1::SetDefaultSumw2();

        // ggH-enriched Category Plots
        histosMap_["leadMuon_pt_onGGH"] = new TH1F(("lead_muon_pt_onGGH_" + systName).c_str(), "Leading Muon p_{T};p_{T}  (GeV);Events", 60, 20, 200);
        histosMap_["leadMuon_eta_onGGH"] = new TH1F(("lead_muon_eta_onGGH_" + systName).c_str(), "Leading Muon \\eta;\\eta;Events", 50, -2.5, 2.5);
        histosMap_["leadMuon_phi_onGGH"] = new TH1F(("lead_muon_phi_onGGH_" + systName).c_str(), "Leading Muon \\phi;\\phi;Events", 64, -3.2, 3.2);
        histosMap_["leadMuon_ptOverM_onGGH"] = new TH1F(("leadMuon_ptOverM_onGGH_" + systName).c_str(), "Leading Muon p_{T} Over Mass ;p_{T}/M  (GeV);Events", 30, 0, 1.5);

        histosMap_["subMuon_pt_onGGH"] = new TH1F(("sub_muon_pt_onGGH_" + systName).c_str(), "Subleading Muon p_{T};p_{T}  (GeV);Events", 60, 20, 200);
        histosMap_["subMuon_eta_onGGH"] = new TH1F(("sub_muon_eta_onGGH_" + systName).c_str(), "Subleading Muon eta;\\eta;Events", 50, -2.5, 2.5);
        histosMap_["subMuon_phi_onGGH"] = new TH1F(("sub_muon_phi_onGGH_" + systName).c_str(), "Subleading Muon \\phi;\\phi;Events", 64, -3.2, 3.2);
        histosMap_["subMuon_ptOverM_onGGH"] = new TH1F(("subMuon_ptOverM_onGGH_" + systName).c_str(), "Subleading Muon p_{T} Over Mass;p_{T}/M (GeV);Events", 30, 0, 1.5);

        histosMap_["dimuon_mass_onGGH"] = new TH1F(("dimuon_mass_onGGH_" + systName).c_str(), "Dimuon Mass;M_{\\mu \\mu}  (Gev);Events", 80, 110, 150);
        histosMap_["dimuon_pt_onGGH"] = new TH1F(("dimuon_pt_onGGH_" + systName).c_str(), "Dimuon p_{T};p_{T}  (GeV);Events", 50, 0, 200);
        histosMap_["dimuon_eta_onGGH"] = new TH1F(("dimuon_eta_onGGH_" + systName).c_str(), "Dimuon \\eta;\\eta;Events", 50, -5.0, 5.0);
        histosMap_["dimuon_phi_onGGH"] = new TH1F(("dimuon_phi_onGGH_" + systName).c_str(), "Dimuon \\phi;\\phi;Events", 36, -3.6, 3.6);
        histosMap_["dimuon_rap_onGGH"] = new TH1F(("dimuon_rap_onGGH_"+ systName).c_str(), "Dimuon rap; rap; Events", 40, -2.5, 2.5);

        histosMap_["dimuon_deta_onGGH"] = new TH1F(("dimuon_deta_onGGH_" + systName).c_str(), "Dimuon deta;deta;Events", 50, 0, 5.0);
        histosMap_["dimuon_dphi_onGGH"] = new TH1F(("dimuon_dphi_onGGH_" + systName).c_str(), "Dimuon dphi;dphi;Events", 18, 0, 3.6);

        histosMap_["num_jets_onGGH"] = new TH1F(("num_jets_onGGH_" + systName).c_str(), "Number of Jets;nJets;Events", 8, 0, 8);
        histosMap_["num_bjets_onGGH"] = new TH1F(("num_bjets_onGGH_" + systName).c_str(), "Number of B Jets;nBJets;Events", 6, 0, 6);

        histosMap_["leadjet_pt_onGGH"] = new TH1F(("leadjet_pt_onGGH_" + systName).c_str(), "Leading Jet p_{T};p_{T}  (GeV);Events", 45, 25, 250);
        histosMap_["leadjet_eta_onGGH"] = new TH1F(("leadjet_eta_onGGH_" + systName).c_str(), "Leading Jet \\eta;\\eta;Events", 24, -4.8, 4.8);

        histosMap_["subjet_pt_onGGH"] = new TH1F(("subjet_pt_onGGH_" + systName).c_str(), "Subleading Jet p_{T};p_{T}  (GeV);Events", 45, 25, 250);
        histosMap_["subjet_eta_onGGH"] = new TH1F(("subjet_eta_onGGH_" + systName).c_str(), "Subleading Jet \\eta;\\eta;Events", 24, -4.8, 4.8);

        histosMap_["dijet_mass_onGGH"] = new TH1F(("dijet_mass_onGGH_" + systName).c_str(), "DiJet Mass;M_{jj}  (GeV);Events", 40, 0, 2000);
        histosMap_["dijet_deta_onGGH"] = new TH1F(("dijet_deta_onGGH_" + systName).c_str(), "DiJet deta;deta;Events", 30, 0, 9);
        histosMap_["dijet_dphi_onGGH"] = new TH1F(("dijet_dphi_onGGH_" + systName).c_str(), "DiJet dphi;dphi;Events", 30, 0, 3.14);

        histosMap_["zeppen_onGGH"] = new TH1F(("zeppen_onGGH_" + systName).c_str(), "; zeppen; ", 20, -10, 10);
        histosMap_["dphimmj_onGGH"] = new TH1F(("dphimmj_onGGH_" + systName).c_str(), "dphimmj;dphimmj;Events", 40, 0, 3.14 );
        histosMap_["detammj_onGGH"] = new TH1F(("detammj_onGGH_" + systName).c_str(), "detammj;detammj;Events", 40, 0, 8);

        
        histosMap_["met_pt_onGGH"] = new TH1F(("met_pt_onGGH_" + systName).c_str(), "MET p_{T};p_{T}  (GeV) ", 40, 0, 140);
        histosMap_["met_phi_onGGH"] = new TH1F(("met_phi_onGGH_" + systName).c_str(), "MET phi", 18, -3.6, 3.6);

        histosMap_["csTheta_onGGH"] = new TH1F(("csTheta_onGGH_" + systName).c_str(), " ; csTheta; ", 10, -1, 1);
        histosMap_["csPhi_onGGH"] = new TH1F(("csPhi_onGGH_" + systName).c_str(), " ; csPhi ; ", 30, -3.14, 3.14);

        histosMap_["bdtScore_onGGH"] = new TH1F(("bdtScore_onGGH_" + systName).c_str(), "; ; ", 40, -0.9, 0.9);

        histosMap_["dimuon_0"] = new TH1F(("dimuon_mass_cat0_" + systName).c_str(), " cat 0 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_1"] = new TH1F(("dimuon_mass_cat1_" + systName).c_str(), " cat 1 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_2"] = new TH1F(("dimuon_mass_cat2_" + systName).c_str(), " cat 2 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_3"] = new TH1F(("dimuon_mass_cat3_" + systName).c_str(), " cat 3 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_4"] = new TH1F(("dimuon_mass_cat4_" + systName).c_str(), " cat 4 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_5"] = new TH1F(("dimuon_mass_cat5_" + systName).c_str(), " cat 5 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_6"] = new TH1F(("dimuon_mass_cat6_" + systName).c_str(), " cat 6 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_7"] = new TH1F(("dimuon_mass_cat7_" + systName).c_str(), " cat 7 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_8"] = new TH1F(("dimuon_mass_cat8_" + systName).c_str(), " cat 8 ; dimuon mass (GeV); Events", 400, 110, 150);
    }

    void fillEntry(eventEntry data, float weight)
    {
        float finalWeight = weight * data.btagSF;
        if (data.higgsCandMass > 110 and data.higgsCandMass < 150)
        {
            if (data.category < 5)
            {
                histosMap_["leadMuon_pt_onGGH"]->Fill(data.muonOnePt, finalWeight);
                histosMap_["leadMuon_eta_onGGH"]->Fill(data.muonOneEta, finalWeight);
                histosMap_["leadMuon_phi_onGGH"]->Fill(data.muonOnePhi, finalWeight);
                histosMap_["leadMuon_ptOverM_onGGH"]->Fill( (data.muonOnePt / data.higgsCandMass) , finalWeight);
                
                histosMap_["subMuon_pt_onGGH"]->Fill(data.muonTwoPt, finalWeight);
                histosMap_["subMuon_eta_onGGH"]->Fill(data.muonTwoEta, finalWeight);
                histosMap_["subMuon_phi_onGGH"]->Fill(data.muonTwoPhi, finalWeight);
                histosMap_["subMuon_ptOverM_onGGH"]->Fill( (data.muonTwoPt / data.higgsCandMass) , finalWeight);

                histosMap_["dimuon_mass_onGGH"]->Fill(data.higgsCandMass, finalWeight);
                histosMap_["dimuon_pt_onGGH"]->Fill(data.higgsCandPt, finalWeight);
                histosMap_["dimuon_eta_onGGH"]->Fill(data.higgsCandEta, finalWeight);
                histosMap_["dimuon_phi_onGGH"]->Fill(data.higgsCandPhi, finalWeight);
                histosMap_["dimuon_rap_onGGH"]->Fill(data.higgsCandRap, finalWeight);

                histosMap_["dimuon_deta_onGGH"]->Fill(data.higgsCandDeltaEta, finalWeight);
                histosMap_["dimuon_dphi_onGGH"]->Fill(data.higgsCandDeltaPhi, finalWeight);

                histosMap_["csTheta_onGGH"]->Fill(data.higgsCandCosThetaCs, finalWeight);
                histosMap_["csPhi_onGGH"]->Fill(data.higgsCandPhiCs, finalWeight);

                histosMap_["num_jets_onGGH"]->Fill(data.nJets, finalWeight);
                histosMap_["num_bjets_onGGH"]->Fill(data.nBtagJetsM, finalWeight);

                if (data.nJets > 0)
                {
                    histosMap_["leadjet_pt_onGGH"]->Fill(data.jetOnePt, finalWeight);
                    histosMap_["leadjet_eta_onGGH"]->Fill(data.jetOneEta, finalWeight);

                    if (data.nJets > 1)
                    {
                        histosMap_["subjet_pt_onGGH"]->Fill(data.jetTwoPt, finalWeight);
                        histosMap_["subjet_eta_onGGH"]->Fill(data.jetTwoEta, finalWeight);

                        histosMap_["dijet_mass_onGGH"]->Fill(data.dijetMass, finalWeight);
                        histosMap_["dijet_deta_onGGH"]->Fill(data.dijetDeltaEta, finalWeight);
                        histosMap_["dijet_dphi_onGGH"]->Fill(data.dijetDeltaPhi, finalWeight);

                        histosMap_["zeppen_onGGH"]->Fill(data.zeppen, finalWeight);
                        histosMap_["dphimmj_onGGH"]->Fill(data.deltaPhiHiggsJet, finalWeight);
                        histosMap_["detammj_onGGH"]->Fill(data.deltaEtaHiggsJet, finalWeight);
                    }
                }

                histosMap_["met_pt_onGGH"]->Fill(data.metPT, finalWeight);
                histosMap_["met_phi_onGGH"]->Fill(data.metPhi, finalWeight);

                histosMap_["bdtScore_onGGH"]->Fill(data.bdtScore, finalWeight);

            }


            histosMap_["dimuon_" + std::to_string(data.category)]->Fill(data.higgsCandMass, finalWeight);
        }
    }

    TH1F *getHistogram(std::string histoName)
    {
        return histosMap_[histoName];
    }

    std::vector<TH1F *> getHistograms()
    {
        std::vector<TH1F *> histosVec;

        for (std::map<std::string, TH1F *>::iterator h = histosMap_.begin(); h != histosMap_.end(); ++h)
        {
            histosVec.push_back(h->second);
        }
        return histosVec;
    }
};