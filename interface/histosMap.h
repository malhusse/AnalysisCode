#include <TH1.h>
#include "eventEntry.h"

class histosMap {
    private:

    std::map<std::string, TH1F*> histosMap_;

    public:

    histosMap(std::string systName) {
        TH1::SetDefaultSumw2();
        histosMap_["leadMuon_pt_onZ"] = new TH1F( ("lead_muon_pt_onZ_" + systName).c_str(), "Leading Muon p_{T};p_{T}  (GeV);Events", 100, 0, 200);
        histosMap_["leadMuon_eta_onZ"] = new TH1F( ("lead_muon_eta_onZ_" + systName).c_str(), "Leading Muon \\eta;\\eta;Events", 50, -2.5, 2.5);
        histosMap_["leadMuon_phi_onZ"] = new TH1F( ("lead_muon_phi_onZ_" + systName).c_str(), "Leading Muon \\phi;\\phi;Events", 64, -3.2, 3.2);
        histosMap_["subMuon_pt_onZ"] = new TH1F( ("sub_muon_pt_onZ_" + systName).c_str(), "Subleading Muon p_{T};p_{T}  (GeV);Events", 100, 0, 200);
        histosMap_["subMuon_eta_onZ"] = new TH1F( ("sub_muon_eta_onZ_" + systName).c_str(), "Subleading Muon eta;\\eta;Events", 50, -2.5, 2.5);
        histosMap_["subMuon_phi_onZ"] = new TH1F( ("sub_muon_phi_onZ_" + systName).c_str(), "Subleading Muon \\phi;\\phi;Events", 64, -3.2, 3.2);
        histosMap_["dimuon_mass_onZ"] = new TH1F( ("dimuon_mass_onZ_" + systName).c_str(), "Dimuon Mass;M_{\\mu \\mu}  (Gev);Events", 80, 70, 110);
        histosMap_["dimuon_pt_onZ"] = new TH1F( ("dimuon_pt_onZ_" + systName).c_str(), "Dimuon p_{T};p_{T}  (GeV);Events", 120, 0, 240);
        histosMap_["dimuon_eta_onZ"] = new TH1F( ("dimuon_eta_onZ_" + systName).c_str(), "Dimuon \\eta;\\eta;Events", 100, -5.0, 5.0);
        histosMap_["dimuon_phi_onZ"] = new TH1F( ("dimuon_phi_onZ_" + systName).c_str(), "Dimuon \\phi;\\phi;Events", 36, -3.6, 3.6);
        histosMap_["dimuon_deta_onZ"] = new TH1F( ("dimuon_deta_onZ_" + systName).c_str(), "Dimuon deta;deta;Events", 50, 0, 5.0);
        histosMap_["dimuon_dphi_onZ"] = new TH1F( ("dimuon_dphi_onZ_" + systName).c_str(), "Dimuon dphi;dphi;Events", 18, 0, 3.6);
        histosMap_["num_jets_onZ"] = new TH1F( ("num_jets_onZ_" + systName).c_str(), "Number of Jets;nJets;Events", 8, 0, 8);
        histosMap_["num_bjets_onZ"] = new TH1F( ("num_bjets_onZ_" + systName).c_str(), "Number of B Jets;nBJets;Events", 6, 0, 6);
        histosMap_["leadjet_pt_onZ"] = new TH1F( ("leadjet_pt_onZ_" + systName).c_str(), "Leading Jet p_{T};p_{T}  (GeV);Events", 300, 0, 300);
        histosMap_["leadjet_eta_onZ"] = new TH1F( ("leadjet_eta_onZ_" + systName).c_str(), "Leading Jet \\eta;\\eta;Events", 94, -4.7, 4.7);
        histosMap_["subjet_pt_onZ"] = new TH1F( ("subjet_pt_onZ_" + systName).c_str(), "Subleading Jet p_{T};p_{T}  (GeV);Events", 250, 0, 250);
        histosMap_["subjet_eta_onZ"] = new TH1F( ("subjet_eta_onZ_" + systName).c_str(), "Subleading Jet \\eta;\\eta;Events", 94, -4.7, 4.7);
        histosMap_["dijet_mass_onZ"] = new TH1F( ("dijet_mass1_onZ_" + systName).c_str(), "DiJet Mass;M_{jj}  (GeV);Events", 60, 0, 600);
        histosMap_["dijet_deta_onZ"] = new TH1F( ("dijet_deta1_onZ_" + systName).c_str(), "DiJet deta;deta;Events", 94, 0, 9.4);
        histosMap_["met_pt_onZ"] = new TH1F( ("met_pt_onZ_" + systName).c_str(), "MET p_{T};p_{T}  (GeV) ",  100, 0, 200);
        histosMap_["met_phi_onZ"] = new TH1F( ("met_phi_onZ_" + systName).c_str(), "MET phi", 36, -3.6, 3.6);
        histosMap_["zeppen_onZ"] = new TH1F( ("zeppen_onZ_" + systName).c_str(), "; zeppen; ",  100, -10, 10);
        histosMap_["csTheta_onZ"] = new TH1F( ("csTheta_onZ_" + systName).c_str(), " ; csTheta; ",  10, -1, 1);
        histosMap_["csPhi_onZ"] = new TH1F( ("csPhi_onZ_" + systName).c_str(), " ; csPhi ; ",  100, -10, 10);
        histosMap_["bdtScore_onZ"] = new TH1F( ("bdtScore_onZ_" + systName).c_str(), "; ; ",  20, -1, 1);   

        // Inclusive (All Categories) Plots
        histosMap_["leadMuon_pt_onH"] = new TH1F( ("lead_muon_pt_onH_" + systName).c_str(), "Leading Muon p_{T};p_{T}  (GeV);Events", 100, 0, 200);
        histosMap_["leadMuon_eta_onH"] = new TH1F( ("lead_muon_eta_onH_" + systName).c_str(), "Leading Muon \\eta;\\eta;Events", 50, -2.5, 2.5);
        histosMap_["leadMuon_phi_onH"] = new TH1F( ("lead_muon_phi_onH_" + systName).c_str(), "Leading Muon \\phi;\\phi;Events", 64, -3.2, 3.2);
        histosMap_["subMuon_pt_onH"] = new TH1F( ("sub_muon_pt_onH_" + systName).c_str(), "Subleading Muon p_{T};p_{T}  (GeV);Events", 100, 0, 200);
        histosMap_["subMuon_eta_onH"] = new TH1F( ("sub_muon_eta_onH_" + systName).c_str(), "Subleading Muon eta;\\eta;Events", 50, -2.5, 2.5);
        histosMap_["subMuon_phi_onH"] = new TH1F( ("sub_muon_phi_onH_" + systName).c_str(), "Subleading Muon \\phi;\\phi;Events", 64, -3.2, 3.2);
        histosMap_["dimuon_mass_onH"] = new TH1F( ("dimuon_mass_onH_" + systName).c_str(), "Dimuon Mass;M_{\\mu \\mu}  (Gev);Events", 80, 110, 150);
        histosMap_["dimuon_pt_onH"] = new TH1F( ("dimuon_pt_onH_" + systName).c_str(), "Dimuon p_{T};p_{T}  (GeV);Events", 200, 0, 400);
        histosMap_["dimuon_eta_onH"] = new TH1F( ("dimuon_eta_onH_" + systName).c_str(), "Dimuon \\eta;\\eta;Events", 100, -5.0, 5.0);
        histosMap_["dimuon_phi_onH"] = new TH1F( ("dimuon_phi_onH_" + systName).c_str(), "Dimuon \\phi;\\phi;Events", 36, -3.6, 3.6);
        histosMap_["dimuon_deta_onH"] = new TH1F( ("dimuon_deta_onH_" + systName).c_str(), "Dimuon deta;deta;Events", 50, 0, 5.0);
        histosMap_["dimuon_dphi_onH"] = new TH1F( ("dimuon_dphi_onH_" + systName).c_str(), "Dimuon dphi;dphi;Events", 18, 0, 3.6);
        histosMap_["num_jets_onH"] = new TH1F( ("num_jets_onH_" + systName).c_str(), "Number of Jets;nJets;Events", 8, 0, 8);
        histosMap_["num_bjets_onH"] = new TH1F( ("num_bjets_onH_" + systName).c_str(), "Number of B Jets;nBJets;Events", 6, 0, 6);
        histosMap_["leadjet_pt_onH"] = new TH1F( ("leadjet_pt_onH_" + systName).c_str(), "Leading Jet p_{T};p_{T}  (GeV);Events", 300, 0, 300);
        histosMap_["leadjet_eta_onH"] = new TH1F( ("leadjet_eta_onH_" + systName).c_str(), "Leading Jet \\eta;\\eta;Events", 94, -4.7, 4.7);
        histosMap_["subjet_pt_onH"] = new TH1F( ("subjet_pt_onH_" + systName).c_str(), "Subleading Jet p_{T};p_{T}  (GeV);Events", 250, 0, 250);
        histosMap_["subjet_eta_onH"] = new TH1F( ("subjet_eta_onH_" + systName).c_str(), "Subleading Jet \\eta;\\eta;Events", 94, -4.7, 4.7);
        histosMap_["dijet_mass_onH"] = new TH1F( ("dijet_mass1_onH_" + systName).c_str(), "DiJet Mass;M_{jj}  (GeV);Events", 60, 0, 600);
        histosMap_["dijet_deta_onH"] = new TH1F( ("dijet_deta1_onH_" + systName).c_str(), "DiJet deta;deta;Events", 94, 0, 9.4);
        histosMap_["met_pt_onH"] = new TH1F( ("met_pt_onH_" + systName).c_str(), "MET p_{T};p_{T}  (GeV) ",  100, 0, 200);
        histosMap_["met_phi_onH"] = new TH1F( ("met_phi_onH_" + systName).c_str(), "MET phi", 36, -3.6, 3.6);
        histosMap_["zeppen_onH"] = new TH1F( ("zeppen_onH_" + systName).c_str(), "; zeppen; ",  100, -10, 10);
        histosMap_["csTheta_onH"] = new TH1F( ("csTheta_onH_" + systName).c_str(), " ; csTheta; ",  10, -1, 1);
        histosMap_["csPhi_onH"] = new TH1F( ("csPhi_onH_" + systName).c_str(), " ; csPhi ; ",  100, -10, 10);
        histosMap_["bdtScore_onH"] = new TH1F( ("bdtScore_onH_" + systName).c_str(), "; ; ",  20, -1, 1);   

        // ggH-enriched Category Plots
        histosMap_["leadMuon_pt_onGGH"] = new TH1F( ("lead_muon_pt_onGGH_" + systName).c_str(), "Leading Muon p_{T};p_{T}  (GeV);Events", 100, 0, 200);
        histosMap_["leadMuon_eta_onGGH"] = new TH1F( ("lead_muon_eta_onGGH_" + systName).c_str(), "Leading Muon \\eta;\\eta;Events", 50, -2.5, 2.5);
        histosMap_["leadMuon_phi_onGGH"] = new TH1F( ("lead_muon_phi_onGGH_" + systName).c_str(), "Leading Muon \\phi;\\phi;Events", 64, -3.2, 3.2);
        histosMap_["subMuon_pt_onGGH"] = new TH1F( ("sub_muon_pt_onGGH_" + systName).c_str(), "Subleading Muon p_{T};p_{T}  (GeV);Events", 100, 0, 200);
        histosMap_["subMuon_eta_onGGH"] = new TH1F( ("sub_muon_eta_onGGH_" + systName).c_str(), "Subleading Muon eta;\\eta;Events", 50, -2.5, 2.5);
        histosMap_["subMuon_phi_onGGH"] = new TH1F( ("sub_muon_phi_onGGH_" + systName).c_str(), "Subleading Muon \\phi;\\phi;Events", 64, -3.2, 3.2);
        histosMap_["dimuon_mass_onGGH"] = new TH1F( ("dimuon_mass_onGGH_" + systName).c_str(), "Dimuon Mass;M_{\\mu \\mu}  (Gev);Events", 80, 110, 150);
        histosMap_["dimuon_pt_onGGH"] = new TH1F( ("dimuon_pt_onGGH_" + systName).c_str(), "Dimuon p_{T};p_{T}  (GeV);Events", 200, 0, 400);
        histosMap_["dimuon_eta_onGGH"] = new TH1F( ("dimuon_eta_onGGH_" + systName).c_str(), "Dimuon \\eta;\\eta;Events", 100, -5.0, 5.0);
        histosMap_["dimuon_phi_onGGH"] = new TH1F( ("dimuon_phi_onGGH_" + systName).c_str(), "Dimuon \\phi;\\phi;Events", 36, -3.6, 3.6);
        histosMap_["dimuon_deta_onGGH"] = new TH1F( ("dimuon_deta_onGGH_" + systName).c_str(), "Dimuon deta;deta;Events", 50, 0, 5.0);
        histosMap_["dimuon_dphi_onGGH"] = new TH1F( ("dimuon_dphi_onGGH_" + systName).c_str(), "Dimuon dphi;dphi;Events", 18, 0, 3.6);
        histosMap_["num_jets_onGGH"] = new TH1F( ("num_jets_onGGH_" + systName).c_str(), "Number of Jets;nJets;Events", 8, 0, 8);
        histosMap_["num_bjets_onGGH"] = new TH1F( ("num_bjets_onGGH_" + systName).c_str(), "Number of B Jets;nBJets;Events", 6, 0, 6);
        histosMap_["leadjet_pt_onGGH"] = new TH1F( ("leadjet_pt_onGGH_" + systName).c_str(), "Leading Jet p_{T};p_{T}  (GeV);Events", 300, 0, 300);
        histosMap_["leadjet_eta_onGGH"] = new TH1F( ("leadjet_eta_onGGH_" + systName).c_str(), "Leading Jet \\eta;\\eta;Events", 94, -4.7, 4.7);
        histosMap_["subjet_pt_onGGH"] = new TH1F( ("subjet_pt_onGGH_" + systName).c_str(), "Subleading Jet p_{T};p_{T}  (GeV);Events", 250, 0, 250);
        histosMap_["subjet_eta_onGGH"] = new TH1F( ("subjet_eta_onGGH_" + systName).c_str(), "Subleading Jet \\eta;\\eta;Events", 94, -4.7, 4.7);
        histosMap_["dijet_mass_onGGH"] = new TH1F( ("dijet_mass1_onGGH_" + systName).c_str(), "DiJet Mass;M_{jj}  (GeV);Events", 60, 0, 600);
        histosMap_["dijet_deta_onGGH"] = new TH1F( ("dijet_deta1_onGGH_" + systName).c_str(), "DiJet deta;deta;Events", 94, 0, 9.4);
        histosMap_["met_pt_onGGH"] = new TH1F( ("met_pt_onGGH_" + systName).c_str(), "MET p_{T};p_{T}  (GeV) ",  100, 0, 200);
        histosMap_["met_phi_onGGH"] = new TH1F( ("met_phi_onGGH_" + systName).c_str(), "MET phi", 36, -3.6, 3.6);
        histosMap_["zeppen_onGGH"] = new TH1F( ("zeppen_onGGH_" + systName).c_str(), "; zeppen; ",  100, -10, 10);
        histosMap_["csTheta_onGGH"] = new TH1F( ("csTheta_onGGH_" + systName).c_str(), " ; csTheta; ",  10, -1, 1);
        histosMap_["csPhi_onGGH"] = new TH1F( ("csPhi_onGGH_" + systName).c_str(), " ; csPhi ; ",  100, -10, 10);
        histosMap_["bdtScore_onGGH"] = new TH1F( ("bdtScore_onGGH_" + systName).c_str(), "; ; ",  20, -1, 1);   

        histosMap_["dimuon_0"] = new TH1F( ("dimuon_mass_cat0_" + systName).c_str(), " cat 0 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_1"] = new TH1F( ("dimuon_mass_cat1_" + systName).c_str(), " cat 1 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_2"] = new TH1F( ("dimuon_mass_cat2_" + systName).c_str(), " cat 2 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_3"] = new TH1F( ("dimuon_mass_cat3_" + systName).c_str(), " cat 3 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_4"] = new TH1F( ("dimuon_mass_cat4_" + systName).c_str(), " cat 4 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_5"] = new TH1F( ("dimuon_mass_cat5_" + systName).c_str(), " cat 5 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_6"] = new TH1F( ("dimuon_mass_cat6_" + systName).c_str(), " cat 6 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_7"] = new TH1F( ("dimuon_mass_cat7_" + systName).c_str(), " cat 7 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_8"] = new TH1F( ("dimuon_mass_cat8_" + systName).c_str(), " cat 8 ; dimuon mass (GeV); Events", 400, 110, 150);
        histosMap_["dimuon_9"] = new TH1F( ("dimuon_mass_cat9_" + systName).c_str(), " cat 9 ; dimuon mass (GeV); Events", 400, 110, 150);
    }

    void fillEntry(eventEntry data, float weight) {
        if (data.higgsCandMass > 70 && data.higgsCandMass < 110) {
            //std::cout << "filling z event " << std::endl;
            histosMap_["leadMuon_pt_onZ"]->Fill( data.muonOnePt ,weight);
            histosMap_["leadMuon_eta_onZ"]->Fill( data.muonOneEta ,weight);
            histosMap_["leadMuon_phi_onZ"]->Fill( data.muonOnePhi ,weight);

            histosMap_["subMuon_pt_onZ"]->Fill( data.muonTwoPt ,weight);
            histosMap_["subMuon_eta_onZ"]->Fill( data.muonTwoEta ,weight);
            histosMap_["subMuon_phi_onZ"]->Fill( data.muonTwoPhi ,weight);

            histosMap_["dimuon_mass_onZ"]->Fill( data.higgsCandMass ,weight);
            histosMap_["dimuon_pt_onZ"]->Fill( data.higgsCandPt ,weight);
            histosMap_["dimuon_eta_onZ"]->Fill( data.higgsCandEta ,weight);
            histosMap_["dimuon_phi_onZ"]->Fill( data.higgsCandPhi ,weight);
            histosMap_["dimuon_deta_onZ"]->Fill( data.higgsCandDeltaEta ,weight);
            histosMap_["dimuon_dphi_onZ"]->Fill( data.higgsCandDeltaPhi ,weight);

            histosMap_["csTheta_onZ"]->Fill( data.higgsCandCosThetaCs ,weight);
            histosMap_["csPhi_onZ"]->Fill( data.higgsCandPhiCs ,weight);

            histosMap_["num_jets_onZ"]->Fill( data.nJets ,weight);
            histosMap_["num_bjets_onZ"]->Fill( data.nBtagJetsM ,weight);
            
            if (data.nJets > 0) {
                histosMap_["leadjet_pt_onZ"]->Fill( data.jetOnePt ,weight);
                histosMap_["leadjet_eta_onZ"]->Fill( data.jetOneEta ,weight);

                if (data.nJets > 1) {
                    histosMap_["subjet_pt_onZ"]->Fill( data.jetTwoPt ,weight);
                    histosMap_["subjet_eta_onZ"]->Fill( data.jetTwoEta ,weight);
                    histosMap_["dijet_mass_onZ"]->Fill( data.dijetMass ,weight);
                    histosMap_["dijet_deta_onZ"]->Fill( data.dijetDeltaEta ,weight);
                    histosMap_["zeppen_onZ"]->Fill( data.zeppen ,weight);
                }
            }

            histosMap_["met_pt_onZ"]->Fill( data.metPT ,weight);
            histosMap_["met_phi_onZ"]->Fill( data.metPhi ,weight);

            histosMap_["bdtScore_onZ"]->Fill( data.bdtScore ,weight);
            //std::cout << "filled z event success" << std::endl;

        }
        
        if (data.higgsCandMass > 110 && data.higgsCandMass < 150) {
            //std::cout << "filling H event " << std::endl;

            histosMap_["leadMuon_pt_onH"]->Fill( data.muonOnePt ,weight);
            histosMap_["leadMuon_eta_onH"]->Fill( data.muonOneEta ,weight);
            histosMap_["leadMuon_phi_onH"]->Fill( data.muonOnePhi ,weight);

            histosMap_["subMuon_pt_onH"]->Fill( data.muonTwoPt ,weight);
            histosMap_["subMuon_eta_onH"]->Fill( data.muonTwoEta ,weight);
            histosMap_["subMuon_phi_onH"]->Fill( data.muonTwoPhi ,weight);

            histosMap_["dimuon_mass_onH"]->Fill( data.higgsCandMass ,weight);
            histosMap_["dimuon_pt_onH"]->Fill( data.higgsCandPt ,weight);
            histosMap_["dimuon_eta_onH"]->Fill( data.higgsCandEta ,weight);
            histosMap_["dimuon_phi_onH"]->Fill( data.higgsCandPhi ,weight);
            histosMap_["dimuon_deta_onH"]->Fill( data.higgsCandDeltaEta ,weight);
            histosMap_["dimuon_dphi_onH"]->Fill( data.higgsCandDeltaPhi ,weight);

            histosMap_["csTheta_onH"]->Fill( data.higgsCandCosThetaCs ,weight);
            histosMap_["csPhi_onH"]->Fill( data.higgsCandPhiCs ,weight);

            histosMap_["num_jets_onH"]->Fill( data.nJets ,weight);
            histosMap_["num_bjets_onH"]->Fill( data.nBtagJetsM ,weight);

            if (data.nJets > 0) {
                histosMap_["leadjet_pt_onH"]->Fill( data.jetOnePt ,weight);
                histosMap_["leadjet_eta_onH"]->Fill( data.jetOneEta ,weight);
            
                if (data.nJets > 1) {
                    histosMap_["subjet_pt_onH"]->Fill( data.jetTwoPt ,weight);
                    histosMap_["subjet_eta_onH"]->Fill( data.jetTwoEta ,weight);
                    histosMap_["dijet_mass_onH"]->Fill( data.dijetMass ,weight);
                    histosMap_["dijet_deta_onH"]->Fill( data.dijetDeltaEta ,weight);
                    histosMap_["zeppen_onH"]->Fill( data.zeppen ,weight);
                }
            }
            
            histosMap_["met_pt_onH"]->Fill( data.metPT ,weight);
            histosMap_["met_phi_onH"]->Fill( data.metPhi ,weight);

            histosMap_["bdtScore_onH"]->Fill( data.bdtScore ,weight);

            //std::cout << "filled H event success " << std::endl;

            if ( data.category < 5)
            {
                //std::cout << "filling ggH event " << std::endl;

                histosMap_["leadMuon_pt_onGGH"]->Fill( data.muonOnePt ,weight);
                histosMap_["leadMuon_eta_onGGH"]->Fill( data.muonOneEta ,weight);
                histosMap_["leadMuon_phi_onGGH"]->Fill( data.muonOnePhi ,weight);

                histosMap_["subMuon_pt_onGGH"]->Fill( data.muonTwoPt ,weight);
                histosMap_["subMuon_eta_onGGH"]->Fill( data.muonTwoEta ,weight);
                histosMap_["subMuon_phi_onGGH"]->Fill( data.muonTwoPhi ,weight);

                histosMap_["dimuon_mass_onGGH"]->Fill( data.higgsCandMass ,weight);
                histosMap_["dimuon_pt_onGGH"]->Fill( data.higgsCandPt ,weight);
                histosMap_["dimuon_eta_onGGH"]->Fill( data.higgsCandEta ,weight);
                histosMap_["dimuon_phi_onGGH"]->Fill( data.higgsCandPhi ,weight);
                histosMap_["dimuon_deta_onGGH"]->Fill( data.higgsCandDeltaEta ,weight);
                histosMap_["dimuon_dphi_onGGH"]->Fill( data.higgsCandDeltaPhi ,weight);

                histosMap_["csTheta_onGGH"]->Fill( data.higgsCandCosThetaCs ,weight);
                histosMap_["csPhi_onGGH"]->Fill( data.higgsCandPhiCs ,weight);

                histosMap_["num_jets_onGGH"]->Fill( data.nJets ,weight);
                histosMap_["num_bjets_onGGH"]->Fill( data.nBtagJetsM ,weight);

                if (data.nJets > 0) {
                    histosMap_["leadjet_pt_onGGH"]->Fill( data.jetOnePt ,weight);
                    histosMap_["leadjet_eta_onGGH"]->Fill( data.jetOneEta ,weight);
                
                    if (data.nJets > 1) {
                        histosMap_["subjet_pt_onGGH"]->Fill( data.jetTwoPt ,weight);
                        histosMap_["subjet_eta_onGGH"]->Fill( data.jetTwoEta ,weight);
                        histosMap_["dijet_mass_onGGH"]->Fill( data.dijetMass ,weight);
                        histosMap_["dijet_deta_onGGH"]->Fill( data.dijetDeltaEta ,weight);
                        histosMap_["zeppen_onGGH"]->Fill( data.zeppen ,weight);
                    }
                }
                
                histosMap_["met_pt_onGGH"]->Fill( data.metPT ,weight);
                histosMap_["met_phi_onGGH"]->Fill( data.metPhi ,weight);

                histosMap_["bdtScore_onGGH"]->Fill( data.bdtScore ,weight);

            //std::cout << "filled ggH event success" << std::endl;

            } 

            //std::cout << "filling dimuon histo " << std::endl;
            
            histosMap_["dimuon_" + std::to_string(data.category)]->Fill(data.higgsCandMass, weight);
            //std::cout << "filled dimuon histo " << std::endl;
            
        }

    }

    std::vector<TH1F*> getHistograms() {
        std::vector<TH1F*> histosVec;
        
        for ( std::map< std::string, TH1F*>::iterator h = histosMap_.begin() ; h != histosMap_.end(); ++h)
        {
            histosVec.push_back( h->second );
        }
        return histosVec;
    }
};