struct eventEntry
{

    // int mcLabel;
    // // Event Weights
    // float eventWeight;
    // float pileupWeight;
    // float zPtWeight;

    // // Scale Factors 
    // float prefireSF;
    // float muonIDSF;
    // float muonIsoSF;
    // float muonTrigSF;
    bool isValid;
    // float totalWeight;
    // Observed Variables //
    
    // Leading Muon
    float muonOnePt;
    float muonOneEta;
    float muonOnePhi;

    // Subleading Muon

    float muonTwoPt;
    float muonTwoEta;
    float muonTwoPhi;

    // Higgs Candidate
    float higgsCandMass;
    float higgsCandPt;
    float higgsCandEta;
    float higgsCandPhi;
    float higgsCandRap;
    float higgsCandDeltaEta;
    float higgsCandDeltaPhi;
    float higgsCandCosThetaCs;
    float higgsCandPhiCs;

    // Jet Quantities
    int nJets;
    int nBtagJetsL;
    int nBtagJetsM;

    // Leading Jet
    float jetOnePt;
    float jetOneEta;
    float jetOnePhi;

    // Subleading Jet

    float jetTwoPt;
    float jetTwoEta;
    float jetTwoPhi;

    // Dijet 
    float dijetMass;
    float dijetDeltaEta;
    float dijetDeltaPhi;

    // Higgs + Jets
    float zeppen;
    float deltaPhiHiggsJet;
    float deltaEtaHiggsJet;

    // MET
    float metPT;
    float metPhi;

    // event cateogrization
    float bdtScore;
    int category;
    float btagSF;

};