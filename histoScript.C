void histoScript(string inputFile, Int_t year)
{
    // std::map<int, float> lumiMap;
    //     lumiMap[2016] = 35545.5 ;
    //     lumiMap[2017] = 41859.5 ;
    //     lumiMap[2018] = 58877.4 ;

    string file = "ntupleFiles/";
    file += std::to_string(year);
    file += "/";
    file += inputFile;

   TProof::Open("workers=4");

    TChain *ch = new TChain("ntupledData");
    cout << file << endl;
    ch->Add(file.c_str());

    gProof->AddInput(new TNamed("outputName", inputFile));
    gProof->SetParameter("getYear",year);

    ch->SetProof();
    ch->Process("plotSelec.C+");

}
