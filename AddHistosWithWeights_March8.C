void AddHistoWithWeight_March8(){
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptTitle(kFALSE);                                                                                                                                                                             
  gStyle->SetLegendBorderSize(0);
  
  float NBins = 100;
  float BinMin = 0;
  float BinMax = 1000;//5000

  TString s_XTitle;

  TString dir = "/uscms_data/d3/3wolf3/ttH/Jan2019/Mar2019/CMSSW_9_4_9/src/Ragnarok/Output_Histograms/";//"/eos/uscms/store/user/3wolf3/";

  TFile *Data_inf = TFile::Open(dir+"Output_Histograms_DoubleEG_MuonEG_DoubleMuon_SingleElectron_SingleMuon_PeriodAll.root");
			//"Output_Histograms_DataAll_PeriodALL_AllDL_AllSF_ptpt.root");

  TString DYJetsToLL_M_10to50 = "Output_Histograms_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root";
  TString DYJetsToLL_M_50 = "Output_Histograms_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root";
  TString ST_s_channel_4f = "Output_Histograms_ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8.root";
  TString ST_t_channel_top_4f = "Output_Histograms_ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root";
  TString ST_tW_top_5f = "Output_Histograms_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root";
  TString TTToHadronic = "Output_Histograms_TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root";
  TString TTToSemiLept = "Output_Histograms_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root";
  TString TTTo2L2Nu = "Output_Histograms_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root";
  TString TTGJets = "Output_Histograms_TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root";
  TString TTWJetsToLNu = "Output_Histograms_TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root";
  TString TTWJetsToQQ = "Output_Histograms_TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root";
  TString TTZToLLNuNu= "Output_Histograms_TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root";
  TString TTZToQQ = "Output_Histograms_TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8.root";
  TString WJetsToLNu = "Output_Histograms_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.root";
  TString ttHTobb = "Output_Histograms_ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8.root";

  TFile *DYJetsToLL_M_10to50_inf = TFile::Open(dir+DYJetsToLL_M_10to50);
  TFile *DYJetsToLL_M_50_inf = TFile::Open(dir+DYJetsToLL_M_50);
  TFile *ST_s_channel_4f_inf = TFile::Open(dir+ST_s_channel_4f);  
  TFile *ST_t_channel_top_4f_inf = TFile::Open(dir+ST_t_channel_top_4f);
  TFile *ST_tW_top_5f_inf = TFile::Open(dir+ST_tW_top_5f);
  TFile *TTToHadronic_inf = TFile::Open(dir+TTToHadronic);
  TFile *TTToSemiLept_inf = TFile::Open(dir+TTToSemiLept);
  TFile *TTTo2L2Nu_inf = TFile::Open(dir+TTTo2L2Nu);
  TFile *TTGJets_inf = TFile::Open(dir+TTGJets);
  TFile *TTWJetsToLNu_inf = TFile::Open(dir+TTWJetsToLNu);
  TFile *TTWJetsToQQ_inf = TFile::Open(dir+TTWJetsToQQ);
  TFile *TTZToLLNuNu_inf = TFile::Open(dir+TTZToLLNuNu);
  TFile *TTZToQQ_inf = TFile::Open(dir+TTZToQQ);
  TFile *WJetsToLNu_inf = TFile::Open(dir+WJetsToLNu);
  TFile *ttHTobb_inf = TFile::Open(dir+ttHTobb);


  float ScaleToFit = 1.0; //Still trying to figure out why NJets is so much off....
  
  
  int PlotLog = 1;
  
  //int withFF = 1;	
  int withFF = 0;
  int NoEffNorm = 0;
  //int NoEffNorm = 1;
//=========================================
  TString Channel = "";
  //TString Channel = "_ee";
  //TString Channel = "_emu";
  //TString Channel = "_mumu";
//========================================
  //TString Var = "h_Jet0DeepCSV_b";
  //TString Var = "h_Jet0PT_SF";
  //TString Var = "h_Jet1PT";

  //TString Var = "h_NumberOfJets";     // WORKS
  //TString Var = "h_NumberOfJets";
  //TString Var = "h_NumberOfJets_SF";
  //TString Var = "h_bdt_score_3j2b_SF";
  //TString Var = "h_bdt_score_3j3b_SF";
  //TString Var = "h_bdt_score_ge4j2b_SF";
  TString Var = "h_bdt_score_ge4j3b_SF";
  //TString Var = "h_MET";


  //TString Var = "h_NumberOfbTags";//WORKS
  //TString Var = "h_NumberOfbTags_SF";//WORKS
  //TString Var = "h_NumberOfPV";//WORKS
  //TString Var = "h_NumberOfLeptons";//WORKS but need to change EventSelection to bin from 1-2
  //TString Var = "h_DeepCSVSF"; //No data in plot
 //TString Var = "h_SumJetPT";//Works but looks  too dense
  //TString Var = "h_MET";//BROKE, something up with TTToHadronic
  //TString Var = "h_Lep0PT";//WORKS
  //TString Var = "h_Lep0PT_SF";
  //TString Var = "h_Lep0Eta";//WORKS
  //TString Var = "h_Lep0Phi";//DOES NOT WORK
  //TString Var = "h_Lep0Iso";//DOES NOT WORK because it's blank
  //TString Var = "h_Lep1PT";
  //TString Var = "h_Lep1Eta";//Didnt have it write out in EventReadout....
  //TString Var = "h_Ele0PT";//Didnt write it out either.......
  //TString Var = "h_Jet0PT_ttlf";//WORKS
  //TString Var = "h_Jet0PT_ttbb";
  //TString Var = "h_Jet0Eta";//WORKS and looks good to data
  //TString Var = "h_Jet0Phi";//Issue with TTToHadronic not having any events???? WORKS and compares well to data
  //TString Var = "h_Jet0DeepCSV_b";//WORKS but seems a bit off....
  //TString Var = "h_Jet1DeepCSV_b";
  //TString Var = "h_Jet0DeepCSV_b_SF";//DOESNT WORK FOR DATA.....WORKS looks good
  //TString Var = "h_Jet1PT";//WORKS
  //TString Var = "h_Jet1Eta";//DOES NOT WORK
  //TString Var = "h_Jet1Phi";//WORKS
  //TString Var = "h_Jet1DeepCSV_b";



  std::string CheckNumberOf("NumberOf");
  std::string CheckbTags("bTags");
  std::string CheckLeading ("0");
  std::string CheckSubLeading ("1");
  std::string CheckLepton ("Lep");
  std::string CheckJet ("Jet");
  std::string CheckTopMass ("TopMass");
  std::string CheckWMass ("WMass");
  std::string CheckPT ("PT");
  std::string CheckHT ("HT");
  std::string CheckBT ("BT");
  std::string CheckMet ("Met");
  std::string CheckEta ("Eta");  
  std::string CheckPhi ("Phi");
  std::string CheckDeltaR ("DeltaR");
  std::string CheckMCMatched ("MCMatched");  
  std::string CheckDeepCSV("DeepCSV");
  std::string Check_ttlf("_ttlf");
  std::string Check_ttb("_ttb");
  std::string Check_ttbb("_ttbb");
  std::string Check_tt2b("_tt2b");
  std::string Check_ttcc("_ttcc");


  std::string CheckVar (Var);
  std::size_t foundNumberOf = CheckVar.find(CheckNumberOf);
  std::size_t foundbTags = CheckVar.find(CheckbTags);
  std::size_t foundLeading = CheckVar.find(CheckLeading);
  std::size_t foundSubLeading = CheckVar.find(CheckSubLeading);
  std::size_t foundLepton = CheckVar.find(CheckLepton);
  std::size_t foundJet = CheckVar.find(CheckJet);
  std::size_t foundTopMass = CheckVar.find(CheckTopMass);
  std::size_t foundWMass = CheckVar.find(CheckWMass);
  std::size_t foundPT = CheckVar.find(CheckPT);
  std::size_t foundHT = CheckVar.find(CheckHT);
  std::size_t foundBT = CheckVar.find(CheckBT);
  std::size_t foundMet = CheckVar.find(CheckMet);
  std::size_t foundEta = CheckVar.find(CheckEta);
  std::size_t foundPhi = CheckVar.find(CheckPhi);
  std::size_t foundDeltaR = CheckVar.find(CheckDeltaR);
  std::size_t foundMCMatched = CheckVar.find(CheckMCMatched);
  std::size_t foundDeepCSV = CheckVar.find(CheckDeepCSV);
  std::size_t found_ttlf = CheckVar.find(Check_ttlf);
  std::size_t found_ttb = CheckVar.find(Check_ttb);
  std::size_t found_ttbb = CheckVar.find(Check_ttbb);
  std::size_t found_tt2b = CheckVar.find(Check_tt2b);
  std::size_t found_ttcc = CheckVar.find(Check_ttcc);

  cout << "++++++++++++++++++++++++++++++++++" << endl;
  cout << "Running over " << Var << endl;
  cout << "++++++++++++++++++++++++++++++++++"<< endl;

  //Making sure if we do MC matching that data uses the corresponding variable, since it doesn't have a truth value...
/*
   TH1F *HtData;
  TString VarData;
  if(foundMCMatched <= 100){
	CheckVar.erase(foundMCMatched,9);
	VarData = CheckVar;
        cout << "$$$$$$$$$$$$$$$$$$$$$$$$$" << endl << VarData << "$$$$$$$$$$$$$$$$$$$$" << endl;
	HtData = (TH1F*)Data_inf->Get(VarData);
  }
  else HtData = (TH1F*)Data_inf->Get(Var);
*/
  TH1F *HtData = (TH1F*)Data_inf->Get(Var+Channel);
/*
  float maxYValue = HtData->GetMaximum();
  if (PlotLog == 1)maxYValue = maxYValue * 50;
  else maxYValue = maxYValue * 1.6;
  cout << "GetMaximum: " << maxYValue << endl;
  HtData->GetYaxis()->SetRange(0,maxYValue*2);
*/
  //int NumberOfBins = HtData->GetSize() - 2;

  BinMin = HtData->GetXaxis()->GetBinLowEdge(1);

  cout << "BinMin: "<< BinMin << endl;

  float BinW = HtData->GetBinWidth(2);
  int NumberBins = HtData->GetNbinsX();
  NBins = NumberBins;//DID THIS WORK FOR Var=NJets?????
  
  cout << "NumberBins: " << NumberBins << endl;

  BinMax = BinMin + (BinW*NumberBins);
  //BinMax = 1000;
  if(foundEta <=100 ){
  	BinMin = -2.4;
  	BinMax = 2.4;
 	NBins = 60;  
   }
  if(foundPhi <= 100){
  	BinMin = -4;
	BinMax = 4;
	NBins = 40;
  }


 TH1F *HtVJets = new TH1F("HtVJets","VJets",NBins,BinMin,BinMax);
 TH1F *HtTTV = new TH1F("HtTTV","TTV",NBins,BinMin,BinMax);
 TH1F *HtST = new TH1F("HtST","ST",NBins,BinMin,BinMax);
 TH1F *HtDiboson = new TH1F("HtDiboson","Diboson",NBins,BinMin,BinMax);
 TH1F *Ht_ttlf = new TH1F("Ht_ttlf","ttlf",NBins,BinMin,BinMax);
 TH1F *Ht_ttb = new TH1F("Ht_ttb","ttb",NBins,BinMin,BinMax);
 TH1F *Ht_ttbb = new TH1F("Ht_ttbb","ttbb",NBins,BinMin,BinMax);
 TH1F *Ht_tt2b = new TH1F("Ht_tt2b","tt2b",NBins,BinMin,BinMax);
 TH1F *Ht_ttcc = new TH1F("Ht_ttcc","ttcc",NBins,BinMin,BinMax);


  TH1F *HtDYJetsToLL_M_10to50 = (TH1F*)DYJetsToLL_M_10to50_inf->Get(Var+Channel);
  TH1F *HtDYJetsToLL_M_50 = (TH1F*)DYJetsToLL_M_50_inf->Get(Var+Channel);
  TH1F *HtST_s_channel_4f = (TH1F*)ST_s_channel_4f_inf->Get(Var+Channel);
  TH1F *HtST_t_channel_top_4f = (TH1F*)ST_t_channel_top_4f_inf->Get(Var+Channel);
  TH1F *HtST_tW_top_5f = (TH1F*)ST_tW_top_5f_inf->Get(Var+Channel);
   TH1F *HtTTToHadronic = (TH1F*)TTToHadronic_inf->Get(Var+Channel);
   TH1F *HtTTToHadronic_ttlf = (TH1F*)TTToHadronic_inf->Get(Var+"_ttlf"+Channel);
  TH1F *HtTTToHadronic_ttb = (TH1F*)TTToHadronic_inf->Get(Var+"_ttb"+Channel);
  TH1F *HtTTToHadronic_ttbb = (TH1F*)TTToHadronic_inf->Get(Var+"_ttbb"+Channel);
   TH1F *HtTTToHadronic_tt2b = (TH1F*)TTToHadronic_inf->Get(Var+"_tt2b"+Channel);
   TH1F *HtTTToHadronic_ttcc = (TH1F*)TTToHadronic_inf->Get(Var+"_ttcc"+Channel);

  TH1F *HtTTToSemiLept = (TH1F*)TTToSemiLept_inf->Get(Var+Channel);
     TH1F *HtTTToSemiLept_ttlf = (TH1F*)TTToSemiLept_inf->Get(Var+"_ttlf"+Channel);
  TH1F *HtTTToSemiLept_ttb = (TH1F*)TTToSemiLept_inf->Get(Var+"_ttb"+Channel);
  TH1F *HtTTToSemiLept_ttbb = (TH1F*)TTToSemiLept_inf->Get(Var+"_ttbb"+Channel);
   TH1F *HtTTToSemiLept_tt2b = (TH1F*)TTToSemiLept_inf->Get(Var+"_tt2b"+Channel);
   TH1F *HtTTToSemiLept_ttcc = (TH1F*)TTToSemiLept_inf->Get(Var+"_ttcc"+Channel);

  TH1F *HtTTTo2L2Nu = (TH1F*)TTTo2L2Nu_inf->Get(Var+Channel);
     TH1F *HtTTTo2L2Nu_ttlf = (TH1F*)TTTo2L2Nu_inf->Get(Var+"_ttlf"+Channel);
  TH1F *HtTTTo2L2Nu_ttb = (TH1F*)TTTo2L2Nu_inf->Get(Var+"_ttb"+Channel);
  TH1F *HtTTTo2L2Nu_ttbb = (TH1F*)TTTo2L2Nu_inf->Get(Var+"_ttbb"+Channel);
   TH1F *HtTTTo2L2Nu_tt2b = (TH1F*)TTTo2L2Nu_inf->Get(Var+"_tt2b"+Channel);
   TH1F *HtTTTo2L2Nu_ttcc = (TH1F*)TTTo2L2Nu_inf->Get(Var+"_ttcc"+Channel);



  TH1F *HtTTGJets = (TH1F*)TTGJets_inf->Get(Var+Channel);
  TH1F *HtTTWJetsToLNu = (TH1F*)TTWJetsToLNu_inf->Get(Var+Channel);
  TH1F *HtTTWJetsToQQ = (TH1F*)TTWJetsToQQ_inf->Get(Var+Channel);
  TH1F *HtTTZToLLNuNu = (TH1F*)TTZToLLNuNu_inf->Get(Var+Channel);
  TH1F *HtTTZToQQ = (TH1F*)TTZToQQ_inf->Get(Var+Channel);
  TH1F *HtWJetsToLNu = (TH1F*)WJetsToLNu_inf->Get(Var+Channel);
  TH1F *HtttHTobb = (TH1F*)ttHTobb_inf->Get(Var+Channel);
  
   
  float BinUnit = (BinMax-BinMin)/NBins;


  //Get number of events from the NJets histogram
  TH1F *cutFlowDYJetsToLL_M_10to50 = (TH1F*)DYJetsToLL_M_10to50_inf->Get("h_cutFlow");
  TH1F *cutFlowDYJetsToLL_M_50 = (TH1F*)DYJetsToLL_M_50_inf->Get("h_cutFlow");
  TH1F *cutFlowST_s_channel_4f = (TH1F*)ST_s_channel_4f_inf->Get("h_cutFlow");
  TH1F *cutFlowST_t_channel_top_4f = (TH1F*)ST_t_channel_top_4f_inf->Get("h_cutFlow");
  TH1F *cutFlowST_tW_top_5f = (TH1F*)ST_tW_top_5f_inf->Get("h_cutFlow");
   TH1F *cutFlowTTToHadronic = (TH1F*)TTToHadronic_inf->Get("h_cutFlow");
  TH1F *cutFlowTTToSemiLept = (TH1F*)TTToSemiLept_inf->Get("h_cutFlow");
  TH1F *cutFlowTTTo2L2Nu = (TH1F*)TTTo2L2Nu_inf->Get("h_cutFlow");
  TH1F *cutFlowTTGJets = (TH1F*)TTGJets_inf->Get("h_cutFlow");
  TH1F *cutFlowTTWJetsToLNu = (TH1F*)TTWJetsToLNu_inf->Get("h_cutFlow");
  TH1F *cutFlowTTWJetsToQQ = (TH1F*)TTWJetsToQQ_inf->Get("h_cutFlow");
  TH1F *cutFlowTTZToLLNuNu = (TH1F*)TTZToLLNuNu_inf->Get("h_cutFlow");
  TH1F *cutFlowTTZToQQ = (TH1F*)TTZToQQ_inf->Get("h_cutFlow");
  TH1F *cutFlowWJetsToLNu = (TH1F*)WJetsToLNu_inf->Get("h_cutFlow");
  TH1F *cutFlowttHTobb = (TH1F*)ttHTobb_inf->Get("h_cutFlow");






//+++++++++++++++++++++++++++++++++++
cout << "================================="<< endl;
cout << "Total Number of Entries (h_cutFlow):"<< endl;
  float NEntriesDYJetsToLL_M_10to50 = cutFlowDYJetsToLL_M_10to50->GetBinContent(1);
  float NEntriesDYJetsToLL_M_50 = cutFlowDYJetsToLL_M_50->GetBinContent(1);
  float NEntriesST_s_channel_4f = cutFlowST_s_channel_4f->GetBinContent(1);
  float NEntriesST_t_channel_top_4f = cutFlowST_t_channel_top_4f->GetBinContent(1);
  float NEntriesST_tW_top_5f = cutFlowST_tW_top_5f->GetBinContent(1);
  float NEntriesTTToHadronic = cutFlowTTToHadronic->GetBinContent(1);
  cout << "NEntriesTTToHadronic = " << NEntriesTTToHadronic << endl;
  float NEntriesTTToSemiLept = cutFlowTTToSemiLept->GetBinContent(1);
  float NEntriesTTTo2L2Nu = cutFlowTTTo2L2Nu->GetBinContent(1);
  float NEntriesTTGJets = cutFlowTTGJets->GetBinContent(1);
  float NEntriesTTWJetsToLNu = cutFlowTTWJetsToLNu->GetBinContent(1);
  float NEntriesTTWJetsToQQ = cutFlowTTWJetsToQQ->GetBinContent(1);
  float NEntriesTTZToLLNuNu = cutFlowTTZToLLNuNu->GetBinContent(1);
  float NEntriesTTZToQQ = cutFlowTTZToQQ->GetBinContent(1);
  float NEntriesWJetsToLNu = cutFlowWJetsToLNu->GetBinContent(1);
  float NEntriesttHTobb = cutFlowttHTobb->GetBinContent(1);




  cout <<"===============================================" << endl;
  cout <<"              PreSelectionEvents" << endl;
  cout <<"==============================================="<< endl;
  cout << "Initial EntriesDYJetsToLL_M_10to50 " << NEntriesDYJetsToLL_M_10to50 << endl;
  cout << "Initial EntriesDYJetsToLL_M_50 " << NEntriesDYJetsToLL_M_50 << endl;
  cout << "Initial EntriesST_s_channel_4f " << NEntriesST_s_channel_4f << endl;
  cout << "Initial EntriesST_t_channel_top_4f " << NEntriesST_t_channel_top_4f << endl;
  cout << "Initial EntriesST_tW_top_5f " << NEntriesST_tW_top_5f << endl;
  cout << "Initial EntriesTTToHadronic " << NEntriesTTToHadronic << endl;
  cout << "Initial EntriesTTToSemiLept " << NEntriesTTToSemiLept << endl;
  cout << "Initial EntriesTTTo2L2Nu " << NEntriesTTTo2L2Nu << endl;
  cout << "Initial EntriesTTGJets " << NEntriesTTGJets << endl;
  cout << "Initial EntriesTTWJetsToLNu " << NEntriesTTWJetsToLNu << endl;
  cout << "Initial EntriesTTWJetsToQQ " << NEntriesTTWJetsToQQ << endl;
  cout << "Initial EntriesTTZToLLNuNu " << NEntriesTTZToLLNuNu << endl;
  cout << "Initial EntriesTTZToQQ " << NEntriesTTZToQQ << endl;
  cout << "Initial EntriesWJetsToLNu " << NEntriesWJetsToLNu << endl;
  cout << "Initial EntriesttHTobb " << NEntriesttHTobb << endl;

//March22 - Chaning to GetEntries from Integral
  float EntriesDYJetsToLL_M_10to50 = HtDYJetsToLL_M_10to50->GetEntries();
  float EntriesDYJetsToLL_M_50 = HtDYJetsToLL_M_50->GetEntries(); 
  float EntriesST_s_channel_4f = HtST_s_channel_4f->GetEntries();
  float EntriesST_t_channel_top_4f = HtST_t_channel_top_4f->GetEntries();
  float EntriesST_tW_top_5f = HtST_tW_top_5f->GetEntries();
  float EntriesTTToHadronic = HtTTToHadronic->GetEntries();
  float EntriesTTToSemiLept = HtTTToSemiLept->GetEntries();
  float EntriesTTTo2L2Nu = HtTTTo2L2Nu->GetEntries();
  float EntriesTTGJets = HtTTGJets->GetEntries();
  float EntriesTTWJetsToLNu = HtTTWJetsToLNu->GetEntries();
  float EntriesTTWJetsToQQ = HtTTWJetsToQQ->GetEntries();
  float EntriesTTZToLLNuNu = HtTTZToLLNuNu->GetEntries();
  float EntriesTTZToQQ = HtTTZToQQ->GetEntries();
  float EntriesWJetsToLNu = HtWJetsToLNu->GetEntries();
  float EntriesttHTobb = HtttHTobb->GetEntries();



  cout <<"===============================================" << endl;
  cout <<"              PostSelectionEvents" << endl;
  cout <<"==============================================="<< endl;
  cout << "Selected EntriesDYJetsToLL_M_10to50 " << EntriesDYJetsToLL_M_10to50 << endl;
  cout << "Selected EntriesDYJetsToLL_M_50 " << EntriesDYJetsToLL_M_50 << endl;
  cout << "Selected EntriesST_s_channel_4f " << EntriesST_s_channel_4f << endl;
  cout << "Selected EntriesST_t_channel_top_4f " << EntriesST_t_channel_top_4f << endl;
  cout << "Selected EntriesST_tW_top_5f " << EntriesST_tW_top_5f << endl;
  cout << "Selected EntriesTTToHadronic " << EntriesTTToHadronic << endl;
  cout << "Selected EntriesTTToSemiLept " << EntriesTTToSemiLept << endl;
  cout << "Selected EntriesTTTo2L2Nu " << EntriesTTTo2L2Nu << endl;
  cout << "Selected EntriesTTGJets " << EntriesTTGJets << endl;
  cout << "Selected EntriesTTWJetsToLNu " << EntriesTTWJetsToLNu << endl;
  cout << "Selected EntriesTTWJetsToQQ " << EntriesTTWJetsToQQ << endl;
  cout << "Selected EntriesTTZToLLNuNu " << EntriesTTZToLLNuNu << endl;
  cout << "Selected EntriesTTZToQQ " << EntriesTTZToQQ << endl;
  cout << "Selected EntriesWJetsToLNu " << EntriesWJetsToLNu << endl;
  cout << "Selected EntriesttHTobb " << EntriesttHTobb << endl;







  float BR_Hbb = 0.5824;
  float BR_Whad = 0.6741;
  float BR_Wlep = 1 - BR_Whad;

  float xsDYJetsToLL_M_10to50 = 18610;
  float xsDYJetsToLL_M_50 = 3*1921.8;
  float xsST_s_channel_4f = 11.36 * BR_Wlep;
  float xsST_t_channel_top_4f = 136.02;
  float xsST_tW_top_5f = 35.85;
  float xsTTToHadronic = 831.76 * BR_Whad * BR_Whad;
  float xsTTToSemiLept = 831.76 * 2 * BR_Whad * BR_Wlep;
  float xsTTTo2L2Nu = 831.76 * BR_Wlep * BR_Wlep; 
  float xsTTGJets = 0.5297;
  float xsTTWJetsToLNu = 0.2043;
  float xsTTWJetsToQQ = 0.4062;
  float xsTTZToLLNuNu = 0.2529;
  float xsTTZToQQ = 0.5297;
  float xsWJetsToLNu = 61526.7;
  float xsttHTobb = 0.5071 * BR_Hbb;



  float eDYJetsToLL_M_10to50;
  float eDYJetsToLL_M_50;
  float eST_s_channel_4f;
  float eST_t_channel_top_4f;
  float eST_tW_top_5f;
  float eTTToHadronic;
  float eTTToSemiLept;
  float eTTTo2L2Nu;
  float eTTGJets;
  float eTTWJetsToLNu;
  float eTTWJetsToQQ;
  float eTTZToLLNuNu;
  float eTTZToQQ;
  float eWJetsToLNu;
  float ettHTobb;

	eDYJetsToLL_M_10to50 = 1/NEntriesDYJetsToLL_M_10to50;
	eDYJetsToLL_M_50 = 1/NEntriesDYJetsToLL_M_50;
	eST_s_channel_4f = 1/NEntriesST_s_channel_4f;
	eST_t_channel_top_4f = 1/NEntriesST_t_channel_top_4f;
	eST_tW_top_5f = 1/NEntriesST_tW_top_5f;
	eTTToHadronic = 1/NEntriesTTToHadronic;
 cout << "eTTToHadronic = " << eTTToHadronic << endl;
	eTTToSemiLept = 1/NEntriesTTToSemiLept;
	eTTTo2L2Nu = 1/NEntriesTTTo2L2Nu;
	eTTGJets = 1/NEntriesTTGJets;
	eTTWJetsToLNu = 1/NEntriesTTWJetsToLNu;
	eTTWJetsToQQ = 1/NEntriesTTWJetsToQQ;
	eTTZToLLNuNu = 1/NEntriesTTZToLLNuNu;
	eTTZToQQ = 1/NEntriesTTZToQQ;
	eWJetsToLNu = 1/NEntriesWJetsToLNu;
	ettHTobb = 1/NEntriesttHTobb;
 


  float wDYJetsToLL_M_10to50 = xsDYJetsToLL_M_10to50 * eDYJetsToLL_M_10to50;
  float wDYJetsToLL_M_50 = xsDYJetsToLL_M_50 * eDYJetsToLL_M_50;
  float wST_s_channel_4f = xsST_s_channel_4f * eST_s_channel_4f;
  float wST_t_channel_top_4f = xsST_t_channel_top_4f * eST_t_channel_top_4f;
  float wST_tW_top_5f = xsST_tW_top_5f * eST_tW_top_5f;
  float wTTToHadronic = xsTTToHadronic*eTTToHadronic;
cout << "wTTToHadronic = " << wTTToHadronic << endl;
  float wTTToSemiLept = xsTTToSemiLept*eTTToSemiLept;
  float wTTTo2L2Nu = xsTTTo2L2Nu*eTTTo2L2Nu;
  float wTTGJets = xsTTGJets * eTTGJets;
  float wTTWJetsToLNu = xsTTWJetsToLNu * eTTWJetsToLNu;
  float wTTWJetsToQQ = xsTTWJetsToQQ * eTTWJetsToQQ;
  float wTTZToLLNuNu = xsTTZToLLNuNu * eTTZToLLNuNu;
  float wTTZToQQ = xsTTZToQQ * eTTZToQQ;
  float wWJetsToLNu = xsWJetsToLNu * eWJetsToLNu;
  float wttHTobb = xsttHTobb * ettHTobb;


  float IntLumi = 41529.0;


  cout <<"===============================================" << endl;
  cout <<"   CrossSection * 1/TotalPreSelectionEvents" << endl;
  cout <<"==============================================="<< endl;
  cout << "DYJetsToLL_M_10to50 " << wDYJetsToLL_M_10to50 << endl;
  cout << "DYJetsToLL_M_50 " << wDYJetsToLL_M_50 << endl;
  cout << "ST_s_channel_4f " << wST_s_channel_4f << endl;
  cout << "ST_t_channel_top_4f " << wST_t_channel_top_4f << endl;
  cout << "ST_tW_top_5f " << wST_tW_top_5f << endl;
  cout << "TTToHadronic " << wTTToHadronic << endl;
  cout << "TTToSemiLept " << wTTToSemiLept << endl;
  cout << "TTTo2L2Nu " << wTTTo2L2Nu << endl;
  cout << "TTGJets " << wTTGJets << endl;
  cout << "TTWJetsToLNu " << wTTWJetsToLNu << endl;
  cout << "TTWJetsToQQ " << wTTWJetsToQQ << endl;
  cout << "TTZToLLNuNu " << wTTZToLLNuNu << endl;
  cout << "TTZToQQ " << wTTZToQQ << endl;
  cout << "WJetsToLNu " << wWJetsToLNu << endl;
  cout << "ttHTobb " << wttHTobb << endl;
  cout <<"===============================================" << endl;
  cout <<"   PostSelectionEvents * IntLumi * CrossSection * 1/TotalPreSelectionEvents" << endl;
  cout <<"==============================================="<< endl;
  cout << "DYJetsToLL_M_10to50 " << EntriesDYJetsToLL_M_10to50 * IntLumi * wDYJetsToLL_M_10to50  << endl;
  cout << "DYJetsToLL_M_50 " << EntriesDYJetsToLL_M_50 * IntLumi * wDYJetsToLL_M_50 << endl;
  cout << "ST_s_channel_4f " << EntriesST_s_channel_4f * IntLumi * wST_s_channel_4f << endl;
  cout << "ST_t_channel_top_4f " << EntriesST_t_channel_top_4f * IntLumi * wST_t_channel_top_4f << endl;
  cout << "ST_tW_top_5f " << EntriesST_tW_top_5f * IntLumi * wST_tW_top_5f << endl;
  cout << "TTToHadronic " << EntriesTTToHadronic * IntLumi * wTTToHadronic << endl;
  cout << "TTToSemiLept " << EntriesTTToSemiLept * IntLumi * wTTToSemiLept << endl;
  cout << "TTTo2L2Nu " << EntriesTTTo2L2Nu * IntLumi * wTTTo2L2Nu << endl;
  cout << "TTGJets " << EntriesTTGJets * IntLumi * wTTGJets << endl;
  cout << "TTWJetsToLNu " << EntriesTTWJetsToLNu * IntLumi * wTTWJetsToLNu << endl;
  cout << "TTWJetsToQQ " << EntriesTTWJetsToQQ * IntLumi * wTTWJetsToQQ << endl;
  cout << "TTZToLLNuNu " << EntriesTTZToLLNuNu * IntLumi * wTTZToLLNuNu << endl;
  cout << "TTZToQQ " << EntriesTTZToQQ * IntLumi * wTTZToQQ << endl;
  cout << "WJetsToLNu " << EntriesWJetsToLNu * IntLumi * wWJetsToLNu << endl;
  cout << "ttHTobb " << EntriesttHTobb * IntLumi * wttHTobb << endl;







  ScaleToFit = 1.0;
  IntLumi = IntLumi * ScaleToFit; 
  

  float kFactor = 0.485351;// kFactor Normalize QCD to data from HT>1000 GeV for only specific cuts 
  kFactor = 1;
 
  cout<<"0000000000000000000000000000000000000000000000000000000000000000"<<endl;
  HtDYJetsToLL_M_10to50->Scale(IntLumi*wDYJetsToLL_M_10to50); 
  HtDYJetsToLL_M_50->Scale(IntLumi*wDYJetsToLL_M_50);
  
  HtST_s_channel_4f->Scale(IntLumi*wST_s_channel_4f);
  HtST_t_channel_top_4f->Scale(IntLumi*wST_t_channel_top_4f);
  HtST_tW_top_5f->Scale(IntLumi*wST_tW_top_5f);
  HtST->Add(HtST_s_channel_4f,HtST_t_channel_top_4f,1,1);
  HtST->Add(HtST_tW_top_5f);


  cout << "Events before scaling in TTToHadronic: "<< HtTTToHadronic->Integral()<< endl;
  HtTTToHadronic->Scale(IntLumi*wTTToHadronic);
  if(Var !="h_Jet0PT_SF" )HtTTToHadronic_ttlf->Scale(IntLumi*wTTToHadronic);
  cout << " TTToHadronic_ttb: IntLumi*wTTToHadronic = "<< IntLumi*wTTToHadronic << endl;
  cout << "HtTTToHadronic_ttb->Integral() = "<< HtTTToHadronic_ttb->Integral() << endl;
  if(HtTTToHadronic_ttb->Integral() >= 1){
  	HtTTToHadronic_ttb->Scale(IntLumi*wTTToHadronic);
  }
  if(Var != "h_Jet0PT_SF")HtTTToHadronic_ttbb->Scale(IntLumi*wTTToHadronic);
  if(Var != "h_Jet0PT_SF")HtTTToHadronic_tt2b->Scale(IntLumi*wTTToHadronic);
  if(Var != "h_Jet0PT_SF")HtTTToHadronic_ttcc->Scale(IntLumi*wTTToHadronic);
  cout << "Events after scaling in TTToHadronic: "<< HtTTToHadronic->Integral()<< endl;
  cout << "Events before scaling in TTToSemiLept: "<< HtTTToSemiLept->Integral()<< endl;
  HtTTToSemiLept->Scale(IntLumi*wTTToSemiLept);
  if(Var != "h_Jet0PT_SF")HtTTToSemiLept_ttlf->Scale(IntLumi*wTTToSemiLept);
  if(Channel != "_emu" && Var !="h_Jet0PT_SF")HtTTToSemiLept_ttb->Scale(IntLumi*wTTToSemiLept);
  if(Var != "h_Jet0PT_SF")HtTTToSemiLept_ttbb->Scale(IntLumi*wTTToSemiLept);
  if(Var != "h_Jet0PT_SF")HtTTToSemiLept_tt2b->Scale(IntLumi*wTTToSemiLept);
  if(Var != "h_Jet0PT_SF")HtTTToSemiLept_ttcc->Scale(IntLumi*wTTToSemiLept);
  cout << "Events after scaling in TTToSemiLept: "<< HtTTToSemiLept->Integral()<< endl;
  cout << "Events before scaling in TTTo2L2Nu: "<< HtTTTo2L2Nu->Integral()<< endl;
  HtTTTo2L2Nu->Scale(IntLumi*wTTTo2L2Nu);
  HtTTTo2L2Nu_ttlf->Scale(IntLumi*wTTTo2L2Nu);
  if(Channel != "_emu")HtTTTo2L2Nu_ttb->Scale(IntLumi*wTTTo2L2Nu);
  HtTTTo2L2Nu_ttbb->Scale(IntLumi*wTTTo2L2Nu);
  HtTTTo2L2Nu_tt2b->Scale(IntLumi*wTTTo2L2Nu);
  HtTTTo2L2Nu_ttcc->Scale(IntLumi*wTTTo2L2Nu);  

  Ht_ttlf->Add(HtTTToHadronic_ttlf,HtTTToSemiLept_ttlf,1,1);
  Ht_ttlf->Add(HtTTTo2L2Nu_ttlf);
  Ht_ttb->Add(HtTTToHadronic_ttb,HtTTToSemiLept_ttb,1,1);
  Ht_ttb->Add(HtTTTo2L2Nu_ttb);
  Ht_ttbb->Add(HtTTToHadronic_ttbb,HtTTToSemiLept_ttbb,1,1);
  Ht_ttbb->Add(HtTTTo2L2Nu_ttbb);
  Ht_tt2b->Add(HtTTToHadronic_tt2b,HtTTToSemiLept_tt2b,1,1);
  Ht_tt2b->Add(HtTTTo2L2Nu_tt2b);
  Ht_ttcc->Add(HtTTToHadronic_ttcc,HtTTToSemiLept_ttcc,1,1);
  Ht_ttcc->Add(HtTTTo2L2Nu_ttcc);



  cout << "Events after scaling in TTTo2L2Nu: "<< HtTTTo2L2Nu->Integral()<< endl;
  
  HtTTGJets->Scale(IntLumi*wTTGJets);
  HtTTWJetsToLNu->Scale(IntLumi*wTTWJetsToLNu);
  HtTTWJetsToQQ->Scale(IntLumi*wTTWJetsToQQ);
  HtTTZToLLNuNu->Scale(IntLumi*wTTZToLLNuNu);
  HtTTZToQQ->Scale(IntLumi*wTTZToQQ);
  HtTTV->Add(HtTTGJets,HtTTWJetsToLNu,1,1);
  HtTTV->Add(HtTTWJetsToQQ);
  HtTTV->Add(HtTTZToLLNuNu);
  HtTTV->Add(HtTTZToQQ);
  cout << "Events before scaling in WJetsToLNu: "<< HtWJetsToLNu->Integral()<< endl;
  
  HtWJetsToLNu->Scale(IntLumi*wWJetsToLNu);
  cout << "Events after scaling in WJetsToLNu: "<< HtWJetsToLNu->Integral()<< endl;
  
  cout << "Events before scaling in ST: "<< HtST->Integral()<< endl;

  cout << "Events before scaling in ttHtobb: "<< HtttHTobb->Integral()<< endl;
  HtVJets->Add(HtDYJetsToLL_M_10to50,HtDYJetsToLL_M_50,1,1);
  HtVJets->Add(HtWJetsToLNu);

  HtttHTobb->Scale(IntLumi*wttHTobb*20);
  cout << "Events after scaling(and x20) in ttHtobb: "<< HtttHTobb->Integral()<< endl;

  cout << "Total Events in DATA: "<< HtData->Integral() << endl;
  cout << "Total Events in MC(-ttH): "<< HtTTToHadronic->Integral() + HtTTToSemiLept->Integral() + HtTTTo2L2Nu->Integral() + HtVJets->Integral() + HtST->Integral() + HtTTV->Integral() << endl;




cout<<"555555555555555555555555555555555555555555555555555555"<<endl;
cout << "Integral of Distributions:" <<endl;
cout << "ttlf: " << Ht_ttlf->Integral()<<endl;
cout << "ttb: " << Ht_ttb->Integral()<<endl;
cout << "ttbb: " << Ht_ttbb->Integral()<<endl;
cout << "tt2b: " << Ht_tt2b->Integral()<<endl;
cout << "ttcc: " << Ht_ttcc->Integral()<<endl;
cout << "Single t: "<< HtST->Integral()<<endl;
cout << "V+Jets: " <<HtVJets->Integral()<<endl;
cout << "tt+V: " << HtTTV->Integral()<<endl;
cout << "DATA: "<<HtData->Integral()<<endl;
cout << "Entries of Distributions:" << endl;
cout << "ttlf: " << Ht_ttlf->GetEntries()<<endl;
cout << "ttb: " << Ht_ttb->GetEntries()<<endl;
cout << "ttbb: " << Ht_ttbb->GetEntries()<<endl;
cout << "tt2b: " << Ht_tt2b->GetEntries()<<endl;
cout << "ttcc: " << Ht_ttcc->GetEntries()<<endl;
cout << "Single t: "<< HtST->GetEntries()<<endl;
cout << "V+Jets: " <<HtVJets->GetEntries()<<endl;
cout << "tt+V: " << HtTTV->GetEntries()<<endl;
cout << "DATA: " << HtData->GetEntries()<<endl;
cout<<"555555555555555555555555555555555555555555555555555555"<<endl;


  cout<<"0000000000000000000000000000000000000000000000000000000000000000"<<endl;
  TH1F *Ht_TTV_VJets = new TH1F("Ht_TTV_VJets","TTV_VJets", NBins, BinMin, BinMax);//Binning needs to be matched with the histogram above
  Ht_TTV_VJets->Add(HtTTV, HtVJets, 1, 1);
  
  TH1F *Ht_TTV_VJets_ST = new TH1F("Ht_TTV_VJets_ST","TTV_VJets_ST",NBins,BinMin,BinMax);
  Ht_TTV_VJets_ST->Add(Ht_TTV_VJets,HtST, 1, 1);
  TH1F *Ht_TTV_VJets_ST_ttbb = new TH1F("Ht_TTV_VJets_ST_ttbb","TTV_VJets_ST_ttbb",NBins,BinMin,BinMax);
  Ht_TTV_VJets_ST_ttbb->Add(Ht_TTV_VJets_ST,Ht_ttbb, 1, 1);
  TH1F *Ht_TTV_VJets_ST_ttbb_tt2b = new TH1F("Ht_TTV_VJets_ST_ttbb_tt2b","TTV_VJets_ST_ttbb_tt2b",NBins,BinMin,BinMax);
  Ht_TTV_VJets_ST_ttbb_tt2b->Add(Ht_TTV_VJets_ST_ttbb,Ht_tt2b, 1, 1);
  TH1F *Ht_TTV_VJets_ST_ttbb_tt2b_ttb= new TH1F("Ht_TTV_VJets_ST_ttbb_tt2b_ttb","TTV_VJets_ST_ttbb_tt2b_ttb",NBins,BinMin,BinMax);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb->Add(Ht_TTV_VJets_ST_ttbb_tt2b,Ht_ttb, 1, 1);
  TH1F *Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc= new TH1F("Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc","TTV_VJets_ST_ttbb_tt2b_ttb_ttcc",NBins,BinMin,BinMax);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc->Add(Ht_TTV_VJets_ST_ttbb_tt2b_ttb,Ht_ttcc, 1, 1);
  TH1F *Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf= new TH1F("Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf","TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf",NBins,BinMin,BinMax);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->Add(Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc,Ht_ttlf, 1, 1);



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //ADD SCALING OF RATIOS FROM MC TO DATA!!!
  float TotalMCEntries = Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->Integral();//(BinMTop_LB, BinMTop_UB);
  float TotalDataEntries = HtData->Integral();//(BinMTop_LB, BinMTop_UB);
  if(withFF == 1) ScaleToFit = TotalDataEntries/TotalMCEntries;
  else ScaleToFit = 1;
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->Scale(ScaleToFit);
    Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc->Scale(ScaleToFit);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb->Scale(ScaleToFit);
   Ht_TTV_VJets_ST_ttbb_tt2b->Scale(ScaleToFit);
    Ht_TTV_VJets_ST_ttbb->Scale(ScaleToFit);
  Ht_TTV_VJets_ST->Scale(ScaleToFit);
  Ht_TTV_VJets->Scale(ScaleToFit);
    HtTTV->Scale(ScaleToFit);






  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
  cout << "Total Data Entries: " << TotalDataEntries << endl << "Total MC Entries: " << TotalMCEntries << endl;
  cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;


  TCanvas *cHT = new TCanvas("cHT", "cHT", 900, 1200);                                                                                                                                         
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0.0);//(0.01)
  pad1->Draw();                                                                                                                                                                                         
  pad1->cd();                                                                                                                                                                                           
  if(PlotLog == 1) pad1->SetLogy();                                                                                                                                                                                      
  pad1->SetTicks();      
  //pad1->SetMaximum(maxYValue*2); 

  float maxYValue = Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->GetMaximum();
  if (PlotLog == 1)maxYValue = maxYValue * 50;
  else maxYValue = maxYValue * 1.6;
  cout << "GetMaximum: " << maxYValue << endl;
  HtData->GetYaxis()->SetRange(0,maxYValue*2);
  

  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->SetLineColor(kBlack-2);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->SetFillColor(kRed-7);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->SetMinimum(0.5);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->SetMaximum(maxYValue);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->Draw("hist");
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->GetXaxis()->SetLabelSize(0);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->SetYTitle("Events");

  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc->SetLineColor(kBlack-2);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc->SetFillColor(kRed+1);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc->SetMinimum(0.5);//1E-1
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc->SetMaximum(maxYValue);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc->Draw("hist same");
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc->GetXaxis()->SetLabelSize(0);
  //  HtZJets_ttH_WJets_QCD_ttJets->SetYTitle("56 Entries/Event");
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc->SetYTitle("Events");
  
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb->SetLineColor(kBlack);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb->SetFillColor(kRed+2);
  Ht_TTV_VJets_ST_ttbb_tt2b_ttb->Draw("hist same");
  Ht_TTV_VJets_ST_ttbb_tt2b->SetLineColor(kBlack);
  Ht_TTV_VJets_ST_ttbb_tt2b->SetFillColor(kRed+3);
  Ht_TTV_VJets_ST_ttbb_tt2b->Draw("hist same");

  Ht_TTV_VJets_ST_ttbb->SetLineColor(kBlack);
  Ht_TTV_VJets_ST_ttbb->SetFillColor(kRed+3);
  Ht_TTV_VJets_ST_ttbb->Draw("hist same");

  Ht_TTV_VJets_ST->SetLineColor(kBlack);
  Ht_TTV_VJets_ST->SetFillColor(kViolet);
  Ht_TTV_VJets_ST->Draw("hist same");

  Ht_TTV_VJets->SetLineColor(kBlack);
  Ht_TTV_VJets->SetFillColor(kGreen-6);
  Ht_TTV_VJets->Draw("hist same");


  HtTTV->SetLineColor(kBlack);
  HtTTV->SetFillColor(kGray);
  HtTTV->Draw("hist same");


  HtttHTobb->SetLineColor(kBlue);
  HtttHTobb->GetYaxis()->SetRange(0,maxYValue);
  HtttHTobb->Draw("hist same");


  const Int_t n = NBins;
  cout << NBins << endl;
  double x[int(NBins)];
  double y[int(NBins)];
  double ex[int(NBins)];
  double ey[int(NBins)];
  double BinWidth = (BinMax - BinMin)/NBins;
  cout<< BinWidth << endl;
  for (int i=0; i<n; i++){
    x[i] = BinMin + i*BinWidth + BinWidth/2;
    y[i] = Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->GetBinContent(i+1);
    //cout << "y value for ge for iteration" << i << " is: " << y[i] << endl;
    ex[i] = BinWidth/2;
    ey[i] = Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->GetBinError(i+1);
    //    cout << "x:y:ex:ey = " << x[i] << " : " <<  y[i] << " : " << ex[i] << " : " << ey[i] << endl;  
  }
 //Draws errors on the stacked histogram 
  TGraphErrors* ge = new TGraphErrors(100, x, y, ex, ey);
  ge->SetFillColor(kBlue+1);
  ge->SetFillStyle(3001);
  //ge->Draw("same E2");

  HtData->SetMarkerStyle(20);
  HtData->SetMarkerSize(1.1);
  TLegend *leg = new TLegend(0.70,0.7,0.88,0.88,"", "brNDC");//0.70,0.6,0.9,0.88,"", "brNDC");
  leg->SetTextFont(62);
  leg->SetTextSize(0.020);//0.025
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillStyle(1001); 
  leg->SetFillColor(10); 
  leg->SetNColumns(2);
  leg->AddEntry(HtData,"#bf{data}","p");
  leg->AddEntry(HtttHTobb,"#bf{t#bar{t}H(Hbb)x20}","l");    //_TuneCP5(POWHEG+Pythia8)}","l");
  //leg->AddEntry(HtZJets_ttH_WJets_QCD_ttJets,"#bf{TTbar_P8M2T4(MG5+Pythia8)}","f");
  //leg->AddEntry(HtZJets_ttH_WJets_QCD_ttJets,"#bf{TTbar_P8M2T4(POWHEG+Pythia8)}","f");
  //if(TTBBFile == "TTBB_SherpaOL_hist_v3.root")leg->AddEntry(HtZJets_ttH_WJets_QCD_ttJets_ttbb,"#bf{t#bar{t}+b#bar{b} (Sherpa+OpenLoops)}","f");
  //else if(TTBBFile == "ttbb_4FS_ckm_amcatnlo_madspin_pythia8_hist_v4_ALL.root")leg->AddEntry(HtZJets_ttH_WJets_QCD_ttJets_ttbb,"#bf{t#bar{t}+b#bar{b} (amc@NLO_Madspin+Pythia8)}","f");
  //if (TTJetsFile == "TTbarIncl_hist_POWHEG_PY8.root")leg->AddEntry(HtZJets_ttH_WJets_QCD_ttJets,"#bf{TTbar_P8M2T4(POWHEG+Pythia8)}","f");
  //else if (TTJetsFile == "TTIncl_ttbbRemoved_POWHEG_PY8_hist_v4_ALLttbbRemoved.root")
  leg->AddEntry(Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf,"#bf{t#bar{t}+lf}","f");
  leg->AddEntry(Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc,"#bf{t#bar{t}+c#bar{c}}","f");
  leg->AddEntry(Ht_TTV_VJets_ST_ttbb_tt2b_ttb,"#bf{t#bar{t}+b}","f");
  leg->AddEntry(Ht_TTV_VJets_ST_ttbb_tt2b,"#bf{t#bar{t}+2b}","f");
  leg->AddEntry(Ht_TTV_VJets_ST_ttbb,"#bf{t#bar{t}+b#bar{b}}","f");
  leg->AddEntry(Ht_TTV_VJets_ST,"#bf{Single t}","f");
  leg->AddEntry(Ht_TTV_VJets,"#bf{V+Jets}","f");
  leg->AddEntry(HtTTV,"#bf{t#bar{t}+V}","f");
  leg->Draw();


  TLatex *   tex = new TLatex(0.15,0.95,"#bf{CMS Preliminary}");//was at 0.13,0.88......
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();

  TString t_Channel;
  if(Channel == "")t_Channel = "l^{+}l^{-}";
  if(Channel == "_ee")t_Channel = "e^{+}e^{-}";
  if(Channel == "_mumu")t_Channel = "#mu^{+}#mu^{-}";
  if(Channel == "_emu")t_Channel = "e^{+}#mu^{-}/e^{-}#mu^{+}";
  TLatex * text = new TLatex(0.15,0.86,t_Channel+", #geq 4 jets,  3 b-tag");//was at 0.13,0.88......
  text->SetNDC();
  text->SetTextAlign(13);
  text->SetTextFont(42);
  text->SetTextSize(0.035);
  text->SetLineWidth(2);
  text->Draw();
  //TLatex *   tex_var = new TLatex(0.35,0.96,Var);
  //tex_var->SetNDC();
  //tex_var->SetTextAlign(13);
  //tex_var->SetTextFont(42);
  //tex_var->SetTextSize(0.04);
  //tex_var->SetLineWidth(2);
  //tex_var->Draw();


  TLatex *   tex_lum = new TLatex(0.73,0.95,"41.5 fb^{-1} (13 TeV)");
  tex_lum->SetNDC();
  tex_lum->SetTextAlign(13);
  tex_lum->SetTextFont(42);
  tex_lum->SetTextSize(0.04);
  tex_lum->SetLineWidth(2);
  tex_lum->Draw();


  TH1F *ratio_ht = new TH1F("ratio_ht", "ratio_ht", NBins, BinMin, BinMax);
  for (int i=0; i<NBins; i++){
    //    if (HtZJets_ttH_WJets_QCD_ttJets->GetBinContent(i)!=0 && HtData->GetBinContent(i) != 0){
    if (Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->GetBinContent(i+1)!=0){
      ratio_ht->SetBinContent(i+1, (HtData->GetBinContent(i+1)/Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->GetBinContent(i+1)));
      ratio_ht->SetBinError(i+1, HtData->GetBinError(i+1)/Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->GetBinContent(i+1));
    }
  }                                                  
  
  
  TH1F *hDataGraph = new TH1F("DataGraph", "", NBins, BinMin, BinMax);
  hDataGraph->SetBarWidth(0);
  hDataGraph->SetMarkerSize(0.8);
  hDataGraph->SetLineWidth(2);
  hDataGraph->SetMarkerStyle(20);

  TH1F *hRatioGraph = new TH1F("DataRatio", "", NBins, BinMin, BinMax);
  for (int i=0; i<NBins; i++){
    float DataContent = 0;  
    DataContent = HtData->GetBinContent(i+1); 
    //    cout << " Bin "  << i+1 << " Data content == "  << DataContent << endl;
    hDataGraph->SetBinContent(i+1, DataContent);
    float RatioContent = 0;
    RatioContent = ratio_ht->GetBinContent(i+1);
    //    cout << " Bin "  << i+1 << " Ratio content == "  << RatioContent << endl;
    if (RatioContent != 0){ hRatioGraph->SetBinContent(i+1, RatioContent);}
    if (RatioContent == 0){ hRatioGraph->SetBinContent(i+1, -1);}//Set the bins without ratio content to -1

  }                                                                                                                                                                                                        
  //hRatioGraph->SetXTitle(Var);



  TGraphAsymmErrors * g = new TGraphAsymmErrors(hDataGraph);
  g->SetMarkerSize(1.5);
  g->SetMarkerStyle (20);
  TGraphAsymmErrors * gRatio = new TGraphAsymmErrors(hRatioGraph);
  gRatio->SetMarkerSize(1.0);
  gRatio->SetMarkerStyle (20);
  const double alpha = 1 - 0.6827;
  //  const double alpha = 1;
  for (int i = 0; i < g->GetN(); ++i) {
    int N = g->GetY()[i];
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
    g->SetPointEYlow(i, N-L);
    g->SetPointEYhigh(i, U-N);
    g->SetPointEXlow(i, 0);
    g->SetPointEXhigh(i, 0);
    double ht_ratio = 0;
    ht_ratio = ratio_ht->GetBinContent(i+1);
    if (ht_ratio <0) ht_ratio = (-1)*ht_ratio;
      gRatio->SetPointEYlow(i, (N-L)/Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->GetBinContent(i));
      gRatio->SetPointEYhigh(i, (U-N)/Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->GetBinContent(i));
      gRatio->SetPointEXlow(i, 0);
      gRatio->SetPointEXhigh(i, 0);
  }
    
  g->Draw("same P");
  g->SetLineColor(1);
  g->SetLineWidth(1);
  g->SetMarkerSize(0.9);
  gPad->RedrawAxis();

  cHT->cd(); 
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);//0.29
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.25);//
  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();
  pad2->SetTicks();



  //Using as a template to draw a graph on
  TH1F *ratio_ht_template = new TH1F("ratio_ht_template", "ratio_ht_template", NBins, BinMin, BinMax);
  ratio_ht_template->GetYaxis()->CenterTitle();
  ratio_ht_template->SetXTitle(Var);
  ratio_ht_template->SetYTitle("data/MC");


/*
  std::string CheckTopMass ("TopMass");
  std::string CheckWMass ("WMass");
  std::string CheckPt ("Pt");
  std::string CheckBT ("BT");
  std::string CheckMet ("Met");


  std::string CheckVar (Var);
  std::size_t foundTopMass = CheckVar.find(CheckTopMass);
  std::size_t foundWMass = CheckVar.find(CheckWMass);
  std::size_t foundBT = CheckVar.find(CheckBT);
  std::size_t foundMet = CheckVar.find(CheckMet);
  std::size_t foundPt = CheckVar.find(CheckPt);
*/

  if(foundNumberOf <= 100) s_XTitle = "Number Of";
  if(foundbTags <= 100)s_XTitle = s_XTitle+" b-Tags";

  if(foundLeading <= 100) s_XTitle = "Leading";
  else if(foundSubLeading <= 100)s_XTitle = "SubLeading";
  
  if(foundJet <= 100)s_XTitle = s_XTitle+" Jet";
  else if(foundLepton <= 100)s_XTitle = s_XTitle+" Lepton";
  
  if (foundPT <= 100)s_XTitle = s_XTitle+" p_{T} (GeV)";
  else if (foundEta <=100)s_XTitle = s_XTitle+" #eta";
  else if (foundPhi <=100)s_XTitle = s_XTitle+" #phi";
  else if (foundDeepCSV <= 100)s_XTitle = s_XTitle+" DeepCSV";


  if (foundMet <=100)ratio_ht_template->SetXTitle("E_{T}^{miss} (GeV)");
  else cout << "Something else went down and I don't know how to name the X-Axis....." << endl;

  if(foundNumberOf <=100 && foundJet <= 100)s_XTitle = s_XTitle+"s";


  if(found_ttlf <=100)s_XTitle = s_XTitle+" (t#bar{t}+lf)";
if(found_ttb <=100)s_XTitle = s_XTitle+" (t#bar{t}+b)";
if(found_ttbb <=100)s_XTitle = s_XTitle+" (t#bar{t}+bb)";
if(found_tt2b <=100)s_XTitle = s_XTitle+" (t#bar{t}+2b)";
if(found_ttcc <=100)s_XTitle = s_XTitle+" (t#bar{t}+cc)";


  ratio_ht_template->SetXTitle(s_XTitle);

//  if(Var =="h_TopMAssBestComb"||Var=="h_TopMassFirstComb"||Var =="h_TopMassFirstComb")ratio_ht_template->SetXTitle("M_{Top} (GeV)");
//  else if (Var == "HTAllJets")ratio_ht_template->SetXTitle("H_{T} (GeV)");
//  else if(Var == "TopPtAllComb" || Var == "h_Pt_PtOrdered_0"|| Var == "h_Pt_PtOrdered_1")ratio_ht_template->SetXTitle("p_{T, top} (GeV)");
//  else ratio_ht_template->SetXTitle("GeV");
  // ratio_ht_template->SetXTitle("Number of jets");
  // ratio_ht_template->SetXTitle("#eta_{top}");
  //  ratio_ht_template->SetXTitle("#phi_{top}");
  ratio_ht_template->SetMaximum(1.75);//2.5
  ratio_ht_template->SetMinimum(0.25);//-0.5
  ratio_ht_template->GetXaxis()->SetLabelSize(0.09);
  ratio_ht_template->GetYaxis()->SetLabelSize(0.09);
  ratio_ht_template->GetXaxis()->SetTitleSize(0.11);
  ratio_ht_template->GetYaxis()->SetTitleSize(0.09);
  ratio_ht_template->GetYaxis()->SetTitleOffset(0.5);
  ratio_ht_template->Draw("PE");

  gRatio->SetLineColor(1);
  gRatio->Draw("same P");


  TLine *line = new TLine(BinMin,1,BinMax,1);
  line->SetLineColor(kBlack);
  line->Draw("same");

  cHT->Modified();

  cout << "TTToHadronic == " << HtTTToHadronic->Integral() << " events" << endl;
  cout << "TTToSemiLept == " << HtTTToSemiLept->Integral() << " events" << endl;   //ERASED Integral(BinMTop_LB, BinMTop_UB) from all of these to get the full integral...
  cout << "TTTo2L2Nu == " << HtTTTo2L2Nu->Integral() << " events" << endl;
  cout << "WJetsToLNu == " << HtWJetsToLNu->Integral() << " events" << endl;
  cout << "ST == " << HtST->Integral() << " events" << endl;
  cout << "ttHtobb == " << HtttHTobb->Integral() << " events" << endl;
  cout << "============ purity MC Base ==========" << endl;
  cout << "Data == " << HtData->Integral() << " events" << endl;
  //cout << "Data Based Purity Mass WindoMass Window== " << 100*(HtTTTo2L2Nu->Integral(BinMTop_LB, BinMTop_UB)/HtData->Integral(BinMTop_LB, BinMTop_UB)) << " %" << endl;
  // cout << "Data Based Purity Full INTEGRAL== " << 100*(HtTTTo2L2Nu->Integral()/HtData->Integral()) << " %" << endl;

  // cout << "MC Based Purity == " << 100*(HtTTTo2L2Nu->Integral(BinMTop_LB, BinMTop_UB)/(	//HtQCD_All->Integral(BinMTop_LB, BinMTop_UB) 
  //								  HtTTToSemiLept->Integral(BinMTop_LB, BinMTop_UB) 
//								  + HtWJets->Integral(BinMTop_LB, BinMTop_UB)
//								  + HtZJets->Integral(BinMTop_LB, BinMTop_UB)
//								  + HtTTTo2L2Nu->Integral(BinMTop_LB,BinMTop_UB)
//								  //								  + HtTTBB->Integral(MTop_LB, MTop_UB)
//								  + HtttHtobb->Integral(BinMTop_LB, BinMTop_UB)
//								  ))<<  " %"<< endl;
  
  cout << "Testing MC Ratio of Purity to Data Purity" << Ht_TTV_VJets_ST_ttbb_tt2b_ttb_ttcc_ttlf->Integral()/HtData->Integral() << endl;;
  cout << "*************************" << endl;
  cout << "Variable:"<< Var << endl;
  cout << "kFactor: " << kFactor << endl;
  cout << "xsTTToHadronic: " << xsTTToHadronic << endl;
  cout << "xsTTToSemiLept: " << xsTTToSemiLept<< endl;
  cout << "xsTTTo2L2Nu: " << xsTTTo2L2Nu << endl;
  cout << "*************************"<< endl;
  //  gRatio->Draw("AP");
  //  TFile *ff = TFile::Open("QCD_Add_Histos.root","RECREATE"); 
  
  cHT->SaveAs("Histogram_"+Var+Channel+"_March8.pdf");

  TFile *ff;
  ff = TFile::Open("Plot_AddHistosWithWeight_"+Var+Channel+"_March8.root","RECREATE");


  cout << ScaleToFit << " is scale to fit." << endl;
  //HtQCD_All->Write();                                                                                                                                                                                   
  HtData->Write();
  ff->Close(); 
  
  HtTTToHadronic->Scale(ScaleToFit);
  HtTTToSemiLept->Scale(ScaleToFit);
  HtTTTo2L2Nu->Scale(ScaleToFit);
  //HtZJets->Scale(ScaleToFit);
  HtWJetsToLNu->Scale(ScaleToFit);
  HtST->Scale(ScaleToFit);
  //HtttHtobb->Scale(ScaleToFit);
  cout << "77777777777777 YIELDS 7777777777777777777" << endl;
  cout << "TTFH: " << HtTTToHadronic->Integral()<<endl;
  cout << "TTSL: " << HtTTToSemiLept->Integral()<<endl;
  cout << "TT2L2Nu: "<<HtTTTo2L2Nu->Integral()<<endl;
  cout << "VJets: " << HtVJets->Integral()<<endl;
  cout << "ST: " << HtST->Integral()<<endl;
  cout << "ttHtobb: "<<HtttHTobb->Integral() << endl;
  cout << "MC Total(-ttH): " << HtTTToHadronic->Integral() +  HtTTToSemiLept->Integral() + HtTTTo2L2Nu->Integral() + HtVJets->Integral() + HtST->Integral() << endl;
  cout << "Data Total: "<< HtData->Integral() << endl;
  cout << "ScaleToFit: "<< ScaleToFit << endl;
  cout << "7777777777777777777777777777777777777777" << endl;
    
  ScaleToFit = TotalDataEntries/TotalMCEntries;



  HtTTToHadronic->Scale(ScaleToFit);
  HtTTToSemiLept->Scale(ScaleToFit);
  HtTTTo2L2Nu->Scale(ScaleToFit);
  HtVJets->Scale(ScaleToFit);
  HtST->Scale(ScaleToFit);
  HtTTV->Scale(ScaleToFit);

  cout << "77777777777777 YIELDS IFF REQUIRED EQUAL EVENTS 7777777777777777777" << endl;
  cout << "TTFH: " << HtTTToHadronic->Integral()<<endl;
  cout << "TTSL: " << HtTTToSemiLept->Integral()<<endl;
  cout << "TT2L2Nu: "<<HtTTTo2L2Nu->Integral()<<endl;
  cout << "WJetsToLNu: " << HtWJetsToLNu->Integral()<<endl;
  cout << "ST: " << HtST->Integral()<<endl;
  cout << "ttHtobb: "<<HtttHTobb->Integral() << endl;
  cout << "MC Total(-ttH): " << HtTTToHadronic->Integral() +  HtTTToSemiLept->Integral() + HtTTTo2L2Nu->Integral() + HtVJets->Integral() + HtTTV->Integral() + HtST->Integral() << endl;
  cout << "Data Total: "<< HtData->Integral() << endl;
  cout << "ScaleToFit(FUDGE FACTOR): "<< ScaleToFit << endl;
  cout << "7777777777777777777777777777777777777777" << endl;

}


