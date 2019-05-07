//v9: Added cut to prevent multiple leptons from passing event selection and added electron iso cut.
//v10: Changed cuts for muons from tighter version as in v9 to looser ones listed on the sync exercise
#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLatex.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TString.h>
#include <TVector.h>
#include <THStack.h>
#include <TLorentzVector.h>
//This will compile in ROOT, still wont compile in g++
//	#include "/uscms_data/d3/3wolf3/ttH/Jan2019/Feb2019/CMSSW_9_4_9/src/TTH/CommonClassifier/interface/DLBDTVars.h"
//	#include "../ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilEventSelection.h"
//	#include "../ttH-LeptonPlusJets/AnalysisCode/interface/YggdrasilEventVars.h"
//	#include "DLBDTClassifier.h"

	//#include "ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilEventSelection.h"
	//#include "../ttH-LeptonPlusJets/YggdrasilTreeMaker/plugins/ttHYggdrasilEventSelection.cc"
//#include "/uscms_data/d3/3wolf3/ttH/Jan2019/Feb2019/CMSSW_9_4_9/src/TTH/CommonClassifier/interface/DLBDTVars.h"
//#include "/uscms_data/d3/3wolf3/ttH/Jan2019/Feb2019/CMSSW_9_4_9/src/TTH/CommonClassifier/interface/DLBDTClassifier.h "
#include "../ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilEventSelection.h"
#include "../ttH-LeptonPlusJets/AnalysisCode/interface/YggdrasilEventVars.h"
#include "../TTH/CommonClassifier/interface/DLBDTVars.h"
#include "../TTH/CommonClassifier/interface/DLBDTClassifier2.h"

//#include "ttHYggdrasilEventSelection.h"
//#include "YggdrasilEventVars.h"
//#include "DLBDTVars.h"
//#include "DLBDTClassifier.h"

#include <ctgmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
using namespace std;

void EventReadOut_MAIN_MarinoBDT_METFilter_PT15ge2j1b_May6(TString var_Sample = "",TString isMC = "", TString singleFile = ""){
//	gSystem->Load("YggdrasilEventVars.h");

		//int main(int argc, char **argv){
//	TString var_Sample = argv[1];
//	TString isMC = argv[2];
//	TString singleFile = argv[3];

		//TString var_Sample="", TString isMC="",TString singleFile=""){
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
	if(singleFile == ""){
		//var_Sample = "SingleElectron_PeriodB";
		//singleFile = "root://cmseos.fnal.gov//store/group/lpctthrun2/UVA/ICHEP2018/MET/DoubleMuon/crab_ICHEP18_postMCsync_v0_DoubleMuon_PeriodB/181113_163516/0000/yggdrasil_treeMaker_ttH_sync_11-12-18_v24_data_1-17.root"; 
		var_Sample = "MuonEG_PeriodB"; 
			//"TTTo2L2Nu";
			//"ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8"; 
			//"MuonEG_PeriodB";
			//"DoubleMuon_PeriodB";
		singleFile = "yggdrasil_treeMaker_ttH_sync_02-25-19_v24_MuonEG_PeriodB-singleSync_data_useJSON.root";
		//"yggdrasil_treeMaker_ttH_sync_MuonEG_PeriodB_10-29-18_v20_data_2017B.root";
		//"yggdrasil_treeMaker_ttH_sync_10-29-18_v20_data_2017B.root";
		//"yggdrasil_treeMaker_ttH_sync_02-22-19_v24_MuonEG_PeriodB-singleSync_data_addMETFilters.root";
		//"/eos/uscms/store/group/lpctthrun2/UVA/ICHEP2018/MC/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v0_ttHTobb/181107_184846/0000/yggdrasil_treeMaker_ttH_sync_11-06-18_v26_recipeTest_44.root"; 
			//"/eos/uscms/store/group/lpctthrun2/UVA/ICHEP2018/MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v0_TTto2L2Nu/181107_184524/0000/yggdrasil_treeMaker_ttH_sync_11-06-18_v26_recipeTest_44.root";
				//"/eos/uscms/store/group/lpctthrun2/UVA/ICHEP2018/MC/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v0_ttHTobb/181107_184846/0000/yggdrasil_treeMaker_ttH_sync_11-06-18_v26_recipeTest_44.root"; 
			//"yggdrasil_treeMaker_ttH_sync_02-22-19_v24_MuonEG_PeriodB-singleSync_data_addMETFilters.root";
			//"Output_Histograms_EventReadOut_MAIN_Feb2_PT15_ge2j1bAllDL_AllSF_ptpt_CONDORnc_11-06-18_v26_recipeTest_4.root"; 

			//"yggdrasil_treeMaker_ttH_sync_02-20-19_v24_MuonEG_PeriodB-singleSync_data_addMETFilters.root";//File with 8/9 MET Filters attached...
				//"yggdrasil_treeMaker_ttH_sync_02-20-19_v24_periodB-emu-singleSync_data_addMETFilters.root";
			//"yggdrasil_treeMaker_ttH_sync_02-06-19_v25_MuonEG_PeriodB-singleSync_data_wOut4JCut.root";//As of Feb21st
			
			//yggdrasil_treeMaker_ttH_sync_10-29-18_v20_data_DATASYNC.root";
			//"/eos/uscms/store/group/lpctthrun2/UVA/ICHEP2018/MuonEG/MuonEG/crab_ICHEP18_postMCsync_v0_MuonEG_PeriodB/181115_162211/0000/yggdrasil_treeMaker_ttH_sync_11-12-18_v24_data_1.root";
			//"/uscms_data/d3/3wolf3/ttH/Jan2019/yggdrasil_treeMaker_ttH_sync_10-29-18_v20_data_DATASYNC.root";
			//"root://cmseos.fnal.gov//store/group/lpctthrun2/UVA/ICHEP2018/MET/DoubleMuon/crab_ICHEP18_postMCsync_v0_DoubleMuon_PeriodB/181113_163516/0000/yggdrasil_treeMaker_ttH_sync_11-12-18_v24_data_9.root";
		isMC = "false";
	}	
	std::string weightpath="./bdt_dl";	//Previous weights used//"./dlbdtweights_v5/";
	DLBDTClassifier bdt(weightpath);

//if(singleFile == "")singleFile = "root://cmseos.fnal.gov//store/group/lpctthrun2/UVA/ICHEP2018/MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v0_TTto2L2Nu/181107_184524/0000/yggdrasil_treeMaker_ttH_sync_11-06-18_v26_recipeTest_18.root"; isMC = "true";
	if(singleFile == ""){
		cout << "Missing third argument for singleFile!!!!" << endl;
		return;
	}
	//"root://cmseos.fnal.gov//store/group/lpctthrun2/UVA/ICHEP2018/MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v0_TTto2L2Nu/181107_184524/0000/yggdrasil_treeMaker_ttH_sync_11-06-18_v26_recipeTest_1-1.root";
	std::string CheckDoubleMuon ("DoubleMuon");
	std::string CheckDoubleEG ("DoubleEG");
	std::string CheckMuonEG ("MuonEG");
	std::string CheckSingleElectron ("SingleElectron");
	std::string CheckSingleMuon ("SingleMuon");
	std::string CheckTTTo ("TTTo");
	std::string ChecksingleFile (singleFile);
	std::string Checkvar_Sample (var_Sample);	
	std::string CheckPeriodB ("PeriodB");
	std::string CheckPeriodC ("PeriodC");
	std::string CheckPeriodD ("PeriodD");
	std::string CheckPeriodE ("PeriodE");
	std::string CheckPeriodF ("PeriodF");

	//TFile* fOut = new TFile("Output_Histograms_EventReadOut_MAIN_BDT_METFilter_PT15ge2j1b_Mar1"+var_Sample+".root", "RECREATE");


	std::size_t foundSingleElectron = ChecksingleFile.find(CheckSingleElectron);
	std::size_t foundSingleMuon = ChecksingleFile.find(CheckSingleMuon);
	std::size_t foundDoubleMuon = ChecksingleFile.find(CheckDoubleMuon);
	std::size_t foundDoubleEG = ChecksingleFile.find(CheckDoubleEG);
	std::size_t foundMuonEG = ChecksingleFile.find(CheckMuonEG);
	std::size_t foundTTTo = ChecksingleFile.find(CheckTTTo);
	std::size_t foundPeriodB = Checkvar_Sample.find(CheckPeriodB);
	std::size_t foundPeriodC = Checkvar_Sample.find(CheckPeriodC);
	std::size_t foundPeriodD = Checkvar_Sample.find(CheckPeriodD);
	std::size_t foundPeriodE = Checkvar_Sample.find(CheckPeriodE);
	std::size_t foundPeriodF = Checkvar_Sample.find(CheckPeriodF);


	int SingleElectronData = 0;
	if(foundSingleElectron <= 200)SingleElectronData = 1;
	int SingleMuonData = 0;
	if(foundSingleMuon <= 200)SingleMuonData = 1;
	int DoubleMuonData = 0;
	if(foundDoubleMuon <= 200)DoubleMuonData = 1;
	int DoubleEGData = 0;
	if(foundDoubleEG <= 200)DoubleEGData = 1;
	int MuonEGData = 0;
	if(foundMuonEG <= 200)MuonEGData = 1;
	int TTbarMC = 0;
	if(foundTTTo <= 200)TTbarMC = 1;

	cout << "DoubleMuonData = " << DoubleMuonData << endl;
	cout << "foundDoubleMuon = " << foundDoubleMuon << endl;
	cout << "FoundTTTo = " << foundTTTo << endl;
	cout << "FoundMuongEGData = " << foundMuonEG << endl;
	cout << "MuonEGData = " << MuonEGData << endl;

	int is_MC = 0;
	if(isMC=="true")is_MC = 1;
	TString var_DL = "AllDL";//"ee";//Only used for naming purposes, and some cout statements for debugging.
	TString var_TrigSFxx = "ptpt";//Used to determine which TrigSF's are applied to events, should be "ptpt"
	if(var_DL==""||var_TrigSFxx==""){
		cout<<"Arguments must be: var_DL={ee,mumu,emu} , var_TrigSFxx=={ptpt, etaeta, 0pt0eta, 1pt1eta}"<<endl;
		return;
	}
	//TString var_DL = "AllDL";
	//TString 
	
//var_DL = "ee";

	//TString var_DL = "mumu";
	//TString var_DL = "emu";
//============================================
	//TString 
	//TString var_SF = "noElMuCSVSF";
	TString var_SF = "AllSF";
	//TString var_SF = "eleSF";
	//TString var_SF = "csvSF";
	//TString var_SF = "muonSF";
//=================================================
//   WHICH TRIGGER SCALE FACTORS SHOULD BE APPLIED:
	//TString var_TrigSFxx = "noTrigSF";
	//TString 
	//var_TrigSFxx = "ptpt";
	//TString var_TrigSFxx = "etaeta";
	//TString var_TrigSFxx = "0pt0eta";
	//TString var_TrigSFxx = "1pt1eta";
//================================================
/*
	TString var_TrigSFee = "DLee_el0pt_el1pt";
	TString var_TrigSFee = "DLee_el0eta_el1eta";
	TString var_TrigSFee = "DLee_el0pt_el1eta";
	TString var_TrigSFee = "DLee_el1pt_el0eta";	

	TString var_TrigSFemu = "DLemu_mu0pt_el0pt";
        TString var_TrigSFemu = "DLemu_mu0eta_el0eta";
        TString var_TrigSFemu = "DLemu_mu0pt_mu0eta";
        TString var_TrigSFemu = "DLemu_el0pt_el0eta"; 

	TString var_TrigSFmumu = "DLmumu_mu0pt_mu1pt";
        TString var_TrigSFmumu = "DLmumu_mu0eta_mu1eta";
        TString var_TrigSFmumu = "DLmumu_mu0pt_mu1eta";
        TString var_TrigSFmumu = "DLmumu_mu1pt_mu0eta";
*/


	//if(var_Sample=="TTTo2L2Nu")

	int EventNumber = 617147596;//
//346108907;//0.00113188
		//368472746;//bdt = 0.315
		//570813386;//bdt 0.28
		//621323193;//0.299
		//604337191;//bdt 0.28
		//271861991;//bdt: 0.18
		//391723;//bdt: 0.37
		//13012599;//bdtscore = 0.55
		//610860840;//Event that SHOULD pass, using to test why the CMSSW compiled version of this script doesnt save it to file....
		//113017724;
		//301904563;
		//280686615;
		//540885024;//DataSync event that has same charge leptons
		//397824342;
		//540885024;//Missing event for DataSync that UVA doesn't have
		
		//5337543;//ee event 2896987;//7641754; //2896987;
	//is_ee:  2896987 caused by mll being 76 and RWTH has 75.9995;
	//is_emu:  3280379 caused by lepton not being tight
	//	6895521
	//	7872745 NO TRIGGERS
	//	7925484 SL ONLY TRIGGER HT150....
	//	6894933
	int RunSingleEvent = 0; //Set to 0 for run over full sample, Set to 1 to run over single event specified by EventNumber
	int DEBUG = 0;// Set to 1 if you want ALL THE MESSAGES
	long n_entries = -1;//1000000000;// Sets the number of events to run over
	//int is_MC = 1;
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
	int ttHSync = 0;
	int EleB = 0;
	int EleC = 0;
	int EleD = 0;
	int EleE = 0;
	int EleF = 0;
	int MuonB = 0;
	int MuonC = 0;
	int MuonD = 0;
	int MuonE = 0;
	int MuonF = 0;
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888	
	int ttHTobb = 0;
	int TTTo2L2Nu = 0;
	int TTToSemiLept = 0;
	int TTWJetsToLNu = 0;//CHECK THE FILE THAT THIS RUNS OVER!!!!!!!
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
	int ST_t_chan_top = 0;
	int ST_t_chan_antitop = 0;
	int ST_tW_top = 0;
	int ST_tW_antitop = 0;
	int ST_s_chan = 0;
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
	int WJetsToLNu = 0;
	int WW = 0;
	int WZ = 0;
	int ZZ = 0;
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
	int DYJetsToLL_M10 = 0;
	int DYJetsToLL_M50 = 0;
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

	if(var_Sample=="TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8")TTTo2L2Nu = 1;
	if(var_Sample=="ttHSync")ttHSync = 1;
	if(var_Sample=="ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8")ttHTobb= 1;
	 if(var_Sample=="TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8")TTToSemiLept = 1;
	 if(var_Sample=="ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8")ST_s_chan = 1;
	 if(var_Sample=="ST_t_chan_antitop")ST_t_chan_antitop = 1;
	 if(var_Sample=="ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8")ST_tW_antitop = 1;
	 if(var_Sample=="ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8")ST_tW_top = 1;
	if(var_Sample=="DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8")DYJetsToLL_M10 = 1;
	if(var_Sample=="DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8")DYJetsToLL_M50 = 1;



//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
	int applyIndCSVsf = 0;//Turn to 1 if you want individual CSV correction factors
	int applyJERsf = 0;
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
////88888888888888888888888888888888888888888888888888888888888888888888888888888888888888

	double TrigSF = 1;
	int xbin_TrigSF = -99;
	int ybin_TrigSF = -99;
	double xValue_TrigSF = -99;
	double yValue_TrigSF = -99;

	double TrigSF_0Eta = -99;
	double TrigSF_0PT = -99;
	double TrigSF_1Eta = -99;
	double TrigSF_1PT = -99; 

	int cutNJets = 2;// was 4 in version 9
	int cutNumbOfGoodJets = 2;

	double cutMinJetPT = 20.0;
	double cutJetPT = 30.0;//35;
	double cutMinJetEta = 2.4;
	double cutJetEta = 2.4;
	double cutDeepCSV = 0.4941;//0.4941 for v9

	double cutMuonRelIso = 0.25;//0.15 for v9

	double cutMinElePT = 15.0;//15.0;
	double cutElePT = 25.0;//Was 38 for v9
	double cutMinEleEta = 2.4;
	double cutEleEta = 2.4;//was 2.1 in v9	
	double cutEleIso = 0.25;

	double cutMinMuonPT = 15;
	double cutMuonPT = 25;//30GeV in v9

	double cutSL_MuonPT = 26;
	double cutSL_ElePT = 30;

	double cutDL_LeadElePT = 0;

	double cutMinMuonIso = 0.25;
	double cutMuonIso = 0.25;
	double cutVetoMuonIso = 0.25;
	double cutVetoMuonPT = 15;
	double cutVetoMuonEta = 2.4;
	//double cutMuonPTLeading = 25.0;
	//double cutMuonPTSubLeading = 15.0;
		//https://gitlab.cern.ch/ttH/reference/blob/ICHEP18/definitions/ICHEP18.md
	double cutMinMuonEta = 2.4;
	double cutMuonEta = 2.4;//2.1 in v9
	
	int cutNbTags = 1;//Was 2 in version 9

//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
	int isData = 0;
	//int isMC = 0;
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
//	TString outSingleFile;
//        outSingleFile = singleFile;
//        int strLength = outSingleFile.Length();
//        outSingleFile.Remove(0,strLength-33);
//        cout << outSingleFile << endl;

	//TH2D *h_EleSF= new TH1F("h_Jet3PT", "Jet3PT",60, 0, 600);




	//	TString inputSLTriggerSF = "SingleEG_JetHT_Trigger_Scale_Factors_ttHbb_Data_MC_v2.0.root";
	//"/uscms_data/d3/3wolf3/ttH/Testing94X/Yggdrasil2018Oct/CMSSW_9_4_0_patch1/src/ttH-LeptonPlusJets/AnalysisCode/macros/SingleEG_JetHT_Trigger_Scale_Factors_ttHbb_Data_MC_v2.0.root";
       // TFile *f_SLTriggerSF = new TFile(inputSLTriggerSF);
       // 	TH2D *h_SLTriggerSF =(TH2D*)f_SLTriggerSF->Get("SFs_ele_pt_ele_sceta_ele28_ht150_OR_ele35_2017BCDEF");


	TString inputDLTriggerSF = "inputData/tth_dileptonic_2DscaleFactors_withSysts_2017BCDEF_11-26-18.root";
	//"/uscms_data/d3/3wolf3/ttH/Testing94X/Yggdrasil2018Oct/CMSSW_9_4_0_patch1/src/ttH-LeptonPlusJets/AnalysisCode/macros/tth_dileptonic_2DscaleFactors_withSysts_2017BCDEF_11-26-18.root";
//tth_dileptonic_2DscaleFactors_10-05-18.root";
        TFile *f_DLTriggerSF = new TFile(inputDLTriggerSF);
		TH2D *h_DLTriggerSFee;
		TH2D *h_DLTriggerSFemu;
		TH2D *h_DLTriggerSFmumu;
		TH2D *h_DLTriggerSFll;	

		TH2D *h_DLTriggerSF_ee_el0pt_el0eta =(TH2D*)f_DLTriggerSF->Get("h_DoubleEl_OR__X__allMET_el0_pt_vs_eta_withSysts");
        	TH2D *h_DLTriggerSF_ee_el0pt_el1pt =(TH2D*)f_DLTriggerSF->Get("h_DoubleEl_OR__X__allMET_el0_pt_vs_el1_pt_withSysts");
		TH2D *h_DLTriggerSF_ee_el0eta_el1eta =(TH2D*)f_DLTriggerSF->Get("h_DoubleEl_OR__X__allMET_el0_eta_vs_el1_eta_withSysts");
                TH2D *h_DLTriggerSF_ee_el1pt_el1eta =(TH2D*)f_DLTriggerSF->Get("h_DoubleEl_OR__X__allMET_el1_pt_vs_eta_withSysts");

		TH2D *h_DLTriggerSF_emu_mu0pt_el0pt =(TH2D*)f_DLTriggerSF->Get("h_EMu_OR__X__allMET_mu0_pt_vs_el0_pt_withSysts");
		TH2D *h_DLTriggerSF_emu_mu0eta_el0eta =(TH2D*)f_DLTriggerSF->Get("h_EMu_OR__X__allMET_mu0_eta_vs_el0_eta_withSysts");
		TH2D *h_DLTriggerSF_emu_mu0pt_mu0eta =(TH2D*)f_DLTriggerSF->Get("h_EMu_OR__X__allMET_mu0_pt_vs_eta_withSysts");
		TH2D *h_DLTriggerSF_emu_el0pt_el0eta =(TH2D*)f_DLTriggerSF->Get("h_EMu_OR__X__allMET_el0_pt_vs_eta_withSysts");
	
		TH2D *h_DLTriggerSF_mumu_mu0pt_mu0eta =(TH2D*)f_DLTriggerSF->Get("h_DoubleMu_OR__X__allMET_mu0_pt_vs_eta_withSysts");
                TH2D *h_DLTriggerSF_mumu_mu1pt_mu1eta =(TH2D*)f_DLTriggerSF->Get("h_DoubleMu_OR__X__allMET_mu1_pt_vs_eta_withSysts");
                TH2D *h_DLTriggerSF_mumu_mu0pt_mu1pt =(TH2D*)f_DLTriggerSF->Get("h_DoubleMu_OR__X__allMET_mu0_pt_vs_mu1_pt_withSysts");
		TH2D *h_DLTriggerSF_mumu_mu0eta_mu1eta =(TH2D*)f_DLTriggerSF->Get("h_DoubleMu_OR__X__allMET_mu0_eta_vs_mu1_eta_withSysts");

	if(var_TrigSFxx=="ptpt"){
		h_DLTriggerSFee = (TH2D*)f_DLTriggerSF->Get("h_DoubleEl_OR__X__allMET_el0_pt_vs_el1_pt_withSysts");
		h_DLTriggerSFemu = (TH2D*)f_DLTriggerSF->Get("h_EMu_OR__X__allMET_mu0_pt_vs_el0_pt_withSysts");
		h_DLTriggerSFmumu = (TH2D*)f_DLTriggerSF->Get("h_DoubleMu_OR__X__allMET_mu0_pt_vs_mu1_pt_withSysts");
	}else if(var_TrigSFxx=="etaeta"){
		h_DLTriggerSFee = (TH2D*)h_DLTriggerSF_ee_el0eta_el1eta->Clone();
                h_DLTriggerSFemu = (TH2D*)h_DLTriggerSF_emu_mu0eta_el0eta->Clone();
                h_DLTriggerSFmumu = (TH2D*)h_DLTriggerSF_mumu_mu0eta_mu1eta->Clone();
	}else if(var_TrigSFxx=="0pt0eta"){
                h_DLTriggerSFee = (TH2D*)h_DLTriggerSF_ee_el0pt_el0eta->Clone();
                h_DLTriggerSFemu = (TH2D*)h_DLTriggerSF_emu_el0pt_el0eta->Clone();
                h_DLTriggerSFmumu = (TH2D*)h_DLTriggerSF_mumu_mu0pt_mu0eta->Clone();
        }else if(var_TrigSFxx=="1pt1eta"){
                h_DLTriggerSFee = (TH2D*)h_DLTriggerSF_ee_el1pt_el1eta->Clone();
                h_DLTriggerSFemu = (TH2D*)h_DLTriggerSF_emu_mu0pt_mu0eta->Clone();
                h_DLTriggerSFmumu = (TH2D*)h_DLTriggerSF_mumu_mu1pt_mu1eta->Clone();
        }



	TString inputEleIsoSF = "inputData/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root";
//"/uscms_data/d3/3wolf3/ttH/Testing94X/Yggdrasil2018Oct/CMSSW_9_4_0_patch1/src/ttH-LeptonPlusJets/AnalysisCode/macros/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root";
	TFile *f_EleIsoSF = new TFile(inputEleIsoSF);
	TH2D *h_EleIsoSF =(TH2D*)f_EleIsoSF->Get("EGamma_SF2D");
	double EleIsoSF_Error = -99;

	TString inputEleIDSF = "inputData/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root";
	//"/uscms_data/d3/3wolf3/ttH/Testing94X/Yggdrasil2018Oct/CMSSW_9_4_0_patch1/src/ttH-LeptonPlusJets/AnalysisCode/macros/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root";
        TFile *f_EleIDSF = new TFile(inputEleIDSF);
        TH2D *h_EleIDSF =(TH2D*)f_EleIDSF->Get("EGamma_SF2D");
        double EleIDSF_Error = -99;

	TString inputMuonIsoSF = "inputData/muon_Run2017-Nov17_RunBCDEF_SF_ISO.root";
	//"/uscms_data/d3/3wolf3/ttH/Testing94X/Yggdrasil2018Oct/CMSSW_9_4_0_patch1/src/ttH-LeptonPlusJets/AnalysisCode/macros/muon_Run2017-Nov17_RunBCDEF_SF_ISO.root";
        TFile *f_MuonIsoSF = new TFile(inputMuonIsoSF);
        TH2D *h_MuonIsoSF =(TH2D*)f_MuonIsoSF->Get("NUM_LooseRelIso_DEN_TightIDandIPCut_pt_abseta");
        double MuonIsoSF_Error = -99;

	TString inputMuonIDSF = "inputData/muon_Run2017-Nov17_RunBCDEF_SF_ID.root";
	//"/uscms_data/d3/3wolf3/ttH/Testing94X/Yggdrasil2018Oct/CMSSW_9_4_0_patch1/src/ttH-LeptonPlusJets/AnalysisCode/macros/muon_Run2017-Nov17_RunBCDEF_SF_ID.root";
        TFile *f_MuonIDSF = new TFile(inputMuonIDSF);
        TH2D *h_MuonIDSF =(TH2D*)f_MuonIDSF->Get("NUM_TightID_DEN_genTracks_pt_abseta");
        double MuonIDSF_Error = -99;



/*
	
	cout << "TESTING EleSF Root File:"<<endl;
	cout << h_EleSF->GetNbinsX() << " should be the number of X bins" << endl;
	cout << h_EleSF->GetNbinsY() << " should be the number of Y bins" << endl;
	cout << h_EleSF->GetBinContent(1,1) << " should be the Bin Content for (1,1) " <<endl;	
	cout << h_EleSF->GetBinContent(2,2) << " should be the Bin Content for (2,2) " <<endl;
	cout << h_EleSF->GetBinContent(3,3) << " should be the Bin Content for (3,3) " <<endl;
	cout << h_EleSF->GetMaximum() << " should be the global maximum " << endl;
	h_EleSF->GetXaxis()->SetRange();
	cout<< h_EleSF->GetMaximum() << " should be the maximum of the Xaxis" << endl;
	h_EleSF->GetYaxis();
	cout<< h_EleSF->GetMaximum() << " should be the maximum of the Yaxis" << endl;
	cout << h_EleSF->FindFirstBinAbove(2.1,2) << "should be the Eta Bin above 2.1 on the Xaxis"<<endl;
	cout << h_EleSF->FindFirstBinAbove(39,1) << "should be the ET Bin above 39 on the Yaxis" <<endl;


	cout << h_EleSF->GetXaxis()->GetBinLowEdge(3) << " is the low edge of bin 3 on Xaxis" << endl;
	int xbin_EleSF = -1;
	int ybin_EleSF = -1;
	double eleSF_scEta = -1.5;
	double eleSF_Et = 80.987;
	for (int xbin = 0; xbin <= h_EleSF->GetNbinsX();xbin++){
		if(eleSF_scEta > h_EleSF->GetXaxis()->GetBinLowEdge(xbin))xbin_EleSF = xbin;
		else continue;
		for (int ybin = 0; ybin <= h_EleSF->GetNbinsY();ybin++){
			if(eleSF_Et > h_EleSF->GetYaxis()->GetBinLowEdge(ybin))ybin_EleSF = ybin;
	                else continue;
		}
	}

*/	

 	TChain *chain = new TChain("ttHTreeMaker/worldTree");

	TString fDir = "/eos/uscms/store/group/lpctthrun2/UVA/ICHEP2018/";
	TString output;
	if  (ttHSync ==1){
	fDir = fDir+"MC/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v0_ttHTobb/181107_184846/0000/";
		//uscms_data/d3/3wolf3/ttH/Testing94X/Yggdrasil2018Oct/CMSSW_9_4_0_patch1/src/ttH-LeptonPlusJets/AnalysisCode/macros/";
		//Old sample //"/uscms_data/d3/3wolf3/ttH/Testing94X/Yggdrasil2018March/CMSSW_9_4_0_patch1/src/ttH-LeptonPlusJets/AnalysisCode/macros/";
	output = "ttHSync";
	isData = 0;//This was set to 1???????????????????
	}
	else if (EleB == 1){
		fDir =fDir+"SingleElectron/DataElB_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_212028/0000/";
		//"SingleElectron/DataElB_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_155726/0000/";
		output = "EleB";
		isData = 1;
	}else if(EleC == 1){
        	fDir =fDir+"SingleElectron/DataElC_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_212254/0000/";//0001
		//SingleElectron/DataElC_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_155820/0000/";
		//fDir="/eos/uscms/store/group/lpctthrun2/SingleElectron/DataElC_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_155820/0001/"
        	output = "EleC";
		isData = 1;
	}else if (EleD == 1){
		fDir = fDir+"SingleElectron/DataElD_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_212343/0000/";
		//SingleElectron/DataElD_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_155911/0000/";
		output = "EleD";
		isData = 1;
	}else if (EleE == 1){
	        fDir = fDir+"SingleElectron/DataElE_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_212430/0000/";//0001
		//SingleElectron/DataElE_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_155957/0000/";
		//SingleElectron/DataElE_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_155957/0001/
	        output = "EleE";
		isData = 1;
	}else if (EleF == 1){
	        fDir = fDir+"SingleElectron/DataElF_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_212520/0000/";//0001
		//SingleElectron/DataElF_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160051/0000/";
		//SingleElectron/DataElF_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160051/0001/
	        output = "EleF";
		isData = 1;
	}else if (MuonB == 1){
	        fDir = fDir+"SingleMuon/DataMuB_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_212716/0000/";//0001
		//SingleMuon/DataMuB_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160150/0000/";
		//SingleMuon/DataMuB_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160150/0001/
	        output = "MuonB";
		isData = 1;
	}else if (MuonC == 1){
	        fDir = fDir+"SingleMuon/DataMuC_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_212803/0000/";
		//SingleMuon/DataMuC_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160242/0000/";
		//SingleMuon/DataMuC_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160242/0001	
	        output = "MuonC";
		isData = 1;
	}else if (MuonD == 1){
	        fDir = fDir+"SingleMuon/DataMuD_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_213006/0000/";
		//SingleMuon/DataMuD_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160334/0000/";
		output = "MuonD";
		isData = 1;
	}else if (MuonE == 1){
	        fDir = fDir+"SingleMuon/DataMuE_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_213057/0000/";//0001
		//SingleMuon/DataMuE_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160423/0000/";
	        output = "MuonE";
		isData = 1;
	}else if (MuonF == 1){
	        fDir = fDir+"SingleMuon/DataMuF_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_213202/0000/";//0001 0002
		//SingleMuon/DataMuF_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160510/0000/";
		//SingleMuon/DataMuF_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160510/0001/
		//SingleMuon/DataMuF_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_160510/0002/
	        output = "MuonF";
		isData = 1;
	}else if (ttHTobb ==1){
		fDir = fDir+"MC/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v0_ttHTobb/181107_184846/0000/";
			//ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/crab_ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8_May11_v2/180511_180918/0000/";
			//"UVA/ICHEP2018/MC/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/crab_ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8_SecondTry/180409_132153/0000/";
		output = "ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8";
	}else if (TTTo2L2Nu == 1){
		fDir = fDir+"MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v0_TTto2L2Nu/181107_184524/0000/";
		//"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/ttto2l2nu_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_214955/0000/";
		//TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/ttto2l2nu_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_161808/0000/";
		output = "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8";
	}else if (TTToSemiLept == 1){
	       	fDir = fDir+"MC/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v0_TTtoSemiLep/181107_184710/0000/";
		//"TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/tttosemilepPSweight_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_215914/0000/";//0001 0002
		//TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/tttosemilep_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_161947/0000/";
	        output = "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8";
	}else if (TTWJetsToLNu == 1){
	        fDir = fDir+"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/WjetIncl_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag_resubmit6/180703_223550/0000/";
		//WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/WjetIncl_Satoshi_jobsubmit_2018_03_13_SubjetsV4/180314_141906/0000/";
		//"UVA/ICHEP2018/MC/TTW/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/180411_204038/0000/";
	        output = "TTWJetsToLNu";
	}else if (ST_t_chan_top == 1){
	        fDir = fDir+"MC/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/crab_ICHEP18_postMCsync_v2_ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/181110_011510/0000/";
		//ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/tchan_top_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_214238/0000/";
		//ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/tchan_top_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_161502/0000/";
	        output = "ST_t-channel_top";
	}else if (ST_t_chan_antitop == 1){
	        fDir = fDir+"ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/tchan_tbar_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_214152/0000/";
		//ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/tchan_tbar_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_161413/0000/";
	        output = "ST_t-channel_antitop";
	}else if (ST_tW_top == 1){
	        fDir = fDir+"UVA/ICHEP2018/MC/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v2_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/181110_011842/0000/";
		//ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/tW_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_214011/0000/";
		//ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/tW_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_161242/0000/";
	        output = "ST_tW_top";
	}else if (ST_tW_antitop == 1){
	        fDir = fDir+"MC/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v2_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/181110_011643/0000/";
		//ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/tbarW_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_214101/0000/";
		//ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/tbarW_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_161326/0000/";
	        output = "ST_tW_antitop";
	}else if (ST_s_chan == 1){
	        fDir = fDir+"MC/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/crab_ICHEP18_postMCsync_v2_ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/181110_011003/0000/";
		//ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/schan_both_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_213920/0000/";
		//ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/schan_both_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_161154/0000/";
	        output = "ST_s-channel";
	}else if (WJetsToLNu == 1){
		fDir = fDir+"UVA/ICHEP2018/MC/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/crab_ICHEP18_postMCsync_v0_WJetsToLNu/181107_185158/0000/";
		//WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/WjetIncl_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag_resubmit6/180703_223550/0000/";
		//WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/WjetIncl_Satoshi_jobsubmit_2018_03_13_SubjetsV4/180314_141906/0000/";
		output = "WJetsToLNu";
	}else if (WW == 1){
	        fDir = fDir+"WW_TuneCP5_13TeV-pythia8/ww_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag_resubmit5/180702_233750/0000/";
		//WW_TuneCP5_13TeV-pythia8/ww_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_162033/0000/";
	        output = "WW";
	}else if (WZ == 1){
	        fDir = fDir+"WZ_TuneCP5_13TeV-pythia8/wz_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag_resubmit5/180702_233837/0000/";
		//WZ_TuneCP5_13TeV-pythia8/wz_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_162129/0000/";
	        output = "WZ";
	}else if (ZZ == 1){
	        fDir = fDir+"ZZ_TuneCP5_13TeV-pythia8/zz_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag_resubmit5/180702_233928/0000/";
		//ZZ_TuneCP5_13TeV-pythia8/zz_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_162226/0000/";
	        output = "ZZ";
	}else if (DYJetsToLL_M10==1){
	        fDir = fDir+"MC/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_ICHEP18_postMCsync_v2_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/181110_002009/0000/";
		//DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/ZjetLowMass_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_213829/0000/";
		//DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/ZjetLowMass_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_161110/0000/";
	        output = "DYJetsToLL_M-10to50";
	}else if (DYJetsToLL_M50==1){
		fDir = fDir+"MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_ICHEP18_postMCsync_v2_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/181110_010811/0000/";
		//DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/ZjetIncl_Satoshi_jobsubmit_2018_06_28_FatJetDoubleBtag/180629_213736/0000/";//0001 0002
		//DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/ZjetIncl_Satoshi_jobsubmit_2018_03_13_SubjetsV2/180313_161015/0000/";
		output = "DYJetsToLL_M-50";
	}


	output = var_Sample;


	


	if(DEBUG == 1)cout << "Output of the files will be: " << output << endl;
	//if(ttHSync == 1) chain->Add(fDir+"yggdrasil_treeMaker_ttHTobb_SYNC_v6_NoJECLocal_NoPUPPI.root");//"yggdrasil_treeMaker_ttHTobb_SYNC_v3.root");
	//if(ttHSync == 1) chain->Add(fDir+"yggdrasil_treeMaker_ttH_sync_10-29-18_v20.root");

			//"yggdrasil_treeMaker_ttH_sync_10-18-18_v17.root");

			//"yggdrasil_treeMaker_ttH_sync_EleIDv1_JEC_MET_ESm.root");

			//"yggdrasil_treeMaker_ttH_sync_10-12-18_v14_full.root");
			//"yggdrasil_treeMaker_ttH_sync_Oct14_BensRecipe.root");


			//"yggdrasil_treeMaker_ttH_sync_EleIDv1_JEC_MET_ESm.root");
			//yggdrasil_treeMaker_ttH_syncEleIDv2_JEC_MET_v14.root
			//"yggdrasil_treeMaker_ttH_sync_10-08-18_v7_elID_full.root");
			//"yggdrasil_treeMaker_ttH_sync_TestingMET_v14.root");

			//"yggdrasil_treeMaker_ttH_sync_10-10-18_v9_full.root");

			//"yggdrasil_treeMaker_ttH_sync_TRIAL_EleIDv1_JEC_Oct10_v14.root"); // Evan's MET recipe
			//"yggdrasil_treeMaker_ttH_sync_10-10-18_v9_full.root ");//Ben's version with correct JetPT
		
		//"yggdrasil_treeMaker_ttH_sync_Oct10_ElectronID_JER_v14.root");
		//"yggdrasil_treeMaker_ttH_sync_10-08-18_v7_elID_full.root");

				//"yggdrasil_treeMaker_ttH_sync_Oct10_ElectronID_JER_v14.root"); Newest Version!
					//"yggdrasil_treeMaker_ttH_sync_10-08-18_v7_elID_full.root");//Oct 9 version without JER
						//"yggdrasil_treeMaker_ttH_sync_10-08-18_v7_elID_full.root"); //Ben's Sample
					//"yggdrasil_treeMaker_ttH_sync_Oct8_ElectronID_v14.root");

	//chain->Add(fDir+singleFile);
	chain->Add(singleFile);
			//"yggdrasil_treeMaker_ttH_sync_11-06-18_v26_recipeTest_1.root");
			//"yggdrasil_treeMaker_1.root");
	stringstream ss;
	string str;
	TString chainFile;
/*	if(ttHSync == 0){
		for (int nFiles=2;nFiles<=999;nFiles++){
			ss << nFiles;
			str = ss.str();
			if(DEBUG == 1)	cout << "This is the string: "<< str << endl;
			//if(ttHTobb==1)chainFile = fDir+"yggdrasil_treeMaker_"+str+".root";
			chainFile = fDir+"yggdrasil_treeMaker_ttH_sync_11-06-18_v26_recipeTest_"+str+".root";
			chain->Add(chainFile);
			if(DEBUG == 1) cout << chainFile << endl;	
			ss.str("");
		}
	}
*/
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
  // *** 1. Define histograms and canvasses
    	//TCanvas *cJetAllPT = new TCanvas("cJetAllPT", "cJetAllPT", 50, 50, 800, 600);
  	
	double sfEvent = -99;
	double sfTrig = -99;
	double sfEleIso = -99;
	double sfEleID = -99;
	double sfMuonIso = -99;
	double sfMuonDO = -99;
	double sfCSV = -99;



	//MCWeights Mar3
	TH1F *h_mcWeight_value = new TH1F("h_mcWeight_value","mcWeight_value",100,-10,10);
	TH1F *h_wgt_generator = new TH1F("h_wgt_generator_","wgt_generator_",200,0,2);
	//USELESS as value is always 25474122  TH1F *h_wgt_nGen_ = new TH1F("h_wgt_nGen_","wgt_nGen_",100,-10,10);
	


	//BDT 2/21
	TH1F *h_bdt_score_ge2jge1b = new TH1F("h_bdt_score_ge2jge1b", "bdt_score_ge2jge1b",50, -1.0, 1.0);
	TH1F *h_bdt_score_3j2b = new TH1F("h_bdt_score_3j2b", "bdt_score_3j2b",50, -1.0, 1.0);
	TH1F *h_bdt_score_3j3b = new TH1F("h_bdt_score_3j3b", "bdt_score_3j3b",10, -1.0, 1.0);
	TH1F *h_bdt_score_ge4j2b = new TH1F("h_bdt_score_ge4j2b", "bdt_score_ge4j2b",40, -1.0, 1.0);
	TH1F *h_bdt_score_ge4j3b = new TH1F("h_bdt_score_ge4j3b", "bdt_score_ge4j3b",20, -1.0, 1.0);
	TH1F *h_bdt_score_ge4jge4b = new TH1F("h_bdt_score_ge4jge4b", "bdt_score_ge4jge4b",10, -1.0, 1.0);


	TH1F *h_bdt_score_ge2jge1b_SF = new TH1F("h_bdt_score_ge2jge1b_SF", "bdt_score_ge2jge1b_SF",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j2b_SF = new TH1F("h_bdt_score_3j2b_SF", "bdt_score_3j2b_SF",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j3b_SF = new TH1F("h_bdt_score_3j3b_SF", "bdt_score_3j3b_SF",10, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j2b_SF = new TH1F("h_bdt_score_ge4j2b_SF", "bdt_score_ge4j2b_SF",40, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j3b_SF = new TH1F("h_bdt_score_ge4j3b_SF", "bdt_score_ge4j3b_SF",20, -1.0, 1.0);
        TH1F *h_bdt_score_ge4jge4b_SF = new TH1F("h_bdt_score_ge4jge4b_SF", "bdt_score_ge4jge4b_SF",10, -1.0, 1.0);

	TH1F *h_bdt_score_ge2jge1b_SF_ttlf = new TH1F("h_bdt_score_ge2jge1b_SF_ttlf", "bdt_score_ge2jge1b_SF",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j2b_SF_ttlf = new TH1F("h_bdt_score_3j2b_SF_ttlf", "bdt_score_3j2b_SF_ttlf",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j3b_SF_ttlf = new TH1F("h_bdt_score_3j3b_SF_ttlf", "bdt_score_3j3b_SF_ttlf",10, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j2b_SF_ttlf = new TH1F("h_bdt_score_ge4j2b_SF_ttlf", "bdt_score_ge4j2b_SF_ttlf",40, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j3b_SF_ttlf = new TH1F("h_bdt_score_ge4j3b_SF_ttlf", "bdt_score_ge4j3b_SF_ttlf",20, -1.0, 1.0);
        TH1F *h_bdt_score_ge4jge4b_SF_ttlf = new TH1F("h_bdt_score_ge4jge4b_SF_ttlf", "bdt_score_ge4jge4b_SF_ttlf",10, -1.0, 1.0);

	TH1F *h_bdt_score_ge2jge1b_SF_ttb = new TH1F("h_bdt_score_ge2jge1b_SF_ttb", "bdt_score_ge2jge1b_SF_ttb",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j2b_SF_ttb = new TH1F("h_bdt_score_3j2b_SF_ttb", "bdt_score_3j2b_SF_ttb",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j3b_SF_ttb = new TH1F("h_bdt_score_3j3b_SF_ttb", "bdt_score_3j3b_SF_ttb",10, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j2b_SF_ttb = new TH1F("h_bdt_score_ge4j2b_SF_ttb", "bdt_score_ge4j2b_SF_ttb",40, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j3b_SF_ttb = new TH1F("h_bdt_score_ge4j3b_SF_ttb", "bdt_score_ge4j3b_SF_ttb",20, -1.0, 1.0);
        TH1F *h_bdt_score_ge4jge4b_SF_ttb = new TH1F("h_bdt_score_ge4jge4b_SF_ttb", "bdt_score_ge4jge4b_SF_ttb",10, -1.0, 1.0);

	TH1F *h_bdt_score_ge2jge1b_SF_ttbb = new TH1F("h_bdt_score_ge2jge1b_SF_ttbb", "bdt_score_ge2jge1b_SF_ttbb",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j2b_SF_ttbb = new TH1F("h_bdt_score_3j2b_SF_ttbb", "bdt_score_3j2b_SF_ttbb",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j3b_SF_ttbb = new TH1F("h_bdt_score_3j3b_SF_ttbb", "bdt_score_3j3b_SF_ttbb",10, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j2b_SF_ttbb = new TH1F("h_bdt_score_ge4j2b_SF_ttbb", "bdt_score_ge4j2b_SF_ttbb",40, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j3b_SF_ttbb = new TH1F("h_bdt_score_ge4j3b_SF_ttbb", "bdt_score_ge4j3b_SF_ttbb",20, -1.0, 1.0);
        TH1F *h_bdt_score_ge4jge4b_SF_ttbb = new TH1F("h_bdt_score_ge4jge4b_SF_ttbb", "bdt_score_ge4jge4b_SF_ttbb",10, -1.0, 1.0);

	TH1F *h_bdt_score_ge2jge1b_SF_tt2b = new TH1F("h_bdt_score_ge2jge1b_SF_tt2b", "bdt_score_ge2jge1b_SF_tt2b",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j2b_SF_tt2b = new TH1F("h_bdt_score_3j2b_SF_tt2b", "bdt_score_3j2b_SF_tt2b",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j3b_SF_tt2b = new TH1F("h_bdt_score_3j3b_SF_tt2b", "bdt_score_3j3b_SF_tt2b",10, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j2b_SF_tt2b = new TH1F("h_bdt_score_ge4j2b_SF_tt2b", "bdt_score_ge4j2b_SF_tt2b",40, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j3b_SF_tt2b = new TH1F("h_bdt_score_ge4j3b_SF_tt2b", "bdt_score_ge4j3b_SF_tt2b",20, -1.0, 1.0);
        TH1F *h_bdt_score_ge4jge4b_SF_tt2b = new TH1F("h_bdt_score_ge4jge4b_SF_tt2b", "bdt_score_ge4jge4b_SF_tt2b",10, -1.0, 1.0);

	TH1F *h_bdt_score_ge2jge1b_SF_ttcc = new TH1F("h_bdt_score_ge2jge1b_SF_ttcc", "bdt_score_ge2jge1b_SF_ttcc",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j2b_SF_ttcc = new TH1F("h_bdt_score_3j2b_SF_ttcc", "bdt_score_3j2b_SF_ttcc",50, -1.0, 1.0);
        TH1F *h_bdt_score_3j3b_SF_ttcc = new TH1F("h_bdt_score_3j3b_SF_ttcc", "bdt_score_3j3b_SF_ttcc",10, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j2b_SF_ttcc = new TH1F("h_bdt_score_ge4j2b_SF_ttcc", "bdt_score_ge4j2b_SF_ttcc",40, -1.0, 1.0);
        TH1F *h_bdt_score_ge4j3b_SF_ttcc = new TH1F("h_bdt_score_ge4j3b_SF_ttcc", "bdt_score_ge4j3b_SF_ttcc",20, -1.0, 1.0);
        TH1F *h_bdt_score_ge4jge4b_SF_ttcc = new TH1F("h_bdt_score_ge4jge4b_SF_ttcc", "bdt_score_ge4jge4b_SF_ttcc",10, -1.0, 1.0);



	TH1F *h_sfEvent = new TH1F("h_sfEvent", "sfEvent",200, 0.5, 1.5);
	TH1F *h_sfTrig = new TH1F("h_sfTrig", "sfTrig",200, 0.5, 1.5);
	TH1F *h_sfEleIso = new TH1F("h_sfEleIso", "sfEleIso",200, 0, 2);
	TH1F *h_sfEleID = new TH1F("h_sfEleID", "sfEleID",200, 0, 2);
	TH1F *h_sfMuonIso = new TH1F("h_sfMuonIso", "sfMuonIso",200, 0, 2);
	TH1F *h_sfMuonID = new TH1F("h_sfMuonID", "sfMuonID",200, 0, 2);
	TH1F *h_sfCSV = new TH1F("h_sfCSV", "sfCSV",200, 0, 2);



	TH1F *h_SumJetPT = new TH1F("h_SumJetPT", "SumJetPT",200, 0, 2000);
  	TH1F *h_SumJetEta = new TH1F("h_SumJetEta", "SumJetEta",50, -5, 5);  
	TH1F *h_MET = new TH1F("h_MET", "MET",60, 0, 600);
	TH1F *h_MET_ee = new TH1F("h_MET_ee", "MET_ee",60, 0, 600);
	TH1F *h_MET_emu = new TH1F("h_MET_emu", "MET_emu",60, 0, 600);
	TH1F *h_MET_mumu = new TH1F("h_MET_mumu", "MET_mumu",60, 0, 600);
	TH1F *h_MET_ttlf = new TH1F("h_MET_ttlf", "MET",60, 0, 600);
	TH1F *h_MET_ttlf_ee = new TH1F("h_MET_ttlf_ee", "MET_ee",60, 0, 600);
	TH1F *h_MET_ttlf_emu = new TH1F("h_MET_ttlf_emu", "MET_emu",60, 0, 600);
	TH1F *h_MET_ttlf_mumu = new TH1F("h_MET_ttlf_mumu", "MET_mumu",60, 0, 600);
	TH1F *h_MET_ttb = new TH1F("h_MET_ttb", "MET",60, 0, 600);
	TH1F *h_MET_ttb_ee = new TH1F("h_MET_ttb_ee", "MET_ee",60, 0, 600);
	TH1F *h_MET_ttb_emu = new TH1F("h_MET_ttb_emu", "MET_emu",60, 0, 600);
	TH1F *h_MET_ttb_mumu = new TH1F("h_MET_ttb_mumu", "MET_mumu",60, 0, 600);
	TH1F *h_MET_ttbb = new TH1F("h_MET_ttbb", "MET",60, 0, 600);
	TH1F *h_MET_ttbb_ee = new TH1F("h_MET_ttbb_ee", "MET_ee",60, 0, 600);
	TH1F *h_MET_ttbb_emu = new TH1F("h_MET_ttbb_emu", "MET_emu",60, 0, 600);
	TH1F *h_MET_ttbb_mumu = new TH1F("h_MET_ttbb_mumu", "MET_mumu",60, 0, 600);
	TH1F *h_MET_tt2b = new TH1F("h_MET_tt2b", "MET",60, 0, 600);
	TH1F *h_MET_tt2b_ee = new TH1F("h_MET_tt2b_ee", "MET_ee",60, 0, 600);
	TH1F *h_MET_tt2b_emu = new TH1F("h_MET_tt2b_emu", "MET_emu",60, 0, 600);
	TH1F *h_MET_tt2b_mumu = new TH1F("h_MET_tt2b_mumu", "MET_mumu",60, 0, 600);
	TH1F *h_MET_ttcc = new TH1F("h_MET_ttcc", "MET",60, 0, 600);
	TH1F *h_MET_ttcc_ee = new TH1F("h_MET_ttcc_ee", "MET_ee",60, 0, 600);
	TH1F *h_MET_ttcc_emu = new TH1F("h_MET_ttcc_emu", "MET_emu",60, 0, 600);
	TH1F *h_MET_ttcc_mumu = new TH1F("h_MET_ttcc_mumu", "MET_mumu",60, 0, 600);

        TH1F *h_MET_SF = new TH1F("h_MET_SF", "MET_SF",60, 0, 600);
        TH1F *h_MET_SF_ee = new TH1F("h_MET_SF_ee", "MET_SF_ee",60, 0, 600);
        TH1F *h_MET_SF_emu = new TH1F("h_MET_SF_emu", "MET_SF_emu",60, 0, 600);
        TH1F *h_MET_SF_mumu = new TH1F("h_MET_SF_mumu", "MET_SF_mumu",60, 0, 600);
        TH1F *h_MET_SF_ttlf = new TH1F("h_MET_SF_ttlf", "MET_SF",60, 0, 600);
        TH1F *h_MET_SF_ttlf_ee = new TH1F("h_MET_SF_ttlf_ee", "MET_SF_ee",60, 0, 600);
        TH1F *h_MET_SF_ttlf_emu = new TH1F("h_MET_SF_ttlf_emu", "MET_SF_emu",60, 0, 600);
        TH1F *h_MET_SF_ttlf_mumu = new TH1F("h_MET_SF_ttlf_mumu", "MET_SF_mumu",60, 0, 600);
        TH1F *h_MET_SF_ttb = new TH1F("h_MET_SF_ttb", "MET_SF",60, 0, 600);
        TH1F *h_MET_SF_ttb_ee = new TH1F("h_MET_SF_ttb_ee", "MET_SF_ee",60, 0, 600);
        TH1F *h_MET_SF_ttb_emu = new TH1F("h_MET_SF_ttb_emu", "MET_SF_emu",60, 0, 600);
        TH1F *h_MET_SF_ttb_mumu = new TH1F("h_MET_SF_ttb_mumu", "MET_SF_mumu",60, 0, 600);
        TH1F *h_MET_SF_ttbb = new TH1F("h_MET_SF_ttbb", "MET_SF",60, 0, 600);
        TH1F *h_MET_SF_ttbb_ee = new TH1F("h_MET_SF_ttbb_ee", "MET_SF_ee",60, 0, 600);
        TH1F *h_MET_SF_ttbb_emu = new TH1F("h_MET_SF_ttbb_emu", "MET_SF_emu",60, 0, 600);
        TH1F *h_MET_SF_ttbb_mumu = new TH1F("h_MET_SF_ttbb_mumu", "MET_SF_mumu",60, 0, 600);
        TH1F *h_MET_SF_tt2b = new TH1F("h_MET_SF_tt2b", "MET_SF",60, 0, 600);
        TH1F *h_MET_SF_tt2b_ee = new TH1F("h_MET_SF_tt2b_ee", "MET_SF_ee",60, 0, 600);
        TH1F *h_MET_SF_tt2b_emu = new TH1F("h_MET_SF_tt2b_emu", "MET_SF_emu",60, 0, 600);
        TH1F *h_MET_SF_tt2b_mumu = new TH1F("h_MET_SF_tt2b_mumu", "MET_SF_mumu",60, 0, 600);
        TH1F *h_MET_SF_ttcc = new TH1F("h_MET_SF_ttcc", "MET_SF",60, 0, 600);
        TH1F *h_MET_SF_ttcc_ee = new TH1F("h_MET_SF_ttcc_ee", "MET_SF_ee",60, 0, 600);
        TH1F *h_MET_SF_ttcc_emu = new TH1F("h_MET_SF_ttcc_emu", "MET_SF_emu",60, 0, 600);
        TH1F *h_MET_SF_ttcc_mumu = new TH1F("h_MET_SF_ttcc_mumu", "MET_SF_mumu",60, 0, 600);


//0000000000000000000000000000000000000000000


	TH1F *h_NumberOfPV = new TH1F("h_NumberOfPV", "NumberOfPV",80, 0, 80);
  
	TH1F *h_JetAllPT = new TH1F("h_JetAllPT", "JetAllPT",60, 0, 600);
  	TH1F *h_Jet0PT= new TH1F("h_Jet0PT", "Jet0PT",60, 0, 600);
	TH1F *h_Jet0PT_ee= new TH1F("h_Jet0PT_ee", "Jet0PT_ee",60, 0, 600);	
	TH1F *h_Jet0PT_emu= new TH1F("h_Jet0PT_emu", "Jet0PT_emu",60, 0, 600);
	TH1F *h_Jet0PT_mumu= new TH1F("h_Jet0PT_mumu", "Jet0PT_mumu",60, 0, 600);
	TH1F *h_Jet0PT_ttlf= new TH1F("h_Jet0PT_ttlf", "Jet0PT",60, 0, 600);
	TH1F *h_Jet0PT_ttlf_ee= new TH1F("h_Jet0PT_ttlf_ee", "Jet0PT",60, 0, 600);
	TH1F *h_Jet0PT_ttlf_emu= new TH1F("h_Jet0PT_ttlf_emu", "Jet0PT",60, 0, 600);
	TH1F *h_Jet0PT_ttlf_mumu= new TH1F("h_Jet0PT_ttlf_mumu", "Jet0PT",60, 0, 600);
	TH1F *h_Jet0PT_ttb= new TH1F("h_Jet0PT_ttb", "Jet0PT",60, 0, 600);
	TH1F *h_Jet0PT_ttb_ee= new TH1F("h_Jet0PT_ttb_ee", "Jet0PT_ee",60, 0, 600);
	TH1F *h_Jet0PT_ttb_emu= new TH1F("h_Jet0PT_ttb_emu", "Jet0PT_emu",60, 0, 600);
	TH1F *h_Jet0PT_ttb_mumu= new TH1F("h_Jet0PT_ttb_mumu", "Jet0PT_mumu",60, 0, 600);
	TH1F *h_Jet0PT_ttbb= new TH1F("h_Jet0PT_ttbb", "Jet0PT",60, 0, 600);
	TH1F *h_Jet0PT_ttbb_ee= new TH1F("h_Jet0PT_ttbb_ee", "Jet0PT_ee",60, 0, 600);
	TH1F *h_Jet0PT_ttbb_emu= new TH1F("h_Jet0PT_ttbb_emu", "Jet0PT_emu",60, 0, 600);
	TH1F *h_Jet0PT_ttbb_mumu= new TH1F("h_Jet0PT_ttbb_mumu", "Jet0PT_mumu",60, 0, 600);
	TH1F *h_Jet0PT_tt2b= new TH1F("h_Jet0PT_tt2b", "Jet0PT",60, 0, 600);
	TH1F *h_Jet0PT_tt2b_ee= new TH1F("h_Jet0PT_tt2b_ee", "Jet0PT_ee",60, 0, 600);
	TH1F *h_Jet0PT_tt2b_emu= new TH1F("h_Jet0PT_tt2b_emu", "Jet0PT_emu",60, 0, 600);
	TH1F *h_Jet0PT_tt2b_mumu= new TH1F("h_Jet0PT_tt2b_mumu", "Jet0PT_mumu",60, 0, 600);
	TH1F *h_Jet0PT_ttcc= new TH1F("h_Jet0PT_ttcc", "Jet0PT",60, 0, 600);
	TH1F *h_Jet0PT_ttcc_ee= new TH1F("h_Jet0PT_ttcc_ee", "Jet0PT_ee",60, 0, 600);
	TH1F *h_Jet0PT_ttcc_emu= new TH1F("h_Jet0PT_ttcc_emu", "Jet0PT_emu",60, 0, 600);
	TH1F *h_Jet0PT_ttcc_mumu= new TH1F("h_Jet0PT_ttcc_mumu", "Jet0PT_mumu",60, 0, 600);

        TH1F *h_Jet0PT_SF= new TH1F("h_Jet0PT_SF", "Jet0PT_SF",60, 0, 600);
        TH1F *h_Jet0PT_SF_ee= new TH1F("h_Jet0PT_SF_ee", "Jet0PT_SF_ee",60, 0, 600);
        TH1F *h_Jet0PT_SF_emu= new TH1F("h_Jet0PT_SF_emu", "Jet0PT_SF_emu",60, 0, 600);
        TH1F *h_Jet0PT_SF_mumu= new TH1F("h_Jet0PT_SF_mumu", "Jet0PT_SF_mumu",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttlf= new TH1F("h_Jet0PT_SF_ttlf", "Jet0PT_SF",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttlf_ee= new TH1F("h_Jet0PT_SF_ttlf_ee", "Jet0PT_SF",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttlf_emu= new TH1F("h_Jet0PT_SF_ttlf_emu", "Jet0PT_SF",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttlf_mumu= new TH1F("h_Jet0PT_SF_ttlf_mumu", "Jet0PT_SF",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttb= new TH1F("h_Jet0PT_SF_ttb", "Jet0PT_SF",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttb_ee= new TH1F("h_Jet0PT_SF_ttb_ee", "Jet0PT_SF_ee",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttb_emu= new TH1F("h_Jet0PT_SF_ttb_emu", "Jet0PT_SF_emu",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttb_mumu= new TH1F("h_Jet0PT_SF_ttb_mumu", "Jet0PT_SF_mumu",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttbb= new TH1F("h_Jet0PT_SF_ttbb", "Jet0PT_SF",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttbb_ee= new TH1F("h_Jet0PT_SF_ttbb_ee", "Jet0PT_SF_ee",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttbb_emu= new TH1F("h_Jet0PT_SF_ttbb_emu", "Jet0PT_SF_emu",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttbb_mumu= new TH1F("h_Jet0PT_SF_ttbb_mumu", "Jet0PT_SF_mumu",60, 0, 600);
        TH1F *h_Jet0PT_SF_tt2b= new TH1F("h_Jet0PT_SF_tt2b", "Jet0PT_SF",60, 0, 600);
        TH1F *h_Jet0PT_SF_tt2b_ee= new TH1F("h_Jet0PT_SF_tt2b_ee", "Jet0PT_SF_ee",60, 0, 600);
        TH1F *h_Jet0PT_SF_tt2b_emu= new TH1F("h_Jet0PT_SF_tt2b_emu", "Jet0PT_SF_emu",60, 0, 600);
        TH1F *h_Jet0PT_SF_tt2b_mumu= new TH1F("h_Jet0PT_SF_tt2b_mumu", "Jet0PT_SF_mumu",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttcc= new TH1F("h_Jet0PT_SF_ttcc", "Jet0PT_SF",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttcc_ee= new TH1F("h_Jet0PT_SF_ttcc_ee", "Jet0PT_SF_ee",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttcc_emu= new TH1F("h_Jet0PT_SF_ttcc_emu", "Jet0PT_SF_emu",60, 0, 600);
        TH1F *h_Jet0PT_SF_ttcc_mumu= new TH1F("h_Jet0PT_SF_ttcc_mumu", "Jet0PT_SF_mumu",60, 0, 600);


//0000000000000000000000



  	TH1F *h_Jet0Eta= new TH1F("h_Jet0Eta", "Jet0Eta",30, -2.4, 2.4); 
	TH1F *h_Jet0Eta_ttlf= new TH1F("h_Jet0Eta_ttlf", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet0Eta_ttb= new TH1F("h_Jet0Eta_ttb", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet0Eta_ttbb= new TH1F("h_Jet0Eta_ttbb", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet0Eta_tt2b= new TH1F("h_Jet0Eta_tt2b", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet0Eta_ttcc= new TH1F("h_Jet0Eta_ttcc", "Jet0Eta",30, -2.4, 2.4);


	TH1F *h_Jet0Eta_SF= new TH1F("h_Jet0Eta_SF", "Jet0Eta",30, -2.4, 2.4);
	TH1F *h_Jet0Eta_SF_ttlf= new TH1F("h_Jet0Eta_SF_ttlf", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet0Eta_SF_ttb= new TH1F("h_Jet0Eta_SF_ttb", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet0Eta_SF_ttbb= new TH1F("h_Jet0Eta_SF_ttbb", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet0Eta_SF_tt2b= new TH1F("h_Jet0Eta_SF_tt2b", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet0Eta_SF_ttcc= new TH1F("h_Jet0Eta_SF_ttcc", "Jet0Eta",30, -2.4, 2.4);






 	TH1F *h_Jet0Phi= new TH1F("h_Jet0Phi", "Jet0Phi",40, -4, 4);
	TH1F *h_Jet0combinedInclusiveSecondaryVertexV2BJetTags = new TH1F("h_Jet0combinedInclusiveSecondaryVertexV2BJetTags","h_Jet0combinedInclusiveSecondaryVertexV2BJetTags",50,0,1);
 

  	TH1F *h_Jet0DeepCSV_b = new TH1F("h_Jet0DeepCSV_b","h_Jet0DeepCSV_b",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ee = new TH1F("h_Jet0DeepCSV_b_ee","h_Jet0DeepCSV_b_ee",50,0,1);
	TH1F *h_Jet0DeepCSV_b_emu = new TH1F("h_Jet0DeepCSV_b_emu","h_Jet0DeepCSV_b_emu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_mumu = new TH1F("h_Jet0DeepCSV_b_mumu","h_Jet0DeepCSV_b_mumu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttlf = new TH1F("h_Jet0DeepCSV_b_ttlf","h_Jet0DeepCSV_b",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttlf_ee = new TH1F("h_Jet0DeepCSV_b_ttlf_ee","h_Jet0DeepCSV_b_ee",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttlf_emu = new TH1F("h_Jet0DeepCSV_b_ttlf_emu","h_Jet0DeepCSV_b_emu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttlf_mumu = new TH1F("h_Jet0DeepCSV_b_ttlf_mumu","h_Jet0DeepCSV_b_mumu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttb = new TH1F("h_Jet0DeepCSV_b_ttb","h_Jet0DeepCSV_b",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttb_ee = new TH1F("h_Jet0DeepCSV_b_ttb_ee","h_Jet0DeepCSV_b_ee",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttb_emu = new TH1F("h_Jet0DeepCSV_b_ttb_emu","h_Jet0DeepCSV_b_emu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttb_mumu = new TH1F("h_Jet0DeepCSV_b_ttb_mumu","h_Jet0DeepCSV_b_mumu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttbb = new TH1F("h_Jet0DeepCSV_b_ttbb","h_Jet0DeepCSV_b",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttbb_ee = new TH1F("h_Jet0DeepCSV_b_ttbb_ee","h_Jet0DeepCSV_b_ee",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttbb_emu = new TH1F("h_Jet0DeepCSV_b_ttbb_emu","h_Jet0DeepCSV_b_emu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttbb_mumu = new TH1F("h_Jet0DeepCSV_b_ttbb_mumu","h_Jet0DeepCSV_b_mumu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_tt2b = new TH1F("h_Jet0DeepCSV_b_tt2b","h_Jet0DeepCSV_b",50,0,1);
	TH1F *h_Jet0DeepCSV_b_tt2b_ee = new TH1F("h_Jet0DeepCSV_b_tt2b_ee","h_Jet0DeepCSV_b_ee",50,0,1);
	TH1F *h_Jet0DeepCSV_b_tt2b_emu = new TH1F("h_Jet0DeepCSV_b_tt2b_emu","h_Jet0DeepCSV_b_emu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_tt2b_mumu = new TH1F("h_Jet0DeepCSV_b_tt2b_mumu","h_Jet0DeepCSV_b_mumu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttcc = new TH1F("h_Jet0DeepCSV_b_ttcc","h_Jet0DeepCSV_b",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttcc_ee = new TH1F("h_Jet0DeepCSV_b_ttcc_ee","h_Jet0DeepCSV_b_ee",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttcc_emu = new TH1F("h_Jet0DeepCSV_b_ttcc_emu","h_Jet0DeepCSV_b_emu",50,0,1);
	TH1F *h_Jet0DeepCSV_b_ttcc_mumu = new TH1F("h_Jet0DeepCSV_b_ttcc_mumu","h_Jet0DeepCSV_b_mumu",50,0,1);

	





	TH1F *h_Jet0DeepCSV = new TH1F("h_Jet0DeepCSV","h_Jet0DeepCSV",50,0,1);
        TH1F *h_Jet0DeepCSV_ttlf = new TH1F("h_Jet0DeepCSV_ttlf","h_Jet0DeepCSV",50,0,1);
        TH1F *h_Jet0DeepCSV_ttb = new TH1F("h_Jet0DeepCSV_ttb","h_Jet0DeepCSV",50,0,1);
        TH1F *h_Jet0DeepCSV_ttbb = new TH1F("h_Jet0DeepCSV_ttbb","h_Jet0DeepCSV",50,0,1);
        TH1F *h_Jet0DeepCSV_tt2b = new TH1F("h_Jet0DeepCSV_tt2b","h_Jet0DeepCSV",50,0,1);
        TH1F *h_Jet0DeepCSV_ttcc = new TH1F("h_Jet0DeepCSV_ttcc","h_Jet0DeepCSV",50,0,1);

	TH1F *h_Jet0DeepCSV_SF = new TH1F("h_Jet0DeepCSV_SF","h_Jet0DeepCSV",50,0,1);
        TH1F *h_Jet0DeepCSV_SF_ttlf = new TH1F("h_Jet0DeepCSV_SF_ttlf","h_Jet0DeepCSV",50,0,1);
        TH1F *h_Jet0DeepCSV_SF_ttb = new TH1F("h_Jet0DeepCSV_SF_ttb","h_Jet0DeepCSV",50,0,1);
        TH1F *h_Jet0DeepCSV_SF_ttbb = new TH1F("h_Jet0DeepCSV_SF_ttbb","h_Jet0DeepCSV",50,0,1);
        TH1F *h_Jet0DeepCSV_SF_tt2b = new TH1F("h_Jet0DeepCSV_SF_tt2b","h_Jet0DeepCSV",50,0,1);
        TH1F *h_Jet0DeepCSV_SF_ttcc = new TH1F("h_Jet0DeepCSV_SF_ttcc","h_Jet0DeepCSV",50,0,1);

	        TH1F *h_Jet1DeepCSV = new TH1F("h_Jet1DeepCSV","h_Jet1DeepCSV",50,0,1);
        TH1F *h_Jet1DeepCSV_ttlf = new TH1F("h_Jet1DeepCSV_ttlf","h_Jet1DeepCSV",50,0,1);
        TH1F *h_Jet1DeepCSV_ttb = new TH1F("h_Jet1DeepCSV_ttb","h_Jet1DeepCSV",50,0,1);
        TH1F *h_Jet1DeepCSV_ttbb = new TH1F("h_Jet1DeepCSV_ttbb","h_Jet1DeepCSV",50,0,1);
        TH1F *h_Jet1DeepCSV_tt2b = new TH1F("h_Jet1DeepCSV_tt2b","h_Jet1DeepCSV",50,0,1);
        TH1F *h_Jet1DeepCSV_ttcc = new TH1F("h_Jet1DeepCSV_ttcc","h_Jet1DeepCSV",50,0,1);

        TH1F *h_Jet1DeepCSV_SF = new TH1F("h_Jet1DeepCSV_SF","h_Jet1DeepCSV",50,0,1);
        TH1F *h_Jet1DeepCSV_SF_ttlf = new TH1F("h_Jet1DeepCSV_SF_ttlf","h_Jet1DeepCSV",50,0,1);
        TH1F *h_Jet1DeepCSV_SF_ttb = new TH1F("h_Jet1DeepCSV_SF_ttb","h_Jet1DeepCSV",50,0,1);
        TH1F *h_Jet1DeepCSV_SF_ttbb = new TH1F("h_Jet1DeepCSV_SF_ttbb","h_Jet1DeepCSV",50,0,1);
        TH1F *h_Jet1DeepCSV_SF_tt2b = new TH1F("h_Jet1DeepCSV_SF_tt2b","h_Jet1DeepCSV",50,0,1);
        TH1F *h_Jet1DeepCSV_SF_ttcc = new TH1F("h_Jet1DeepCSV_SF_ttcc","h_Jet1DeepCSV",50,0,1);







	TH1F *h_Jet0DeepCSV_b_SF = new TH1F("h_Jet0DeepCSV_b_SF","h_Jet0DeepCSV_b",50,0,1);
	TH1F *h_Jet0DeepCSV_b_SF_ttlf = new TH1F("h_Jet0DeepCSV_b_SF_ttlf","h_Jet0DeepCSV_b",50,0,1);
        TH1F *h_Jet0DeepCSV_b_SF_ttb = new TH1F("h_Jet0DeepCSV_b_SF_ttb","h_Jet0DeepCSV_b",50,0,1);
        TH1F *h_Jet0DeepCSV_b_SF_ttbb = new TH1F("h_Jet0DeepCSV_b_SF_ttbb","h_Jet0DeepCSV_b",50,0,1);
        TH1F *h_Jet0DeepCSV_b_SF_tt2b = new TH1F("h_Jet0DeepCSV_b_SF_tt2b","h_Jet0DeepCSV_b",50,0,1);
        TH1F *h_Jet0DeepCSV_b_SF_ttcc = new TH1F("h_Jet0DeepCSV_b_SF_ttcc","h_Jet0DeepCSV_b",50,0,1);

   	TH1F *h_Jet0SFDeepCSV= new TH1F("h_Jet0SFDeepCSV", "Jet0SFDeepCSV",20, 0.5, 1.5);
	TH1F *h_Jet0DeepCSV_bb = new TH1F("h_Jet0DeepCSV_bb","h_Jet0DeepCSV_bb",50,0,1);
	TH1F *h_Jet0DeepCSV_bb_ttlf = new TH1F("h_Jet0DeepCSV_bb_ttlf","h_Jet0DeepCSV_bb",50,0,1);
	TH1F *h_Jet0DeepCSV_bb_ttb = new TH1F("h_Jet0DeepCSV_bb_ttb","h_Jet0DeepCSV_bb",50,0,1);
	TH1F *h_Jet0DeepCSV_bb_ttbb = new TH1F("h_Jet0DeepCSV_bb_ttbb","h_Jet0DeepCSV_bb",50,0,1);
	TH1F *h_Jet0DeepCSV_bb_tt2b = new TH1F("h_Jet0DeepCSV_bb_tt2b","h_Jet0DeepCSV_bb",50,0,1);
	TH1F *h_Jet0DeepCSV_bb_ttcc = new TH1F("h_Jet0DeepCSV_bb_ttcc","h_Jet0DeepCSV_bb",50,0,1);
	


        TH1F *h_Jet0DeepCSV_bb_SF = new TH1F("h_Jet0DeepCSV_bb_SF","h_Jet0DeepCSV_bb",50,0,1);
        TH1F *h_Jet0DeepCSV_bb_SF_ttlf = new TH1F("h_Jet0DeepCSV_bb_SF_ttlf","h_Jet0DeepCSV_bb",50,0,1);
        TH1F *h_Jet0DeepCSV_bb_SF_ttb = new TH1F("h_Jet0DeepCSV_bb_SF_ttb","h_Jet0DeepCSV_bb",50,0,1);
        TH1F *h_Jet0DeepCSV_bb_SF_ttbb = new TH1F("h_Jet0DeepCSV_bb_SF_ttbb","h_Jet0DeepCSV_bb",50,0,1);
        TH1F *h_Jet0DeepCSV_bb_SF_tt2b = new TH1F("h_Jet0DeepCSV_bb_SF_tt2b","h_Jet0DeepCSV_bb",50,0,1);
        TH1F *h_Jet0DeepCSV_bb_SF_ttcc = new TH1F("h_Jet0DeepCSV_bb_SF_ttcc","h_Jet0DeepCSV_bb",50,0,1);





  	TH1F *h_Jet1PT= new TH1F("h_Jet1PT", "Jet1PT",60, 0, 600);
  	TH1F *h_Jet1PT_ee= new TH1F("h_Jet1PT_ee", "Jet1PT_ee",60, 0, 600);
	TH1F *h_Jet1PT_emu= new TH1F("h_Jet1PT_emu", "Jet1PT_emu",60, 0, 600);
	TH1F *h_Jet1PT_mumu= new TH1F("h_Jet1PT_mumu", "Jet1PT_mumu",60, 0, 600);
	TH1F *h_Jet1PT_ttlf= new TH1F("h_Jet1PT_ttlf", "Jet1PT",60, 0, 600);
	TH1F *h_Jet1PT_ttlf_ee= new TH1F("h_Jet1PT_ttlf_ee", "Jet1PT_ee",60, 0, 600);
	TH1F *h_Jet1PT_ttlf_emu= new TH1F("h_Jet1PT_ttlf_emu", "Jet1PT_emu",60, 0, 600);
	TH1F *h_Jet1PT_ttlf_mumu= new TH1F("h_Jet1PT_ttlf_mumu", "Jet1PT_mumu",60, 0, 600);
	TH1F *h_Jet1PT_ttb= new TH1F("h_Jet1PT_ttb", "Jet1PT",60, 0, 600);
	TH1F *h_Jet1PT_ttb_ee= new TH1F("h_Jet1PT_ttb_ee", "Jet1PT_ee",60, 0, 600);	
	TH1F *h_Jet1PT_ttb_emu= new TH1F("h_Jet1PT_ttb_emu", "Jet1PT_emu",60, 0, 600);
	TH1F *h_Jet1PT_ttb_mumu= new TH1F("h_Jet1PT_ttb_mumu", "Jet1PT_mumu",60, 0, 600);
	TH1F *h_Jet1PT_ttbb= new TH1F("h_Jet1PT_ttbb", "Jet1PT",60, 0, 600);
	TH1F *h_Jet1PT_ttbb_ee= new TH1F("h_Jet1PT_ttbb_ee", "Jet1PT_ee",60, 0, 600);
	TH1F *h_Jet1PT_ttbb_emu= new TH1F("h_Jet1PT_ttbb_emu", "Jet1PT_emu",60, 0, 600);
	TH1F *h_Jet1PT_ttbb_mumu= new TH1F("h_Jet1PT_ttbb_mumu", "Jet1PT_mumu",60, 0, 600);
	TH1F *h_Jet1PT_tt2b= new TH1F("h_Jet1PT_tt2b", "Jet1PT",60, 0, 600);
	TH1F *h_Jet1PT_tt2b_ee= new TH1F("h_Jet1PT_tt2b_ee", "Jet1PT_ee",60, 0, 600);
	TH1F *h_Jet1PT_tt2b_emu= new TH1F("h_Jet1PT_tt2b_emu", "Jet1PT_emu",60, 0, 600);
	TH1F *h_Jet1PT_tt2b_mumu= new TH1F("h_Jet1PT_tt2b_mumu", "Jet1PT_mumu",60, 0, 600);
	TH1F *h_Jet1PT_ttcc= new TH1F("h_Jet1PT_ttcc", "Jet1PT",60, 0, 600);
	TH1F *h_Jet1PT_ttcc_ee= new TH1F("h_Jet1PT_ttcc_ee", "Jet1PT_ee",60, 0, 600);
	TH1F *h_Jet1PT_ttcc_emu= new TH1F("h_Jet1PT_ttcc_emu", "Jet1PT_emu",60, 0, 600);
	TH1F *h_Jet1PT_ttcc_mumu= new TH1F("h_Jet1PT_ttcc_mumu", "Jet1PT_mumu",60, 0, 600);

        TH1F *h_Jet1PT_SF = new TH1F("h_Jet1PT_SF", "Jet1PT_SF",60, 0, 600);
        TH1F *h_Jet1PT_SF_ee= new TH1F("h_Jet1PT_SF_ee", "Jet1PT_SF_ee",60, 0, 600);
        TH1F *h_Jet1PT_SF_emu= new TH1F("h_Jet1PT_SF_emu", "Jet1PT_SF_emu",60, 0, 600);
        TH1F *h_Jet1PT_SF_mumu= new TH1F("h_Jet1PT_SF_mumu", "Jet1PT_SF_mumu",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttlf= new TH1F("h_Jet1PT_SF_ttlf", "Jet1PT_SF",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttlf_ee= new TH1F("h_Jet1PT_SF_ttlf_ee", "Jet1PT_SF_ee",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttlf_emu= new TH1F("h_Jet1PT_SF_ttlf_emu", "Jet1PT_SF_emu",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttlf_mumu= new TH1F("h_Jet1PT_SF_ttlf_mumu", "Jet1PT_SF_mumu",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttb= new TH1F("h_Jet1PT_SF_ttb", "Jet1PT_SF",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttb_ee= new TH1F("h_Jet1PT_SF_ttb_ee", "Jet1PT_SF_ee",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttb_emu= new TH1F("h_Jet1PT_SF_ttb_emu", "Jet1PT_SF_emu",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttb_mumu= new TH1F("h_Jet1PT_SF_ttb_mumu", "Jet1PT_SF_mumu",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttbb= new TH1F("h_Jet1PT_SF_ttbb", "Jet1PT_SF",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttbb_ee= new TH1F("h_Jet1PT_SF_ttbb_ee", "Jet1PT_SF_ee",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttbb_emu= new TH1F("h_Jet1PT_SF_ttbb_emu", "Jet1PT_SF_emu",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttbb_mumu= new TH1F("h_Jet1PT_SF_ttbb_mumu", "Jet1PT_SF_mumu",60, 0, 600);
        TH1F *h_Jet1PT_SF_tt2b= new TH1F("h_Jet1PT_SF_tt2b", "Jet1PT_SF",60, 0, 600);
        TH1F *h_Jet1PT_SF_tt2b_ee= new TH1F("h_Jet1PT_SF_tt2b_ee", "Jet1PT_SF_ee",60, 0, 600);
        TH1F *h_Jet1PT_SF_tt2b_emu= new TH1F("h_Jet1PT_SF_tt2b_emu", "Jet1PT_SF_emu",60, 0, 600);
        TH1F *h_Jet1PT_SF_tt2b_mumu= new TH1F("h_Jet1PT_SF_tt2b_mumu", "Jet1PT_SF_mumu",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttcc= new TH1F("h_Jet1PT_SF_ttcc", "Jet1PT_SF",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttcc_ee= new TH1F("h_Jet1PT_SF_ttcc_ee", "Jet1PT_SF_ee",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttcc_emu= new TH1F("h_Jet1PT_SF_ttcc_emu", "Jet1PT_SF_emu",60, 0, 600);
        TH1F *h_Jet1PT_SF_ttcc_mumu= new TH1F("h_Jet1PT_SF_ttcc_mumu", "Jet1PT_SF_mumu",60, 0, 600);



//000000000000000000000000000000000
	TH1F *h_Jet1Eta= new TH1F("h_Jet1Eta", "Jet1Eta",30, -2.4, 2.4);
	TH1F *h_Jet1Eta_ttlf= new TH1F("h_Jet1Eta_ttlf", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet1Eta_ttb= new TH1F("h_Jet1Eta_ttb", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet1Eta_ttbb= new TH1F("h_Jet1Eta_ttbb", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet1Eta_tt2b= new TH1F("h_Jet1Eta_tt2b", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet1Eta_ttcc= new TH1F("h_Jet1Eta_ttcc", "Jet0Eta",30, -2.4, 2.4);

	TH1F *h_Jet1Eta_SF= new TH1F("h_Jet1Eta_SF", "Jet1Eta",30, -2.4, 2.4);
        TH1F *h_Jet1Eta_SF_ttlf= new TH1F("h_Jet1Eta_SF_ttlf", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet1Eta_SF_ttb= new TH1F("h_Jet1Eta_SF_ttb", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet1Eta_SF_ttbb= new TH1F("h_Jet1Eta_SF_ttbb", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet1Eta_SF_tt2b= new TH1F("h_Jet1Eta_SF_tt2b", "Jet0Eta",30, -2.4, 2.4);
        TH1F *h_Jet1Eta_SF_ttcc= new TH1F("h_Jet1Eta_SF_ttcc", "Jet0Eta",30, -2.4, 2.4);




 	TH1F *h_Jet1Phi= new TH1F("h_Jet1Phi", "Jet1Phi",40, -4, 4); 
  	TH1F *h_Jet1combinedInclusiveSecondaryVertexV2BJetTags = new TH1F("h_Jet1combinedInclusiveSecondaryVertexV2BJetTags","h_Jet1combinedInclusiveSecondaryVertexV2BJetTags",50,0,1);
   	TH1F *h_Jet1DeepCSV_b = new TH1F("h_Jet1DeepCSV_b","h_Jet1DeepCSV_b",50,0,1);
	TH1F *h_Jet1DeepCSV_b_ttlf = new TH1F("h_Jet1DeepCSV_b_ttlf","h_Jet1DeepCSV_b",50,0,1);
	TH1F *h_Jet1DeepCSV_b_ttb = new TH1F("h_Jet1DeepCSV_b_ttb","h_Jet1DeepCSV_b",50,0,1);
	TH1F *h_Jet1DeepCSV_b_ttbb = new TH1F("h_Jet1DeepCSV_b_ttbb","h_Jet1DeepCSV_b",50,0,1);
	TH1F *h_Jet1DeepCSV_b_tt2b = new TH1F("h_Jet1DeepCSV_b_tt2b","h_Jet1DeepCSV_b",50,0,1);
	TH1F *h_Jet1DeepCSV_b_ttcc = new TH1F("h_Jet1DeepCSV_b_ttcc","h_Jet1DeepCSV_b",50,0,1);

	TH1F *h_Jet1DeepCSV_b_SF = new TH1F("h_Jet1DeepCSV_b_SF","h_Jet1DeepCSV_b",50,0,1);
        TH1F *h_Jet1DeepCSV_b_SF_ttlf = new TH1F("h_Jet1DeepCSV_b_SF_ttlf","h_Jet1DeepCSV_b",50,0,1);
        TH1F *h_Jet1DeepCSV_b_SF_ttb = new TH1F("h_Jet1DeepCSV_b_SF_ttb","h_Jet1DeepCSV_b",50,0,1);
        TH1F *h_Jet1DeepCSV_b_SF_ttbb = new TH1F("h_Jet1DeepCSV_b_SF_ttbb","h_Jet1DeepCSV_b",50,0,1);
        TH1F *h_Jet1DeepCSV_b_SF_tt2b = new TH1F("h_Jet1DeepCSV_b_SF_tt2b","h_Jet1DeepCSV_b",50,0,1);
        TH1F *h_Jet1DeepCSV_b_SF_ttcc = new TH1F("h_Jet1DeepCSV_b_SF_ttcc","h_Jet1DeepCSV_b",50,0,1);

	TH1F *h_Jet1SFDeepCSV= new TH1F("h_Jet1SFDeepCSV", "Jet1SFDeepCSV",20, 0.5, 1.5);
   	TH1F *h_Jet1DeepCSV_bb = new TH1F("h_Jet1DeepCSV_bb","h_Jet1DeepCSV_bb",50,0,1);
	 TH1F *h_Jet1DeepCSV_bb_ttlf = new TH1F("h_Jet1DeepCSV_bb_ttlf","h_Jet1DeepCSV_bb",50,0,1);
	 TH1F *h_Jet1DeepCSV_bb_ttb = new TH1F("h_Jet1DeepCSV_bb_ttb","h_Jet1DeepCSV_bb",50,0,1);
	 TH1F *h_Jet1DeepCSV_bb_ttbb = new TH1F("h_Jet1DeepCSV_bb_ttbb","h_Jet1DeepCSV_bb",50,0,1);
	 TH1F *h_Jet1DeepCSV_bb_tt2b = new TH1F("h_Jet1DeepCSV_bb_tt2b","h_Jet1DeepCSV_bb",50,0,1);
	 TH1F *h_Jet1DeepCSV_bb_ttcc = new TH1F("h_Jet1DeepCSV_bb_ttcc","h_Jet1DeepCSV_bb",50,0,1);

        TH1F *h_Jet1DeepCSV_bb_SF = new TH1F("h_Jet1DeepCSV_bb_SF","h_Jet1DeepCSV_bb",50,0,1);
         TH1F *h_Jet1DeepCSV_bb_SF_ttlf = new TH1F("h_Jet1DeepCSV_bb_SF_ttlf","h_Jet1DeepCSV_bb",50,0,1);
         TH1F *h_Jet1DeepCSV_bb_SF_ttb = new TH1F("h_Jet1DeepCSV_bb_SF_ttb","h_Jet1DeepCSV_bb",50,0,1);
         TH1F *h_Jet1DeepCSV_bb_SF_ttbb = new TH1F("h_Jet1DeepCSV_bb_SF_ttbb","h_Jet1DeepCSV_bb",50,0,1);
         TH1F *h_Jet1DeepCSV_bb_SF_tt2b = new TH1F("h_Jet1DeepCSV_bb_SF_tt2b","h_Jet1DeepCSV_bb",50,0,1);
         TH1F *h_Jet1DeepCSV_bb_SF_ttcc = new TH1F("h_Jet1DeepCSV_bb_SF_ttcc","h_Jet1DeepCSV_bb",50,0,1);





  	TH1F *h_Jet2PT= new TH1F("h_Jet2PT", "Jet2PT",60, 0, 600);
  	TH1F *h_Jet2Eta= new TH1F("h_Jet2Eta", "Jet2Eta",30, -2.4, 2.4);
 	TH1F *h_Jet2Phi= new TH1F("h_Jet2Phi", "Jet2Phi",40, -4, 4);
   	TH1F *h_Jet2combinedInclusiveSecondaryVertexV2BJetTags = new TH1F("h_Jet2combinedInclusiveSecondaryVertexV2BJetTags","h_Jet2combinedInclusiveSecondaryVertexV2BJetTags",50,0,1);
   	TH1F *h_Jet2DeepCSV_b = new TH1F("h_Jet2DeepCSV_b","h_Jet2DeepCSV_b",50,0,1);
   	TH1F *h_Jet2SFDeepCSV= new TH1F("h_Jet2SFDeepCSV", "Jet2SFDeepCSV",20, 0.5, 1.5);
	TH1F *h_Jet2DeepCSV_bb = new TH1F("h_Jet2DeepCSV_bb","h_Jet2DeepCSV_bb",50,0,1);

  	TH1F *h_Jet3PT= new TH1F("h_Jet3PT", "Jet3PT",60, 0, 600);
  	TH1F *h_Jet3Eta= new TH1F("h_Jet3Eta", "Jet3Eta",30, -2.4, 2.4);
 	TH1F *h_Jet3Phi= new TH1F("h_Jet3Phi", "Jet3Phi",40, -4, 4);
    	TH1F *h_Jet3combinedInclusiveSecondaryVertexV2BJetTags = new TH1F("h_Jet3combinedInclusiveSecondaryVertexV2BJetTags","h_Jet3combinedInclusiveSecondaryVertexV2BJetTags",50,0,1);
   	TH1F *h_Jet3DeepCSV_b = new TH1F("h_Jet3DeepCSV_b","h_Jet3DeepCSV_b",50,0,1);
   	TH1F *h_Jet3SFDeepCSV= new TH1F("h_Jet3SFDeepCSV", "Jet3SFDeepCSV",20, 0.5, 1.5);
	TH1F *h_Jet3DeepCSV_bb = new TH1F("h_Jet3DeepCSV_bb","h_Jet3DeepCSV_bb",50,0,1);

  	TH1F *h_Jet4PT= new TH1F("h_Jet4PT", "Jet4PT",60, 0, 600);
  	TH1F *h_Jet4Eta= new TH1F("h_Jet4Eta", "Jet4Eta",30, -2.4, 2.4);
 	TH1F *h_Jet4Phi= new TH1F("h_Jet4Phi", "Jet4Phi",30, -3, 3);
    	TH1F *h_Jet4combinedInclusiveSecondaryVertexV2BJetTags = new TH1F("h_Jet4combinedInclusiveSecondaryVertexV2BJetTags","h_Jet4combinedInclusiveSecondaryVertexV2BJetTags",50,0,1);
   	TH1F *h_Jet4DeepCSV_b = new TH1F("h_Jet4DeepCSV_b","h_Jet4DeepCSV_b",50,0,1);
   	TH1F *h_Jet4DeepCSV_bb = new TH1F("h_Jet4DeepCSV_bb","h_Jet4DeepCSV_bb",50,0,1);

  	TH1F *h_Jet5PT= new TH1F("h_Jet5PT", "Jet5PT",60, 0, 600);
  	TH1F *h_Jet5Eta= new TH1F("h_Jet5Eta", "Jet5Eta",30, -2.4, 2.4);
 	TH1F *h_Jet5Phi= new TH1F("h_Jet5Phi", "Jet5Phi",40, -4, 4);
    	TH1F *h_Jet5combinedInclusiveSecondaryVertexV2BJetTags = new TH1F("h_Jet5combinedInclusiveSecondaryVertexV2BJetTags","h_Jet5combinedInclusiveSecondaryVertexV2BJetTags",50,0,1);
   	TH1F *h_Jet5DeepCSV_b = new TH1F("h_Jet5DeepCSV_b","h_Jet5DeepCSV_b",50,0,1);
   	TH1F *h_Jet5DeepCSV_bb = new TH1F("h_Jet5DeepCSV_bb","h_Jet5DeepCSV_bb",50,0,1);

  	TH1F *h_Jet6PT= new TH1F("h_Jet6PT", "Jet6PT",60, 0, 600);
    	TH1F *h_Jet6combinedInclusiveSecondaryVertexV2BJetTags = new TH1F("h_Jet6combinedInclusiveSecondaryVertexV2BJetTags","h_Jet6combinedInclusiveSecondaryVertexV2BJetTags",50,0,1);
  	TH1F *h_Jet6Eta= new TH1F("h_Jet6Eta", "Jet6Eta",30, -2.4, 2.4);
 	TH1F *h_Jet6Phi= new TH1F("h_Jet6Phi", "Jet6Phi",40, -4, 4);
   	TH1F *h_Jet6DeepCSV_b = new TH1F("h_Jet6DeepCSV_b","h_Jet6DeepCSV_b",50,0,1);
   	TH1F *h_Jet6DeepCSV_bb = new TH1F("h_Jet6DeepCSV_bb","h_Jet6DeepCSV_bb",50,0,1);

	double Muon0PT = -99;double Muon0Eta = -99;

  	TH1F *h_Muon0PT = new TH1F("h_Muon0PT", "Muon0PT",60, 0, 300);
  	TH1F *h_Muon0PT_SF = new TH1F("h_Muon0PT_SF", "Muon0PT_SF",60, 0, 300);
	TH1F *h_Muon0Eta = new TH1F("h_Muon0Eta", "Muon0Eta",30, -2.4, 2.4);
  	TH1F *h_Muon0Eta_SF = new TH1F("h_Muon0Eta_SF", "Muon0Eta_SF",30, -2.4, 2.4);
	TH1F *h_Muon0Phi = new TH1F("h_Muon0Phi", "Muon0Phi",40, -4, 4);  
	TH1F *h_Muon0Iso = new TH1F("h_Muon0Iso", "Muon0Iso",30, 0, 0.3);	

	double Ele0PT = -99;double Ele0Eta = -99;

  	TH1F *h_Ele0PT = new TH1F("h_Ele0PT", "Ele0PT",60, 0, 300);
  	TH1F *h_Ele0PT_SF = new TH1F("h_Ele0PT_SF", "Ele0PT_SF",60, 0, 300);
	TH1F *h_Ele0Eta = new TH1F("h_Ele0Eta", "Ele0Eta",30, -2.4, 2.4);
	TH1F *h_Ele0Eta_SF = new TH1F("h_Ele0Eta_SF", "Ele0Eta_SF",30, -2.4, 2.4);
	TH1F *h_Ele0Phi = new TH1F("h_Ele0Phi", "Ele0Phi",40, -4, 4);
	TH1F *h_Ele0Iso = new TH1F("h_Ele0Iso", "Ele0Iso",30, 0, 0.3);

	TH1F *h_Lep0PT = new TH1F("h_Lep0PT", "Lep0PT",60, 0, 300);
	TH1F *h_Lep0PT_ee = new TH1F("h_Lep0PT_ee", "Lep0PT_ee",60, 0, 300);
	TH1F *h_Lep0PT_emu = new TH1F("h_Lep0PT_emu", "Lep0PT_emu",60, 0, 300);
	TH1F *h_Lep0PT_mumu = new TH1F("h_Lep0PT_mumu", "Lep0PT_mumu",60, 0, 300);
	TH1F *h_Lep0PT_ttlf = new TH1F("h_Lep0PT_ttlf", "Lep0PT",60, 0, 300);
	TH1F *h_Lep0PT_ttlf_ee = new TH1F("h_Lep0PT_ttlf_ee", "Lep0PT_ee",60, 0, 300);
	TH1F *h_Lep0PT_ttlf_emu = new TH1F("h_Lep0PT_ttlf_emu", "Lep0PT_emu",60, 0, 300);
	TH1F *h_Lep0PT_ttlf_mumu = new TH1F("h_Lep0PT_ttlf_mumu", "Lep0PT_mumu",60, 0, 300);
	TH1F *h_Lep0PT_ttb = new TH1F("h_Lep0PT_ttb", "Lep0PT",60, 0, 300);
	TH1F *h_Lep0PT_ttb_ee = new TH1F("h_Lep0PT_ttb_ee", "Lep0PT_ee",60, 0, 300);
	TH1F *h_Lep0PT_ttb_emu = new TH1F("h_Lep0PT_ttb_emu", "Lep0PT_emu",60, 0, 300);
	TH1F *h_Lep0PT_ttb_mumu = new TH1F("h_Lep0PT_ttb_mumu", "Lep0PT_mumu",60, 0, 300);
	TH1F *h_Lep0PT_ttbb = new TH1F("h_Lep0PT_ttbb", "Lep0PT",60, 0, 300);
	TH1F *h_Lep0PT_ttbb_ee = new TH1F("h_Lep0PT_ttbb_ee", "Lep0PT_ee",60, 0, 300);
	TH1F *h_Lep0PT_ttbb_emu = new TH1F("h_Lep0PT_ttbb_emu", "Lep0PT_emu",60, 0, 300);
	TH1F *h_Lep0PT_ttbb_mumu = new TH1F("h_Lep0PT_ttbb_mumu", "Lep0PT_mumu",60, 0, 300);
	TH1F *h_Lep0PT_tt2b = new TH1F("h_Lep0PT_tt2b", "Lep0PT",60, 0, 300);
	TH1F *h_Lep0PT_tt2b_ee = new TH1F("h_Lep0PT_tt2b_ee", "Lep0PT_ee",60, 0, 300);
	TH1F *h_Lep0PT_tt2b_emu = new TH1F("h_Lep0PT_tt2b_emu", "Lep0PT_emu",60, 0, 300);
	TH1F *h_Lep0PT_tt2b_mumu = new TH1F("h_Lep0PT_tt2b_mumu", "Lep0PT_mumu",60, 0, 300);
	TH1F *h_Lep0PT_ttcc = new TH1F("h_Lep0PT_ttcc", "Lep0PT",60, 0, 300);
	TH1F *h_Lep0PT_ttcc_ee = new TH1F("h_Lep0PT_ttcc_ee", "Lep0PT_ee",60, 0, 300);
	TH1F *h_Lep0PT_ttcc_emu = new TH1F("h_Lep0PT_ttcc_emu", "Lep0PT_emu",60, 0, 300);
	TH1F *h_Lep0PT_ttcc_mumu = new TH1F("h_Lep0PT_ttcc_mumu", "Lep0PT_mumu",60, 0, 300);
        TH1F *h_Lep0PT_SF = new TH1F("h_Lep0PT_SF", "Lep0PT_SF",60, 0, 300);
	TH1F *h_Lep0PT_SF_ee = new TH1F("h_Lep0PT_SF_ee", "Lep0PT_SF_ee",60, 0, 300);
	TH1F *h_Lep0PT_SF_emu = new TH1F("h_Lep0PT_SF_emu", "Lep0PT_SF_emu",60, 0, 300);
	TH1F *h_Lep0PT_SF_mumu = new TH1F("h_Lep0PT_SF_mumu", "Lep0PT_SF_mumu",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttlf = new TH1F("h_Lep0PT_SF_ttlf", "Lep0PT_SF",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttlf_ee = new TH1F("h_Lep0PT_SF_ttlf_ee", "Lep0PT_SF_ee",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttlf_emu = new TH1F("h_Lep0PT_SF_ttlf_emu", "Lep0PT_SF_emu",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttlf_mumu = new TH1F("h_Lep0PT_SF_ttlf_mumu", "Lep0PT_SF_mumu",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttb = new TH1F("h_Lep0PT_SF_ttb", "Lep0PT_SF",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttb_ee = new TH1F("h_Lep0PT_SF_ttb_ee", "Lep0PT_SF_ee",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttb_emu = new TH1F("h_Lep0PT_SF_ttb_emu", "Lep0PT_SF_emu",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttb_mumu = new TH1F("h_Lep0PT_SF_ttb_mumu", "Lep0PT_SF_mumu",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttbb = new TH1F("h_Lep0PT_SF_ttbb", "Lep0PT_SF",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttbb_ee = new TH1F("h_Lep0PT_SF_ttbb_ee", "Lep0PT_SF_ee",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttbb_emu = new TH1F("h_Lep0PT_SF_ttbb_emu", "Lep0PT_SF_emu",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttbb_mumu = new TH1F("h_Lep0PT_SF_ttbb_mumu", "Lep0PT_SF_mumu",60, 0, 300);
	TH1F *h_Lep0PT_SF_tt2b = new TH1F("h_Lep0PT_SF_tt2b", "Lep0PT_SF",60, 0, 300);
	TH1F *h_Lep0PT_SF_tt2b_ee = new TH1F("h_Lep0PT_SF_tt2b_ee", "Lep0PT_SF_ee",60, 0, 300);
	TH1F *h_Lep0PT_SF_tt2b_emu = new TH1F("h_Lep0PT_SF_tt2b_emu", "Lep0PT_SF_emu",60, 0, 300);
	TH1F *h_Lep0PT_SF_tt2b_mumu = new TH1F("h_Lep0PT_SF_tt2b_mumu", "Lep0PT_SF_mumu",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttcc = new TH1F("h_Lep0PT_SF_ttcc", "Lep0PT_SF",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttcc_ee = new TH1F("h_Lep0PT_SF_ttcc_ee", "Lep0PT_SF_ee",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttcc_emu = new TH1F("h_Lep0PT_SF_ttcc_emu", "Lep0PT_SF_emu",60, 0, 300);
	TH1F *h_Lep0PT_SF_ttcc_mumu = new TH1F("h_Lep0PT_SF_ttcc_mumu", "Lep0PT_SF_mumu",60, 0, 300);





	TH1F *h_Lep0Eta = new TH1F("h_Lep0Eta", "Lep0Eta",30, -2.4, 2.4);
	TH1F *h_Lep0Eta_ttlf = new TH1F("h_Lep0Eta_ttlf", "Lep0Eta",30, -2.4, 2.4);        
	TH1F *h_Lep0Eta_ttb = new TH1F("h_Lep0Eta_ttb", "Lep0Eta",30, -2.4, 2.4);
	TH1F *h_Lep0Eta_ttbb = new TH1F("h_Lep0Eta_ttbb", "Lep0Eta",30, -2.4, 2.4);
	TH1F *h_Lep0Eta_tt2b = new TH1F("h_Lep0Eta_tt2b", "Lep0Eta",30, -2.4, 2.4);
	TH1F *h_Lep0Eta_ttcc = new TH1F("h_Lep0Eta_ttcc", "Lep0Eta",30, -2.4, 2.4);

	TH1F *h_Lep0Eta_SF = new TH1F("h_Lep0Eta_SF", "Lep0Eta_SF",30, -2.4, 2.4);
        TH1F *h_Lep0Eta_SF_ttlf = new TH1F("h_Lep0Eta_SF_ttlf", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep0Eta_SF_ttb = new TH1F("h_Lep0Eta_SF_ttb", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep0Eta_SF_ttbb = new TH1F("h_Lep0Eta_SF_ttbb", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep0Eta_SF_tt2b = new TH1F("h_Lep0Eta_SF_tt2b", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep0Eta_SF_ttcc = new TH1F("h_Lep0Eta_SF_ttcc", "Lep0Eta",30, -2.4, 2.4);





	TH1F *h_Lep0Phi = new TH1F("h_Lep0Phi", "Lep0Phi",40, -4, 4);
        TH1F *h_Lep0Iso = new TH1F("h_Lep0Iso", "Lep0Iso",30, 0, 0.3);
	TH1F *h_Lep0IsoSF = new TH1F("h_Lep0IsoSF", "Lep0IsoSF",20, 0.9, 1.1);	
	TH1F *h_Lep0IDSF = new TH1F("h_Lep0IDSF", "Lep0IDSF",20, 0.9, 1.1);

        TH1F *h_Lep1PT = new TH1F("h_Lep1PT", "Lep1PT",60, 0, 300);
	TH1F *h_Lep1PT_ee = new TH1F("h_Lep1PT_ee", "Lep1PT_ee",60, 0, 300);
	TH1F *h_Lep1PT_emu = new TH1F("h_Lep1PT_emu", "Lep1PT_emu",60, 0, 300);
	TH1F *h_Lep1PT_mumu = new TH1F("h_Lep1PT_mumu", "Lep1PT_mumu",60, 0, 300);
	TH1F *h_Lep1PT_ttlf = new TH1F("h_Lep1PT_ttlf", "Lep1PT",60, 0, 300);
	TH1F *h_Lep1PT_ttlf_ee = new TH1F("h_Lep1PT_ttlf_ee", "Lep1PT_ee",60, 0, 300);
	TH1F *h_Lep1PT_ttlf_emu = new TH1F("h_Lep1PT_ttlf_emu", "Lep1PT_emu",60, 0, 300);
	TH1F *h_Lep1PT_ttlf_mumu = new TH1F("h_Lep1PT_ttlf_mumu", "Lep1PT_mumu",60, 0, 300);
	TH1F *h_Lep1PT_ttb = new TH1F("h_Lep1PT_ttb", "Lep1PT",60, 0, 300);
	TH1F *h_Lep1PT_ttb_ee = new TH1F("h_Lep1PT_ttb_ee", "Lep1PT_ee",60, 0, 300);
	TH1F *h_Lep1PT_ttb_emu = new TH1F("h_Lep1PT_ttb_emu", "Lep1PT_emu",60, 0, 300);
	TH1F *h_Lep1PT_ttb_mumu = new TH1F("h_Lep1PT_ttb_mumu", "Lep1PT_mumu",60, 0, 300);
	TH1F *h_Lep1PT_ttbb = new TH1F("h_Lep1PT_ttbb", "Lep1PT",60, 0, 300);
	TH1F *h_Lep1PT_ttbb_ee = new TH1F("h_Lep1PT_ttbb_ee", "Lep1PT_ee",60, 0, 300);
	TH1F *h_Lep1PT_ttbb_emu = new TH1F("h_Lep1PT_ttbb_emu", "Lep1PT_emu",60, 0, 300);
	TH1F *h_Lep1PT_ttbb_mumu = new TH1F("h_Lep1PT_ttbb_mumu", "Lep1PT_mumu",60, 0, 300);
	TH1F *h_Lep1PT_tt2b = new TH1F("h_Lep1PT_tt2b", "Lep1PT",60, 0, 300);
	TH1F *h_Lep1PT_tt2b_ee = new TH1F("h_Lep1PT_tt2b_ee", "Lep1PT_ee",60, 0, 300);
	TH1F *h_Lep1PT_tt2b_emu = new TH1F("h_Lep1PT_tt2b_emu", "Lep1PT_emu",60, 0, 300);
	TH1F *h_Lep1PT_tt2b_mumu = new TH1F("h_Lep1PT_tt2b_mumu", "Lep1PT_mumu",60, 0, 300);
	TH1F *h_Lep1PT_ttcc = new TH1F("h_Lep1PT_ttcc", "Lep1PT",60, 0, 300);
	TH1F *h_Lep1PT_ttcc_ee = new TH1F("h_Lep1PT_ttcc_ee", "Lep1PT_ee",60, 0, 300);
	TH1F *h_Lep1PT_ttcc_emu = new TH1F("h_Lep1PT_ttcc_emu", "Lep1PT_emu",60, 0, 300);
	TH1F *h_Lep1PT_ttcc_mumu = new TH1F("h_Lep1PT_ttcc_mumu", "Lep1PT_mumu",60, 0, 300);



        TH1F *h_Lep1PT_SF = new TH1F("h_Lep1PT_SF", "Lep1PT_SF",60, 0, 300);
	TH1F *h_Lep1PT_SF_ee = new TH1F("h_Lep1PT_SF_ee", "Lep1PT_SF_ee",60, 0, 300);
	TH1F *h_Lep1PT_SF_emu = new TH1F("h_Lep1PT_SF_emu", "Lep1PT_SF_emu",60, 0, 300);
	TH1F *h_Lep1PT_SF_mumu = new TH1F("h_Lep1PT_SF_mumu", "Lep1PT_SF_mumu",60, 0, 300);
 	TH1F *h_Lep1PT_SF_ttlf = new TH1F("h_Lep1PT_SF_ttlf", "Lep1PT_SF_ttlf",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttlf_ee = new TH1F("h_Lep1PT_SF_ttlf_ee", "Lep1PT_SF_ttlf_ee",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttlf_emu = new TH1F("h_Lep1PT_SF_ttlf_emu", "Lep1PT_SF_ttlf_emu",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttlf_mumu = new TH1F("h_Lep1PT_SF_ttlf_mumu", "Lep1PT_SF_ttlf_mumu",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttb = new TH1F("h_Lep1PT_SF_ttb", "Lep1PT_SF_ttb",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttb_ee = new TH1F("h_Lep1PT_SF_ttb_ee", "Lep1PT_SF_ttb_ee",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttb_emu = new TH1F("h_Lep1PT_SF_ttb_emu", "Lep1PT_SF_ttb_emu",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttb_mumu = new TH1F("h_Lep1PT_SF_ttb_mumu", "Lep1PT_SF_ttb_mumu",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttbb = new TH1F("h_Lep1PT_SF_ttbb", "Lep1PT_SF_ttbb",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttbb_ee = new TH1F("h_Lep1PT_SF_ttbb_ee", "Lep1PT_SF_ttbb_ee",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttbb_emu = new TH1F("h_Lep1PT_SF_ttbb_emu", "Lep1PT_SF_ttbb_emu",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttbb_mumu = new TH1F("h_Lep1PT_SF_ttbb_mumu", "Lep1PT_SF_ttbb_mumu",60, 0, 300);
	TH1F *h_Lep1PT_SF_tt2b = new TH1F("h_Lep1PT_SF_tt2b", "Lep1PT_SF_tt2b",60, 0, 300);
	TH1F *h_Lep1PT_SF_tt2b_ee = new TH1F("h_Lep1PT_SF_tt2b_ee", "Lep1PT_SF_tt2b_ee",60, 0, 300);
	TH1F *h_Lep1PT_SF_tt2b_emu = new TH1F("h_Lep1PT_SF_tt2b_emu", "Lep1PT_SF_tt2b_emu",60, 0, 300);
	TH1F *h_Lep1PT_SF_tt2b_mumu = new TH1F("h_Lep1PT_SF_tt2b_mumu", "Lep1PT_SF_tt2b_mumu",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttcc = new TH1F("h_Lep1PT_SF_ttcc", "Lep1PT_SF_ttcc",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttcc_ee = new TH1F("h_Lep1PT_SF_ttcc_ee", "Lep1PT_SF_ttcc_ee",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttcc_emu = new TH1F("h_Lep1PT_SF_ttcc_emu", "Lep1PT_SF_ttcc_emu",60, 0, 300);
	TH1F *h_Lep1PT_SF_ttcc_mumu = new TH1F("h_Lep1PT_SF_ttcc_mumu", "Lep1PT_SF_ttcc_mumu",60, 0, 300);


        TH1F *h_Lep1Eta = new TH1F("h_Lep1Eta", "Lep1Eta",30, -2.4, 2.4);
        TH1F *h_Lep1Eta_ttlf = new TH1F("h_Lep1Eta_ttlf", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep1Eta_ttb = new TH1F("h_Lep1Eta_ttb", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep1Eta_ttbb = new TH1F("h_Lep1Eta_ttbb", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep1Eta_tt2b = new TH1F("h_Lep1Eta_tt2b", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep1Eta_ttcc = new TH1F("h_Lep1Eta_ttcc", "Lep0Eta",30, -2.4, 2.4);
	TH1F *h_Lep1Eta_SF = new TH1F("h_Lep1Eta_SF", "Lep1Eta_SF",30, -2.4, 2.4);
	TH1F *h_Lep1Eta_SF_ttlf = new TH1F("h_Lep1Eta_SF_ttlf", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep1Eta_SF_ttb = new TH1F("h_Lep1Eta_SF_ttb", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep1Eta_SF_ttbb = new TH1F("h_Lep1Eta_SF_ttbb", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep1Eta_SF_tt2b = new TH1F("h_Lep1Eta_SF_tt2b", "Lep0Eta",30, -2.4, 2.4);
        TH1F *h_Lep1Eta_SF_ttcc = new TH1F("h_Lep1Eta_SF_ttcc", "Lep0Eta",30, -2.4, 2.4);


        TH1F *h_Lep1Phi = new TH1F("h_Lep1Phi", "Lep1Phi",40, -4, 4);
        TH1F *h_Lep1Iso = new TH1F("h_Lep1Iso", "Lep1Iso",30, 0, 0.3);
	TH1F *h_Lep1IsoSF = new TH1F("h_Lep1IsoSF", "Lep1IsoSF",20, 0.9, 1.1);
	TH1F *h_Lep1IDSF = new TH1F("h_Lep1IDSF", "Lep1IDSF",20, 0.9, 1.1);



	TH1F *h_DeepCSVSF = new TH1F("h_DeepCSVSF","DeepCSVSF",40, 0 ,1);
	TH1F *h_DeepCSVSF_ttlf = new TH1F("h_DeepCSVSF_ttlf","DeepCSVSF",40, 0 ,1);
	TH1F *h_DeepCSVSF_ttb = new TH1F("h_DeepCSVSF_ttb","DeepCSVSF",40, 0 ,1);
	TH1F *h_DeepCSVSF_ttbb = new TH1F("h_DeepCSVSF_ttbb","DeepCSVSF",40, 0 ,1);
	TH1F *h_DeepCSVSF_tt2b = new TH1F("h_DeepCSVSF_tt2b","DeepCSVSF",40, 0 ,1);
	TH1F *h_DeepCSVSF_ttcc = new TH1F("h_DeepCSVSF_ttcc","DeepCSVSF",40, 0 ,1);




  	TH1F *h_NumberOfJets = new TH1F("h_NumberOfJets", "NumberOfJets",8, 2 , 10);
	TH1F *h_NumberOfJets_ee = new TH1F("h_NumberOfJets_ee", "NumberOfJets_ee",8, 2 , 10);	
	TH1F *h_NumberOfJets_emu = new TH1F("h_NumberOfJets_emu", "NumberOfJets_emu",8, 2 , 10);
	TH1F *h_NumberOfJets_mumu = new TH1F("h_NumberOfJets_mumu", "NumberOfJets_mumu",8, 2 , 10);
	TH1F *h_NumberOfJets_ttlf = new TH1F("h_NumberOfJets_ttlf", "NumberOfJets",8, 2, 10);
	TH1F *h_NumberOfJets_ttlf_ee = new TH1F("h_NumberOfJets_ttlf_ee", "NumberOfJets_ee",8, 2, 10);
	TH1F *h_NumberOfJets_ttlf_emu = new TH1F("h_NumberOfJets_ttlf_emu", "NumberOfJets_emu",8, 2, 10);
	TH1F *h_NumberOfJets_ttlf_mumu = new TH1F("h_NumberOfJets_ttlf_mumu", "NumberOfJets_mumu",8, 2, 10);
	TH1F *h_NumberOfJets_ttb = new TH1F("h_NumberOfJets_ttb", "NumberOfJets",8, 2, 10);
	TH1F *h_NumberOfJets_ttb_ee = new TH1F("h_NumberOfJets_ttb_ee", "NumberOfJets_ee",8, 2, 10);
	TH1F *h_NumberOfJets_ttb_emu = new TH1F("h_NumberOfJets_ttb_emu", "NumberOfJets_emu",8, 2, 10);
	TH1F *h_NumberOfJets_ttb_mumu = new TH1F("h_NumberOfJets_ttb_mumu", "NumberOfJets_mumu",8, 2, 10);
	TH1F *h_NumberOfJets_ttbb = new TH1F("h_NumberOfJets_ttbb", "NumberOfJets",8, 2, 10);
	TH1F *h_NumberOfJets_ttbb_ee = new TH1F("h_NumberOfJets_ttbb_ee", "NumberOfJets_ee",8, 2, 10);
	TH1F *h_NumberOfJets_ttbb_emu = new TH1F("h_NumberOfJets_ttbb_emu", "NumberOfJets_emu",8, 2, 10);
	TH1F *h_NumberOfJets_ttbb_mumu = new TH1F("h_NumberOfJets_ttbb_mumu", "NumberOfJets_mumu",8, 2, 10);
	TH1F *h_NumberOfJets_tt2b = new TH1F("h_NumberOfJets_tt2b", "NumberOfJets",8, 2, 10);
	TH1F *h_NumberOfJets_tt2b_ee = new TH1F("h_NumberOfJets_tt2b_ee", "NumberOfJets_ee",8, 2, 10);
	TH1F *h_NumberOfJets_tt2b_emu = new TH1F("h_NumberOfJets_tt2b_emu", "NumberOfJets_emu",8, 2, 10);
	TH1F *h_NumberOfJets_tt2b_mumu = new TH1F("h_NumberOfJets_tt2b_mumu", "NumberOfJets_mumu",8, 2, 10);
	TH1F *h_NumberOfJets_ttcc = new TH1F("h_NumberOfJets_ttcc", "NumberOfJets",8, 2, 10);
	TH1F *h_NumberOfJets_ttcc_ee = new TH1F("h_NumberOfJets_ttcc_ee", "NumberOfJets_ee",8, 2, 10);
	TH1F *h_NumberOfJets_ttcc_emu = new TH1F("h_NumberOfJets_ttcc_emu", "NumberOfJets_emu",8, 2, 10);
	TH1F *h_NumberOfJets_ttcc_mumu = new TH1F("h_NumberOfJets_ttcc_mumu", "NumberOfJets_mumu",8, 2, 10);

	TH1F *h_NumberOfJets_SF = new TH1F("h_NumberOfJets_SF", "NumberOfJets_SF",8, 2, 10);  
	TH1F *h_NumberOfJets_SF_ee = new TH1F("h_NumberOfJets_SF_ee", "NumberOfJets_SF_ee",8, 2, 10);
	TH1F *h_NumberOfJets_SF_emu = new TH1F("h_NumberOfJets_SF_emu", "NumberOfJets_SF_emu",8, 2, 10);
	TH1F *h_NumberOfJets_SF_mumu = new TH1F("h_NumberOfJets_SF_mumu", "NumberOfJets_SF_mumu",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttlf = new TH1F("h_NumberOfJets_SF_ttlf", "NumberOfJets",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttlf_ee = new TH1F("h_NumberOfJets_SF_ttlf_ee", "NumberOfJets_ee",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttlf_emu = new TH1F("h_NumberOfJets_SF_ttlf_emu", "NumberOfJets_emu",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttlf_mumu = new TH1F("h_NumberOfJets_SF_ttlf_mumu", "NumberOfJets_mumu",8, 2, 10);	



        TH1F *h_NumberOfJets_SF_ttb = new TH1F("h_NumberOfJets_SF_ttb", "NumberOfJets",8, 2, 10);
        TH1F *h_NumberOfJets_SF_ttb_ee = new TH1F("h_NumberOfJets_SF_ttb_ee", "NumberOfJets_ee",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttb_emu = new TH1F("h_NumberOfJets_SF_ttb_emu", "NumberOfJets_emu",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttb_mumu = new TH1F("h_NumberOfJets_SF_ttb_mumu", "NumberOfJets_mumu",8, 2, 10);
  	TH1F *h_NumberOfJets_SF_ttbb = new TH1F("h_NumberOfJets_SF_ttbb", "NumberOfJets",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttbb_ee = new TH1F("h_NumberOfJets_SF_ttbb_ee", "NumberOfJets_ee",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttbb_emu = new TH1F("h_NumberOfJets_SF_ttbb_emu", "NumberOfJets_emu",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttbb_mumu = new TH1F("h_NumberOfJets_SF_ttbb_mumu", "NumberOfJets_mumu",8, 2, 10);
        TH1F *h_NumberOfJets_SF_tt2b = new TH1F("h_NumberOfJets_SF_tt2b", "NumberOfJets",8, 2, 10);
	TH1F *h_NumberOfJets_SF_tt2b_ee = new TH1F("h_NumberOfJets_SF_tt2b_ee", "NumberOfJets_ee",8, 2, 10);
	TH1F *h_NumberOfJets_SF_tt2b_emu = new TH1F("h_NumberOfJets_SF_tt2b_emu", "NumberOfJets_emu",8, 2, 10);
	TH1F *h_NumberOfJets_SF_tt2b_mumu = new TH1F("h_NumberOfJets_SF_tt2b_mumu", "NumberOfJets_mumu",8, 2, 10);
        TH1F *h_NumberOfJets_SF_ttcc = new TH1F("h_NumberOfJets_SF_ttcc", "NumberOfJets",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttcc_ee = new TH1F("h_NumberOfJets_SF_ttcc_ee", "NumberOfJets_ee",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttcc_emu = new TH1F("h_NumberOfJets_SF_ttcc_emu", "NumberOfJets_emu",8, 2, 10);
	TH1F *h_NumberOfJets_SF_ttcc_mumu = new TH1F("h_NumberOfJets_SF_ttcc_mumu", "NumberOfJets_mumu",8, 2, 10);

	TH1F *h_NumberOfbTags = new TH1F("h_NumberOfbTags", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_ee = new TH1F("h_NumberOfbTags_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_emu = new TH1F("h_NumberOfbTags_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_mumu = new TH1F("h_NumberOfbTags_mumu", "NumberOfbTags_mumu",5, 1,6);
	TH1F *h_NumberOfbTags_ttlf = new TH1F("h_NumberOfbTags_ttlf", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_ttlf_ee = new TH1F("h_NumberOfbTags_ttlf_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_ttlf_emu = new TH1F("h_NumberOfbTags_ttlf_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_ttlf_mumu = new TH1F("h_NumberOfbTags_ttlf_mumu", "NumberOfbTags_mumu",5, 1,6);
	TH1F *h_NumberOfbTags_ttb = new TH1F("h_NumberOfbTags_ttb", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_ttb_ee = new TH1F("h_NumberOfbTags_ttb_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_ttb_emu = new TH1F("h_NumberOfbTags_ttb_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_ttb_mumu = new TH1F("h_NumberOfbTags_ttb_mumu", "NumberOfbTags_mumu",5, 1,6);
	TH1F *h_NumberOfbTags_ttbb = new TH1F("h_NumberOfbTags_ttbb", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_ttbb_ee = new TH1F("h_NumberOfbTags_ttbb_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_ttbb_emu = new TH1F("h_NumberOfbTags_ttbb_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_ttbb_mumu = new TH1F("h_NumberOfbTags_ttbb_mumu", "NumberOfbTags_mumu",5, 1,6);
	TH1F *h_NumberOfbTags_tt2b = new TH1F("h_NumberOfbTags_tt2b", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_tt2b_ee = new TH1F("h_NumberOfbTags_tt2b_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_tt2b_emu = new TH1F("h_NumberOfbTags_tt2b_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_tt2b_mumu = new TH1F("h_NumberOfbTags_tt2b_mumu", "NumberOfbTags_mumu",5, 1,6);
	TH1F *h_NumberOfbTags_ttcc = new TH1F("h_NumberOfbTags_ttcc", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_ttcc_ee = new TH1F("h_NumberOfbTags_ttcc_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_ttcc_emu = new TH1F("h_NumberOfbTags_ttcc_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_ttcc_mumu = new TH1F("h_NumberOfbTags_ttcc_mumu", "NumberOfbTags_mumu",5, 1,6);


	TH1F *h_NumberOfbTags_SF = new TH1F("h_NumberOfbTags_SF", "NumberOfbTags_SF",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ee = new TH1F("h_NumberOfbTags_SF_ee", "NumberOfbTags_SF_ee",5, 1,6);
	TH1F *h_NumberOfbTags_SF_emu = new TH1F("h_NumberOfbTags_SF_emu", "NumberOfbTags_SF_emu",5, 1,6);
	TH1F *h_NumberOfbTags_SF_mumu = new TH1F("h_NumberOfbTags_SF_mumu", "NumberOfbTags_SF_mumu",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttlf = new TH1F("h_NumberOfbTags_SF_ttlf", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttlf_ee = new TH1F("h_NumberOfbTags_SF_ttlf_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttlf_emu = new TH1F("h_NumberOfbTags_SF_ttlf_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttlf_mumu = new TH1F("h_NumberOfbTags_SF_ttlf_mumu", "NumberOfbTags_mumu",5, 1,6);
        TH1F *h_NumberOfbTags_SF_ttb = new TH1F("h_NumberOfbTags_SF_ttb", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttb_ee = new TH1F("h_NumberOfbTags_SF_ttb_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttb_emu = new TH1F("h_NumberOfbTags_SF_ttb_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttb_mumu = new TH1F("h_NumberOfbTags_SF_ttb_mumu", "NumberOfbTags_mumu",5, 1,6);
        TH1F *h_NumberOfbTags_SF_ttbb = new TH1F("h_NumberOfbTags_SF_ttbb", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttbb_ee = new TH1F("h_NumberOfbTags_SF_ttbb_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttbb_emu = new TH1F("h_NumberOfbTags_SF_ttbb_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttbb_mumu = new TH1F("h_NumberOfbTags_SF_ttbb_mumu", "NumberOfbTags_mumu",5, 1,6);
        TH1F *h_NumberOfbTags_SF_tt2b = new TH1F("h_NumberOfbTags_SF_tt2b", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_SF_tt2b_ee = new TH1F("h_NumberOfbTags_SF_tt2b_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_SF_tt2b_emu = new TH1F("h_NumberOfbTags_SF_tt2b_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_SF_tt2b_mumu = new TH1F("h_NumberOfbTags_SF_tt2b_mumu", "NumberOfbTags_mumu",5, 1,6);
        TH1F *h_NumberOfbTags_SF_ttcc = new TH1F("h_NumberOfbTags_SF_ttcc", "NumberOfbTags",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttcc_ee = new TH1F("h_NumberOfbTags_SF_ttcc_ee", "NumberOfbTags_ee",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttcc_emu = new TH1F("h_NumberOfbTags_SF_ttcc_emu", "NumberOfbTags_emu",5, 1,6);
	TH1F *h_NumberOfbTags_SF_ttcc_mumu = new TH1F("h_NumberOfbTags_SF_ttcc_mumu", "NumberOfbTags_mumu",5, 1,6);
  	TH1F *h_NumberOfLeptons = new TH1F("h_NumberOfLeptons", "NumberOfLeptons",15, 0, 15);
  	TH1F *h_NumberOfElectrons = new TH1F("h_NumberOfElectrons", "NumberOfElectrons",15, 0, 15);
    	TH1F *h_NumberOfMuons = new TH1F("h_NumberOfMuons", "NumberOfMuons",15, 0, 15);

	TH1F *h_additionalJetEventId = new TH1F("h_additionalJetEventId","AdditionalJetEventId",21000,0,21000);
	TH1F *h_cutFlow = new TH1F("h_cutFlow","CutFlow",10,0,10);

  	THStack *hs_JetAllPT = new THStack("hs_JetAllPT","");
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
  // *** 2. Set tree structure and variables to read
  	yggdrasilEventVars *eve=0;


	int xbin_EleIsoSF = -1;
        int ybin_EleIsoSF = -1;
        double eleIDSF_scEta = -99;
        double eleIDSF_Pt = -99;

        int xbin_EleIDSF = -1;
        int ybin_EleIDSF = -1;
        double eleIsoSF_scEta = -99;
       	double eleIsoSF_Pt = -99;

	int xbin_MuonIsoSF = -1;
        int ybin_MuonIsoSF = -1;
        double muonIDSF_absEta = -99;
        double muonIDSF_Pt = -99;

        int xbin_MuonIDSF = -1;
        int ybin_MuonIDSF = -1;
        double muonIsoSF_absEta = -99;
        double muonIsoSF_Pt = -99;

	double wgt_generator_ = 0.0;
	double mcWeight_value = 0.0;



	double evtSF = 1.0;
	int bTags = 0;
  	double SumJetPT = 0.0;
  	double SumJetEta = 0.0;
  	double MET = 0.0;
  	double METphi = 0.0;
	int NumberOfPV = 0;


	double PotJetPT = 0;
	int JetNumber = 0;
  	double JetPT=0.0;
  	double JetEta=0.0;
	double JetCSV_b = 0.0;
	double JetCSV_bb = 0.0;
	double JetCSVsf = 0.0;
	double CSVsf = 0.0;

  	double MuonPT=0.0;
  	double MuonEta=0.0;
  	double MuonPhi=0.0;
	double MuonIso=0;
	double MuonIsoSF = 0;
	double MuonIDSF = 0;


  	double ElePT=0.0;
  	double EleEta=0.0;
  	double ElePhi=0.0; 
	double EleIso = 0;
	double EleIDSF = 0;
	double EleIsoSF = 0;

  	int NumberOfJets = 0;
	int NumberOfLooseJets = 0;
	int NumberOfBadJets = 0;
  	int NumberOfLeptons = 0;
	int NumberOfLooseLeptons = 0;
	int NumberOfBadLeptons = 0;
  	int NumberOfElectrons=0;
	int NumberOfLooseElectrons = 0;
	int NumberOfBadElectrons = 0;
  	int NumberOfMuons=0;
	int NumberOfGoodMuons = 0;
	int NumberOfLooseMuons = 0;
	int NumberOfBadMuons = 0;

	int passedSLeTrigOnly = 0;
	int passedSLeTrig = 0;
	//SL_Electron
	//	B,C,D,E,F
  	int passedHLT_Ele35_WPTight_Gsf_v_=0;
	int passedHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_ = 0;

	int passedSLmuTrig = 0;
	//SL:Muon
	//	B,C,D
	int passedHLT_IsoMu24_eta2p1_v_ = 0;
  	int passedHLT_IsoMu27_v_=0;

	int passedDLeeTrig = 0;
  	//DL:ee
  	//	B,C,D,E,F
	int passedHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = 0;
	int passedHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0;

	int passedDLemuTrig = 0;
	//DL:emu
	//	B,C,D,E,F
	int passedHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = 0;
	int passedHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0;
	int passedHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0;
	int passedHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0;


	int passedDLmumuTrig = 0;
	//DL:mumu
	//	B
	int passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_ = 0;
	//DL:mumu
	//	C,D,E,F
	int passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_ = 0;
	// 	C,D,E,F
	int passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_ = 0;






	double EtaJets = 0.0;  

  	chain->SetBranchAddress("eve.", &eve );

  	long total_entries = chain->GetEntries();
  	cout << "Chain entries: " << total_entries << endl;
	if(n_entries < 0 )n_entries = total_entries;
	
	cout << "Will be running over " << n_entries << " entries..." << endl;

	int NumberOfGoodJets = 0;
	//int NumberOfLooseJets = 0;
	//int NumberOfBadJets = 0;
	int NumberOfGoodElectrons = 0;
	int NumberOfGoodLeptons = 0;
	//int NumberOfBadElectrons = 0;
	//int NumberOfLooseElectrons = 0;
	//int NumberOfGoodMuons = 0;
	//int NumberOfBadMuons = 0;	
	//int NumberOfLooseMuons = 0;

	//BDT 2/21
        std::vector<TLorentzVector> selectedLeptonP4;
        std::vector<double> selectedLeptonCharge;
        std::vector<TLorentzVector> selectedJetP4;
        std::vector<double> selectedJetCSV;
        TLorentzVector metP4;

	TLorentzVector Lepton1;
	TLorentzVector Lepton2;
	TLorentzVector Lepton3;
	TLorentzVector Lepton4;
	TLorentzVector DiLepton;

	TLorentzVector JetLorentzVector;
        double bdt_score;

	TLorentzVector Jet1;
        TLorentzVector Jet2;
        TLorentzVector Jet3;
        TLorentzVector Jet4;
        TLorentzVector Jet5;
        TLorentzVector Jet6;
        TLorentzVector Jet7;
        TLorentzVector Jet8;
        TLorentzVector Jet9;
        TLorentzVector Jet10;
        TLorentzVector Jet11;
        TLorentzVector Jet12;




	std::vector<int> WhichGoodJets;
	std::vector<double> PotJetsPT;
	std::vector<double> GoodJetsE;
	std::vector<double> GoodJetsPT;
	std::vector<double> GoodJetsEta;
	std::vector<double> GoodJetsPhi;
        std::vector<double> GoodJetsM;	
	std::vector<double> GoodJetsCSVv2;
	std::vector<double> GoodJetsDeepCSV_b;
	std::vector<double> GoodJetsSFDeepCSV;
	std::vector<double> GoodJetsDeepCSV_bb;
	std::vector<double> GoodJetsCSVSF;
	std::vector<double> PotJetsJERsf;
	std::vector<double> GoodJetsJERsf;

	std::vector<double> WhichBadJets;
	std::vector<double> BadJetsE;
        std::vector<double> BadJetsPT;
        std::vector<double> BadJetsEta;
        std::vector<double> BadJetsPhi;
        std::vector<double> BadJetsDeepCSV_b;

        std::vector<double> WhichLooseJets;
	std::vector<double> LooseJetsE;
        std::vector<double> LooseJetsPT;
        std::vector<double> LooseJetsEta;
        std::vector<double> LooseJetsPhi;
        std::vector<double> LooseJetsDeepCSV_b;
	std::vector<double> LooseJetsDeepCSV_bb;


	std::vector<int> WhichGoodLeptons;
	std::vector<int> GoodLeptonsCharge;
	std::vector<double> GoodLeptonsE;
	std::vector<double> GoodLeptonsPT;
        std::vector<double> GoodLeptonsEta;
        std::vector<double> GoodLeptonsPhi;	
	std::vector<double> GoodLeptonsIso;

        std::vector<int> WhichBadLeptons;
        std::vector<int> BadLeptonsCharge;
	std::vector<double> BadLeptonsE;
        std::vector<double> BadLeptonsPT;
        std::vector<double> BadLeptonsEta;
        std::vector<double> BadLeptonsPhi;
        std::vector<double> BadLeptonsIso;



        std::vector<int> WhichLooseLeptons;
	std::vector<int> LooseLeptonsCharge;
	std::vector<double> LooseLeptonsE;
        std::vector<double> LooseLeptonsPT;
        std::vector<double> LooseLeptonsEta;
        std::vector<double> LooseLeptonsPhi;
        std::vector<double> LooseLeptonsIso;



	std::vector<int> WhichGoodElectrons;
	std::vector<int> GoodElectronsCharge;
	std::vector<double> GoodElectronsE;
	std::vector<double> GoodElectronsPT;
        std::vector<double> GoodElectronsEta;
	std::vector<double> GoodElectronsPhi;	
	std::vector<double> GoodElectronsIso;
	std::vector<double> GoodElectronsIsoSF;
	std::vector<double> GoodElectronsIDSF;

	std::vector<int> WhichBadElectrons;
        std::vector<int> BadElectronsCharge;
	std::vector<double> BadElectronsE;
        std::vector<double> BadElectronsPT;
        std::vector<double> BadElectronsEta;
        std::vector<double> BadElectronsPhi;
        std::vector<double> BadElectronsIso;	



        std::vector<int> WhichLooseElectrons;
	std::vector<int> LooseElectronsCharge;
	std::vector<double> LooseElectronsE;
        std::vector<double> LooseElectronsPT;
        std::vector<double> LooseElectronsEta;
        std::vector<double> LooseElectronsPhi;
        std::vector<double> LooseElectronsIso;	
        std::vector<double> LooseElectronsIsoSF;
        std::vector<double> LooseElectronsIDSF;


	std::vector<int> WhichGoodMuons;	
	std::vector<int> GoodMuonsCharge;
	std::vector<double> GoodMuonsE;
	std::vector<double> GoodMuonsPT;
        std::vector<double> GoodMuonsEta;
	std::vector<double> GoodMuonsPhi;
	std::vector<double> GoodMuonsIso;
	std::vector<double> GoodMuonsIsoSF;
	std::vector<double> GoodMuonsIDSF;


	std::vector<int> WhichBadMuons;
        std::vector<int> BadMuonsCharge;
        std::vector<double> BadMuonsE;
	std::vector<double> BadMuonsPT;
        std::vector<double> BadMuonsEta;
        std::vector<double> BadMuonsPhi;
        std::vector<double> BadMuonsIso;


        std::vector<int> WhichLooseMuons;
	std::vector<int> LooseMuonsCharge;
   	std::vector<double> LooseMuonsE;
     	std::vector<double> LooseMuonsPT;
        std::vector<double> LooseMuonsEta;
        std::vector<double> LooseMuonsPhi;
        std::vector<double> LooseMuonsIso;
        std::vector<double> LooseMuonsIsoSF;
        std::vector<double> LooseMuonsIDSF;

	std::vector<int> indexLeptons;
	std::vector<int> indexBadLeptons;

	std::vector<int> indexLooseLeptons;
	std::vector<int> indexElectrons;
	std::vector<int> indexLooseElectrons;
	std::vector<int> indexBadElectrons;
	std::vector<int> indexMuons;
	std::vector<int> indexBadMuons;
	std::vector<int> indexLooseMuons;
	std::vector<int> indexJets;
	std::vector<int> indexBadJets;
	std::vector<int> indexLooseJets;



	int Lep0Charge = 0;
	int Lep1Charge = 0;

	double Lep0E = 0;
	double Lep0PT = 0;
	double Lep0Eta = 0;
	double Lep0Phi = 0;
	double Lep0Iso = 0;
	double Lep0IsoSF = 0;
	double Lep0IDSF = 0;

	double Lep1E = 0;
	double Lep1PT = 0;
	double Lep1Eta = 0;
	double Lep1Phi = 0;
	double Lep1Iso = 0;
	double Lep1IsoSF = 0;
	double Lep1IDSF = 0;

	double lep1_isoSF = -99;

	double csvSF = 0;
	double mll = 0;

	int is_e = 0;
	int is_mu = 0;
	int is_ee = 0;
	int is_emu = 0;
	int is_mumu = 0;
	
	int is_ttlf = 0;
	int is_ttb = 0;
	int is_ttbb = 0;
	int is_tt2b = 0;
	int is_ttcc = 0;
	int additionalJetEventId = -99;



	int cfTotal = 0;
	int cfTrig = 0;
	int cfTrigLep = 0;
	int cfTrigLepJet = 0;
	int cfTrigLepJetBtag = 0;






	double deltaEta = 0;
	double deltaPhi = 0;

	double R_jetVeto = 0;
	int flagJetDR = 0;

	double sfJet0 = 0.0;
	double sfJet1 = 0.0;
	double sfJet2 = 0.0;
	double sfJet3 = 0.0;

	int flagEventEle = 0;
	int flagEventMuon = 0;


	int failedEvent = 0;
/*
        ofstream fCSV;
	//ofstream failedCSV;
	TString outputCSV;
	//TString failedOutputCSV;
	//failedOutputCSV = "failedEvents_EventReadOut_v22_Dec6"+outSingleFile+".csv";
	//if(applyJERsf ==1 && applyIndCSVsf == 1)outputCSV = "syncExercise_"+output+"_EventReadOut_v15_BensRecipe_Oct17_v14.csv";
	//else if (applyJERsf == 1 && applyIndCSVsf == 0)outputCSV = "syncExercise_"+output+"_EventReadOut_v15_BensRecipe_Oct17.csv";
	//else if (applyJERsf == 0 && applyIndCSVsf == 1)outputCSV = "syncExercise_"+output+"_EventReadOut_v15_BensRecipe_Oct17_v14.csv";
	outputCSV = "syncExercise_"+output+"_EventReadOut_BDT_METFilter_Feb22_PT15_2j1b_"+outSingleFile+".csv";
        fCSV.open(outputCSV);
	//failedCSV.open(failedOutputCSV);
        fCSV<<"run,lumi,event,is_e,is_mu,is_ee,is_emu,is_mumu,n_jets,n_btags,";
        fCSV<<"lep1_pt,lep1_eta,lep1_iso,lep1_pdgId,lep1_idSF,lep1_isoSF,";
        fCSV<<"lep2_pt,lep2_eta,lep2_iso,lep2_pdgId,lep2_idSF,lep2_isoSF,";
        fCSV<<"jet1_pt,jet1_eta,jet1_phi,jet1_jesSF,jet1_jesSF_up,jet1_jesSF_down,jet1_jesSF_PileUpDataMC_down,jet1_jesSF_RelativeFSR_up,jet1_jerSF_nominal,jet1_csv,jet1_PUJetId,jet1_PUJetDiscriminant,";
        fCSV<<"jet2_pt,jet2_eta,jet2_phi,jet2_jesSF,jet2_jesSF_up,jet2_jesSF_down,jet2_jesSF_PileUpDataMC_down,jet2_jesSF_RelativeFSR_up,jet2_jerSF_nominal,jet2_csv,jet2_PUJetId,jet2_PUJetDiscriminant,";
        fCSV<<"MET_pt,MET_phi,mll,ttHFCategory,ttHFGenFilterTag,n_interactions,puWeight,csvSF,csvSF_lf_up,csvSF_hf_down,csvSF_cErr1_down"<<endl;
*/

	//failedCSV<<"run,lumi,event,is_e,is_mu,is_ee,is_emu,is_mumu,n_jets,n_btags,";
        //failedCSV<<"lep1_pt,lep1_eta,lep1_iso,lep1_pdgId,lep1_idSF,lep1_isoSF,";
        //failedCSV<<"lep2_pt,lep2_eta,lep2_iso,lep2_pdgId,lep2_idSF,lep2_isoSF,";
        //failedCSV<<"jet1_pt,jet1_eta,jet1_phi,jet1_jesSF,jet1_jesSF_up,jet1_jesSF_down,jet1_jesSF_PileUpDataMC_down,jet1_jesSF_RelativeFSR_up,jet1_jerSF_nominal,jet1_csv,jet1_PUJetId,jet1_PUJetDiscriminant,";
        //failedCSV<<"jet2_pt,jet2_eta,jet2_phi,jet2_jesSF,jet2_jesSF_up,jet2_jesSF_down,jet2_jesSF_PileUpDataMC_down,jet2_jesSF_RelativeFSR_up,jet2_jerSF_nominal,jet2_csv,jet2_PUJetId,jet2_PUJetDiscriminant,";
        //failedCSV<<"MET_pt,MET_phi,mll,ttHFCategory,ttHFGenFilterTag,n_interactions,puWeight,csvSF,csvSF_lf_up,csvSF_hf_down,csvSF_cErr1_down"<<endl;


//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
// *** 3. Start Running Over Events
  	for(int i = 0; i < n_entries; i++) {
    		
		bdt_score = 0;
		evtSF = 1.0;
		wgt_generator_ = 0.0;
		mcWeight_value = 0.0;		



    		chain->GetEntry(i);
	    
		if(is_MC == 1){
			wgt_generator_ = eve->wgt_generator_;
			//mcWeight_value = eve->mcWeight_value;
		}

		cfTotal = cfTotal + 1;

		WhichGoodJets.clear();
		PotJetsPT.clear();
        	GoodJetsPT.clear();
        	GoodJetsEta.clear();
        	GoodJetsCSVv2.clear();
        	GoodJetsDeepCSV_b.clear();
        	GoodJetsSFDeepCSV.clear();
		GoodJetsDeepCSV_bb.clear();
		GoodJetsCSVSF.clear();
		PotJetsJERsf.clear();
		GoodJetsJERsf.clear();

		WhichGoodMuons.clear();
		GoodMuonsCharge.clear();
	        GoodMuonsE.clear();
		GoodMuonsPT.clear();
                GoodMuonsEta.clear();
                GoodMuonsPhi.clear();
                GoodMuonsIso.clear();
		GoodMuonsIsoSF.clear();
		GoodMuonsIDSF.clear();


		WhichLooseMuons.clear();
		LooseMuonsCharge.clear(); 
       		LooseMuonsE.clear();
		LooseMuonsPT.clear(); 
	      	LooseMuonsEta.clear();
		LooseMuonsPhi.clear();
		LooseMuonsIso.clear();
		LooseMuonsIsoSF.clear();
		LooseMuonsIDSF.clear();


		WhichBadMuons.clear();


		WhichGoodElectrons.clear();
		GoodElectronsE.clear();
		GoodElectronsPT.clear();
                GoodElectronsEta.clear();
                GoodElectronsPhi.clear();
                GoodElectronsIso.clear();
		GoodElectronsIsoSF.clear();	        
		GoodElectronsIDSF.clear();

		WhichLooseElectrons.clear();
		LooseElectronsCharge.clear();
                LooseElectronsE.clear();
                LooseElectronsPT.clear();
                LooseElectronsEta.clear();
                LooseElectronsPhi.clear();
                LooseElectronsIso.clear();
		LooseElectronsIsoSF.clear();
		LooseElectronsIDSF.clear();

		WhichBadElectrons.clear();
                BadElectronsCharge.clear();
                BadElectronsE.clear();
                BadElectronsPT.clear();
                BadElectronsEta.clear();
                BadElectronsPhi.clear();
                BadElectronsIso.clear();



		WhichGoodElectrons.clear();
		WhichGoodMuons.clear();


//0000000000000000000000000000000000000000
		if(RunSingleEvent == 1 && eve->evt_ != EventNumber)continue;
		//EventNumber = eve->evt_;
//000000000000000000000000000000000000000


		if(DEBUG==1)EventNumber = eve->evt_;
		if(eve->evt_ == EventNumber){// || eve->evt_%1000 == 0){
			//failedEvent = 1;
			
 			/*cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
			for(int u = 0; u<99;u++){
				cout << "JetPT [0][X]" << endl;
				cout<< u<< " ... " << eve->jet_pt_[0][u] << endl;
				cout << "JetPT [1][X]" << endl;
				cout << u << " ___ " << eve->jet_pt_[1][u] << endl;
 				cout << "JetPT [2][X]" << endl;
                                cout << u << " ___ " << eve->jet_pt_[2][u] << endl;
			}
			*/
			cout<<"========================================================"<<endl;
			cout <<"   EVENT DETAILS BEFORE ANY CORRECTIONS"<<endl;
			cout<<"========================================================"<<endl;
			cout<< "EventNumber: " <<eve->evt_<<endl;
			cout << "Event Run Number: "<< eve->run_ << endl;
			cout << "Event Lumi Section: " << eve->lumi_ << endl;
			if(is_MC == 1){
				cout<<"eve->wgt_generator_: "<<eve->wgt_generator_<<endl;
				//cout<<"eve->mcWeight_value: "<<eve->mcWeight_value << endl;
			}
			cout<<"eve->passHLT_Ele35_WPTight_Gsf_v_: "<<eve->passHLT_Ele35_WPTight_Gsf_v_<<endl;
                        cout<< "eve->passHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_: "<<eve->passHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_<<endl;
                        cout<<"eve->passHLT_IsoMu24_eta2p1_v_: "<<eve->passHLT_IsoMu24_eta2p1_v_<<endl;
                        cout<<"eve->passHLT_IsoMu27_v_: "<<eve->passHLT_IsoMu27_v_<<endl;
                        cout<<"eve->passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_: "<<eve->passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_<<endl;
                        cout<<"eve->passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_: "<<eve->passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_<<endl;
                        cout<<"eve->passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_: "<<eve->passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_<<endl;
                        cout<<"eve->passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_: "<<eve->passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_<<endl;
                        cout<<"eve->passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_: "<<eve->passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_<<endl;
                        cout<<"eve->passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_: "<<eve->passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_<<endl;
                        cout<<"eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_: "<<eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_<<endl;
                        cout<<"eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_: "<<eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_<<endl;
                        cout<<"eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_: "<<eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_<<endl;
                        cout<<"===================================================================="<<endl;
                        cout<<" MET: " <<eve->MET_[0]<<endl;
			cout<<" MET_Type1xy: " <<eve->MET_Type1xy_[0]<<endl;
			cout<<" Number Of PV: "<< eve->numPVs_<<endl;
                        cout<<"==================================================================="<<endl;
		/*	cout<<"eve->passMETFilter_Flag_goodVertices_v_: "<<eve->passMETFilter_Flag_goodVertices_v_<<endl;
			cout<<"eve->passMETFilter_Flag_globalTightHalo2016Filter_v_: "<< eve->passMETFilter_Flag_globalTightHalo2016Filter_v_ << endl;
			cout<<"eve->passMETFilter_Flag_globalSuperTightHalo2016Filter_v_: "<< eve->passMETFilter_Flag_globalSuperTightHalo2016Filter_v_ << endl;
			cout<<"eve->passMETFilter_Flag_HBHENoiseFilter_v_: "<< eve->passMETFilter_Flag_HBHENoiseFilter_v_ << endl;
			cout<<"eve->passMETFilter_Flag_HBHENoiseIsoFilter_v_: "<< eve->passMETFilter_Flag_HBHENoiseIsoFilter_v_ << endl;	
			cout<<"eve->passMETFilter_Flag_EcalDeadCellTriggerPrimitiveFilter_v_: "<< eve->passMETFilter_Flag_EcalDeadCellTriggerPrimitiveFilter_v_ << endl;
			cout<<"eve->passMETFilter_Flag_BadPFMuonFilter_v_: "<< eve->passMETFilter_Flag_BadPFMuonFilter_v_ << endl;	
			cout<<"eve->passMETFilter_Flag_BadChargedCandidateFilter_v_: "<< eve->passMETFilter_Flag_BadChargedCandidateFilter_v_ << endl;
			cout<<"eve->passMETFilter_Flag_ecalBadCalibFilter_v_: "<< eve->passMETFilter_Flag_ecalBadCalibFilter_v_ << endl;
			cout<<"eve->passMETFilter_Flag_eeBadScFilter_v_: "<< eve->passMETFilter_Flag_eeBadScFilter_v_ << endl;
		*/	cout<<"==================================================================="<<endl;
                        cout<< "NumberOfJets: " <<((eve->jet_pt_[0]).size()) <<endl;
                        NumberOfJets = ((eve->jet_pt_[0]).size()) + 5  ;
                        cout<<"---------------------------"<<endl;
                        cout<< "jet0_pt: "<< eve->jet_pt_[0][0] << endl;
                        cout << "jet0_eta: "<< eve->jet_eta_[0][0] << endl;
        		cout << "jet0_phi: "<<eve->jet_phi_[0][0] << endl; 
	               cout << "jet0_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][0] << endl;
                	cout << "jet0_DeepCSV_bb: " << eve->jet_DeepCSV_bb_[0][0] << endl;
			cout << "jet0_puid: "<< eve->jet_puid_[0][0] << endl; 
		       cout<<"---------------------------"<<endl;
			if(NumberOfJets > 1){
				cout<< "jet1_pt: " <<eve->jet_pt_[0][1] << endl;
                        	cout << "jet1_eta: "<< eve->jet_eta_[0][1] << endl;
                        	cout << "jet1_phi: "<< eve->jet_phi_[0][1] << endl;
				cout << "jet1_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][1] << endl;
				cout << "jet1_DeepCSV_bb: " << eve->jet_DeepCSV_bb_[0][1] << endl;
				cout << "jet1_puid: "<< eve->jet_puid_[0][1] << endl;
			}
			cout<<"---------------------------"<<endl;
			if(NumberOfJets > 2){
                                cout<< "jet2_pt: "<<eve->jet_pt_[0][2] << endl;
                                cout << "jet2_eta: "<<eve->jet_eta_[0][2] << endl;
				cout << "jet2_phi: "<<eve->jet_phi_[0][2] << endl;
                                cout << "jet2_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][2] << endl;
                        	cout << "jet2_DeepCSV_bb: " << eve->jet_DeepCSV_bb_[0][2] << endl;
				cout << "jet2_puid: "<< eve->jet_puid_[0][2] << endl;
			}
                        cout<<"---------------------------"<<endl;
                        if(NumberOfJets > 3){
                                cout<< "jet3_pt: "<<eve->jet_pt_[0][3] << endl;
                                cout << "jet3_eta: " << eve->jet_eta_[0][3] << endl;
                                cout << "jet3_phi: " << eve->jet_phi_[0][3] << endl;
                                cout << "jet3_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][3] << endl;
                        	cout << "jet3_DeepCSV_bb: " << eve->jet_DeepCSV_bb_[0][3] << endl;
				cout << "jet3_puid: "<< eve->jet_puid_[0][4] << endl;
			}
                        cout<<"---------------------------"<<endl;
                        if(NumberOfJets > 4){
                                cout<< "jet4_pt: "<<eve->jet_pt_[0][4] << endl;
                                cout << "jet4_eta: " << eve->jet_eta_[0][4] << endl;
                                cout << "jet4_phi: " << eve->jet_phi_[0][4] << endl;
                                cout << "jet4_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][4] << endl;
                        	cout << "jet4_puid: "<< eve->jet_puid_[0][4] << endl;
			}
                        cout<<"---------------------------"<<endl;
                        if(NumberOfJets > 5){
                                cout<< "jet5_pt: "<<eve->jet_pt_[0][5] << endl;
                                cout << "jet5_eta: " << eve->jet_eta_[0][5] << endl;
                                cout << "jet5_phi: " << eve->jet_phi_[0][5] << endl;
                                cout << "jet5_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][5] << endl;
                        	cout << "jet5_puid: "<< eve->jet_puid_[0][5] << endl;
			}
                        cout<<"---------------------------"<<endl;
                        if(NumberOfJets > 6){
                                cout<< "jet6_pt: "<<eve->jet_pt_[0][6] << endl;
                                cout << "jet6_eta: " << eve->jet_eta_[0][6] << endl;
                                cout << "jet6_phi: " << eve->jet_phi_[0][6] << endl;
                                cout << "jet6_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][6] << endl;
                        	cout << "jet6_puid: "<< eve->jet_puid_[0][6] << endl;
			}
                        cout<<"---------------------------"<<endl;
                        if(NumberOfJets > 7){
                                cout<< "jet7_pt: "<<eve->jet_pt_[0][7] << endl;
                                cout << "jet7_eta: " << eve->jet_eta_[0][7] << endl;
                                cout << "jet7_phi: " << eve->jet_phi_[0][7] << endl;
                                cout << "jet7_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][7] << endl;
                        	cout << "jet7_puid: "<< eve->jet_puid_[0][7] << endl;
			}
                        if(NumberOfJets > 8){
                                cout<< "jet8_pt: "<<eve->jet_pt_[0][8] << endl;
                                cout << "jet8_eta: " << eve->jet_eta_[0][8] << endl;
                                cout << "jet8_phi: " << eve->jet_phi_[0][8] << endl;
                                cout << "jet8_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][8] << endl;
                        	cout << "jet8_puid: "<< eve->jet_puid_[0][8] << endl;
			}
                        if(NumberOfJets > 9){
                                cout<< "jet9_pt: "<<eve->jet_pt_[0][9] << endl;
                                cout << "jet9_eta: " << eve->jet_eta_[0][9] << endl;
                                cout << "jet9_phi: " << eve->jet_phi_[0][9] << endl;
                                cout << "jet9_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][9] << endl;
                        	cout << "jet9_puid: "<< eve->jet_puid_[0][9] << endl;
			}
                        if(NumberOfJets > 10){
                                cout<< "jet10_pt: "<<eve->jet_pt_[0][10] << endl;
                                cout << "jet10_eta: " << eve->jet_eta_[0][10] << endl;
                                cout << "jet10_phi: " << eve->jet_phi_[0][10] << endl;
                                cout << "jet10_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][10] << endl;
                        }
                        if(NumberOfJets > 11){
                                cout<< "jet11_pt: "<<eve->jet_pt_[0][11] << endl;
                                cout << "jet11_eta: " << eve->jet_eta_[0][11] << endl;
                                cout << "jet11_phi: " << eve->jet_phi_[0][11] << endl;
                                cout << "jet11_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][11] << endl;
                        }
                        if(NumberOfJets > 12){
                                cout<< "jet12_pt: "<<eve->jet_pt_[0][12] << endl;
                                cout << "jet12_eta: " << eve->jet_eta_[0][12] << endl;
                                cout << "jet12_phi: " << eve->jet_phi_[0][12] << endl;
                                cout << "jet12_DeepCSV_b: " << eve->jet_DeepCSV_b_[0][12] << endl;
                        }
                        cout<<"======================================================="<<endl;
                        cout<< "NumberOfLeptons: " << (eve->lepton_pt_).size() << endl;
			if((eve->lepton_pt_).size() > 0){
				cout << "lep0_pt: " << eve->lepton_pt_[0] << endl;
				cout << "lep0_pt_1: "<< eve->lepton_pt_[1] << endl;
				cout << "lep0_pt_2: "<< eve->lepton_pt_[2] << endl;
                        	cout<< "lep0_eta: " << eve->lepton_eta_[0] << endl;
                        	cout<< "lep0_phi: " << eve->lepton_phi_[0] << endl;
			        cout<<"lep0_scEta: "<< eve->lepton_scEta_[0]<<endl;
				cout<<"lep0_relIso: "<< eve->lepton_relIso_[0]<<endl;
				cout<<"THIS IS WHERE IT BREAKS?!?!?!!?!?!?"<<endl;
				cout<<"lepton_isTight_: "<< eve->lepton_isTight_[0] <<endl;
				 cout<<"========================================================"<<endl;
                        }
			if((eve->lepton_pt_).size() > 1){
                                cout<< "lep1_pt: " << eve->lepton_pt_[1] << endl;
                                cout<< "lep1_eta: " << eve->lepton_eta_[1] << endl;
                                cout<< "lep1_phi: " << eve->lepton_phi_[1] << endl;
			        cout<<"lep1_scEta: "<<eve->lepton_scEta_[1]<<endl;
                                cout<<"lep1_relIso: "<< eve->lepton_relIso_[1]<<endl;
				cout<<"lepton_isTight_: "<<eve->lepton_isTight_[1]<<endl;
				cout<<"========================================================"<<endl;
                        }
                        if((eve->lepton_pt_).size() > 2){
                                cout<< "lep2_pt: " << eve->lepton_pt_[2] << endl;
                                cout<< "lep2_eta: " << eve->lepton_eta_[2] << endl;
                                cout<< "lep2_phi: " << eve->lepton_phi_[2] << endl;
                                cout<<"lep2_scEta: "<<eve->lepton_scEta_[2]<<endl;
                                cout<<"lep2_relIso: "<< eve->lepton_relIso_[2]<<endl;
				cout<<"lepton_isTight_: "<<eve->lepton_isTight_[2]<<endl;
				cout<<"========================================================"<<endl;
			}
                        if((eve->lepton_pt_).size() > 3){
                                cout<< "lep3_pt: " << eve->lepton_pt_[3] << endl;
                                cout<< "lep3_eta: " << eve->lepton_eta_[3] << endl;
                                cout<< "lep3_phi: " << eve->lepton_phi_[3] << endl;
                                cout<<"lep3_scEta: "<<eve->lepton_scEta_[3]<<endl;
                              
				cout<<"lep3_relIso: "<< eve->lepton_relIso_[3]<<endl;
				cout<<"========================================================"<<endl;
			}

			cout<<"========================================================"<<endl;
			cout<<"========================================================"<<endl;


		}


		evtSF = 1.0;

        	NumberOfGoodJets = 0;
        	NumberOfLooseJets = 0;
		NumberOfBadJets = 0;
		NumberOfJets = 0;
        	NumberOfLeptons = 0;
        	NumberOfMuons = 0;
        	NumberOfElectrons = 0;


		is_ttlf = 0;
		is_ttb = 0;
		is_ttbb = 0;
		is_tt2b = 0;
		is_ttcc = 0;







		if (is_MC == 1 && TTbarMC == 1){
			additionalJetEventId = (eve->additionalJetEventId_)%100;
			if(additionalJetEventId == 0)is_ttlf = 1;
			else if (additionalJetEventId == 51)is_ttb = 1;
			else if (additionalJetEventId == 52)is_tt2b = 1;
			else if (additionalJetEventId == 53 || additionalJetEventId == 54 || additionalJetEventId == 55)is_ttbb = 1;
			else if (additionalJetEventId == 41 || additionalJetEventId == 42 || additionalJetEventId == 43 || additionalJetEventId == 45)is_ttcc = 1;
		}


		if(DEBUG == 1){
			cout << "is_MC = " << is_MC << endl;
			cout << "is_ttlf = " << is_ttlf << endl;
		}

        	MuonPT = 0;
		MuonIso = 0;
        	ElePT = 0;
		EleIso=0;
		EleIDSF=0;

        	MET = 0;
        	METphi = 0;
		NumberOfPV = 0;


		PotJetPT = 0;
        	JetPT = 0;
		JetCSV_b = 0;
		JetCSV_bb = 0;
		JetCSVsf = 0;
		CSVsf = 1;
        	SumJetPT = 0.0;
        	JetEta = 0;
        	SumJetEta = 0.0;

		sfJet0 = 0;
		sfJet1 = 0;
		sfJet2 = 0;
		sfJet3 = 0;


		flagEventEle = 0;
		flagEventMuon = 0;

		//additionalJetEventId = 0;

		passedSLeTrigOnly = 0;
		passedSLeTrig = 0;
		passedSLmuTrig = 0;
		passedDLeeTrig = 0;
		passedDLemuTrig = 0;
		passedDLmumuTrig = 0;


		passedHLT_Ele35_WPTight_Gsf_v_ = 0;
		passedHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_ = 0;
		passedHLT_IsoMu24_eta2p1_v_ = 0;
		passedHLT_IsoMu27_v_ = 0;
		passedHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = 0;
		passedHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0;
		passedHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = 0;
		passedHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0;
		passedHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0;
		passedHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = 0;
		passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_ = 0;
		passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_ = 0;
		passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_ = 0;







		if(i%5000 == 0){
			cout<<"Running over "<<i<<"th Event..."<<endl;
		}
		//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		// *** 3.1 Cuts on Trigger and NJets
		NumberOfPV = eve->numPVs_;
// CUT 1: Number of PV > = 1
		if(NumberOfPV == 0)continue;//cut
		passedHLT_Ele35_WPTight_Gsf_v_=eve->passHLT_Ele35_WPTight_Gsf_v_;
		passedHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_=eve->passHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_;
		passedHLT_IsoMu24_eta2p1_v_=eve->passHLT_IsoMu24_eta2p1_v_;
		passedHLT_IsoMu27_v_=eve->passHLT_IsoMu27_v_;
		passedHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_=eve->passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_;
		passedHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_=eve->passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_;
		passedHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_=eve->passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_;
		passedHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_=eve->passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_;
		passedHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_=eve->passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_;
		passedHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_=eve->passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_;
		passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_=eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_;
		passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_=eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_;
		passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_=eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_;

// CUT 2: Must Pass Single Lepton Trigger
 
		if(passedHLT_Ele35_WPTight_Gsf_v_  == 1) passedSLeTrig = 1;
		if(passedHLT_IsoMu24_eta2p1_v_ == 1 && ((eve->evt_ >= 297046 && eve->evt_ <= 303434) || is_MC == 1))passedSLmuTrig = 1;
		if(passedHLT_IsoMu27_v_ ==1)passedSLmuTrig = 1;
		if(passedHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_ == 1)passedDLeeTrig = 1;
		if(passedHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ == 1)passedDLeeTrig = 1;
		if(passedHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ == 1)passedDLemuTrig = 1;
		if(passedHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ == 1)passedDLemuTrig = 1;
		if(passedHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ == 1)passedDLemuTrig = 1;
		if(passedHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ == 1)passedDLemuTrig = 1;
		if ((foundPeriodB <= 100) && passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_ == 1)passedDLmumuTrig = 1;
		//if(passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_ == 1 && ((eve->evt_ >= 297046  && eve->evt_ <= 299329) || is_MC == 1))passedDLmumuTrig = 1;
		if(passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_ == 1 && (foundPeriodC <= 100 || foundPeriodD <= 100 || foundPeriodE <= 100 || foundPeriodF <= 100) )passedDLmumuTrig = 1;
		//if(passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_ == 1 && ((eve->evt_ >= 299368 && eve->evt_ <= 306462) || is_MC ==1 ))passedDLmumuTrig = 1;
		if(passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_ == 1 && (foundPeriodC <= 100 || foundPeriodD <= 100 || foundPeriodE <= 100 || foundPeriodF <= 100) )passedDLmumuTrig = 1;
		//if(passedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_ == 1 && ((eve->evt_ >= 299368 && eve->evt_ <= 306462) || is_MC ==1 ) )passedDLmumuTrig = 1;
		

		if(passedHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_ == 1 ){
			passedSLeTrigOnly = 1;
			//passedSLeTrig = 1;
					
			//Trigger not used for the DL channel.....
		}


		//Trigger Rule from gitlab ICHEP18 page
		if(DoubleMuonData == 1 && passedDLmumuTrig == 0)continue;
		if(DoubleEGData == 1 && passedDLeeTrig == 0 )continue;
		if(MuonEGData == 1 && passedDLemuTrig == 0 )continue;

                if(SingleElectronData == 1 && passedDLeeTrig == 1 && passedSLeTrig == 1)continue;
		if(SingleElectronData == 1 && passedDLeeTrig == 1 && passedSLeTrig == 0)continue;
		if(SingleElectronData == 1 && passedDLeeTrig == 0 && passedSLeTrig == 0)continue;
				

		if(SingleMuonData == 1 && passedDLmumuTrig == 1 && passedSLmuTrig == 1)continue;
		if(SingleMuonData == 1 && passedDLmumuTrig == 1 && passedSLmuTrig == 0)continue;
                if(SingleMuonData == 1 && passedDLmumuTrig == 0 && passedSLmuTrig == 0)continue;


/*		if(is_MC == 0){
			if(eve->passMETFilter_Flag_goodVertices_v_ == 0) continue;
                	//if(eve->passMETFilter_Flag_globalTightHalo2016Filter_v_ == 0) continue;
                	if(eve->passMETFilter_Flag_globalSuperTightHalo2016Filter_v_ == 0) continue;
                	if(eve->passMETFilter_Flag_HBHENoiseFilter_v_ == 0) continue;
                	if(eve->passMETFilter_Flag_HBHENoiseIsoFilter_v_ == 0) continue;
                	if(eve->passMETFilter_Flag_EcalDeadCellTriggerPrimitiveFilter_v_ == 0) continue;
                	if(eve->passMETFilter_Flag_BadPFMuonFilter_v_ == 0) continue;
                	if(eve->passMETFilter_Flag_BadChargedCandidateFilter_v_ == 0) continue;
                	if(eve->passMETFilter_Flag_eeBadScFilter_v_ == 0 && is_MC == 0 )continue;
			if(eve->passMETFilter_Flag_ecalBadCalibFilter_v_ == 0) continue;
		}
*/


 		cfTrig = cfTrig + 1;
		NumberOfJets = ((eve->jet_pt_[0]).size());     
		if(DEBUG == 1)cout << "Number of Jets: "<< NumberOfJets << endl;
		//if(NumberOfJets < 4)continue;
		//cfJets = cfNJets + 1;
  		NumberOfLeptons = (eve->lepton_pt_).size();
		if(DEBUG == 1) cout << "Number of Leptons: "<< NumberOfLeptons << endl;
		NumberOfMuons = 0;
		NumberOfBadMuons = 0;
		NumberOfLooseMuons = 0;
		NumberOfGoodMuons = 0;
		NumberOfElectrons = 0;
		NumberOfLooseElectrons = 0;
		NumberOfGoodElectrons = 0;
		NumberOfBadElectrons = 0;

		NumberOfLooseLeptons = 0;

		WhichGoodElectrons.clear();
		WhichBadElectrons.clear();
		WhichLooseElectrons.clear();

		WhichGoodMuons.clear();
		WhichBadMuons.clear();
		WhichLooseMuons.clear();

		WhichGoodLeptons.clear();
		WhichLooseLeptons.clear();
		WhichBadLeptons.clear();
        	//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		// ***  3.2 Find which leptons pass the PT and Eta cuts
		if(DEBUG == 1) cout << "Just before the Lepton Loop " << endl;

		for (int nLept = 0; nLept < NumberOfLeptons;nLept++){
               		if(NumberOfLeptons < 1)break;
			//if(eve->lepton_isTight_[nLept] != 1){
			//	cout<<"=========================="<<endl;
			//	cout<<nLept<<" Lepton is not tight!"<<endl;
			//	WhichBadLeptons.push_back(nLept);
			//	NumberOfBadLeptons = NumberOfBadLeptons + 1;
			//	continue;
			//}
		 	if(eve->lepton_isMuon_[nLept] == 1){
                        	MuonPT = eve->lepton_pt_[nLept];
				MuonEta = eve->lepton_eta_[nLept];
				MuonIso = eve->lepton_relIso_[nLept];
//CUT 5 Loose Leptons: Muon
				if(eve->lepton_isTight_[nLept] != 1 || MuonPT < cutMinMuonPT || abs(MuonEta) > cutMinMuonEta || MuonIso > cutMinMuonIso){
					NumberOfBadMuons = NumberOfBadMuons + 1;
					NumberOfBadLeptons = NumberOfBadLeptons + 1;
					WhichBadLeptons.push_back(nLept);
					WhichBadMuons.push_back(nLept);
					continue;//cut
				}
				//NumberOfLooseMuons = NumberOfLooseMuons + 1;
				//WhichLooseMuons.push_back(nLept);
//CUT 3 Muon: PT > 30 GeV       
//CUT 3 Muon: |ETA| < 2.1			
				if(abs(MuonEta) > cutMuonEta){
                                        NumberOfBadLeptons = NumberOfBadLeptons + 1;
                                        WhichBadLeptons.push_back(nLept);				
					WhichBadMuons.push_back(nLept);
					NumberOfBadMuons = NumberOfBadMuons + 1;
					continue;//cut
				}
//CUT 3 Muon: Rel Iso < 0.15
				if(MuonIso > cutMinMuonIso){
			                NumberOfBadLeptons = NumberOfBadLeptons + 1;
                                        WhichBadLeptons.push_back(nLept);
					WhichBadMuons.push_back(nLept);
					NumberOfBadMuons = NumberOfBadMuons + 1;
					continue;//cut
				}


	                           if(MuonPT < cutMuonPT){
                                        WhichLooseMuons.push_back(nLept);
                                        NumberOfLooseMuons = NumberOfLooseMuons + 1;
                                       NumberOfLooseLeptons = NumberOfLooseLeptons + 1;
                                        WhichLooseLeptons.push_back(nLept);
                                        continue;
                                }

                                WhichLooseMuons.push_back(nLept);
                                NumberOfLooseMuons = NumberOfLooseMuons + 1; 
                                NumberOfLooseLeptons = NumberOfLooseLeptons + 1;
                                WhichLooseLeptons.push_back(nLept);


				NumberOfGoodMuons = NumberOfGoodMuons +1;
				WhichGoodMuons.push_back(nLept);
				NumberOfGoodLeptons = NumberOfGoodLeptons + 1;
				WhichGoodLeptons.push_back(nLept);
                	}
               		 else if(eve->lepton_isMuon_[nLept] == 0){
                        	NumberOfElectrons=NumberOfElectrons+1;
                        	ElePT = eve->lepton_pt_[nLept];
				EleEta = eve->lepton_eta_[nLept];
				EleIso = eve->lepton_relIso_[nLept];				



//CUT 5 Loose Leptons: Ele				
	/*
				if(abs(eve->lepton_scEta_[nLept])<=1.479 && EleIso > (0.0287+(0.506/ElePT)) ){
					if(RunSingleEvent==1)cout<<"FAILED DUE TO EleIsoCut"<<endl<< "EleIso = "<<EleIso<<endl<<"EleIsoCut = "<< 0.0287+(0.506/ElePT) << endl;
					NumberOfBadLeptons = NumberOfBadLeptons + 1;
                                        WhichBadLeptons.push_back(nLept);
                                        WhichBadElectrons.push_back(nLept);
                                        NumberOfBadElectrons = NumberOfBadElectrons + 1;
                                        continue;
				}
				if(abs(eve->lepton_scEta_[nLept])>1.479 && EleIso > (0.0445+(0.963/ElePT)) ){
                              		if(RunSingleEvent==1)cout<<"FAILED DUE TO EleIsoCut"<<endl<< "EleIso = "<<EleIso<<endl<<"EleIsoCut = "<< 0.0445+(0.963/ElePT) << endl;
				          NumberOfBadLeptons = NumberOfBadLeptons + 1;
                                        WhichBadLeptons.push_back(nLept);
                                        WhichBadElectrons.push_back(nLept);
                                        NumberOfBadElectrons = NumberOfBadElectrons + 1;
                                        continue;
                                }
	*/
				if(eve->lepton_isTight_[nLept] != 1 || ElePT < cutMinElePT || abs(EleEta) > cutMinEleEta ){
                                        NumberOfBadLeptons = NumberOfBadLeptons + 1;
                                        WhichBadLeptons.push_back(nLept);			
					WhichBadElectrons.push_back(nLept);
					NumberOfBadElectrons = NumberOfBadElectrons + 1;
					continue;//cut
				}
				if(abs(eve->lepton_scEta_[nLept]) >=1.4442 && abs(eve->lepton_scEta_[nLept]) <= 1.5660){
                                        NumberOfBadLeptons = NumberOfBadLeptons + 1;
                                        WhichBadLeptons.push_back(nLept);
					WhichBadElectrons.push_back(nLept);
					NumberOfBadElectrons = NumberOfBadElectrons + 1;
					 continue;//cut
				}
				//NumberOfLooseLeptons = NumberOfLooseLeptons + 1;
				//WhichLooseLeptons.push_back(nLept);
				//NumberOfLooseElectrons = NumberOfLooseElectrons + 1;
				//WhichLooseElectrons.push_back(nLept);
//CUT 4 Ele: PT > 38 GeV
				
				if(ElePT < cutElePT){
                                        NumberOfLooseLeptons = NumberOfLooseLeptons + 1;
                                        WhichLooseLeptons.push_back(nLept);				
					WhichLooseElectrons.push_back(nLept);
					NumberOfLooseElectrons = NumberOfLooseElectrons + 1;
					continue;//cut
				}
				//if(EleIso > cutEleIso)continue;//cut
//CUT 4 Ele: |ETA| < 2.1
                        	if(abs(EleEta) > cutEleEta){
                                        NumberOfBadLeptons = NumberOfBadLeptons + 1;
                                        WhichBadLeptons.push_back(nLept);				
					WhichBadElectrons.push_back(nLept);
					NumberOfBadElectrons = NumberOfBadElectrons + 1;				
					continue;//cut
                      		}
			
                                        NumberOfLooseLeptons = NumberOfLooseLeptons + 1;
                                        WhichLooseLeptons.push_back(nLept);
                                        WhichLooseElectrons.push_back(nLept);
                                        NumberOfLooseElectrons = NumberOfLooseElectrons + 1;

				NumberOfGoodLeptons = NumberOfGoodLeptons + 1;
				NumberOfGoodElectrons = NumberOfGoodElectrons + 1;
				WhichGoodElectrons.push_back(nLept);
				WhichGoodLeptons.push_back(nLept);
                	}
        	}
		//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
//CUT 3,4 Muon, Ele: Must be at least one good lepton
		//if(NumberOfGoodElectrons == 0 && NumberOfGoodMuons == 0)continue;//cut
		//if(WhichGoodLeptons.size() == 0)continue;//cut
//CUT 5 Loose Leptons: Can only be one good electron
		//if(NumberOfLooseElectrons > 1) continue;//cut
		//if(NumberOfLooseElectrons == 1 && NumberOfGoodElectrons == 0)continue;//cut
//CUT 5 Loose Leptons: Can only be one good muon
		//if(NumberOfGoodElectrons > 2 || NumberOfGoodMuons > 2)continue;//cut
		//if(NumberOfLooseMuons > 2) continue;//cut
		//if(NumberOfLooseElectrons > 2)continue;
		
		if(DEBUG == 1){ 
			cout << "Out of the Lepton Loop " << endl;
			cout << "NumberOfGoodLeptons: "<< NumberOfGoodLeptons << endl;
			cout << "NumberOfLooseLeptons: "<< NumberOfLooseLeptons << endl;
			cout << "NumberOfBadLeptons: " << NumberOfBadLeptons << endl;
		}
		
		if(NumberOfLooseLeptons == 0)continue;
		//if(NumberOfLooseMuons == 1 && NumberOfGoodMuons == 2)continue;	//cut	
		//if(NumberOfLooseMuons >= 1 && NumberOfLooseElectrons >= 1)continue;//cut

		if(NumberOfGoodElectrons == 1) flagEventEle = 1;
		else if (NumberOfGoodMuons == 1) flagEventMuon = 1;

		cfTrigLep = cfTrigLep + 1;		

		GoodLeptonsE.clear();
		GoodLeptonsPT.clear();
		LooseLeptonsCharge.clear();
		LooseLeptonsE.clear();
		LooseLeptonsPT.clear();
		GoodLeptonsEta.clear();
		LooseLeptonsEta.clear();
		GoodLeptonsIso.clear();
		LooseLeptonsIso.clear();

		GoodMuonsCharge.clear();
		GoodMuonsE.clear();
		GoodMuonsPT.clear();

		LooseMuonsCharge.clear();
		LooseMuonsE.clear();
		LooseMuonsPT.clear();

		BadMuonsCharge.clear();
		BadMuonsE.clear();
		BadMuonsPT.clear();
		BadMuonsEta.clear();
		BadMuonsPhi.clear();


		GoodMuonsEta.clear();
		LooseMuonsEta.clear();
		GoodMuonsPhi.clear();
		LooseMuonsPhi.clear();
		GoodMuonsIso.clear();
		LooseMuonsIso.clear();

		GoodElectronsCharge.clear();
		GoodElectronsPT.clear();
		GoodElectronsE.clear();
		GoodElectronsIsoSF.clear();
		GoodElectronsIDSF.clear();

		LooseElectronsPT.clear();
		GoodElectronsEta.clear();
		LooseElectronsEta.clear();
		GoodElectronsPhi.clear();
		LooseElectronsPhi.clear();
		GoodElectronsIso.clear();
		LooseElectronsIso.clear();


		BadElectronsCharge.clear();
		BadElectronsE.clear();
		BadElectronsPT.clear();
		BadElectronsEta.clear();
		BadElectronsPhi.clear();
		BadElectronsIso.clear();




//===============
                for(int bl = 0; bl<WhichBadLeptons.size();bl++){
                        BadLeptonsCharge.push_back(eve->lepton_charge_[WhichBadLeptons.at(bl)]);
                        BadLeptonsPT.push_back(eve->lepton_pt_[WhichBadLeptons.at(bl)]);
                        BadLeptonsEta.push_back(eve->lepton_eta_[WhichBadLeptons.at(bl)]);
                        BadLeptonsPhi.push_back(eve->lepton_phi_[WhichBadLeptons.at(bl)]);
                        BadLeptonsIso.push_back(eve->lepton_relIso_[WhichBadLeptons.at(bl)]);

                }
                if(WhichBadMuons.size() != 0){
                        for(int b=0; b<WhichBadMuons.size();b++){
                                BadMuonsCharge.push_back(eve->lepton_charge_[WhichBadMuons.at(b)]);
                                BadMuonsPT.push_back(eve->lepton_pt_[WhichBadMuons.at(b)]);
                                BadMuonsEta.push_back(eve->lepton_eta_[WhichBadMuons.at(b)]);
                                BadMuonsPhi.push_back(eve->lepton_phi_[WhichBadMuons.at(b)]);
                                BadMuonsIso.push_back(eve->lepton_relIso_[WhichBadMuons.at(b)]);
                        }
                }
                if(WhichBadElectrons.size() != 0){
                        for(int b=0;b<WhichBadElectrons.size();b++){
                                BadElectronsCharge.push_back(eve->lepton_charge_[WhichBadElectrons.at(b)]);
                                BadElectronsPT.push_back(eve->lepton_pt_[WhichBadElectrons.at(b)]);
                                BadElectronsEta.push_back(eve->lepton_eta_[WhichBadElectrons.at(b)]);
                                BadElectronsPhi.push_back(eve->lepton_phi_[WhichBadElectrons.at(b)]);
                                BadElectronsIso.push_back(eve->lepton_relIso_[WhichBadElectrons.at(b)]);
                        }
                }
//================




		for(int ll = 0; ll<WhichLooseLeptons.size();ll++){
			LooseLeptonsCharge.push_back(eve->lepton_charge_[WhichLooseLeptons.at(ll)]);
			LooseLeptonsE.push_back(eve->lepton_e_[WhichLooseLeptons.at(ll)]);
			LooseLeptonsPT.push_back(eve->lepton_pt_[WhichLooseLeptons.at(ll)]);
			LooseLeptonsEta.push_back(eve->lepton_eta_[WhichLooseLeptons.at(ll)]);
                        LooseLeptonsPhi.push_back(eve->lepton_phi_[WhichLooseLeptons.at(ll)]);
                        LooseLeptonsIso.push_back(eve->lepton_relIso_[WhichLooseLeptons.at(ll)]);

		}
                if(WhichLooseMuons.size() != 0){
                        for(int l=0; l<WhichLooseMuons.size();l++){
                       		LooseMuonsCharge.push_back(eve->lepton_charge_[WhichLooseMuons.at(l)]);
			        LooseMuonsE.push_back(eve->lepton_e_[WhichLooseMuons.at(l)]);
				LooseMuonsPT.push_back(eve->lepton_pt_[WhichLooseMuons.at(l)]);
                                LooseMuonsEta.push_back(eve->lepton_eta_[WhichLooseMuons.at(l)]);
                                LooseMuonsPhi.push_back(eve->lepton_phi_[WhichLooseMuons.at(l)]);
                                LooseMuonsIso.push_back(eve->lepton_relIso_[WhichLooseMuons.at(l)]);
                        
                                        xbin_MuonIsoSF = -1;
                                        ybin_MuonIsoSF = -1;
                                        muonIsoSF_absEta = abs(eve->lepton_eta_[WhichLooseMuons.at(l)]);
                                        muonIsoSF_Pt = eve->lepton_pt_[WhichLooseMuons.at(l)];
                                for (int ybin = 0; ybin <= h_MuonIsoSF->GetNbinsY();ybin++){
                                        if(muonIsoSF_absEta > h_MuonIsoSF->GetYaxis()->GetBinLowEdge(ybin))ybin_MuonIsoSF = ybin;
                                        else continue;
                                        for (int xbin = 0; xbin <= h_MuonIsoSF->GetNbinsX();xbin++){
                                                if(muonIsoSF_Pt > h_MuonIsoSF->GetXaxis()->GetBinLowEdge(xbin))xbin_MuonIsoSF = xbin;
                                                else continue;
                                        }
                                }

                                LooseMuonsIsoSF.push_back(h_MuonIsoSF->GetBinContent(xbin_MuonIsoSF,ybin_MuonIsoSF));

                                if(var_SF == "muonSF" || var_SF =="AllSF"){
                                        //evtSF=evtSF*(h_MuonIsoSF->GetBinContent(xbin_MuonIsoSF,ybin_MuonIsoSF));
                                        h_sfMuonIso->Fill(h_MuonIsoSF->GetBinContent(xbin_MuonIsoSF,ybin_MuonIsoSF));
                                }

					xbin_MuonIDSF = -1;
                                        ybin_MuonIDSF = -1;
                                        muonIDSF_absEta = abs(eve->lepton_eta_[WhichLooseMuons.at(l)]);
                                        muonIDSF_Pt = eve->lepton_pt_[WhichLooseMuons.at(l)];
                                for (int ybin = 0; ybin <= h_MuonIDSF->GetNbinsY();ybin++){
                                        if(muonIDSF_absEta > h_MuonIDSF->GetYaxis()->GetBinLowEdge(ybin))ybin_MuonIDSF = ybin;
                                        else continue;
                                        for (int xbin = 0; xbin <= h_MuonIDSF->GetNbinsX();xbin++){
                                                if(muonIDSF_Pt > h_MuonIDSF->GetXaxis()->GetBinLowEdge(xbin))xbin_MuonIDSF = xbin;
                                                else continue;
                                        }
                                }
                                LooseMuonsIDSF.push_back(h_MuonIDSF->GetBinContent(xbin_MuonIDSF,ybin_MuonIDSF));

			}
                }                                                                                                                                                                                                
                if(WhichLooseElectrons.size() != 0){
                        for(int l=0;l<WhichLooseElectrons.size();l++){
                                LooseElectronsCharge.push_back(eve->lepton_charge_[WhichLooseElectrons.at(l)]);
				LooseElectronsE.push_back(eve->lepton_e_[WhichLooseElectrons.at(l)]);
				LooseElectronsPT.push_back(eve->lepton_pt_[WhichLooseElectrons.at(l)]);
                                LooseElectronsEta.push_back(eve->lepton_eta_[WhichLooseElectrons.at(l)]);
                                LooseElectronsPhi.push_back(eve->lepton_phi_[WhichLooseElectrons.at(l)]);
                                LooseElectronsIso.push_back(eve->lepton_relIso_[WhichLooseElectrons.at(l)]);
                        
					xbin_EleIsoSF = -1;
                                        ybin_EleIsoSF = -1;
                                        eleIsoSF_scEta = eve->lepton_scEta_[WhichLooseElectrons.at(l)];
                                        eleIsoSF_Pt = eve->lepton_pt_[WhichLooseElectrons.at(l)];
                                for (int xbin = 0; xbin <= h_EleIsoSF->GetNbinsX();xbin++){
                                        if(eleIsoSF_scEta > h_EleIsoSF->GetXaxis()->GetBinLowEdge(xbin))xbin_EleIsoSF = xbin;
                                        else continue;
                                        for (int ybin = 0; ybin <= h_EleIsoSF->GetNbinsY();ybin++){
                                                if(eleIsoSF_Pt > h_EleIsoSF->GetYaxis()->GetBinLowEdge(ybin))ybin_EleIsoSF = ybin;
                                                else continue;
                                        }
                                }
                                LooseElectronsIsoSF.push_back(h_EleIsoSF->GetBinContent(xbin_EleIsoSF,ybin_EleIsoSF));
                                if(var_SF == "eleSF" || var_SF =="AllSF"){
                                        //evtSF = evtSF*(h_EleIsoSF->GetBinContent(xbin_EleIsoSF,ybin_EleIsoSF));
                                        h_sfEleIso->Fill(h_EleIsoSF->GetBinContent(xbin_EleIsoSF,ybin_EleIsoSF));
                                }

					xbin_EleIDSF = -1;
                                        ybin_EleIDSF = -1;
                                        eleIDSF_scEta = eve->lepton_scEta_[WhichLooseElectrons.at(l)];
                                        eleIDSF_Pt = eve->lepton_pt_[WhichLooseElectrons.at(l)];
                                for (int xbin = 0; xbin <= h_EleIDSF->GetNbinsX();xbin++){
                                        if(eleIDSF_scEta > h_EleIDSF->GetXaxis()->GetBinLowEdge(xbin))xbin_EleIDSF = xbin;
                                        else continue;
                                        for (int ybin = 0; ybin <= h_EleIDSF->GetNbinsY();ybin++){
                                                if(eleIDSF_Pt > h_EleIDSF->GetYaxis()->GetBinLowEdge(ybin))ybin_EleIDSF = ybin;
                                                else continue;
                                        }
                                }
                                LooseElectronsIDSF.push_back(h_EleIDSF->GetBinContent(xbin_EleIDSF,ybin_EleIDSF));
                                if(var_SF == "eleSF" || var_SF =="AllSF"){
                                        //evtSF=evtSF*(h_EleIDSF->GetBinContent(xbin_EleIDSF,ybin_EleIDSF));
                                        h_sfEleID->Fill(h_EleIDSF->GetBinContent(xbin_EleIDSF,ybin_EleIDSF));
                                }
                                EleIDSF_Error = h_EleIDSF->GetBinError(xbin_EleIDSF,ybin_EleIDSF);
		

			}
                }



//==========
                indexBadLeptons.clear();
                if(BadLeptonsPT.size() > 0){
                        for (int i = 0 ; i != BadLeptonsPT.size() ; i++) {
                                indexBadLeptons.push_back(i);
                        }
                        sort(indexBadLeptons.begin(), indexBadLeptons.end(),[&](const int& a,const int&b){
                                return (BadLeptonsPT[a] > BadLeptonsPT[b]);
                        }
                        );
                        if(DEBUG == 2){
                                for (int i = 0 ; i != indexBadLeptons.size() ; i++) {
                                        cout << "Index: " << indexBadLeptons[i] << "  BadLeptonsPT: " << BadLeptonsPT.at(indexBadLeptons[i]) << endl;
                                }
                        }

                }

                indexBadMuons.clear();
                if(BadMuonsPT.size() > 0){
                        for (int i = 0 ; i != BadMuonsPT.size() ; i++) {
                                indexBadMuons.push_back(i);
                        }
                        sort(indexBadMuons.begin(), indexBadMuons.end(),[&](const int& a,const int&b){
                                return (BadMuonsPT[a] > BadMuonsPT[b]);
                        }
                        );
                        if(DEBUG == 2){
                                for (int i = 0 ; i != indexBadMuons.size() ; i++) {
                                        cout << "Index: " << indexBadMuons[i] << "  BadMuonsPT: " << BadMuonsPT.at(indexBadMuons[i]) << endl;
                                }
                        }

                }

                indexBadElectrons.clear();
                if(BadElectronsPT.size() > 0){
                        for (int i = 0 ; i != BadElectronsPT.size() ; i++) {
                                indexBadElectrons.push_back(i);
                        }
                        sort(indexBadElectrons.begin(), indexBadElectrons.end(),[&](const int& a,const int&b){
                                return (BadElectronsPT[a] > BadElectronsPT[b]);
                        }
                        );
                        if(DEBUG ==2){
                                for (int i = 0 ; i != indexBadElectrons.size() ; i++) {
                                        cout << "Index: " << indexBadElectrons[i] << "  BadElectronsPT: " << BadElectronsPT.at(indexBadElectrons[i]) << endl;
                                }
                        }

                }


//==========
                indexLooseLeptons.clear();
                if(LooseLeptonsPT.size() > 0){
                        for (int i = 0 ; i != LooseLeptonsPT.size() ; i++) {
                                indexLooseLeptons.push_back(i);
                        }
                        sort(indexLooseLeptons.begin(), indexLooseLeptons.end(),[&](const int& a,const int&b){
                                return (LooseLeptonsPT[a] > LooseLeptonsPT[b]);
                        }
                        );
                        if(DEBUG == 2){
                                for (int i = 0 ; i != indexLooseLeptons.size() ; i++) {
                                        cout << "Index: " << indexLooseLeptons[i] << "  LooseLeptonsPT: " << LooseLeptonsPT.at(indexLooseLeptons[i]) << endl;
                                }
                        }

                }
                indexLooseMuons.clear();
                if(LooseMuonsPT.size() > 0){
                        for (int i = 0 ; i != LooseMuonsPT.size() ; i++) {
                                indexLooseMuons.push_back(i);
                        }
                        sort(indexLooseMuons.begin(), indexLooseMuons.end(),[&](const int& a,const int&b){
                                return (LooseMuonsPT[a] > LooseMuonsPT[b]);
                        }
                        );
                        if(DEBUG == 2){
                                for (int i = 0 ; i != indexLooseMuons.size() ; i++) {
                                        cout << "Index: " << indexLooseMuons[i] << "  LooseMuonsPT: " << LooseMuonsPT.at(indexLooseMuons[i]) << endl;
                                }
                        }

                }

                indexLooseElectrons.clear();
                if(LooseElectronsPT.size() > 0){
                        for (int i = 0 ; i != LooseElectronsPT.size() ; i++) {
                                indexLooseElectrons.push_back(i);
                        }
                        sort(indexLooseElectrons.begin(), indexLooseElectrons.end(),[&](const int& a,const int&b){
                                return (LooseElectronsPT[a] > LooseElectronsPT[b]);
                        }
                        );
                        if(DEBUG ==2){
                                for (int i = 0 ; i != indexLooseElectrons.size() ; i++) {
                                        cout << "Index: " << indexLooseElectrons[i] << "  LooseElectronsPT: " << LooseElectronsPT.at(indexLooseElectrons[i]) << endl;
                                }
                        }

                }







		//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		// *** 3.3 Push the PT,Eta, and Phi variables from good leptons to vectors
		if(WhichGoodLeptons.size() != 0){
			for(int l = 0;l<WhichGoodLeptons.size();l++){
					GoodLeptonsCharge.push_back(eve->lepton_charge_[WhichGoodLeptons.at(l)]);
				 	GoodLeptonsE.push_back(eve->lepton_e_[WhichGoodLeptons.at(l)]);
					GoodLeptonsPT.push_back(eve->lepton_pt_[WhichGoodLeptons.at(l)]);
                        	        GoodLeptonsEta.push_back(eve->lepton_eta_[WhichGoodLeptons.at(l)]);
                        	        GoodLeptonsPhi.push_back(eve->lepton_phi_[WhichGoodLeptons.at(l)]);
					GoodLeptonsIso.push_back(eve->lepton_relIso_[WhichGoodLeptons.at(l)]);
			}	
		}

		xbin_MuonIsoSF = -1;
                ybin_MuonIsoSF = -1;
                muonIDSF_absEta = -99;
                muonIDSF_Pt = -99;

                xbin_MuonIDSF = -1;
                ybin_MuonIDSF = -1;
                muonIsoSF_absEta = -99;
                muonIsoSF_Pt = -99;




		if(WhichGoodMuons.size() != 0){
			for(int l=0; l<WhichGoodMuons.size();l++){
				GoodMuonsCharge.push_back(eve->lepton_charge_[WhichGoodMuons.at(l)]);
				GoodMuonsE.push_back(eve->lepton_e_[WhichGoodMuons.at(l)]);
				GoodMuonsPT.push_back(eve->lepton_pt_[WhichGoodMuons.at(l)]);
				GoodMuonsEta.push_back(eve->lepton_eta_[WhichGoodMuons.at(l)]);
				GoodMuonsPhi.push_back(eve->lepton_phi_[WhichGoodMuons.at(l)]);
				GoodMuonsIso.push_back(eve->lepton_relIso_[WhichGoodMuons.at(l)]);
			
				        xbin_MuonIsoSF = -1;
                                        ybin_MuonIsoSF = -1;
                                        muonIsoSF_absEta = abs(eve->lepton_eta_[WhichGoodMuons.at(l)]);
                                        muonIsoSF_Pt = eve->lepton_pt_[WhichGoodMuons.at(l)];
                                for (int ybin = 0; ybin <= h_MuonIsoSF->GetNbinsY();ybin++){
                                        if(muonIsoSF_absEta > h_MuonIsoSF->GetYaxis()->GetBinLowEdge(ybin))ybin_MuonIsoSF = ybin;
                                        else continue;
                                        for (int xbin = 0; xbin <= h_MuonIsoSF->GetNbinsX();xbin++){
                                                if(muonIsoSF_Pt > h_MuonIsoSF->GetXaxis()->GetBinLowEdge(xbin))xbin_MuonIsoSF = xbin;
                                                else continue;
                                        }
                                }

				GoodMuonsIsoSF.push_back(h_MuonIsoSF->GetBinContent(xbin_MuonIsoSF,ybin_MuonIsoSF));
                                
				if(var_SF == "muonSF" || var_SF =="AllSF"){
					//evtSF=evtSF*(h_MuonIsoSF->GetBinContent(xbin_MuonIsoSF,ybin_MuonIsoSF));
					h_sfMuonIso->Fill(h_MuonIsoSF->GetBinContent(xbin_MuonIsoSF,ybin_MuonIsoSF));
				}
				MuonIsoSF_Error = h_MuonIsoSF->GetBinError(xbin_MuonIsoSF,ybin_MuonIsoSF);
                                if(DEBUG==1){
					cout<<"&&&&&&&&&&&&&&&&&&"<<endl;
					cout<<"xbin_MuonIsoSF = "<< xbin_MuonIsoSF<<endl;
					cout<<"ybin_MuonIsoSF = "<<ybin_MuonIsoSF<<endl;
                                	cout<<"MuonIsoSF  =  "<<h_MuonIsoSF->GetBinContent(xbin_MuonIsoSF,ybin_MuonIsoSF)<<endl;
                                	cout<<"MuonIsoSF_Error = " << MuonIsoSF_Error <<endl;
                                	cout<<"&&&&&&&&&&&&&&&&&&&"<<endl;
				}

                                        xbin_MuonIDSF = -1;
                                        ybin_MuonIDSF = -1;
                                        muonIDSF_absEta = abs(eve->lepton_eta_[WhichGoodMuons.at(l)]);
                                        muonIDSF_Pt = eve->lepton_pt_[WhichGoodMuons.at(l)];
                                for (int ybin = 0; ybin <= h_MuonIDSF->GetNbinsY();ybin++){
                                        if(muonIDSF_absEta > h_MuonIDSF->GetYaxis()->GetBinLowEdge(ybin))ybin_MuonIDSF = ybin;
                                        else continue;
                                        for (int xbin = 0; xbin <= h_MuonIDSF->GetNbinsX();xbin++){
                                                if(muonIDSF_Pt > h_MuonIDSF->GetXaxis()->GetBinLowEdge(xbin))xbin_MuonIDSF = xbin;
                                                else continue;
                                        }
                                }
                                GoodMuonsIDSF.push_back(h_MuonIDSF->GetBinContent(xbin_MuonIDSF,ybin_MuonIDSF));
                                if(var_SF == "muonSF" || var_SF =="AllSF"){
					//evtSF=evtSF*(h_MuonIDSF->GetBinContent(xbin_MuonIDSF,ybin_MuonIDSF));
					h_sfMuonID->Fill(h_MuonIDSF->GetBinContent(xbin_MuonIDSF,ybin_MuonIDSF));
				}
				MuonIDSF_Error = h_MuonIDSF->GetBinError(xbin_MuonIDSF,ybin_MuonIDSF);
                                if(DEBUG==1){
					cout<<"&&&&&&&&&&&&&&&&&&"<<endl;
					cout<<"xbin_MuonIDSF = "<< xbin_MuonIDSF <<endl;
					cout<<"ybin_MuonIDSF = "<< ybin_MuonIDSF << endl;
                                	cout<<"MuonIDSF  =  "<<h_MuonIDSF->GetBinContent(xbin_MuonIDSF,ybin_MuonIDSF)<<endl;
                                	cout<<"MuonIDSF_Error = " << MuonIDSF_Error <<endl;
                                	cout<<"&&&&&&&&&&&&&&&&&&&"<<endl;
				}

			}
		}

		xbin_EleIsoSF = -1;
		ybin_EleIsoSF = -1;
		eleIDSF_scEta = -99;
		eleIDSF_Pt = -99;
		
		xbin_EleIDSF = -1;
		ybin_EleIDSF = -1;
		eleIsoSF_scEta = -99;
		eleIsoSF_Pt = -99;




		if(WhichGoodElectrons.size() != 0){
			for(int l=0;l<WhichGoodElectrons.size();l++){
                                GoodElectronsCharge.push_back(eve->lepton_charge_[WhichGoodElectrons.at(l)]);
				GoodElectronsE.push_back(eve->lepton_e_[WhichGoodElectrons.at(l)]);
				GoodElectronsPT.push_back(eve->lepton_pt_[WhichGoodElectrons.at(l)]);
                        	GoodElectronsEta.push_back(eve->lepton_eta_[WhichGoodElectrons.at(l)]);
				GoodElectronsPhi.push_back(eve->lepton_phi_[WhichGoodElectrons.at(l)]);
				GoodElectronsIso.push_back(eve->lepton_relIso_[WhichGoodElectrons.at(l)]);
				        xbin_EleIsoSF = -1;
        				ybin_EleIsoSF = -1;
        				eleIsoSF_scEta = eve->lepton_scEta_[WhichGoodElectrons.at(l)];
        				eleIsoSF_Pt = eve->lepton_pt_[WhichGoodElectrons.at(l)];
        			for (int xbin = 0; xbin <= h_EleIsoSF->GetNbinsX();xbin++){
                			if(eleIsoSF_scEta > h_EleIsoSF->GetXaxis()->GetBinLowEdge(xbin))xbin_EleIsoSF = xbin;
                			else continue;
                			for (int ybin = 0; ybin <= h_EleIsoSF->GetNbinsY();ybin++){
                        			if(eleIsoSF_Pt > h_EleIsoSF->GetYaxis()->GetBinLowEdge(ybin))ybin_EleIsoSF = ybin;
                        			else continue;
                			}
        			}
				GoodElectronsIsoSF.push_back(h_EleIsoSF->GetBinContent(xbin_EleIsoSF,ybin_EleIsoSF));
				if(var_SF == "eleSF" || var_SF =="AllSF"){
					//evtSF = evtSF*(h_EleIsoSF->GetBinContent(xbin_EleIsoSF,ybin_EleIsoSF));
					h_sfEleIso->Fill(h_EleIsoSF->GetBinContent(xbin_EleIsoSF,ybin_EleIsoSF));
				}
				EleIsoSF_Error = h_EleIsoSF->GetBinError(xbin_EleIsoSF,ybin_EleIsoSF);
				if(DEBUG==1){
					cout<<"&&&&&&&&&&&&&&&&&&"<<endl;
					cout<<"EleIsoSF  =  "<<h_EleIsoSF->GetBinContent(xbin_EleIsoSF,ybin_EleIsoSF)<<endl;
					cout<<"EleIsoSF_Error = " << EleIsoSF_Error <<endl;
					cout<<"&&&&&&&&&&&&&&&&&&&"<<endl;
				}

					xbin_EleIDSF = -1;
                                        ybin_EleIDSF = -1;
                                        eleIDSF_scEta = eve->lepton_scEta_[WhichGoodElectrons.at(l)];
                                        eleIDSF_Pt = eve->lepton_pt_[WhichGoodElectrons.at(l)];
                                for (int xbin = 0; xbin <= h_EleIDSF->GetNbinsX();xbin++){
                                        if(eleIDSF_scEta > h_EleIDSF->GetXaxis()->GetBinLowEdge(xbin))xbin_EleIDSF = xbin;
                                        else continue;
                                        for (int ybin = 0; ybin <= h_EleIDSF->GetNbinsY();ybin++){
                                                if(eleIDSF_Pt > h_EleIDSF->GetYaxis()->GetBinLowEdge(ybin))ybin_EleIDSF = ybin;
                                                else continue;
                                        }
                                }
                                GoodElectronsIDSF.push_back(h_EleIDSF->GetBinContent(xbin_EleIDSF,ybin_EleIDSF));
                                if(var_SF == "eleSF" || var_SF =="AllSF"){
					//evtSF=evtSF*(h_EleIDSF->GetBinContent(xbin_EleIDSF,ybin_EleIDSF));
					h_sfEleID->Fill(h_EleIDSF->GetBinContent(xbin_EleIDSF,ybin_EleIDSF));
				}
				EleIDSF_Error = h_EleIDSF->GetBinError(xbin_EleIDSF,ybin_EleIDSF);
                 		if(DEBUG==1){
			               cout<<"&&&&&&&&&&&&&&&&&&"<<endl;
                                	cout<<"EleIDSF  =  "<<h_EleIDSF->GetBinContent(xbin_EleIDSF,ybin_EleIDSF)<<endl;
                                	cout<<"EleIDSF_Error = " << EleIDSF_Error <<endl;
                                	cout<<"&&&&&&&&&&&&&&&&&&&"<<endl;
				}
			}
		//cout<<"EleSF_Error = " << EleSF_Error <<endl;
		
		}
		
		indexLeptons.clear();
                if(GoodLeptonsPT.size() > 0){
                        for (int i = 0 ; i != GoodLeptonsPT.size() ; i++) {
                                indexLeptons.push_back(i);
                        }
                        sort(indexLeptons.begin(), indexLeptons.end(),[&](const int& a,const int&b){
                                return (GoodLeptonsPT[a] > GoodLeptonsPT[b]);
                        }
                        );
                        if(DEBUG == 2){
                                for (int i = 0 ; i != indexLeptons.size() ; i++) {
                                        cout << "Index: " << indexLeptons[i] << "  GoodLeptonsPT: " << GoodLeptonsPT.at(indexLeptons[i]) << endl;
                                }
                        }

                }



		// *** Arrange Muon distributions py PT
                indexMuons.clear();
                if(GoodMuonsPT.size() > 0){
                        for (int i = 0 ; i != GoodMuonsPT.size() ; i++) {
                                indexMuons.push_back(i);
                        }
                        sort(indexMuons.begin(), indexMuons.end(),[&](const int& a,const int&b){
                                return (GoodMuonsPT[a] > GoodMuonsPT[b]);
                        }
                        );
                        if(DEBUG == 2){
				for (int i = 0 ; i != indexMuons.size() ; i++) {
                                	cout << "Index: " << indexMuons[i] << "  GoodMuonsPT: " << GoodMuonsPT.at(indexMuons[i]) << endl;
				}
                        }

                }
		// *** Arrange Electron distributions by PT
                indexElectrons.clear();
                if(GoodElectronsPT.size() > 0){
                        for (int i = 0 ; i != GoodElectronsPT.size() ; i++) {
                                indexElectrons.push_back(i);
                        }
                        sort(indexElectrons.begin(), indexElectrons.end(),[&](const int& a,const int&b){
                                return (GoodElectronsPT[a] > GoodElectronsPT[b]);
                        }
                        );
			if(DEBUG ==2){
				for (int i = 0 ; i != indexElectrons.size() ; i++) {
                                	cout << "Index: " << indexElectrons[i] << "  GoodElectronsPT: " << GoodElectronsPT.at(indexElectrons[i]) << endl;
                        	}
			}
			
                }

		//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		// *** 3.4 Looping over Jets
		SumJetPT = 0;
		SumJetEta = 0;
		NumberOfGoodJets = 0;
		NumberOfLooseJets = 0;
		NumberOfBadJets = 0;
	
		WhichGoodJets.clear();
		WhichLooseJets.clear();
		WhichBadJets.clear();

		if(DEBUG == 1) cout << "Starting the Jet Loop now over NumberOfJets="<<NumberOfJets << endl;

		//if(NumberOfJets < 2)continue;		
		for( int nObjA = 0;nObjA < NumberOfJets ;nObjA++){
			PotJetsPT.push_back(-1);
			PotJetsJERsf.push_back(1);

			// *** 4.A) Selection of Objects
			//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Checking for overlap with jets and leptons
			flagJetDR = 0;
			R_jetVeto = 0;
			deltaPhi = 0;
			deltaEta = 0;
			if(DEBUG==1)cout << "NumberOfLooseLeptons:           " << NumberOfLooseLeptons << endl;
			if(NumberOfLooseLeptons > 0){
			for(int r = 0;r < NumberOfLooseLeptons; r++){
				
				if(DEBUG==1){cout << "r is: " << r<< endl;
						cout<<"NumberOfLooseLeptons==="<<NumberOfLooseLeptons<<endl;
						cout<<"NumberOfLooseMuons==="<<NumberOfLooseMuons<<endl;
						cout<<"NumberOfLooseElectrons=="<<NumberOfLooseElectrons<<endl;
						cout<<"WhichLooseMuons.size()==== "<<WhichLooseMuons.size()<<endl;
						cout<<"WhichLooseElectrons.size()=== " <<WhichLooseElectrons.size()<<endl;
						cout<<"WhichLooseLeptons.size() === " << WhichLooseLeptons.size() << endl;
						cout << "indexLooseLeptons.size() === " << indexLooseLeptons.size() << endl;
						cout<<"indexLooseLeptons.at(r) ==== " << indexLooseLeptons.at(r) << endl;
						cout<< "WhichLooseLeptons.at(r) ==== " << WhichLooseLeptons.at(r) << endl;
				
				}
				deltaEta = eve->lepton_eta_[WhichLooseLeptons.at(r)] - eve->jet_eta_[0][nObjA];
				if(DEBUG==1)cout << "deltaEta: "<< deltaEta << endl;
				deltaPhi = abs(eve->lepton_phi_[WhichLooseLeptons.at(r)]-eve->jet_phi_[0][nObjA]);
				if(deltaPhi > 3.14159)deltaPhi = 2*3.14159 - deltaPhi;
				if(DEBUG==1)cout << "deltaPhi: " << deltaPhi << endl;
				R_jetVeto = sqrt(abs(pow(deltaEta,2) + pow(deltaPhi,2)));
				//sqrt(pow((eve->lepton_eta_[r] - eve->jet_eta_[0][nObjA]),2) - pow((eve->lepton_phi_[r]-eve->jet_phi_[0][nObjA]),2));
				if (R_jetVeto < 0.4){ 
					flagJetDR = 1;
					if(RunSingleEvent == 1){ cout<<"Jet Number: "<< nObjA << " was flagged for Lepton Overlap Veto with deltaR: " << R_jetVeto << endl;
				
						cout<<"############################################################"<<endl;
						cout << "Checking Jet "<< nObjA <<" against Loose Lepton "<< r<<endl;
						cout << "R_jetVeto: "<< R_jetVeto << endl;	
						cout << "Lepton Eta: " << eve->lepton_eta_[WhichLooseLeptons.at(r)] << endl;
						cout << "Jet Eta: " << eve->jet_eta_[0][nObjA] << endl;
						cout << "Lepton Phi: " << eve->lepton_phi_[WhichLooseLeptons.at(r)] << endl;
						cout << "Jet Phi: " << eve->jet_phi_[0][nObjA] << endl;
						cout<< "##########################################################"<<endl;
					}
				}
				//cout << "End of Loop Number: " << r << endl;
			}
			}
			if(flagJetDR == 1){
				NumberOfBadJets = NumberOfBadJets + 1;
				WhichBadJets.push_back(nObjA);		
				continue;
			}
			if(abs(eve->jet_eta_[0][nObjA]) > cutMinJetEta){
				NumberOfBadJets = NumberOfBadJets + 1;
				WhichBadJets.push_back(nObjA);
				continue;
			}
			if(DEBUG==1)cout<<"Looping over Jet Number:"<<nObjA<<endl;
			PotJetsPT.at(nObjA) = eve->jet_pt_[0][nObjA];
			//PotJetsJERsf.push_back(1);
			//if(abs(eve->jet_eta_[0][nObjA]) > cutMinJetEta)continue;//cut
			if(DEBUG==1)cout << "PotJetsPT being filled..." <<endl;
			//PotJetsPT.push_back(eve->jet_pt_[0][nObjA]);
	       		if(DEBUG==1)cout << "PotJetsPT just filled at iteration:"<<nObjA <<endl;
		        if(abs(eve->jet_eta_[0][nObjA]) < 0.522) PotJetsJERsf.at(nObjA) = 1.1595;
                        else if(0.522 > abs(eve->jet_eta_[0][nObjA]) && abs(eve->jet_eta_[0][nObjA]) < 0.783) PotJetsJERsf.at(nObjA) = 1.1948;
                        else if(0.738 > abs(eve->jet_eta_[0][nObjA]) && abs(eve->jet_eta_[0][nObjA]) < 1.131) PotJetsJERsf.at(nObjA) = 1.1464;
                        else if(1.131 > abs(eve->jet_eta_[0][nObjA]) && abs(eve->jet_eta_[0][nObjA]) < 1.305) PotJetsJERsf.at(nObjA) = 1.1609;
                        else if(1.305 > abs(eve->jet_eta_[0][nObjA]) && abs(eve->jet_eta_[0][nObjA]) < 1.740) PotJetsJERsf.at(nObjA) = 1.1278;
                        else if(1.740 > abs(eve->jet_eta_[0][nObjA]) && abs(eve->jet_eta_[0][nObjA]) < 1.930) PotJetsJERsf.at(nObjA) = 1.1000;
                        else if(1.930 > abs(eve->jet_eta_[0][nObjA]) && abs(eve->jet_eta_[0][nObjA]) < 2.043) PotJetsJERsf.at(nObjA) = 1.1426;
                        else if(2.043 > abs(eve->jet_eta_[0][nObjA]) && abs(eve->jet_eta_[0][nObjA]) < 2.322) PotJetsJERsf.at(nObjA) = 1.1512;
                        else if(2.322 > abs(eve->jet_eta_[0][nObjA]) && abs(eve->jet_eta_[0][nObjA]) < 2.500) PotJetsJERsf.at(nObjA) = 1.2963;
                        else if(2.500 > abs(eve->jet_eta_[0][nObjA]) && abs(eve->jet_eta_[0][nObjA]) < 2.853) PotJetsJERsf.at(nObjA) = 1.3418;
                        else if(2.853 > abs(eve->jet_eta_[0][nObjA])){
                                cout << "A JET JUST GOT THROUGH THAT WAS OUTSIDE AN ETA CUT..."<< endl;
                                cout << "JetEta: "<< eve->jet_eta_[0][nObjA] << endl;
                        }
	        
		        if(applyJERsf == 1){
                       		if(DEBUG == 1 && i%5000){
                                        cout<<"JetPT before JER SF: "<< PotJetsPT.at(nObjA) <<endl;
                                }
                                PotJetsPT.at(nObjA) = PotJetsPT.at(nObjA)*PotJetsJERsf.at(nObjA);
                                if(DEBUG == 1 && i%5000){
                                        cout<<"JetPT after JER SF: "<< PotJetsPT.at(nObjA) <<endl;
                                }
                        }
			if(DEBUG==1)cout<<"Just before the call on jet_eta_ cuts..."<<endl;		
			//PotJetPT = eve->jet_pt[0][j] * PotJetsJERsf.at(nObjA);
			if(eve->jet_eta_[0][nObjA] > cutJetEta || eve->jet_pt_[0][nObjA] < cutMinJetPT || eve->jet_puid_[0][nObjA]!=7 ){
				NumberOfBadJets = NumberOfBadJets + 1;
                                WhichBadJets.push_back(nObjA);
                                continue;
			}
			if(eve->jet_pt_[0][nObjA] < cutJetPT && eve->jet_pt_[0][nObjA] > cutMinJetPT && NumberOfGoodLeptons >= 1 && NumberOfLooseLeptons >= 2){
			//	NumberOfGoodJets = NumberOfGoodJets + 1;
			//	WhichGoodJets.push_back(nObjA);

				NumberOfLooseJets = NumberOfLooseJets + 1;
				WhichLooseJets.push_back(nObjA);

				continue;
			}	
			if(DEBUG==1)cout<<"Might be the PotsJetsPT not being filled...."<<endl;
			//if(PotJetsPT.at(nObjA) < cutJetPT){
                                //NumberOfLooseJets = NumberOfLooseJets + 1;
                                //WhichLooseJets.push_back(nObjA);
				//NumberOfBadJets = NumberOfBadJets + 1;
				//WhichBadJets.push_back(nObjA);
				//continue;//cut
			//}
				if(DEBUG == 1)cout << "Out of Range???" << endl;	
			NumberOfGoodJets = NumberOfGoodJets+1;
				if(DEBUG==1)cout << "Just before the WhichGoodJets cut..." << endl;
			WhichGoodJets.push_back(nObjA);
				if(DEBUG==1)cout<<"Survived the WhichGoodJets cut at nObjA="<<nObjA<<endl;
		}
		if (DEBUG == 1)cout << "Before the number of good jets cut..." << endl;
		//if(NumberOfGoodJets < cutNumbOfGoodJets){
		//	cout<<NumberOfGoodJets<<" Good Jets but require "<<cutNumbOfGoodJets<<" jets!!!!!!" << endl;
		 	//continue;//cut
		//}
		if (DEBUG == 1)cout<<"Made it past the number of good jets cut..."<<endl;
		cfTrigLepJet = cfTrigLepJet + 1;

		SumJetPT = 0;
		SumJetEta = 0;
		bTags = 0;
		JetNumber = 0;

		for(int d = 0;d<WhichBadJets.size();d++){
			JetNumber = WhichBadJets.at(d);
			BadJetsPT.push_back(eve->jet_pt_[0][JetNumber]);
                        BadJetsEta.push_back(eve->jet_eta_[0][JetNumber]);
                        BadJetsPhi.push_back(eve->jet_phi_[0][JetNumber]);
                        BadJetsDeepCSV_b.push_back(eve->jet_DeepCSV_b_[0][JetNumber]);			
		}
		JetNumber = 0;
		for(int e = 0;e<WhichLooseJets.size();e++){
                        JetNumber = WhichLooseJets.at(e);
                        LooseJetsPT.push_back(eve->jet_pt_[0][JetNumber]);
                        LooseJetsEta.push_back(eve->jet_eta_[0][JetNumber]);
                        LooseJetsPhi.push_back(eve->jet_phi_[0][JetNumber]);
                        LooseJetsDeepCSV_b.push_back(eve->jet_DeepCSV_b_[0][JetNumber]);
			LooseJetsDeepCSV_bb.push_back(eve->jet_DeepCSV_bb_[0][JetNumber]);
                }
                JetNumber = 0;

								if(DEBUG == 1)cout << "WhichGoodJets.size() = " << WhichGoodJets.size() << endl;
		for(int j = 0; j<WhichGoodJets.size(); j++){
			 		if(DEBUG == 1) cout << "Just started inside the WhichGoodJets loop...."<< endl;
			JetNumber = WhichGoodJets.at(j);
			if(DEBUG == 1) cout << "Assigned the j from WhichGoodJets"<< endl;
			GoodJetsPT.push_back(eve->jet_pt_[0][JetNumber]);
			GoodJetsCSVSF.push_back(2.22144*((1.+(0.540134*(eve->jet_pt_[0][JetNumber])))/(1.+(1.30246*(eve->jet_pt_[0][JetNumber])))) );
                        //0.9201*((1.+(0.0115429*x))/(1.+(0.0119144*x)))
				//Line 237//(0.972902+0.000201811*(eve->jet_pt_[0][JetNumber])+3.96396e-08*(eve->jet_pt_[0][JetNumber])*(eve->jet_pt_[0][JetNumber])+-4.53965e-10*(eve->jet_pt_[0][JetNumber])*(eve->jet_pt_[0][JetNumber])*(eve->jet_pt_[0][JetNumber]) );
				//Line 116//(1.00303*((1.+(0.0263384*(eve->jet_pt_[0][JetNumber])))/(1.+(0.0294513*(eve->jet_pt_[0][JetNumber])))) );
			//Line 78  Medium WP, combined..//(2.22144*((1.+(0.540134*(eve->jet_pt_[0][JetNumber])))/(1.+(1.30246*(eve->jet_pt_[0][JetNumber])))) );
				//Line 192//(0.562751*((1.+(0.509404*(eve->jet_pt_[0][JetNumber])))/(1.+(0.315111*(eve->jet_pt_[0][JetNumber])))) );
				//Line 244//0.744235+(0.959064/sqrt(eve->jet_pt_[0][JetNumber])));
				//Line 154//(0.9201*((1.+(0.0115429*(eve->jet_pt_[0][JetNumber])))/(1.+(0.0119144*(eve->jet_pt_[0][JetNumber])))));
			//0.9201*((1.+(0.0115429*x))/(1.+(0.0119144*x))) 
			if(var_SF == "csvSF" || var_SF =="AllSF"){
				//evtSF = evtSF * GoodJetsCSVSF.at(j);
				//CSVsf = CSVsf*GoodJetsCSVSF.at(j);
				CSVsf = CSVsf * eve->jet_DeepCSV_SF_[0][JetNumber];
				h_sfCSV->Fill(GoodJetsCSVSF.at(j));
			}
			SumJetPT = SumJetPT + eve->jet_pt_[0][JetNumber];
			GoodJetsJERsf.push_back(-1);
			GoodJetsEta.push_back(eve->jet_eta_[0][JetNumber]);	
//-------------------------------------
//***** MAY WANT TO MOVE THIS TO BEFORE ANY CUTS ARE MADE..... Aug 15
//---------------------------------
			if(DEBUG == 1) cout << "Calculating the GoodJetsJERsf from eta regions..."<< endl;
			if(abs(GoodJetsEta.at(j)) < 0.522) GoodJetsJERsf.push_back(1.1595);
			else if(0.522 > abs(GoodJetsEta.at(j)) && abs(GoodJetsEta.at(j)) < 0.783) GoodJetsJERsf.at(j) = 1.1948;
                        else if(0.738 > abs(GoodJetsEta.at(j)) && abs(GoodJetsEta.at(j)) < 1.131) GoodJetsJERsf.at(j) = 1.1464;
			else if(1.131 > abs(GoodJetsEta.at(j)) && abs(GoodJetsEta.at(j)) < 1.305) GoodJetsJERsf.at(j) = 1.1609;
                        else if(1.305 > abs(GoodJetsEta.at(j)) && abs(GoodJetsEta.at(j)) < 1.740) GoodJetsJERsf.at(j) = 1.1278;
                        else if(1.740 > abs(GoodJetsEta.at(j)) && abs(GoodJetsEta.at(j)) < 1.930) GoodJetsJERsf.at(j) = 1.1000;
                        else if(1.930 > abs(GoodJetsEta.at(j)) && abs(GoodJetsEta.at(j)) < 2.043) GoodJetsJERsf.at(j) = 1.1426;
                        else if(2.043 > abs(GoodJetsEta.at(j)) && abs(GoodJetsEta.at(j)) < 2.322) GoodJetsJERsf.at(j) = 1.1512;
                        else if(2.322 > abs(GoodJetsEta.at(j)) && abs(GoodJetsEta.at(j)) < 2.500) GoodJetsJERsf.at(j) = 1.2963;
			else if(2.500 > abs(GoodJetsEta.at(j)) && abs(GoodJetsEta.at(j)) < 2.853) GoodJetsJERsf.at(j) = 1.3418;
			else if(2.853 > abs(GoodJetsEta.at(j))){
				cout << "A JET JUST GOT THROUGH THAT WAS OUTSIDE AN ETA CUT..."<< endl;
				cout << "JetEta: "<< GoodJetsEta.at(j) << endl;
				break;
			}




			if(DEBUG == 1) cout << "Survived the JER filling..."<< endl;
			GoodJetsPhi.push_back(eve->jet_phi_[0][JetNumber]);
			GoodJetsM.push_back(eve->jet_m_[0][JetNumber]);
			SumJetEta = SumJetEta + eve->jet_eta_[0][JetNumber];
			GoodJetsCSVv2.push_back(eve->jet_combinedInclusiveSecondaryVertexV2BJetTags_[0][JetNumber]); 
			
			
			//GoodJetsDeepCSV_b.push_back(eve->jet_DeepCSV_b_[0][j]);
			//if statement for CSV scale factor HERE
			
			//if(output == "EleB"||output == "MuonB"){
			//		
			//}
			JetPT = eve->jet_pt_[0][JetNumber];
			JetCSV_b = eve->jet_DeepCSV_b_[0][JetNumber];
			JetCSV_bb = eve->jet_DeepCSV_bb_[0][JetNumber];
			JetCSVsf = eve->jet_DeepCSV_SF_[0][JetNumber];
/*
			if(foundPeriodB <= 50){
				JetCSVsf = 1.02475*((1.+(0.0059128*JetPT))/(1.+(0.00790218*JetPT)));
			}else if(foundPeriodC <= 50 || foundPeriodD <= 50 || foundPeriodE <= 50){
				JetCSVsf = 1.2493*((1.+(0.13089*JetPT))/(1.+(0.180274*JetPT)));
			}else if(foundPeriodF <= 50){
				JetCSVsf = 33.8749*((1.+(0.335857*JetPT))/(1.+(13.1409*JetPT)));
			}
*/
/*			else{//Use these for MC samples
				//JetCSVsf = 1.00303*((1.+(0.0263384*JetPT))/(1.+(0.0294513*JetPT)));//One from previous versions
				JetCSVsf=0.9201*((1.+(0.0115429*JetPT))/(1.+(0.0119144*JetPT)));//Using the DeepCSV_94XSF_v3_B_F.csv from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X for tightWP
			}
*/	

			if(applyJERsf == 1)GoodJetsPT.at(j) = GoodJetsPT.at(j)*GoodJetsJERsf.at(j);

		
	//=============Applying the CSV correction factors=========================================================
			if(DEBUG == 1) cout << "Filling in the DeepCSV_b scores for good jets..."<< endl;
			if (applyIndCSVsf == 1){
				GoodJetsDeepCSV_b.push_back(JetCSV_b*JetCSVsf);//eve->jet_DeepCSV_b_[0][j]);
				GoodJetsDeepCSV_bb.push_back(JetCSV_bb*JetCSVsf);
				GoodJetsSFDeepCSV.push_back(JetCSVsf);
				//if(CSVsf > 0)CSVsf = CSVsf*JetCSVsf;//Only calculates the JetCSVsf and the nominal CSVsf = Jet1csvsf*Jet2CSVsf.... if the flag is set
				//else CSVsf = JetCSVsf;
			}else{
				GoodJetsDeepCSV_b.push_back(JetCSV_b);
				GoodJetsDeepCSV_bb.push_back(JetCSV_bb);
				GoodJetsSFDeepCSV.push_back(1.00);
			}




		//==================================================
			if(DEBUG == 1) cout << "About to count b-tags according to DeepCSV cut..."<< endl;
			//if( (GoodJetsDeepCSV_b.at(j) + GoodJetsDeepCSV_bb.at(j)  ) > cutDeepCSV) bTags = bTags+1;
			//if( ( NumberOfGoodLeptons == 2 || (NumberOfLooseLeptons == 2 && NumberOfGoodLeptons == 1)) && ((LooseJetsDeepCSV_b.at(j) + LooseJetsDeepCSV_bb.at(j)) > cutDeepCSV) ) bTags = bTags+1; 
			//GoodJetsDeepCSV_b.push_back(eve->jet_DeepCSV_b_[0][JetNumber]);
		}

		//for(int lJets = 0; lJets < NumberOfLooseJets; lJets++){
		//	if( (LooseJetsDeepCSV_b.at(lJets) + LooseJetsDeepCSV_bb.at(lJets)) > cutDeepCSV) bTags = bTags + 1;
		//}

		if(DEBUG == 1) cout << "About to organize the jets according to PT..."<< endl;
		// *** Arrange Jets Distributions by PT
                indexJets.clear();
                if(GoodJetsPT.size() > 0){
		       for (int i = 0 ; i != NumberOfGoodJets; i++){
                                indexJets.push_back(i);
                        }
                        sort(indexJets.begin(), indexJets.end(),[&](const int& a,const int&b){
                                return (GoodJetsPT.at(a) > GoodJetsPT.at(b));
                        }
                        );
                        
                        if(DEBUG == 2){
				for (int i = 0 ; i != indexJets.size() ; i++) {
                                	cout << "Index: " << indexJets[i] << "  GoodJetssPT: " << GoodJetsPT.at(indexJets[i]) << endl;

                        	}
			}
                }
//===========
                indexLooseJets.clear();
                if(LooseJetsPT.size() > 0){
                       for (int i = 0 ; i != NumberOfLooseJets; i++){
                                indexLooseJets.push_back(i);
                        }
                        sort(indexLooseJets.begin(), indexLooseJets.end(),[&](const int& a,const int&b){
                                return (LooseJetsPT.at(a) > LooseJetsPT.at(b));
                        }
                        );

                        if(DEBUG == 2){
                                for (int i = 0 ; i != indexLooseJets.size() ; i++) {
                                        cout << "Index: " << indexLooseJets[i] << "  LooseJetsPT: " << LooseJetsPT.at(indexLooseJets[i]) << endl;

                                }
                        }
                }

//==============
		indexBadJets.clear();
                if(BadJetsPT.size() > 0){
                       for (int i = 0 ; i != NumberOfBadJets; i++){
                                indexBadJets.push_back(i);
                        }
                        sort(indexBadJets.begin(), indexBadJets.end(),[&](const int& a,const int&b){
                                return (BadJetsPT.at(a) > BadJetsPT.at(b));
                        }
                        );
                        
                        if(DEBUG == 2){
                                for (int i = 0 ; i != indexBadJets.size() ; i++) {
                                        cout << "Index: " << indexBadJets[i] << "  BadJetssPT: " << BadJetsPT.at(indexBadJets[i]) << endl;

                                }
                        }
                }


//============




		if(DEBUG == 1) cout << "Finished organizing the jets according to PT..."<< endl;
		//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		// *** 3.5  Final Checks 
		
		
		if(DEBUG ==1) cout<< "bTags: " << bTags << endl;

		//if(NumberOfGoodJets < cutNumbOfGoodJets)continue;//cut
		if(NumberOfGoodJets != WhichGoodJets.size()) cout<< "NumberOfGoodJets:"<<NumberOfGoodJets<< "  doesnt equal size of WhichGoodJets:"<< WhichGoodJets.size()<<endl;
		//if(bTags < cutNbTags) continue;//cut

		cfTrigLepJetBtag = cfTrigLepJetBtag + 1;
		

		if(DEBUG == 1)cout << "Filling Histograms now...." << endl;



		//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		// *** 3.6 Fill out distributions 

		//if(GoodMuonsPT.size() == 0 && GoodElectronsPT.size() == 0){
        	//        cout << endl<< endl << "NO ELECTRONS OR MUONS FOUND IN THE EVENT!!!!!!!!!!!" << endl << endl;
        	        //break;
        	//}else 
	/*
		if(GoodElectronsPT.size() > 0 && GoodMuonsPT.size() == 0){ 
			h_Ele0PT->Fill(GoodElectronsPT.at(indexElectrons[0]));
			h_Ele0PT_SF->Fill(GoodElectronsPT.at(indexElectrons[0]),evtSF);
			h_Ele0Eta->Fill(GoodElectronsEta.at(indexElectrons[0]));
			h_Ele0Eta_SF->Fill(GoodElectronsEta.at(indexElectrons[0]),evtSF);
			h_Ele0Phi->Fill(GoodElectronsPhi.at(indexElectrons[0]));
			h_Ele0Iso->Fill(GoodElectronsIso.at(indexElectrons[0]));
			if(DEBUG==1)cout<<"Good Electrons Iso: " << GoodElectronsIso.at(indexElectrons[0]) << endl;
		}
		else if(GoodMuonsPT.size() > 0 && GoodElectronsPT.size() == 0){
			h_Muon0PT->Fill(GoodMuonsPT.at(indexMuons[0]));
			h_Muon0PT_SF->Fill(GoodMuonsPT.at(indexMuons[0]),evtSF);
			h_Muon0Eta->Fill(GoodMuonsEta.at(indexMuons[0]));
			h_Muon0Eta_SF->Fill(GoodMuonsEta.at(indexMuons[0]),evtSF);
			h_Muon0Phi->Fill(GoodMuonsPhi.at(indexMuons[0]));
			h_Muon0Iso->Fill(GoodMuonsIso.at(indexMuons[0]));
			if(DEBUG==1)cout<<"Good Muons Iso: " << GoodMuonsIso.at(indexMuons[0]) << endl;	
		}
	*/
		//888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		// *** Filling leading Lepton histograms here, though should be the same as below with leading Electrons and Muons...
		if(DEBUG==1)cout<<"Filling Lepton Histograms..."<<endl;
                if(GoodLeptonsPT.size() > 0){
			if(DEBUG == 1){
				cout<< "indexLeptons[0]: "<< indexLeptons[0]<<endl;
				cout<< "GoodLeptonsPT: "  << GoodLeptonsPT.at(indexLeptons[0]) <<endl;
			}
			//h_Lep0PT->Fill(GoodLeptonsPT.at(indexLeptons[0]));
       	        	//h_Lep0Eta->Fill(GoodLeptonsEta.at(indexLeptons[0]));
                	//h_Lep0Phi->Fill(GoodLeptonsPhi.at(indexLeptons[0]));
                	//h_Lep0Iso->Fill(GoodLeptonsIso.at(indexLeptons[0]));		
		}

//888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
                if((NumberOfGoodLeptons >= 1 && NumberOfLooseLeptons >= 2  ) && NumberOfGoodJets >= 2  && NumberOfLooseJets > 0 ){
/*                        for(int lJets = 0; lJets < NumberOfLooseJets ;lJets++){
                                GoodJetsPT.push_back(LooseJetsPT.at(lJets));
                                GoodJetsEta.push_back(LooseJetsEta.at(lJets));
                                GoodJetsPhi.push_back(LooseJetsPhi.at(lJets));
                                GoodJetsDeepCSV_b.push_back(LooseJetsDeepCSV_b.at(lJets));
                                GoodJetsDeepCSV_bb.push_back(LooseJetsDeepCSV_bb.at(lJets));
                                WhichGoodJets.push_back(NumberOfGoodJets+lJets);
                        }
                        NumberOfGoodJets = NumberOfGoodJets + 1;
                }
*/
                JetNumber = 0;
                for(int e = 0;e<WhichLooseJets.size();e++){  //[ ] WHICH GOOD JETS     or indexJets?????
                        JetNumber = WhichLooseJets.at(e);
                        GoodJetsPT.push_back(eve->jet_pt_[0][JetNumber]);
                        GoodJetsEta.push_back(eve->jet_eta_[0][JetNumber]);
                        GoodJetsPhi.push_back(eve->jet_phi_[0][JetNumber]);
                        GoodJetsM.push_back(eve->jet_m_[0][JetNumber]);
			GoodJetsDeepCSV_b.push_back(eve->jet_DeepCSV_b_[0][JetNumber]);
                        GoodJetsDeepCSV_bb.push_back(eve->jet_DeepCSV_bb_[0][JetNumber]);
                	NumberOfGoodJets = NumberOfGoodJets + 1;
		}
		//NumberOfGoodJets = NumberOfGoodJets + 1;
                indexJets.clear();
                if(GoodJetsPT.size() > 0){
                       for (int i = 0 ; i != NumberOfGoodJets; i++){
                                indexJets.push_back(i);
                        }
                        sort(indexJets.begin(), indexJets.end(),[&](const int& a,const int&b){
                                return (GoodJetsPT.at(a) > GoodJetsPT.at(b));
                        }
                        );

                        if(DEBUG == 1){
                                for (int i = 0 ; i != indexJets.size() ; i++) {
                                        cout << "Index: " << indexJets[i] << "  GoodJetssPT: " << GoodJetsPT.at(indexJets[i]) << endl;

                                }
                        }
                }
	}//END EVENT RUNNING LOOP......
	for(int j=0;j<NumberOfGoodJets;j++){
		if( (GoodJetsDeepCSV_b.at(j) + GoodJetsDeepCSV_bb.at(j)  ) > cutDeepCSV) bTags = bTags+1;
	}
//888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888






		MET = eve->MET_Type1xy_sync_[0];
		METphi = eve->MET_Type1xy_phi_sync_[0];
		//if(MET<20)continue;
		//NumberOfPV = eve->numPVs_;		


		is_e = 0;
		is_mu = 0;
		is_ee = 0;
		is_emu = 0;
		is_mumu = 0;
		

		Lep0Charge = 0;
		Lep1Charge = 0;
	
		Ele0PT = -99;
		Ele0Eta = -99;

		Muon0PT = -99;
		Muon0Eta = -99;

		Lep0E = 0;
		Lep0PT = -99;
                Lep0Eta = -99;
		Lep0Phi = 0;
                Lep0Iso = 0;
		Lep0IsoSF = 0;
		Lep0IDSF = 0;


		Lep1E = 0;
                Lep1PT = 0;
                Lep1Eta = 0;
		Lep1Phi = 0;
                Lep1Iso = 0;
		Lep1IsoSF = 0;        
		Lep1IDSF = 0;

		csvSF = 0;
	        mll = 0;

		if(DEBUG==1){
			cout<<"About to check for e/mu/ee/mumu/emu conditions..."<<endl;		
			cout<<"NumberOfGoodJets: "<< NumberOfGoodJets <<endl;
			cout<<"NumberOfGoodElectrons: "<< NumberOfGoodElectrons <<endl;
			cout<<"NumberOfLooseElectrons: "<< NumberOfLooseElectrons <<endl;
			cout<<"NumberOfGoodMuons: "<< NumberOfGoodMuons <<endl;
			cout<<"NumberOfLooseMuons: "<< NumberOfLooseMuons <<endl;
			cout<<"NumberOfBtags: "<< bTags <<endl;
		}
		


		if(MET > 20 && NumberOfGoodJets >= 4 && NumberOfGoodElectrons == 1 && NumberOfGoodMuons == 0 && NumberOfLooseElectrons == 1 && NumberOfLooseMuons == 0  && bTags >= 2 && (passedSLeTrig == 1 || passedSLeTrigOnly ==1) ){
			if(eve->evt_ == EventNumber)cout<<"Event Flagged as SL(e)..."<<endl;
			is_e = 1;
			Lep0E = GoodElectronsE.at(indexElectrons[0]);
			Lep0PT = GoodElectronsPT.at(indexElectrons[0]);
			Lep0Eta = GoodElectronsEta.at(indexElectrons[0]);
			Lep0Phi = GoodElectronsPhi.at(indexElectrons[0]);
			Lep0Iso = GoodElectronsIso.at(indexElectrons[0]);
			Lep0IsoSF = GoodElectronsIsoSF.at(indexElectrons[0]);
			Lep0IDSF = GoodElectronsIDSF.at(indexElectrons[0]);

			Ele0PT = Lep0PT;Ele0Eta=Lep0Eta;

		}else if(MET > 20 && NumberOfGoodJets >= 4 && NumberOfGoodElectrons == 0 && NumberOfGoodMuons == 1 && NumberOfLooseElectrons == 0 && NumberOfLooseMuons == 1 && passedSLmuTrig == 1 && bTags >= 2 ){
			if(eve->evt_ == EventNumber)cout<<"Event Flagged as SL(mu)..."<<endl;
			is_mu = 1;
			Lep0E = GoodMuonsE.at(indexMuons[0]);
			Lep0PT = GoodMuonsPT.at(indexMuons[0]);
                        Lep0Eta = GoodMuonsEta.at(indexMuons[0]);
                        Lep0Phi = GoodMuonsPhi.at(indexMuons[0]);
			Lep0Iso = GoodMuonsIso.at(indexMuons[0]);
			Lep0IsoSF = GoodMuonsIsoSF.at(indexMuons[0]);
			Lep0IDSF = GoodMuonsIDSF.at(indexMuons[0]);


			Muon0PT = Lep0PT;Muon0Eta=Lep0Eta;


		}else if( MET > 40 && (NumberOfGoodElectrons == 2 || (NumberOfGoodElectrons == 1 && NumberOfLooseElectrons == 2)) && NumberOfGoodMuons == 0 && NumberOfLooseElectrons == 2  && NumberOfLooseMuons == 0 &&  bTags >= 1 && NumberOfGoodJets >= 2 && ((is_MC == 0 && passedDLeeTrig == 1 ) || (is_MC == 1  && (passedDLeeTrig == 1 || passedSLeTrig==1)))){

			if(eve->evt_ == EventNumber)cout<<"Event Flagged as DL(ee)..."<<endl;
			is_ee = 1;
			Lep0Charge = GoodElectronsCharge.at(indexElectrons[0]);
			Lep0E = GoodElectronsE.at(indexElectrons[0]);
			Lep0PT = GoodElectronsPT.at(indexElectrons[0]);
                        Lep0Eta = GoodElectronsEta.at(indexElectrons[0]);
                        Lep0Phi = GoodElectronsPhi.at(indexElectrons[0]);
			Lep0Iso = GoodElectronsIso.at(indexElectrons[0]);
			Lep0IsoSF = GoodElectronsIsoSF.at(indexElectrons[0]);			
			Lep0IDSF = GoodElectronsIDSF.at(indexElectrons[0]);

			Ele0PT = Lep0PT;Ele0Eta=Lep0Eta;


			if(NumberOfGoodElectrons==2){
				Lep1Charge = GoodElectronsCharge.at(indexElectrons[1]);
				Lep1E = GoodElectronsE.at(indexElectrons[1]);
				Lep1PT = GoodElectronsPT.at(indexElectrons[1]);
                        	Lep1Eta = GoodElectronsEta.at(indexElectrons[1]);
                        	Lep1Phi = GoodElectronsPhi.at(indexElectrons[1]);
				Lep1Iso = GoodElectronsIso.at(indexElectrons[1]);
				Lep1IsoSF = GoodElectronsIsoSF.at(indexElectrons[1]);
				Lep1IDSF = GoodElectronsIDSF.at(indexElectrons[1]);
			}
			if(NumberOfGoodElectrons == 1 && NumberOfLooseElectrons == 2){
				Lep1Charge = LooseElectronsCharge.at(indexLooseElectrons[1]);
				Lep1E = LooseElectronsE.at(indexLooseElectrons[1]);
				Lep1PT = LooseElectronsPT.at(indexLooseElectrons[1]);
                                Lep1Eta = LooseElectronsEta.at(indexLooseElectrons[1]);
                             	Lep1Phi = LooseElectronsPhi.at(indexLooseElectrons[1]);
				Lep1Iso = LooseElectronsIso.at(indexLooseElectrons[1]);
				Lep1IsoSF = LooseElectronsIsoSF.at(indexLooseElectrons[1]);
                        	Lep1IDSF = LooseElectronsIDSF.at(indexLooseElectrons[1]);
			}

			Lepton1.SetPtEtaPhiE(Lep0PT,Lep0Eta,Lep0Phi,Lep0E);
			Lepton2.SetPtEtaPhiE(Lep1PT,Lep1Eta,Lep1Phi,Lep1E);
			DiLepton = Lepton1 + Lepton2;			

			mll = DiLepton.M();
			if((mll > 76 && mll < 106) || mll < 20  )is_ee = 0;
			if(Lep0Charge == Lep1Charge)is_ee = 0;

/*
			if( AllSF==1 || var_TrigSFee == "DLee_el0pt_el1pt"){
				h_DLTriggerSF_ee_el0pt_el0eta
			}
        		else if(AllSF==1 || var_TrigSFee == "DLee_el0eta_el1eta"){

			}else if(AllSF == 1 || var_TrigSFee == "DLee_el0pt_el1eta"){

			}else if(AllSF == 1 || var_TrigSFee == "DLee_el1pt_el0eta"){

			}else cout << "No DLee SF applied" <<endl;
*/


		}else if(((NumberOfGoodElectrons == 1 && NumberOfGoodMuons == 1 && NumberOfLooseElectrons == 1 && NumberOfLooseMuons == 1 ) || (NumberOfGoodElectrons == 1 && NumberOfLooseElectrons == 1 && NumberOfLooseMuons == 1)||(NumberOfGoodMuons==1 && NumberOfLooseMuons == 1 && NumberOfLooseElectrons==1) ) && NumberOfGoodJets >= 2 && bTags >= 1  && (((is_MC == 1 && (passedDLemuTrig == 1 || passedSLeTrig == 1 || passedSLmuTrig == 1)) || ((MuonEGData == 1 && passedDLemuTrig == 1) || (SingleElectronData==1 && passedDLemuTrig == 0 && passedSLeTrig == 1 && passedSLmuTrig == 0) || (SingleMuonData == 1 && passedDLemuTrig == 0 && passedSLeTrig == 0 && passedSLmuTrig == 1 ))))){

			if(eve->evt_ == EventNumber)cout<<"&&&&&&&&   Event Flagged as DL(emu)...   &&&&&&&&&&&&"<<endl;
			is_emu = 1;
			if(NumberOfGoodElectrons == 1 && NumberOfGoodMuons == 1){
				if(eve->evt_ == EventNumber)cout<<" Good Electron and Good Muon...."<<endl;
				if(GoodElectronsPT.at(indexElectrons[0]) > GoodMuonsPT.at(indexMuons[0])){
					Lep0Charge = GoodElectronsCharge.at(indexElectrons[0]);
					 Lep0E = GoodElectronsE.at(indexElectrons[0]);
					Lep0PT = GoodElectronsPT.at(indexElectrons[0]);
                        		Lep0Eta = GoodElectronsEta.at(indexElectrons[0]);
                        		Lep0Phi = GoodElectronsPhi.at(indexElectrons[0]);
					Lep0Iso = GoodElectronsIso.at(indexElectrons[0]);			
					Lep0IsoSF = GoodElectronsIsoSF.at(indexElectrons[0]);				
					Lep0IDSF = GoodElectronsIDSF.at(indexElectrons[0]);
					
					Ele0PT = Lep0PT;Ele0Eta=Lep0Eta;
					
					Lep1Charge = GoodMuonsCharge.at(indexMuons[0]);
					Lep1E = GoodMuonsE.at(indexMuons[0]);
					Lep1PT = GoodMuonsPT.at(indexMuons[0]);
                	        	Lep1Eta = GoodMuonsEta.at(indexMuons[0]);
                	        	Lep1Phi = GoodMuonsPhi.at(indexMuons[0]);
					Lep1Iso = GoodMuonsIso.at(indexMuons[0]);
					Lep1IsoSF = GoodMuonsIsoSF.at(indexMuons[0]);
                        		Lep1IDSF = GoodMuonsIDSF.at(indexMuons[0]);
				
					Muon0PT = Lep1PT;Muon0Eta=Lep0Eta;

					 Lepton1.SetPtEtaPhiE(Lep0PT,Lep0Eta,Lep0Phi,Lep0E);
                        		Lepton2.SetPtEtaPhiE(Lep1PT,Lep1Eta,Lep1Phi,Lep1E);
                        		DiLepton = Lepton1 + Lepton2;
                        		mll = DiLepton.M();
					if(eve->evt_ == EventNumber)cout<<"&&&&&&&&   About to check the charge   &&&&&&&&&&&&"<<endl;
                        		if(mll < 20 ) is_emu = 0;
                        		if( Lep1Charge == Lep0Charge){
                                		is_emu = 0;
 			                 	if(eve->evt_ == EventNumber) cout<<"Same Charge for Leptons............Event will have the is_emu flag removed!"<<endl;
                        		}

	
					//mll = Lep0PT + Lep1PT;
			//		Lepton1.SetPtEtaPhiE(Lep0PT,Lep0Eta,Lep0Phi,Lep0E);
                        //		Lepton2.SetPtEtaPhiE(Lep1PT,Lep1Eta,Lep1Phi,Lep1E);
                        //		DiLepton = Lepton1 + Lepton2;

                        //		mll = DiLepton.M();
				}else{
					Lep0Charge = GoodMuonsCharge.at(indexMuons[0]);
					Lep0E = GoodMuonsE.at(indexMuons[0]);
					Lep0PT = GoodMuonsPT.at(indexMuons[0]);
                	                Lep0Eta = GoodMuonsEta.at(indexMuons[0]);
                	                Lep0Phi = GoodMuonsPhi.at(indexMuons[0]);
					Lep0Iso = GoodMuonsIso.at(indexMuons[0]);
					Lep0IsoSF = GoodMuonsIsoSF.at(indexMuons[0]);
                        		Lep0IDSF = GoodMuonsIDSF.at(indexMuons[0]);					
				
					Muon0PT = Lep0PT;Muon0Eta=Lep0Eta;
	
					Lep1Charge = GoodElectronsCharge.at(indexElectrons[0]);
					Lep1E = GoodElectronsE.at(indexElectrons[0]);
					Lep1PT = GoodElectronsPT.at(indexElectrons[0]);
                	                Lep1Eta = GoodElectronsEta.at(indexElectrons[0]);
                	                Lep1Phi = GoodElectronsPhi.at(indexElectrons[0]);
					Lep1Iso = GoodElectronsIso.at(indexElectrons[0]);
					Lep1IsoSF = GoodElectronsIsoSF.at(indexElectrons[0]);
					Lep1IDSF = GoodElectronsIDSF.at(indexElectrons[0]);

					Ele0PT = Lep1PT;Ele0Eta=Lep1Eta;


					Lepton1.SetPtEtaPhiE(Lep0PT,Lep0Eta,Lep0Phi,Lep0E);
                        		Lepton2.SetPtEtaPhiE(Lep1PT,Lep1Eta,Lep1Phi,Lep1E);
                        		DiLepton = Lepton1 + Lepton2;
                        		mll = DiLepton.M();
					if(eve->evt_ == EventNumber)cout<<"&&&&&&&&   MuonPT > ElectronPT   &&&&&&&&&&&&"<<endl;
                        		if(mll < 20 ) is_emu = 0;
                        		if( Lep1Charge == Lep0Charge){
                                		is_emu = 0;
                                		if(eve->evt_ == EventNumber)cout<<"Same Charge for Leptons............Event will have the is_emu flag removed!"<<endl;
                        		}
					if(eve->evt_ == EventNumber)cout<<"&&&&&&&&   Passed the    &&&&&&&&&&&&"<<endl;


					//mll = Lep0PT + Lep1PT;
			//		Lepton1.SetPtEtaPhiE(Lep0PT,Lep0Eta,Lep0Phi,Lep0E);
                        //		Lepton2.SetPtEtaPhiE(Lep1PT,Lep1Eta,Lep1Phi,Lep1E);
                        //		DiLepton = Lepton1 + Lepton2;

                        //		mll = DiLepton.M();

				}
			}else if(NumberOfGoodElectrons == 1 && NumberOfGoodMuons == 0){
					Lep0Charge = GoodElectronsCharge.at(indexElectrons[0]);    
				    	Lep0E = GoodElectronsE.at(indexElectrons[0]);
					Lep0PT = GoodElectronsPT.at(indexElectrons[0]);
                                        Lep0Eta = GoodElectronsEta.at(indexElectrons[0]);
                                      	Lep0Phi = GoodElectronsPhi.at(indexElectrons[0]);
					Lep0Iso = GoodElectronsIso.at(indexElectrons[0]);
					Lep0IsoSF = GoodElectronsIsoSF.at(indexElectrons[0]);
					Lep0IDSF = GoodElectronsIDSF.at(indexElectrons[0]);

					Ele0PT = Lep0PT;Ele0Eta=Lep0Eta;

					Lep1Charge = LooseMuonsCharge.at(indexLooseMuons[0]);
                                 	Lep1E = LooseMuonsE.at(indexLooseMuons[0]);
			       		Lep1PT = LooseMuonsPT.at(indexLooseMuons[0]);
                                        Lep1Eta = LooseMuonsEta.at(indexLooseMuons[0]);
                                        Lep1Phi = LooseMuonsPhi.at(indexLooseMuons[0]);
					Lep1Iso = LooseMuonsIso.at(indexLooseMuons[0]);
					Lep1IsoSF = LooseMuonsIsoSF.at(indexLooseMuons[0]);
                        		Lep1IDSF = LooseMuonsIDSF.at(indexLooseMuons[0]);



					Muon0PT = Lep1PT;Muon0Eta=Lep1Eta;
					

			 		Lepton1.SetPtEtaPhiE(Lep0PT,Lep0Eta,Lep0Phi,Lep0E);
                        		Lepton2.SetPtEtaPhiE(Lep1PT,Lep1Eta,Lep1Phi,Lep1E);
                        		DiLepton = Lepton1 + Lepton2;
                        		mll = DiLepton.M();

                        		if(mll < 20 ) is_emu = 0;
                        		if( Lep1Charge == Lep0Charge){
                                		is_emu = 0;
                                		if(eve->evt_ == EventNumber)cout<<"Same Charge for Leptons............Event will have the is_emu flag removed!"<<endl;
                        		}


			}else if(NumberOfGoodElectrons == 0 && NumberOfGoodMuons == 1){
                                        Lep0Charge = GoodMuonsCharge.at(indexMuons[0]);
					Lep0E = GoodMuonsE.at(indexMuons[0]);
					Lep0PT = GoodMuonsPT.at(indexMuons[0]);
                                        Lep0Eta = GoodMuonsEta.at(indexMuons[0]);
                                        Lep0Phi = GoodMuonsPhi.at(indexMuons[0]);
					Lep0Iso = GoodMuonsIso.at(indexMuons[0]);
					Lep0IsoSF = GoodMuonsIsoSF.at(indexMuons[0]);
                        		Lep0IDSF = GoodMuonsIDSF.at(indexMuons[0]);

					Muon0PT = Lep0PT;Muon0Eta=Lep0Eta;

					Lep1Charge = LooseElectronsCharge.at(indexLooseElectrons[0]);
                                        Lep1E = LooseElectronsE.at(indexLooseElectrons[0]);
					Lep1PT = LooseElectronsPT.at(indexLooseElectrons[0]);
                                        Lep1Eta = LooseElectronsEta.at(indexLooseElectrons[0]);
                                        Lep1Phi = LooseElectronsPhi.at(indexLooseElectrons[0]);
					Lep1Iso = LooseElectronsIso.at(indexLooseElectrons[0]);
					Lep1IsoSF = LooseElectronsIsoSF.at(indexLooseElectrons[0]);
                        		Lep1IDSF = LooseElectronsIDSF.at(indexLooseElectrons[0]);




					Ele0PT = Lep1PT;Ele0Eta=Lep1Eta;

					 Lepton1.SetPtEtaPhiE(Lep0PT,Lep0Eta,Lep0Phi,Lep0E);
                        		Lepton2.SetPtEtaPhiE(Lep1PT,Lep1Eta,Lep1Phi,Lep1E);
                        		DiLepton = Lepton1 + Lepton2;
                        		mll = DiLepton.M();

                        		if(mll < 20 ) is_emu = 0;
                        		if( Lep1Charge == Lep0Charge){
                                		is_emu = 0;
                                		if(eve->evt_ == EventNumber) cout<<"Same Charge for Leptons............Event will have the is_emu flag removed!"<<endl;
                        		}


			}else{
				cout<<"SOMETHING WENT WRONG WITH THE DL(emu) selection for event: "<< eve->evt_ <<endl;
			}
			//mll = Lep0PT + Lep1PT;
                        

		        //Lepton1.SetPtEtaPhiE(Lep0PT,Lep0Eta,Lep0Phi,Lep0E);
                        //Lepton2.SetPtEtaPhiE(Lep1PT,Lep1Eta,Lep1Phi,Lep1E);
                        //DiLepton = Lepton1 + Lepton2;
                        //mll = DiLepton.M();

                        if(mll < 20 ) is_emu = 0;
                        if( Lep1Charge == Lep0Charge){
                                is_emu = 0;
                                //cout<<"Same Charge for Leptons............Event will have the is_emu flag removed!"<<endl;
			}

		}else if(MET > 40 && NumberOfGoodElectrons == 0 && NumberOfLooseElectrons == 0 && (NumberOfGoodMuons == 2 || (NumberOfGoodMuons==1 && NumberOfLooseMuons==2)) && NumberOfLooseMuons == 2  && ((is_MC==1 && (passedDLmumuTrig == 1 || passedSLmuTrig == 1)) || (DoubleMuonData == 1 && passedDLmumuTrig == 1) || (SingleMuonData == 1 && passedDLmumuTrig == 0 && passedSLmuTrig == 1)) && NumberOfGoodJets >= 2 && bTags >=1 ){
			if(eve->evt_ == EventNumber)cout<<"Event Flagged as DL(mumu)..."<<endl;
			is_mumu = 1;
			if(NumberOfGoodMuons == 2){
				 Lep0Charge = GoodMuonsCharge.at(indexMuons[0]);
				Lep0E = GoodMuonsE.at(indexMuons[0]);
				Lep0PT = GoodMuonsPT.at(indexMuons[0]);
                        	Lep0Eta = GoodMuonsEta.at(indexMuons[0]);
                       		Lep0Phi = GoodMuonsPhi.at(indexMuons[0]);
			 	Lep0Iso = GoodMuonsIso.at(indexMuons[0]);
				Lep0IsoSF = GoodMuonsIsoSF.at(indexMuons[0]);
                        	Lep0IDSF = GoodMuonsIDSF.at(indexMuons[0]);

				Muon0PT = Lep0PT;Muon0Eta=Lep0Eta;

				Lep1Charge = GoodMuonsCharge.at(indexMuons[1]);
				Lep1E = GoodMuonsE.at(indexMuons[1]);
				Lep1PT = GoodMuonsPT.at(indexMuons[1]);
                        	Lep1Eta = GoodMuonsEta.at(indexMuons[1]);
                        	Lep1Phi = GoodMuonsPhi.at(indexMuons[1]);
				Lep1Iso = GoodMuonsIso.at(indexMuons[1]);
				Lep1IsoSF = GoodMuonsIsoSF.at(indexMuons[1]);
                        	Lep1IDSF = GoodMuonsIDSF.at(indexMuons[1]);

			}else if(NumberOfGoodMuons == 1 && NumberOfLooseMuons == 2){
				Lep0Charge = LooseMuonsCharge.at(indexLooseMuons[0]);
				Lep0E = LooseMuonsE.at(indexLooseMuons[0]);
				Lep0PT = LooseMuonsPT.at(indexLooseMuons[0]);
                                Lep0Eta = LooseMuonsEta.at(indexLooseMuons[0]);
                                Lep0Phi = LooseMuonsPhi.at(indexLooseMuons[0]);
				Lep0Iso = LooseMuonsIso.at(indexLooseMuons[0]);
				Lep0IsoSF = LooseMuonsIsoSF.at(indexLooseMuons[0]);
                        	Lep0IDSF = LooseMuonsIDSF.at(indexLooseMuons[0]);				


				Muon0PT = Lep0PT;Muon0Eta=Lep0Eta;

				Lep1Charge = LooseMuonsCharge.at(indexLooseMuons[1]);
                                Lep1E = LooseMuonsE.at(indexLooseMuons[1]);
				Lep1PT = LooseMuonsPT.at(indexLooseMuons[1]);
                                Lep1Eta = LooseMuonsEta.at(indexLooseMuons[1]);
                                Lep1Phi = LooseMuonsPhi.at(indexLooseMuons[1]);
				Lep1Iso = LooseMuonsIso.at(indexLooseMuons[1]);
				Lep1IsoSF = LooseMuonsIsoSF.at(indexLooseMuons[1]);
                        	Lep1IDSF = LooseMuonsIDSF.at(indexLooseMuons[1]);


			}
			Lepton1.SetPtEtaPhiE(Lep0PT,Lep0Eta,Lep0Phi,Lep0E);
                        Lepton2.SetPtEtaPhiE(Lep1PT,Lep1Eta,Lep1Phi,Lep1E);
                        DiLepton = Lepton1 + Lepton2;

                        mll = DiLepton.M();

		//	mll = Lep0PT + Lep1PT;
                        if((mll > 76 && mll < 106) || mll < 20){
				if(eve->evt_ == EventNumber)cout<< "EVENT FLAGGED FOR DILEPTON MASS"<<endl<<"  mll = "<<mll <<endl; //<<is_mumu = 0;
				mll = 0;
				is_mumu = 0;
			}
			if(Lep0Charge == Lep1Charge){
				if(eve->evt_ == EventNumber)cout<< "EVENT FLAGGED FOR SAME SIGN DILEPTONS"<<endl;
				is_mumu = 0;
			}
		}
		else if(eve->evt_==EventNumber)cout<<"EVENT WAS NOT FLAGGED FOR ANY CHANNEL!"<<endl;


		if(is_e == 1 & is_ee == 1){
			cout<<"Event: "<< eve->evt_ << " has been flagged as both SL(e) and DL(ee)"<<endl;
			break;
		}

                if(is_mu == 1 & is_mumu == 1){
                        cout<<"Event: "<< eve->evt_ << " has been flagged as both SL(mu) and DL(mumu)"<<endl;
                        break;
                }

                if(is_e == 1 & is_emu == 1){
                        cout<<"Event: "<< eve->evt_ << " has been flagged as both SL(e) and DL(emu)"<<endl;
                        break;
                }

                if(is_mu == 1 & is_emu == 1){
                        cout<<"Event: "<< eve->evt_ << " has been flagged as both SL(mu) and DL(emu)"<<endl;
                        break;
                }
		if(var_DL=="ee" && is_ee != 1){
			if(DEBUG==1)cout<<"EVENT WAS FLAGGED AS NOT AN ee EVENT, BREAKING OUT OF LOOP!"<<endl;
			continue;
		}
		if(var_DL=="emu" && is_emu != 1){
                        if(DEBUG==1)cout<<"EVENT WAS FLAGGED AS NOT AN emu EVENT, BREAKING OUT OF LOOP!"<<endl;
                        continue;
                }
		if(var_DL=="mumu" && is_mumu != 1){
                        if(DEBUG==1)cout<<"EVENT WAS FLAGGED AS NOT AN mumu EVENT, BREAKING OUT OF LOOP!"<<endl;
                        continue;
                }


		if(MuonEGData == 1 && (is_ee == 1 || is_mumu == 1)){
			if(DEBUG==1)cout<<"EVENT WAS FLAGGED AS NOT emu event coming from MuonEGData"<<endl;
			continue;
		}

					if(DEBUG==1)cout<<"GOING OVER TRIGGER SCALE FACTORS..."<<endl;
		TrigSF = 1;
                xbin_TrigSF = -1;
                ybin_TrigSF = -1;
		xValue_TrigSF = -99;
		yValue_TrigSF = -99;
		TrigSF_0Eta = Lep0Eta;
                TrigSF_0PT = Lep0PT;
		TrigSF_1Eta = Lep1Eta;
                TrigSF_1PT = Lep1PT;
		if(TrigSF_0PT < 20){
			TrigSF_0PT = 20.5;
			//cout<<endl<<endl<<"Leading Lepton PT is below 20GeV...."<<endl<<endl;
			//cout<< "Event Number: " << eve->evt_ << endl;
			continue;//break;//This was set to break for the full sync version//CHECKTHIS
		}	
		if(TrigSF_1PT < 20)TrigSF_1PT = 20.5;
		if(var_TrigSFxx=="ptpt"){
			xValue_TrigSF = TrigSF_0PT;yValue_TrigSF = TrigSF_1PT;
		}else if(var_TrigSFxx=="etaeta"){
			xValue_TrigSF = abs(Lep0Eta);yValue_TrigSF = abs(Lep1Eta);

		}else if(var_TrigSFxx=="0pt0eta"){
                        xValue_TrigSF = TrigSF_0PT;yValue_TrigSF = abs(Lep0Eta);

                }else if(var_TrigSFxx=="1pt1eta"){
                        xValue_TrigSF = TrigSF_1PT;yValue_TrigSF = abs(Lep1Eta);

                }
		//if(is_ee==1)h_DLTriggerSFll = (TH2D*)h_DLTriggerSFee->Clone();
		//else if(is_emu==1)h_DLTriggerSFll = (TH2D*)h_DLTriggerSFemu->Clone();
		//else if(is_mumu==1)h_DLTriggerSFll = (TH2D*)h_DLTriggerSFmumu->Clone();

					if(DEBUG==1)cout<<"RUNNING OVER TRIG 2D HISTOGRAMS...."<<endl;
	if(var_TrigSFxx != "noTrigSF"){
		if(is_ee==1){
						if(DEBUG==1)cout<<"JUST ENTERED THE LOOP FOR TRIG 2d histograms...."<<endl;
			for (int xbin = 0; xbin <= h_DLTriggerSFee->GetNbinsX();xbin++){
							//if(DEBUG==1)cout<<"WELCOME TO HELL!"<<endl;
                        	if(xValue_TrigSF > h_DLTriggerSFee->GetXaxis()->GetBinLowEdge(xbin))xbin_TrigSF = xbin;
                        		//if(DEBUG==1)cout<<"IN THE XAXIS LOOP!!!!"<<endl;
				else continue;
			}
                        for (int ybin = 0; ybin <= h_DLTriggerSFee->GetNbinsY();ybin++){
                                			//if(DEBUG==1)cout<<"IN THE YAXIS LOOP!!!!!"<<endl;
				if(yValue_TrigSF > h_DLTriggerSFee->GetYaxis()->GetBinLowEdge(ybin))ybin_TrigSF = ybin;
                                			//if(DEBUG==1)cout<<"RUNNING OVER yValue_TrigSF now...."<<endl;
				else continue;
                        }
                			if(DEBUG==1)cout<<"FINISHED RUNNING OVER YAXIS IN TRIG HISTOGRAM..."<<endl;
			
			TrigSF = h_DLTriggerSFee->GetBinContent(xbin_TrigSF,ybin_TrigSF);
			
			//cout<<endl<<endl<<endl<<endl;cout<<"TRIG SF = " << TrigSF << endl << xbin_TrigSF <<endl<<ybin_TrigSF<<endl; cout<<endl<<endl<<endl<<endl;
			
		//	h_Ele0PT->Fill(Lep0PT);
                  //      h_Ele0PT_SF->Fill(Lep0PT,evtSF);
                  //      h_Ele0Eta->Fill(Lep0Eta);
                  //      h_Ele0Eta_SF->Fill(Lep0Eta,evtSF);
                  //      h_Ele0Phi->Fill(Lep0Phi);
                  //      h_Ele0Iso->Fill(Lep0Iso);

		}

		if(is_emu==1){
                                                if(DEBUG==1)cout<<"JUST ENTERED THE LOOP FOR TRIG 2d histograms...."<<endl;
                        for (int xbin = 0; xbin <= h_DLTriggerSFemu->GetNbinsX();xbin++){
                                                        //if(DEBUG==1)cout<<"WELCOME TO HELL!"<<endl;
                                if(xValue_TrigSF > h_DLTriggerSFemu->GetXaxis()->GetBinLowEdge(xbin))xbin_TrigSF = xbin;
                                        //if(DEBUG==1)cout<<"IN THE XAXIS LOOP!!!!"<<endl;
                                else continue;
                        }
			for (int ybin = 0; ybin <= h_DLTriggerSFemu->GetNbinsY();ybin++){
                                                       // if(DEBUG==1)cout<<"IN THE YAXIS LOOP!!!!!"<<endl;
                                        if(yValue_TrigSF > h_DLTriggerSFemu->GetYaxis()->GetBinLowEdge(ybin))ybin_TrigSF = ybin;
                                                        //if(DEBUG==1)cout<<"RUNNING OVER yValue_TrigSF now...."<<endl;
                                        else continue;
                                }
                                        if(DEBUG==1)cout<<"FINISHED RUNNING OVER YAXIS IN TRIG HISTOGRAM..."<<endl;

			TrigSF = h_DLTriggerSFemu->GetBinContent(xbin_TrigSF,ybin_TrigSF);
                }
		if(is_mumu==1){
                                            //    if(DEBUG==1)cout<<"JUST ENTERED THE LOOP FOR TRIG 2d histograms...."<<endl;
                        for (int xbin = 0; xbin <= h_DLTriggerSFmumu->GetNbinsX();xbin++){
                                          //              if(DEBUG==1)cout<<"WELCOME TO HELL!"<<endl;
                                if(xValue_TrigSF > h_DLTriggerSFmumu->GetXaxis()->GetBinLowEdge(xbin))xbin_TrigSF = xbin;
                                        //if(DEBUG==1)cout<<"IN THE XAXIS LOOP!!!!"<<endl;
                                else continue;
                        }        
			for (int ybin = 0; ybin <= h_DLTriggerSFmumu->GetNbinsY();ybin++){
                                              //          if(DEBUG==1)cout<<"IN THE YAXIS LOOP!!!!!"<<endl;
                                if(yValue_TrigSF > h_DLTriggerSFmumu->GetYaxis()->GetBinLowEdge(ybin))ybin_TrigSF = ybin;
                                                //        if(DEBUG==1)cout<<"RUNNING OVER yValue_TrigSF now...."<<endl;
                                else continue;
                        }
                                        if(DEBUG==1)cout<<"FINISHED RUNNING OVER YAXIS IN TRIG HISTOGRAM..."<<endl;
                        
		TrigSF = h_DLTriggerSFmumu->GetBinContent(xbin_TrigSF,ybin_TrigSF);
                }
	}



		//if(TrigSF<=0)TrigSF = 1.0;
		evtSF=1.0;
		//TrigSF = h_DLTriggerSFll->GetBinContent(xbin_TrigSF,ybin_TrigSF);
		if((is_ee==1||is_emu==1||is_mumu==1) && (var_TrigSFxx=="ptpt"|| var_TrigSFxx=="etaeta" || var_TrigSFxx=="0pt0eta"||var_TrigSFxx=="1pt1eta")){
			evtSF= evtSF * TrigSF;
			//cout<<"evtSF: "<<evtSF<<endl;
			//cout<<"TrigSF: "<<TrigSF<<endl;
			h_sfTrig->Fill(TrigSF);
		}
		if (var_SF == "AllSF")evtSF = evtSF * CSVsf;

		evtSF = evtSF * Lep0IDSF * Lep0IsoSF * Lep1IDSF * Lep1IsoSF;	

		if(is_MC == 1 && wgt_generator_ > 0.001 && wgt_generator_ < 1.9999)evtSF = evtSF * wgt_generator_;


		//evtSF=1;
		//if(Lep0PT>200)cout<<"44444444444444444444444444444"<<endl<<"Leading Lepton PT: "<<Lep0PT<<" for event: "<<eve->evt_<<endl<< "has TrigSF: "<<TrigSF<<"44444444444444444444444444444"<<endl;
		
		if(eve->evt_ == EventNumber || eve->evt_%7777719==0){	//%10==0){
			cout << endl;
			cout<<"====================================="<<endl;
                        cout<<"EVENT DETAILS AFTER SELECTION"<<endl;
        		cout << "EVENT NUMBER: " << eve->evt_ << endl; 
			cout<< "wgt_generator_: " << wgt_generator_ << endl;
	               cout<<"====================================="<<endl;
			cout<<"is_e: " << is_e <<endl;
			cout<<"is_mu: " << is_mu <<endl;
			cout<<"is_ee: " << is_ee <<endl;
			cout<<"is_mumu: " << is_mumu <<endl;
			cout<<"is_emu: " << is_emu <<endl;
			cout<<"=====================================" << endl;
			cout<<"MET = MET_Type1xy_sync_ should be same as below: "<< MET<<endl;
			cout<<"MET_Type1xy: "<< eve->MET_Type1xy_[0]<<endl;
			cout<<"Number of PVs: "<<NumberOfPV<<endl;
			cout<<"Number of bTags: "<<bTags<<endl;
			cout<<"TrigSF: "<<TrigSF<<endl;
			cout<<"evtSF: "<<evtSF<<endl;
			cout<<"mll : "<<mll<<endl;
			cout<<"------------------------------------"<<endl;
			cout<<"Number of Good Jets: "<<NumberOfGoodJets<<endl;
			cout<<"0st Jet Statistics__________________"<<endl;
			if(NumberOfGoodJets >= 1){
			cout<<"GoodJetPT: "<<GoodJetsPT.at(indexJets[0])<<endl;
			cout<<"GoodJetEta: "<<GoodJetsEta.at(indexJets[0])<<endl;
			cout<<"GoodJetPhi: "<<GoodJetsPhi.at(indexJets[0])<<endl;
			cout<<"GoodJetDeepCSV_b: "<<GoodJetsDeepCSV_b.at(indexJets[0])<<endl;
			cout<<"GoodJetsCSVSF: "<<GoodJetsCSVSF.at(indexJets[0])<<endl;
			}
			if(NumberOfGoodJets >= 2){
                        cout<<"1st Jet Statistics__________________"<<endl;
                        cout<<"GoodJetPT: "<<GoodJetsPT.at(indexJets[1])<<endl;
                        cout<<"GoodJetEta: "<<GoodJetsEta.at(indexJets[1])<<endl;
                        cout<<"GoodJetPhi: "<<GoodJetsPhi.at(indexJets[1])<<endl;
                        cout<<"GoodJetDeepCSV_b: "<<GoodJetsDeepCSV_b.at(indexJets[1])<<endl;
			cout<<"GoodJetsCSVSF: "<<GoodJetsCSVSF.at(indexJets[1])<<endl;
			}
			if(NumberOfGoodJets >= 3){
                        cout<<"2 Jet Statistics__________________"<<endl;
                        cout<<"GoodJetPT: "<<GoodJetsPT.at(indexJets[2])<<endl;
                        cout<<"GoodJetEta: "<<GoodJetsEta.at(indexJets[2])<<endl;
                        cout<<"GoodJetPhi: "<<GoodJetsPhi.at(indexJets[2])<<endl;
                        cout<<"GoodJetDeepCSV_b: "<<GoodJetsDeepCSV_b.at(indexJets[2])<<endl;
			}  
                        if(NumberOfGoodJets >= 4){
                        cout<<"3 Jet Statistics__________________"<<endl;
                        cout<<"GoodJetPT: "<<GoodJetsPT.at(indexJets[3])<<endl;
                        cout<<"GoodJetEta: "<<GoodJetsEta.at(indexJets[3])<<endl;
                        cout<<"GoodJetPhi: "<<GoodJetsPhi.at(indexJets[3])<<endl;
                        cout<<"GoodJetDeepCSV_b: "<<GoodJetsDeepCSV_b.at(indexJets[3])<<endl;
                        }
                        if(NumberOfGoodJets >= 5){
                        cout<<"4 Jet Statistics__________________"<<endl;
                        cout<<"GoodJetPT: "<<GoodJetsPT.at(indexJets[4])<<endl;
                        cout<<"GoodJetEta: "<<GoodJetsEta.at(indexJets[4])<<endl;
                        cout<<"GoodJetPhi: "<<GoodJetsPhi.at(indexJets[4])<<endl;
                        cout<<"GoodJetDeepCSV_b: "<<GoodJetsDeepCSV_b.at(indexJets[4])<<endl;
                        }
                        if(NumberOfGoodJets >= 6){
                        cout<<"5 Jet Statistics__________________"<<endl;
                        cout<<"GoodJetPT: "<<GoodJetsPT.at(indexJets[5])<<endl;
                        cout<<"GoodJetEta: "<<GoodJetsEta.at(indexJets[5])<<endl;
                        cout<<"GoodJetPhi: "<<GoodJetsPhi.at(indexJets[5])<<endl;
                        cout<<"GoodJetDeepCSV_b: "<<GoodJetsDeepCSV_b.at(indexJets[5])<<endl;
                        }
                        if(NumberOfGoodJets >= 7){
                        cout<<"6 Jet Statistics__________________"<<endl;
                        cout<<"GoodJetPT: "<<GoodJetsPT.at(indexJets[6])<<endl;
                        cout<<"GoodJetEta: "<<GoodJetsEta.at(indexJets[6])<<endl;
                        cout<<"GoodJetPhi: "<<GoodJetsPhi.at(indexJets[6])<<endl;
                        cout<<"GoodJetDeepCSV_b: "<<GoodJetsDeepCSV_b.at(indexJets[6])<<endl;
                        }
                        if(NumberOfGoodJets >= 8){
                        cout<<"7 Jet Statistics__________________"<<endl;
                        cout<<"GoodJetPT: "<<GoodJetsPT.at(indexJets[7])<<endl;
                        cout<<"GoodJetEta: "<<GoodJetsEta.at(indexJets[7])<<endl;
                        cout<<"GoodJetPhi: "<<GoodJetsPhi.at(indexJets[7])<<endl;
                        cout<<"GoodJetDeepCSV_b: "<<GoodJetsDeepCSV_b.at(indexJets[7])<<endl;
                        }

			cout<<"========================================"<<endl;
			cout<<"Number of Loose Jets: "<<NumberOfLooseJets<<endl;
                        if(NumberOfLooseJets >= 1){
			cout<<"0st Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<LooseJetsPT.at(indexLooseJets[0])<<endl;
                        cout<<"JetEta: "<<LooseJetsEta.at(indexLooseJets[0])<<endl;
                        cout<<"JetPhi: "<<LooseJetsPhi.at(indexLooseJets[0])<<endl;
                        cout<<"JetDeepCSV_b: "<<LooseJetsDeepCSV_b.at(indexLooseJets[0])<<endl;
                        }
                        if(NumberOfLooseJets >= 2){
                        cout<<"1 Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<LooseJetsPT.at(indexLooseJets[1])<<endl;
                        cout<<"JetEta: "<<LooseJetsEta.at(indexLooseJets[1])<<endl;
                        cout<<"JetPhi: "<<LooseJetsPhi.at(indexLooseJets[1])<<endl;
                        cout<<"JetDeepCSV_b: "<<LooseJetsDeepCSV_b.at(indexLooseJets[1])<<endl;
                        }
                        if(NumberOfLooseJets >= 3){
                        cout<<"2 Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<LooseJetsPT.at(indexLooseJets[2])<<endl;
                        cout<<"JetEta: "<<LooseJetsEta.at(indexLooseJets[2])<<endl;
                        cout<<"JetPhi: "<<LooseJetsPhi.at(indexLooseJets[2])<<endl;
                        cout<<"JetDeepCSV_b: "<<LooseJetsDeepCSV_b.at(indexLooseJets[2])<<endl;
                        }
                        if(NumberOfLooseJets >= 4){
                        cout<<"3 Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<LooseJetsPT.at(indexLooseJets[3])<<endl;
                        cout<<"JetEta: "<<LooseJetsEta.at(indexLooseJets[3])<<endl;
                        cout<<"JetPhi: "<<LooseJetsPhi.at(indexLooseJets[3])<<endl;
                        cout<<"JetDeepCSV_b: "<<LooseJetsDeepCSV_b.at(indexLooseJets[3])<<endl;
                        }
			cout<<"=============================================="<<endl;
			cout<<"Number of Bad Jets: "<<NumberOfBadJets<<endl;
                        if(NumberOfBadJets >= 1){
                        cout<<"0st Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<BadJetsPT.at(indexBadJets[0])<<endl;
                        cout<<"JetEta: "<<BadJetsEta.at(indexBadJets[0])<<endl;
                        cout<<"JetPhi: "<<BadJetsPhi.at(indexBadJets[0])<<endl;
                        cout<<"JetDeepCSV_b: "<<BadJetsDeepCSV_b.at(indexBadJets[0])<<endl;
                        }
			if(NumberOfBadJets >= 2){
                        cout<<"1 Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<BadJetsPT.at(indexBadJets[1])<<endl;
                        cout<<"JetEta: "<<BadJetsEta.at(indexBadJets[1])<<endl;
                        cout<<"JetPhi: "<<BadJetsPhi.at(indexBadJets[1])<<endl;
                        cout<<"JetDeepCSV_b: "<<BadJetsDeepCSV_b.at(indexBadJets[1])<<endl;
                        }
                        if(NumberOfBadJets >= 3){
                        cout<<"2 Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<BadJetsPT.at(indexBadJets[2])<<endl;
                        cout<<"JetEta: "<<BadJetsEta.at(indexBadJets[2])<<endl;
                        cout<<"JetPhi: "<<BadJetsPhi.at(indexBadJets[2])<<endl;
                        cout<<"JetDeepCSV_b: "<<BadJetsDeepCSV_b.at(indexBadJets[2])<<endl;
                        }
                        if(NumberOfBadJets >= 4){
                        cout<<"3 Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<BadJetsPT.at(indexBadJets[3])<<endl;
                        cout<<"JetEta: "<<BadJetsEta.at(indexBadJets[3])<<endl;
                        cout<<"JetPhi: "<<BadJetsPhi.at(indexBadJets[3])<<endl;
                        cout<<"JetDeepCSV_b: "<<BadJetsDeepCSV_b.at(indexBadJets[3])<<endl;
                        }
                        if(NumberOfBadJets >= 5){
                        cout<<"4 Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<BadJetsPT.at(indexBadJets[4])<<endl;
                        cout<<"JetEta: "<<BadJetsEta.at(indexBadJets[4])<<endl;
                        cout<<"JetPhi: "<<BadJetsPhi.at(indexBadJets[4])<<endl;
                        cout<<"JetDeepCSV_b: "<<BadJetsDeepCSV_b.at(indexBadJets[4])<<endl;
                        }
                        if(NumberOfBadJets >= 6){
                        cout<<"5 Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<BadJetsPT.at(indexBadJets[5])<<endl;
                        cout<<"JetEta: "<<BadJetsEta.at(indexBadJets[5])<<endl;
                        cout<<"JetPhi: "<<BadJetsPhi.at(indexBadJets[5])<<endl;
                        cout<<"JetDeepCSV_b: "<<BadJetsDeepCSV_b.at(indexBadJets[5])<<endl;
                        }
                        if(NumberOfBadJets >= 7){
                        cout<<"6th Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<BadJetsPT.at(indexBadJets[6])<<endl;
                        cout<<"JetEta: "<<BadJetsEta.at(indexBadJets[6])<<endl;
                        cout<<"JetPhi: "<<BadJetsPhi.at(indexBadJets[6])<<endl;
                        cout<<"JetDeepCSV_b: "<<BadJetsDeepCSV_b.at(indexBadJets[6])<<endl;
                        }
                        if(NumberOfBadJets >= 8){
                        cout<<"7th Jet Statistics__________________"<<endl;
                        cout<<"JetPT: "<<BadJetsPT.at(indexBadJets[7])<<endl;
                        cout<<"JetEta: "<<BadJetsEta.at(indexBadJets[7])<<endl;
                        cout<<"JetPhi: "<<BadJetsPhi.at(indexBadJets[7])<<endl;
                        cout<<"JetDeepCSV_b: "<<BadJetsDeepCSV_b.at(indexBadJets[7])<<endl;
                        }
			cout<<"========================================="<<endl;
			cout<<"Number of Good Electrons: "<<NumberOfGoodElectrons<<endl;
			if(NumberOfGoodElectrons >= 1){
                        cout<<"0st Electron Statistics__________________"<<endl;
                        cout<<"EleCharge: "<<GoodElectronsCharge.at(indexElectrons[0])<<endl;
			cout<<"EleE: " <<GoodElectronsE.at(indexElectrons[0])<<endl;
			cout<<"ElePT: "<<GoodElectronsPT.at(indexElectrons[0])<<endl;
                        cout<<"EleEta: "<<GoodElectronsEta.at(indexElectrons[0])<<endl;
                        cout<<"ElePhi: "<<GoodElectronsPhi.at(indexElectrons[0])<<endl;
                        cout<<"EleIso: "<<GoodElectronsIso.at(indexElectrons[0])<<endl;
                        cout<<"EleIsoSF: "<<GoodElectronsIsoSF.at(indexElectrons[0])<<endl;
			cout<<"EleIDSF: "<<GoodElectronsIDSF.at(indexElectrons[0])<<endl;
			}
			if(NumberOfGoodElectrons >= 2){
                        cout<<"1 Electron Statistics__________________"<<endl;
                        cout<<"EleCharge: "<<GoodElectronsCharge.at(indexElectrons[1])<<endl;
			cout<<"EleE: "<<GoodElectronsE.at(indexElectrons[1])<<endl;
			cout<<"ElePT: "<<GoodElectronsPT.at(indexElectrons[1])<<endl;
                        cout<<"EleEta: "<<GoodElectronsEta.at(indexElectrons[1])<<endl;
                        cout<<"ElePhi: "<<GoodElectronsPhi.at(indexElectrons[1])<<endl;
                        cout<<"EleIso: "<<GoodElectronsIso.at(indexElectrons[1])<<endl;
                        cout<<"EleIsoSF: "<<GoodElectronsIsoSF.at(indexElectrons[1])<<endl;
			cout<<"EleIDSF: "<<GoodElectronsIDSF.at(indexElectrons[1])<<endl;
			}
			if(NumberOfGoodElectrons >= 3){
                        cout<<"2 Electron Statistics__________________"<<endl;
                        cout<<"EleCharge: "<<GoodElectronsCharge.at(indexElectrons[2])<<endl;
			cout<<"ElePT: "<<GoodElectronsPT.at(indexElectrons[2])<<endl;
                        cout<<"EleEta: "<<GoodElectronsEta.at(indexElectrons[2])<<endl;
                        cout<<"ElePhi: "<<GoodElectronsPhi.at(indexElectrons[2])<<endl;
                        cout<<"EleIso: "<<GoodElectronsIso.at(indexElectrons[2])<<endl;
                        }
			cout<<"--------------------------------------------------"<<endl;
			cout<<"Number of Loose Electrons: "<<NumberOfLooseElectrons<<endl;
			if(NumberOfLooseElectrons >= 1){
                        cout<<"0st Electron Statistics__________________"<<endl;
                        cout<<"EleCharge: "<<LooseElectronsCharge.at(indexLooseElectrons[0])<<endl;
			cout<<"EleE: "<<LooseElectronsE.at(indexLooseElectrons[0])<<endl;
			cout<<"ElePT: "<<LooseElectronsPT.at(indexLooseElectrons[0])<<endl;
                        cout<<"EleEta: "<<LooseElectronsEta.at(indexLooseElectrons[0])<<endl;
                        cout<<"ElePhi: "<<LooseElectronsPhi.at(indexLooseElectrons[0])<<endl;
                        cout<<"EleIso: "<<LooseElectronsIso.at(indexLooseElectrons[0])<<endl;
			}
			if(NumberOfLooseElectrons >= 2){
                        cout<<"1 Electron Statistics__________________"<<endl;
                        cout<<"EleCharge: "<<LooseElectronsCharge.at(indexLooseElectrons[1])<<endl;
			cout<<"EleE: "<<LooseElectronsE.at(indexLooseElectrons[1])<<endl;
			cout<<"ElePT: "<<LooseElectronsPT.at(indexLooseElectrons[1])<<endl;
                        cout<<"EleEta: "<<LooseElectronsEta.at(indexLooseElectrons[1])<<endl;
                        cout<<"ElePhi: "<<LooseElectronsPhi.at(indexLooseElectrons[1])<<endl;
                        cout<<"EleIso: "<<LooseElectronsIso.at(indexLooseElectrons[1])<<endl;
                        }
                        if(NumberOfLooseElectrons >= 3){
                        cout<<"2 Electron Statistics__________________"<<endl;
                        cout<<"EleCharge: "<<LooseElectronsCharge.at(indexLooseElectrons[2])<<endl;
			cout<<"ElePT: "<<LooseElectronsPT.at(indexLooseElectrons[2])<<endl;
                        cout<<"EleEta: "<<LooseElectronsEta.at(indexLooseElectrons[2])<<endl;
                        cout<<"ElePhi: "<<LooseElectronsPhi.at(indexLooseElectrons[2])<<endl;
                        cout<<"EleIso: "<<LooseElectronsIso.at(indexLooseElectrons[2])<<endl;
                        }
			cout<<"-------------------------------------------------"<<endl;
			cout<<"Number of Bad Electrons: "<<NumberOfBadElectrons<<endl;
			if(NumberOfBadElectrons >= 1){
                        cout<<"0st Electron Statistics__________________"<<endl;
			cout<<"EleCharge: "<<BadElectronsCharge.at(indexBadElectrons[0])<<endl;
                        cout<<"ElePT: "<<BadElectronsPT.at(indexBadElectrons[0])<<endl;
                        cout<<"EleEta: "<<BadElectronsEta.at(indexBadElectrons[0])<<endl;
                        cout<<"ElePhi: "<<BadElectronsPhi.at(indexBadElectrons[0])<<endl;
                        cout<<"EleIso: "<<BadElectronsIso.at(indexBadElectrons[0])<<endl;
                        }
			if(NumberOfBadElectrons >= 2){
                        cout<<"1 Electron Statistics__________________"<<endl;
                         cout<<"EleCharge: "<<BadElectronsCharge.at(indexBadElectrons[1])<<endl;
			cout<<"ElePT: "<<BadElectronsPT.at(indexBadElectrons[1])<<endl;
                        cout<<"EleEta: "<<BadElectronsEta.at(indexBadElectrons[1])<<endl;
                        cout<<"ElePhi: "<<BadElectronsPhi.at(indexBadElectrons[1])<<endl;
                        cout<<"EleIso: "<<BadElectronsIso.at(indexBadElectrons[1])<<endl;
                        }
	                if(NumberOfBadElectrons >= 3){
                        cout<<"2 Electron Statistics__________________"<<endl;
                         cout<<"EleCharge: "<<BadElectronsCharge.at(indexBadElectrons[2])<<endl;
			cout<<"ElePT: "<<BadElectronsPT.at(indexBadElectrons[2])<<endl;
                        cout<<"EleEta: "<<BadElectronsEta.at(indexBadElectrons[2])<<endl;
                        cout<<"ElePhi: "<<BadElectronsPhi.at(indexBadElectrons[2])<<endl;
                        cout<<"EleIso: "<<BadElectronsIso.at(indexBadElectrons[2])<<endl;
                        }
			cout<<"=================================================="<<endl;
                        cout<<"Number of Good Muons: "<<NumberOfGoodMuons<<endl;
			if(NumberOfGoodMuons >= 1){
                        cout<<"0st Muon Statistics__________________"<<endl;
                         cout<<"MuonCharge: "<<GoodMuonsCharge.at(indexMuons[0])<<endl;
			cout<<"MuonPT: "<<GoodMuonsPT.at(indexMuons[0])<<endl;
                        cout<<"MuonEta: "<<GoodMuonsEta.at(indexMuons[0])<<endl;
                        cout<<"MuonPhi: "<<GoodMuonsPhi.at(indexMuons[0])<<endl;
                        cout<<"MuonIso: "<<GoodMuonsIso.at(indexMuons[0])<<endl;
                        cout<<"MuonIsoSF: "<<GoodMuonsIsoSF.at(indexMuons[0])<<endl;
			cout<<"MuonIDSF: "<<GoodMuonsIDSF.at(indexMuons[0])<<endl;
			}
                        if(NumberOfGoodMuons >= 2){
                        cout<<"1 Muon Statistics__________________"<<endl;
                         cout<<"MuonCharge: "<<GoodMuonsCharge.at(indexMuons[1])<<endl;
			cout<<"MuonPT: "<<GoodMuonsPT.at(indexMuons[1])<<endl;
                        cout<<"MuonEta: "<<GoodMuonsEta.at(indexMuons[1])<<endl;
                        cout<<"MuonPhi: "<<GoodMuonsPhi.at(indexMuons[1])<<endl;
                        cout<<"MuonIso: "<<GoodMuonsIso.at(indexMuons[1])<<endl;
			cout<<"MuonIsoSF: "<<GoodMuonsIsoSF.at(indexMuons[1])<<endl;
			cout<<"MuonIDSF: "<<GoodMuonsIDSF.at(indexMuons[1])<<endl;
                        }
                        if(NumberOfGoodMuons >= 3){
                        cout<<"2 Muon Statistics__________________"<<endl;
			 cout<<"MuonCharge: "<<GoodMuonsCharge.at(indexMuons[2])<<endl;
                        cout<<"MuonPT: "<<GoodMuonsPT.at(indexMuons[2])<<endl;
                        cout<<"MuonEta: "<<GoodMuonsEta.at(indexMuons[2])<<endl;
                        cout<<"MuonPhi: "<<GoodMuonsPhi.at(indexMuons[2])<<endl;
                        cout<<"MuonIso: "<<GoodMuonsIso.at(indexMuons[2])<<endl;
                        }
                        if(NumberOfGoodMuons >= 4){
                        cout<<"3 Muon Statistics__________________"<<endl;
			 cout<<"MuonCharge: "<<GoodMuonsCharge.at(indexMuons[3])<<endl;
                        cout<<"MuonPT: "<<GoodMuonsPT.at(indexMuons[3])<<endl;
                        cout<<"MuonEta: "<<GoodMuonsEta.at(indexMuons[3])<<endl;
                        cout<<"MuonPhi: "<<GoodMuonsPhi.at(indexMuons[3])<<endl;
                        cout<<"MuonIso: "<<GoodMuonsIso.at(indexMuons[3])<<endl;
                        }
			cout<<"------------------------------------------------"<<endl;
                        cout<<"Number of Loose Muons: "<<NumberOfLooseMuons<<endl;
			if(NumberOfLooseMuons >= 1){
                        cout<<"0 Muon Statistics__________________"<<endl;
                         cout<<"MuonCharge: "<<LooseMuonsCharge.at(indexLooseMuons[0])<<endl;
			cout<<"MuonPT: "<<LooseMuonsPT.at(indexLooseMuons[0])<<endl;
                        cout<<"MuonEta: "<<LooseMuonsEta.at(indexLooseMuons[0])<<endl;
                        cout<<"MuonPhi: "<<LooseMuonsPhi.at(indexLooseMuons[0])<<endl;
                        cout<<"MuonIso: "<<LooseMuonsIso.at(indexLooseMuons[0])<<endl;
                        }
			if(NumberOfLooseMuons >= 2){
                        cout<<"1 Muon Statistics__________________"<<endl;
                        cout<<"MuonCharge: "<<LooseMuonsCharge.at(indexLooseMuons[1])<<endl;
			cout<<"MuonPT: "<<LooseMuonsPT.at(indexLooseMuons[1])<<endl;
                        cout<<"MuonEta: "<<LooseMuonsEta.at(indexLooseMuons[1])<<endl;
                        cout<<"MuonPhi: "<<LooseMuonsPhi.at(indexLooseMuons[1])<<endl;
                        cout<<"MuonIso: "<<LooseMuonsIso.at(indexLooseMuons[1])<<endl;
                        }
			if(NumberOfLooseMuons >= 3){
                        cout<<"2 Muon Statistics__________________"<<endl;
                        cout<<"MuonCharge: "<<LooseMuonsCharge.at(indexLooseMuons[2])<<endl;
			cout<<"MuonPT: "<<LooseMuonsPT.at(indexLooseMuons[2])<<endl;
                        cout<<"MuonEta: "<<LooseMuonsEta.at(indexLooseMuons[2])<<endl;
                        cout<<"MuonPhi: "<<LooseMuonsPhi.at(indexLooseMuons[2])<<endl;
                        cout<<"MuonIso: "<<LooseMuonsIso.at(indexLooseMuons[2])<<endl;
                        }
			if(NumberOfLooseMuons >= 4){
                        cout<<"3 Muon Statistics__________________"<<endl;
                        cout<<"MuonCharge: "<<LooseMuonsCharge.at(indexLooseMuons[3])<<endl;
			cout<<"MuonPT: "<<LooseMuonsPT.at(indexLooseMuons[3])<<endl;
                        cout<<"MuonEta: "<<LooseMuonsEta.at(indexLooseMuons[3])<<endl;
                        cout<<"MuonPhi: "<<LooseMuonsPhi.at(indexLooseMuons[3])<<endl;
                        cout<<"MuonIso: "<<LooseMuonsIso.at(indexLooseMuons[3])<<endl;
                        }
			cout<<"--------------------------------------------"<<endl;
                        cout<<"Number of Bad Muons: "<<NumberOfBadMuons<<endl;
			if(NumberOfBadMuons >= 1){
                        cout<<"0st Muon Statistics__________________"<<endl;
                        cout<<"MuonCharge: "<<BadMuonsCharge.at(indexBadMuons[0])<<endl;
			cout<<"MuonPT: "<<BadMuonsPT.at(indexBadMuons[0])<<endl;
                        cout<<"MuonEta: "<<BadMuonsEta.at(indexBadMuons[0])<<endl;
                        cout<<"MuonPhi: "<<BadMuonsPhi.at(indexBadMuons[0])<<endl;
                        cout<<"MuonIso: "<<BadMuonsIso.at(indexBadMuons[0])<<endl;
                        }
			if(NumberOfBadMuons >= 2){
                        cout<<"1 Muon Statistics__________________"<<endl;
                        cout<<"MuonCharge: "<<BadMuonsCharge.at(indexBadMuons[1])<<endl;
			cout<<"MuonPT: "<<BadMuonsPT.at(indexBadMuons[1])<<endl;
                        cout<<"MuonEta: "<<BadMuonsEta.at(indexBadMuons[1])<<endl;
                        cout<<"MuonPhi: "<<BadMuonsPhi.at(indexBadMuons[1])<<endl;
                        cout<<"MuonIso: "<<BadMuonsIso.at(indexBadMuons[1])<<endl;
                        }
			if(NumberOfBadMuons >= 3){
                        cout<<"2 Muon Statistics__________________"<<endl;
                        cout<<"MuonCharge: "<<BadMuonsCharge.at(indexBadMuons[2])<<endl;
			cout<<"MuonPT: "<<BadMuonsPT.at(indexBadMuons[2])<<endl;
                        cout<<"MuonEta: "<<BadMuonsEta.at(indexBadMuons[2])<<endl;
                        cout<<"MuonPhi: "<<BadMuonsPhi.at(indexBadMuons[2])<<endl;
                        cout<<"MuonIso: "<<BadMuonsIso.at(indexBadMuons[2])<<endl;
                        }
                        if(NumberOfBadMuons >= 4){
                        cout<<"3 Muon Statistics__________________"<<endl;
                        cout<<"MuonCharge: "<<BadMuonsCharge.at(indexBadMuons[3])<<endl;
			cout<<"MuonPT: "<<BadMuonsPT.at(indexBadMuons[3])<<endl;
                        cout<<"MuonEta: "<<BadMuonsEta.at(indexBadMuons[3])<<endl;
                        cout<<"MuonPhi: "<<BadMuonsPhi.at(indexBadMuons[3])<<endl;
                        cout<<"MuonIso: "<<BadMuonsIso.at(indexBadMuons[3])<<endl;
                        }
                        cout<<"------------------------------------"<<endl;
	
//=================================================	
			cout<< Lep0E << "   "<<Lep1E <<endl;
			cout << Lep0PT << "   "<<Lep1PT<<endl;
			cout << Lep0Eta << "   "<<Lep1Eta<<endl;
			cout << "M  " <<DiLepton.M() << endl;
			cout <<"Pt  "<< DiLepton.Pt() << endl;
			cout << "M2  "<<DiLepton.M2() << endl;
			cout << "E  " <<DiLepton.E() << endl;
			cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
			cout << endl;


		}	

		if(is_e != 1 && is_mu != 1 && is_emu != 1 && is_ee != 1 && is_mumu != 1){
			if(DEBUG==1)cout<<"Event: "<<eve->evt_ << " has not been flagged as an event for SL nor DL"<<endl;				
			continue;
		}


		//EMW Jan 6 2019

		if(Lep1PT < 15 ){
			is_emu = 0;
			is_ee = 0;
			is_mumu = 0;
			if(DEBUG==1)cout<<"Lep1PT less than 15GeV"<<endl;
			continue;
		}
		if(NumberOfGoodJets < cutNumbOfGoodJets){
			if(DEBUG==1)cout<<"Number of Good Jets less than cutNumbOfGoodJets"<<endl;
			continue;		
		}
		if(bTags < cutNbTags)continue;


		//BDT 2/21 EMW TODO

                selectedLeptonCharge.clear();
                selectedJetCSV.clear();

                selectedLeptonP4.push_back(Lepton1);
                selectedLeptonP4.push_back(Lepton2);
                selectedLeptonCharge.push_back(Lep0Charge);
                selectedLeptonCharge.push_back(Lep1Charge);
                //JetLorentzVector.clear();
                	if(DEBUG==1){
				cout<<"Starting the JetLorentzVector filling..."<<endl;
				cout << "NumberOfGoodJets: " << NumberOfGoodJets << endl;
			}
		for(int i = 0; i < NumberOfGoodJets;i++){
        		if(DEBUG==1){
				cout<<"i: "<< i << endl;
				cout<<"indexJets[i]: " << indexJets[i] << endl;
				cout <<"GoodJetsPT.at(indexJets[i]) = "<< GoodJetsPT.at(indexJets[i]) << endl;
				cout << "GoodJetsEta.at(indexJets[i]) = " << GoodJetsEta.at(indexJets[i]) << endl;
				cout << "GoodJetsPhi.at(indexJets[i]) = " << GoodJetsPhi.at(indexJets[i]) << endl;
				cout << "GoodJetsM.at(indexJets[i]) = "<< GoodJetsM.at(indexJets[i]) << endl;
			}
	                JetLorentzVector.SetPtEtaPhiM(GoodJetsPT.at(indexJets[i]),GoodJetsEta.at(indexJets[i]),GoodJetsPhi.at(indexJets[i]),GoodJetsM.at(indexJets[i]));
        
	                selectedJetP4.push_back(JetLorentzVector);
                        selectedJetCSV.push_back(GoodJetsDeepCSV_b.at(indexJets[i]) + GoodJetsDeepCSV_bb.at(indexJets[i]));
		}
		if(DEBUG==1)cout<<"Survived the Filling of the JetLorentzVector..."<<endl;
		metP4.SetPtEtaPhiM(MET,0,METphi,0);
	//mar1
	bdt_score = bdt.GetBDTOutput(selectedLeptonP4, selectedLeptonCharge, selectedJetP4, selectedJetCSV, metP4);

		if(DEBUG ==1 && bdt_score > -50 ){
			cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
			cout<<"Event Number: "<< eve->evt_ << endl;
			cout<<"NJets: " << NumberOfGoodJets << endl;
			cout<<"NbTags: " << bTags << endl;
			cout<<"bdt_score: " << bdt_score << endl;
			cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		}

		selectedLeptonP4.clear();
                selectedLeptonCharge.clear();
                selectedJetP4.clear();
                selectedJetCSV.clear();





		//if(isData == 0)additionalJetEventId = eve->additionalJetEventId_;
		if(isData == 1)additionalJetEventId = -9999;
//		if(DEBUG ==1)cout<<"Writing values to CSV file...."<<endl;

/*
		fCSV<< eve->run_<<","<<eve->lumi_<<","<<eve->evt_<<","<< is_e<<","<<is_mu<<","<<is_ee<<","<<is_emu<<","<<is_mumu<<","<<NumberOfGoodJets<<","<<bTags<<",";
		//if(failedEvent)failedCSV<< eve->run_<<","<<eve->lumi_<<","<<eve->evt_<<","<< is_e<<","<<is_mu<<","<<is_ee<<","<<is_emu<<","<<is_mumu<<","<<NumberOfGoodJets<<","<<bTags<<",";
		if(DEBUG ==1)cout<<"Writing values to CSV file 2...."<<endl;
	
		if(DEBUG ==1)cout<<"Writing values to CSV file 2.5...."<<endl;
		fCSV << Lep0PT << "," << Lep0Eta << ","<< Lep0Iso <<","<<"lep1_pdgId,"<<Lep0IDSF<<","<<Lep0IsoSF<<",";
                //if(failedEvent)failedCSV << Lep0PT << "," << Lep0Eta << ","<< Lep0Iso <<","<<"lep1_pdgId,lep1_idSF,lep1_isoSF,";	
		if(DEBUG ==1)cout<<"Writing values to CSV file 3...."<<endl;
		fCSV<< Lep1PT <<","<< Lep1Eta <<","<< Lep1Iso <<","<<"lep2_pdgId"<<","<<Lep1IDSF<<","<<Lep1IsoSF<<",";		     
                //if(failedEvent)failedCSV<< Lep1PT <<","<< Lep1Eta <<","<< Lep1Iso <<","<<"lep2_pdgId"<<","<<"lep2_idSF"<<","<<"lep2_isoSF"<<",";		
        	if(DEBUG ==1)cout<<"Writing values to CSV file 4...."<<endl;
		fCSV<< GoodJetsPT.at(indexJets[0]) << "," << GoodJetsEta.at(indexJets[0]) << "," << GoodJetsPhi.at(indexJets[0])<<","<<"jet1_jesSF"<<","<<"jet1_jesSF_up,"<<"jet1_jesSF_down,"<<"jet1_jesSF_PileUpDataMC_down,"<<"jet1_jesSF_RelativeFSR_up,"<<"jet1_jerSF_nominal,"<<GoodJetsDeepCSV_b.at(indexJets[0])<<","<<"jet1_PUJetId,"<<"jet1_PUJetDiscriminant,";
                //if(failedEvent == 1)failedCSV<< GoodJetsPT.at(indexJets[0]) << "," << GoodJetsEta.at(indexJets[0]) << "," << GoodJetsPhi.at(indexJets[0])<<","<<"jet1_jesSF"<<","<<"jet1_jesSF_up,"<<"jet1_jesSF_down,"<<"jet1_jesSF_PileUpDataMC_down,"<<"jet1_jesSF_RelativeFSR_up,"<<"jet1_jerSF_nominal,"<<GoodJetsDeepCSV_b.at(indexJets[0])<<","<<"jet1_PUJetId,"<<"jet1_PUJetDiscriminant,";

		if(DEBUG ==1)cout<<"Writing values to CSV file 5...."<<endl;
		fCSV<<GoodJetsPT.at(indexJets[1])<<","<<GoodJetsEta.at(indexJets[1])<<","<<GoodJetsPhi.at(indexJets[1])<<","<<"jet2_jesSF"<<","<<"jet2_jesSF_up,"<<"jet2_jesSF_down,"<<"jet2_jesSF_PileUpDataMC_down,"<<"jet2_jesSF_RelativeFSR_up,"<<"jet2_jerSF_nominal,"<<GoodJetsDeepCSV_b.at(indexJets[1])<<","<<"jet2_PUJetId,"<<"jet2_PUJetDiscriminant,";
                //if(failedEvent==1)failedCSV<<GoodJetsPT.at(indexJets[1])<<","<<GoodJetsEta.at(indexJets[1])<<","<<GoodJetsPhi.at(indexJets[1])<<","<<"jet2_jesSF"<<","<<"jet2_jesSF_up,"<<"jet2_jesSF_down,"<<"jet2_jesSF_PileUpDataMC_down,"<<"jet2_jesSF_RelativeFSR_up,"<<"jet2_jerSF_nominal,"<<GoodJetsDeepCSV_b.at(indexJets[1])<<","<<"jet2_PUJetId,"<<"jet2_PUJetDiscriminant,";

		if(DEBUG ==1)cout<<"Writing values to CSV file 6...."<<endl;
		fCSV<<MET<<","<<eve->MET_phi_[0]<<","<< mll << "," << additionalJetEventId << "," <<"ttHFGenFilterTag,"<<NumberOfPV<<","<<"puWeight,"<< CSVsf <<",csvSF_lf_up,csvSF_hf_down,csvSF_cErr1_down"<<endl;
		//if(failedEvent==1)failedCSV<<MET<<","<<eve->MET_phi_[0]<<","<< mll << "," << additionalJetEventId << "," <<"ttHFGenFilterTag,"<<NumberOfPV<<","<<"puWeight,"<< CSVsf <<",csvSF_lf_up,csvSF_hf_down,csvSF_cErr1_down"<<endl;

*/
		if(NumberOfGoodJets >= 2 && bTags >=1 )h_bdt_score_ge2jge1b->Fill(bdt_score);
		if(NumberOfGoodJets == 3 && bTags == 2)h_bdt_score_3j2b->Fill(bdt_score);
		if(NumberOfGoodJets == 3 && bTags == 3)h_bdt_score_3j3b->Fill(bdt_score);		
		if(NumberOfGoodJets >= 4 && bTags == 2)h_bdt_score_ge4j2b->Fill(bdt_score);
		if(NumberOfGoodJets >= 4 && bTags == 3)h_bdt_score_ge4j3b->Fill(bdt_score);
		if(NumberOfGoodJets >= 4 && bTags >= 4)h_bdt_score_ge4jge4b->Fill(bdt_score);

		if(NumberOfGoodJets >= 2 && bTags >=1 )h_bdt_score_ge2jge1b_SF->Fill(bdt_score,evtSF);
                if(NumberOfGoodJets == 3 && bTags == 2)h_bdt_score_3j2b_SF->Fill(bdt_score,evtSF);
                if(NumberOfGoodJets == 3 && bTags == 3)h_bdt_score_3j3b_SF->Fill(bdt_score,evtSF);                                
		if(NumberOfGoodJets >= 4 && bTags == 2)h_bdt_score_ge4j2b_SF->Fill(bdt_score,evtSF);
                if(NumberOfGoodJets >= 4 && bTags == 3)h_bdt_score_ge4j3b_SF->Fill(bdt_score,evtSF);
                if(NumberOfGoodJets >= 4 && bTags >= 4)h_bdt_score_ge4jge4b_SF->Fill(bdt_score,evtSF);




		h_MET->Fill(MET);
		h_MET_SF->Fill(MET,evtSF);
		h_NumberOfPV->Fill(NumberOfPV);
	
/*
		if((is_ee==1||is_emu==1||is_mumu==1) && (var_TrigSFxx=="ptpt"|| var_TrigSFxx=="etaeta" || var_TrigSFxx=="0pt0eta"||var_TrigSFxx=="1pt1eta")){
                        evtSF=evtSF*TrigSF;
                        h_sfTrig->Fill(TrigSF);
                }
*/

		h_Lep0PT->Fill(Lep0PT);
 		h_Lep0PT_SF->Fill(Lep0PT,evtSF);
		h_Lep0Eta->Fill(Lep0Eta);
		h_Lep0Eta_SF->Fill(Lep0Eta,evtSF);
		h_Lep0Phi->Fill(Lep0Phi);
		h_Lep0Iso->Fill(Lep0Iso);
		h_Lep0IsoSF->Fill(Lep0IsoSF);
		h_Lep0IDSF->Fill(Lep0IDSF);


		h_Lep1PT->Fill(Lep1PT);
		h_Lep1PT_SF->Fill(Lep1PT,evtSF);
		h_Lep1Eta->Fill(Lep1Eta);
		h_Lep1Eta_SF->Fill(Lep1Eta,evtSF);
		h_Lep1Phi->Fill(Lep1Phi);
		h_Lep1Iso->Fill(Lep1Iso);
		h_Lep1IsoSF->Fill(Lep1IsoSF);
		h_Lep1IDSF->Fill(Lep1IDSF);



		h_Ele0PT->Fill(Ele0PT);
                h_Ele0PT_SF->Fill(Ele0PT,evtSF);
                h_Ele0Eta->Fill(Ele0Eta);
                h_Ele0Eta_SF->Fill(Ele0Eta,evtSF);
		
		h_Muon0PT->Fill(Muon0PT);
		//cout<<"EventSF just before filling Muon0PT: "<<evtSF<<endl;
                //if(TrigSF==0){
		//	cout<<"TrigSF = 0  for event: "<<eve->evt_;
		//	cout<<"Muon0PT is "<<Muon0PT<<endl;
		//}
		h_Muon0PT_SF->Fill(Muon0PT,evtSF);
                h_Muon0Eta->Fill(Muon0Eta);
                h_Muon0Eta_SF->Fill(Muon0Eta,evtSF);
		h_sfEvent->Fill(evtSF);
		



		if(DEBUG==1)cout<<"Finished writing to CSV file...."<<endl;

		h_Jet0PT->Fill(GoodJetsPT.at(indexJets[0]));//0]);//LeadingJetPT);
		h_Jet0PT_SF->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
		h_Jet1PT->Fill(GoodJetsPT.at(indexJets[1]));//
		h_Jet1PT_SF->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
		//h_Jet2PT->Fill(GoodJetsPT.at(indexJets[2]));//2]);
		//h_Jet3PT->Fill(GoodJetsPT.at(indexJets[3]));//3]);//Jet3PT);
		                if(DEBUG==1)cout<<"Filling Jet 1 and 2 ETA histograms...."<<endl;
        	h_Jet0Eta->Fill(GoodJetsEta.at(indexJets[0]));//0]);//LeadingJetEta);
        	h_Jet1Eta->Fill(GoodJetsEta.at(indexJets[1]));//1]);
        	h_Jet0Eta_SF->Fill(GoodJetsEta.at(indexJets[0]),evtSF);
		h_Jet1Eta_SF->Fill(GoodJetsEta.at(indexJets[1]),evtSF);
		
		//h_Jet2Eta->Fill(GoodJetsEta.at(indexJets[2]));//2]);
        	//h_Jet3Eta->Fill(GoodJetsEta.at(indexJets[3]));//3]);//Jet3Eta);
		if(DEBUG==1)cout<<"Filling Jet 1 and 2 Phi histograms...."<<endl;

		h_Jet0Phi->Fill(GoodJetsPhi.at(indexJets[0]));//0]);//LeadingJetEta);
        	h_Jet1Phi->Fill(GoodJetsPhi.at(indexJets[1]));//1]);
        	//h_Jet2Phi->Fill(GoodJetsPhi.at(indexJets[2]));//2]);
        	//h_Jet3Phi->Fill(GoodJetsPhi.at(indexJets[3]));//3]);


		//h_Jet0combinedInclusiveSecondaryVertexV2BJetTags->Fill(GoodJetsCSVv2.at(indexJets[0]));//0]);//Jet0combinedInclusiveSecondaryVertexV2BJetTags);
		//h_Jet1combinedInclusiveSecondaryVertexV2BJetTags->Fill(GoodJetsCSVv2.at(indexJets[1]));//1]);//Jet1combinedInclusiveSecondaryVertexV2BJetTags);
		//h_Jet2combinedInclusiveSecondaryVertexV2BJetTags->Fill(GoodJetsCSVv2.at(indexJets[2]));//2]);//Jet2combinedInclusiveSecondaryVertexV2BJetTags);
		//h_Jet3combinedInclusiveSecondaryVertexV2BJetTags->Fill(GoodJetsCSVv2.at(indexJets[3]));//3]);//Jet3combinedInclusiveSecondaryVertexV2BJetTags);
		if(DEBUG==1)cout<<"Filling Jet 1 and 2 DeepCSV_b histograms...."<<endl;
	 	h_Jet0DeepCSV_b->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
	 	h_Jet1DeepCSV_b->Fill(GoodJetsDeepCSV_b.at(indexJets[1]));    	
		h_Jet0DeepCSV_b_SF->Fill(GoodJetsDeepCSV_b.at(indexJets[0]),evtSF);
                h_Jet1DeepCSV_b_SF->Fill(GoodJetsDeepCSV_b.at(indexJets[1]),evtSF);


		//h_Jet2DeepCSV_b->Fill(GoodJetsDeepCSV_b.at(indexJets[2]));//2]);//Jet2DeepCSV_b);
          	//h_Jet3DeepCSV_b->Fill(GoodJetsDeepCSV_b.at(indexJets[3]));//3]);//Jet3DeepCSV_b);
		if(DEBUG==1)cout<<"Filling Jet 1 and 2 DeepCSV histograms...."<<endl;
		h_Jet0SFDeepCSV->Fill(GoodJetsSFDeepCSV.at(indexJets[0]));
		h_Jet1SFDeepCSV->Fill(GoodJetsSFDeepCSV.at(indexJets[1]));
                //h_Jet2SFDeepCSV->Fill(GoodJetsSFDeepCSV.at(indexJets[2]));
                //h_Jet3SFDeepCSV->Fill(GoodJetsSFDeepCSV.at(indexJets[3]));

		if(DEBUG==1)cout<<"Filling Jet 1 and 2 DeepCSV_bb histograms...."<<endl;
		h_Jet0DeepCSV_bb->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]));
	 	h_Jet1DeepCSV_bb->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]));         	
		h_Jet0DeepCSV_bb_SF->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                h_Jet1DeepCSV_bb_SF->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);


		h_Jet0DeepCSV->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]));
		h_Jet1DeepCSV->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]));
		 h_Jet0DeepCSV_SF->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                h_Jet1DeepCSV_SF->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);



		//h_Jet2DeepCSV_bb->Fill(GoodJetsDeepCSV_bb.at(indexJets[2]));//2]);//Jet2DeepCSV_bb);
         	//h_Jet3DeepCSV_bb->Fill(GoodJetsDeepCSV_bb.at(indexJets[3]));//3]);//Jet3DeepCSV_bb);

		if(DEBUG==1)cout<<"Moving into number of jets greater than 3..."<<endl;

		if(NumberOfGoodJets > 4){
			//h_Jet4combinedInclusiveSecondaryVertexV2BJetTags->Fill(GoodJetsCSVv2.at(indexJets[4]));//4]);//Jet4combinedInclusiveSecondaryVertexV2BJetTags);
			h_Jet4DeepCSV_b->Fill(GoodJetsDeepCSV_b.at(indexJets[4]));//4]);
			h_Jet4DeepCSV_bb->Fill(GoodJetsDeepCSV_bb.at(indexJets[4]));//4]);
			h_Jet4PT->Fill(GoodJetsPT.at(indexJets[4]));//4]);
			h_Jet4Eta->Fill(GoodJetsEta.at(indexJets[4]));//4]);
			h_Jet4Phi->Fill(GoodJetsEta.at(indexJets[4]));//4]);
		}


		if(NumberOfGoodJets > 5){
        	//	h_Jet5combinedInclusiveSecondaryVertexV2BJetTags->Fill(GoodJetsCSVv2[indexJets[5]]);//5]);//Jet4combinedInclusiveSecondaryVertexV2BJetTags);
        		h_Jet5DeepCSV_b->Fill(GoodJetsDeepCSV_b[indexJets[5]]);//5]);
        		h_Jet5DeepCSV_bb->Fill(GoodJetsDeepCSV_bb[indexJets[5]]);//5]);
        		h_Jet5PT->Fill(GoodJetsPT[indexJets[5]]);//5]);
			h_Jet5Eta->Fill(GoodJetsEta[indexJets[5]]);//5]);
        		h_Jet5Phi->Fill(GoodJetsPhi[indexJets[5]]);//5]);
		}
		if(NumberOfGoodJets > 6){
        	//	h_Jet6combinedInclusiveSecondaryVertexV2BJetTags->Fill(GoodJetsCSVv2[indexJets[6]]);//6]);//Jet4combinedInclusiveSecondaryVertexV2BJetTags);
        		h_Jet6DeepCSV_b->Fill(GoodJetsDeepCSV_b[indexJets[6]]);//6]);
        		h_Jet6DeepCSV_bb->Fill(GoodJetsDeepCSV_bb[indexJets[6]]);//6]);
        		h_Jet6PT->Fill(GoodJetsPT[indexJets[6]]);//6]);
			if(GoodJetsPT.at(indexJets[6]) < cutMinJetPT && DEBUG==1  )cout<<endl<<i<<endl<<"GoodJetsPT[6]: "<<GoodJetsPT.at(indexJets.at(6))<<endl;
        		h_Jet6Eta->Fill(GoodJetsEta[indexJets[6]]);//6]);
			h_Jet6Phi->Fill(GoodJetsPhi[indexJets[6]]);//6]);

		}		
		if(DEBUG ==1)cout << "Finished looping over entries...." << endl;


		h_SumJetPT->Fill(SumJetPT);
        	h_SumJetEta->Fill(SumJetEta);

		h_DeepCSVSF->Fill(CSVsf);

		h_NumberOfJets->Fill(NumberOfGoodJets);//WhichGoodJets.size())
		h_NumberOfJets_SF->Fill(NumberOfGoodJets,evtSF);		







		h_NumberOfbTags->Fill(bTags);
		h_NumberOfbTags_SF->Fill(bTags,evtSF);



		h_NumberOfElectrons->Fill(WhichGoodElectrons.size());
  		h_NumberOfMuons->Fill(WhichGoodMuons.size());
		h_NumberOfLeptons->Fill(WhichGoodLeptons.size());

		h_additionalJetEventId->Fill(additionalJetEventId);


		if(is_ee==1){
			h_Jet0PT_ee->Fill(GoodJetsPT.at(indexJets[0]));
                        h_Jet0PT_SF_ee->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                        h_Jet1PT_ee->Fill(GoodJetsPT.at(indexJets[1]));
                        h_Jet1PT_SF_ee->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                        h_NumberOfJets_ee->Fill(NumberOfGoodJets);
                        h_NumberOfJets_SF_ee->Fill(NumberOfGoodJets,evtSF);
                        h_NumberOfbTags_ee->Fill(bTags);
                        h_NumberOfbTags_SF_ee->Fill(bTags,evtSF);
                        h_Lep0PT_ee->Fill(Lep0PT);
                        h_Lep0PT_SF_ee->Fill(Lep0PT,evtSF);
                        h_Lep1PT_ee->Fill(Lep1PT);
                        h_Lep1PT_SF_ee->Fill(Lep1PT,evtSF);
                        h_MET_ee->Fill(MET);
                        h_MET_SF_ee->Fill(MET,evtSF);
                        h_Jet0DeepCSV_b_ee->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));	
		}
		if(is_emu==1){
                        h_Jet0PT_emu->Fill(GoodJetsPT.at(indexJets[0]));
                        h_Jet0PT_SF_emu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                        h_Jet1PT_emu->Fill(GoodJetsPT.at(indexJets[1]));
                        h_Jet1PT_SF_emu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                        h_NumberOfJets_emu->Fill(NumberOfGoodJets);
                        h_NumberOfJets_SF_emu->Fill(NumberOfGoodJets,evtSF);
                        h_NumberOfbTags_emu->Fill(bTags);
                        h_NumberOfbTags_SF_emu->Fill(bTags,evtSF);
                        h_Lep0PT_emu->Fill(Lep0PT);
                        h_Lep0PT_SF_emu->Fill(Lep0PT,evtSF);
                        h_Lep1PT_emu->Fill(Lep1PT);
                        h_Lep1PT_SF_emu->Fill(Lep1PT,evtSF);
                        h_MET_emu->Fill(MET);
                        h_MET_SF_emu->Fill(MET,evtSF);
                        h_Jet0DeepCSV_b_emu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));                     
                }
		if(is_mumu==1){
                        h_Jet0PT_mumu->Fill(GoodJetsPT.at(indexJets[0]));
                        h_Jet0PT_SF_mumu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                        h_Jet1PT_mumu->Fill(GoodJetsPT.at(indexJets[1]));
                        h_Jet1PT_SF_mumu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                        h_NumberOfJets_mumu->Fill(NumberOfGoodJets);
                        h_NumberOfJets_SF_mumu->Fill(NumberOfGoodJets,evtSF);
                        h_NumberOfbTags_mumu->Fill(bTags);
                        h_NumberOfbTags_SF_mumu->Fill(bTags,evtSF);
                        h_Lep0PT_mumu->Fill(Lep0PT);
                        h_Lep0PT_SF_mumu->Fill(Lep0PT,evtSF);
                        h_Lep1PT_mumu->Fill(Lep1PT);
                        h_Lep1PT_SF_mumu->Fill(Lep1PT,evtSF);
                        h_MET_mumu->Fill(MET);
                        h_MET_SF_mumu->Fill(MET,evtSF);
                        h_Jet0DeepCSV_b_mumu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                }

		if(is_MC == 1 && TTbarMC == 1){
			if(is_ttlf == 1){

				if(NumberOfGoodJets >= 2 && bTags >=1 )h_bdt_score_ge2jge1b_SF_ttlf->Fill(bdt_score,evtSF);
		                if(NumberOfGoodJets == 3 && bTags == 2)h_bdt_score_3j2b_SF_ttlf->Fill(bdt_score,evtSF);
		                if(NumberOfGoodJets == 3 && bTags == 3)h_bdt_score_3j3b_SF_ttlf->Fill(bdt_score,evtSF);
		                if(NumberOfGoodJets >= 4 && bTags == 2)h_bdt_score_ge4j2b_SF_ttlf->Fill(bdt_score,evtSF);
		                if(NumberOfGoodJets >= 4 && bTags == 3)h_bdt_score_ge4j3b_SF_ttlf->Fill(bdt_score,evtSF);
		                if(NumberOfGoodJets >= 4 && bTags >= 4)h_bdt_score_ge4jge4b_SF_ttlf->Fill(bdt_score,evtSF);



				h_MET_ttlf->Fill(MET);
				h_Lep0Eta_ttlf->Fill(Lep0Eta);
				h_Lep1Eta_ttlf->Fill(Lep1Eta);
				h_DeepCSVSF_ttlf->Fill(CSVsf);
				h_NumberOfJets_ttlf->Fill(NumberOfGoodJets);
				h_NumberOfJets_SF_ttlf->Fill(NumberOfGoodJets,evtSF);
				h_NumberOfbTags_ttlf->Fill(bTags);
				h_NumberOfbTags_SF_ttlf->Fill(bTags,evtSF);
				h_Jet0Eta_ttlf->Fill(GoodJetsEta.at(indexJets[0]));
				h_Jet0Eta_SF_ttlf->Fill(GoodJetsEta.at(indexJets[0]),evtSF);
				h_Jet0PT_ttlf->Fill(GoodJetsPT.at(indexJets[0]));
				h_Jet0PT_SF_ttlf->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
				h_Jet1Eta_ttlf->Fill(GoodJetsEta.at(indexJets[1]));
				h_Jet1Eta_SF_ttlf->Fill(GoodJetsEta.at(indexJets[1]),evtSF);
				h_Jet1PT_ttlf->Fill(GoodJetsPT.at(indexJets[1]));
				h_Jet1PT_SF_ttlf->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
				h_Jet0DeepCSV_b_ttlf->Fill(GoodJetsDeepCSV_b.at(indexJets[0])); 
				h_Jet1DeepCSV_b_ttlf->Fill(GoodJetsDeepCSV_b.at(indexJets[1]));
				h_Jet0DeepCSV_b_SF_ttlf->Fill(GoodJetsDeepCSV_b.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_b_SF_ttlf->Fill(GoodJetsDeepCSV_b.at(indexJets[1]),evtSF);
				h_Jet0DeepCSV_bb_ttlf->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]));
				h_Jet1DeepCSV_bb_ttlf->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]));
				h_Jet0DeepCSV_bb_SF_ttlf->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                               h_Jet1DeepCSV_bb_SF_ttlf->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);


				h_Jet0DeepCSV_ttlf->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]));
               	 		h_Jet1DeepCSV_ttlf->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]));
                 		h_Jet0DeepCSV_SF_ttlf->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                		h_Jet1DeepCSV_SF_ttlf->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);

				h_Lep0PT_ttlf->Fill(Lep0PT);				
				h_Lep0PT_SF_ttlf->Fill(Lep0PT,evtSF);
				h_Lep1PT_ttlf->Fill(Lep1PT);
				h_Lep1PT_SF_ttlf->Fill(Lep1PT,evtSF);
				h_Lep0Eta_SF_ttlf->Fill(Lep0Eta,evtSF);
				h_Lep1Eta_SF_ttlf->Fill(Lep1Eta,evtSF);
				if(is_ee==1){
					h_Jet0PT_ttlf_ee->Fill(GoodJetsPT.at(indexJets[0]));
					h_Jet0PT_SF_ttlf_ee->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
					h_Jet1PT_ttlf_ee->Fill(GoodJetsPT.at(indexJets[1]));
                                	h_Jet1PT_SF_ttlf_ee->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
					h_NumberOfJets_ttlf_ee->Fill(NumberOfGoodJets);
                                	h_NumberOfJets_SF_ttlf_ee->Fill(NumberOfGoodJets,evtSF);
					h_NumberOfbTags_ttlf_ee->Fill(bTags);
                                	h_NumberOfbTags_SF_ttlf_ee->Fill(bTags,evtSF);
					h_Lep0PT_ttlf_ee->Fill(Lep0PT);
					h_Lep0PT_SF_ttlf_ee->Fill(Lep0PT,evtSF);
                                	h_Lep1PT_ttlf_ee->Fill(Lep1PT);
					h_Lep1PT_SF_ttlf_ee->Fill(Lep1PT,evtSF);
					h_MET_ttlf_ee->Fill(MET);
					h_MET_SF_ttlf_ee->Fill(MET,evtSF);
					h_Jet0DeepCSV_b_ttlf_ee->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
				}

				if(is_emu==1){
                                        h_Jet0PT_ttlf_emu->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttlf_emu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttlf_emu->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttlf_emu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttlf_emu->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttlf_emu->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttlf_emu->Fill(bTags);
                                        h_NumberOfbTags_SF_ttlf_emu->Fill(bTags,evtSF);
                                        h_Lep0PT_ttlf_emu->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttlf_emu->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttlf_emu->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttlf_emu->Fill(Lep1PT,evtSF);
                                        h_MET_ttlf_emu->Fill(MET);
                                        h_MET_SF_ttlf_emu->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttlf_emu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }


				if(is_mumu==1){
                                        h_Jet0PT_ttlf_mumu->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttlf_mumu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttlf_mumu->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttlf_mumu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttlf_mumu->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttlf_mumu->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttlf_mumu->Fill(bTags);
                                        h_NumberOfbTags_SF_ttlf_mumu->Fill(bTags,evtSF);
                                        h_Lep0PT_ttlf_mumu->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttlf_mumu->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttlf_mumu->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttlf_mumu->Fill(Lep1PT,evtSF);
                                        h_MET_ttlf_mumu->Fill(MET);
                                        h_MET_SF_ttlf_mumu->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttlf_mumu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }




			}else if(is_ttb == 1){


				if(NumberOfGoodJets >= 2 && bTags >=1 )h_bdt_score_ge2jge1b_SF_ttb->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets == 3 && bTags == 2)h_bdt_score_3j2b_SF_ttb->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets == 3 && bTags == 3)h_bdt_score_3j3b_SF_ttb->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags == 2)h_bdt_score_ge4j2b_SF_ttb->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags == 3)h_bdt_score_ge4j3b_SF_ttb->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags >= 4)h_bdt_score_ge4jge4b_SF_ttb->Fill(bdt_score,evtSF);



                                h_MET_ttb->Fill(MET);
				h_Lep0Eta_ttb->Fill(Lep0Eta);
				h_Lep1Eta_ttb->Fill(Lep1Eta);
				h_DeepCSVSF_ttb->Fill(CSVsf);
				h_NumberOfJets_ttb->Fill(NumberOfGoodJets);
                                h_NumberOfJets_SF_ttb->Fill(NumberOfGoodJets,evtSF);
				h_NumberOfbTags_ttb->Fill(bTags);
                                h_NumberOfbTags_SF_ttb->Fill(bTags,evtSF);
				h_Jet0Eta_ttb->Fill(GoodJetsEta.at(indexJets[0]));
				h_Jet0Eta_SF_ttb->Fill(GoodJetsEta.at(indexJets[0]),evtSF);
				h_Jet0PT_ttb->Fill(GoodJetsPT.at(indexJets[0]));
                                h_Jet1Eta_ttb->Fill(GoodJetsEta.at(indexJets[1]));
				h_Jet1Eta_SF_ttb->Fill(GoodJetsEta.at(indexJets[1]),evtSF);
				h_Jet1PT_ttb->Fill(GoodJetsPT.at(indexJets[1]));
                                h_Jet0DeepCSV_b_ttb->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
				h_Jet1DeepCSV_b_ttb->Fill(GoodJetsDeepCSV_b.at(indexJets[1]));                                
				h_Jet0DeepCSV_b_SF_ttb->Fill(GoodJetsDeepCSV_b.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_b_SF_ttb->Fill(GoodJetsDeepCSV_b.at(indexJets[1]),evtSF);
				h_Jet0DeepCSV_bb_ttb->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]));
				h_Jet1DeepCSV_bb_ttb->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]));
				h_Jet0DeepCSV_bb_SF_ttb->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_bb_SF_ttb->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);


				h_Jet0DeepCSV_ttb->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]));
                                h_Jet1DeepCSV_ttb->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]));
                                h_Jet0DeepCSV_SF_ttb->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_SF_ttb->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);



				h_Lep0PT_ttb->Fill(Lep0PT);
                                h_Lep1PT_ttb->Fill(Lep1PT);
				h_Lep0PT_SF_ttb->Fill(Lep0PT,evtSF);
				h_Lep0Eta_SF_ttb->Fill(Lep0Eta,evtSF);
				h_Lep1PT_SF_ttb->Fill(Lep1PT,evtSF);
                                h_Lep1Eta_SF_ttb->Fill(Lep1Eta,evtSF);

				if(is_ee==1){
                                        h_Jet0PT_ttb_ee->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttb_ee->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttb_ee->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttb_ee->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttb_ee->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttb_ee->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttb_ee->Fill(bTags);
                                        h_NumberOfbTags_SF_ttb_ee->Fill(bTags,evtSF);
                                        h_Lep0PT_ttb_ee->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttb_ee->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttb_ee->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttb_ee->Fill(Lep1PT,evtSF);
                                        h_MET_ttb_ee->Fill(MET);
                                        h_MET_SF_ttb_ee->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttb_ee->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }
				if(is_emu==1){
                                        h_Jet0PT_ttb_emu->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttb_emu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttb_emu->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttb_emu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttb_emu->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttb_emu->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttb_emu->Fill(bTags);
                                        h_NumberOfbTags_SF_ttb_emu->Fill(bTags,evtSF);
                                        h_Lep0PT_ttb_emu->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttb_emu->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttb_emu->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttb_emu->Fill(Lep1PT,evtSF);
                                        h_MET_ttb_emu->Fill(MET);
                                        h_MET_SF_ttb_emu->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttb_emu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }
				if(is_mumu==1){
                                        h_Jet0PT_ttb_mumu->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttb_mumu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttb_mumu->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttb_mumu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttb_mumu->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttb_mumu->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttb_mumu->Fill(bTags);
                                        h_NumberOfbTags_SF_ttb_mumu->Fill(bTags,evtSF);
                                        h_Lep0PT_ttb_mumu->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttb_mumu->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttb_mumu->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttb_mumu->Fill(Lep1PT,evtSF);
                                        h_MET_ttb_mumu->Fill(MET);
                                        h_MET_SF_ttb_mumu->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttb_mumu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }
				



                        }else if(is_ttbb == 1){
 

				if(NumberOfGoodJets >= 2 && bTags >=1 )h_bdt_score_ge2jge1b_SF_ttbb->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets == 3 && bTags == 2)h_bdt_score_3j2b_SF_ttbb->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets == 3 && bTags == 3)h_bdt_score_3j3b_SF_ttbb->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags == 2)h_bdt_score_ge4j2b_SF_ttbb->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags == 3)h_bdt_score_ge4j3b_SF_ttbb->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags >= 4)h_bdt_score_ge4jge4b_SF_ttbb->Fill(bdt_score,evtSF);





                               h_MET_ttbb->Fill(MET);
				h_Lep0Eta_ttbb->Fill(Lep0Eta);
				h_Lep1Eta_ttbb->Fill(Lep1Eta);
				h_DeepCSVSF_ttbb->Fill(CSVsf);
				h_NumberOfJets_ttbb->Fill(NumberOfGoodJets);
                                h_NumberOfJets_SF_ttbb->Fill(NumberOfGoodJets,evtSF);
				h_NumberOfbTags_ttbb->Fill(bTags);
                                h_NumberOfbTags_SF_ttbb->Fill(bTags,evtSF);
				h_Jet0Eta_ttbb->Fill(GoodJetsEta.at(indexJets[0]));
				h_Jet0Eta_SF_ttbb->Fill(GoodJetsEta.at(indexJets[0]),evtSF);
				h_Jet0PT_ttbb->Fill(GoodJetsPT.at(indexJets[0]));
                                h_Jet1Eta_ttbb->Fill(GoodJetsEta.at(indexJets[1]));
				h_Jet1Eta_SF_ttbb->Fill(GoodJetsEta.at(indexJets[1]),evtSF);
				h_Jet1PT_ttbb->Fill(GoodJetsPT.at(indexJets[1]));
                                h_Jet0DeepCSV_b_ttbb->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
				h_Jet1DeepCSV_b_ttbb->Fill(GoodJetsDeepCSV_b.at(indexJets[1]));
                                h_Jet0DeepCSV_b_SF_ttbb->Fill(GoodJetsDeepCSV_b.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_b_SF_ttbb->Fill(GoodJetsDeepCSV_b.at(indexJets[1]),evtSF);
				h_Jet0DeepCSV_bb_ttbb->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]));
				h_Jet1DeepCSV_bb_ttbb->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]));
				h_Jet0DeepCSV_bb_SF_ttbb->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_bb_SF_ttbb->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);
				
				h_Jet0DeepCSV_ttbb->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]));
                                h_Jet1DeepCSV_ttbb->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]));
                                h_Jet0DeepCSV_SF_ttbb->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_SF_ttbb->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);


				h_Lep0PT_ttbb->Fill(Lep0PT);
                                h_Lep0PT_SF_ttbb->Fill(Lep0PT,evtSF);
				h_Lep1PT_ttbb->Fill(Lep1PT);
				h_Lep1PT_SF_ttbb->Fill(Lep1PT,evtSF);
                                h_Lep0Eta_SF_ttbb->Fill(Lep0Eta,evtSF);
                                h_Lep1Eta_SF_ttbb->Fill(Lep1Eta,evtSF);

				if(is_ee==1){
                                        h_Jet0PT_ttbb_ee->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttbb_ee->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttbb_ee->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttbb_ee->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttbb_ee->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttbb_ee->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttbb_ee->Fill(bTags);
                                        h_NumberOfbTags_SF_ttbb_ee->Fill(bTags,evtSF);
                                        h_Lep0PT_ttbb_ee->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttbb_ee->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttbb_ee->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttbb_ee->Fill(Lep1PT,evtSF);
                                        h_MET_ttbb_ee->Fill(MET);
                                        h_MET_SF_ttbb_ee->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttbb_ee->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }
                                if(is_emu==1){
                                        h_Jet0PT_ttbb_emu->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttbb_emu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttbb_emu->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttbb_emu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttbb_emu->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttbb_emu->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttbb_emu->Fill(bTags);
                                        h_NumberOfbTags_SF_ttbb_emu->Fill(bTags,evtSF);
                                        h_Lep0PT_ttbb_emu->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttbb_emu->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttbb_emu->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttbb_emu->Fill(Lep1PT,evtSF);
                                        h_MET_ttbb_emu->Fill(MET);
                                        h_MET_SF_ttbb_emu->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttbb_emu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }
                                if(is_mumu==1){
                                        h_Jet0PT_ttbb_mumu->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttbb_mumu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttbb_mumu->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttbb_mumu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttbb_mumu->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttbb_mumu->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttbb_mumu->Fill(bTags);
                                        h_NumberOfbTags_SF_ttbb_mumu->Fill(bTags,evtSF);
                                        h_Lep0PT_ttbb_mumu->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttbb_mumu->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttbb_mumu->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttbb_mumu->Fill(Lep1PT,evtSF);
                                        h_MET_ttbb_mumu->Fill(MET);
                                        h_MET_SF_ttbb_mumu->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttbb_mumu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }
				

                        }else if(is_tt2b == 1){



				if(NumberOfGoodJets >= 2 && bTags >=1 )h_bdt_score_ge2jge1b_SF_tt2b->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets == 3 && bTags == 2)h_bdt_score_3j2b_SF_tt2b->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets == 3 && bTags == 3)h_bdt_score_3j3b_SF_tt2b->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags == 2)h_bdt_score_ge4j2b_SF_tt2b->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags == 3)h_bdt_score_ge4j3b_SF_tt2b->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags >= 4)h_bdt_score_ge4jge4b_SF_tt2b->Fill(bdt_score,evtSF);



                                h_MET_tt2b->Fill(MET);
				h_Lep0Eta_tt2b->Fill(Lep0Eta);
				h_Lep1Eta_tt2b->Fill(Lep1Eta);
				h_DeepCSVSF_tt2b->Fill(CSVsf);
				h_NumberOfJets_tt2b->Fill(NumberOfGoodJets);
                                h_NumberOfJets_SF_tt2b->Fill(NumberOfGoodJets,evtSF);
				h_NumberOfbTags_tt2b->Fill(bTags);
				h_NumberOfbTags_SF_tt2b->Fill(bTags,evtSF);
                                h_Jet0Eta_tt2b->Fill(GoodJetsEta.at(indexJets[0]));
				h_Jet0Eta_SF_tt2b->Fill(GoodJetsEta.at(indexJets[0]),evtSF);
				h_Jet0PT_tt2b->Fill(GoodJetsPT.at(indexJets[0]));
                                h_Jet1Eta_tt2b->Fill(GoodJetsEta.at(indexJets[1]));
				h_Jet1Eta_SF_tt2b->Fill(GoodJetsEta.at(indexJets[1]),evtSF);
				h_Jet1PT_tt2b->Fill(GoodJetsPT.at(indexJets[1]));
                                h_Jet0DeepCSV_b_tt2b->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                        	h_Jet1DeepCSV_b_tt2b->Fill(GoodJetsDeepCSV_b.at(indexJets[1]));        
				h_Jet0DeepCSV_b_SF_tt2b->Fill(GoodJetsDeepCSV_b.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_b_SF_tt2b->Fill(GoodJetsDeepCSV_b.at(indexJets[1]),evtSF);
				h_Jet0DeepCSV_bb_tt2b->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]));
                        	h_Jet1DeepCSV_bb_tt2b->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]));
				h_Jet0DeepCSV_bb_SF_tt2b->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_bb_SF_tt2b->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);
				

				h_Jet0DeepCSV_tt2b->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]));
                                h_Jet1DeepCSV_tt2b->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]));
                                h_Jet0DeepCSV_SF_tt2b->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_SF_tt2b->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);


				h_Lep0PT_tt2b->Fill(Lep0PT);
                                h_Lep1PT_tt2b->Fill(Lep1PT);
				h_Lep0PT_SF_tt2b->Fill(Lep0PT,evtSF);
                                h_Lep0Eta_SF_tt2b->Fill(Lep0Eta,evtSF);
                                h_Lep1PT_SF_tt2b->Fill(Lep1PT,evtSF);
                                h_Lep1Eta_SF_tt2b->Fill(Lep1Eta,evtSF);

				if(is_ee==1){
                                        h_Jet0PT_tt2b_ee->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_tt2b_ee->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_tt2b_ee->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_tt2b_ee->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_tt2b_ee->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_tt2b_ee->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_tt2b_ee->Fill(bTags);
                                        h_NumberOfbTags_SF_tt2b_ee->Fill(bTags,evtSF);
                                        h_Lep0PT_tt2b_ee->Fill(Lep0PT);
                                        h_Lep0PT_SF_tt2b_ee->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_tt2b_ee->Fill(Lep1PT);
                                        h_Lep1PT_SF_tt2b_ee->Fill(Lep1PT,evtSF);
                                        h_MET_tt2b_ee->Fill(MET);
                                        h_MET_SF_tt2b_ee->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_tt2b_ee->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }
                                if(is_emu==1){
                                        h_Jet0PT_tt2b_emu->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_tt2b_emu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_tt2b_emu->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_tt2b_emu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_tt2b_emu->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_tt2b_emu->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_tt2b_emu->Fill(bTags);
                                        h_NumberOfbTags_SF_tt2b_emu->Fill(bTags,evtSF);
                                        h_Lep0PT_tt2b_emu->Fill(Lep0PT);
                                        h_Lep0PT_SF_tt2b_emu->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_tt2b_emu->Fill(Lep1PT);
                                        h_Lep1PT_SF_tt2b_emu->Fill(Lep1PT,evtSF);
                                        h_MET_tt2b_emu->Fill(MET);
                                        h_MET_SF_tt2b_emu->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_tt2b_emu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }
                                if(is_mumu==1){
                                        h_Jet0PT_tt2b_mumu->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_tt2b_mumu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_tt2b_mumu->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_tt2b_mumu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_tt2b_mumu->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_tt2b_mumu->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_tt2b_mumu->Fill(bTags);
                                        h_NumberOfbTags_SF_tt2b_mumu->Fill(bTags,evtSF);
                                        h_Lep0PT_tt2b_mumu->Fill(Lep0PT);
                                        h_Lep0PT_SF_tt2b_mumu->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_tt2b_mumu->Fill(Lep1PT);
                                        h_Lep1PT_SF_tt2b_mumu->Fill(Lep1PT,evtSF);
                                        h_MET_tt2b_mumu->Fill(MET);
                                        h_MET_SF_tt2b_mumu->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_tt2b_mumu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }






			}else if(is_ttcc == 1){



				if(NumberOfGoodJets >= 2 && bTags >=1 )h_bdt_score_ge2jge1b_SF_ttcc->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets == 3 && bTags == 2)h_bdt_score_3j2b_SF_ttcc->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets == 3 && bTags == 3)h_bdt_score_3j3b_SF_ttcc->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags == 2)h_bdt_score_ge4j2b_SF_ttcc->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags == 3)h_bdt_score_ge4j3b_SF_ttcc->Fill(bdt_score,evtSF);
                                if(NumberOfGoodJets >= 4 && bTags >= 4)h_bdt_score_ge4jge4b_SF_ttcc->Fill(bdt_score,evtSF);







                                h_MET_ttcc->Fill(MET);
				h_Lep0Eta_ttcc->Fill(Lep0Eta);
				h_Lep1Eta_ttcc->Fill(Lep1Eta);
				h_DeepCSVSF_ttcc->Fill(CSVsf);
				h_NumberOfJets_ttcc->Fill(NumberOfGoodJets);
                                h_NumberOfJets_SF_ttcc->Fill(NumberOfGoodJets,evtSF);
				h_NumberOfbTags_ttcc->Fill(bTags);
				h_NumberOfbTags_SF_ttcc->Fill(bTags,evtSF);
                                h_Jet0Eta_ttcc->Fill(GoodJetsEta.at(indexJets[0]));
				h_Jet0Eta_SF_ttcc->Fill(GoodJetsEta.at(indexJets[0]),evtSF);
				h_Jet0PT_ttcc->Fill(GoodJetsPT.at(indexJets[0]));
                                h_Jet1Eta_ttcc->Fill(GoodJetsEta.at(indexJets[1]));
				h_Jet1Eta_SF_ttcc->Fill(GoodJetsEta.at(indexJets[1]),evtSF);
				h_Jet1PT_ttcc->Fill(GoodJetsPT.at(indexJets[1]));
                                h_Jet0DeepCSV_b_ttcc->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                        	h_Jet1DeepCSV_b_ttcc->Fill(GoodJetsDeepCSV_b.at(indexJets[1]));        
				h_Jet0DeepCSV_b_SF_ttcc->Fill(GoodJetsDeepCSV_b.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_b_SF_ttcc->Fill(GoodJetsDeepCSV_b.at(indexJets[1]),evtSF);
				h_Jet0DeepCSV_bb_ttcc->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]));
                        	h_Jet1DeepCSV_bb_ttcc->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]));
				h_Jet0DeepCSV_bb_SF_ttcc->Fill(GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_bb_SF_ttcc->Fill(GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);

				h_Jet0DeepCSV_ttcc->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]));
                                h_Jet1DeepCSV_ttcc->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]));
                                h_Jet0DeepCSV_SF_ttcc->Fill(GoodJetsDeepCSV_b.at(indexJets[0]) + GoodJetsDeepCSV_bb.at(indexJets[0]),evtSF);
                                h_Jet1DeepCSV_SF_ttcc->Fill(GoodJetsDeepCSV_b.at(indexJets[1]) + GoodJetsDeepCSV_bb.at(indexJets[1]),evtSF);


				h_Lep0PT_ttcc->Fill(Lep0PT);
                                h_Lep1PT_ttcc->Fill(Lep1PT);
				h_Lep0PT_SF_ttcc->Fill(Lep0PT,evtSF);
                                h_Lep0Eta_SF_ttcc->Fill(Lep0Eta,evtSF);
                                h_Lep1PT_SF_ttcc->Fill(Lep1PT,evtSF);
                                h_Lep1Eta_SF_ttcc->Fill(Lep1Eta,evtSF);

				if(is_ee==1){
                                        h_Jet0PT_ttcc_ee->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttcc_ee->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttcc_ee->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttcc_ee->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttcc_ee->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttcc_ee->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttcc_ee->Fill(bTags);
                                        h_NumberOfbTags_SF_ttcc_ee->Fill(bTags,evtSF);
                                        h_Lep0PT_ttcc_ee->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttcc_ee->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttcc_ee->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttcc_ee->Fill(Lep1PT,evtSF);
                                        h_MET_ttcc_ee->Fill(MET);
                                        h_MET_SF_ttcc_ee->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttcc_ee->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }
                                if(is_emu==1){
                                        h_Jet0PT_ttcc_emu->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttcc_emu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttcc_emu->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttcc_emu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttcc_emu->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttcc_emu->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttcc_emu->Fill(bTags);
                                        h_NumberOfbTags_SF_ttcc_emu->Fill(bTags,evtSF);
                                        h_Lep0PT_ttcc_emu->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttcc_emu->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttcc_emu->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttcc_emu->Fill(Lep1PT,evtSF);
                                        h_MET_ttcc_emu->Fill(MET);
                                        h_MET_SF_ttcc_emu->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttcc_emu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }
                                if(is_mumu==1){
                                        h_Jet0PT_ttcc_mumu->Fill(GoodJetsPT.at(indexJets[0]));
                                        h_Jet0PT_SF_ttcc_mumu->Fill(GoodJetsPT.at(indexJets[0]),evtSF);
                                        h_Jet1PT_ttcc_mumu->Fill(GoodJetsPT.at(indexJets[1]));
                                        h_Jet1PT_SF_ttcc_mumu->Fill(GoodJetsPT.at(indexJets[1]),evtSF);
                                        h_NumberOfJets_ttcc_mumu->Fill(NumberOfGoodJets);
                                        h_NumberOfJets_SF_ttcc_mumu->Fill(NumberOfGoodJets,evtSF);
                                        h_NumberOfbTags_ttcc_mumu->Fill(bTags);
                                        h_NumberOfbTags_SF_ttcc_mumu->Fill(bTags,evtSF);
                                        h_Lep0PT_ttcc_mumu->Fill(Lep0PT);
                                        h_Lep0PT_SF_ttcc_mumu->Fill(Lep0PT,evtSF);
                                        h_Lep1PT_ttcc_mumu->Fill(Lep1PT);
                                        h_Lep1PT_SF_ttcc_mumu->Fill(Lep1PT,evtSF);
                                        h_MET_ttcc_mumu->Fill(MET);
                                        h_MET_SF_ttcc_mumu->Fill(MET,evtSF);
                                        h_Jet0DeepCSV_b_ttcc_mumu->Fill(GoodJetsDeepCSV_b.at(indexJets[0]));
                                }











			}

		}






		//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
		// *** 3.6 Flush out all vectors for next iteration in event loop 

		bdt_score = 0.0;

		WhichGoodJets.clear();
		PotJetsPT.clear();
		GoodJetsPT.clear();
		GoodJetsEta.clear();
		GoodJetsPhi.clear();
		GoodJetsCSVv2.clear();
		GoodJetsDeepCSV_b.clear();
		GoodJetsDeepCSV_bb.clear();

		GoodLeptonsCharge.clear();
		GoodLeptonsE.clear();
		GoodLeptonsPT.clear();
		GoodLeptonsEta.clear();
		GoodLeptonsPhi.clear();
		GoodLeptonsIso.clear();
		

		GoodElectronsE.clear();
		GoodElectronsPT.clear();
		GoodElectronsEta.clear();

		GoodMuonsE.clear();
		GoodMuonsPT.clear();
		GoodMuonsEta.clear();
		GoodMuonsPhi.clear();
		GoodMuonsIso.clear();

		NumberOfGoodJets = 0;
		NumberOfJets = 0;
		NumberOfLeptons = 0;
		NumberOfMuons = 0;
		NumberOfElectrons = 0;
	
		additionalJetEventId = 0;

		flagEventEle = 0;
		flagEventMuon = 0;


		bTags = 0;

		MuonPT = 0;
		ElePT = 0;

		MET = 0;
		NumberOfPV = 0;

		JetPT = 0;
		SumJetPT = 0.0;
		JetEta = 0;
		SumJetEta = 0.0;
  	}// *** END EVENT LOOP
	

	//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
	// *** 4. Fill out histograms and write to Output_Histograms*.root file



/*
	fCSV.close();
	//failedCSV.close();
*/
	h_cutFlow->SetBinContent(1,cfTotal);
	h_cutFlow->GetXaxis()->SetBinLabel(1,"Total");
	h_cutFlow->SetBinContent(2,cfTrig);
        h_cutFlow->GetXaxis()->SetBinLabel(2,"Trigger");
	h_cutFlow->SetBinContent(3,cfTrigLep);
        h_cutFlow->GetXaxis()->SetBinLabel(3,"Lepton");
	h_cutFlow->SetBinContent(4,cfTrigLepJet);
	h_cutFlow->GetXaxis()->SetBinLabel(4,"Jets");
	h_cutFlow->SetBinContent(5,cfTrigLepJetBtag);
        h_cutFlow->GetXaxis()->SetBinLabel(5,"bTags");


	TString outSingleFile;
	outSingleFile = singleFile;
	int strLength = outSingleFile.Length();
  	outSingleFile.Remove(0,strLength-33);
	cout << outSingleFile << endl;


	TFile* fOut = new TFile("Output_Histograms_EventReadOut_MAIN_MarinoBDT_METFilter_PT15ge2j1b_May6_CONDOR_" + outSingleFile, "RECREATE");
	//TFile* fOut = new TFile("Output_Histograms_EventReadOut_MAIN_BDT_METFilter_PT15ge2j1b_Mar1"+var_DL+"_"+var_SF+"_"+var_TrigSFxx + "_CONDOR" + outSingleFile, "RECREATE");
  	//TFile* fOut = new TFile("Output_Histograms_EventReadOut_v22_"+output+"_"+var_DL+"_"+var_SF+"_"+var_TrigSFxx + "_CONDOR" + outSingleFile, "RECREATE");
  	//Output_Histograms_"+output+"_EventReadOut_v21_"+var_DL+"_"+var_SF+"_"+var_TrigSFxx + singleFile+".root", "RECREATE");
	//TLatex *   tex_JetAllPT = new TLatex(0.10,0.96,"#bf{CMS Preliminary}");
  	//plotPretty(h_JetAllPT,cJetAllPT,tex_JetAllPT,"PT","Events");
  	//plotStacker(hs_JetAllPT,h_JetAllPT,h_Jet0PT,h_Jet1PT,h_Jet2PT,h_Jet3PT,h_Jet4PT,h_Jet5PT);
	//gDirectory->GetFile();

	fOut->cd();


	h_cutFlow->Write();
	h_bdt_score_ge2jge1b->Write();
	h_bdt_score_3j2b->Write();
	h_bdt_score_3j3b->Write();
	h_bdt_score_ge4j2b->Write();
	h_bdt_score_ge4j3b->Write();
	h_bdt_score_ge4jge4b->Write();

	h_bdt_score_ge2jge1b_SF->Write();
        h_bdt_score_3j2b_SF->Write();
        h_bdt_score_3j3b_SF->Write();
        h_bdt_score_ge4j2b_SF->Write();
        h_bdt_score_ge4j3b_SF->Write();
        h_bdt_score_ge4jge4b_SF->Write();

	h_bdt_score_ge2jge1b_SF_ttlf->Write();
        h_bdt_score_3j2b_SF_ttlf->Write();
        h_bdt_score_3j3b_SF_ttlf->Write();
        h_bdt_score_ge4j2b_SF_ttlf->Write();
        h_bdt_score_ge4j3b_SF_ttlf->Write();
        h_bdt_score_ge4jge4b_SF_ttlf->Write();

	h_bdt_score_ge2jge1b_SF_ttb->Write();
        h_bdt_score_3j2b_SF_ttb->Write();
        h_bdt_score_3j3b_SF_ttb->Write();
        h_bdt_score_ge4j2b_SF_ttb->Write();
        h_bdt_score_ge4j3b_SF_ttb->Write();
        h_bdt_score_ge4jge4b_SF_ttb->Write();

	h_bdt_score_ge2jge1b_SF_ttbb->Write();
        h_bdt_score_3j2b_SF_ttbb->Write();
        h_bdt_score_3j3b_SF_ttbb->Write();
        h_bdt_score_ge4j2b_SF_ttbb->Write();
        h_bdt_score_ge4j3b_SF_ttbb->Write();
        h_bdt_score_ge4jge4b_SF_ttbb->Write();

	h_bdt_score_ge2jge1b_SF_tt2b->Write();
        h_bdt_score_3j2b_SF_tt2b->Write();
        h_bdt_score_3j3b_SF_tt2b->Write();
        h_bdt_score_ge4j2b_SF_tt2b->Write();
        h_bdt_score_ge4j3b_SF_tt2b->Write();
        h_bdt_score_ge4jge4b_SF_tt2b->Write();

	h_bdt_score_ge2jge1b_SF_ttcc->Write();
        h_bdt_score_3j2b_SF_ttcc->Write();
        h_bdt_score_3j3b_SF_ttcc->Write();
        h_bdt_score_ge4j2b_SF_ttcc->Write();
        h_bdt_score_ge4j3b_SF_ttcc->Write();
        h_bdt_score_ge4jge4b_SF_ttcc->Write();


	h_wgt_generator->Fill(wgt_generator_);




	h_sfEvent->Write();
	h_sfMuonID->Write();
	h_sfMuonIso->Write();
	h_sfEleID->Write();
	h_sfEleIso->Write();
	h_sfCSV->Write();
	h_sfTrig->Write();

	h_Jet0PT_ee->Write();
        h_Jet0PT_SF_ee->Write();
    	h_Jet1PT_ee->Write();
   	h_Jet1PT_SF_ee->Write();
   	h_NumberOfJets_ee->Write();
   	h_NumberOfJets_SF_ee->Write();
   	h_NumberOfbTags_ee->Write();
     	h_NumberOfbTags_SF_ee->Write();
	h_Lep0PT_ee->Write();
     	h_Lep0PT_SF_ee->Write();
  	h_Lep1PT_ee->Write();
     	h_Lep1PT_SF_ee->Write();
     	h_MET_ee->Write();
   	h_MET_SF_ee->Write();
    	h_Jet0DeepCSV_b_ee->Write();

	h_Jet0PT_emu->Write();
	h_Jet0PT_SF_emu->Write();
        h_Jet1PT_emu->Write();
        h_Jet1PT_SF_emu->Write();
        h_NumberOfJets_emu->Write();
        h_NumberOfJets_SF_emu->Write();
        h_NumberOfbTags_emu->Write();
        h_NumberOfbTags_SF_emu->Write();
        h_Lep0PT_emu->Write();
        h_Lep0PT_SF_emu->Write();
        h_Lep1PT_emu->Write();
        h_Lep1PT_SF_emu->Write();
        h_MET_emu->Write();
        h_MET_SF_emu->Write();
        h_Jet0DeepCSV_b_emu->Write();

	h_Jet0PT_mumu->Write();
        h_Jet0PT_SF_mumu->Write();
        h_Jet1PT_mumu->Write();
        h_Jet1PT_SF_mumu->Write();
        h_NumberOfJets_mumu->Write();
        h_NumberOfJets_SF_mumu->Write();
        h_NumberOfbTags_mumu->Write();
        h_NumberOfbTags_SF_mumu->Write();
        h_Lep0PT_mumu->Write();
        h_Lep0PT_SF_mumu->Write();
        h_Lep1PT_mumu->Write();
        h_Lep1PT_SF_mumu->Write();
        h_MET_mumu->Write();
        h_MET_SF_mumu->Write();
        h_Jet0DeepCSV_b_mumu->Write();

//00000000000000000

	h_Jet0PT_ttlf_ee->Write();
        h_Jet0PT_SF_ttlf_ee->Write();
        h_Jet1PT_ttlf_ee->Write();
        h_Jet1PT_SF_ttlf_ee->Write();
        h_NumberOfJets_ttlf_ee->Write();
        h_NumberOfJets_SF_ttlf_ee->Write();
        h_NumberOfbTags_ttlf_ee->Write();
        h_NumberOfbTags_SF_ttlf_ee->Write();
        h_Lep0PT_ttlf_ee->Write();
        h_Lep0PT_SF_ttlf_ee->Write();
        h_Lep1PT_ttlf_ee->Write();
        h_Lep1PT_SF_ttlf_ee->Write();
        h_MET_ttlf_ee->Write();
        h_MET_SF_ttlf_ee->Write();
        h_Jet0DeepCSV_b_ttlf_ee->Write();

        h_Jet0PT_ttlf_emu->Write();
        h_Jet0PT_SF_ttlf_emu->Write();
        h_Jet1PT_ttlf_emu->Write();
        h_Jet1PT_SF_ttlf_emu->Write();
        h_NumberOfJets_ttlf_emu->Write();
        h_NumberOfJets_SF_ttlf_emu->Write();
        h_NumberOfbTags_ttlf_emu->Write();
        h_NumberOfbTags_SF_ttlf_emu->Write();
        h_Lep0PT_ttlf_emu->Write();
        h_Lep0PT_SF_ttlf_emu->Write();
        h_Lep1PT_ttlf_emu->Write();
        h_Lep1PT_SF_ttlf_emu->Write();
        h_MET_ttlf_emu->Write();
        h_MET_SF_ttlf_emu->Write();
        h_Jet0DeepCSV_b_ttlf_emu->Write();
        
        h_Jet0PT_ttlf_mumu->Write();
        h_Jet0PT_SF_ttlf_mumu->Write();
        h_Jet1PT_ttlf_mumu->Write();
        h_Jet1PT_SF_ttlf_mumu->Write();
        h_NumberOfJets_ttlf_mumu->Write();
        h_NumberOfJets_SF_ttlf_mumu->Write();
        h_NumberOfbTags_ttlf_mumu->Write();
        h_NumberOfbTags_SF_ttlf_mumu->Write();
        h_Lep0PT_ttlf_mumu->Write();
        h_Lep0PT_SF_ttlf_mumu->Write();
        h_Lep1PT_ttlf_mumu->Write();
        h_Lep1PT_SF_ttlf_mumu->Write();
        h_MET_ttlf_mumu->Write();
        h_MET_SF_ttlf_mumu->Write();
        h_Jet0DeepCSV_b_ttlf_mumu->Write();

//000000000000000000


	h_Jet0PT_ttb_ee->Write();
        h_Jet0PT_SF_ttb_ee->Write();
        h_Jet1PT_ttb_ee->Write();
        h_Jet1PT_SF_ttb_ee->Write();
        h_NumberOfJets_ttb_ee->Write();
        h_NumberOfJets_SF_ttb_ee->Write();
        h_NumberOfbTags_ttb_ee->Write();
        h_NumberOfbTags_SF_ttb_ee->Write();
        h_Lep0PT_ttb_ee->Write();
        h_Lep0PT_SF_ttb_ee->Write();
        h_Lep1PT_ttb_ee->Write();
        h_Lep1PT_SF_ttb_ee->Write();
        h_MET_ttb_ee->Write();
        h_MET_SF_ttb_ee->Write();
        h_Jet0DeepCSV_b_ttb_ee->Write();

        h_Jet0PT_ttb_emu->Write();
        h_Jet0PT_SF_ttb_emu->Write();
        h_Jet1PT_ttb_emu->Write();
        h_Jet1PT_SF_ttb_emu->Write();
        h_NumberOfJets_ttb_emu->Write();
        h_NumberOfJets_SF_ttb_emu->Write();
        h_NumberOfbTags_ttb_emu->Write();
        h_NumberOfbTags_SF_ttb_emu->Write();
        h_Lep0PT_ttb_emu->Write();
        h_Lep0PT_SF_ttb_emu->Write();
        h_Lep1PT_ttb_emu->Write();
        h_Lep1PT_SF_ttb_emu->Write();
        h_MET_ttb_emu->Write();
        h_MET_SF_ttb_emu->Write();
        h_Jet0DeepCSV_b_ttb_emu->Write();

        h_Jet0PT_ttb_mumu->Write();
        h_Jet0PT_SF_ttb_mumu->Write();
        h_Jet1PT_ttb_mumu->Write();
        h_Jet1PT_SF_ttb_mumu->Write();
        h_NumberOfJets_ttb_mumu->Write();
        h_NumberOfJets_SF_ttb_mumu->Write();
        h_NumberOfbTags_ttb_mumu->Write();
        h_NumberOfbTags_SF_ttb_mumu->Write();
        h_Lep0PT_ttb_mumu->Write();
        h_Lep0PT_SF_ttb_mumu->Write();
        h_Lep1PT_ttb_mumu->Write();
        h_Lep1PT_SF_ttb_mumu->Write();
        h_MET_ttb_mumu->Write();
        h_MET_SF_ttb_mumu->Write();
        h_Jet0DeepCSV_b_ttb_mumu->Write();
//0000000000000000000000
 	h_Jet0PT_ttbb_ee->Write();
        h_Jet0PT_SF_ttbb_ee->Write();
        h_Jet1PT_ttbb_ee->Write();
        h_Jet1PT_SF_ttbb_ee->Write();
        h_NumberOfJets_ttbb_ee->Write();
        h_NumberOfJets_SF_ttbb_ee->Write();
        h_NumberOfbTags_ttbb_ee->Write();
        h_NumberOfbTags_SF_ttbb_ee->Write();
        h_Lep0PT_ttbb_ee->Write();
        h_Lep0PT_SF_ttbb_ee->Write();
        h_Lep1PT_ttbb_ee->Write();
        h_Lep1PT_SF_ttbb_ee->Write();
        h_MET_ttbb_ee->Write();
        h_MET_SF_ttbb_ee->Write();
        h_Jet0DeepCSV_b_ttbb_ee->Write();

        h_Jet0PT_ttbb_emu->Write();
        h_Jet0PT_SF_ttbb_emu->Write();
        h_Jet1PT_ttbb_emu->Write();
        h_Jet1PT_SF_ttbb_emu->Write();
        h_NumberOfJets_ttbb_emu->Write();
        h_NumberOfJets_SF_ttbb_emu->Write();
        h_NumberOfbTags_ttbb_emu->Write();
        h_NumberOfbTags_SF_ttbb_emu->Write();
        h_Lep0PT_ttbb_emu->Write();
        h_Lep0PT_SF_ttbb_emu->Write();
        h_Lep1PT_ttbb_emu->Write();
        h_Lep1PT_SF_ttbb_emu->Write();
        h_MET_ttbb_emu->Write();
        h_MET_SF_ttbb_emu->Write();
        h_Jet0DeepCSV_b_ttbb_emu->Write();

        h_Jet0PT_ttbb_mumu->Write();
        h_Jet0PT_SF_ttbb_mumu->Write();
        h_Jet1PT_ttbb_mumu->Write();
        h_Jet1PT_SF_ttbb_mumu->Write();
        h_NumberOfJets_ttbb_mumu->Write();
        h_NumberOfJets_SF_ttbb_mumu->Write();
        h_NumberOfbTags_ttbb_mumu->Write();
        h_NumberOfbTags_SF_ttbb_mumu->Write();
        h_Lep0PT_ttbb_mumu->Write();
        h_Lep0PT_SF_ttbb_mumu->Write();
        h_Lep1PT_ttbb_mumu->Write();
        h_Lep1PT_SF_ttbb_mumu->Write();
        h_MET_ttbb_mumu->Write();
        h_MET_SF_ttbb_mumu->Write();
        h_Jet0DeepCSV_b_ttbb_mumu->Write();


//00000000000000000000
 h_Jet0PT_tt2b_ee->Write();
        h_Jet0PT_SF_tt2b_ee->Write();
        h_Jet1PT_tt2b_ee->Write();
        h_Jet1PT_SF_tt2b_ee->Write();
        h_NumberOfJets_tt2b_ee->Write();
        h_NumberOfJets_SF_tt2b_ee->Write();
        h_NumberOfbTags_tt2b_ee->Write();
        h_NumberOfbTags_SF_tt2b_ee->Write();
        h_Lep0PT_tt2b_ee->Write();
        h_Lep0PT_SF_tt2b_ee->Write();
        h_Lep1PT_tt2b_ee->Write();
        h_Lep1PT_SF_tt2b_ee->Write();
        h_MET_tt2b_ee->Write();
        h_MET_SF_tt2b_ee->Write();
        h_Jet0DeepCSV_b_tt2b_ee->Write();

        h_Jet0PT_tt2b_emu->Write();
        h_Jet0PT_SF_tt2b_emu->Write();
        h_Jet1PT_tt2b_emu->Write();
        h_Jet1PT_SF_tt2b_emu->Write();
        h_NumberOfJets_tt2b_emu->Write();
        h_NumberOfJets_SF_tt2b_emu->Write();
        h_NumberOfbTags_tt2b_emu->Write();
        h_NumberOfbTags_SF_tt2b_emu->Write();
        h_Lep0PT_tt2b_emu->Write();
        h_Lep0PT_SF_tt2b_emu->Write();
        h_Lep1PT_tt2b_emu->Write();
        h_Lep1PT_SF_tt2b_emu->Write();
        h_MET_tt2b_emu->Write();
        h_MET_SF_tt2b_emu->Write();
        h_Jet0DeepCSV_b_tt2b_emu->Write();

        h_Jet0PT_tt2b_mumu->Write();
        h_Jet0PT_SF_tt2b_mumu->Write();
        h_Jet1PT_tt2b_mumu->Write();
        h_Jet1PT_SF_tt2b_mumu->Write();
        h_NumberOfJets_tt2b_mumu->Write();
        h_NumberOfJets_SF_tt2b_mumu->Write();
        h_NumberOfbTags_tt2b_mumu->Write();
        h_NumberOfbTags_SF_tt2b_mumu->Write();
        h_Lep0PT_tt2b_mumu->Write();
        h_Lep0PT_SF_tt2b_mumu->Write();
        h_Lep1PT_tt2b_mumu->Write();
        h_Lep1PT_SF_tt2b_mumu->Write();
        h_MET_tt2b_mumu->Write();
        h_MET_SF_tt2b_mumu->Write();
        h_Jet0DeepCSV_b_tt2b_mumu->Write();


//00000000000000000000

 h_Jet0PT_ttcc_ee->Write();
        h_Jet0PT_SF_ttcc_ee->Write();
        h_Jet1PT_ttcc_ee->Write();
        h_Jet1PT_SF_ttcc_ee->Write();
        h_NumberOfJets_ttcc_ee->Write();
        h_NumberOfJets_SF_ttcc_ee->Write();
        h_NumberOfbTags_ttcc_ee->Write();
        h_NumberOfbTags_SF_ttcc_ee->Write();
        h_Lep0PT_ttcc_ee->Write();
        h_Lep0PT_SF_ttcc_ee->Write();
        h_Lep1PT_ttcc_ee->Write();
        h_Lep1PT_SF_ttcc_ee->Write();
        h_MET_ttcc_ee->Write();
        h_MET_SF_ttcc_ee->Write();
        h_Jet0DeepCSV_b_ttcc_ee->Write();

        h_Jet0PT_ttcc_emu->Write();
        h_Jet0PT_SF_ttcc_emu->Write();
        h_Jet1PT_ttcc_emu->Write();
        h_Jet1PT_SF_ttcc_emu->Write();
        h_NumberOfJets_ttcc_emu->Write();
        h_NumberOfJets_SF_ttcc_emu->Write();
        h_NumberOfbTags_ttcc_emu->Write();
        h_NumberOfbTags_SF_ttcc_emu->Write();
        h_Lep0PT_ttcc_emu->Write();
        h_Lep0PT_SF_ttcc_emu->Write();
        h_Lep1PT_ttcc_emu->Write();
        h_Lep1PT_SF_ttcc_emu->Write();
        h_MET_ttcc_emu->Write();
        h_MET_SF_ttcc_emu->Write();
        h_Jet0DeepCSV_b_ttcc_emu->Write();

        h_Jet0PT_ttcc_mumu->Write();
        h_Jet0PT_SF_ttcc_mumu->Write();
        h_Jet1PT_ttcc_mumu->Write();
        h_Jet1PT_SF_ttcc_mumu->Write();
        h_NumberOfJets_ttcc_mumu->Write();
        h_NumberOfJets_SF_ttcc_mumu->Write();
        h_NumberOfbTags_ttcc_mumu->Write();
        h_NumberOfbTags_SF_ttcc_mumu->Write();
        h_Lep0PT_ttcc_mumu->Write();
        h_Lep0PT_SF_ttcc_mumu->Write();
        h_Lep1PT_ttcc_mumu->Write();
        h_Lep1PT_SF_ttcc_mumu->Write();
        h_MET_ttcc_mumu->Write();
        h_MET_SF_ttcc_mumu->Write();
        h_Jet0DeepCSV_b_ttcc_mumu->Write();


//00000000000000000000
	h_NumberOfJets_ttlf->Write();
        h_NumberOfJets_SF_ttlf->Write();
	h_NumberOfbTags_ttlf->Write();
        h_NumberOfbTags_SF_ttlf->Write();
	h_Lep0PT_SF_ttlf->Write();
	h_Lep0Eta_SF_ttlf->Write();
	h_Lep1PT_SF_ttlf->Write();
	h_Lep1Eta_SF_ttlf->Write();	
	h_Jet0PT_ttlf->Write();
	h_Jet0PT_SF_ttlf->Write();
        h_Jet1PT_ttlf->Write();
	h_Jet1PT_SF_ttlf->Write();
        h_Jet0DeepCSV_b_ttlf->Write();
        h_Jet0DeepCSV_bb_ttlf->Write();
	h_Jet1DeepCSV_b_ttlf->Write();
        h_Jet1DeepCSV_bb_ttlf->Write();

	h_Jet0DeepCSV_b_SF_ttlf->Write();
        h_Jet0DeepCSV_bb_SF_ttlf->Write();
        h_Jet1DeepCSV_b_SF_ttlf->Write();
        h_Jet1DeepCSV_bb_SF_ttlf->Write();



	h_NumberOfJets_ttb->Write();
        h_NumberOfJets_SF_ttb->Write();
	h_NumberOfbTags_ttb->Write();
        h_NumberOfbTags_SF_ttb->Write();
	h_Lep0PT_SF_ttb->Write();
        h_Lep0Eta_SF_ttb->Write();
        h_Lep1PT_SF_ttb->Write();
        h_Lep1Eta_SF_ttb->Write();
	h_Jet0PT_ttb->Write();
        h_Jet1PT_ttb->Write();
	h_Jet0PT_SF_ttb->Write();
        h_Jet1PT_SF_ttb->Write();
        h_Jet0DeepCSV_b_ttb->Write();
        h_Jet0DeepCSV_bb_ttb->Write();
	h_Jet1DeepCSV_b_ttb->Write();
        h_Jet1DeepCSV_bb_ttb->Write();

	h_Jet0DeepCSV_b_SF_ttb->Write();
        h_Jet0DeepCSV_bb_SF_ttb->Write();
        h_Jet1DeepCSV_b_SF_ttb->Write();
        h_Jet1DeepCSV_bb_SF_ttb->Write();



	h_NumberOfJets_ttbb->Write();
        h_NumberOfJets_SF_ttbb->Write();
	h_NumberOfbTags_ttbb->Write();
        h_NumberOfbTags_SF_ttbb->Write();
        h_Lep0PT_SF_ttbb->Write();
        h_Lep0Eta_SF_ttbb->Write();
        h_Lep1PT_SF_ttbb->Write();
        h_Lep1Eta_SF_ttbb->Write();	
	h_Jet0PT_ttbb->Write();
        h_Jet1PT_ttbb->Write();
	h_Jet0PT_SF_ttbb->Write();
        h_Jet1PT_SF_ttbb->Write();
        h_Jet0DeepCSV_b_ttbb->Write();
        h_Jet0DeepCSV_bb_ttbb->Write();
        h_Jet1DeepCSV_b_ttbb->Write();
        h_Jet1DeepCSV_bb_ttbb->Write();

	h_Jet0DeepCSV_b_SF_ttbb->Write();
        h_Jet0DeepCSV_bb_SF_ttbb->Write();
        h_Jet1DeepCSV_b_SF_ttbb->Write();
        h_Jet1DeepCSV_bb_SF_ttbb->Write();



	h_NumberOfJets_tt2b->Write();
        h_NumberOfJets_SF_tt2b->Write();
	h_NumberOfbTags_tt2b->Write();
        h_NumberOfbTags_SF_tt2b->Write();
	h_Lep0PT_SF_tt2b->Write();
        h_Lep0Eta_SF_tt2b->Write();
        h_Lep1PT_SF_tt2b->Write();
        h_Lep1Eta_SF_tt2b->Write();
	h_Jet0PT_tt2b->Write();
        h_Jet1PT_tt2b->Write();
	h_Jet0PT_SF_tt2b->Write();
        h_Jet1PT_SF_tt2b->Write();
        h_Jet0DeepCSV_b_tt2b->Write();
        h_Jet0DeepCSV_bb_tt2b->Write();
        h_Jet1DeepCSV_b_tt2b->Write();
        h_Jet1DeepCSV_bb_tt2b->Write();

	h_Jet0DeepCSV_b_SF_tt2b->Write();
        h_Jet0DeepCSV_bb_SF_tt2b->Write();
        h_Jet1DeepCSV_b_SF_tt2b->Write();
        h_Jet1DeepCSV_bb_SF_tt2b->Write();



	h_NumberOfJets_ttcc->Write();
	h_NumberOfJets_SF_ttcc->Write();
        h_NumberOfbTags_ttcc->Write();
        h_NumberOfbTags_SF_ttcc->Write();
        h_Lep0PT_SF_ttcc->Write();
        h_Lep0Eta_SF_ttcc->Write();
        h_Lep1PT_SF_ttcc->Write();
        h_Lep1Eta_SF_ttcc->Write();
	h_Jet0PT_ttcc->Write();
        h_Jet1PT_ttcc->Write();
	h_Jet0PT_SF_ttcc->Write();
        h_Jet1PT_SF_ttcc->Write();
        h_Jet0DeepCSV_b_ttcc->Write();
        h_Jet0DeepCSV_bb_ttcc->Write();
        h_Jet1DeepCSV_b_ttcc->Write();
        h_Jet1DeepCSV_bb_ttcc->Write();

	h_Jet0DeepCSV_b_SF_ttcc->Write();
        h_Jet0DeepCSV_bb_SF_ttcc->Write();
        h_Jet1DeepCSV_b_SF_ttcc->Write();
        h_Jet1DeepCSV_bb_SF_ttcc->Write();



 	h_SumJetPT->Write();
 	h_SumJetEta->Write();
        h_Lep0PT->Write();
	h_Lep0PT_ttlf->Write();
	h_Lep0PT_ttb->Write();
	h_Lep0PT_ttbb->Write();
	h_Lep0PT_tt2b->Write();
	h_Lep0PT_ttcc->Write();

	h_Lep1PT->Write();
        h_Lep1PT_SF->Write();
	h_Lep1PT_ttlf->Write();
        h_Lep1PT_ttb->Write();
        h_Lep1PT_ttbb->Write();
        h_Lep1PT_tt2b->Write();
        h_Lep1PT_ttcc->Write();

	h_Lep1Eta->Write();


	h_Lep0PT_SF->Write();
	h_Lep0Eta->Write();
        h_Lep0Eta_ttlf->Write();
	h_Lep0Eta_ttb->Write();
	h_Lep0Eta_ttbb->Write();
	h_Lep0Eta_tt2b->Write();
	h_Lep0Eta_ttcc->Write();
	h_Lep0Eta_SF->Write();
        h_Lep0Phi->Write();
        h_Lep0Iso->Write(); 	
	h_Lep0IsoSF->Write();

	h_Muon0PT->Write();
  	h_Muon0PT_SF->Write();
	h_Muon0Eta->Write();
  	h_Muon0Eta_SF->Write();
	h_Muon0Phi->Write();
	h_Muon0Iso->Write();
    	h_Ele0PT->Write();
  	h_Ele0PT_SF->Write();
	h_Ele0Eta->Write();
  	h_Ele0Eta_SF->Write();
	h_Ele0Phi->Write();
  	h_Ele0Iso->Write();

	h_additionalJetEventId->Write();	


	h_DeepCSVSF->Write();
	h_DeepCSVSF_ttlf->Write();
	h_DeepCSVSF_ttb->Write();
	h_DeepCSVSF_ttbb->Write();
	h_DeepCSVSF_tt2b->Write();
	h_DeepCSVSF_ttcc->Write();




	h_NumberOfJets->Write();
 	h_NumberOfJets_SF->Write();
	h_NumberOfbTags->Write();
	h_NumberOfbTags_SF->Write();
	
	h_NumberOfLeptons->Write();
  	h_NumberOfMuons->Write();
  	h_NumberOfElectrons->Write();
  	h_NumberOfPV->Write();
	h_MET->Write();
	h_MET_ttlf->Write();
	h_MET_ttb->Write();
	h_MET_ttbb->Write();
	h_MET_tt2b->Write();
	h_MET_ttcc->Write();

  	h_Jet0combinedInclusiveSecondaryVertexV2BJetTags->Write(); 
  	h_Jet1combinedInclusiveSecondaryVertexV2BJetTags->Write();
  	h_Jet2combinedInclusiveSecondaryVertexV2BJetTags->Write();
  	h_Jet3combinedInclusiveSecondaryVertexV2BJetTags->Write();
  	h_Jet4combinedInclusiveSecondaryVertexV2BJetTags->Write();
  	h_Jet5combinedInclusiveSecondaryVertexV2BJetTags->Write();
  	h_Jet6combinedInclusiveSecondaryVertexV2BJetTags->Write();

  	h_Jet0PT->Write();
  	h_Jet0PT_SF->Write();
	h_Jet1PT->Write();
	h_Jet1PT_SF->Write();  
	h_Jet2PT->Write();
  	h_Jet3PT->Write();
  	h_Jet4PT->Write();
  	h_Jet5PT->Write();

   	h_Jet0Eta->Write();  
   	h_Jet0Eta_ttlf->Write();
	h_Jet0Eta_ttb->Write();	
	h_Jet0Eta_ttbb->Write();
	h_Jet0Eta_tt2b->Write();
	h_Jet0Eta_ttcc->Write();

	h_Jet0Eta_SF->Write();
	h_Jet0Eta_SF_ttlf->Write();
        h_Jet0Eta_SF_ttb->Write();
        h_Jet0Eta_SF_ttbb->Write();
        h_Jet0Eta_SF_tt2b->Write();
        h_Jet0Eta_SF_ttcc->Write();


	h_Jet1Eta->Write();
   	h_Jet1Eta_ttlf->Write();
	h_Jet1Eta_ttb->Write();
        h_Jet1Eta_ttbb->Write();
        h_Jet1Eta_tt2b->Write();
        h_Jet1Eta_ttcc->Write();
	h_Jet1Eta_SF->Write();
	h_Jet1Eta_SF_ttlf->Write();
        h_Jet1Eta_SF_ttb->Write();
        h_Jet1Eta_SF_ttbb->Write();
        h_Jet1Eta_SF_tt2b->Write();
        h_Jet1Eta_SF_ttcc->Write();

	h_Jet2Eta->Write();
   	h_Jet3Eta->Write();
   	h_Jet4Eta->Write();
   	h_Jet5Eta->Write();

   	h_Jet0Phi->Write();
   	h_Jet1Phi->Write();
   	h_Jet2Phi->Write();
   	h_Jet3Phi->Write();
   	h_Jet4Phi->Write();
   	h_Jet5Phi->Write();


	h_Jet0DeepCSV->Write();
        h_Jet1DeepCSV->Write();
        h_Jet0DeepCSV_SF->Write();
        h_Jet1DeepCSV_SF->Write();
	h_Jet0DeepCSV_ttlf->Write();
        h_Jet1DeepCSV_ttlf->Write();
     	h_Jet0DeepCSV_SF_ttlf->Write();
        h_Jet1DeepCSV_SF_ttlf->Write();

	h_Jet0DeepCSV_ttb->Write();
        h_Jet1DeepCSV_ttb->Write();
        h_Jet0DeepCSV_SF_ttb->Write();
        h_Jet1DeepCSV_SF_ttb->Write();

	h_Jet0DeepCSV_ttbb->Write();
        h_Jet1DeepCSV_ttbb->Write();
        h_Jet0DeepCSV_SF_ttbb->Write();
        h_Jet1DeepCSV_SF_ttbb->Write();

	h_Jet0DeepCSV_tt2b->Write();
        h_Jet1DeepCSV_tt2b->Write();
        h_Jet0DeepCSV_SF_tt2b->Write();
        h_Jet1DeepCSV_SF_tt2b->Write();


	h_Jet0DeepCSV_ttcc->Write();
        h_Jet1DeepCSV_ttcc->Write();
        h_Jet0DeepCSV_SF_ttcc->Write();
        h_Jet1DeepCSV_SF_ttcc->Write();



	h_Jet0DeepCSV_b->Write();
	h_Jet1DeepCSV_b->Write();
	h_Jet0DeepCSV_b_SF->Write();
        h_Jet1DeepCSV_b_SF->Write();

	h_Jet2DeepCSV_b->Write();
	h_Jet3DeepCSV_b->Write();

        h_Jet0SFDeepCSV->Write();
        h_Jet1SFDeepCSV->Write();
        h_Jet2SFDeepCSV->Write();
        h_Jet3SFDeepCSV->Write();

	h_Jet0DeepCSV_bb->Write();
	h_Jet1DeepCSV_bb->Write();
	h_Jet0DeepCSV_bb_SF->Write();
        h_Jet1DeepCSV_bb_SF->Write();

	h_Jet2DeepCSV_bb->Write();
	h_Jet3DeepCSV_bb->Write();

	
	//fOut->Close();
	//cout << "Closing sync file now... " << outputCSV << endl;
	cout << "Closing root file now... "<< fOut->GetName() <<endl;
	fOut->Close();
	return;


} 

int main(int argc, char **argv){
//      TString var_Sample = argv[1];
//      //      TString isMC = argv[2];
//      //      TString singleFile = argv[3];
	//int k = 0, var_Sample, isMC, singleFile){
	EventReadOut_MAIN_MarinoBDT_METFilter_PT15ge2j1b_May6(argv[1],argv[2],argv[3]);	
//EventReadOut_v22("TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8","true","root://cmseos.fnal.gov//store/group/lpctthrun2/UVA/ICHEP2018/MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_ICHEP18_postMCsync_v0_TTto2L2Nu/181107_184524/0000/yggdrasil_treeMaker_ttH_sync_11-06-18_v26_recipeTest_1-1.root");

//(TString var_Sample="", TString isMC="",TString singleFile="")

	return 0;
}







