#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TChain.h"
#include "TAttMarker.h"
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <sstream>

// Following headers help decode L1T ntuples
#include "../L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "../L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "../L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "../L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "../L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloTowerDataFormat.h"
#include "../L1Trigger/L1TNtuples/interface/L1AnalysisRecoJetDataFormat.h"
#include "../L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"
#include "../L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"

#define N_IETA_BINS 80
#define N_IPHI_BINS 72
#define N_ETT_BINS 50
#define ETT_RANGE 1000
#define N_INDIV_EVENTS 20

// Number of events that pass any cuts, for normalising later
int nPassing = 0;

// 1d formatter
void formatPlot1D(TH1D* plot1d, int colour){
    plot1d->GetXaxis()->SetTitleOffset(1.2);
    plot1d->GetYaxis()->SetTitleOffset(1.4);
    plot1d->SetMinimum(0.1);
    plot1d->SetLineColor(colour);
    plot1d->SetLineWidth(2);
    // plot1d->Scale(1. / (double) nPassing);
    plot1d->Draw("HIST");
    plot1d->SetStats(false);
}

// 2d formatter
void formatPlot2D(TH2D* plot2d){
    plot2d->GetXaxis()->SetTitleOffset(1.2);
    plot2d->GetYaxis()->SetTitleOffset(1.4);
    plot2d->Draw("colz");
    plot2d->SetStats(false);
}

int doHOvEStudy(TString inFileName,TString outFileName, TString puSelStr, TString singleJetStr,TString genStr, TString genMatchStr, TString subSystem){

    std::cout << "Running doHOvEStudy.cxx..." << std::endl;
    uint puSel = puSelStr.Atoi();
    bool rateAboveSingleJet = (singleJetStr.Atoi() > 0);
    bool genDecayInHCAL = (genStr.Atoi() > 0);
    bool genMatch = (genMatchStr.Atoi() > 0);
    std::vector<TString> towers;
    // std::vector<TString> jets;
    // std::vector<TString> gen;

    bool isTest = false;
    bool isData = !genDecayInHCAL;
    int maxTowerBarrel = 16;
    int maxTowerEndcap = 28;

    int minTowerForHOvE = 0;
    int maxTowerForHOvE = maxTowerBarrel;
    if (subSystem == "barrel"){
    minTowerForHOvE = 0;
    maxTowerForHOvE = maxTowerBarrel;
    }
    else if (subSystem == "endcap"){
	 minTowerForHOvE = maxTowerBarrel+1;
	 maxTowerForHOvE = maxTowerEndcap;
    }
    else if (subSystem == "both") {
	minTowerForHOvE = 0;
	maxTowerForHOvE = maxTowerEndcap;
    }
    else{
	std::cout << "bad SUBSYSTEM string - opts: barrel, endcap, both" << std::endl;
	return 0;
    }

    // Load files
    TChain * treeL1Gen;

    std::cout << "Loading up the TChain..." << std::endl;
    TChain * eventTree = new TChain("l1EventTree/L1EventTree");
    TChain * treeL1Towemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
    TChain * treeL1CaloTPemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
    TChain * treeL1emu = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
    TChain * treeL1reco = new TChain("l1RecoTree/RecoTree");

    towers.push_back(inFileName);
    if (!isData){
        treeL1Gen = new TChain("l1GeneratorTree/L1GenTree");
    }
    //     towers.push_back("/home/users/mcitron/triggerStudies/CMSSW_10_3_1/src/L1Ntuple.root");
    //     jets.push_back("/home/users/mcitron/triggerStudies/CMSSW_10_3_1/src/L1Ntuple.root");
    // else{
    //     towers.push_back("/hadoop/cms/store/user/mcitron/ZeroBias/zbMETa/181212_124635/0000/L1Ntuple_10.root");
    //     jets.push_back("/hadoop/cms/store/user/mcitron/ZeroBias/zbMETa/181212_124635/0000/L1Ntuple_10.root");
    //
    // }

    uint minFiles = towers.size();

    for(uint i = 0; i < minFiles; ++i) {
	eventTree->Add(towers[i]);
	treeL1Towemu->Add(towers[i]);
	treeL1CaloTPemu->Add(towers[i]);
	treeL1emu->Add(towers[i]);
	treeL1reco->Add(towers[i]);
	if (!isData){
	    treeL1Gen->Add(towers[i]);
	}
    }

    L1Analysis::L1AnalysisEventDataFormat           *event_ = new L1Analysis::L1AnalysisEventDataFormat();
    L1Analysis::L1AnalysisL1CaloTowerDataFormat     *l1Towemu_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
    L1Analysis::L1AnalysisCaloTPDataFormat     *l1CaloTPemu_ = new L1Analysis::L1AnalysisCaloTPDataFormat();
    L1Analysis::L1AnalysisL1UpgradeDataFormat       *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
    L1Analysis::L1AnalysisGeneratorDataFormat       *l1gen_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
    L1Analysis::L1AnalysisRecoVertexDataFormat       *l1reco_ = new L1Analysis::L1AnalysisRecoVertexDataFormat();
    // L1Analysis::L1AnalysisL1UpgradeDataFormat       *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();

    eventTree->SetBranchAddress("Event", &event_);
    treeL1Towemu->SetBranchAddress("L1CaloTower", &l1Towemu_);
    treeL1emu->SetBranchAddress("L1Upgrade", &l1emu_);
    treeL1reco->SetBranchAddress("Vertex", &l1reco_);
    treeL1CaloTPemu->SetBranchAddress("CaloTP", &l1CaloTPemu_);
    if (!isData){
	treeL1Gen->SetBranchAddress("Generator", &l1gen_);
    }

    // get number of entries
    Long64_t nentriesTowers;

    nentriesTowers = treeL1emu->GetEntries();
    // nentriesHCAL = treeL1HCALemu->GetEntries();
    // nentriesECAL = treeL1ECALemu->GetEntries();
    //

    // Initialise histograms
    TH1D* nTotHist = new TH1D("nTot",";;nPass",3,0,3);

    // Tower hists
    TH1D* hAllTowEt = new TH1D("towerEt", ";Tower E_{T}; # Towers", 40, -0.5, 39.5);
    TH1D* hAllTowEta = new TH1D("towerEta", "Towers vs iEta;iEta; # Towers", N_IETA_BINS, -40., 40.);
    TH1D* hTowTPETEta = new TH1D("towerTPETEta", "Average tower TP E_{T} vs. iEta; iEta; Average tower TP E_{T}", N_IETA_BINS, -40., 40.);

    TH1D* hTowPhiB = new TH1D("towerPhiB", "Towers vs iPhi in Barrel;iPhi; # Towers", N_IPHI_BINS, 0., 72.);
    TH1D* hTowPhiE = new TH1D("towerPhiE", "Towers vs iPhi in End cap;iPhi; # Towers", N_IPHI_BINS, 0., 72.);

    TH1D* hTowTPETphiB = new TH1D("towerTPETPhiB", "Average tower TP E_{T} vs. iPhi in Barrel;iPhi; Average tower TP E_{T}", N_IPHI_BINS, 0., 72.);
    TH1D* hTowTPETphiE = new TH1D("towerTPETPhiE", "Average tower TP E_{T} vs. iPhi in End cap;iPhi; Average tower TP E_{T}", N_IPHI_BINS, 0., 72.);
    // Jet hists
    TH1D * hJetEt = new TH1D("jetET",";ET;",100,0,1000);
    // TH1D * hMaxJetEt = new TH1D("maxJetEt",";ET;",100,0,2000);
    // TH1D * hMaxJetSeedEt = new TH1D("maxJetSeedEt",";ET;",100,0,1000);
    // TH1D * hMaxJetSeedEm = new TH1D("maxJetSeedEm",";ET;",100,0,1000);
    // TH1D * hMaxJetSeedHad = new TH1D("maxJetSeedHad",";ET;",100,0,1000);

    TH1D * hMaxTowerEt = new TH1D("maxTowerEt",";ET;",100,0,1000);
    TH1D * hMaxTowerEm = new TH1D("maxTowerEm",";ET;",100,0,1000);
    TH1D * hMaxTowerHad = new TH1D("maxTowerHad",";ET;",100,0,1000);

    std::vector<TString> ratioStrings = {"HOvE","HOvE3","HOvE9","H3OvE9","H3OvE3","H9OvE9"};
    std::map<const TString, TH2D*> maxJet2DHists;
    std::map<const TString, TH1D*> maxTower1DHists;
    std::map<const TString, TH2D*> maxTower2DHists;
    // std::map<const TString, TH1D*> maxJet1DHists;
    std::map<const TString, TH1D*> allTower1DHists;
    std::map<const TString, TH1D*> maxJet1DHists;

    std::map<const TString, double> hOvEThresholds;
    hOvEThresholds["None"] = -100;
    hOvEThresholds["5"] = 0.699;
    hOvEThresholds["10"] = 1;
    hOvEThresholds["20"] = 1.30;
    hOvEThresholds["100"] = 2;
    hOvEThresholds["Inf"] = 4;
    std::map<const TString, double> hadThresholds;
    hadThresholds["None"] = -100;
    hadThresholds["10"] = 10;
    hadThresholds["20"] = 20;
    hadThresholds["30"] = 30;
    hadThresholds["50"] = 50;
    TH2D * maxTowerEtaPhi = new TH2D("maxTowerEtaPhi",";eta;phi",200,-100.,100.,200,-100,100);
    TH2D * maxJetEtaPhi = new TH2D("maxJetEtaPhi",";eta;phi",200,-100.,100.,200,-100,100);

    std::map<const TString, double> hadVariablesSumJets;
    for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
	maxTower2DHists["maxTower"+*ratioIt] = new TH2D("maxTower"+*ratioIt,";ET;log(H/E)",100,0,500,110,-11,11);
	for (auto hOvEIt = hOvEThresholds.begin(); hOvEIt != hOvEThresholds.end(); hOvEIt++){
	    for (auto hadIt = hadThresholds.begin(); hadIt != hadThresholds.end(); hadIt++){
		maxJet1DHists["maxJetHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxJetHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";ET;",140,-12,12);
		maxTower1DHists["maxTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";ET;",140,-12,12);
		maxJet1DHists["maxJetPt"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxJetPt"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";ET;",1000,0,1000);
		maxTower1DHists["maxTowerHad"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxTowerHad"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";ET;",500,0,500);
		maxJet1DHists["maxJetDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxJetDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;",7,-0.5,6.5);
		maxTower1DHists["maxTowerDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxTowerDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;",7,-0.5,6.5);
		maxTower2DHists["maxTowerDepth2D"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH2D("maxTowerDepth2D"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;Energy",7,-0.5,6.5,100,0,50);
		maxJet1DHists["maxJetDepthRatio"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxJetDepthRatio"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;",100,0,1);
		maxJet2DHists["maxJetDepth1Vs3"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH2D("maxJetDepth1Vs3"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Energy 1; Energy 3",100,0,50,100,0,50);
		maxJet2DHists["maxJetDepth1Vs2"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH2D("maxJetDepth1Vs2"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Energy 1; Energy 2",100,0,50,100,0,50);
		maxJet2DHists["maxJetDepth2Vs3"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH2D("maxJetDepth2Vs3"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Energy 2; Energy 3",100,0,50,100,0,50);
		maxJet2DHists["maxJetDepth2D"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH2D("maxJetDepth2D"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;Energy",7,-0.5,6.5,100,0,50);
		maxJet1DHists["maxJetDenom"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxJetDenom"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;",100,0,100);
		maxJet1DHists["maxJetDenomIt"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxJetDenomIt"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;",100,0,10);
		maxJet1DHists["maxJetDenomHigh"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxJetDenomHigh"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;",100,0,100);
		maxJet1DHists["maxJetTotalHad"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxJetTotalHad"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;",100,0,100);
		maxJet1DHists["minJetEcal"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("minJetEcal"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;",100,0,100);
		maxJet1DHists["minJetLowDepth"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("minJetLowDepth"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;",100,0,100);
		maxJet1DHists["nJetPassSel"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("nJetPassSel"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";Depth;",10,0,10);
		hadVariablesSumJets[hOvEIt->first+"Had"+hadIt->first+*ratioIt] = 0;
	    }
	}
    }
    std::map<const TString, double> hadVariablesMaxTower;
    std::map<const TString, std::vector<double>> depthVariablesAllTowers;
    std::map<const TString, double> emVariablesMaxTower;


    TH1D* hPartVr = new TH1D("partVr", ";V_{R};", 100, 0., 1000.);
    TH1D* hPartVz = new TH1D("partVz", ";V_{Z};", 100, 0., 1000.);
    // Histogram arrays for storing individual event information

    // Main loop

    // Number of entries to run over
    Long64_t nEvents;
    if (isTest == true) nEvents = 1000;
    else nEvents = 1500000;

    nEvents = std::min(nentriesTowers,nEvents);//std::min(std::min(nentriesTowers,nentriesHCAL),nentriesECAL);
    for (Long64_t jentry = 0; jentry < nEvents; ++jentry) {
	nTotHist->Fill(0.5);
	std::map<const TString, std::vector<double> > hadVariablesAllJets;
	std::map<const TString, std::vector<std::vector<double>> > depthVariablesAllJets;
	std::map<const TString, std::vector<double> > emVariablesAllJets;
	// initialise some variables
	uint nHCALTPemu(0), nECALTPemu(0), nTowemu(0), nPart(0),nJetemu(0);
	int maxTowerIPhi(-1),maxTowerIEta(-1);
	double hcalTPEtEm(0), ecalTPEtEm(0), towEtemu(0), towHademu(0),  towEmemu(0),
	       maxTowerEt(0), maxTowerHad(0), maxTowerEm(0), maxTower3x3Em(0), 
	       maxTower3x3Had(0), maxTower9x9Em(0), maxTower9x9Had(0),ratioTower(0),ratioJet(0);

	for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
	    hadVariablesMaxTower[*ratioIt] = 0;
	    depthVariablesAllTowers[*ratioIt].clear();
	    emVariablesMaxTower[*ratioIt] = 0;
	}

	int hcalTPEtaEm(0), ecalTPEtaEm(0), towEtaemu(0);
	int hcalTPPhiEm(0), ecalTPPhiEm(0), towPhiemu(0);

	double l1MetEmu(0.);
	double l1MetHCALEmu(0.);
	double l1MetECALEmu(0.);

	double l1EttEmu(0.);
	double l1EttHCALEmu(0.);
	double l1EttECALEmu(0.);

	double l1MetPhiEmu(-1.);
	double l1MetPhiHCALEmu(-1.);
	double l1MetPhiECALEmu(-1.);

	int seedTowerIEta(-1);
	int seedTowerIPhi(-1);

	// run number
	int run(306042);

	//counter
	if( (jentry % 1000) == 0 ) std::cout << "Done " << jentry << " events of " << nEvents << std::endl;
	if( (jentry % 100) == 0 ) std::cout << "." << std::flush;

	eventTree->GetEntry(jentry);

	run = event_->run;
	int lumi = event_->lumi;
	int event = event_->event;

	treeL1Towemu->GetEntry(jentry);
	treeL1CaloTPemu->GetEntry(jentry);
	treeL1reco->GetEntry(jentry);
	if (!isData){
	    treeL1Gen->GetEntry(jentry);
	}
	treeL1emu->GetEntry(jentry);
	uint nVtxReco = l1reco_->nVtx;
	if (nVtxReco < puSel) continue;
	nTowemu = l1Towemu_->nTower;
	nECALTPemu = l1CaloTPemu_->nECALTP;
	nHCALTPemu = l1CaloTPemu_->nHCALTP;
	nJetemu = l1emu_->nJets;
	nTotHist->Fill(1.5);
	// // gen stuff
	bool decayInHCAL = false;
	std::vector<int> genMatchId;
	if (genDecayInHCAL){
	    nPart = l1gen_->nPart;
	    for(uint partIt = 0; partIt < nPart; ++partIt){
		if (abs(l1gen_->partId[partIt]) == 5 && (abs(l1gen_->partParent[partIt]) == 6000113 || abs(l1gen_->partParent[partIt]) == 9000007) || abs(l1gen_->partParent[partIt]) == 9000006) {
		    double hcalRadius = TMath::Sqrt(l1gen_->partVx[partIt]*l1gen_->partVx[partIt]+l1gen_->partVy[partIt]*l1gen_->partVy[partIt]);
		    if (minTowerForHOvE <= maxTowerBarrel && hcalRadius > 180 && hcalRadius < 270 && abs(l1gen_->partVz[partIt]) < 390) decayInHCAL = true;
		    if (maxTowerForHOvE >= maxTowerBarrel+1 && hcalRadius > 44 && hcalRadius < 270 && abs(l1gen_->partVz[partIt]) > 390 && abs(l1gen_->partVz[partIt]) < 568) decayInHCAL = true;
		    genMatchId.push_back(partIt);
		}
	    }
	}
	if (!decayInHCAL && genDecayInHCAL) continue;
	nTotHist->Fill(2.5);
	if (rateAboveSingleJet){
	    double maxPt = -1;
	    double secondPt = -1;
	    for (uint jetIt = 0; jetIt < nJetemu; ++jetIt){
		seedTowerIPhi = l1emu_->jetTowerIPhi[jetIt];
		// std::cout << seedTowerIPhi << std::endl;
		seedTowerIEta = l1emu_->jetTowerIEta[jetIt];
		// if (abs(seedTowerIEta) <= maxTowerForHOvE && abs(seedTowerIEta) >= minTowerForHOvE){
		if (abs(seedTowerIEta) <= maxTowerEndcap){
		    if (secondPt < l1emu_->jetEt[jetIt]) secondPt = l1emu_->jetEt[jetIt];
		    if (secondPt > maxPt) {
			secondPt = maxPt;
			maxPt = l1emu_->jetEt[jetIt];
		    }
		}
	    }
	    if (maxPt > 180 || secondPt > 150){
		continue;
	    }
	}


	for (uint jetIt = 0; jetIt < nJetemu; ++jetIt){
	    hJetEt->Fill(l1emu_->jetEt[jetIt]);
	    seedTowerIPhi = l1emu_->jetTowerIPhi[jetIt];
	    seedTowerIEta = l1emu_->jetTowerIEta[jetIt];
	    double seedTowerHad(0), seedTowerEm(0), seedTower3x3Em(0), 
		   seedTower3x3Had(0), seedTower9x9Em(0), seedTower9x9Had(0);
	    std::vector<double> depthVarHad(7,0);
	    // std::vector<double> timeVarHad(7,0);
	    std::vector<double> depthVar3x3Had(7,0); 
	    // std::vector<double> timeVar3x3Had(7,0); 
	    std::vector<double> depthVar9x9Had(7,0);
	    // std::vector<double> timeVar9x9Had(7,0);

	    for(uint HcalTPIt = 0; HcalTPIt < nHCALTPemu; ++HcalTPIt){
		towEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // use for HB HE restrictions                                 
		towPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt];  // hcalTPiphi[HcalTPIt]; // use for deltaR calculation
		towHademu = l1CaloTPemu_->hcalTPet[HcalTPIt]; // used for energy normalization in the energy depth plots    
		int nDepths = l1CaloTPemu_->hcalTPnDepths[HcalTPIt];
		if (abs(towEtaemu) >= minTowerForHOvE && abs(towEtaemu) <= maxTowerForHOvE){
		    if (nDepths == 0) continue;
		    std::vector<double> hcalTPDepth(7,0); 
		    // std::vector<double hcalTPTiming(7,0);
		    hcalTPDepth[0] = l1CaloTPemu_->hcalTPDepth1[HcalTPIt];
		    hcalTPDepth[1] = l1CaloTPemu_->hcalTPDepth2[HcalTPIt];
		    hcalTPDepth[2] = l1CaloTPemu_->hcalTPDepth3[HcalTPIt];
		    hcalTPDepth[3] = l1CaloTPemu_->hcalTPDepth4[HcalTPIt];
		    hcalTPDepth[4] = l1CaloTPemu_->hcalTPDepth5[HcalTPIt];
		    hcalTPDepth[5] = l1CaloTPemu_->hcalTPDepth6[HcalTPIt];
		    hcalTPDepth[6] = l1CaloTPemu_->hcalTPDepth7[HcalTPIt];

		    // hcalTPtiming[0] = l1CaloTPemu_->hcalTPtiming1[HcalTPIt];// - exp_time;
		    // hcalTPtiming[1] = l1CaloTPemu_->hcalTPtiming2[HcalTPIt];// - exp_time;
		    // hcalTPtiming[2] = l1CaloTPemu_->hcalTPtiming3[HcalTPIt];// - exp_time;
		    // hcalTPtiming[3] = l1CaloTPemu_->hcalTPtiming4[HcalTPIt];// - exp_time;
		    // hcalTPtiming[4] = l1CaloTPemu_->hcalTPtiming5[HcalTPIt];// - exp_time;
		    // hcalTPtiming[5] = l1CaloTPemu_->hcalTPtiming6[HcalTPIt];// - exp_time;
		    // hcalTPtiming[6] = l1CaloTPemu_->hcalTPtiming7[HcalTPIt];// - exp_time;
		    // std::cout << towHademu << " " <<  hcalTPDepth[0]+hcalTPDepth[1]+hcalTPDepth[2]+hcalTPDepth[3]+hcalTPDepth[4]+hcalTPDepth[5]+hcalTPDepth[6] << std::endl;
		    double weightedTime = 0;
		    if (towEtaemu == seedTowerIEta && towPhiemu == seedTowerIPhi){
			seedTowerHad = towHademu;
			for (int i = 0; i < nDepths-1; i++) depthVarHad[i] = hcalTPDepth[i];
			// for (int i = 0; i < nDepths-1; i++) timeVarHad[i] = hcalTPtiming[i];
		    }
		    for (int iSeedTowerIEta = -4; iSeedTowerIEta <= 4; ++iSeedTowerIEta){
			for (int iSeedTowerIPhi = -4; iSeedTowerIPhi <= 4; ++iSeedTowerIPhi){
			    int wrappedIPhi = seedTowerIPhi+iSeedTowerIPhi;
			    if (wrappedIPhi > 72) wrappedIPhi -= 72;
			    if (wrappedIPhi < 0) wrappedIPhi += 72;
			    if (towEtaemu == seedTowerIEta+iSeedTowerIEta && towPhiemu == wrappedIPhi){
				seedTower9x9Had += towHademu;
				for (int i = 0; i < nDepths-1; i++) depthVar9x9Had[i] += hcalTPDepth[i];
				if (abs(iSeedTowerIPhi) <= 1 && abs(iSeedTowerIEta) <= 1){
				    seedTower3x3Had += towHademu;
				    for (int i = 0; i < nDepths-1; i++) depthVar3x3Had[i] += hcalTPDepth[i];
				}
			    }
			}
		    }
		}
	    }
	    for(uint EcalTPIt = 0; EcalTPIt < nECALTPemu; ++EcalTPIt){
		towEtaemu = l1CaloTPemu_->ecalTPieta[EcalTPIt]; // use for HB HE restrictions                                 
		towPhiemu = l1CaloTPemu_->ecalTPCaliphi[EcalTPIt];  // hcalTPiphi[HcalTPIt]; // use for deltaR calculation
		towEmemu = l1CaloTPemu_->ecalTPet[EcalTPIt]; // used for energy normalization in the energy depth plots    
		if (abs(towEtaemu) >= minTowerForHOvE && abs(towEtaemu) <= maxTowerForHOvE){
		    if (towEtaemu == seedTowerIEta && towPhiemu == seedTowerIPhi){
			seedTowerEm = towEmemu;
		    }
		    for (int iSeedTowerIEta = -4; iSeedTowerIEta <= 4; ++iSeedTowerIEta){
			for (int iSeedTowerIPhi = -4; iSeedTowerIPhi <= 4; ++iSeedTowerIPhi){
			    int wrappedIPhi = seedTowerIPhi+iSeedTowerIPhi;
			    if (wrappedIPhi > 72) wrappedIPhi -= 72;
			    if (wrappedIPhi < 0) wrappedIPhi += 72;
			    if (towEtaemu == seedTowerIEta+iSeedTowerIEta && towPhiemu == wrappedIPhi){
				seedTower9x9Em += towEmemu;
				if (abs(iSeedTowerIPhi) <= 1 && abs(iSeedTowerIEta) <= 1){
				    seedTower3x3Em += towEmemu;
				}
			    }
			}
		    }
		}
	    }


	    depthVariablesAllJets["HOvE"].push_back(depthVarHad);
	    depthVariablesAllJets["HOvE3"].push_back(depthVarHad);
	    depthVariablesAllJets["HOvE9"].push_back(depthVarHad);
	    depthVariablesAllJets["H3OvE3"].push_back(depthVar3x3Had);
	    depthVariablesAllJets["H3OvE9"].push_back(depthVar3x3Had);
	    depthVariablesAllJets["H9OvE9"].push_back(depthVar9x9Had);

	    hadVariablesAllJets["HOvE"].push_back(seedTowerHad);
	    hadVariablesAllJets["HOvE3"].push_back(seedTowerHad);
	    hadVariablesAllJets["HOvE9"].push_back(seedTowerHad);
	    hadVariablesAllJets["H3OvE3"].push_back(seedTower3x3Had);
	    hadVariablesAllJets["H3OvE9"].push_back(seedTower3x3Had);
	    hadVariablesAllJets["H9OvE9"].push_back(seedTower9x9Had);

	    emVariablesAllJets["HOvE"].push_back(seedTowerEm);
	    emVariablesAllJets["HOvE3"].push_back(seedTower3x3Em);
	    emVariablesAllJets["HOvE9"].push_back(seedTower9x9Em);
	    emVariablesAllJets["H3OvE3"].push_back(seedTower3x3Em);
	    emVariablesAllJets["H3OvE9"].push_back(seedTower9x9Em);
	    emVariablesAllJets["H9OvE9"].push_back(seedTower9x9Em);
	    // std::cout << seedTower3x3Had << " " <<seedTower3x3Em << std::endl;
	}

	// if (decayInHCAL){
	//     hMaxJetEt->Fill(l1emu_->jetEt[0]);
	// }

	// Max tower stuff
	std::vector<double> depthVarTowerHad(7,0);
	std::vector<double> depthVarTower3x3Had(7,0); 
	std::vector<double> depthVarTower9x9Had(7,0);
	bool foundMax = false;
	for(uint HcalTPIt = 0; HcalTPIt < nHCALTPemu; ++HcalTPIt){
	    towEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // use for HB HE restrictions                                 
	    towPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt];  // hcalTPiphi[HcalTPIt]; // use for deltaR calculation
	    towHademu = l1CaloTPemu_->hcalTPet[HcalTPIt]; // used for energy normalization in the energy depth plots    
	    int nDepths = l1CaloTPemu_->hcalTPnDepths[HcalTPIt];
	    if (nDepths == 0) continue;
	    std::vector<double> hcalTPDepth(7,0); 
	    hcalTPDepth[0] = l1CaloTPemu_->hcalTPDepth1[HcalTPIt];
	    hcalTPDepth[1] = l1CaloTPemu_->hcalTPDepth2[HcalTPIt];
	    hcalTPDepth[2] = l1CaloTPemu_->hcalTPDepth3[HcalTPIt];
	    hcalTPDepth[3] = l1CaloTPemu_->hcalTPDepth4[HcalTPIt];
	    hcalTPDepth[4] = l1CaloTPemu_->hcalTPDepth5[HcalTPIt];
	    hcalTPDepth[5] = l1CaloTPemu_->hcalTPDepth6[HcalTPIt];
	    hcalTPDepth[6] = l1CaloTPemu_->hcalTPDepth7[HcalTPIt];
	    if ((hcalTPDepth[2] + hcalTPDepth[3] + hcalTPDepth[4] + hcalTPDepth[5] + hcalTPDepth[6]) < 10) continue;

	    if (abs(towEtaemu) <= maxTowerForHOvE && abs(towEtaemu) >= minTowerForHOvE){
		if (towHademu > maxTowerHad){
		    maxTowerIPhi = towPhiemu;
		    maxTowerIEta = towEtaemu;
		    maxTowerHad = towHademu;
		    foundMax = true;
		}
	    }
	}

	for(uint HcalTPIt = 0; HcalTPIt < nHCALTPemu; ++HcalTPIt){
	    towEtaemu = l1CaloTPemu_->hcalTPieta[HcalTPIt]; // use for HB HE restrictions                                 
	    towPhiemu = l1CaloTPemu_->hcalTPCaliphi[HcalTPIt];  // hcalTPiphi[HcalTPIt]; // use for deltaR calculation
	    towHademu = l1CaloTPemu_->hcalTPet[HcalTPIt]; // used for energy normalization in the energy depth plots    
	    int nDepths = l1CaloTPemu_->hcalTPnDepths[HcalTPIt];
	    std::vector<double> hcalTPDepth(7,0); 
	    if (nDepths == 0) continue;
	    hcalTPDepth[0] = l1CaloTPemu_->hcalTPDepth1[HcalTPIt];
	    hcalTPDepth[1] = l1CaloTPemu_->hcalTPDepth2[HcalTPIt];
	    hcalTPDepth[2] = l1CaloTPemu_->hcalTPDepth3[HcalTPIt];
	    hcalTPDepth[3] = l1CaloTPemu_->hcalTPDepth4[HcalTPIt];
	    hcalTPDepth[4] = l1CaloTPemu_->hcalTPDepth5[HcalTPIt];
	    hcalTPDepth[5] = l1CaloTPemu_->hcalTPDepth6[HcalTPIt];
	    hcalTPDepth[6] = l1CaloTPemu_->hcalTPDepth7[HcalTPIt];

	    if (abs(towEtaemu) <= maxTowerForHOvE && abs(towEtaemu) >= minTowerForHOvE){
		for (int iTopTowerIEta = -4; iTopTowerIEta <= 4; ++iTopTowerIEta){
		    for (int iTopTowerIPhi = -4; iTopTowerIPhi <= 4; ++iTopTowerIPhi){
			int wrappedIPhi = maxTowerIPhi + iTopTowerIPhi;
			if (wrappedIPhi > 72) wrappedIPhi -= 72;
			if (wrappedIPhi < 0) wrappedIPhi += 72;
			if (towEtaemu == maxTowerIEta+iTopTowerIEta && towPhiemu == wrappedIPhi){
			    maxTower9x9Had += towHademu;
			    for (int i = 0; i < nDepths-1; i++) depthVarTower9x9Had[i] = hcalTPDepth[i];
			    if (abs(iTopTowerIPhi) <= 1 && abs(iTopTowerIEta) <= 1){
				maxTower3x3Had += towHademu;
				for (int i = 0; i < nDepths-1; i++) depthVarTower3x3Had[i] += hcalTPDepth[i];
				if (abs(iTopTowerIPhi) == 0 && abs(iTopTowerIEta) == 0){
				    for (int i = 0; i < nDepths-1; i++) depthVarTowerHad[i] += hcalTPDepth[i];
				}
			    }
			}
		    }
		}
	    }
	}
	if (maxTowerHad > 0){
	    for(uint EcalTPIt = 0; EcalTPIt < nECALTPemu; ++EcalTPIt){
		towEtaemu = l1CaloTPemu_->ecalTPieta[EcalTPIt]; // use for HB HE restrictions                                 
		towPhiemu = l1CaloTPemu_->ecalTPCaliphi[EcalTPIt];  // ecalTPiphi[EcalTPIt]; // use for deltaR calculation
		towEmemu = l1CaloTPemu_->ecalTPet[EcalTPIt]; // used for energy normalization in the energy depth plots    

		if (abs(towEtaemu) <= maxTowerForHOvE && abs(towEtaemu) >= minTowerForHOvE){
		    for (int iTopTowerIEta = -4; iTopTowerIEta <= 4; ++iTopTowerIEta){
			for (int iTopTowerIPhi = -4; iTopTowerIPhi <= 4; ++iTopTowerIPhi){
			    int wrappedIPhi = maxTowerIPhi + iTopTowerIPhi;
			    if (wrappedIPhi > 72) wrappedIPhi -= 72;
			    if (wrappedIPhi < 0) wrappedIPhi += 72;
			    if (towEtaemu == maxTowerIEta+iTopTowerIEta && towPhiemu == wrappedIPhi){
				maxTower9x9Em += towEmemu;
				if (abs(iTopTowerIPhi) == 0 && abs(iTopTowerIEta) == 0){
				    maxTowerEm = towEmemu;
				}
				if (abs(iTopTowerIPhi) <= 1 && abs(iTopTowerIEta) <= 1){
				    maxTower3x3Em += towEmemu;
				}
			    }
			}
		    }
		}
	    }
	}

	depthVariablesAllTowers["HOvE"] = depthVarTowerHad;
	depthVariablesAllTowers["HOvE3"] = depthVarTowerHad;
	depthVariablesAllTowers["HOvE9"] = depthVarTowerHad;
	depthVariablesAllTowers["H3OvE3"] = depthVarTower3x3Had;
	depthVariablesAllTowers["H3OvE9"] = depthVarTower3x3Had;
	depthVariablesAllTowers["H9OvE9"] = depthVarTower9x9Had;

	hadVariablesMaxTower["HOvE"] = maxTowerHad;
	hadVariablesMaxTower["HOvE3"] = maxTowerHad;
	hadVariablesMaxTower["HOvE9"] = maxTowerHad;
	hadVariablesMaxTower["H3OvE3"] = maxTower3x3Had;
	hadVariablesMaxTower["H3OvE9"] = maxTower3x3Had;
	hadVariablesMaxTower["H9OvE9"] = maxTower9x9Had;

	emVariablesMaxTower["HOvE"] = maxTowerEm;
	emVariablesMaxTower["HOvE3"] = maxTower3x3Em;
	emVariablesMaxTower["HOvE9"] = maxTower9x9Em;
	emVariablesMaxTower["H3OvE3"] = maxTower3x3Em;
	emVariablesMaxTower["H3OvE9"] = maxTower9x9Em;
	emVariablesMaxTower["H9OvE9"] = maxTower9x9Em;

	hMaxTowerEt->Fill(maxTowerEt);
	hMaxTowerHad->Fill(maxTowerHad);
	hMaxTowerEm->Fill(maxTowerEm);

	for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
	    if (hadVariablesMaxTower[*ratioIt] == 0) ratioTower = -10.1;
	    else if (emVariablesMaxTower[*ratioIt] == 0) ratioTower = 10.1;
	    if (hadVariablesMaxTower[*ratioIt] > 0) ratioTower = TMath::Log10(hadVariablesMaxTower[*ratioIt]/emVariablesMaxTower[*ratioIt]);

	    maxTower2DHists["maxTower"+*ratioIt]->Fill(hadVariablesMaxTower[*ratioIt],ratioTower);
	    for (auto hOvEIt = hOvEThresholds.begin(); hOvEIt != hOvEThresholds.end(); hOvEIt++){
		for (auto hadIt = hadThresholds.begin(); hadIt != hadThresholds.end(); hadIt++){
		    if (foundMax){
			if (ratioTower >= hOvEIt->second && hadVariablesMaxTower[*ratioIt] >= hadIt->second){
			    maxTower1DHists["maxTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(ratioTower);
			    maxTower1DHists["maxTowerHad"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(hadVariablesMaxTower[*ratioIt]);
			    if (hadIt->first == "30" && hOvEIt->first == "Inf") maxTowerEtaPhi->Fill(maxTowerIEta,maxTowerIPhi);
			    maxTower1DHists["maxTowerDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(hadVariablesMaxTower[*ratioIt]);
			    for (int i = 0; i < 7; i++) {
				maxTower1DHists["maxTowerDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(i,depthVariablesAllTowers[*ratioIt][i]);
				maxTower2DHists["maxTowerDepth2D"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(i,depthVariablesAllTowers[*ratioIt][i]);
			    }
			}
		    }
		    double maxPt = -1;
		    double maxRatio = -100;
		    int maxEta = -1;
		    int maxPhi = -1;
		    int foundJet = -1;
		    double maxDenom = -1;
		    double maxDenomHigh = -1;
		    double maxHad = -1;
		    double ratioJet = -100;
		    double maxDepthRatio = 0;
		    double minEcal = 1E6;
		    double minLowDepth = 1E6;
		    int nJetPassSel = 0;
		    for (uint jetIt = 0; jetIt < nJetemu; ++jetIt){
			seedTowerIPhi = l1emu_->jetTowerIPhi[jetIt];
			seedTowerIEta = l1emu_->jetTowerIEta[jetIt];
			if (hadVariablesAllJets[*ratioIt][jetIt] == 0) ratioJet = -10.1;
			else if (emVariablesAllJets[*ratioIt][jetIt] == 0) ratioJet = 10.1;
			else ratioJet = TMath::Log10(hadVariablesAllJets[*ratioIt][jetIt]/emVariablesAllJets[*ratioIt][jetIt]);
			double denom = depthVariablesAllJets[*ratioIt][jetIt][2]+depthVariablesAllJets[*ratioIt][jetIt][3]+depthVariablesAllJets[*ratioIt][jetIt][4]+depthVariablesAllJets[*ratioIt][jetIt][5]+depthVariablesAllJets[*ratioIt][jetIt][6];
			double denomHigh = depthVariablesAllJets[*ratioIt][jetIt][3]+depthVariablesAllJets[*ratioIt][jetIt][4]+depthVariablesAllJets[*ratioIt][jetIt][5]+depthVariablesAllJets[*ratioIt][jetIt][6];
			double total = denom + depthVariablesAllJets[*ratioIt][jetIt][0] + depthVariablesAllJets[*ratioIt][jetIt][1];
			double depthRatio;

			if (total == 0) depthRatio = 0;
			else depthRatio = denom/total;

			double ecal = emVariablesAllJets[*ratioIt][jetIt];

			if (abs(seedTowerIEta) <= maxTowerForHOvE && abs(seedTowerIEta) >= minTowerForHOvE && ratioJet >= hOvEIt->second && denom >= hadIt->second){
			// if (abs(seedTowerIEta) <= maxTowerForHOvE && abs(seedTowerIEta) >= minTowerForHOvE && ratioJet >= hOvEIt->second && hadVariablesAllJets[*ratioIt][jetIt] >= hadIt->second){
			    nJetPassSel++;
			    if (maxPt < l1emu_->jetEt[jetIt]) {
				maxPt = l1emu_->jetEt[jetIt];
				maxPhi = seedTowerIPhi;
				maxEta = seedTowerIEta;
			    }
			    if (maxRatio < ratioJet){
				maxRatio = ratioJet;
			    }
			    if (maxDenomHigh < denomHigh){
				maxDenomHigh = denomHigh;
			    }
			    if (maxDenom < denom){
				foundJet = jetIt;
				maxDenom = denom;
			    }
			    if (maxHad < total){
				maxHad = total;
			    }
			    if (maxDepthRatio < depthRatio){
				maxDepthRatio = depthRatio;
			    }
			    if (minLowDepth > (total-denom)){
				minLowDepth = total-denom;
			    }
			    if (minEcal > ecal){
				minEcal = ecal;
			    }
			}
		    }
			maxJet1DHists["nJetPassSel"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(nJetPassSel);
		    if (foundJet >= 0){
			// if (maxDenom > 30 && hadIt->first == "None" && hOvEIt->first == "None" && maxPt > 0 && *ratioIt=="H3OvE3"){
			//     std::cout << "Decay in HCAL " << decayInHCAL << std::endl;
			//     for (auto partIt :genMatchId){
			// 	std::cout << "Radius: " << TMath::Sqrt(l1gen_->partVx[partIt]*l1gen_->partVx[partIt]+l1gen_->partVy[partIt]*l1gen_->partVy[partIt]) << "Z: " << l1gen_->partVz[partIt] << std::endl;
			//     }
			// 	std::cout << "HCAL Barrel radius/Z: 180-270, Z < 390 HCAL Endcap radius/Z: 44-270, Z > 390" << std::endl;
			//     std::cout << std::endl;
			// }
			maxJet1DHists["maxJetDenom"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(maxDenom);
			maxJet1DHists["maxJetDenomIt"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(foundJet);
			maxJet1DHists["maxJetDenomHigh"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(maxDenomHigh);
			maxJet1DHists["maxJetTotalHad"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(maxHad);
			maxJet1DHists["minJetEcal"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(minEcal);
			maxJet1DHists["minJetLowDepth"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(minLowDepth);
			for (int i = 0; i < 7; i++) {
			    maxJet1DHists["maxJetDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(i,depthVariablesAllJets[*ratioIt][foundJet][i]);
			    maxJet2DHists["maxJetDepth2D"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(i,depthVariablesAllJets[*ratioIt][foundJet][i]);
			}
			maxJet1DHists["maxJetDepthRatio"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(maxDepthRatio);
			if (abs(maxEta) <= maxTowerBarrel) {
			    maxJet2DHists["maxJetDepth1Vs3"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(depthVariablesAllJets[*ratioIt][foundJet][1],depthVariablesAllJets[*ratioIt][foundJet][3]);
			    maxJet2DHists["maxJetDepth2Vs3"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(depthVariablesAllJets[*ratioIt][foundJet][2],depthVariablesAllJets[*ratioIt][foundJet][3]);
			    maxJet2DHists["maxJetDepth1Vs2"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(depthVariablesAllJets[*ratioIt][foundJet][1],depthVariablesAllJets[*ratioIt][foundJet][2]);
			}
			hadVariablesSumJets[hOvEIt->first+"Had"+hadIt->first+*ratioIt] += hadVariablesAllJets[*ratioIt][foundJet];
			if (maxPt > 0) maxJet1DHists["maxJetPt"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(maxPt);
			if (maxRatio > -100) maxJet1DHists["maxJetHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(maxRatio);
			if (hadIt->first == "30" && hOvEIt->first == "Inf" && maxPt > 0) maxJetEtaPhi->Fill(maxEta,maxPhi);
		    }

		}
	    }
	}
	}

	// End event loop, now plot histos
	// TCanvas* canvas = new TCanvas("canvas","",750,700);
	//
	int ecalColour = 2;
	int hcalColour = 4;
	int towerColour = 3;
	int jetColour = 3;

	TFile outputFile = TFile(outFileName,"RECREATE");
	outputFile.cd();
	nTotHist->Write();
	hMaxTowerHad->Write();
	maxTowerEtaPhi->Write();
	maxJetEtaPhi->Write();
	TDirectory * maxTower2DDir = outputFile.mkdir("maxTower2DHists");
	maxTower2DDir->cd();
	for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
	    formatPlot2D(maxTower2DHists["maxTower"+*ratioIt]);
	    maxTower2DHists["maxTower"+*ratioIt]->Write();
	}
	outputFile.cd();
	TDirectory * ptDirJet = outputFile.mkdir("ptHistsJet");
	for (auto hOvEIt = hOvEThresholds.begin(); hOvEIt != hOvEThresholds.end(); hOvEIt++){
	    TDirectory * hOvEDir = ptDirJet->mkdir("hOvE"+hOvEIt->first);
	    hOvEDir->cd();
	    for (auto hadIt = hadThresholds.begin(); hadIt != hadThresholds.end(); hadIt++){
		TDirectory * hadDir = hOvEDir->mkdir("had"+hadIt->first);
		hadDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxJet1DHists["maxJetHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    formatPlot1D(maxJet1DHists["maxJetPt"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxJet1DHists["maxJetHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		    maxJet1DHists["maxJetPt"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
		TDirectory * depthDir = hadDir->mkdir("depth");
		depthDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxJet1DHists["maxJetDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxJet1DHists["maxJetDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Scale(1./maxJet1DHists["maxJetDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->GetEntries());
		    maxJet1DHists["maxJetDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		    maxJet2DHists["maxJetDepth2D"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		    maxJet2DHists["maxJetDepth1Vs3"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		    maxJet2DHists["maxJetDepth1Vs2"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		    maxJet2DHists["maxJetDepth2Vs3"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
		TDirectory * depthRatioDir = hadDir->mkdir("depthRatio");
		depthRatioDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxJet1DHists["maxJetDepthRatio"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxJet1DHists["maxJetDepthRatio"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
		TDirectory * depthDenomDir = hadDir->mkdir("depthDenom");
		depthDenomDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxJet1DHists["maxJetDenom"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxJet1DHists["maxJetDenom"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
		TDirectory * depthDenomItDir = hadDir->mkdir("depthDenomIt");
		depthDenomItDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxJet1DHists["maxJetDenomIt"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxJet1DHists["maxJetDenomIt"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
		TDirectory * depthDenomHighDir = hadDir->mkdir("depthDenomHigh");
		depthDenomHighDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxJet1DHists["maxJetDenomHigh"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxJet1DHists["maxJetDenomHigh"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
		TDirectory * totalHadDir = hadDir->mkdir("totalHad");
		totalHadDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxJet1DHists["maxJetTotalHad"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxJet1DHists["maxJetTotalHad"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
		TDirectory * minEcalDir = hadDir->mkdir("minEcal");
		minEcalDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxJet1DHists["minJetEcal"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxJet1DHists["minJetEcal"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
		TDirectory * minLowDepthDir = hadDir->mkdir("minLowDepth");
		minLowDepthDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxJet1DHists["minJetLowDepth"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxJet1DHists["minJetLowDepth"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
		TDirectory * nJetPassSelDir = hadDir->mkdir("nJetPassSel");
		nJetPassSelDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxJet1DHists["nJetPassSel"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxJet1DHists["nJetPassSel"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
	    }
	}
	TDirectory * ptDirTower = outputFile.mkdir("ptHistsTower");
	for (auto hOvEIt = hOvEThresholds.begin(); hOvEIt != hOvEThresholds.end(); hOvEIt++){
	    TDirectory * hOvEDir = ptDirTower->mkdir("hOvE"+hOvEIt->first);
	    hOvEDir->cd();
	    for (auto hadIt = hadThresholds.begin(); hadIt != hadThresholds.end(); hadIt++){
		TDirectory * hadDir = hOvEDir->mkdir("had"+hadIt->first);
		hadDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxTower1DHists["maxTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    formatPlot1D(maxTower1DHists["maxTowerHad"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxTower1DHists["maxTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		    maxTower1DHists["maxTowerHad"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
		TDirectory * depthDir = hadDir->mkdir("depth");
		depthDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(maxTower1DHists["maxTowerDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    maxTower1DHists["maxTowerDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Scale(1./maxTower1DHists["maxTowerDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->GetEntries());
		    maxTower1DHists["maxTowerDepthProfile"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		    maxTower2DHists["maxTowerDepth2D"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
	    }
	}
	return 0;
    }
    bool checkRootFileOk(TFile * tf){
	gErrorIgnoreLevel=kSysError;
	// TFile * tf = new TFile(fileName,"READ");
	bool isOK = false;
	if (!tf) isOK = false;
	else if (!(tf->IsZombie() or tf->TestBit(TFile::kRecovered))) isOK = true;
	// tf->Close();
	// delete tf;
	return isOK;
    }
    //main plotting function
    int main(int argc, char *argv[]){
	TString inFileName;
	TString outFileName;
	TString puSelStr;
	TString singleJetStr = "0";
	TString genDecayInHCALStr = "0";
	TString genMatchStr = "0";
	TString subSystem = "barrel";
	if (argc >= 3){
	    inFileName = argv[1];
	    outFileName = argv[2];
	    if (argc == 4){
		puSelStr = argv[3];
	    }
	    else if (argc == 5){
		puSelStr = argv[3];
		singleJetStr = argv[4];
	    }
	    else if (argc == 6){
		puSelStr = argv[3];
		singleJetStr = argv[4];
		genDecayInHCALStr = argv[5];
	    }
	    else if (argc == 7){
		puSelStr = argv[3];
		singleJetStr = argv[4];
		genDecayInHCALStr = argv[5];
		genMatchStr = argv[6];
	    }
	    else if (argc == 8){
		puSelStr = argv[3];
		singleJetStr = argv[4];
		genDecayInHCALStr = argv[5];
		genMatchStr = argv[6];
		subSystem = argv[7];
	    }
	    else {
		std::cout << "Usage: <inFile> <outFile> <min pu> <rate above singleJet> <gen decay in HCAL> <subsystem>" << std::endl;
		return 0;
	    }
	}
	else{
	    std::cout << "Usage: <inFile> <outFile> <min pu> <rate above singleJet> <gen decay in HCAL> <subsystem>" << std::endl;
	    return 0;
	}
	TFile * inputFile = TFile::Open(argv[1],"READ");
	if (!checkRootFileOk(inputFile)) 
	{
	    std::cout << "File " << argv[1] << " is corrupt or doesn't exist!" << std::endl;
	    return 0;
	}
	inputFile->Close();
	doHOvEStudy(argv[1],argv[2],puSelStr,singleJetStr,genDecayInHCALStr,genMatchStr,subSystem);
	return 0;
    }
