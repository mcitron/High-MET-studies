#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TChain.h"
#include "TAttMarker.h"
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

// Following headers help decode L1T ntuples
#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisCaloTPDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1CaloTowerDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoJetDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"

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
    plot1d->SetMinimum(0.);
    plot1d->SetLineColor(colour);
    plot1d->SetLineWidth(2);
    plot1d->Scale(1. / (double) nPassing);
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

//main plotting function
void doHOvEStudy(){

    std::cout << "Running doHOvEStudy.cxx..." << std::endl;

    vector<string> towers;
    vector<string> jets;
    vector<string> gen;

    bool isTest = false;
    bool isData = true;

    // Load files
    TChain * treeL1Gen;

    cout << "Loading up the TChain..." << endl;
    TChain * eventTree = new TChain("l1EventTree/L1EventTree");
    TChain * treeL1Towemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
    TChain * treeL1emu = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
    if (!isData){
	    towers.push_back("/home/users/mcitron/triggerStudies/CMSSW_10_3_1/src/L1Ntuple.root");
	    jets.push_back("/home/users/mcitron/triggerStudies/CMSSW_10_3_1/src/L1Ntuple.root");
	    treeL1Gen = new TChain("l1GeneratorTree/L1GenTree");
    }
    else{
	    towers.push_back("/hadoop/cms/store/user/mcitron/ZeroBias/zbMETa/181212_124635/0000/L1Ntuple_10.root");
	    jets.push_back("/hadoop/cms/store/user/mcitron/ZeroBias/zbMETa/181212_124635/0000/L1Ntuple_10.root");

    }

    int minFiles = std::min( towers.size(), jets.size() );

    for(uint i = 0; i < minFiles; ++i) {
	eventTree->Add(towers[i].c_str());
	treeL1Towemu->Add(towers[i].c_str());
	treeL1emu->Add(towers[i].c_str());
	if (!isData){
	treeL1Gen->Add(towers[i].c_str());
	}
    }

    L1Analysis::L1AnalysisEventDataFormat           *event_ = new L1Analysis::L1AnalysisEventDataFormat();
    L1Analysis::L1AnalysisL1CaloTowerDataFormat     *l1Towemu_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
    L1Analysis::L1AnalysisL1UpgradeDataFormat       *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
    L1Analysis::L1AnalysisGeneratorDataFormat       *l1gen_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
    // L1Analysis::L1AnalysisL1UpgradeDataFormat       *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();

    eventTree->SetBranchAddress("Event", &event_);
    treeL1Towemu->SetBranchAddress("L1CaloTower", &l1Towemu_);
    treeL1emu->SetBranchAddress("L1Upgrade", &l1emu_);
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

    std::vector<TString> ratioStrings = {"HOvE","HOvE3","HOvE9","H3OvE3","H9OvE9"};
    std::map<const TString, TH2D*> maxTower2DHists;
    // std::map<const TString, TH2D*> maxJet2DHists;
    std::map<const TString, TH1D*> maxTower1DHists;
    // std::map<const TString, TH1D*> maxJet1DHists;
    std::map<const TString, TH1D*> allTower1DHists;
    std::map<const TString, TH1D*> allJet1DHists;
    std::map<const TString, double> hOvEThresholds;
    hOvEThresholds["HOvE5"] = 0.699;
    hOvEThresholds["HOvE10"] = 1;
    hOvEThresholds["HOvE20"] = 1.30;
    hOvEThresholds["HOvE100"] = 2;
    hOvEThresholds["HOvEInf"] = 4;
    for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
	maxTower2DHists["maxTower"+*ratioIt] = new TH2D("maxTower"+*ratioIt,";ET;log(H/E)",100,0,500,110,-11,11);
	// maxJet2DHists["maxJet"+*ratioIt] = new TH2D("maxJet"+*ratioIt,";ET;log(H/E)",110,-11,11,100,0,500);
	for (auto hOvEIt = hOvEThresholds.begin(); hOvEIt != hOvEThresholds.end(); hOvEIt++){
	    allTower1DHists["allTower"+hOvEIt->first+*ratioIt] = new TH1D("allTower"+hOvEIt->first+*ratioIt,";ET;",500,0,500);
	    allJet1DHists["allJet"+hOvEIt->first+*ratioIt] = new TH1D("allJet"+hOvEIt->first+*ratioIt,";ET;",1000,0,1000);
	    maxTower1DHists["maxTower"+hOvEIt->first+*ratioIt] = new TH1D("maxTower"+hOvEIt->first+*ratioIt,";ET;",500,0,500);
	    // maxJet1DHists["maxJet"+hOvEIt->first+*ratioIt] = new TH1D("maxJet"+hOvEIt->first+*ratioIt,";ET;",1000,0,1000);
	}
    }
    std::map<const TString, double> hadVariables;
    std::map<const TString, double> emVariables;
    // Gen hists 
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
	++nPassing;
	// initialise some variables
	int nHCALTPemu(0), nECALTPemu(0), nTowemu(-1), nPart(-1),nJetemu(-1),maxTowerIPhi(-1),maxTowerIEta(-1);
	double hcalTPEtEm(0), ecalTPEtEm(0), towEtemu(0), towHademu(0),  towEmemu(0),
	       maxTowerEt(0), maxTowerHad(0), maxTowerEm(0), maxTower3x3Em(0), 
	       maxTower3x3Had(0), maxTower9x9Em(0), maxTower9x9Had(0),ratio(0);

	for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
	    hadVariables[*ratioIt] = 0;
	    emVariables[*ratioIt] = 0;
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
	if( (jentry % 10000) == 0 ) cout << "Done " << jentry << " events of " << nEvents << endl;
	if( (jentry % 1000) == 0 ) cout << "." << flush;

	eventTree->GetEntry(jentry);

	run = event_->run;
	int lumi = event_->lumi;
	int event = event_->event;

	treeL1Towemu->GetEntry(jentry);
	if (!isData){
	    treeL1Gen->GetEntry(jentry);
	}
	treeL1emu->GetEntry(jentry);

	nTowemu = l1Towemu_->nTower;
	nJetemu = l1emu_->nJets;
	for (uint jetIt = 0; jetIt < nJetemu; ++jetIt){
	    hJetEt->Fill(l1emu_->jetEt[jetIt]);
	}
	if (nJetemu){
	    seedTowerIPhi = l1emu_->jetTowerIPhi[0];
	    seedTowerIEta = l1emu_->jetTowerIEta[0];
	}
	bool decayInHCAL = true;
	// nPart = l1gen_->nPart;
	// for(uint partIt = 0; partIt < nPart; ++partIt){
	//     if (abs(l1gen_->partId[partIt]) == 1000039) {
	// 	double hcalRadius = TMath::Sqrt(l1gen_->partVx[partIt]*l1gen_->partVx[partIt]+l1gen_->partVy[partIt]*l1gen_->partVy[partIt]);
	// 	hPartVr->Fill(hcalRadius);
	// 	hPartVz->Fill(abs(l1gen_->partVz[partIt]));
	// 	if (hcalRadius > 180 && hcalRadius < 270 && abs(l1gen_->partVz[partIt]) < 390) decayInHCAL = true;
	//     }
	// }
	// if (decayInHCAL){
	//     hMaxJetEt->Fill(l1emu_->jetEt[0]);
	// }
	//
	// Retrieve tower objects from the emulator tree
	for(uint towIt = 0; towIt < nTowemu; ++towIt){
	    towEtemu  = l1Towemu_->iet[towIt];
	    towHademu  = l1Towemu_->ihad[towIt];
	    towEmemu  = l1Towemu_->iem[towIt];
	    towEtaemu = l1Towemu_->ieta[towIt];
	    towPhiemu = l1Towemu_->iphi[towIt];
	    if (towEtemu > maxTowerEt && abs(towEtaemu) <= 16){
		maxTowerIPhi = towPhiemu;
		maxTowerIEta = towEtaemu;
		maxTowerEt = towEtemu;
		maxTowerHad = towHademu;
		maxTowerEm = towEmemu;
	    }
	}

	for(uint towIt = 0; towIt < nTowemu; ++towIt){
	    towEtemu  = l1Towemu_->iet[towIt];
	    towHademu  = l1Towemu_->ihad[towIt];
	    towEmemu  = l1Towemu_->iem[towIt];
	    towEtaemu = l1Towemu_->ieta[towIt];
	    towPhiemu = l1Towemu_->iphi[towIt];

	    if (abs(towEtaemu) <= 16) {
		hTowPhiB->Fill(towPhiemu);
		hTowTPETphiB->Fill(towPhiemu, towEtemu);
	    }

	    else if (abs(towEtaemu) > 16 && abs(towEtaemu) <= 28) {
		hTowPhiE->Fill(towPhiemu);
		hTowTPETphiE->Fill(towPhiemu, towEtemu);
	    } // Ignoring HF
	    if (abs(towEtaemu) <= 16){
		hAllTowEt->Fill(towEtemu);
		hAllTowEta->Fill(towEtaemu);
		hTowTPETEta->Fill(towEtaemu, towEtemu);
		for (int iTopTowerIEta = -4; iTopTowerIEta <= 4; ++iTopTowerIEta){
		    for (int iTopTowerIPhi = -4; iTopTowerIPhi <= 4; ++iTopTowerIPhi){
			if (towEtaemu == maxTowerIEta+iTopTowerIEta && towPhiemu == maxTowerIPhi+iTopTowerIPhi){
			    maxTower9x9Em += towEmemu;
			    maxTower9x9Had += towHademu;
			    if (abs(iTopTowerIPhi) <= 1 && abs(iTopTowerIEta) <= 1){
				maxTower3x3Em += towEmemu;
				maxTower3x3Had += towHademu;
			    }
			}
		    }
		}
	    }
	}

	hadVariables["HOvE"] = maxTowerHad;
	hadVariables["HOvE3"] = maxTowerHad;
	hadVariables["HOvE9"] = maxTowerHad;
	hadVariables["H3OvE3"] = maxTower3x3Had;
	hadVariables["H9OvE9"] = maxTower9x9Had;

	emVariables["HOvE"] = maxTowerEm;
	emVariables["HOvE3"] = maxTower3x3Em;
	emVariables["HOvE9"] = maxTower9x9Em;
	emVariables["H3OvE3"] = maxTower3x3Em;
	emVariables["H9OvE9"] = maxTower9x9Em;

	hMaxTowerEt->Fill(maxTowerEt);
	hMaxTowerHad->Fill(maxTowerHad);
	hMaxTowerEm->Fill(maxTowerEm);

	for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
	    if (emVariables[*ratioIt] == 0) ratio = 10.1;
	    else if (hadVariables[*ratioIt] > 0) ratio = TMath::Log10(hadVariables[*ratioIt]/emVariables[*ratioIt]);
	    // std::cout << std::endl;
	    // std::cout << hadVariables[*ratioIt] << " "<< ratio << std::endl;
	    maxTower2DHists["maxTower"+*ratioIt]->Fill(hadVariables[*ratioIt],ratio);
	    for (auto hOvEIt = hOvEThresholds.begin(); hOvEIt != hOvEThresholds.end(); hOvEIt++){
		if (ratio >= hOvEIt->second){
		    maxTower1DHists["maxTower"+hOvEIt->first+*ratioIt]->Fill(hadVariables[*ratioIt]);
		}
	    }
	}

    }

    // End event loop, now plot histos
    TCanvas* canvas = new TCanvas("canvas","",750,700);

    int ecalColour = 2;
    int hcalColour = 4;
    int towerColour = 3;
    int jetColour = 3;

    // Plot TP ET vs iEta for ECAL, HCAL and towers
    // Need to do divide before plotting occupancy otherwise normalisation is wrong
    hTowTPETEta->Divide(hAllTowEta);
    hTowTPETEta->Scale(nPassing);
    formatPlot1D(hTowTPETEta, towerColour);
    canvas->SaveAs("./Plots/Towers/TowTPETEta.pdf");

    formatPlot1D(hAllTowEta, towerColour);
    canvas->SaveAs("./Plots/Towers/TowEta.pdf");

    hTowTPETphiB->Divide(hTowPhiB);
    hTowTPETphiB->Scale(nPassing);
    formatPlot1D(hTowTPETphiB, towerColour);
    canvas->SaveAs("./Plots/Towers/TowTPETphiB.pdf");

    hTowTPETphiE->Divide(hTowPhiE);
    hTowTPETphiE->Scale(nPassing);
    canvas->SetLogy();
    canvas->SetLogz();
    formatPlot1D(hTowTPETphiE, towerColour);
    canvas->SaveAs("./Plots/Towers/TowTPETphiE.pdf");

    formatPlot1D(hTowPhiB, towerColour);
    canvas->SaveAs("./Plots/Towers/TowPhiBarrel.pdf");

    formatPlot1D(hTowPhiE, towerColour);
    canvas->SaveAs("./Plots/Towers/TowPhiEndcap.pdf");

    formatPlot1D(hPartVr, towerColour);
    canvas->SaveAs("./Plots/Gen/PartVr.pdf");

    formatPlot1D(hPartVz, towerColour);
    canvas->SaveAs("./Plots/Gen/PartVz.pdf");

    formatPlot1D(hMaxTowerEt, jetColour);
    canvas->SaveAs("./Plots/Jets/maxTowerEt.pdf");

    formatPlot1D(hMaxTowerHad, jetColour);
    canvas->SaveAs("./Plots/Jets/maxTowerHad.pdf");

    formatPlot1D(hMaxTowerEm, jetColour);
    canvas->SaveAs("./Plots/Jets/maxTowerEm.pdf");

    for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
	canvas->SetLogy(0);
	formatPlot2D(maxTower2DHists["maxTower"+*ratioIt]);
	canvas->SaveAs("./Plots/Towers/maxTower"+*ratioIt+".pdf");
	canvas->SetLogy();
	for (auto hOvEIt = hOvEThresholds.begin(); hOvEIt != hOvEThresholds.end(); hOvEIt++){
	    formatPlot1D(maxTower1DHists["maxTower"+hOvEIt->first+*ratioIt],towerColour);
	    canvas->SaveAs("./Plots/Towers/maxTower"+hOvEIt->first+*ratioIt+".pdf");
	}
    }
}
