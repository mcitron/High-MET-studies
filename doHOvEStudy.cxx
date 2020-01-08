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
    plot1d->SetMinimum(0.);
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

int doHOvEStudy(TString inFileName,TString outFileName, TString puSelStr, TString singleJetStr,TString genStr){

    std::cout << "Running doHOvEStudy.cxx..." << std::endl;
    uint puSel = puSelStr.Atoi();
    bool rateAboveSingleJet = (singleJetStr.Atoi() > 0);
    bool genDecayInHCAL = (genStr.Atoi() > 0);
    std::vector<TString> towers;
    // std::vector<TString> jets;
    // std::vector<TString> gen;

    bool isTest = false;
    bool isData = !genDecayInHCAL;
    int maxTowerBarrel = 16;
    int maxTowerEndcap = 28;

    int minTowerForHOvE = maxTowerBarrel+1;
    int maxTowerForHOvE = maxTowerEndcap;

    // Load files
    TChain * treeL1Gen;

    std::cout << "Loading up the TChain..." << std::endl;
    TChain * eventTree = new TChain("l1EventTree/L1EventTree");
    TChain * treeL1Towemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
    TChain * treeL1emu = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
    TChain * treeL1reco = new TChain("l1RecoTree/RecoTree");

    towers.push_back(inFileName);
    // if (!isData){
    //     towers.push_back("/home/users/mcitron/triggerStudies/CMSSW_10_3_1/src/L1Ntuple.root");
    //     jets.push_back("/home/users/mcitron/triggerStudies/CMSSW_10_3_1/src/L1Ntuple.root");
    //     treeL1Gen = new TChain("l1GeneratorTree/L1GenTree");
    // }
    // else{
    //     towers.push_back("/hadoop/cms/store/user/mcitron/ZeroBias/zbMETa/181212_124635/0000/L1Ntuple_10.root");
    //     jets.push_back("/hadoop/cms/store/user/mcitron/ZeroBias/zbMETa/181212_124635/0000/L1Ntuple_10.root");
    //
    // }

    uint minFiles = towers.size();

    for(uint i = 0; i < minFiles; ++i) {
	eventTree->Add(towers[i]);
	treeL1Towemu->Add(towers[i]);
	treeL1emu->Add(towers[i]);
	treeL1reco->Add(towers[i]);
	if (!isData){
	    treeL1Gen->Add(towers[i]);
	}
    }

    L1Analysis::L1AnalysisEventDataFormat           *event_ = new L1Analysis::L1AnalysisEventDataFormat();
    L1Analysis::L1AnalysisL1CaloTowerDataFormat     *l1Towemu_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
    L1Analysis::L1AnalysisL1UpgradeDataFormat       *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
    L1Analysis::L1AnalysisGeneratorDataFormat       *l1gen_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
    L1Analysis::L1AnalysisRecoVertexDataFormat       *l1reco_ = new L1Analysis::L1AnalysisRecoVertexDataFormat();
    // L1Analysis::L1AnalysisL1UpgradeDataFormat       *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();

    eventTree->SetBranchAddress("Event", &event_);
    treeL1Towemu->SetBranchAddress("L1CaloTower", &l1Towemu_);
    treeL1emu->SetBranchAddress("L1Upgrade", &l1emu_);
    treeL1reco->SetBranchAddress("Vertex", &l1reco_);
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
    TH1D* nTotHist = new TH1D("nTot",";;nPass",2,0,2);

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
    TH2D * hMaxTowerVsSecondTowerHad = new TH2D("maxTowerHadVsSecondTowerHad",";ET;",100,0,200,100,0,200);

    std::vector<TString> ratioStrings = {"HOvE","HOvE3","HOvE9","H3OvE3","H9OvE9"};
    std::map<const TString, TH2D*> maxTower2DHists;
    // std::map<const TString, TH2D*> maxJet2DHists;
    std::map<const TString, TH1D*> maxTower1DHists;
    std::map<const TString, TH1D*> secondTower1DHists;
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

    for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
	maxTower2DHists["maxTower"+*ratioIt] = new TH2D("maxTower"+*ratioIt,";ET;log(H/E)",100,0,500,110,-11,11);
	for (auto hOvEIt = hOvEThresholds.begin(); hOvEIt != hOvEThresholds.end(); hOvEIt++){
	    for (auto hadIt = hadThresholds.begin(); hadIt != hadThresholds.end(); hadIt++){
		maxJet1DHists["maxJetHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxJetHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";ET;",1000,0,1000);
		maxTower1DHists["maxTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("maxTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";ET;",500,0,500);
		secondTower1DHists["secondTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt] = new TH1D("secondTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt,";ET;",500,0,500);
	    }
	}
    }
    std::map<const TString, double> hadVariablesMaxTower;
    std::map<const TString, double> emVariablesMaxTower;
    std::map<const TString, double> hadVariablesSecondTower;
    std::map<const TString, double> emVariablesSecondTower;


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
	std::map<const TString, std::vector<double> > emVariablesAllJets;
	// initialise some variables
	uint nHCALTPemu(0), nECALTPemu(0), nTowemu(0), nPart(0),nJetemu(0);
	int maxTowerIPhi(-1),maxTowerIEta(-1),secondTowerIPhi(0),secondTowerIEta(0);
	double hcalTPEtEm(0), ecalTPEtEm(0), towEtemu(0), towHademu(0),  towEmemu(0),
	       maxTowerEt(0), maxTowerHad(0), maxTowerEm(0), maxTower3x3Em(0), 
	       maxTower3x3Had(0), maxTower9x9Em(0), maxTower9x9Had(0),ratioTower(0),ratioJet(0),
	       secondTowerEt(0), secondTowerHad(0), secondTowerEm(0), secondTower3x3Em(0), 
	       secondTower3x3Had(0), secondTower9x9Em(0), secondTower9x9Had(0),ratioSecondTower(0);

	for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
	    hadVariablesMaxTower[*ratioIt] = 0;
	    emVariablesMaxTower[*ratioIt] = 0;
	    hadVariablesSecondTower[*ratioIt] = 0;
	    emVariablesSecondTower[*ratioIt] = 0;
	    // hadVariablesAllJets[*ratioIt] = std::vector<double>;
	    // emVariablesMaxTower[*ratioIt] = std::vector<double>;
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
	treeL1reco->GetEntry(jentry);
	if (!isData){
	    treeL1Gen->GetEntry(jentry);
	}
	treeL1emu->GetEntry(jentry);
	uint nVtxReco = l1reco_->nVtx;
	if (nVtxReco < puSel) continue;
	nTowemu = l1Towemu_->nTower;
	nJetemu = l1emu_->nJets;
	nTotHist->Fill(1.5);
	// // gen stuff
	if (genDecayInHCAL){
	    bool decayInHCAL = false;
	    nPart = l1gen_->nPart;
	    for(uint partIt = 0; partIt < nPart; ++partIt){
		if (abs(l1gen_->partId[partIt]) == 1000039) {
		    double hcalRadius = TMath::Sqrt(l1gen_->partVx[partIt]*l1gen_->partVx[partIt]+l1gen_->partVy[partIt]*l1gen_->partVy[partIt]);
		    if (hcalRadius > 180 && hcalRadius < 270 && abs(l1gen_->partVz[partIt]) < 390) decayInHCAL = true;
		    if (hcalRadius > 44 && hcalRadius < 270 && abs(l1gen_->partVz[partIt]) > 390 && abs(l1gen_->partVz[partIt]) < 568) decayInHCAL = true;
		}
	    }
	}
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
		for(uint towIt = 0; towIt < nTowemu; ++towIt){
		    towEtemu  = l1Towemu_->iet[towIt];
		    towHademu = l1Towemu_->ihad[towIt];
		    towEmemu  = l1Towemu_->iem[towIt];
		    towEtaemu = l1Towemu_->ieta[towIt];
		    towPhiemu = l1Towemu_->iphi[towIt];
		    if (abs(towEtaemu) >= minTowerForHOvE && abs(towEtaemu) <= maxTowerForHOvE){
			if (towEtaemu == seedTowerIEta && towPhiemu == seedTowerIPhi){
			    seedTowerHad = towHademu;
			    seedTowerEm = towEmemu;
			}
			for (int iSeedTowerIEta = -4; iSeedTowerIEta <= 4; ++iSeedTowerIEta){
			    for (int iSeedTowerIPhi = -4; iSeedTowerIPhi <= 4; ++iSeedTowerIPhi){
				int wrappedIPhi = seedTowerIPhi+iSeedTowerIPhi;
				if (wrappedIPhi > 72) wrappedIPhi -= 72;
				if (wrappedIPhi < 0) wrappedIPhi += 72;
				if (towEtaemu == seedTowerIEta+iSeedTowerIEta && towPhiemu == wrappedIPhi){
				    seedTower9x9Em += towEmemu;
				    seedTower9x9Had += towHademu;
				    if (abs(iSeedTowerIPhi) <= 1 && abs(iSeedTowerIEta) <= 1){
					seedTower3x3Em += towEmemu;
					seedTower3x3Had += towHademu;
				    }
				}
			    }
			}
		    }
		}

		hadVariablesAllJets["HOvE"].push_back(seedTowerHad);
		hadVariablesAllJets["HOvE3"].push_back(seedTowerHad);
		hadVariablesAllJets["HOvE9"].push_back(seedTowerHad);
		hadVariablesAllJets["H3OvE3"].push_back(seedTower3x3Had);
		hadVariablesAllJets["H9OvE9"].push_back(seedTower9x9Had);

		emVariablesAllJets["HOvE"].push_back(seedTowerEm);
		emVariablesAllJets["HOvE3"].push_back(seedTower3x3Em);
		emVariablesAllJets["HOvE9"].push_back(seedTower9x9Em);
		emVariablesAllJets["H3OvE3"].push_back(seedTower3x3Em);
		emVariablesAllJets["H9OvE9"].push_back(seedTower9x9Em);
	    }

	    // if (decayInHCAL){
	    //     hMaxJetEt->Fill(l1emu_->jetEt[0]);
	    // }

	    // Max tower stuff
	    for(uint towIt = 0; towIt < nTowemu; ++towIt){
		towEtemu  = l1Towemu_->iet[towIt];
		towHademu  = l1Towemu_->ihad[towIt];
		towEmemu  = l1Towemu_->iem[towIt];
		towEtaemu = l1Towemu_->ieta[towIt];
		towPhiemu = l1Towemu_->iphi[towIt];
		if (towEtemu > secondTowerEt && abs(towEtaemu) <= maxTowerForHOvE && abs(towEtaemu) >= minTowerForHOvE){
		    if (towEtemu > maxTowerEt){
			secondTowerIPhi = maxTowerIPhi;
			secondTowerIEta = maxTowerIEta; 
			secondTowerEt = maxTowerEt;
			secondTowerHad = maxTowerHad;
			secondTowerEm = maxTowerEm;

			maxTowerIPhi = towPhiemu;
			maxTowerIEta = towEtaemu;
			maxTowerEt = towEtemu;
			maxTowerHad = towHademu;
			maxTowerEm = towEmemu;
		    }
		    else{
			secondTowerIPhi = towPhiemu;
			secondTowerIEta = towEtaemu;
			secondTowerEt = towEtemu;
			secondTowerHad = towHademu;
			secondTowerEm = towEmemu;
		    }
		}
	    }

	    for(uint towIt = 0; towIt < nTowemu; ++towIt){
		towEtemu  = l1Towemu_->iet[towIt];
		towHademu  = l1Towemu_->ihad[towIt];
		towEmemu  = l1Towemu_->iem[towIt];
		towEtaemu = l1Towemu_->ieta[towIt];
		towPhiemu = l1Towemu_->iphi[towIt];

		if (abs(towEtaemu) <= maxTowerBarrel) {
		    hTowPhiB->Fill(towPhiemu);
		    hTowTPETphiB->Fill(towPhiemu, towEtemu);
		}

		else if (abs(towEtaemu) > maxTowerBarrel && abs(towEtaemu) <= maxTowerEndcap) {
		    hTowPhiE->Fill(towPhiemu);
		    hTowTPETphiE->Fill(towPhiemu, towEtemu);
		} // Ignoring HF
		if (abs(towEtaemu) <= maxTowerForHOvE && abs(towEtaemu) >= minTowerForHOvE){
		    hAllTowEt->Fill(towEtemu);
		    hAllTowEta->Fill(towEtaemu);
		    hTowTPETEta->Fill(towEtaemu, towEtemu);
		    for (int iTopTowerIEta = -4; iTopTowerIEta <= 4; ++iTopTowerIEta){
			for (int iTopTowerIPhi = -4; iTopTowerIPhi <= 4; ++iTopTowerIPhi){
			    int wrappedIPhi = maxTowerIPhi + iTopTowerIPhi;
			    if (wrappedIPhi > 72) wrappedIPhi -= 72;
			    if (wrappedIPhi < 0) wrappedIPhi += 72;
			    if (towEtaemu == maxTowerIEta+iTopTowerIEta && towPhiemu == wrappedIPhi){
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
		if (abs(towEtaemu) <= maxTowerForHOvE && abs(towEtaemu) >= minTowerForHOvE){
		    hAllTowEt->Fill(towEtemu);
		    hAllTowEta->Fill(towEtaemu);
		    hTowTPETEta->Fill(towEtaemu, towEtemu);
		    for (int iTopTowerIEta = -4; iTopTowerIEta <= 4; ++iTopTowerIEta){
			for (int iTopTowerIPhi = -4; iTopTowerIPhi <= 4; ++iTopTowerIPhi){
			    int wrappedIPhi = maxTowerIPhi + iTopTowerIPhi;
			    if (wrappedIPhi > 72) wrappedIPhi -= 72;
			    if (wrappedIPhi < 0) wrappedIPhi += 72;
			    if (towEtaemu == secondTowerIEta+iTopTowerIEta && towPhiemu == wrappedIPhi){
				secondTower9x9Em += towEmemu;
				secondTower9x9Had += towHademu;
				if (abs(iTopTowerIPhi) <= 1 && abs(iTopTowerIEta) <= 1){
				    secondTower3x3Em += towEmemu;
				    secondTower3x3Had += towHademu;
				}
			    }
			}
		    }
		}
	    }
	    hadVariablesSecondTower["HOvE"] = secondTowerHad;
	    hadVariablesSecondTower["HOvE3"] = secondTowerHad;
	    hadVariablesSecondTower["HOvE9"] = secondTowerHad;
	    hadVariablesSecondTower["H3OvE3"] = secondTower3x3Had;
	    hadVariablesSecondTower["H9OvE9"] = secondTower9x9Had;

	    hadVariablesMaxTower["HOvE"] = maxTowerHad;
	    hadVariablesMaxTower["HOvE3"] = maxTowerHad;
	    hadVariablesMaxTower["HOvE9"] = maxTowerHad;
	    hadVariablesMaxTower["H3OvE3"] = maxTower3x3Had;
	    hadVariablesMaxTower["H9OvE9"] = maxTower9x9Had;

	    emVariablesMaxTower["HOvE"] = maxTowerEm;
	    emVariablesMaxTower["HOvE3"] = maxTower3x3Em;
	    emVariablesMaxTower["HOvE9"] = maxTower9x9Em;
	    emVariablesMaxTower["H3OvE3"] = maxTower3x3Em;
	    emVariablesMaxTower["H9OvE9"] = maxTower9x9Em;

	    hMaxTowerEt->Fill(maxTowerEt);
	    hMaxTowerHad->Fill(maxTowerHad);
	    hMaxTowerVsSecondTowerHad->Fill(maxTowerHad,secondTowerHad);
	    hMaxTowerEm->Fill(maxTowerEm);

	    for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		if (emVariablesMaxTower[*ratioIt] == 0) ratioTower = 10.1;
		else if (hadVariablesMaxTower[*ratioIt] > 0) ratioTower = TMath::Log10(hadVariablesMaxTower[*ratioIt]/emVariablesMaxTower[*ratioIt]);
		else ratioTower = -10.1;
		if (emVariablesSecondTower[*ratioIt] == 0) ratioSecondTower = 10.1;
		else if (hadVariablesMaxTower[*ratioIt] > 0) ratioSecondTower = TMath::Log10(hadVariablesMaxTower[*ratioIt]/emVariablesMaxTower[*ratioIt]);
		else ratioSecondTower = -10.1;

		maxTower2DHists["maxTower"+*ratioIt]->Fill(hadVariablesMaxTower[*ratioIt],ratioTower);
		for (auto hOvEIt = hOvEThresholds.begin(); hOvEIt != hOvEThresholds.end(); hOvEIt++){
		    for (auto hadIt = hadThresholds.begin(); hadIt != hadThresholds.end(); hadIt++){
			if (ratioTower >= hOvEIt->second && hadVariablesMaxTower[*ratioIt] >= hadIt->second){
			    maxTower1DHists["maxTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(hadVariablesMaxTower[*ratioIt]);
			    if (hadIt->first == "30" && hOvEIt->first == "Inf") maxJetEtaPhi->Fill(maxTowerIEta,maxTowerIPhi);
			    if (ratioSecondTower >= hOvEIt->second && hadVariablesSecondTower[*ratioIt] >= hadIt->second){
				secondTower1DHists["secondTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(hadVariablesMaxTower[*ratioIt]);
			    }
			    else if (hadVariablesSecondTower[*ratioIt] <= hadIt->second){
				secondTower1DHists["secondTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(400);
			    }
			}
			double maxPt = -1;
			int maxEta = -1;
			int maxPhi = -1;
			for (uint jetIt = 0; jetIt < nJetemu; ++jetIt){
			    seedTowerIPhi = l1emu_->jetTowerIPhi[jetIt];
			    seedTowerIEta = l1emu_->jetTowerIEta[jetIt];
			    if (emVariablesAllJets[*ratioIt][jetIt] == 0) ratioJet = 10.1;
			    else if (hadVariablesAllJets[*ratioIt][jetIt]> 0) ratioJet = TMath::Log10(hadVariablesAllJets[*ratioIt][jetIt]/emVariablesAllJets[*ratioIt][jetIt]);
			    else ratioJet = -10.1;
			    if (abs(seedTowerIEta) <= maxTowerForHOvE && abs(seedTowerIEta) >= minTowerForHOvE && ratioJet >= hOvEIt->second && hadVariablesAllJets[*ratioIt][jetIt] >= hadIt->second){
				if (maxPt < l1emu_->jetEt[jetIt]) {
				    maxPt = l1emu_->jetEt[jetIt];
				    maxPhi = seedTowerIPhi;
				    maxEta = seedTowerIEta;
				}
			    }
			}
			if (maxPt > 0) maxJet1DHists["maxJetHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Fill(maxPt);
			if (hadIt->first == "30" && hOvEIt->first == "Inf" && maxPt > 0) maxJetEtaPhi->Fill(maxEta,maxPhi);
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
	hMaxTowerVsSecondTowerHad->Write();
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
		    maxJet1DHists["maxJetHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
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
		    maxTower1DHists["maxTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
	    }
	}
	TDirectory * ptDirSecondTower = outputFile.mkdir("ptHistsSecondTower");
	for (auto hOvEIt = hOvEThresholds.begin(); hOvEIt != hOvEThresholds.end(); hOvEIt++){
	    TDirectory * hOvEDir = ptDirSecondTower->mkdir("hOvE"+hOvEIt->first);
	    hOvEDir->cd();
	    for (auto hadIt = hadThresholds.begin(); hadIt != hadThresholds.end(); hadIt++){
		TDirectory * hadDir = hOvEDir->mkdir("had"+hadIt->first);
		hadDir->cd();
		for (auto ratioIt = ratioStrings.begin(); ratioIt != ratioStrings.end(); ratioIt++){
		    formatPlot1D(secondTower1DHists["secondTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt],jetColour);
		    secondTower1DHists["secondTowerHOvE"+hOvEIt->first+"Had"+hadIt->first+*ratioIt]->Write();
		}
	    }
	}
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
	TString singleJetStr;
	TString genDecayInHCALStr;
	if (argc >= 3){
	    inFileName = argv[1];
	    outFileName = argv[2];
	    if (argc == 4){
		puSelStr = argv[3];
		singleJetStr = "0";
		genDecayInHCALStr = "0";
	    }
	    else if (argc == 5){
		puSelStr = argv[3];
		singleJetStr = argv[4];
		genDecayInHCALStr = "0";
	    }
	    else if (argc == 6){
		puSelStr = argv[3];
		singleJetStr = argv[4];
		genDecayInHCALStr = argv[5];
	    }
	    else if (argc > 6){
		std::cout << "Usage: <inFile> <outFile> <min pu>" << std::endl;
		return 0;
	    }
	    else{
		puSelStr = "0";
		singleJetStr = "0";
	    }
	}
	else{
	    std::cout << "Usage: <inFile> <outFile> <min pu>" << std::endl;
	    return 0;
	}
	TFile * inputFile = TFile::Open(argv[1],"READ");
	if (!checkRootFileOk(inputFile)) 
	{
	    std::cout << "File " << argv[1] << " is corrupt or doesn't exist!" << std::endl;
	    return 0;
	}
	inputFile->Close();
	doHOvEStudy(argv[1],argv[2],puSelStr,singleJetStr,genDecayInHCALStr);
	return 0;
    }
