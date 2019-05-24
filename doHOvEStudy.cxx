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

    // Load files
    towers.push_back("/home/users/mcitron/triggerStudies/CMSSW_10_3_1/src/L1Ntuple.root");
    jets.push_back("/home/users/mcitron/triggerStudies/CMSSW_10_3_1/src/L1Ntuple.root");

    cout << "Loading up the TChain..." << endl;
    TChain * eventTree = new TChain("l1EventTree/L1EventTree");
    TChain * treeL1Towemu = new TChain("l1CaloTowerEmuTree/L1CaloTowerTree");
    TChain * treeL1emu = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
    TChain * treeL1Gen = new TChain("l1GeneratorTree/L1GenTree");

    int minFiles = std::min( towers.size(), jets.size() );

    for(uint i = 0; i < minFiles; ++i) {
	eventTree->Add(towers[i].c_str());
	treeL1Towemu->Add(towers[i].c_str());
	treeL1emu->Add(towers[i].c_str());
	treeL1Gen->Add(towers[i].c_str());
    }

    L1Analysis::L1AnalysisEventDataFormat           *event_ = new L1Analysis::L1AnalysisEventDataFormat();
    L1Analysis::L1AnalysisL1CaloTowerDataFormat     *l1Towemu_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
    L1Analysis::L1AnalysisL1UpgradeDataFormat       *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
    L1Analysis::L1AnalysisGeneratorDataFormat       *l1gen_ = new L1Analysis::L1AnalysisGeneratorDataFormat();
    // L1Analysis::L1AnalysisL1UpgradeDataFormat       *l1emu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();

    eventTree->SetBranchAddress("Event", &event_);
    treeL1Towemu->SetBranchAddress("L1CaloTower", &l1Towemu_);
    treeL1emu->SetBranchAddress("L1Upgrade", &l1emu_);
    treeL1Gen->SetBranchAddress("Generator", &l1gen_);

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
    TH1D * hMaxJetEtDecayInHcal = new TH1D("maxJetEtDecayInHcal",";ET;",100,0,2000);
    TH1D * hMaxJetSeedEtDecayInHcal = new TH1D("maxJetSeedEtDecayInHcal",";ET;",100,0,1000);
    TH1D * hMaxJetSeedEmDecayInHcal = new TH1D("maxJetSeedEmDecayInHcal",";ET;",100,0,1000);
    TH1D * hMaxJetSeedHadDecayInHcal = new TH1D("maxJetSeedHadDecayInHcal",";ET;",100,0,1000);
    TH1D * hMaxJetSeedHOvEDecayInHcal = new TH1D("maxJetSeedHOvEDecayInHcal",";ET;",110,-11,11);
    TH1D * hMaxJetSeedHOvE9x9DecayInHcal = new TH1D("maxJetSeedHOvE9x9DecayInHcal",";ET;",110,-11,11);
    TH1D * hMaxJetSeedHOvE3x3DecayInHcal = new TH1D("maxJetSeedHOvE3x3DecayInHcal",";ET;",110,-11,11);

    TH1D * hMaxTopTowerEtDecayInHcal = new TH1D("maxTopTowerEtDecayInHcal",";ET;",100,0,1000);
    TH1D * hMaxTopTowerEmDecayInHcal = new TH1D("maxTopTowerEmDecayInHcal",";ET;",100,0,1000);
    TH1D * hMaxTopTowerHadDecayInHcal = new TH1D("maxTopTowerHadDecayInHcal",";ET;",100,0,1000);
    TH1D * hMaxTopTowerHOvEDecayInHcal = new TH1D("maxTopTowerHOvEDecayInHcal",";ET;",110,-11,11);
    TH1D * hMaxTopTowerHOvE9x9DecayInHcal = new TH1D("maxTopTowerHOvE9x9DecayInHcal",";ET;",110,-11,11);
    TH1D * hMaxTopTowerHOvE3x3DecayInHcal = new TH1D("maxTopTowerHOvE3x3DecayInHcal",";ET;",110,-11,11);
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
	int nHCALTPemu(0), nECALTPemu(0), nTowemu(-1), nPart(-1),nJetemu(-1),topTowerIPhi(-1),topTowerIEta(-1);
	double hcalTPEtEm(0), ecalTPEtEm(0), towEtemu(0), towHademu(0),
	       towEmemu(0), seedEt(0),seedEm(0),seedHad(0), 
	       seed3x3Em(0),seed3x3Had(0),seed9x9Em(0),seed9x9Had(0),
	       topTowerEt(0), topTowerHad(0), topTowerEm(0), topTower3x3Em(0), 
	       topTower3x3Had(0), topTower9x9Em(0), topTower9x9Had(0);
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
	treeL1Gen->GetEntry(jentry);
	treeL1emu->GetEntry(jentry);

	nTowemu = l1Towemu_->nTower;
	nPart = l1gen_->nPart;
	nJetemu = l1emu_->nJets;
	for (uint jetIt = 0; jetIt < nJetemu; ++jetIt){
	    hJetEt->Fill(l1emu_->jetEt[jetIt]);
	}
	if (nJetemu){
	    seedTowerIPhi = l1emu_->jetTowerIPhi[0];
	    seedTowerIEta = l1emu_->jetTowerIEta[0];
	}
	bool decayInHCAL = false;
	for(uint partIt = 0; partIt < nPart; ++partIt){
	    if (abs(l1gen_->partId[partIt]) == 1000039) {
		double hcalRadius = TMath::Sqrt(l1gen_->partVx[partIt]*l1gen_->partVx[partIt]+l1gen_->partVy[partIt]*l1gen_->partVy[partIt]);
		hPartVr->Fill(hcalRadius);
		hPartVz->Fill(abs(l1gen_->partVz[partIt]));
		if (hcalRadius > 180 && hcalRadius < 270 && abs(l1gen_->partVz[partIt]) < 390) decayInHCAL = true;
	    }
	}
	if (decayInHCAL){
	    hMaxJetEtDecayInHcal->Fill(l1emu_->jetEt[0]);
	}

	// Retrieve tower objects from the emulator tree
	for(uint towIt = 0; towIt < nTowemu; ++towIt){
	    towEtemu  = l1Towemu_->iet[towIt];
	    towHademu  = l1Towemu_->ihad[towIt];
	    towEmemu  = l1Towemu_->iem[towIt];
	    towEtaemu = l1Towemu_->ieta[towIt];
	    towPhiemu = l1Towemu_->iphi[towIt];
	    if (towEtemu > topTowerEt){
		topTowerIPhi = towPhiemu;
		topTowerIEta = towEtaemu;
		topTowerEt = towEtemu;
		topTowerHad = towHademu;
		topTowerEm = towEmemu;
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

	    hAllTowEt->Fill(towEtemu);
	    hAllTowEta->Fill(towEtaemu);
	    hTowTPETEta->Fill(towEtaemu, towEtemu);
	    for (int iTopTowerIEta = -4; iTopTowerIEta <= 4; ++iTopTowerIEta){
		for (int iTopTowerIPhi = -4; iTopTowerIPhi <= 4; ++iTopTowerIPhi){
		    if (towEtaemu == topTowerIEta+iTopTowerIEta && towPhiemu == topTowerIPhi+iTopTowerIPhi){
			topTower9x9Em += towEmemu;
			topTower9x9Had += towHademu;
			if (abs(iTopTowerIPhi) <= 1 && abs(iTopTowerIEta) <= 1){
			    topTower3x3Em += towEmemu;
			    topTower3x3Had += towHademu;
			}
		    }
		}
	    }
	    for (int iSeedIEta = -4; iSeedIEta <= 4; ++iSeedIEta){
		for (int iSeedIPhi = -4; iSeedIPhi <= 4; ++iSeedIPhi){
		    if (towEtaemu == seedTowerIEta+iSeedIEta && towPhiemu == seedTowerIPhi+iSeedIPhi){
			seed9x9Had += towHademu;
			seed9x9Em += towEmemu;
			if (abs(iSeedIPhi) <= 1 && abs(iSeedIEta) <= 1){
			    seed3x3Had += towHademu;
			    seed3x3Em += towEmemu;
			}
		    }
		}
	    }
	    if (towEtaemu == seedTowerIEta && towPhiemu == seedTowerIPhi){
		seedEt = towEtemu;
		seedHad = towHademu;
		seedEm = towEmemu;
	    }
	}
	if (decayInHCAL){
	    hMaxJetSeedEtDecayInHcal->Fill(seedEt);
	    hMaxJetSeedHadDecayInHcal->Fill(seedHad);
	    hMaxJetSeedEmDecayInHcal->Fill(seedEm);
	    if (seedHad > 30){
		if (seedEm == 0){
		    hMaxJetSeedHOvEDecayInHcal->Fill(10);
		}
		else{
		    hMaxJetSeedHOvEDecayInHcal->Fill(TMath::Log10(seedHad/seedEm));
		}
		if (seed3x3Em == 0){
		    hMaxJetSeedHOvE3x3DecayInHcal->Fill(10);
		}
		else{
		    hMaxJetSeedHOvE3x3DecayInHcal->Fill(TMath::Log10(seedHad/seed3x3Em));
		}
		if (seed9x9Em == 0){
		    hMaxJetSeedHOvE9x9DecayInHcal->Fill(10);
		}
		else{
		    hMaxJetSeedHOvE9x9DecayInHcal->Fill(TMath::Log10(seedHad/seed9x9Em));
		}
	    }
	    hMaxTopTowerEtDecayInHcal->Fill(topTowerEt);
	    hMaxTopTowerHadDecayInHcal->Fill(topTowerHad);
	    hMaxTopTowerEmDecayInHcal->Fill(topTowerEm);
	    // std::cout << topTowerEm << " "<< topTower3x3Em << " " << topTower9x9Em << std::endl;
	    if (topTowerHad > 30){
		if (topTowerEm == 0){
		    hMaxTopTowerHOvEDecayInHcal->Fill(10);
		}
		else{
		    hMaxTopTowerHOvEDecayInHcal->Fill(TMath::Log10(topTowerHad/topTowerEm));
		}
		if (topTower3x3Em == 0){
		    hMaxTopTowerHOvE3x3DecayInHcal->Fill(10);
		}
		else{
		    hMaxTopTowerHOvE3x3DecayInHcal->Fill(TMath::Log10(topTowerHad/topTower3x3Em));
		}
		if (topTower9x9Em == 0){
		    hMaxTopTowerHOvE9x9DecayInHcal->Fill(10);
		}
		else{
		    hMaxTopTowerHOvE9x9DecayInHcal->Fill(TMath::Log10(topTowerHad/topTower9x9Em));
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

    formatPlot1D(hMaxJetEtDecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxJetEtDecayInHcal.pdf");

    formatPlot1D(hMaxJetSeedEtDecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxJetSeedEtDecayInHcal.pdf");

    formatPlot1D(hMaxJetSeedHadDecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxJetSeedHadDecayInHcal.pdf");

    formatPlot1D(hMaxJetSeedEmDecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxJetSeedEmDecayInHcal.pdf");

    formatPlot1D(hMaxJetSeedHOvEDecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxJetSeedHOvEDecayInHcal.pdf");

    formatPlot1D(hMaxJetSeedHOvE3x3DecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxJetSeedHOvE3x3DecayInHcal.pdf");

    formatPlot1D(hMaxJetSeedHOvE9x9DecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxJetSeedHOvE9x9DecayInHcal.pdf");

    formatPlot1D(hMaxTopTowerEtDecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxTopTowerEtDecayInHcal.pdf");

    formatPlot1D(hMaxTopTowerHadDecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxTopTowerHadDecayInHcal.pdf");

    formatPlot1D(hMaxTopTowerEmDecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxTopTowerEmDecayInHcal.pdf");

    formatPlot1D(hMaxTopTowerHOvEDecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxTopTowerHOvEDecayInHcal.pdf");

    formatPlot1D(hMaxTopTowerHOvE3x3DecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxTopTowerHOvE3x3DecayInHcal.pdf");

    formatPlot1D(hMaxTopTowerHOvE9x9DecayInHcal, jetColour);
    canvas->SaveAs("./Plots/Jets/maxTopTowerHOvE9x9DecayInHcal.pdf");
}
