
// -*- C++ -*-
//
// Package:    DisplacedJets/TrackAnalyzer
// Class:      TrackAnalyzer
// 
/**\class TrackAnalyzer TrackAnalyzer.cc TrackAnalyzer/TrackAnalyzer/plugins/TrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Robert Hardenbrook
//         Created:  Sat, 15 Nov 2014 15:18:18 GMT
//
//

#define DEBUG

//root includes
#include "TFile.h"  
#include "TH2F.h"   
#include "TH1F.h"   
#include "TTree.h"                      

// system include files                                                                                                                                                 
#include <vector> 
#include <string> 
#include <cmath>
#include <memory>
#include <iostream>
#include <algorithm>

// user include files                                                                                                                                                   
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//geometry
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

//formats
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPData.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

//mesages
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//parameters
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class TrackAnalyzer : public edm::EDAnalyzer {

public:
  explicit TrackAnalyzer(const edm::ParameterSet&);
  ~TrackAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:

  //methods
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //file configuration tags
  std::string outputFileName_;
  TFile *outputFile_;
  bool isMC_;

  //tracking tags
  edm::InputTag tag_generalTracks_;
  edm::InputTag tag_trackIPTagInfoCollection_;
  edm::EDGetTokenT<reco::TrackIPTagInfoCollection> token_trackIPTagInfoCollection_;
  
  //jet tags
  edm::InputTag tag_ak4CaloJets_;
  edm::InputTag tag_ak5CaloJets_;
  
  //gen info
  edm::InputTag tag_genParticles_;
  edm::InputTag tag_ak5GenJets_;
  edm::InputTag tag_genMetCalo_;

  //cuts
  Double_t cut_jetPt;
  Double_t cut_jetEta;
  
  //output related
  TTree* trackTree_;   

  static const Int_t MAX_TRACKS = 9999;
  static const Int_t MAX_JETS = 999;
  
  //tree variables
  Int_t nTracks = MAX_TRACKS;
  Int_t nCaloJets = MAX_JETS;

  // track kinematics
  Float_t trackPt[MAX_TRACKS];
  Float_t trackEta[MAX_TRACKS];
  Float_t trackPhi[MAX_TRACKS];

  // jet kinematics
  Float_t caloJetPt[MAX_JETS];
  Float_t caloJetEta[MAX_JETS];
  Float_t caloJetPhi[MAX_JETS];

  // size
  Float_t caloJetN90[MAX_JETS];
  Float_t caloJetN60[MAX_JETS];
  Float_t caloJetTowerArea[MAX_JETS];

  // energy contribution
  Float_t caloJetHfrac[MAX_JETS];
  Float_t caloJetEfrac[MAX_JETS];


  //output histograms

  
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig)
{
  //output configuration
  outputFileName_    =   iConfig.getUntrackedParameter<std::string>("outputFileName");
  isMC_    =   iConfig.getUntrackedParameter<bool>("isMC");

  //tags
  tag_generalTracks_ = iConfig.getUntrackedParameter<edm::InputTag>("generalTracks");
  tag_ak5CaloJets_ = iConfig.getUntrackedParameter<edm::InputTag>("ak5CaloJets");
  //  tag_trackIPTagInfoCollection_ = iConfig.getUntrackedParameter<edm::InputTag>("trackIPTagInfoCollection");
  token_trackIPTagInfoCollection_ = consumes<reco::TrackIPTagInfoCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackIPTagInfoCollection"));


  //cuts 
  cut_jetPt = iConfig.getUntrackedParameter<double>("jetPt");
  cut_jetEta = iConfig.getUntrackedParameter<double>("jetEta");


  //mc tags
  if(isMC_) {
    tag_ak5GenJets_ = iConfig.getUntrackedParameter<edm::InputTag>("ak5GenJets");
    tag_genMetCalo_ = iConfig.getUntrackedParameter<edm::InputTag>("genMetCalo");
    tag_genParticles_ = iConfig.getUntrackedParameter<edm::InputTag>("genParticles");
  }

}


TrackAnalyzer::~TrackAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void 
TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   //////////////////////
   // Extract Collections 
   //////////////////////

   // AOD Compatible
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(tag_generalTracks_, tracks); 
  
  edm::Handle<reco::CaloJetCollection> ak5CaloJets;
  iEvent.getByLabel(tag_ak5CaloJets_, ak5CaloJets);

  edm::Handle<reco::TrackIPTagInfoCollection> trackIPTagInfoCollection;
  iEvent.getByToken(token_trackIPTagInfoCollection_, trackIPTagInfoCollection);
  //iEvent.getByLabel(tag_trackIPTagInfoCollection_, trackIPTagInfoCollection);


  //SIM Compatible 
  if(isMC_) {
    edm::Handle<reco::GenParticleCollection > genParticle;
    iEvent.getByLabel(tag_genParticles_, genParticle);
    
    // edm::Handle<std::vector<int> > genParticleID;
    // iEvent.getByLabel(tag_genParticle_, genParticleID);
    
    // edm::Handle<std::vector<float> > genMetCalo;
    // iEvent.getByLabel(tag_genMetCalo_, genMetCalo);
    
    // edm::Handle<std::vector<float> > ak5GenJets;
    // iEvent.getByLabel(tag_ak5GenJets_, ak5GenJets);
  }
  

  // RECO Compatible
  
  // RAW Compatible

  // All Formats Compatible
  
  
  //////////////////////////////////
  // Calculate Variables
  //////////////////////////////////
  

  /////////////////////////////////
  // Fill Trees
  /////////////////////////////////
  
  nTracks = 0;
  nCaloJets = 0;

  //
  Int_t jj = 0;
  for(reco::CaloJetCollection::const_iterator jet = ak5CaloJets->begin(); jet != ak5CaloJets->end(); ++jet, jj++){

    //cuts 
    if (jet->pt() < cut_jetPt) continue;
    if (fabs(jet->eta()) > cut_jetEta) continue;

    caloJetPt[jj] = jet->pt();
    caloJetEta[jj] = jet->eta();
    caloJetPhi[jj] = jet->phi();        
    
    //area 
    //caloJetN90[jj] = jet->n90();
    //    caloJetN60[jj] = jet->n60();
    caloJetTowerArea[jj] = jet->towersArea();        

    //energy
    caloJetHfrac[jj] = jet->energyFractionHadronic();
    caloJetEfrac[jj] = jet->emEnergyFraction();

    nCaloJets++;
  }



#ifdef DEBUG
  std::cout << "[DEBUG] Begin Extracting IP Info" << std::endl;
#endif

  //extract the collection from the handle
  const reco::TrackIPTagInfoCollection & ip = *(trackIPTagInfoCollection.product()); 
  std::cout << "Found " << ip.size() << " TagInfo" << std::endl;

  //iterate over the ip info
  reco::TrackIPTagInfoCollection::const_iterator ipinfo = ip.begin(); 
  for(; ipinfo != ip.end(); ++ipinfo){    

    std::cout << "Jet pt: " << ipinfo->jet()->pt() << std::endl;
    std::cout << "Tot tracks: " << ipinfo->tracks().size() << std::endl;    

    reco::TrackRefVector selTracks = ipinfo->selectedTracks();
    int n_seltracks = selTracks.size();

    GlobalPoint pv(ipinfo->primaryVertex()->position().x(),ipinfo->primaryVertex()->position().y(),ipinfo->primaryVertex()->position().z());
    //    std::cout << pv << " vs " << ipinfo->primaryVertex()->position()   << std::endl;

#ifdef DEBUG
      std::cout << "[DEBUG] N Selected Tracks: " << n_seltracks << std::endl; 
      std::cout << "[DEBUG] Extracting Significances from Info" << std::endl; 
#endif 

    for(int tt = 0; tt < n_seltracks; tt++){      
      
      reco::btag::TrackIPData data = ipinfo->impactParameterData()[tt];  
      
      std::cout << selTracks[tt]->pt() << "\t";
      std::cout << ipinfo->probabilities(0)[tt] << "\t";
      std::cout << ipinfo->probabilities(1)[tt] << "\t";
      std::cout << data.ip3d.value() << "\t";
      std::cout << data.ip3d.significance() << "\t";
      std::cout << data.distanceToJetAxis.value() << "\t";
      std::cout << data.distanceToJetAxis.significance() << "\t";
      std::cout << data.distanceToGhostTrack.value() << "\t";
      std::cout << data.distanceToGhostTrack.significance() << "\t";
      std::cout << data.closestToJetAxis << "\t";
      std::cout << (data.closestToJetAxis - pv).mag() << "\t";
      std::cout << data.closestToGhostTrack << "\t";
      std::cout << (data.closestToGhostTrack - pv).mag() << "\t";
      std::cout << data.ip2d.value() << "\t";
      std::cout << data.ip2d.significance() <<  std::endl;     

    }
  }

  trackTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackAnalyzer::beginJob()
{

#ifdef DEBUG
  std::cout << "[DEBUG] Setting Up Output File And Tree" << std::endl;
#endif 

  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
  trackTree_  = new TTree("tree","track tree");

  trackTree_->Branch("nCaloJets", &nCaloJets, "nCaloJets/I");
  trackTree_->Branch("nTracks", &nTracks, "nTracks/I");

  trackTree_->Branch("caloJetPt", &caloJetPt, "caloJetPt[nCaloJets]/F");
  trackTree_->Branch("caloJetPhi", &caloJetPhi, "caloJetPhi[nCaloJets]/F");
  trackTree_->Branch("caloJetEta", &caloJetEta, "caloJetEta[nCaloJets]/F");

  trackTree_->Branch("caloJetn90", &caloJetN90, "caloJetN90[nCaloJets]/F");
  trackTree_->Branch("caloJetn60", &caloJetN60, "caloJetN60[nCaloJets]/F");
  trackTree_->Branch("caloJetTowerArea", &caloJetTowerArea, "caloJetTowerArea[nCaloJets]/F");
  trackTree_->Branch("caloJetHfrac", &caloJetHfrac, "caloJetHfrac[nCaloJets]/F");
  trackTree_->Branch("caloJetEfrac", &caloJetEfrac, "caloJetEFrac[nCaloJets]/F");



  trackTree_->Branch("trackPt", &trackPt, "trackPt[nTracks]/F");
  trackTree_->Branch("trackPhi", &trackPhi, "trackPhi[nTracks]/F");
  trackTree_->Branch("trackEta", &trackEta, "trackEta[nTracks]/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackAnalyzer::endJob() 
{
  trackTree_->Write();
  outputFile_->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
TrackAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
TrackAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TrackAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TrackAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
