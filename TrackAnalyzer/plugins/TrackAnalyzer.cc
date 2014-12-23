
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

//#define DEBUG

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
#include <assert.h> 

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
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/GeometrySurface/interface/Line.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPData.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "DataFormats/BTauReco/interface/VertexTypes.h"

//mesages
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//parameters
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//kineamtics
#include <RecoBTag/SecondaryVertex/interface/TrackKinematics.h>

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
  virtual Float_t getJetMedian(Float_t daArray[], int size, int jetid, bool is_signed);
  virtual Float_t getJetVariance(Float_t values[], Float_t mean, int size, int jetid, bool is_signed);
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //file configuration tags
  std::string outputFileName_;
  TFile *outputFile_;
  bool isMC_;
  bool isSignalMC_;

  //tracking tags
  edm::InputTag tag_generalTracks_;
  edm::InputTag tag_trackIPTagInfoCollection_;
  edm::InputTag tag_lifetimeIPTagInfo_;
  edm::InputTag tag_secondaryVertexTagInfo_; 

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
  TTree* jetTree_;   

  static const Int_t MAX_TRACKS = 9999;
  static const Int_t MAX_JETS = 999;
  static const Int_t MAX_VTX = 999;
  static const Int_t debug = 0; 

#ifdef DEBUG
  debug = 2; 
#endif 
 
  //bookkeeping
  int evNum = 0;
  Int_t jetid = 0;

  //tree variables
  Int_t nLiTracks = 0;
  Int_t nCaloJets = 0;

  ////////////////// CALO JETS ////////////

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

  ///////////////////// GEN MATCHED ////////////////////

  Int_t genMatch[MAX_JETS];
  Int_t genMatchTrack[MAX_TRACKS];

  Float_t genPt[MAX_JETS];
  Float_t genEta[MAX_JETS];
  Float_t genPhi[MAX_JETS];
  Float_t genM[MAX_JETS];

  //////////////////// JET TAG /////////////////

  Int_t nTracks = 0;

  //track infoTags
  Float_t jetTrackDR[MAX_TRACKS];
  Float_t trackIP3D[MAX_TRACKS];
  Float_t trackIP2D[MAX_TRACKS];
  Float_t trackIPSig3D[MAX_TRACKS];
  Float_t trackIPSig2D[MAX_TRACKS];
  Float_t trackDistanceJetAxis[MAX_TRACKS];
  Float_t trackDistanceJetAxisSig[MAX_TRACKS];
  Float_t trackChi2[MAX_TRACKS];
  Int_t trackJetID[MAX_TRACKS];
  
  // track kinematics
  Float_t trackPt[MAX_TRACKS];
  Float_t trackEta[MAX_TRACKS];
  Float_t trackPhi[MAX_TRACKS];

  // tag jets kinematics
  Int_t nTagJets = 0;
  Int_t jetJetID[MAX_JETS];
  Float_t tagJetPt[MAX_JETS];
  Float_t tagJetEta[MAX_JETS];
  Float_t tagJetPhi[MAX_JETS];
  // size
  Float_t tagJetN90[MAX_JETS];
  Float_t tagJetN60[MAX_JETS];
  Float_t tagJetTowerArea[MAX_JETS];

  // energy contribution
  Float_t tagJetHfrac[MAX_JETS];
  Float_t tagJetEfrac[MAX_JETS];
  Int_t tagJetNSelTracks[MAX_JETS];

  //////////////////// LIFETIME TAG /////////////////

  // track kinematics
  Float_t liTrackPt[MAX_TRACKS];
  Float_t liTrackEta[MAX_TRACKS];
  Float_t liTrackPhi[MAX_TRACKS];
  
  // lifetime jet kinematics
  Int_t nLiJets = 0;
  Int_t liJetID[MAX_JETS];
  Float_t liJetPt[MAX_JETS];
  Float_t liJetEta[MAX_JETS];
  Float_t liJetPhi[MAX_JETS];
  Int_t liJetNSelTracks[MAX_JETS];

  // lifetime track infoTags
  Float_t liTrackIP3D[MAX_TRACKS];
  Float_t liTrackIP2D[MAX_TRACKS];
  Float_t liTrackIPSig3D[MAX_TRACKS];
  Float_t liTrackIPSig2D[MAX_TRACKS];
  Float_t liTrackDistanceJetAxis[MAX_TRACKS];
  Float_t liTrackDistanceJetAxisSig[MAX_TRACKS];
  Float_t liTrackChi2[MAX_TRACKS];
  
  //tracks to jet comparisons
  Int_t liTrackJetID[MAX_TRACKS];
  Float_t liJetTrackDR[MAX_TRACKS];

  //////////////////// SV TAG ////////////////

  Int_t nSV = 0;

  Int_t svVertexJetID[MAX_VTX]; 
  Int_t nSvJets = 0;
  Float_t svJetPt[MAX_JETS];
  Float_t svJetEta[MAX_JETS];
  Float_t svJetPhi[MAX_JETS];
  Float_t svJetID[MAX_JETS];

  //kinematics
  Float_t svMass[MAX_VTX];
  Float_t svPx[MAX_VTX];
  Float_t svPy[MAX_VTX];
  Float_t svPt[MAX_VTX];
  Float_t svEta[MAX_VTX];
  Float_t svPhi[MAX_VTX];

  //position
  Float_t svX[MAX_VTX];
  Float_t svY[MAX_VTX];
  Float_t svZ[MAX_VTX];

  //position error
  Float_t svXErr[MAX_VTX];
  Float_t svYErr[MAX_VTX];
  Float_t svZErr[MAX_VTX];

  //quality
  Float_t svChi2[MAX_VTX];
  Float_t svNChi2[MAX_VTX];
  Float_t svNDof[MAX_VTX];
  Int_t svIsValid[MAX_VTX];
  
  //tracking
  Int_t svNTracks[MAX_VTX];
  //  Int_t svTrackVertexID[MAX_VTX];   
  Float_t svTotalCharge[MAX_VTX];

  // Flight Information
  Float_t svFlight[MAX_VTX];
  Float_t svFlightErr[MAX_VTX];
  Float_t svFlight2D[MAX_VTX];
  Float_t svFlight2DErr[MAX_VTX];

  // various DR 
  Float_t svDRFlightJet[MAX_VTX];
  Float_t svDRTrackJet[MAX_VTX];
  Float_t svDRTrackFlight[MAX_VTX];

  ///////////////// JET TREE SPECIFIC MEMBERS /////////////////

  //book keeping
  Int_t nJetWithSv; 

  // significance and aboslute IP weighted track energy
  Float_t jetEIPSig2D[MAX_JETS];
  Float_t jetEIP2D[MAX_JETS];
  Float_t jetEIPSig3D[MAX_JETS];
  Float_t jetEIP3D[MAX_JETS];

  // significance log IP weighted track energy
  Float_t jetELogIPSig2D[MAX_JETS];
  // Float_t jetELogIP2D[MAX_JETS];
  Float_t jetELogIPSig3D[MAX_JETS];
  //Float_t jetELogIP3D[MAX_JETS];
  
  // signed weighted pt 
  Float_t jetEIPSignedSig2D[MAX_JETS];
  Float_t jetEIPSigned2D[MAX_JETS];
  Float_t jetEIPSignedSig3D[MAX_JETS];
  Float_t jetEIPSigned3D[MAX_JETS];

  // absolute IP sums
  Float_t jetIPSum2D[MAX_JETS];
  Float_t jetIPSum3D[MAX_JETS];
  Float_t jetIPSignedSum2D[MAX_JETS];
  Float_t jetIPSignedSum3D[MAX_JETS];

  // IP significance sums
  Float_t jetIPSigSum2D[MAX_JETS];
  Float_t jetIPSigSum3D[MAX_JETS];
  Float_t jetIPSigInvSum2D[MAX_JETS];
  Float_t jetIPSigInvSum3D[MAX_JETS];

  // IP significance log sums
  Float_t jetIPSigLogSum2D[MAX_JETS];
  Float_t jetIPSigLogSum3D[MAX_JETS];

  // IP signed significance sums
  Float_t jetIPSignedSigSum2D[MAX_JETS];
  Float_t jetIPSignedSigSum3D[MAX_JETS];
  Float_t jetIPSignedSigInvSum2D[MAX_JETS];
  Float_t jetIPSignedSigInvSum3D[MAX_JETS];

  // IP significance averages
  Float_t jetMeanIPSig2D[MAX_JETS];
  Float_t jetMeanIPSig3D[MAX_JETS];
  Float_t jetMedianIPSig2D[MAX_JETS];
  Float_t jetMedianIPSig3D[MAX_JETS];
  Float_t jetVarianceIPSig2D[MAX_JETS];
  Float_t jetVarianceIPSig3D[MAX_JETS];

  // signed averages 
  Float_t jetMeanIPSignedSig2D[MAX_JETS];
  Float_t jetMeanIPSignedSig3D[MAX_JETS];
  Float_t jetMedianIPSignedSig2D[MAX_JETS];
  Float_t jetMedianIPSignedSig3D[MAX_JETS];
  Float_t jetVarianceIPSignedSig2D[MAX_JETS];
  Float_t jetVarianceIPSignedSig3D[MAX_JETS];

  // Absolute IP averages
  Float_t jetMeanIP2D[MAX_JETS];
  Float_t jetMeanIP3D[MAX_JETS];
  Float_t jetMedianIP2D[MAX_JETS];
  Float_t jetMedianIP3D[MAX_JETS];

  // Absolute Signed IP averages
  Float_t jetMeanIPSigned2D[MAX_JETS];
  Float_t jetMeanIPSigned3D[MAX_JETS];
  Float_t jetMedianIPSigned2D[MAX_JETS];
  Float_t jetMedianIPSigned3D[MAX_JETS];

  // SV information
  Int_t jetNSv[MAX_JETS];
  Float_t jetSvMass[MAX_JETS];
  Float_t jetSvLxy[MAX_JETS];
  Float_t jetSvLxySig[MAX_JETS];
  Float_t jetSvLxyz[MAX_JETS];
  Float_t jetSvLxyzSig[MAX_JETS];
  Int_t jetSvNTrack[MAX_JETS]; //vertex track multiplicty

  // SV position
  Float_t jetSvX[MAX_JETS];
  Float_t jetSvY[MAX_JETS];
  Float_t jetSvZ[MAX_JETS];

  Float_t jetSvZErr[MAX_JETS];
  Float_t jetSvYErr[MAX_JETS];
  Float_t jetSvXErr[MAX_JETS];

  //SV quality
  Float_t jetSvChi2[MAX_JETS];
  Float_t jetSvNChi2[MAX_JETS];
  Float_t jetSvNDof[MAX_JETS];
  Int_t jetSvIsValid[MAX_JETS];

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
  isSignalMC_    =   iConfig.getUntrackedParameter<bool>("isSignalMC");

  //tags
  tag_generalTracks_ = iConfig.getUntrackedParameter<edm::InputTag>("generalTracks");
  tag_ak5CaloJets_ = iConfig.getUntrackedParameter<edm::InputTag>("ak5CaloJets");
  tag_secondaryVertexTagInfo_ = iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertexTagInfo");  
  tag_lifetimeIPTagInfo_ = iConfig.getUntrackedParameter<edm::InputTag>("lifetimeIPTagInfo");  

//  tag_trackIPTagInfoCollection_ = iConfig.getUntrackedParameter<edm::InputTag>("trackIPTagInfoCollection");
  token_trackIPTagInfoCollection_ = consumes<reco::TrackIPTagInfoCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackIPTagInfoCollection"));
  
  
  //cuts 
  cut_jetPt = iConfig.getUntrackedParameter<double>("jetPt");
  cut_jetEta = iConfig.getUntrackedParameter<double>("jetEta");


  //mc tags
  if(isMC_) {
    //tag_ak5GenJets_ = iConfig.getUntrackedParameter<edm::InputTag>("ak5GenJets");
    //tag_genMetCalo_ = iConfig.getUntrackedParameter<edm::InputTag>("genMetCalo");
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

  edm::Handle<reco::TrackIPTagInfoCollection> lifetimeIPTagInfo;
  iEvent.getByLabel(tag_lifetimeIPTagInfo_, lifetimeIPTagInfo);

  edm::Handle<reco::SecondaryVertexTagInfoCollection> secondaryVertexTagInfo;
  iEvent.getByLabel(tag_secondaryVertexTagInfo_, secondaryVertexTagInfo);  

  //SIM Compatible 
  if(isMC_) {
    //edm::Handle<reco::GenParticleCollection > genParticle;
    //iEvent.getByLabel(tag_genParticles_, genParticle);
    
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
     
  //
  Int_t jj = 0;
  nCaloJets = 0;
  for(reco::CaloJetCollection::const_iterator jet = ak5CaloJets->begin(); jet != ak5CaloJets->end(); ++jet, jj++){

    //cuts 
    if (jet->pt() < cut_jetPt) continue;
    //if (fabs(jet->eta()) > cut_jetEta) continue;
    
    // initialize to zero in case there isnt any match
    genMatch[jj] = 0; 
    genPt[jj] = 0;
    genEta[jj] = 0;
    genPhi[jj] = 0;
    genM[jj] = 0;

    caloJetPt[jj] = jet->pt();
    caloJetEta[jj] = jet->eta();
    caloJetPhi[jj] = jet->phi();        
    
    //area 
    //caloJetN90[jj] = jet->n90();
    //caloJetN60[jj] = jet->n60();
    caloJetTowerArea[jj] = jet->towersArea();        

    //energy
    caloJetHfrac[jj] = jet->energyFractionHadronic();
    caloJetEfrac[jj] = jet->emEnergyFraction();

    nCaloJets++;
  } // end loop over calojets

  edm::Handle<reco::GenParticleCollection> genParticles;

  // generator matching to truely displaced jets
  if (isSignalMC_) {
    iEvent.getByLabel("genParticles", genParticles);    

    for(size_t pp = 0; pp < genParticles->size(); ++pp) {
      const reco::GenParticle & part = (*genParticles)[pp];
      int id = part.pdgId();
      int st = part.status();  

      if (st != 3 || fabs(id) > 6 ) continue;     

      // const reco::Candidate * mom = part.mother();
      double genpt = part.pt(), geneta = part.eta(), genphi = part.phi(), genmass = part.mass();
      
      //      if (debug > 1 ) 
      //      std::cout << "[GEN CAND] id " << id << " status " << st << " pt " << genpt << " eta " << geneta << " phi " << genphi  << std::endl;      

      // check each jet if it matches
      for(int jj = 0; jj < nCaloJets; jj++ ) {

	float calopt = caloJetPt[jj],  caloeta = caloJetEta[jj], calophi = caloJetPhi[jj];

	//	if (calopt > 100) {std::cout << "[CALO JET] pt " << calopt << " eta " << caloeta << " phi " << calophi << std::endl;}

	float dr = reco::deltaR( geneta, genphi, caloeta, calophi);
	float dpt = fabs(calopt - genpt) / genpt;

	// found a match
	if (dr < .5 && dpt < .5) {
	  std::cout << "[GEN MATCHED] id " << id << " status " << st << " pt " << genpt << " eta " << geneta  <<  " phi " << genphi << std::endl;      
	  genMatch[jj]++; 
	  genPt[jj] = genpt;
	  genEta[jj] = geneta;
	  genPhi[jj] = geneta;
	  genM[jj] = genmass;
	}
      } //end loop over caloejts
    }
  }

  if(debug > 1 ) std::cout << "[DEBUG] Begin Extracting Tag Info" << std::endl;

  //extract the collection from the handle
  const reco::TrackIPTagInfoCollection & ip = *(trackIPTagInfoCollection.product()); 
  const reco::TrackIPTagInfoCollection & lifetime = *(lifetimeIPTagInfo.product()); 
  const reco::SecondaryVertexTagInfoCollection & sv = *(secondaryVertexTagInfo.product()); 

  //  std::auto_ptr<reco::VertexCollection> displacedSvVertexCollectionResult(new reco::VertexCollection);
  //  reco::VertexCollection displacedSvVertexCollection;
  if (debug > 1) {
    std::cout << "-- Found " << ip.size() << " IP TagInfo" << std::endl;
    std::cout << "-- Found " << lifetime.size() << " Lifetime TagInfo" << std::endl;
    std::cout << "-- Found " << sv.size() << " Secondary Vertex TagInfo" << std::endl;
  }

  nTracks = 0;
  nTagJets = 0;

  //iterate over the ip info
  reco::TrackIPTagInfoCollection::const_iterator ipinfo = ip.begin(); 
  jj = 0;
  for(; ipinfo != ip.end(); ++ipinfo, jj++){    
    
    // redudancey for tree output
    jetJetID[nTagJets] = jetid; 

    //skip low pt jets
    if (ipinfo->jet()->pt() < cut_jetPt) continue;


    if(debug > 1){
      std::cout << "JET ID: " << jj << std::endl;
      std::cout << "Jet pt: " << ipinfo->jet()->pt() << std::endl;
      std::cout << "Tot tracks: " << ipinfo->tracks().size() << std::endl;    
    }

    
    //kinematics
    tagJetPt[nTagJets] = ipinfo->jet()->pt();
    tagJetEta[nTagJets] = ipinfo->jet()->eta();
    tagJetPhi[nTagJets] = ipinfo->jet()->phi();        

    //area 
    // tagJetN90[jj] = jet->n90();
    // tagJetN60[jj] = jet->n60();
    // tagJetTowerArea[jj] = ipinfo->jet()->towersArea();        

    // energy type
    //tagJetHfrac[jj] = ipinfo->jet()->energyFractionHadronic();
    //tagJetEfrac[jj] = ipinfo->jet()->emEnergyFraction();

    reco::TrackRefVector selTracks = ipinfo->selectedTracks();
    tagJetNSelTracks[nTagJets] = selTracks.size();
    

    if(debug > 1){
      std::cout << "[DEBUG] N Selected Tracks: " << tagJetNSelTracks[nTagJets] << std::endl; 
      std::cout << "[DEBUG] Extracting Significances from Info" << std::endl; 
    }


   // Loop over the tracks and fill tree information
    
    for(int tt = 0; tt < tagJetNSelTracks[nTagJets]; tt++){      
      
      reco::btag::TrackIPData data = ipinfo->impactParameterData()[tt];  
      const reco::Track ptrack = *(selTracks[tt]);
      
      if (ptrack.pt() < .01) {
	continue;  //minimal track pt
      }
            
      //kinematics
      trackPt[nTracks] = ptrack.pt();
      trackPhi[nTracks] = ptrack.phi();
      trackEta[nTracks] = ptrack.eta();

      //jet track association  
      trackJetID[nTracks] = jetid;
      jetTrackDR[nTracks] = reco::deltaR( ptrack.eta(), ptrack.phi(),
				   ipinfo->jet()->eta(), ipinfo->jet()->phi());

      //ip info tags
      trackIP3D[nTracks] = data.ip3d.value();
      trackIPSig3D[nTracks] = data.ip3d.significance();
      trackIP2D[nTracks] = data.ip2d.value();
      trackIPSig2D[nTracks] = data.ip2d.significance();
      trackDistanceJetAxis[nTracks] = data.distanceToJetAxis.value();
      trackDistanceJetAxisSig[nTracks] = data.distanceToJetAxis.significance();

      //track quality
      trackChi2[nTracks] = ptrack.normalizedChi2();

      nTracks++; // number of tracks in event
    } // End Matched Track Loop

    nTagJets++; // index for the array being filed
    jetid++; // global id for the jet
  } // End Jet Loop (loop over ip tags)


  //iterate over the lifetime ip info
  reco::TrackIPTagInfoCollection::const_iterator liinfo = lifetime.begin(); 

  jj = 0;
  jetid -= nTagJets; //reset the jetid counter
  nLiJets = 0;
  nLiTracks = 0;
  
  for(; liinfo != lifetime.end(); ++liinfo, jj++){        
    //skip low pt jets
    if (liinfo->jet()->pt() < cut_jetPt) continue;    

    reco::TrackRefVector liTracks = liinfo->selectedTracks();    

    // track selection  
    int n_liTracks = liTracks.size();

    liJetNSelTracks[nLiJets] = n_liTracks;

    // redundancy for tree output
    liJetID[nLiJets] = jetid; 

    //kinematics
    liJetPt[nLiJets] = liinfo->jet()->pt();
    liJetEta[nLiJets] = liinfo->jet()->eta();
    liJetPhi[nLiJets] = liinfo->jet()->phi();        

    bool is_matched = genMatch[nLiJets];

    for(int tt = 0; tt < n_liTracks; tt++){            
      reco::btag::TrackIPData data = liinfo->impactParameterData()[tt];  
      const reco::Track litrack = *(liTracks[tt]);
      
      if (litrack.pt() < .01) {
	continue;  //minimal track pt
      }
            
      //kinematics
      liTrackPt[nLiTracks] = litrack.pt();
      liTrackPhi[nLiTracks] = litrack.phi();
      liTrackEta[nLiTracks] = litrack.eta();

      //jet track association  
      liTrackJetID[nLiTracks] = jetid;
      liJetTrackDR[nLiTracks] = reco::deltaR( litrack.eta(), litrack.phi(),
				   liinfo->jet()->eta(), liinfo->jet()->phi());

      // mark the track as part of a gen matched jet
      genMatchTrack[nLiTracks] = 0;
      if (is_matched) genMatchTrack[nLiTracks] = 1;

      //ip info tags
      liTrackIP3D[nLiTracks] = data.ip3d.value();
      liTrackIPSig3D[nLiTracks] = data.ip3d.significance();
      liTrackIP2D[nLiTracks] = data.ip2d.value();
      liTrackIPSig2D[nLiTracks] = data.ip2d.significance();
      liTrackDistanceJetAxis[nLiTracks] = data.distanceToJetAxis.value();
      liTrackDistanceJetAxisSig[nLiTracks] = data.distanceToJetAxis.significance();

      //track quality
      liTrackChi2[nLiTracks] = litrack.normalizedChi2();

      nLiTracks++; // number of tracks in event
    } // End Matched Track Loop

    jetid++; //global jet id
    nLiJets++; //index in filled array
  }// end ip info jet loop

  //clear out the sv collection each event
  //displacedSVertexCollection->clear();

  reco::SecondaryVertexTagInfoCollection::const_iterator svinfo = sv.begin();
  jj = 0;
  jetid -= nLiJets; //reset the jetid counter
  nSvJets = 0;
  nSV = 0; // global event counter fo SV
  for(; svinfo != sv.end(); ++svinfo, jj++){    

    //skip low pt jets
    if (svinfo->jet()->pt() < cut_jetPt) continue;    

    reco::TrackRefVector svTracks = svinfo->selectedTracks();    

    // global jet identifier 
    svJetID[nSvJets] = jetid; 
    
    // jet kinematics
    svJetPt[nSvJets] = svinfo->jet()->pt();
    svJetEta[nSvJets] = svinfo->jet()->eta();
    svJetPhi[nSvJets] = svinfo->jet()->phi();        

    // number of secondary vertices 
    int nSvVertex = 0;
    int thisJetnSV = svinfo->nVertices();
    for(int vv = 0; vv < thisJetnSV; vv++) {

      if(debug > 1)      std::cout << "[Tracking] jet " << jetid << " nSV " << nSV << " nSvVertex " << nSvVertex << " thisJetnSV " << thisJetnSV << std::endl;

      svVertexJetID[nSV] = jetid;
      reco::Vertex svVertex = svinfo->secondaryVertex(vv);  

      // vertex kinematics
      svMass[nSV] = svinfo->secondaryVertex(vv).p4().mass();
      svPx[nSV] = svVertex.p4().px();
      svPy[nSV] = svVertex.p4().py();
      svPt[nSV] = svVertex.p4().pt();
      svEta[nSV] = svVertex.p4().eta();
      svPhi[nSV] = svVertex.p4().phi();
      
      // quality
      svChi2[nSV] = svVertex.chi2();
      svNChi2[nSV] = svVertex.normalizedChi2();
      svNDof[nSV] = svVertex.ndof();
      svIsValid[nSV] = svVertex.isValid(); 	

      // positional space
      svX[nSV] = svVertex.x();
      svXErr[nSV] = svVertex.xError();
      svY[nSV] = svVertex.y();   
      svYErr[nSV] = svVertex.yError();   
      svZ[nSV] = svVertex.z();      
      svZErr[nSV] = svVertex.zError();      
      
      // flight
      svFlight[nSV] = svinfo->flightDistance(vv).value();
      svFlightErr[nSV] = svinfo->flightDistance(vv).error();
      svFlight2D[nSV] = svinfo->flightDistance(vv, true).value();
      svFlight2DErr[nSV] = svinfo->flightDistance(vv, true).error();

      // charge
      svTotalCharge[nSV] = 0;

      //tracking 
      svNTracks[nSV] = svVertex.nTracks();
      reco::TrackKinematics vertexKinematics;

      Bool_t hasRefittedTracks = svVertex.hasRefittedTracks();      

      // loop over the tracks associated with the Secondary Vertex to compute total charge
      for(int tt = 0; tt < svNTracks[nSV]; tt++) {		
	//svTrackVertexID[tt] = vv; 		

	reco::Track track  = *svinfo->track(tt);//->get(); //get the track from TrackRef 
	Double_t trackWeight = svinfo->trackWeight(vv, tt);	

	// only keep high weight tracks
        if(trackWeight < 0.5){
	  continue;
	}

	// use the refitted tracks if they exist
        if (hasRefittedTracks) {
	  reco::Track refitTrack = svVertex.refittedTracks()[vv];

          vertexKinematics.add(refitTrack, trackWeight);
          svTotalCharge[nSV] += refitTrack.charge();
        }
        else {
          vertexKinematics.add(track, trackWeight);
          svTotalCharge[nSV] += track.charge();
        }
      } // end sv tracking look

      // vectors for deltaR calculation
      math::XYZTLorentzVector vertexSum = vertexKinematics.weightedVectorSum();
      math::XYZVector jetDir = svinfo->jet()->momentum().Unit();
      GlobalVector flightDir = svinfo->flightDirection(vv);

      // deltaR variables
      svDRFlightJet[nSV] = reco::deltaR(flightDir, jetDir);
      svDRTrackJet[nSV] = reco::deltaR(vertexSum, jetDir); 
      svDRTrackFlight[nSV] = reco::deltaR(vertexSum, flightDir);       
      
      nSvVertex++; // index in array corresponding to a single jet for vertex
      nSV++; // event sv index  (includes all jets)
    } // end sv vertex loop

    nSvJets++; //index in array for jet
    jetid++; //global jet identifier
  } // end sv tag (jet) loop
  

  /////////////////////////////// JET TREE CALCULATIONS ////////////////////////////

  // loop over all calo jets
  for(int jj = 0; jj < nCaloJets; jj++) {
    
    //double check that all the pt's match up between the tags
    bool pt_equal = caloJetPt[jj] == liJetPt[jj] && liJetPt[jj] == svJetPt[jj]; 
    bool eta_equal = caloJetEta[jj] == liJetEta[jj] && liJetEta[jj] == svJetEta[jj]; 
    bool phi_equal = caloJetPhi[jj] == liJetPhi[jj] && liJetPhi[jj] == svJetPhi[jj]; 
    bool all_match  = pt_equal && eta_equal && phi_equal;
    assert(all_match);
    
    // significance and aboslute IP weighted track energy
    jetEIPSig2D[jj] = 0;
    jetEIP2D[jj] = 0;
    jetEIPSig3D[jj] = 0;
    jetEIP3D[jj] = 0;

    //significance log weighted
    jetELogIPSig2D[jj] = 0;
    //jetELogIP2D[jj] = 0;
    jetELogIPSig3D[jj] = 0;
    //jetELogIP3D[jj] = 0;

    // signed versions of track weighting
    jetEIPSignedSig2D[jj] = 0;
    jetEIPSigned2D[jj] = 0;
    jetEIPSignedSig3D[jj] = 0;
    jetEIPSigned3D[jj] = 0;

    // IP significance sums
    jetIPSigSum2D[jj] = 0;
    jetIPSigSum3D[jj] = 0;
    jetIPSigInvSum2D[jj] = 0;
    jetIPSigInvSum3D[jj] = 0;

    // IP significance log sums
    jetIPSigLogSum2D[jj] = 0;
    jetIPSigLogSum3D[jj] = 0;

    // IP signed significance sums
    jetIPSignedSigSum2D[jj] = 0;
    jetIPSignedSigSum3D[jj] = 0;
    jetIPSignedSigInvSum2D[jj] = 0;
    jetIPSignedSigInvSum3D[jj] = 0;

    // IP significance averages
    jetMeanIPSig2D[jj] = 0;
    jetMeanIPSig3D[jj] = 0;
    jetMedianIPSig2D[jj] = 0;
    jetMedianIPSig3D[jj] = 0;
    jetVarianceIPSig2D[jj] = 0;
    jetVarianceIPSig3D[jj] = 0;

    // signed averages 
    jetMeanIPSignedSig2D[jj] = 0;
    jetMeanIPSignedSig3D[jj] = 0;
    jetMedianIPSignedSig2D[jj] = 0;
    jetMedianIPSignedSig3D[jj] = 0;
    jetVarianceIPSignedSig2D[jj] = 0;
    jetVarianceIPSignedSig3D[jj] = 0;

    //  ip tag related
    for(int tt = 0; tt < nLiTracks; tt++){            

      // only look at tracks that correspond to the current jet
      if( liTrackJetID[tt] != jetJetID[jj]) { continue; }

      float pt = liTrackPt[tt];
      float ip2d = liTrackIP2D[tt];
      float ip3d = liTrackIP3D[tt];
      float ip2ds = liTrackIPSig2D[tt];
      float ip3ds = liTrackIPSig3D[tt];

      // ip weighted track pt sums
      jetEIPSig2D[jj] += fabs(ip2ds) * pt;
      jetEIPSig3D[jj] += fabs(ip3ds) * pt;
      jetEIP2D[jj] += fabs(ip2d) * pt;
      jetEIP3D[jj] += fabs(ip3d) * pt;

      // ip log weighted track pt sums
      jetELogIPSig2D[jj] += (ip2ds ? log(fabs(ip2ds)) * pt : 0);
      jetELogIPSig3D[jj] += (ip3ds ? log(fabs(ip3ds)) * pt : 0);
      //jetELogIP2D[jj] += (ip2d ? log(fabs(ip2d)) * pt : 0);
      //jetELogIP3D[jj] += (ip3d ? log(fabs(ip3d)) * pt : 0);

      jetEIPSignedSig2D[jj] += (ip2ds) * pt;
      jetEIPSignedSig3D[jj] += (ip3ds) * pt;
      jetEIPSigned2D[jj] += (ip2d) * pt;
      jetEIPSigned3D[jj] += (ip3d) * pt; 

      // absolute ip averages
      // unsigned
      jetIPSum2D[jj] += fabs(ip2d);
      jetIPSum3D[jj] += fabs(ip3d);
      //signed
      jetIPSignedSum2D[jj] += ip2d;
      jetIPSignedSum3D[jj] += ip3d;

      // ip significance sums
      // unsigned 
      jetIPSigSum2D[jj] += fabs(ip2ds);
      jetIPSigSum3D[jj] += fabs(ip3ds);
      jetIPSigInvSum2D[jj] += (ip2ds ? 1.0 / fabs(ip2ds): 0);
      jetIPSigInvSum3D[jj] += (ip3ds ? 1.0 / fabs(ip3ds): 0);

      // signed
      jetIPSignedSigSum2D[jj] += (ip2ds);
      jetIPSignedSigSum3D[jj] += (ip3ds);
      jetIPSignedSigInvSum2D[jj] += 1.0 / (ip2ds);
      jetIPSignedSigInvSum3D[jj] += 1.0 / (ip3ds);

      // ip significance log sum
      jetIPSigLogSum2D[jj] += (ip2ds ? log(fabs(ip2ds)) : 0);
      jetIPSigLogSum3D[jj] += (ip3ds ? log(fabs(ip3ds)) : 0);

      // unsigned averages
      jetMeanIPSig2D[jj] += fabs(ip2ds) / float(liJetNSelTracks[jj]);
      jetMeanIPSig3D[jj] += fabs(ip3ds) / float(liJetNSelTracks[jj]);          

      // unsigned averages
      jetMeanIPSignedSig2D[jj] += ip2ds / float(liJetNSelTracks[jj]);
      jetMeanIPSignedSig3D[jj] += ip3ds / float(liJetNSelTracks[jj]);
    }

    bool has_track = liJetNSelTracks[jj] > 0;

    // do the average for the IP weighted
    jetEIPSig2D[jj] /= (has_track ? jetIPSigSum2D[jj] : 1);
    jetEIPSig3D[jj] /= (has_track ? jetIPSigSum3D[jj] : 1);
    jetEIP2D[jj] /= (has_track ?  jetIPSum2D[jj]: 1);
    jetEIP3D[jj] /= (has_track ?  jetIPSum3D[jj]: 1);

    // do the average for the IP weighted
    jetELogIPSig2D[jj] /= (has_track ? jetIPSigLogSum2D[jj] : 1);
    jetELogIPSig3D[jj] /= (has_track ? jetIPSigLogSum3D[jj] : 1);
    //jetELogIP2D[jj] /= (has_track ?  jetIPLogSum2D[jj]: 1);
    //jetELogIP3D[jj] /= (has_track ?  jetIPLogSum3D[jj]: 1);

    // averages for the signed values as well
    jetEIPSignedSig2D[jj] /= (has_track ?  jetIPSignedSum2D[jj]: 1);
    jetEIPSignedSig3D[jj] /= (has_track ?  jetIPSignedSum3D[jj]: 1);
    jetEIPSigned2D[jj] /= (has_track ?  jetIPSignedSigSum2D[jj]: 1);
    jetEIPSigned3D[jj] /= (has_track ?  jetIPSignedSigSum3D[jj]: 1);

    // compute combination variables
    jetMedianIPSignedSig2D[jj] = TrackAnalyzer::getJetMedian(liTrackIP2D, nLiTracks, jetJetID[jj], true);
    jetMedianIPSignedSig3D[jj] = TrackAnalyzer::getJetMedian(liTrackIP3D, nLiTracks, jetJetID[jj], true);
    jetVarianceIPSignedSig2D[jj] = TrackAnalyzer::getJetVariance(liTrackIP2D, jetMeanIPSignedSig2D[jj], nLiTracks, jetJetID[jj], true);
    jetVarianceIPSignedSig3D[jj] = TrackAnalyzer::getJetVariance(liTrackIP3D, jetMeanIPSignedSig3D[jj], nLiTracks, jetJetID[jj], true);    

    jetMedianIPSig2D[jj] = TrackAnalyzer::getJetMedian(liTrackIP2D, nLiTracks, jetJetID[jj], false);
    jetMedianIPSig3D[jj] = TrackAnalyzer::getJetMedian(liTrackIP3D, nLiTracks, jetJetID[jj], false);
    jetVarianceIPSig2D[jj] = TrackAnalyzer::getJetVariance(liTrackIP2D, jetMeanIPSig2D[jj], nLiTracks, jetJetID[jj], false);
    jetVarianceIPSig3D[jj] = TrackAnalyzer::getJetVariance(liTrackIP3D, jetMeanIPSig3D[jj], nLiTracks, jetJetID[jj], false);    

    // loop over all secondary vertices in the event
    for(int vv = 0; vv < nSV; vv++){    

      //skip low pt jets
      //if (svinfo2->jet()->pt() < cut_jetPt) continue;    
      // make sure the svinfo jet matches the jet we are looping on
      //if(debug > 1) std::cout << "[JETS] --------  svVertexJetID: " << svJetID[vv] << " jetjetid: " << jetJetID[jj] << std::endl;
      if (svVertexJetID[vv] != jetJetID[jj]) continue;
      
      // SV information
      jetSvLxy[jj] = 0;
      jetSvLxySig[jj] = 0;
      jetSvLxyz[jj] = 0;
      jetSvLxyzSig[jj] = 0;
      jetSvNTrack[jj] = 0; //vertex track multiplicty

      // SV position
      jetSvX[jj] = 0;
      jetSvY[jj] = 0;
      jetSvZ[jj] = 0;
      // error
      jetSvZErr[jj] = 0;
      jetSvYErr[jj] = 0;
      jetSvXErr[jj] = 0;

      //SV quality
      jetSvChi2[jj] = 0;
      jetSvNChi2[jj] = 0;
      jetSvNDof[jj] = 0;
      jetSvIsValid[jj] = 0;
           
      // find the vertex with the highest sum(p_t) in the full event list corresponding to the current jet
      jetNSv[jj] = 0;

      int bestSV = 0;
      int highestSum = -1;	
      for(int vvv = 0; vvv < nSV; vvv++) {
	// check the jet ids match between the SV and current jet	
	if (svVertexJetID[vvv] != jetJetID[jj]) continue;

	// counter the number of SV collected
	jetNSv[jj]++;
	//float sv_pt = svPt[vvv];
	float sv_ntracks = svNTracks[vvv];
	
	if (sv_ntracks > highestSum) {
	  bestSV = vv;
	  highestSum = sv_ntracks;
	  if(debug > 1) std::cout << "[JETS] -------- highest sum Pt PV: " << highestSum << " index " << vv << std::endl;
	}
      } // end vertex loop

      // no SV, nothing to do
      if (jetNSv[jj] == 0) continue;

      // SV information
      jetSvMass[jj] = svMass[bestSV];
      jetSvLxy[jj] = svFlight2D[bestSV];
      jetSvLxySig[jj] = svFlight2D[bestSV] / svFlight2DErr[bestSV];
      jetSvLxyz[jj] = svFlight[bestSV];
      jetSvLxyzSig[jj] =  svFlight[bestSV] / svFlightErr[bestSV];
      jetSvNTrack[jj] = svNTracks[bestSV]; //vertex track multiplicty

      // SV position and error
      jetSvX[jj] = svX[bestSV];
      jetSvY[jj] = svY[bestSV];
      jetSvZ[jj] = svZ[bestSV];

      jetSvZErr[jj] = svXErr[bestSV];
      jetSvYErr[jj] = svYErr[bestSV];
      jetSvXErr[jj] = svZErr[bestSV];

      //SV quality
      jetSvChi2[jj] = svChi2[bestSV];
      jetSvNChi2[jj] = svNChi2[bestSV];
      jetSvNDof[jj] = svNDof[bestSV];
      jetSvIsValid[jj] = svIsValid[bestSV];            
      
    } // end svtag (jet) loop

  } // end loop over jets

  evNum++;
  trackTree_->Fill();
  jetTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackAnalyzer::beginJob()
{

  if(debug > 1) std::cout << "[DEBUG] Setting Up Output File And Tree" << std::endl;
  
  // storage 
  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
  trackTree_  = new TTree("tracks","track, vertex,jet index tree");
  jetTree_    = new TTree("jets", "jet indexed tree");

  // book-keeping
  trackTree_->Branch("nCaloJets", &nCaloJets, "nCaloJets/I");
  trackTree_->Branch("nTracks", &nTracks, "nTracks/I");
  trackTree_->Branch("nTagJets", &nTagJets, "nTagJets/I");  
  trackTree_->Branch("nSvJets", &nSvJets, "nSvJets/I");
  trackTree_->Branch("nSV", &nSV, "nSV/I");
  trackTree_->Branch("nLiJets", &nLiJets, "nLiJets/I");

  trackTree_->Branch("evNum", &evNum, "evNum/I");

  ////////////////////////////// Calo Jet Information////////////////////////
  
  //jet kinematics
  trackTree_->Branch("caloJetPt", &caloJetPt, "caloJetPt[nCaloJets]/F");
  trackTree_->Branch("caloJetPhi", &caloJetPhi, "caloJetPhi[nCaloJets]/F");
  trackTree_->Branch("caloJetEta", &caloJetEta, "caloJetEta[nCaloJets]/F");

  // jet area
  trackTree_->Branch("caloJetn90", &caloJetN90, "caloJetN90[nCaloJets]/F");
  trackTree_->Branch("caloJetn60", &caloJetN60, "caloJetN60[nCaloJets]/F");
  trackTree_->Branch("caloJetTowerArea", &caloJetTowerArea, "caloJetTowerArea[nCaloJets]/F");

  // tag jet energy fraction
  trackTree_->Branch("caloJetHfrac", &caloJetHfrac, "caloJetHfrac[nCaloJets]/F");
  trackTree_->Branch("caloJetEfrac", &caloJetEfrac, "caloJetEFrac[nCaloJets]/F");

  ////////////////////////////// GEN MATCHING /////////////////////

  trackTree_->Branch("genMatch", &genMatch, "genMatch[nCaloJets]/I");  
  trackTree_->Branch("genMatchTrack", &genMatchTrack, "genMatchTrack[nTracks]/I");  
  trackTree_->Branch("genPt", &genPt, "genPt[nCaloJets]/I");  
  trackTree_->Branch("genEta", &genEta, "genEta[nCaloJets]/I");  
  trackTree_->Branch("genPhi", &genPhi, "genPhi[nCaloJets]/I");  
  trackTree_->Branch("genM", &genM, "genM[nCaloJets]/I");  

  ////////////////////////////// IP Jet tags////////////////////////

  trackTree_->Branch("jetJetID", &jetJetID, "jetJetID[nTagJets]/I");

  //tag jet kinematics
  trackTree_->Branch("tagJetPt", &tagJetPt, "tagJetPt[nTagJets]/F");
  trackTree_->Branch("tagJetPhi", &tagJetPhi, "tagJetPhi[nTagJets]/F");
  trackTree_->Branch("tagJetEta", &tagJetEta, "tagJetEta[nTagJets]/F");

  // tag jet area
  // trackTree_->Branch("tagJetn90", &tagJetN90, "tagJetN90[nTagJets]/F");
  // trackTree_->Branch("tagJetn60", &tagJetN60, "tagJetN60[nTagJets]/F");
  // trackTree_->Branch("tagJetTowerArea", &tagJetTowerArea, "tagJetTowerArea[nTagJets]/F");

  // tag jet energy fraction
  // trackTree_->Branch("tagJetHfrac", &tagJetHfrac, "tagJetHfrac[nTagJets]/F");
  // trackTree_->Branch("tagJetEfrac", &tagJetEfrac, "tagJetEFrac[nTagJets]/F");

  // tag jet track content
  trackTree_->Branch("tagJetNSelTracks", &tagJetNSelTracks, "tagJetNSelTracks[nTagJets]/F");  

  //track kinematics
  trackTree_->Branch("trackPt", &trackPt, "trackPt[nTracks]/F");
  trackTree_->Branch("trackPhi", &trackPhi, "trackPhi[nTracks]/F");
  trackTree_->Branch("trackEta", &trackEta, "trackEta[nTracks]/F");

  //track quality
  trackTree_->Branch("trackChi2", &trackChi2, "trackChi2[nTracks]/F");

  //ip tag info
  trackTree_->Branch("trackIP2D", &trackIP2D, "trackIP2D[nTracks]/F");
  trackTree_->Branch("trackIPSig2D", &trackIPSig2D, "trackIPSig2D[nTracks]/F");
  trackTree_->Branch("trackIP3D", &trackIP3D, "trackIP3D[nTracks]/F");
  trackTree_->Branch("trackIPSig3D", &trackIPSig3D, "trackIPSig3D[nTracks]/F");
  trackTree_->Branch("trackDistanceJetAxis", &trackDistanceJetAxis, "trackDistanceJetAxis[nTracks]/F");
  trackTree_->Branch("trackDistanceJetAxisSig", &trackDistanceJetAxisSig, "trackDistanceJetAxisSig[nTracks]/F");

  //track jet association
  trackTree_->Branch("trackJetID", &trackJetID, "trackJetID[nTracks]/I");  
  trackTree_->Branch("jetTrackDR", &jetTrackDR, "jetTrackDR[nTracks]/F");  

  ////////////////////////////// LIFETIME Jet tags////////////////////////

  // lifetime jets
  trackTree_->Branch("liJetID", &liJetID, "liJetID[nLiJets]/I");
  trackTree_->Branch("liJetPt", &liJetPt, "liJetPt[nLiJets]/F");
  trackTree_->Branch("liJetPhi", &liJetPhi, "liJetPhi[nLiJets]/F");
  trackTree_->Branch("liJetEta", &liJetEta, "liJetEta[nLiJets]/F");
  trackTree_->Branch("liJetNSelTracks", &liJetNSelTracks, "liJetNSelTracks[nLiJets]/I");  

  // lifetime track info
  trackTree_->Branch("nLiTracks", &nLiTracks, "nLiTracks/I");
  trackTree_->Branch("liTrackEta", &liTrackEta, "liTrackEta[nLiTracks]/F");
  trackTree_->Branch("liTrackPhi", &liTrackPhi, "liTrackPhi[nLiTracks]/F");
  trackTree_->Branch("liTrackPt", &liTrackPt, "liTrackPt[nLiTracks]/F");
  trackTree_->Branch("liTrackJetID", &liTrackJetID, "liTrackJetID[nLiTracks]/I");  
  trackTree_->Branch("liJetTrackDR", &liJetTrackDR, "liJetTrackDR[nLiTracks]/F");  

  // lifetime ip tag info
  trackTree_->Branch("liTrackIP2D", &liTrackIP2D, "liTrackIP2D[nLiTracks]/F");
  trackTree_->Branch("liTrackIPSig2D", &liTrackIPSig2D, "liTrackIPSig2D[nLiTracks]/F");
  trackTree_->Branch("liTrackIP3D", &liTrackIP3D, "liTrackIP3D[nLiTracks]/F");
  trackTree_->Branch("liTrackIPSig3D", &liTrackIPSig3D, "liTrackIPSig3D[nLiTracks]/F");
  trackTree_->Branch("liTrackDistanceJetAxis", &liTrackDistanceJetAxis, "liTrackDistanceJetAxis[nLiTracks]/F");
  trackTree_->Branch("liTrackDistanceJetAxisSig", &liTrackDistanceJetAxisSig, "liTrackDistanceJetAxisSig[nLiTracks]/F");


  //////////////////////////////SV Jet tags////////////////////////

  //sv jets
  trackTree_->Branch("svJetPt", &svJetPt, "svJetPt[nSvJets]/F");
  trackTree_->Branch("svJetEta", &svJetEta, "svJetEta[nSvJets]/F");
  trackTree_->Branch("svJetPhi", &svJetPhi, "svJetPhi[nSvJets]/F");
  trackTree_->Branch("svJetID", &svJetID, "svJetID[nSvJets]/I");

  // vertex
  trackTree_->Branch("svVertexJetID", &svVertexJetID, "svVertexJetID[nSV]/I");

  // quality
  trackTree_->Branch("svChi2", &svChi2, "svChi2[nSV]/F");
  trackTree_->Branch("svNChi2", &svNChi2, "svNChi2[nSV]/F");
  trackTree_->Branch("svNDof", &svNDof, "svNDof[nSV]/F");
  trackTree_->Branch("svIsValid", &svIsValid, "svIsValid[nSV]/I");

  // kinematics
  trackTree_->Branch("svMass", &svMass, "svMass[nSV]/F");
  trackTree_->Branch("svPt", &svPt, "svPt[nSV]/F");
  trackTree_->Branch("svPx", &svPx, "svPx[nSV]/F");
  trackTree_->Branch("svPy", &svPy, "svPy[nSV]/F");
  trackTree_->Branch("svEta", &svEta, "svEta[nSV]/F");
  trackTree_->Branch("svPhi", &svPhi, "svPhi[nSV]/F");

  // position and err
  trackTree_->Branch("svX", &svX, "svX[nSV]/F");
  trackTree_->Branch("svY", &svY, "svY[nSV]/F");
  trackTree_->Branch("svZ", &svZ, "svZ[nSV]/F");
  trackTree_->Branch("svXErr", &svXErr, "svXErr[nSV]/F");
  trackTree_->Branch("svYErr", &svYErr, "svYErr[nSV]/F");
  trackTree_->Branch("svZErr", &svZErr, "svZErr[nSV]/F");

  // vertex sum track charge
  trackTree_->Branch("svTotalCharge", &svTotalCharge, "svTotalCharge[nSV]/F");

  // flight
  trackTree_->Branch("svFlight", &svFlight, "svFlight[nSV]/F");
  trackTree_->Branch("svFlightErr", &svFlightErr, "svFlightErr[nSV]/F");
  trackTree_->Branch("svFlight2D", &svFlight2D, "svFlight2D[nSV]/F");
  trackTree_->Branch("svFlight2DErr", &svFlight2DErr, "svFlight2DErr[nSV]/F");

  // DR quantities
  trackTree_->Branch("svDRFlightJet", &svDRFlightJet, "svDRFlightJet[nSV]/F");
  trackTree_->Branch("svDRTrackJet", &svDRTrackJet, "svDRTrackJet[nSV]/F");
  trackTree_->Branch("svDRTrackFlight", &svDRTrackFlight, "svDRTrackFlight[nSV]/F");

  ///////////  ///////////  ///////////  ///////////  ///////////  ///////////  ////
  //////////////////////////////// JET TREE QUANITIES //////////////////////////////
  ///////////  ///////////  ///////////  ///////////  ///////////  /////////// /////
  //////// Everything is either a flat number or indexed by nCaloJets //////////////

  //book keeping
  jetTree_->Branch("nCaloJets", &nCaloJets, "nCaloJets/I");
  jetTree_->Branch("nTracks", &nTracks, "nTracks/I");
  jetTree_->Branch("evNum", &evNum, "evNum/I");
  jetTree_->Branch("jetID", &jetJetID, "jetID[nCaloJets]/I");

  jetTree_->Branch("nJetWithSv", &nJetWithSv, "nJetWithSv/I"); //number of SV in the event

  //////////////// GEN MATCHING ///////////////

  jetTree_->Branch("genMatch", &genMatch, "genMatch[nCaloJets]/I");  
  jetTree_->Branch("genPt", &genPt, "genPt[nCaloJets]/I");  
  jetTree_->Branch("genEta", &genEta, "genEta[nCaloJets]/I");  
  jetTree_->Branch("genPhi", &genPhi, "genPhi[nCaloJets]/I");  
  jetTree_->Branch("genM", &genM, "genM[nCaloJets]/I");  

  //////////////// CALO JETS ///////////////////
  
  //jet kinematics
  jetTree_->Branch("caloJetPt", &caloJetPt, "caloJetPt[nCaloJets]/F");
  jetTree_->Branch("caloJetPhi", &caloJetPhi, "caloJetPhi[nCaloJets]/F");
  jetTree_->Branch("caloJetEta", &caloJetEta, "caloJetEta[nCaloJets]/F");

  // jet area
  jetTree_->Branch("caloJetn90", &caloJetN90, "caloJetN90[nCaloJets]/F");
  jetTree_->Branch("caloJetn60", &caloJetN60, "caloJetN60[nCaloJets]/F");
  jetTree_->Branch("caloJetTowerArea", &caloJetTowerArea, "caloJetTowerArea[nCaloJets]/F");

  // tag jet energy fraction
  jetTree_->Branch("caloJetHfrac", &caloJetHfrac, "caloJetHfrac[nCaloJets]/F");
  jetTree_->Branch("caloJetEfrac", &caloJetEfrac, "caloJetEFrac[nCaloJets]/F");

  //////////////// IP TAG INFORMATION ////////////

  //number of selected tracks
  jetTree_->Branch("jetNTracks", &liJetNSelTracks, "jetNTracks[nCaloJets]/I");  

  // significance and absolute IP weighted track energy
  // unsigned 
  jetTree_->Branch("jetEIPSig2D", &jetEIPSig2D, "jetEIPSig2D[nCaloJets]/F");
  jetTree_->Branch("jetEIP2D", &jetEIP2D, "jetEIP2D[nCaloJets]/F");
  jetTree_->Branch("jetEIPSig3D", &jetEIPSig3D, "jetEIPSig3D[nCaloJets]/F");
  jetTree_->Branch("jetEIP3D", &jetEIP3D, "jetEIP3D[nCaloJets]/F");

  // signed
  jetTree_->Branch("jetEIPSignedSig2D", &jetEIPSignedSig2D, "jetEIPSignedSig2D[nCaloJets]/F");
  jetTree_->Branch("jetEIPSigned2D", &jetEIPSigned2D, "jetEIPSigned2D[nCaloJets]/F");
  jetTree_->Branch("jetEIPSignedSig3D", &jetEIPSignedSig3D, "jetEIPSignedSig3D[nCaloJets]/F");
  jetTree_->Branch("jetEIPSigned3D", &jetEIPSigned3D, "jetEIPSigned3D[nCaloJets]/F");

  // log weighted track pt 
  jetTree_->Branch("jetELogIPSig2D", &jetELogIPSig2D, "jetELogIPSig2D[nCaloJets]/F");
  //jetTree_->Branch("jetELogIP2D", &jetELogIP2D, "jetELogIP2D[nCaloJets]/F");
  jetTree_->Branch("jetELogIPSig3D", &jetELogIPSig3D, "jetELogIPSig3D[nCaloJets]/F");
  //jetTree_->Branch("jetELogIP3D", &jetELogIP3D, "jetELogIP3D[nCaloJets]/F");

  //ip significance sums -- sum(|IPsig|)
  jetTree_->Branch("jetIPSigSum2D", &jetIPSigSum2D, "jetIPSigSum2D[nCaloJets]/F");
  jetTree_->Branch("jetIPSigSum3D", &jetIPSigSum3D, "jetIPSigSum3D[nCaloJets]/F");

  //ip sig inverse sums -- sum(1 / |IPsig|)
  jetTree_->Branch("jetIPSigInvSum2D", &jetIPSigInvSum2D, "jetIPSigInvSum2D[nCaloJets]/F");
  jetTree_->Branch("jetIPSigInvSum3D", &jetIPSigInvSum3D, "jetIPSigInvSum3D[nCaloJets]/F");

  //signed ip significance sums  -- sum(IPsig)
  jetTree_->Branch("jetIPSignedSigSum2D", &jetIPSignedSigSum2D, "jetIPSignedSigSum2D[nCaloJets]/F");
  jetTree_->Branch("jetIPSignedSigSum3D", &jetIPSignedSigSum3D, "jetIPSignedSigSum3D[nCaloJets]/F");

  //signed ip sig inverse sums -- sum(1 / IPsig)
  jetTree_->Branch("jetIPSignedSigInvSum2D", &jetIPSignedSigInvSum2D, "jetIPSignedSigInvSum2D[nCaloJets]/F");
  jetTree_->Branch("jetIPSignedSigInvSum3D", &jetIPSignedSigInvSum3D, "jetIPSigInvSum3D[nCaloJets]/F");

  //ip sig log sums  -- sum(log(|IPsig|))
  jetTree_->Branch("jetIPSigLogSum2D", &jetIPSigLogSum2D, "jetIPSigLogSum2D[nCaloJets]/F");
  jetTree_->Branch("jetIPSigLogSum3D", &jetIPSigLogSum3D, "jetIPSigLogSum3D[nCaloJets]/F");

  // ip sig averages 
  jetTree_->Branch("jetMeanIPSig2D", &jetMeanIPSig2D, "jetMeanIPSig2D[nCaloJets]/F");
  jetTree_->Branch("jetMeanIPSig3D", &jetMeanIPSig3D, "jetMeanIPSig3D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPSig2D", &jetMedianIPSig2D, "jetMedianIPSig2D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPSig3D", &jetMedianIPSig3D, "jetMedianIPSig3D[nCaloJets]/F");
  jetTree_->Branch("jetVarianceIPSig2D", &jetVarianceIPSig2D, "jetVarianceIPSig2D[nCaloJets]/F");
  jetTree_->Branch("jetVarianceIPSig3D", &jetVarianceIPSig3D, "jetVarianceIPSig3D[nCaloJets]/F");

  // ip sig signed averages
  jetTree_->Branch("jetMeanIPSignedSig2D", &jetMeanIPSignedSig2D, "jetMeanIPSignedSig2D[nCaloJets]/F");
  jetTree_->Branch("jetMeanIPSignedSig3D", &jetMeanIPSignedSig3D, "jetMeanIPSignedSig3D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPSignedSig2D", &jetMedianIPSignedSig2D, "jetMedianIPSignedSig2D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPSignedSig3D", &jetMedianIPSignedSig3D, "jetMedianIPSignedSig3D[nCaloJets]/F");
  jetTree_->Branch("jetVarianceIPSignedSig2D", &jetVarianceIPSignedSig2D, "jetVarianceIPSignedSig2D[nCaloJets]/F");
  jetTree_->Branch("jetVarianceIPSignedSig3D", &jetVarianceIPSignedSig3D, "jetVarianceIPSignedSig3D[nCaloJets]/F");

  // ip value averages
  jetTree_->Branch("jetMeanIP2D", &jetMeanIP2D, "jetMeanIP2D[nCaloJets]/F");
  jetTree_->Branch("jetMeanIP3D", &jetMeanIP3D, "jetMeanIP3D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIP2D", &jetMedianIP2D, "jetMedianIP2D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIP3D", &jetMedianIP3D, "jetMedianIP3D[nCaloJets]/F");

  // signed ip value averages
  jetTree_->Branch("jetMeanIPSigned2D", &jetMeanIPSigned2D, "jetMeanIPSigned2D[nCaloJets]/F");
  jetTree_->Branch("jetMeanIPSigned3D", &jetMeanIPSigned3D, "jetMeanIPSigned3D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPSigned2D", &jetMedianIPSigned2D, "jetMedianIPSigned2D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPSigned3D", &jetMedianIPSigned3D, "jetMedianIPSigned3D[nCaloJets]/F");

  //////////////SECONDARY VTX INFORMATION //////////////

  // SV Information
  jetTree_->Branch("jetNSv", &jetNSv, "jetNSv[nCaloJets]/I");
  jetTree_->Branch("jetSvNTrack", &jetSvNTrack, "jetSvNTrack[nCaloJets]/I");
  jetTree_->Branch("jetSvMass", &jetSvMass, "jetSvMass[nCaloJets]/F");
  jetTree_->Branch("jetSvLxy", &jetSvLxy, "jetSvLxy[nCaloJets]/F");
  jetTree_->Branch("jetSvLxySig", &jetSvLxySig, "jetSvLxySig[nCaloJets]/F");
  jetTree_->Branch("jetSvLxyz", &jetSvLxyz, "jetSvLxyz[nCaloJets]/F");
  jetTree_->Branch("jetSvLxyzSig", &jetSvLxyzSig, "jetSvLxyzSig[nCaloJets]/F");

  // SV position
  jetTree_->Branch("jetSvX", &jetSvX, "jetSvX[nCaloJets]/F");
  jetTree_->Branch("jetSvY", &jetSvY, "jetSvY[nCaloJets]/F");
  jetTree_->Branch("jetSvZ", &jetSvZ, "jetSvZ[nCaloJets]/F");
  jetTree_->Branch("jetSvXErr", &jetSvXErr, "jetSvXErr[nCaloJets]/F");
  jetTree_->Branch("jetSvYErr", &jetSvYErr, "jetSvYErr[nCaloJets]/F");
  jetTree_->Branch("jetSvZErr", &jetSvZErr, "jetSvZErr[nCaloJets]/F");

  // SV Quality
  jetTree_->Branch("jetSvChi2", &jetSvChi2, "jetSvChi2[nCaloJets]/F");
  jetTree_->Branch("jetSvNChi2", &jetSvNChi2, "jetSvNChi2[nCaloJets]/F");
  jetTree_->Branch("jetSvNDof", &jetSvNDof, "jetSvNDof[nCaloJets]/F");  
  jetTree_->Branch("jetSvIsValid", &jetSvIsValid, "jetSvIsValid[nCaloJets]/I");  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackAnalyzer::endJob() 
{
  jetTree_->Write();
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

Float_t TrackAnalyzer::getJetMedian(Float_t values[], int size, int jetid, bool is_signed) {
  // Allocate an array of the same size and sort it.
  
  //the real size is only the tracks corresponding to the specific jet
  Int_t true_size = 0;
  for (int i = 0; i < size; ++i) {
    if (liTrackJetID[i] == jetid) {
      true_size++; 
    }
  }

  //loop over the full array but only fill with the corresponding jetid
  Float_t* sorted = new Float_t[true_size];
  int true_index = 0;
  for (int i = 0; true_index < true_size; ++i) {
    if (liTrackJetID[i] == jetid) {
      sorted[true_index] = is_signed ? values[i] : fabs(values[i]);
      true_index++; 
    }
  }

  // sort the array
  for (int i = true_size - 1; i > 0; --i) {
    for (int j = 0; j < i; ++j) {
      if (sorted[j] > sorted[j+1]) {
	Float_t dTemp = sorted[j];
	sorted[j] = sorted[j+1];
	sorted[j+1] = dTemp;
      }
    }
  }

  if(true_size == 0) {
    return 0; //no tracks
  }

  // Middle or average of middle values in the sorted array.
  Float_t dMedian = 0.0;
  if ((true_size % 2) == 0) {
    dMedian = (sorted[true_size/2] + sorted[(true_size/2) - 1])/2.0;
  } else {
    dMedian = sorted[true_size/2];
  }
  delete [] sorted;

  //safety check
  if(!is_signed ) { 
    if (dMedian < 0 ) std::cout << "ASSERT FAIL: "  << dMedian << std::endl;
    assert( dMedian >= 0 ); 
  }

  return dMedian;
}

Float_t TrackAnalyzer::getJetVariance(Float_t values[], Float_t mean, int n, int jetid, bool is_signed) {
  Float_t sum = 0;
  int true_n = 0;

  for (int i = 0; i < n; i++)  {
    if (liTrackJetID[n] == jetid) {
      Float_t val = is_signed ? values[i] : fabs(values[i]);
      sum += (val - mean) * (val - mean);
      true_n++;
    }
  }

  //safety check
  if(!is_signed ) { assert( sum >= 0 ); }

  return sum / float((true_n - 1));
}


//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
