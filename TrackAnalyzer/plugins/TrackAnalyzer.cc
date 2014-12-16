
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

  static const Int_t MAX_TRACKS = 9999;
  static const Int_t MAX_JETS = 999;
  static const Int_t MAX_VTX = 999;

  //bookkeeping
  int evNum = 0;
  
  //every jet has its own id
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

  //  auto_ptr<std::vector<reco::Vertex> > displacedSVertexCollection(new std::vector<reco::Vertex>(0,0));


  Int_t svVertexJetID[MAX_VTX]; 
  Int_t nSvJets = 0;
  Float_t svJetPt[MAX_JETS];
  Float_t svJetEta[MAX_JETS];
  Float_t svJetPhi[MAX_JETS];
  Float_t svJetID[MAX_JETS];

  Int_t nSvVertex = 0;

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
  Bool_t svIsValid[MAX_VTX];
  
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
  tag_secondaryVertexTagInfo_ = iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertexTagInfo");  
  tag_lifetimeIPTagInfo_ = iConfig.getUntrackedParameter<edm::InputTag>("lifetimeIPTagInfo");  

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

  edm::Handle<reco::TrackIPTagInfoCollection> lifetimeIPTagInfo;
  iEvent.getByLabel(tag_lifetimeIPTagInfo_, lifetimeIPTagInfo);

  edm::Handle<reco::SecondaryVertexTagInfoCollection> secondaryVertexTagInfo;
  iEvent.getByLabel(tag_secondaryVertexTagInfo_, secondaryVertexTagInfo);  

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
     
  //
  Int_t jj = 0;
  nCaloJets = 0;
  for(reco::CaloJetCollection::const_iterator jet = ak5CaloJets->begin(); jet != ak5CaloJets->end(); ++jet, jj++){

    //cuts 
    if (jet->pt() < cut_jetPt) continue;
    //if (fabs(jet->eta()) > cut_jetEta) continue;

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

#ifdef DEBUG
  std::cout << "[DEBUG] Begin Extracting Tag Info" << std::endl;
#endif

  //extract the collection from the handle
  const reco::TrackIPTagInfoCollection & ip = *(trackIPTagInfoCollection.product()); 
  const reco::TrackIPTagInfoCollection & lifetime = *(lifetimeIPTagInfo.product()); 
  const reco::SecondaryVertexTagInfoCollection & sv = *(secondaryVertexTagInfo.product()); 

  std::cout << "-- Found " << ip.size() << " IP TagInfo" << std::endl;
  std::cout << "-- Found " << lifetime.size() << " Lifetime TagInfo" << std::endl;
  std::cout << "-- Found " << sv.size() << " Secondary Vertex TagInfo" << std::endl;

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

#ifdef DEBUG
    std::cout << "JET ID: " << jj << std::endl;
    std::cout << "Jet pt: " << ipinfo->jet()->pt() << std::endl;
    std::cout << "Tot tracks: " << ipinfo->tracks().size() << std::endl;    
#endif
    
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
    
#ifdef DEBUG
      std::cout << "[DEBUG] N Selected Tracks: " << tagJetNSelTracks[nTagJets] << std::endl; 
      std::cout << "[DEBUG] Extracting Significances from Info" << std::endl; 
#endif 

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
    liJetNSelTracks[nLiJets] = liTracks.size();

    // redundancy for tree output
    liJetID[nLiJets] = jetid; 

    //kinematics
    liJetPt[nLiJets] = liinfo->jet()->pt();
    liJetEta[nLiJets] = liinfo->jet()->eta();
    liJetPhi[nLiJets] = liinfo->jet()->phi();        

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
  nSvVertex = 0;
  for(; svinfo != sv.end(); ++svinfo, jj++){    

    //skip low pt jets
    if (svinfo->jet()->pt() < cut_jetPt) continue;    

    //grab the tagging variables
    //    TaggingVariableList svVars = svinfo->taggingVariables();        

    reco::TrackRefVector svTracks = svinfo->selectedTracks();    
    int nSV = svinfo->nVertices();     

    // global jet identifier 
    svJetID[nSvJets] = jetid; 

    // jet kinematics
    svJetPt[nSvJets] = svinfo->jet()->pt();
    svJetEta[nSvJets] = svinfo->jet()->eta();
    svJetPhi[nSvJets] = svinfo->jet()->phi();        

    // number of secondary vertices 
    for(int vv = 0; vv < nSV; vv++) {
      svVertexJetID[vv] = jetid;

      reco::Vertex svVertex = svinfo->secondaryVertex(vv);  

      //place it in the collection to view later 
      //displacedSVertexCollection->push_back(svVertex);

      // vertex kinematics
      svMass[nSvVertex] = svinfo->secondaryVertex(vv).p4().mass();
      svPx[nSvVertex] = svVertex.p4().px();
      svPy[nSvVertex] = svVertex.p4().py();
      svPt[nSvVertex] = svVertex.p4().pt();
      svEta[nSvVertex] = svVertex.p4().eta();
      svPhi[nSvVertex] = svVertex.p4().phi();
      
      // quality
      svChi2[nSvVertex] = svVertex.chi2();
      svNChi2[nSvVertex] = svVertex.normalizedChi2();
      svNDof[nSvVertex] = svVertex.ndof();
      svIsValid[nSvVertex] = svVertex.isValid(); 	

      // positional space
      svX[nSvVertex] = svVertex.x();
      svXErr[nSvVertex] = svVertex.xError();
      svY[nSvVertex] = svVertex.y();   
      svYErr[nSvVertex] = svVertex.yError();   
      svZ[nSvVertex] = svVertex.z();      
      svZErr[nSvVertex] = svVertex.zError();      

      // std::vector<float> tagValList = svVars.getList(reco::btau::vertexJetDeltaR,false);    
      // std::copy( tagValList.begin(), tagValList.end(), &svFlight);
      // tagValList = svVars.getList(reco::btau::flightDistance2dVal,false);
      // std::copy( tagValList.begin(), tagValList.end(), &svFlight);
      // tagValList = svVars.getList(reco::btau::flightDistance2dSig,false);
      // std::copy( tagValList.begin(), tagValList.end(), );
      // tagValList = svVars.getList(reco::btau::flightDistance3dVal,false);
      // std::copy( tagValList.begin(), tagValList.end(), );
      // tagValList = svVars.getList(reco::btau::flightDistance3dSig,false);
      
      // flight
      svFlight[nSvVertex] = svinfo->flightDistance(vv).value();
      svFlightErr[nSvVertex] = svinfo->flightDistance(vv).error();
      svFlight2D[nSvVertex] = svinfo->flightDistance(vv, true).value();
      svFlight2DErr[nSvVertex] = svinfo->flightDistance(vv, true).error();

      // charge
      svTotalCharge[nSvVertex] = 0;

      //tracking 
      svNTracks[vv] = svVertex.nTracks();
      reco::TrackKinematics vertexKinematics;

      Bool_t hasRefittedTracks = svVertex.hasRefittedTracks();      

      // loop over the tracks associated with the Secondary Vertex to compute total charge
      for(int tt = 0; tt < svNTracks[nSvVertex]; tt++) {		
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
          svTotalCharge[nSvVertex] += refitTrack.charge();
        }
        else {
          vertexKinematics.add(track, trackWeight);
          svTotalCharge[nSvVertex] += track.charge();
        }
      } // end sv tracking look

      // vectors for deltaR calculation
      math::XYZTLorentzVector vertexSum = vertexKinematics.weightedVectorSum();
      math::XYZVector jetDir = svinfo->jet()->momentum().Unit();
      GlobalVector flightDir = svinfo->flightDirection(vv);

      // deltaR variables
      svDRFlightJet[nSvVertex] = reco::deltaR(flightDir, jetDir);
      svDRTrackJet[nSvVertex] = reco::deltaR(vertexSum, jetDir); 
      svDRTrackFlight[nSvVertex] = reco::deltaR(vertexSum, flightDir);       
      
      nSvVertex++; // index in array for vertex
    } // end sv vertex loop

    nSvJets++; //index in array for jet
    jetid++; //global jet identifier
  } // end sv tag (jet) loop

  //iEvent.put(displacedSVertexCollection, "displacedSVertexCollection");

  evNum++;
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
  trackTree_  = new TTree("tree","tag tree");
  
  // book-keeping
  trackTree_->Branch("nCaloJets", &nCaloJets, "nCaloJets/I");
  trackTree_->Branch("nTracks", &nTracks, "nTracks/I");

  trackTree_->Branch("jetJetID", &jetJetID, "jetJetID[nTagJets]/I");
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

  ////////////////////////////// IP Jet tags////////////////////////

  //jet tags
  trackTree_->Branch("nTagJets", &nTagJets, "nTagJets/I");

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

  //lifetime jets
  trackTree_->Branch("nLiJets", &nLiJets, "nLiJets/I");
  trackTree_->Branch("liJetID", &liJetID, "liJetID[nLiJets]/I");
  trackTree_->Branch("liJetPt", &liJetPt, "liJetPt[nLiJets]/F");
  trackTree_->Branch("liJetPhi", &liJetPhi, "liJetPhi[nLiJets]/F");
  trackTree_->Branch("liJetEta", &liJetEta, "liJetEta[nLiJets]/F");
  trackTree_->Branch("liJetNSelTracks", &liJetNSelTracks, "liJetNSelTracks[nLiJets]/F");  

  //lifetime ip tag info
  trackTree_->Branch("liTrackIP2D", &liTrackIP2D, "liTrackIP2D[nTracks]/F");
  trackTree_->Branch("liTrackIPSig2D", &liTrackIPSig2D, "liTrackIPSig2D[nTracks]/F");
  trackTree_->Branch("liTrackIP3D", &liTrackIP3D, "liTrackIP3D[nTracks]/F");
  trackTree_->Branch("liTrackIPSig3D", &liTrackIPSig3D, "liTrackIPSig3D[nTracks]/F");
  trackTree_->Branch("liTrackDistanceJetAxis", &liTrackDistanceJetAxis, "liTrackDistanceJetAxis[nTracks]/F");
  trackTree_->Branch("liTrackDistanceJetAxisSig", &liTrackDistanceJetAxisSig, "liTrackDistanceJetAxisSig[nTracks]/F");

  trackTree_->Branch("liTrackJetID", &liTrackJetID, "liTrackJetID[nTracks]/I");  
  trackTree_->Branch("liJetTrackDR", &liJetTrackDR, "liJetTrackDR[nTracks]/F");  


  //////////////////////////////SV Jet tags////////////////////////

  //sv jets
  trackTree_->Branch("nSvJets", &nSvJets, "nSVJets/I");
  trackTree_->Branch("svJetPt", &svJetPt, "svJetPt[nSvJets]/F");
  trackTree_->Branch("svJetEta", &svJetEta, "svJetEta[nSvJets]/F");
  trackTree_->Branch("svJetPhi", &svJetPhi, "svJetPhi[nSvJets]/F");
  trackTree_->Branch("svJetID", &svJetID, "svJetID[nSvJets]/F");

  // vertex
  trackTree_->Branch("nSvVertex", &nSvVertex, "svVertexJetID/I");
  trackTree_->Branch("svVertexJetID", &svVertexJetID, "svVertexJetID[nSvVertex]/F");

  // quality
  trackTree_->Branch("svChi2", &svChi2, "svChi2[nSvVertex]/F");
  trackTree_->Branch("svNChi2", &svNChi2, "svNChi2[nSvVertex]/F");
  trackTree_->Branch("svNDof", &svNDof, "svNDof[nSvVertex]/F");
  trackTree_->Branch("svIsValid", &svIsValid, "svIsValid[nSvVertex]/F");

  // kinematics
  trackTree_->Branch("svMass", &svMass, "svMass[nSvVertex]/F");
  trackTree_->Branch("svPt", &svPt, "svPt[nSvVertex]/F");
  trackTree_->Branch("svPx", &svPx, "svPx[nSvVertex]/F");
  trackTree_->Branch("svPy", &svPy, "svPy[nSvVertex]/F");
  trackTree_->Branch("svEta", &svEta, "svEta[nSvVertex]/F");
  trackTree_->Branch("svPhi", &svPhi, "svPhi[nSvVertex]/F");

  // position and err
  trackTree_->Branch("svX", &svX, "svX[nSvVertex]/F");
  trackTree_->Branch("svY", &svY, "svY[nSvVertex]/F");
  trackTree_->Branch("svZ", &svZ, "svZ[nSvVertex]/F");
  trackTree_->Branch("svXErr", &svXErr, "svXErr[nSvVertex]/F");
  trackTree_->Branch("svYErr", &svYErr, "svYErr[nSvVertex]/F");
  trackTree_->Branch("svZErr", &svZErr, "svZErr[nSvVertex]/F");

  // vertex sum track charge
  trackTree_->Branch("svTotalCharge", &svTotalCharge, "svTotalCharge[nSvVertex]/F");

  // flight
  trackTree_->Branch("svFlight", &svFlight, "svFlight[nSvVertex]/F");
  trackTree_->Branch("svFlightErr", &svFlightErr, "svFlightErr[nSvVertex]/F");
  trackTree_->Branch("svFlight2D", &svFlight2D, "svFlight2D[nSvVertex]/F");
  trackTree_->Branch("svFlight2DErr", &svFlight2DErr, "svFlight2DErr[nSvVertex]/F");

  // DR quantities
  trackTree_->Branch("svDRFlightJet", &svDRFlightJet, "svDRFlightJet[nSvVertex]/F");
  trackTree_->Branch("svDRTrackJet", &svDRTrackJet, "svDRTrackJet[nSvVertex]/F");
  trackTree_->Branch("svDRTrackFlight", &svDRTrackFlight, "svDRTrackFlight[nSvVertex]/F");

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
