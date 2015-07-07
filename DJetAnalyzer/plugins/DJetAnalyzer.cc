// -*- C++ -*-
//
// Package:    DisplacedJets/DJetAnalyzer
// Class:      DJetAnalyzer
// 
/**\class DJetAnalyzer DJetAnalyzer.cc DJetAnalyzer/DJetAnalyzer/plugins/DJetAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Robert Hardenbrook
//         Created:  Sat, 15 Nov 2014 15:18:18 GMT
//
//

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

// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// trigger
#include "DataFormats/Common/interface/TriggerResults.h"

// C++ EDM Replacements
#include "DataFormats/Common/interface/Ref.h"

// geometry
#include "DataFormats/GeometrySurface/interface/Line.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

// gen information
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// tracking
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

//vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BTauReco/interface/VertexTypes.h"

// sim information
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

// jets
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"

// btau
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPData.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"

// messages
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// kinematics
#include <RecoBTag/SecondaryVertex/interface/TrackKinematics.h>
#include "DataFormats/Math/interface/deltaR.h"

// user defined includes
#include "DisplacedJets/DisplacedJetSVAssociator/interface/JetVertexAssociation.h"
#include "DisplacedJets/DisplacedJet/interface/DisplacedJet.h"
#include "DisplacedJets/DisplacedJet/interface/DisplacedJetEvent.h"
#include "DisplacedJets/DJetAnalyzer/interface/DJetAnalyzer.h"
//
// class declaration
//
DJetAnalyzer::DJetAnalyzer(const edm::ParameterSet& iConfig)
{
  // output configuration
  debug		  = iConfig.getUntrackedParameter<int>("debugLevel");  
  outputFileName_ = iConfig.getUntrackedParameter<std::string>("outputFileName");
  jetTreeName_	  = iConfig.getUntrackedParameter<std::string>("jetTreeName");
  trackTreeName_  = iConfig.getUntrackedParameter<std::string>("trackTreeName");
  vertexTreeName_ = iConfig.getUntrackedParameter<std::string>("vertexTreeName");
  genTreeName_	  = iConfig.getUntrackedParameter<std::string>("genTreeName");

  // sample information
  isMC_	      = iConfig.getUntrackedParameter<bool>("isMC");
  isSignalMC_ = iConfig.getUntrackedParameter<bool>("isSignalMC");
  
  // tag classification
  shortTagThresDist  = iConfig.getUntrackedParameter<double>("shortTagThreshold");
  mediumTagThresDist = iConfig.getUntrackedParameter<double>("mediumTagThreshold");
  longTagThresDist   = iConfig.getUntrackedParameter<double>("longTagThreshold");
  
  //DHT working point
  dHTWorkingPoint  = iConfig.getUntrackedParameter<int>("dHTWorkingPoint");

  // analysis todos
  doGenMatch_		  = iConfig.getUntrackedParameter<bool>("doGenMatch");
  doSimMatch_		  = iConfig.getUntrackedParameter<bool>("doSimMatch");
  applyEventPreSelection_ = iConfig.getUntrackedParameter<bool>("applyEventPreSelection");
  applyJetPreSelection_	  = iConfig.getUntrackedParameter<bool>("applyJetPreSelection");

  // collection tags
  tag_generalTracks_		  = iConfig.getUntrackedParameter<edm::InputTag>("generalTracks");
  tag_ak4CaloJets_		  = iConfig.getUntrackedParameter<edm::InputTag>("ak4CaloJets");
  tag_secondaryVertexTagInfo_	  = iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertexTagInfo");  
  tag_lifetimeIPTagInfo_	  = iConfig.getUntrackedParameter<edm::InputTag>("lifetimeIPTagInfo"); 

  // vertex tags
  tag_secondaryVertices_	  = iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertex"); 
  tag_inclusiveVertexCandidates_  = iConfig.getUntrackedParameter<edm::InputTag>("inclusiveVertexCand"); 
  tag_inclusiveSecondaryVertices_ = iConfig.getUntrackedParameter<edm::InputTag>("inclusiveVertexSecondary"); 
  tag_offlinePrimaryVertices_	  = iConfig.getUntrackedParameter<edm::InputTag>("offlinePrimaryVertices"); 
  
  //cuts 
  cut_jetPt  = iConfig.getUntrackedParameter<double>("jetPt");
  cut_jetEta = iConfig.getUntrackedParameter<double>("jetEta");

  //mc tags
  if(isMC_) {
    //tag_ak5GenJets_ = iConfig.getUntrackedParameter<edm::InputTag>("ak5GenJets");
    //tag_genMetCalo_ = iConfig.getUntrackedParameter<edm::InputTag>("genMetCalo");
    tag_genParticles_ = iConfig.getUntrackedParameter<edm::InputTag>("genParticles");
    tag_simVertex_    = iConfig.getUntrackedParameter<edm::InputTag>("simVertices");
  } 
}

DJetAnalyzer::~DJetAnalyzer(){ }

void DJetAnalyzer::fillHandles(const edm::Event & iEvent ) {

  // AOD Compatible
  iEvent.getByLabel(tag_generalTracks_, gTracks); 
  iEvent.getByLabel(tag_ak4CaloJets_, ak4CaloJets);

  // tag info
  iEvent.getByLabel(tag_lifetimeIPTagInfo_, lifetimeIPTagInfo);
  iEvent.getByLabel(tag_secondaryVertexTagInfo_, secondaryVertexTagInfo);  

  // vertex info
  iEvent.getByLabel(tag_secondaryVertices_, secondaryVertices);  
  iEvent.getByLabel(tag_inclusiveVertexCandidates_, inclusiveVertexCandidates);  
  iEvent.getByLabel(tag_inclusiveSecondaryVertices_, inclusiveSecondaryVertices);  
  iEvent.getByLabel(tag_offlinePrimaryVertices_, offlinePrimaryVertices);  

  // and sim matching quantities related to MC
  if(isMC_) {    
    if(doGenMatch_) iEvent.getByLabel(tag_genParticles_, genParticles);
    if(doSimMatch_) iEvent.getByLabel(tag_simVertex_, simVertices);
  }  
}

void 
DJetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(debug > 0) std::cout << "[----------------- ANALYZE EVENT DEBUG LEVEL: " << debug  << " --------------------]" << std::endl;
   
  /////////////////////////////////
  // Event Setup
  /////////////////////////////////

  // fill the handles before grabbing the products from them
  fillHandles(iEvent);

  // collection products
  const reco::CaloJetCollection &		    caloJets	     = *(ak4CaloJets.product());
  const reco::VertexCollection &                    pvCollection     = *(offlinePrimaryVertices.product());
  const reco::TrackIPTagInfoCollection &	    lifetimeTagInfo  = *(lifetimeIPTagInfo.product()); ;
  const reco::SecondaryVertexTagInfoCollection &    svTagInfo	     = *(secondaryVertexTagInfo.product()); 
  const reco::VertexCollection &		    inc		     = *(inclusiveVertexCandidates.product());
  const reco::VertexCollection &		    incSV	     = *(inclusiveSecondaryVertices.product());
  const reco::GenParticleCollection &		    genCollection    = *(genParticles.product());     
  const edm::SimVertexContainer &		    simVtxCollection = *(simVertices.product()); 
  const reco::TrackCollection &                     generalTracks    = *(gTracks.product());

  // build the displaced event from the calo jet collection and kinematic cuts
  DisplacedJetEvent djEvent(isMC_, caloJets, pvCollection, cut_jetPt, cut_jetEta, debug);

  // pull out the first primary vertex in the collection (the default PV)
  const reco::Vertex & firstPV = *pvCollection.begin(); 
   
  /////////////////////////////////
  // Fill Trees
  /////////////////////////////////

  // event information 
  run	= iEvent.id().run();
  lumi	= iEvent.id().luminosityBlock();
  event = iEvent.id().event();      

  // ncalo jets indexes every jet branch
  nCaloJets = djEvent.getNJets();   

  // merge in the event info
  djEvent.mergeCaloIPTagInfo(lifetimeTagInfo); // add the ip info built from the JTA

  // dump information related to the preselection
  dumpPreSelection(djEvent);

  evNum++;   
  if(debug > 1) std::cout << "[DEBUG] Fill Run Stat Tree" << std::endl;
  runStatTree_->Fill(); // always fill the run stats  

  // dont keep events not passing preselection 
  if (applyEventPreSelection_ && eventPassEventPreSelection == 0) {
    return;
  }

  djEvent.mergeSVTagInfo(svTagInfo); // add the secondary vertexer from btag
  djEvent.addIVFVertices(incSV); // add the inclusive secondary vertices (includes matching and calculations)

  // mc matching (gen vertex and gen particle to calo jet matching) 
  // (genparticles, particle matching, vtx matching, vtx id matching, threshold for vtx match)
  if(isMC_ && doGenMatch_) djEvent.doGenMatching(genCollection, true, true, true, 0.3, 0.7, 0.05);    
  // only dump sim information for matching
  if(isMC_ && doSimMatch_) dumpSimInfo(simVtxCollection);  

  // only do tagging after all information has been added / merged
  // currently using the same thresholds for each category
  const std::vector<float> thresholds{0.0, 1.0, 2.0, 3.0, 4.0};

  // tagging thresholds (nvtx, short, medium, long)
  // ThresDist = Threshold Distance in cm for the vertex used to categorize the jets
  djEvent.doJetTagging(thresholds, thresholds, thresholds, thresholds,
		       shortTagThresDist, mediumTagThresDist, longTagThresDist, dHTWorkingPoint);
    
  // dump the displaced jet info into the corresponding branches by event
  dumpCaloInfo(djEvent);
  dumpIPInfo(djEvent);
  dumpIVFInfo(djEvent);  
  dumpSVTagInfo(djEvent);
  dumpDJTags(djEvent);

  // dump the track information
  nTracks = 0; 
  dumpTrackInfo(djEvent, generalTracks, 0);

  //dumpDTrackInfo(djEvent);

  // dump the vertex info in the event TODO
  dumpPVInfo(djEvent, pvCollection);

  if(debug > 1) std::cout << "[DEBUG] Fill Event Tree" << std::endl;
  eventTree_->Fill();
  if(debug > 1) std::cout << "[DEBUG] Fill Track Tree" << std::endl;
  trackTree_->Fill();
  if(debug > 1) std::cout << "[DEBUG] Fill Jet Tree" << std::endl;
  jetTree_->Fill();
  if(debug > 1) std::cout << "[DEBUG] Fill VTX Tree" << std::endl; 
  vertexTree_->Fill();
  //   if(debug > 1) std::cout << "[DEBUG] Fill GEN Tree" << std::endl; TODO
  //  genTree_->Fill();
}

void 
DJetAnalyzer::beginJob()
{

  if(debug > 1) std::cout << "[DEBUG 1] Setting Up Output File And Tree" << std::endl;
  
  // storage 
  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
  trackTree_  = new TTree(trackTreeName_.c_str(), "track, vertex, jet index tree");
  jetTree_    = new TTree(jetTreeName_.c_str(), "jet indexed tree");
  vertexTree_ = new TTree(vertexTreeName_.c_str(), "vertex indexed tree");
  genTree_    = new TTree(genTreeName_.c_str(), "Gen Particle Info tree");
  eventTree_  = new TTree("eventInfo", "Event Information Tree");
  runStatTree_  = new TTree("runStats", "Run Statistics");

  ///////////  ///////////  ///////////  ///////////  ///////////  ///////////  ////
  //////////////////////////////// EVENT TREE QUANITIES ////////////////////////////
  ///////////  ///////////  ///////////  ///////////  ///////////  /////////// /////

  // global book keeping
  runStatTree_->Branch("run", &run, "run/I");
  runStatTree_->Branch("lumi", &lumi, "lumi/I");
  runStatTree_->Branch("event", &event, "event/I");

  // local event book keeping 
  runStatTree_->Branch("evNum", &evNum, "evNum/I");
  // analysis information
  runStatTree_->Branch("eventPassEventPreSelection", &eventPassEventPreSelection, "eventPassEventPreSelection/I");

  ///////////  ///////////  ///////////  ///////////  ///////////  ///////////  ////
  //////////////////////////////// EVENT TREE QUANITIES ////////////////////////////
  ///////////  ///////////  ///////////  ///////////  ///////////  /////////// /////

  // global book keeping
  eventTree_->Branch("run", &run, "run/I");
  eventTree_->Branch("lumi", &lumi, "lumi/I");
  eventTree_->Branch("event", &event, "event/I");

  // local event book keeping 
  eventTree_->Branch("evNum", &evNum, "evNum/I");
  eventTree_->Branch("nCaloJets", &nCaloJets, "nCaloJets/I");
  eventTree_->Branch("eventNWP", &nWP, "eventNWP/I");

  // analysis information
  eventTree_->Branch("eventPassEventPreSel", &eventPassEventPreSelection, "eventPassEventPreSel/I");
  eventTree_->Branch("eventCaloHT", &eventCaloHT, "eventCaloHT/F");
  eventTree_->Branch("eventCaloDHT", &eventCaloDHT, "eventCaloDHT[eventNWP]/F");
  eventTree_->Branch("eventCaloMET", &eventCaloMET, "eventCaloMET/F");

  // ivf related
  eventTree_->Branch("eventNIVFReco", &eventNIVFReco, "eventNIVFReco/I");
  eventTree_->Branch("eventNIVFRecoGenMatch", &eventNIVFRecoGenMatch, "eventNIVFRecoGenMatch/I");

  // event tags
  eventTree_->Branch("eventNNoVertexTags", &eventNNoVertexTags, "eventNNoVertexTags[eventNWP]/I");
  eventTree_->Branch("eventNShortTags", &eventNShortTags, "eventNShortTags[eventNWP]/I");
  eventTree_->Branch("eventNMediumTags", &eventNMediumTags, "eventNMediumTags[eventNWP]/I");
  eventTree_->Branch("eventNLongTags", &eventNLongTags, "eventNLongTags[eventNWP]/I");
  eventTree_->Branch("eventNTotalTags", &eventNTotalTags, "eventNTotalTags[eventNWP]/I");

  ///////////  ///////////  ///////////  ///////////  ///////////  ///////////  ////
  //////////////////////////////// JET TREE QUANITIES //////////////////////////////
  ///////////  ///////////  ///////////  ///////////  ///////////  /////////// /////
  //////// Everything is either a flat number or indexed by nCaloJets //////////////
  //////// Exceptions are made for IVF gen matching quantities for moms and sons ///

  // global book keeping
  jetTree_->Branch("run", &run, "run/I");
  jetTree_->Branch("lumi", &lumi, "lumi/I");
  jetTree_->Branch("event", &event, "event/I");

  // branch indices must be defined first
  jetTree_->Branch("nCaloJets", &nCaloJets, "nCaloJets/I");
  jetTree_->Branch("nJetWithSv", &nJetWithSv, "nJetWithSv/I"); //number of SV in the event

  // local book keeping
  jetTree_->Branch("evNum", &evNum, "evNum/I");
  jetTree_->Branch("jetID", &caloJetID, "jetID[nCaloJets]/I");

  // analysis book keeping
  jetTree_->Branch("eventPassEventPreSelection", &eventPassEventPreSelection, "eventPassEventPreSelection/I");
  jetTree_->Branch("jetPassPreSelection", &jetPassPreSelection, "jetPassPreSelection[nCaloJets]/I");

  //////////////// CALO JETS ///////////////////
  
  //jet kinematics
  jetTree_->Branch("caloJetPt", &caloJetPt, "caloJetPt[nCaloJets]/F");
  jetTree_->Branch("caloJetPhi", &caloJetPhi, "caloJetPhi[nCaloJets]/F");
  jetTree_->Branch("caloJetEta", &caloJetEta, "caloJetEta[nCaloJets]/F");

  // jet size 
  jetTree_->Branch("caloJetn90", &caloJetN90, "caloJetN90[nCaloJets]/F");
  jetTree_->Branch("caloJetn60", &caloJetN60, "caloJetN60[nCaloJets]/F");
  jetTree_->Branch("caloJetTowerArea", &caloJetTowerArea, "caloJetTowerArea[nCaloJets]/F");

  // tag jet energy fraction
  jetTree_->Branch("caloJetHfrac", &caloJetHfrac, "caloJetHfrac[nCaloJets]/F");
  jetTree_->Branch("caloJetEfrac", &caloJetEfrac, "caloJetEFrac[nCaloJets]/F");

  jetTree_->Branch("caloGenMatch", &caloGenMatch, "caloGenMatch[nCaloJets]/I");  
  jetTree_->Branch("caloGenPt", &caloGenPt, "caloGenPt[nCaloJets]/F");  
  jetTree_->Branch("caloGenEta", &caloGenEta, "caloGenEta[nCaloJets]/F");  
  jetTree_->Branch("caloGenPhi", &caloGenPhi, "caloGenPhi[nCaloJets]/F");  
  jetTree_->Branch("caloGenM", &caloGenM, "caloGenM[nCaloJets]/F");  

  //////////////// IP TAG INFORMATION ////////////

  //number of selected tracks
  jetTree_->Branch("jetNTracks", &jetNTracks, "jetNTracks[nCaloJets]/I");  

  // significance and absolute IP weighted track energy
  // unsigned 
  jetTree_->Branch("jetEIPSig2D", &jetEIPSig2D, "jetEIPSig2D[nCaloJets]/F");
  jetTree_->Branch("jetEIPSigLog2D", &jetEIPSigLog2D, "jetEIPSigLog2D[nCaloJets]/F");
  jetTree_->Branch("jetEIPSig3D", &jetEIPSig3D, "jetEIPSig3D[nCaloJets]/F");
  jetTree_->Branch("jetEIPSigLog3D", &jetEIPSigLog3D, "jetEIPSigLog3D[nCaloJets]/F");

  // log weighted track pt 
  jetTree_->Branch("jetELogIPSig2D", &jetELogIPSig2D, "jetELogIPSig2D[nCaloJets]/F");
  jetTree_->Branch("jetELogIPSig3D", &jetELogIPSig3D, "jetELogIPSig3D[nCaloJets]/F");

  //ip significance sums -- sum(|IPsig|)
  // jetTree_->Branch("jetIPSigSum2D", &jetIPSigSum2D, "jetIPSigSum2D[nCaloJets]/F");
  // jetTree_->Branch("jetIPSigSum3D", &jetIPSigSum3D, "jetIPSigSum3D[nCaloJets]/F");

  // ip sig log sums  -- sum(log(|IPsig|))
  jetTree_->Branch("jetIPSigLogSum2D", &jetIPSigLogSum2D, "jetIPSigLogSum2D[nCaloJets]/F");
  jetTree_->Branch("jetIPSigLogSum3D", &jetIPSigLogSum3D, "jetIPSigLogSum3D[nCaloJets]/F");

  jetTree_->Branch("jetDistSigLogSum", &jetDistSigLogSum, "jetDistSigLogSum[nCaloJets]/F");
  jetTree_->Branch("jetDistLogSum", &jetDistLogSum, "jetDistLogSum[nCaloJets]/F");

  //jetTree_->Branch("jetIPLogSum2D", &jetIPLogSum2D, "jetIPLogSum2D[nCaloJets]/F");
  //jetTree_->Branch("jetIPLogSum3D", &jetIPLogSum3D, "jetIPLogSum3D[nCaloJets]/F");

  // ip sig averages 
  jetTree_->Branch("jetMeanIPSig2D", &jetMeanIPSig2D, "jetMeanIPSig2D[nCaloJets]/F");
  jetTree_->Branch("jetMeanIPSig3D", &jetMeanIPSig3D, "jetMeanIPSig3D[nCaloJets]/F");
  jetTree_->Branch("jetMeanIPLogSig2D", &jetMeanIPLogSig2D, "jetMeanIPLogSig2D[nCaloJets]/F");
  jetTree_->Branch("jetMeanIPLogSig3D", &jetMeanIPLogSig3D, "jetMeanIPLogSig3D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPSig2D", &jetMedianIPSig2D, "jetMedianIPSig2D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPSig3D", &jetMedianIPSig3D, "jetMedianIPSig3D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPLogSig2D", &jetMedianIPLogSig2D, "jetMedianIPLogSig2D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPLogSig3D", &jetMedianIPLogSig3D, "jetMedianIPLogSig3D[nCaloJets]/F");
  jetTree_->Branch("jetVarianceIPSig2D", &jetVarianceIPSig2D, "jetVarianceIPSig2D[nCaloJets]/F");
  jetTree_->Branch("jetVarianceIPSig3D", &jetVarianceIPSig3D, "jetVarianceIPSig3D[nCaloJets]/F");
  jetTree_->Branch("jetVarianceIPLogSig2D", &jetVarianceIPLogSig2D, "jetVarianceIPLogSig2D[nCaloJets]/F");
  jetTree_->Branch("jetVarianceIPLogSig3D", &jetVarianceIPLogSig3D, "jetVarianceIPLogSig3D[nCaloJets]/F");
  jetTree_->Branch("jetVarianceJetDist", &jetVarianceJetDist, "jetVarianceJetDist[nCaloJets]/F");
  jetTree_->Branch("jetVarianceJetDistSig", &jetVarianceJetDistSig, "jetVarianceJetDistSig[nCaloJets]/F");

  // ip value averages
  jetTree_->Branch("jetMeanIP2D", &jetMeanIP2D, "jetMeanIP2D[nCaloJets]/F");
  jetTree_->Branch("jetMeanIP3D", &jetMeanIP3D, "jetMeanIP3D[nCaloJets]/F");
  jetTree_->Branch("jetMeanIPLog2D", &jetMeanIPLog2D, "jetMeanIPLog2D[nCaloJets]/F");
  jetTree_->Branch("jetMeanIPLog3D", &jetMeanIPLog3D, "jetMeanIPLog3D[nCaloJets]/F");
  jetTree_->Branch("jetMeanJetDist", &jetMeanJetDist, "jetMeanJetDist[nCaloJets]/F");
  jetTree_->Branch("jetMedianIP2D", &jetMedianIP2D, "jetMedianIP2D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIP3D", &jetMedianIP3D, "jetMedianIP3D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPLog2D", &jetMedianIPLog2D, "jetMedianIPLog2D[nCaloJets]/F");
  jetTree_->Branch("jetMedianIPLog3D", &jetMedianIPLog3D, "jetMedianIPLog3D[nCaloJets]/F");
  jetTree_->Branch("jetMedianJetDist", &jetMedianJetDist, "jetMedianJetDist[nCaloJets]/F");

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
  jetTree_->Branch("jetSvNDof", &jetSvNDof, "jetSvNDof[nCaloJets]/F");  
  // jetTree_->Branch("jetSvNChi2", &jetSvNChi2, "jetSvNChi2[nCaloJets]/F");
  // jetTree_->Branch("jetSvIsValid", &jetSvIsValid, "jetSvIsValid[nCaloJets]/I");  

  //////////////INCLUSIVE SECONDARY VTX INFORMATION //////////////

  // IVF Information
  jetTree_->Branch("jetIVFNTrack", &jetIVFNTrack, "jetIVFNTrack[nCaloJets]/I");
  jetTree_->Branch("jetIVFMass", &jetIVFMass, "jetIVFMass[nCaloJets]/F");
  jetTree_->Branch("jetIVFLxy", &jetIVFLxy, "jetIVFLxy[nCaloJets]/F");
  jetTree_->Branch("jetIVFLxySig", &jetIVFLxySig, "jetIVFLxySig[nCaloJets]/F");
  jetTree_->Branch("jetIVFLxyz", &jetIVFLxyz, "jetIVFLxyz[nCaloJets]/F");
  jetTree_->Branch("jetIVFLxyzSig", &jetIVFLxyzSig, "jetIVFLxyzSig[nCaloJets]/F");

  // IVF Position
  jetTree_->Branch("jetIVFX", &jetIVFX, "jetIVFX[nCaloJets]/F");
  jetTree_->Branch("jetIVFY", &jetIVFY, "jetIVFY[nCaloJets]/F");
  jetTree_->Branch("jetIVFZ", &jetIVFZ, "jetIVFZ[nCaloJets]/F");
  jetTree_->Branch("jetIVFXErr", &jetIVFXErr, "jetIVFXErr[nCaloJets]/F");
  jetTree_->Branch("jetIVFYErr", &jetIVFYErr, "jetIVFYErr[nCaloJets]/F");
  jetTree_->Branch("jetIVFZErr", &jetIVFZErr, "jetIVFZErr[nCaloJets]/F");

  // IVF Vertex Picking Metric
  jetTree_->Branch("jetIVFMatchingScore", &jetIVFMatchingScore, "jetIVFMatchingScore[nCaloJets]/F");  

  // IVF Gen Matching
  jetTree_->Branch("jetIVFGenVertexMatched", &jetIVFGenVertexMatched, "jetIVFGenVertexMatched[nCaloJets]/I");  
  jetTree_->Branch("jetIVFGenVertexMatchMetric", &jetIVFGenVertexMatchMetric, "jetIVFGenVertexMatchMetric[nCaloJets]/F");  
  // jetTree_->Branch("jetIVFSimVertexMatched", &jetIVFSimVertexMatched, "jetIVFSimVertexMatched[nCaloJets]/I");  
  // jetTree_->Branch("jetIVFSimVertexMatchMetric", &jetIVFSimVertexMatchMetric, "jetIVFSimVertexMatchMetric[nCaloJets]/F");  

  // IVF Vertex ID
  // mom
  jetTree_->Branch("jetIVFVertexIDNMom", &jetIVFVertexIDNMom, "jetIVFVertexIDNMom/I");
  jetTree_->Branch("jetIVFVertexIDMom", &jetIVFVertexIDMom, "jetIVFVertexIDMom[jetIVFVertexIDNMom]/I");
  jetTree_->Branch("jetIVFVertexIDMomPt", &jetIVFVertexIDMomPt, "jetIVFVertexIDMomPt[jetIVFVertexIDNMom]/F");
  jetTree_->Branch("jetIVFVertexIDMomJetID", &jetIVFVertexIDMomJetID, "jetIVFVertexIDMomJetID[jetIVFVertexIDNMom]/I");
  // highest pt mom (jet indexed)
  jetTree_->Branch("jetIVFVertexIDMomHighestPtID", &jetIVFVertexIDMomHighestPtID, "jetIVFVertexIDMomHighestPtID[nCaloJets]/I");
  jetTree_->Branch("jetIVFVertexIDMomHighestPt", &jetIVFVertexIDMomHighestPt, "jetIVFVertexIDMomHighestPt[nCaloJets]/F");
  // son
  jetTree_->Branch("jetIVFVertexIDNSon", &jetIVFVertexIDNSon, "jetIVFVertexIDNSon/I");
  jetTree_->Branch("jetIVFVertexIDSon", &jetIVFVertexIDSon, "jetIVFVertexIDSon[jetIVFVertexIDNSon]/I");
  jetTree_->Branch("jetIVFVertexIDSonPt", &jetIVFVertexIDSonPt, "jetIVFVertexIDSonPt[jetIVFVertexIDNSon]/F");
  jetTree_->Branch("jetIVFVertexIDSonJetID", &jetIVFVertexIDSonJetID, "jetIVFVertexIDSonJetID[jetIVFVertexIDNSon]/I");

  // Gen Matching
  // jetTree_->Branch("jetSvNGenMatched", &jetSvNGenMatched, "jetSvNGenMatched/I");  
  // jetTree_->Branch("jetSvNGenFake", &jetSvNGenFake, "jetSvNGenFake/I");  
  // gen
  jetTree_->Branch("jetSvGenVertexMatched", &jetSvGenVertexMatched, "jetSvGenVertexMatched[nCaloJets]/I");  
  jetTree_->Branch("jetSvGenVertexMatchMetric", &jetSvGenVertexMatchMetric, "jetSvGenVertexMatchMetric[nCaloJets]/F");  
  // sim
  // jetTree_->Branch("jetSvNSimMatched", &jetSvNSimMatched, "jetSvNSimMatched/I");  
  // jetTree_->Branch("jetSvNSimFake", &jetSvNSimFake, "jetSvNSimFake/I");  
  // jetTree_->Branch("jetSvSimMatched", &jetSvSimVertexMatched, "jetSvSimMatched[nCaloJets]/I");  
  // jetTree_->Branch("jetSvSimMatchMetric", &jetSvSimVertexMatchMetric, "jetSvSimMatchMetric[nCaloJets]/F");  

  ///////////  ///////////  ///////////  ///////////  ///////////  ///////////  ////
  //////////////////////////////// TRACK TREE QUANITIES ////////////////////////////
  ///////////  ///////////  ///////////  ///////////  ///////////  /////////// /////
  //////// Everything is either a flat number or indexed by nCaloJets //////////////


  // indices which index the branches 
  trackTree_->Branch("nCaloJets", &nCaloJets, "nCaloJets/I");
  trackTree_->Branch("nLiTracks", &nLiTracks, "nLiTracks/I");
  trackTree_->Branch("nSvJets", &nSvJets, "nSvJets/I");
  trackTree_->Branch("nSV", &nSV, "nSV/I");
  trackTree_->Branch("nLiJets", &nLiJets, "nLiJets/I");

  // analyzer event number
  trackTree_->Branch("evNum", &evNum, "evNum/I");

  // file run numbers
  trackTree_->Branch("run", &run, "run/I");
  trackTree_->Branch("lumi", &lumi, "lumi/I");
  trackTree_->Branch("event", &event, "event/I");

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
  trackTree_->Branch("genMatchTrack", &genMatchTrack, "genMatchTrack[nLiTracks]/I");  
  trackTree_->Branch("caloGenPt", &caloGenPt, "caloGenPt[nCaloJets]/F");  
  trackTree_->Branch("caloGenEta", &caloGenEta, "caloGenEta[nCaloJets]/F");  
  trackTree_->Branch("caloGenPhi", &caloGenPhi, "caloGenPhi[nCaloJets]/F");  
  trackTree_->Branch("caloGenM", &caloGenM, "caloGenM[nCaloJets]/F");  

  ////////////////////////////// LIFETIME Jet tags////////////////////////

  // // lifetime jets
  // trackTree_->Branch("liJetID", &liJetID, "liJetID[nLiJets]/I");
  // trackTree_->Branch("liJetPt", &liJetPt, "liJetPt[nLiJets]/F");
  // trackTree_->Branch("liJetPhi", &liJetPhi, "liJetPhi[nLiJets]/F");
  // trackTree_->Branch("liJetEta", &liJetEta, "liJetEta[nLiJets]/F");
  // trackTree_->Branch("liJetNSelTracks", &liJetNSelTracks, "liJetNSelTracks[nLiJets]/I");  

  // // lifetime track info
  // trackTree_->Branch("liTrackEta", &liTrackEta, "liTrackEta[nLiTracks]/F");
  // trackTree_->Branch("liTrackPhi", &liTrackPhi, "liTrackPhi[nLiTracks]/F");
  // trackTree_->Branch("liTrackPt", &liTrackPt, "liTrackPt[nLiTracks]/F");
  // trackTree_->Branch("liTrackJetID", &liTrackJetID, "liTrackJetID[nLiTracks]/I");  
  // trackTree_->Branch("liJetTrackDR", &liJetTrackDR, "liJetTrackDR[nLiTracks]/F");  

  // // lifetime ip tag info
  // trackTree_->Branch("liTrackIP2D", &liTrackIP2D, "liTrackIP2D[nLiTracks]/F");
  // trackTree_->Branch("liTrackIPSig2D", &liTrackIPSig2D, "liTrackIPSig2D[nLiTracks]/F");
  // trackTree_->Branch("liTrackIP3D", &liTrackIP3D, "liTrackIP3D[nLiTracks]/F");
  // trackTree_->Branch("liTrackIPSig3D", &liTrackIPSig3D, "liTrackIPSig3D[nLiTracks]/F");
  // trackTree_->Branch("liTrackDistanceJetAxis", &liTrackDistanceJetAxis, "liTrackDistanceJetAxis[nLiTracks]/F");
  // trackTree_->Branch("liTrackDistanceJetAxisSig", &liTrackDistanceJetAxisSig, "liTrackDistanceJetAxisSig[nLiTracks]/F");

  //////////////////////////////SV Jet tags////////////////////////

  //sv jets
  // trackTree_->Branch("svJetPt", &svJetPt, "svJetPt[nSvJets]/F");
  // trackTree_->Branch("svJetEta", &svJetEta, "svJetEta[nSvJets]/F");
  // trackTree_->Branch("svJetPhi", &svJetPhi, "svJetPhi[nSvJets]/F");
  // trackTree_->Branch("svJetID", &svJetID, "svJetID[nSvJets]/I");

  // quality
  trackTree_->Branch("svChi2", &svChi2, "svChi2[nSV]/F");
  //  trackTree_->Branch("svNChi2", &svNChi2, "svNChi2[nSV]/F"); // what is this?
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
  //////////////////////////////// VTX TREE QUANITIES //////////////////////////////
  ///////////  ///////////  ///////////  ///////////  ///////////  /////////// /////
  /// indices of branches are NOT related between vertex types //////
  
  // arbitrary evNum from running the plugin onces
  vertexTree_->Branch("evNum", &evNum, "evNum/I");

  // file run numbers
  vertexTree_->Branch("run", &run, "run/I");
  vertexTree_->Branch("lumi", &lumi, "lumi/I");
  vertexTree_->Branch("event", &event, "event/I");

  // primary vertices
  vertexTree_->Branch("pvN", &pvN, "pvN/I");
  vertexTree_->Branch("pvMass", &pvMass, "pvMass[pvN]/F");
  vertexTree_->Branch("pvNTrack", &pvNTrack, "pvNTrack[pvN]/I");
  vertexTree_->Branch("pvSumPtSq", &pvSumPtSq, "pvSumPtSq[pvN]/F");

  // position measurements
  vertexTree_->Branch("pvX", &pvX, "pvX[pvN]/F");
  vertexTree_->Branch("pvY", &pvY, "pvY[pvN]/F");
  vertexTree_->Branch("pvZ", &pvZ, "pvZ[pvN]/F");
  vertexTree_->Branch("pvXErr", &pvXErr, "pvXErr[pvN]/F");
  vertexTree_->Branch("pvYErr", &pvYErr, "pvYErr[pvN]/F");
  vertexTree_->Branch("pvZErr", &pvZErr, "pvZErr[pvN]/F");

  // quality
  vertexTree_->Branch("pvChi2", &pvMass, "pvMass[pvN]/F");

  // Standalone Vertexing
  vertexTree_->Branch("vtxN", &vtxN, "vtxN/I");
  vertexTree_->Branch("vtxIsFake", &vtxIsFake, "vtxIsFake[vtxN]/I");
  vertexTree_->Branch("vtxNTracks", &vtxNTracks, "vtxNTracks[vtxN]/I");
  vertexTree_->Branch("vtxChi2", &vtxChi2, "vtxChi2[vtxN]/F");
  vertexTree_->Branch("vtxNDof", &vtxNDof, "vtxNDof[vtxN]/I");
  vertexTree_->Branch("vtxX", &vtxX, "vtxX[vtxN]/F");
  vertexTree_->Branch("vtxY", &vtxY, "vtxY[vtxN]/F");
  vertexTree_->Branch("vtxZ", &vtxZ, "vtxZ[vtxN]/F");
  vertexTree_->Branch("vtxLxy", &vtxLxy, "vtxLxy[vtxN]/F");
  vertexTree_->Branch("vtxLxyz", &vtxLxyz, "vtxLxyz[vtxN]/F");
  vertexTree_->Branch("vtxXsig", &vtxXSig, "vtxXSig[vtxN]/F");
  vertexTree_->Branch("vtxYsig", &vtxYSig, "vtxYSig[vtxN]/F");
  vertexTree_->Branch("vtxZsig", &vtxZSig, "vtxZSig[vtxN]/F");
  vertexTree_->Branch("vtxLxySig", &vtxLxySig, "vtxLxySig[vtxN]/F");

  // Inclusive candidate vertices
  vertexTree_->Branch("vtxIncCandN", &vtxIncCandN, "vtxIncCandN/I");
  vertexTree_->Branch("vtxIncCandIsFake", &vtxIncCandIsFake, "vtxIncCandIsFake[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandNTracks", &vtxIncCandNTracks, "vtxIncCandNTracks[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandChi2", &vtxIncCandChi2, "vtxIncCandChi2[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandNDof", &vtxIncCandNDof, "vtxIncCandNDof[vtxIncCandN]/I");
  vertexTree_->Branch("vtxIncCandX", &vtxIncCandX, "vtxIncCandX[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandY", &vtxIncCandY, "vtxIncCandY[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandZ", &vtxIncCandZ, "vtxIncCandZ[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandLxy", &vtxIncCandLxy, "vtxIncCandLxy[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandLxyz", &vtxIncCandLxyz, "vtxIncCandLxyz[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandXsig", &vtxIncCandXSig, "vtxIncCandXSig[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandYsig", &vtxIncCandYSig, "vtxIncCandYSig[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandZsig", &vtxIncCandZSig, "vtxIncCandZSig[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandLxySig", &vtxIncCandLxySig, "vtxIncCandLxySig[vtxIncCandN]/F");

  vertexTree_->Branch("vtxIncCandNGenMatched", &vtxIncCandNGenMatched, "vtxIncCandNGenMatched/I");
  vertexTree_->Branch("vtxIncCandNGenFake", &vtxIncCandNGenFake, "vtxIncCandNGenFake/I");
  vertexTree_->Branch("vtxIncCandGenMatched", &vtxIncCandGenMatched, "vtxIncCandGenMatched[vtxIncCandN]/I");
  vertexTree_->Branch("vtxIncCandGenMatchMetric", &vtxIncCandGenMatchMetric, "vtxIncCandGenMatchMetric[vtxIncCandN]/F");
  vertexTree_->Branch("vtxIncCandNSimMatched", &vtxIncCandNSimMatched, "vtxIncCandNSimMatched/I");
  vertexTree_->Branch("vtxIncCandNSimFake", &vtxIncCandNSimFake, "vtxIncCandNSimFake/I");
  //  vertexTree_->Branch("vtxIncCandSimMatched", &vtxIncCandSimMatched, "vtxIncCandSimMatched[vtxIncCandN]/I");
  // vertexTree_->Branch("vtxIncCandSimMatchMetric", &vtxIncCandSimMatchMetric, "vtxIncCandSimMatchMetric[vtxIncCandN]/F");

  //inclusive secondary vertices (after merging and arbitration) 
  vertexTree_->Branch("vtxIncSecN", &vtxIncSecN, "vtxIncSecN/I");
  vertexTree_->Branch("vtxIncSecIsFake", &vtxIncSecIsFake, "vtxIncSecIsFake[vtxIncSecN]/F");
  vertexTree_->Branch("vtxIncSecNTracks", &vtxIncSecNTracks, "vtxIncSecNTracks[vtxIncSecN]/I");
  vertexTree_->Branch("vtxIncSecChi2", &vtxIncSecChi2, "vtxIncSecChi2[vtxIncSecN]/F");
  vertexTree_->Branch("vtxIncSecNDof", &vtxIncSecNDof, "vtxIncSecNDof[vtxIncSecN]/I");
  vertexTree_->Branch("vtxIncSecX", &vtxIncSecX, "vtxIncSecX[vtxIncSecN]/F");
  vertexTree_->Branch("vtxIncSecY", &vtxIncSecY, "vtxIncSecY[vtxIncSecN]/F");
  vertexTree_->Branch("vtxIncSecZ", &vtxIncSecZ, "vtxIncSecZ[vtxIncSecN]/F");
  vertexTree_->Branch("vtxIncSecLxy", &vtxIncSecLxy, "vtxIncSecLxy[vtxIncSecN]/F");
  vertexTree_->Branch("vtxIncSecLxyz", &vtxIncSecLxyz, "vtxIncSecLxyz[vtxIncSecN]/F");
  vertexTree_->Branch("vtxIncSecXsig", &vtxIncSecXSig, "vtxIncSecXSig[vtxIncSecN]/F");
  vertexTree_->Branch("vtxIncSecYsig", &vtxIncSecYSig, "vtxIncSecYSig[vtxIncSecN]/F");
  vertexTree_->Branch("vtxIncSecZsig", &vtxIncSecZSig, "vtxIncSecZSig[vtxIncSecN]/F");
  vertexTree_->Branch("vtxIncSecLxySig", &vtxIncSecLxySig, "vtxIncSecLxySig[vtxIncSecN]/F");    

  vertexTree_->Branch("vtxIncSecNGenMatched", &vtxIncSecNGenMatched, "vtxIncSecNGenMatched/I");
  vertexTree_->Branch("vtxIncSecNGenFake", &vtxIncSecNGenFake, "vtxIncSecNGenFake/I");
  vertexTree_->Branch("vtxIncSecGenMatched", &vtxIncSecGenMatched, "vtxIncSecGenMatched[vtxIncSecN]/I");
  vertexTree_->Branch("vtxIncSecGenMatchMetric", &vtxIncSecGenMatchMetric, "vtxIncSecGenMatchMetric[vtxIncSecN]/F");
  // vertexTree_->Branch("vtxIncSecNSimMatched", &vtxIncSecNSimMatched, "vtxIncSecNSimMatched/I");
  // vertexTree_->Branch("vtxIncSecNSimFake", &vtxIncSecNSimFake, "vtxIncSecNSimFake/I");
  // vertexTree_->Branch("vtxIncSecSimMatched", &vtxIncSecSimMatched, "vtxIncSecSimMatched[vtxIncSecN]/I");
  // vertexTree_->Branch("vtxIncSecSimMatchMetric", &vtxIncSecSimMatchMetric, "vtxIncSecSimMatchMetric[vtxIncSecN]/F");


  ///////////  ///////////  ///////////  ///////////  ///////////  ///////////  ////
  //////////////////////////////// SIMULATION QUANITIES ////////////////////////////
  ///////////  ///////////  ///////////  ///////////  ///////////  /////////// /////

  // file run numbers
  genTree_->Branch("evNum",  &evNum,  "evNum/I");
  genTree_->Branch("run",  &run,  "run/I");
  genTree_->Branch("lumi",  &lumi,  "lumi/I");
  genTree_->Branch("event",  &event,  "event/I");

  // GEN Particle Information
  genTree_->Branch("genPartN", &genPartN, "genPartN/I");
  genTree_->Branch("genPartPID", &genPartPID, "genPartPID[genPartN]/I");
  genTree_->Branch("genPartStatus", &genPartStatus, "genPartStatus[genPartN]/I");
  genTree_->Branch("genPartPt", &genPartPt, "genPartPt[genPartN]/F");
  genTree_->Branch("genPartEta", &genPartEta, "genPartEta[genPartN]/F");
  genTree_->Branch("genPartPhi", &genPartPhi, "genPartPhi[genPartN]/F");
  genTree_->Branch("genPartVX", &genPartVX, "genPartVX[genPartN]/F");
  genTree_->Branch("genPartVY", &genPartVY, "genPartVY[genPartN]/F");
  genTree_->Branch("genPartVZ", &genPartVZ, "genPartVZ[genPartN]/F");
  genTree_->Branch("genPartVLxy", &genPartVLxy, "genPartVLxy[genPartN]/F");
  genTree_->Branch("genPartVLxyz", &genPartVLxyz, "genPartVLxyz[genPartN]/F");

  // SIM Vertex Quantites
  genTree_->Branch("simVtxN", &simVtxN, "simVtxN/I");
  genTree_->Branch("simVtxProcType", &simVtxProcType, "simVtxProcType[simVtxN]/I");
  genTree_->Branch("simVtxID", &simVtxID, "simVtxID[simVtxN]/I");
  genTree_->Branch("simVtxTOF", &simVtxTOF, "simVtxTOF[simVtxN]/F");
  genTree_->Branch("simVtxX", &simVtxX, "simVtxX[simVtxN]/F");
  genTree_->Branch("simVtxY", &simVtxY, "simVtxY[simVtxN]/F");
  genTree_->Branch("simVtxZ", &simVtxZ, "simVtxZ[simVtxN]/F");
  genTree_->Branch("simVtxLxy", &simVtxLxy, "simVtxLxy[simVtxN]/F");
  genTree_->Branch("simVtxLxyz", &simVtxLxyz, "simVtxLxyz[simVtxN]/F");    


  ///////////  ///////////  ///////////  ///////////  ///////////  ///////////  ////
  //////////////////////////////// TRACKING QUANITIES // ////////////////////////////
  ///////////  ///////////  ///////////  ///////////  ///////////  /////////// /////

  // nominal kinematics
  trackTree_->Branch("trCharge", &trCharge, "trCharge[nTracks]/F");
  trackTree_->Branch("trQOverP", &trQOverP, "trQOverP[nTracks]/F");
  trackTree_->Branch("trPt", &trPt, "trPt[nTracks]/F");
  trackTree_->Branch("trPtError", &trPtError, "trPtError[nTracks]/F");
  trackTree_->Branch("trEta", &trEta, "trEta[nTracks]/F");
  trackTree_->Branch("trEtaError", &trEtaError, "trEtaError[nTracks]/F");
  trackTree_->Branch("trPhi", &trPhi, "trPhi[nTracks]/F");
  trackTree_->Branch("trPhiError", &trPhiError, "trPhiError[nTracks]/F");

  // tracking angles 
  trackTree_->Branch("trTheta", &trTheta, "trTheta[nTracks]/F");
  trackTree_->Branch("trThetaError", &trThetaError, "trThetaError[nTracks]/F");
  trackTree_->Branch("trThetaSig", &trThetaSig, "trThetaSig[nTracks]/F");
  trackTree_->Branch("trLambda", &trLambda, "trLambda[nTracks]/F");
  trackTree_->Branch("trLambdaError", &trLambdaError, "trLambdaError[nTracks]/F");
  trackTree_->Branch("trLambdaSig", &trLambdaSig, "trLambdaSig[nTracks]/F");

  // impact parameter proxies
  trackTree_->Branch("trDxy", &trDxy, "trDxy[nTracks]/F");
  trackTree_->Branch("trDxyError", &trDxyError, "trDxyError[nTracks]/F");
  trackTree_->Branch("trDxySig", &trDxySig, "trDxySig[nTracks]/F");
  trackTree_->Branch("trDz", &trDz, "trDz[nTracks]/F");
  trackTree_->Branch("trDzError", &trDzError, "trDzError[nTracks]/F");
  trackTree_->Branch("trDzSig", &trDzSig, "trDzSig[nTracks]/F");
  trackTree_->Branch("trDsz", &trDsz, "trDsz[nTracks]/F");
  trackTree_->Branch("trDszError", &trDszError, "trDszError[nTracks]/F");
  trackTree_->Branch("trDszSig", &trDszSig, "trDszSig[nTracks]/F");

  // reference positions
  trackTree_->Branch("trRefX", &trRefX, "trRefX[nTracks]/F");
  trackTree_->Branch("trRefY", &trRefY, "trRefY[nTracks]/F");
  trackTree_->Branch("trRefZ", &trRefZ, "trRefZ[nTracks]/F");

  // inner hit positions
  trackTree_->Branch("trInnerX", &trInnerX, "trInnerX[nTracks]/F");
  trackTree_->Branch("trInnerY", &trInnerY, "trInnerY[nTracks]/F");
  trackTree_->Branch("trInnerZ", &trInnerZ, "trInnerZ[nTracks]/F");

  // inner hit positions
  trackTree_->Branch("trInnerX", &trInnerX, "trInnerX[nTracks]/F");
  trackTree_->Branch("trInnerY", &trInnerY, "trInnerY[nTracks]/F");
  trackTree_->Branch("trInnerZ", &trInnerZ, "trInnerZ[nTracks]/F");
  trackTree_->Branch("trInnerEta", &trInnerEta, "trInnerEta[nTracks]/F");
  trackTree_->Branch("trInnerPhi", &trInnerPhi, "trInnerPhi[nTracks]/F");
  trackTree_->Branch("trInnerPt", &trInnerPt, "trInnerPt[nTracks]/F");
  trackTree_->Branch("trInnerPx", &trInnerPt, "trInnerPx[nTracks]/F");
  trackTree_->Branch("trInnerPz", &trInnerPt, "trInnerPy[nTracks]/F");
  trackTree_->Branch("trInnerP", &trInnerP, "trInnerP[nTracks]/F");

  // outer hit kinematics
  trackTree_->Branch("trOuterX", &trOuterX, "trOuterX[nTracks]/F");
  trackTree_->Branch("trOuterY", &trOuterY, "trOuterY[nTracks]/F");
  trackTree_->Branch("trOuterZ", &trOuterZ, "trOuterZ[nTracks]/F");
  trackTree_->Branch("trOuterEta", &trOuterEta, "trOuterEta[nTracks]/F");
  trackTree_->Branch("trOuterPhi", &trOuterPhi, "trOuterPhi[nTracks]/F");

  // outer hit kinematics
  trackTree_->Branch("trOuterPt", &trOuterPt, "trOuterPt[nTracks]/F");
  trackTree_->Branch("trOuterPx", &trOuterPt, "trOuterPx[nTracks]/F");
  trackTree_->Branch("trOuterPz", &trOuterPt, "trOuterPy[nTracks]/F");
  trackTree_->Branch("trOuterP", &trOuterP, "trOuterP[nTracks]/F");

  // track quality
  trackTree_->Branch("trChi2", &trChi2, "trChi2[nTracks]/F");
  trackTree_->Branch("trNDof", &trNDoF, "trNDoF[nTracks]/F");
  trackTree_->Branch("trNChi2", &trNChi2, "trNChi2[nTracks]/F");
  trackTree_->Branch("trValidFraction", &trValidFraction, "trValidFraction[nTracks]/F");
  trackTree_->Branch("trNLost", &trNLost, "trNLost[nTracks]/I");
  trackTree_->Branch("trNFound", &trNFound, "trNFound[nTracks]/I");
  //trackTree_->Branch("trAlgo", &trAlgo, "trAlgo[nTracks]/F");
  trackTree_->Branch("trAlgoInt", &trAlgoInt, "trAlgoInt[nTracks]/I");

}

void DJetAnalyzer::endJob() 
{
  outputFile_->cd();
  runStatTree_->Write();
  jetTree_->Write();
  trackTree_->Write();
  eventTree_->Write();
  // genTree_->Write(); TODO
  vertexTree_->Write(); 
  outputFile_->Close();
}

void DJetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void DJetAnalyzer::dumpPreSelection(DisplacedJetEvent & djEvent) {
  if(debug > 1 ) std:: cout << "[DEBUG] PV Vertex Dumping" << std::endl;
  eventPassEventPreSelection = djEvent.doesPassPreSelection() ? 1 : 0; 

  const DisplacedJetCollection djetCollection = djEvent.getDisplacedJets();
  DisplacedJetCollection::const_iterator djet = djetCollection.begin();
  int jj = 0;
  for(; djet != djetCollection.end(); ++djet, ++jj) {        
    jetPassPreSelection[jj]	= djet->doesPassPreSelection() ? 1 : 0;
    eventNJetsPassPreSelection += djet->doesPassPreSelection() ? 1 : 0;
  }
}

// dump all information related to primary vertices into the event
void DJetAnalyzer::dumpPVInfo(DisplacedJetEvent & djEvent, const reco::VertexCollection & pv) {

  if(debug > 1 ) std:: cout << "[DEBUG] PV Vertex Dumping" << std::endl;
  // dump the PV information
  reco::VertexCollection::const_iterator iterPV = pv.begin();
  pvN = 0;
  for(; iterPV != pv.end(); ++iterPV, ++pvN) {
    float x  = iterPV->x(), y = iterPV->y(), z = iterPV->z();
    float xE = iterPV->xError(), yE = iterPV->yError(), zE = iterPV->zError();

    // absolute position
    pvX[pvN]	= x;
    pvY[pvN]	= y;
    pvZ[pvN]	= z;
    pvXErr[pvN] = xE;
    pvYErr[pvN] = yE;
    pvZErr[pvN] = zE;    

    // position variables
    pvMass[pvN]	   = 0;
    pvLxy[pvN]	   = std::sqrt( x * x + y * y );
    pvLxyz[pvN]	   = std::sqrt( x * x + y * y + z * z) ;
    pvLxySig[pvN]  = std::sqrt( x * x + y * y ) / std::sqrt(xE * xE + yE * yE);
    pvLxyzSig[pvN] = std::sqrt( x * x + y * y + z * z / std::sqrt(xE * xE + yE * yE + zE * zE));

    // quality
    pvChi2[pvN]   = iterPV->chi2();   

    // track sums
    reco::Vertex::trackRef_iterator vtxIter = iterPV->tracks_begin();
    pvSumPtSq[pvN] = 0;
    pvNTrack[pvN] = 0;
    for(; vtxIter != iterPV->tracks_end(); ++vtxIter) {
      pvNTrack[pvN]++;
      pvSumPtSq[pvN] += (*vtxIter)->pt() * (*vtxIter)->pt(); 
    }
  }
}

void DJetAnalyzer::dumpCaloInfo(DisplacedJetEvent & djEvent) {
  if(debug > 1 ) std:: cout << "[DEBUG] Dumping Calo Info" << std::endl;

  // calo HT
  eventCaloHT  = djEvent.caloHT; 
  // calo MET
  eventCaloMET = djEvent.caloMET;
  

  const DisplacedJetCollection djetCollection = djEvent.getDisplacedJets();
  DisplacedJetCollection::const_iterator djet = djetCollection.begin();
  int jj = 0;
  for(; djet != djetCollection.end(); ++djet, ++jj) {        
    // id
    caloJetID[jj]	 = djet->jetID;
    // kinematics
    caloJetPt[jj]	 = djet->caloPt;
    caloJetEta[jj]	 = djet->caloEta;
    caloJetPhi[jj]	 = djet->caloPhi;        
    // energy fractions
    caloJetHfrac[jj]	 = djet->caloEMEnergyFrac;
    caloJetEfrac[jj]	 = djet->caloHadEnergyFrac;
    // jet size
    caloJetN60[jj]	 = djet->caloN60;
    caloJetN90[jj]	 = djet->caloN90;
    caloJetTowerArea[jj] = djet->caloTowerArea;
    // gen matching to particle quantities
    caloGenMatch[jj]	 = djet->isCaloGenMatched ? 1 : 0;       
    caloGenPt[jj]	 = djet->caloGenPt;
    caloGenEta[jj]	 = djet->caloGenEta;
    caloGenPhi[jj]	 = djet->caloGenPhi;
  }
}

void DJetAnalyzer::dumpSVTagInfo(DisplacedJetEvent & djEvent) {
  if(debug > 1 ) std:: cout << "[DEBUG] Dumping SV Info" << std::endl;
  const DisplacedJetCollection djetCollection = djEvent.getDisplacedJets();
  DisplacedJetCollection::const_iterator djet = djetCollection.begin();
  int jj = 0;
  for(; djet != djetCollection.end(); ++djet, ++jj) {    
    // positional space    
    jetNSv[jj]                    = djet->svNVertex; 
    jetSvX[jj]			  = djet->svX;
    jetSvY[jj]			  = djet->svY;
    jetSvZ[jj]			  = djet->svZ;
    jetSvXErr[jj]		  = djet->svXError;
    jetSvYErr[jj]		  = djet->svYError;
    jetSvZErr[jj]		  = djet->svZError;
    // qualities
    jetSvMass[jj]		  = djet->svMass;
    jetSvLxy[jj]		  = djet->svLxy;
    jetSvLxySig[jj]		  = djet->svLxySig;
    jetSvLxyz[jj]		  = djet->svLxyz;
    jetSvLxyzSig[jj]		  = djet->svLxyzSig;
    jetSvNTrack[jj]		  = djet->svNTracks;    
    jetSvChi2[jj]                 = djet->svChi2;
    jetSvNDof[jj]                 = djet->svNDof;
    // matching
    jetSvGenVertexMatched[jj]	  = djet->svIsGenMatched;
    jetSvGenVertexMatchMetric[jj] = djet->svGenVertexMatchMetric;
  } // djet loop  
}

void DJetAnalyzer::dumpIVFInfo(DisplacedJetEvent & djEvent) {
  if(debug > 1 ) std:: cout << "[DEBUG] Dumping IVF Info" << std::endl;

  jetIVFVertexIDNMom  = 0;
  jetIVFVertexIDNSon  = 0;

  eventNIVFReco		= djEvent.nIVFReconstructed;
  eventNIVFRecoGenMatch = djEvent.nIVFGenMatched;

  const DisplacedJetCollection djetCollection = djEvent.getDisplacedJets();
  DisplacedJetCollection::const_iterator djet = djetCollection.begin();
  int jj = 0, jjmm = 0, jjss = 0;
  for(; djet != djetCollection.end(); ++djet, ++jj) {    
    // IVF position
    jetIVFX[jj]			     = djet->ivfX;
    jetIVFY[jj]			     = djet->ivfY;
    jetIVFZ[jj]			     = djet->ivfZ;
    jetIVFXErr[jj]		     = djet->ivfXError;
    jetIVFYErr[jj]		     = djet->ivfYError;
    jetIVFZErr[jj]		     = djet->ivfZError;    
    // qualities
    jetIVFNTrack[jj]		     = djet->ivfNTracks;
    jetIVFMass[jj]		     = djet->ivfMass;
    jetIVFLxySig[jj]		     = djet->ivfLxySig;
    jetIVFLxyzSig[jj]		     = djet->ivfLxyzSig;
    jetIVFLxy[jj]		     = djet->ivfLxy;
    jetIVFLxyz[jj]		     = djet->ivfLxyz;    
    // track based matching score of vertex to the jet
    jetIVFMatchingScore[jj]	     = djet->ivfMatchingScore;

    // gen matching of the vertex position to the IVF chosen
    jetIVFGenVertexMatched[jj]	     = djet->ivfIsGenMatched;
    jetIVFGenVertexMatchMetric[jj]   = djet->ivfGenVertexMatchMetric;      
    // sim matching to vertex from sim
    jetIVFSimVertexMatched[jj]	     = djet->ivfIsSimMatched;

    // initializers for IVF ID
    jetIVFVertexIDMomHighestPtID[jj] = 0;

    // highest mom pt and id
    jetIVFVertexIDMomHighestPtID[jj] = djet->ivfHighestPtMomID;
    jetIVFVertexIDMomHighestPt[jj]   = djet->ivfHighestPtMomPt;

    // loop over the vertex ID information
    // moms
    std::vector<std::pair<const int, const float>> moms = djet->genMomVector;
    int momsize = moms.size();
    jetIVFVertexIDNMom += momsize;    
    for(int mm = 0; mm < momsize; ++mm, ++jjmm) { 
      if(debug > 5 ) std:: cout << "[DEBUG] WRITING IVF VERTEX INFO jjmm: " << jjmm << " mm " << mm <<  std::endl;
      jetIVFVertexIDMom[jjmm]   = moms[mm].first;
      jetIVFVertexIDMomPt[jjmm] = moms[mm].second;
      jetIVFVertexIDMomJetID[jjmm] = djet->jetID;
    }
    // sons
    std::vector<std::pair<const int, const float>> sons = djet->genSonVector;
    int sonsize = sons.size();
    jetIVFVertexIDNSon += sonsize;
    for(int ss = 0; ss < sonsize; ++ss, ++jjss) {
      if(debug > 5 ) std:: cout << "[DEBUG] WRITING IVF VERTEX INFO jjss: " << jjss << " ss " << ss <<  std::endl;
      jetIVFVertexIDSon[jjss]	   = sons[ss].first;
      jetIVFVertexIDSonPt[jjss]	   = sons[ss].second;
      jetIVFVertexIDSonJetID[jjss] = djet->jetID;
    }
  }
}

void DJetAnalyzer::dumpIPInfo(DisplacedJetEvent & djEvent) {
  if(debug > 1 ) std:: cout << "[DEBUG] Dumping IP Info" << std::endl;

  const DisplacedJetCollection djetCollection = djEvent.getDisplacedJets();
  DisplacedJetCollection::const_iterator djet = djetCollection.begin();
  int jj = 0;
  for(; djet != djetCollection.end(); ++djet, ++jj) {    

    if(debug > 4 ) std:: cout << "[DEBUG 4] djet->ipsiglogsum2d: " << djet->ipSigLogSum2D << std::endl;

    // jet track association
    jetNTracks[jj]	      = djet->nTracks;
    // IP significance log sums
    jetIPSigLogSum2D[jj]      = djet->ipSigLogSum2D;
    jetIPSigLogSum3D[jj]      = djet->ipSigLogSum3D;    
    // distance from jet axis
    jetDistLogSum[jj]	      = djet->jetDistLogSum;
    jetDistSigLogSum[jj]      = djet->jetDistSigLogSum;
    //2d eip
    jetEIPSig2D[jj]	      = djet->eipSigSum2D;
    jetEIPSigLog2D[jj]	      = djet->eipSigLogSum2D;
    // 3d eip
    jetEIPSig3D[jj]	      = djet->eipSigSum3D;
    jetEIPSigLog3D[jj]	      = djet->eipSigLogSum3D;
    // IP significance averages
    jetMeanIPSig2D[jj]	      = djet->meanIPSig2D;
    jetMeanIPSig3D[jj]	      = djet->meanIPSig3D;
    jetMeanIPLogSig2D[jj]     = djet->meanIPLogSig2D;
    jetMeanIPLogSig3D[jj]     = djet->meanIPLogSig3D;
    // edian
    jetMedianIPSig2D[jj]      = djet->medianIPSig2D;
    jetMedianIPSig3D[jj]      = djet->medianIPSig3D;
    jetMedianIPLogSig2D[jj]   = djet->medianIPLogSig2D;
    jetMedianIPLogSig3D[jj]   = djet->medianIPLogSig3D;
    // variance
    jetVarianceIPSig2D[jj]    = djet->varianceIPSig2D;
    jetVarianceIPSig3D[jj]    = djet->varianceIPSig3D;
    jetVarianceIPLogSig2D[jj] = djet->varianceIPLogSig2D;
    jetVarianceIPLogSig3D[jj] = djet->varianceIPLogSig3D;
    jetVarianceJetDistSig[jj] = djet->varianceJetDistSig;    

    // IP value averages
    // mean
    jetMeanIP2D[jj]	   = djet->meanIP2D;
    jetMeanIP3D[jj]	   = djet->meanIP3D;
    jetMeanIPLog2D[jj]	   = djet->meanIPLog2D;
    jetMeanIPLog3D[jj]	   = djet->meanIPLog3D;
    jetMeanJetDist[jj]     = djet->meanJetDist;
    // median
    jetMedianIP2D[jj]	   = djet->medianIP2D;
    jetMedianIP3D[jj]	   = djet->medianIP3D;
    jetMedianIPLog2D[jj]   = djet->medianIPLog2D;
    jetMedianIPLog3D[jj]   = djet->medianIPLog3D;
    jetMedianJetDist[jj]   = djet->medianJetDist;
    // variance
    jetVarianceIP2D[jj]	   = djet->varianceIP2D;
    jetVarianceIP3D[jj]	   = djet->varianceIP3D;
    jetVarianceIPLog2D[jj] = djet->varianceIPLog2D;
    jetVarianceIPLog3D[jj] = djet->varianceIPLog3D;
    jetVarianceJetDist[jj] = djet->varianceJetDist;
  }
}


void
DJetAnalyzer::dumpSimInfo(const edm::SimVertexContainer & simVtx) {

  if(debug > 1 ) std::cout << "[DEBUG] Sim Vertex Dumping" << std::endl;

  simVtxN = 0;
  edm::SimVertexContainer::const_iterator	iterSimVtx = simVtx.begin();
  for(; iterSimVtx != simVtx.end(); ++iterSimVtx){    
	
    Int_t	proc_type = iterSimVtx->processType();
	
    // only keep vertices from the generator level
    if (proc_type != 0) continue;
	
    simVtxProcType[simVtxN] = proc_type;
    simVtxID[simVtxN]	    = iterSimVtx->vertexId();
	
    const math::XYZTLorentzVectorD &    pos = iterSimVtx->position();      
    simVtxTOF[simVtxN]			= pos.t();
    simVtxX[simVtxN]			= pos.x();
    simVtxY[simVtxN]			= pos.y();
    simVtxZ[simVtxN]			= pos.z();

    simVtxLxy[simVtxN]  = std::sqrt(pos.x() * pos.x()  + pos.y() * pos.y() );
    simVtxLxyz[simVtxN] = std::sqrt(pos.x() * pos.x()  + pos.y() * pos.y() + pos.z() * pos.z());
	
    simVtxN++;
  }
}

// dump all of the track info available at AOD into the branches and label with a collection ID
void DJetAnalyzer::dumpTrackInfo(DisplacedJetEvent djEvent, reco::TrackCollection tracks, const int & collectionID) {
  reco::TrackCollection::const_iterator tt = tracks.begin();
  for(; tt != tracks.end(); ++tt) {
    // nominal kinematics
    trCollectionID[nTracks]  = collectionID; 
    trCharge[nTracks]	     = tt->charge();
    trQOverP[nTracks]	     = tt->qoverp();
    trPt[nTracks]	     = tt->pt();
    trEta[nTracks]	     = tt->eta();
    trPhi[nTracks]	     = tt->phi();
    //tracking angles
    trTheta[nTracks]	     = 0;
    trThetaError[nTracks]    = 0;
    trThetaSig[nTracks]	     = 0;
    trLambda[nTracks]	     = 0;
    trLambdaError[nTracks]   = 0;
    trLambdaSig[nTracks]     = 0;

    // impact parameter proxies
    trDxy[nTracks]	     = tt->dxy();
    trDxyError[nTracks]	     = 0;
    trDxySig[nTracks]	     = 0;
    trDz[nTracks]	     = 0;
    trDzError[nTracks]	     = 0;
    trDzSig[nTracks]	     = 0;
    trDsz[nTracks]	     = 0;
    trDszError[nTracks]	     = 0;
    trDszSig[nTracks]	     = 0;

    // reference point
    trRefX[nTracks]	     = tt->vx();
    trRefY[nTracks]	     = tt->vy();
    trRefZ[nTracks]	     = tt->vz();

    // inner positions
    trInnerX[nTracks]	     = 0;
    trInnerY[nTracks]	     = 0;
    trInnerZ[nTracks]	     = 0;
    trInnerEta[nTracks]	     = 0;
    trInnerPhi[nTracks]	     = 0;

    // inne momentum
    trInnerPt[nTracks]	     = 0;
    trInnerPx[nTracks]	     = 0;
    trInnerPy[nTracks]	     = 0;
    trInnerPz[nTracks]	     = 0;
    trInnerP[nTracks]	     = 0;

    // outer position
    trOuterX[nTracks]	     = 0;
    trOuterY[nTracks]	     = 0;
    trOuterZ[nTracks]	     = 0;
    trOuterEta[nTracks]	     = 0;
    trOuterPhi[nTracks]	     = 0;
    trOuterRadius[nTracks]   = 0;

    // outer momentum
    trOuterPt[nTracks]	     = 0;
    trOuterPx[nTracks]	     = 0;
    trOuterPy[nTracks]	     = 0;
    trOuterPz[nTracks]	     = 0;
    trOuterP[nTracks]	     = 0;

    // quality
    trChi2[nTracks]	     = tt->chi2();
    trNDoF[nTracks]	     = tt->ndof();
    trNChi2[nTracks]	     = tt->normalizedChi2();
    trValidFraction[nTracks] = 0;
    trNLost[nTracks]	     = 0;
    trNFound[nTracks]	     = 0;

    //    trAlgo[nTracks]	     = ;
    trAlgoInt[nTracks]	     = tt->algo();

    nTracks++;
  }
    
  
}

void DJetAnalyzer::dumpGenInfo(const reco::GenParticleCollection & gen) {
  if(debug > 1 ) std::cout << "[DEBUG] Gen Particle Dumping" << std::endl;

  reco::GenParticleCollection::const_iterator iterGenParticle = gen.begin();
  genPartN = 0;
  for(; iterGenParticle != gen.end(); ++iterGenParticle){    
    float vx = iterGenParticle->vx(), vy = iterGenParticle->vy(), vz = iterGenParticle->vz();
    // gen id and status
    genPartPID[genPartN]    = iterGenParticle->pdgId();
    genPartStatus[genPartN] = iterGenParticle->status();
    // gen kinematics
    genPartPt[genPartN]  = iterGenParticle->pt();
    genPartEta[genPartN] = iterGenParticle->eta();
    genPartPhi[genPartN] = iterGenParticle->phi();
    // vertex position
    genPartVX[genPartN] = vx;
    genPartVY[genPartN] = vy;
    genPartVZ[genPartN] = vz;
    // vertex flight
    genPartVLxy[genPartN] = std::sqrt( vx * vx + vy * vy );
    genPartVLxyz[genPartN] = std::sqrt( vx * vx + vy * vy + vz * vz);
    
    genPartN++;    
  }
}

// dumps the djevent tags into trees
void DJetAnalyzer::dumpDJTags(DisplacedJetEvent & djEvent) { 
  if(debug > 1 ) std::cout << "[DEBUG] Dumping Tag INfo" << std::endl;
  nWP = djEvent.nNoVertexTagsVector.size();  
  for(int wp = 0; wp < nWP; ++wp) {
    eventNNoVertexTags[wp] = djEvent.nNoVertexTagsVector[wp];
    eventNShortTags[wp]	   = djEvent.nShortTagsVector[wp];
    eventNMediumTags[wp]   = djEvent.nMediumTagsVector[wp];
    eventNLongTags[wp]     = djEvent.nLongTagsVector[wp];
    eventNTotalTags[wp]    = djEvent.nNoVertexTagsVector[wp] + djEvent.nShortTagsVector[wp] + 
      djEvent.nMediumTagsVector[wp] + djEvent.nLongTagsVector[wp];

    eventCaloDHT[wp]       = djEvent.caloDHT[wp];
  } // loop: working points  
}


//define this as a plug-in
DEFINE_FWK_MODULE(DJetAnalyzer);
