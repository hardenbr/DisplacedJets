// \class JetTracksAssociatorAtInnerHitProducer JetTracksAssociatorAtInnerHitProducer.cc 
//
// Original Author:  Joshua Hardenbrook

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"

#include "JetTracksAssociatorAtInnerHitProducer.h"

JetTracksAssociatorAtInnerHitProducer::JetTracksAssociatorAtInnerHitProducer(const edm::ParameterSet& fConfig)
  : mAssociator (fConfig.getParameter<double> ("coneSize")){

  mJets = consumes<edm::View <reco::Jet> >(fConfig.getParameter<edm::InputTag> ("jets"));
  mTracks = consumes<reco::TrackCollection>(fConfig.getParameter<edm::InputTag>("tracks"));

  //pvSrc = consumes<reco::VertexCollection>(fConfig.getParameter<edm::InputTag> ("pvSrc"));

  produces<reco::JetTracksAssociation::Container> ();
}

JetTracksAssociatorAtInnerHitProducer::~JetTracksAssociatorAtInnerHitProducer() {}

void JetTracksAssociatorAtInnerHitProducer::produce(edm::Event& fEvent, const edm::EventSetup& fSetup) {
  edm::Handle <edm::View <reco::Jet> > jets_h;
  fEvent.getByToken (mJets, jets_h);
  edm::Handle <reco::TrackCollection> tracks_h;
  fEvent.getByToken (mTracks, tracks_h);
  
  std::auto_ptr<reco::JetTracksAssociation::Container> jetTracks (new reco::JetTracksAssociation::Container (reco::JetRefBaseProd(jets_h)));

  // format inputs
  std::vector <edm::RefToBase<reco::Jet> > allJets;
  allJets.reserve (jets_h->size());
  for (unsigned i = 0; i < jets_h->size(); ++i) allJets.push_back (jets_h->refAt(i));
  std::vector <reco::TrackRef> allTracks;
  //reco::TrackRefVector allTracks;
  //  allTracks.reserve (tracks_h->size());

  // add all the tracks
  for (unsigned i = 0; i < tracks_h->size(); ++i) {
    allTracks.push_back (reco::TrackRef (tracks_h, i));
  }
  
  //run the associator
  mAssociator.produce (&*jetTracks, allJets, allTracks, fSetup);

  // store output
  fEvent.put (jetTracks);
}


