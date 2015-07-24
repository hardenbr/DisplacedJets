// \class JetTracksAssociatorAtInnerHitProducer JetTracksAssociatorAtInnerHitProducer.cc 
//
// Original Author:  Andrea Rizzi
//         Created:  Wed Apr 12 11:12:49 CEST 2006
// Accommodated for Jet Package by: Fedor Ratnikov Jul. 30, 2007
//
//
#ifndef JetTracksAssociatorAtInnerHitProducer_h
#define JetTracksAssociatorAtInnerHitProducer_h

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/Common/interface/EDProductfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "JetTracksAssociationInnerHit.h"
//#include "RecoJets/JetAssociationAlgorithms/interface/JetTracksAssociationDRVertex.h"
//#include "RecoJets/JetAssociationAlgorithms/interface/JetTracksAssociationDRVertexAssigned.h"

class JetTracksAssociatorAtInnerHitProducer : public edm::stream::EDProducer<> {
   public:
      JetTracksAssociatorAtInnerHitProducer(const edm::ParameterSet&);
      virtual ~JetTracksAssociatorAtInnerHitProducer();

      virtual void produce(edm::Event&, const edm::EventSetup&);

   private:
     edm::EDGetTokenT<edm::View <reco::Jet>> mJets;
     edm::EDGetTokenT<reco::TrackCollection> mTracks;

     JetTracksAssociationInnerHit mAssociator;
     //JetTracksAssociationDRVertexAssigned mAssociatorAssigned;

     //     edm::EDGetTokenT<reco::VertexCollection> pvSrc; /// if useAssigned, will read this PV collection. 
};

DEFINE_FWK_MODULE(JetTracksAssociatorAtInnerHitProducer);

#endif

