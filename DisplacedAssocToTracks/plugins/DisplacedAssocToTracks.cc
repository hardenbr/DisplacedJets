// -*- C++ -*-
//
// Package:    DisplacedJets/DisplacedAssocToTracks
// Class:      DisplacedAssocToTracks
// 
/**\class DisplacedAssocToTracks DisplacedAssocToTracks.cc DisplacedJets/DisplacedAssocToTracks/plugins/DisplacedAssocToTracks.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Robert Hardenbrook
//         Created:  Thu, 26 Mar 2015 07:39:53 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducer.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/BTauReco/interface/VertexTypes.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


//
// class declaration
//

class DisplacedAssocToTracks : public edm::EDProducer {
   public:
      explicit DisplacedAssocToTracks(const edm::ParameterSet&);
      ~DisplacedAssocToTracks();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  
  edm::InputTag tag_jetTracksAssociation_;
  edm::InputTag tag_genParticles_;   
  std::string outputLabel_;
  float jetPtCut_;
  int debug_;
  bool isSignalMC_;
  bool doGenMatch_;


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
DisplacedAssocToTracks::DisplacedAssocToTracks(const edm::ParameterSet& iConfig)
{

  tag_jetTracksAssociation_ = iConfig.getUntrackedParameter<edm::InputTag>("jetTracksAssociation");
  tag_genParticles_ = iConfig.getUntrackedParameter<edm::InputTag>("genParticleTag");
  jetPtCut_ = iConfig.getUntrackedParameter<double>("jetPtCut");
  outputLabel_ = iConfig.getUntrackedParameter<std::string>("outputLabel");
  //isSignalMC_    =   iConfig.getUntrackedParameter<bool>("isSignalMC");
  //doGenMatch_    =   iConfig.getUntrackedParameter<bool>("doGenMatch");
  debug_    =   iConfig.getUntrackedParameter<int>("debug");

  produces<reco::TrackCollection >(outputLabel_);    
}


DisplacedAssocToTracks::~DisplacedAssocToTracks()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DisplacedAssocToTracks::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   // Handle<reco::SecondaryVertexTagInfoCollection> secondaryVertexTagInfo;
   // iEvent.getByLabel(tag_secondaryVertexTagInfo_, secondaryVertexTagInfo);

   edm::Handle<reco::JetTracksAssociationCollection> jetTracksAssociation;
   iEvent.getByLabel(tag_jetTracksAssociation_, jetTracksAssociation);

   std::auto_ptr<reco::TrackCollection > jetTrackCollectionResult(new reco::TrackCollection);
   reco::TrackCollection jetTrackCollection;                                                                                      

   if (debug_ > 0) std::cout << "[DAssoc2Tracks] starting jet tracks iterator" << std::endl;

   for(reco::JetTracksAssociationCollection::const_iterator it = jetTracksAssociation->begin();
       it != jetTracksAssociation->end(); ++it){
     const reco::Jet * jet = &*it->first; 
     
     //expect only calo jets for track assocations
     if (debug_ > 0) std::cout << "[DAssoc2Tracks] Casting jet to calojet" << std::endl;

     reco::CaloJet const * pCaloJet = dynamic_cast<reco::CaloJet const *>(jet);
     if ( pCaloJet == 0 ) {
       throw cms::Exception("InvalidInput") << "Expecting calo jets only in JetTracksAssociationXtrpCalo";
     }

     float pt = pCaloJet->detectorP4().pt();
     if (debug_ > 0) std::cout << "[DAssoc2Tracks] pt cut" << std::endl;
     if (debug_ > 0) std::cout << "[DAssoc2Tracks] jet pt " << pt << std::endl;
     if (pt < jetPtCut_ ) continue;

     if (it->second.isNull()) { continue; }
     std::vector<reco::Track> tracks = *it->second.product();     
     if (debug_ > 0) std::cout << "[DAssoc2Tracks] Track Iteration" << std::endl;
     for(std::vector<reco::Track>::const_iterator tt = tracks.begin(); tt != tracks.end(); tt++ ){
       jetTrackCollection.push_back(*tt);
     }
   }
   
   if (debug_ > 0) std::cout << "[DAssoc2Tracks] Cast" << std::endl;
   *jetTrackCollectionResult = jetTrackCollection;

   if (debug_ > 0) std::cout << "[DAssoc2Tracks] Final Put" << std::endl;
   iEvent.put(jetTrackCollectionResult, outputLabel_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
DisplacedAssocToTracks::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DisplacedAssocToTracks::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
DisplacedAssocToTracks::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DisplacedAssocToTracks::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DisplacedAssocToTracks::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DisplacedAssocToTracks::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisplacedAssocToTracks::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedAssocToTracks);
