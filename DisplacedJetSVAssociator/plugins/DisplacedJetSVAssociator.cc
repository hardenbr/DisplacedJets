// -*- C++ -*-
//
// Package:    DisplacedJets/DisplacedJetSVAssociator
// Class:      DisplacedJetSVAssociator
// 
/**\class DisplacedJetSVAssociator DisplacedJetSVAssociator.cc DisplacedJets/DisplacedJetSVAssociator/plugins/DisplacedJetSVAssociator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Robert Hardenbrook
//         Created:  Sun, 26 Apr 2015 14:08:28 GMT
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
#include "DisplacedJets/DisplacedJetSVAssociator/interface/JetVertexAssociation.h"

//
// class declaration
//

class DisplacedJetSVAssociator : public edm::EDProducer {
public:
  explicit DisplacedJetSVAssociator(const edm::ParameterSet&);
  ~DisplacedJetSVAssociator();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  typedef std::vector<reco::JetVertexAssociation> jetVertexAssocationCollection; 
  
private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  

  edm::InputTag tag_caloJets_;
  edm::InputTag tag_secondaryVertices_;
  edm::InputTag tag_primaryVertices_;
  std::string outputLabel_;


  int debug_;
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
DisplacedJetSVAssociator::DisplacedJetSVAssociator(const edm::ParameterSet& iConfig) {
  tag_caloJets_ = iConfig.getUntrackedParameter<edm::InputTag>("caloJets");
  tag_secondaryVertices_ = iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertices");
  tag_primaryVertices_ = iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices");
  tag_algoName_ = iConfig.getUntrackedParameter<std::string>("algoName");
  outputLabel_ = iConfig.getUntrackedParameter<std::string>("outputLabel");

  produces<jetVertexAssociationCollection>(outputLabel_);

}


DisplacedJetSVAssociator::~DisplacedJetSVAssociator()
{

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DisplacedJetSVAssociator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   std::auto_ptr<jetVertexAssociationCollection> jetVertexAssociationCollection(new jetVertexAssociationCollection);
      
   edm::Handle<reco::CaloJetCollection> caloJets;
   edm::Handle<reco::VertexCollection>	secondaryVertices;
   edm::Handle<reco::VertexCollection>	primaryVertices;

   //get the products from the tags
   iEvent.GetByLabel(tag_caloJets_, caloJets);
   iEvent.GetByLabel(tag_secondaryVertices_, secondaryVertices);
   iEvent.GetByLabel(tag_primaryVertices_, primaryVertices);
   iEvent.GetByLabel(tag_algoName_, algoName);

   const reco::CaloJetCollection &  jets = *(caloJets.product());
   const reco::VertexCollection &   sv	 = *(secondaryVertices.product());
   const reco::VertexCollection &   pv	 = *(primaryVertices.product());

   //use the first primary vertex for now
   const reco::Vertex primaryVertex = pv[0];
   
   // loop over each jet
   reco::CaloJetCollection::const_iterator jj = jets.begin();
   for(; jj != jets.end() ; ++jj){

     edm::LorentzVector p4 = jj->detectorP4();
     float jet_pt = p4.pt();
     float jet_eta = p4.eta();
     float jet_phi = jj.phi();

     //start building the vertex scores
     edm::AssociationMap<edm::OneToValue<reco::Vetex, float>> vertexScores; 
     float bestScore = 9999;
     reco::Vertex bestVertex;

     //loop over each of the vertices to be scored
     reco::VertexCollection::const_terator ss = sv.begin();
     for(; ss != sv.end() ; ++ss) {

       float score = 0;             
       math::XYZTLorentzVectorD kin_sum= ss->p4(0, .5);
             
       //loop over each track in the vertex       	
       reco::Vertex::const_iterator tt = ss.tracks_begin();
       for(; tt != ss.tracks_end(); ++tt) {

	 float track_outerEta = tt->outerEta();
	 float track_outerPhi = tt->outerPhi();

	 float dR = reco:deltaR(jet_eta, jet_phi, track_outerEta, track_outerPhi);

	 if (dR > 1.0) continue;
	 else score += 1.0 / dR;	 	  	 
       }

       vertexScores->insert(*ss, score);   
     }

     // put together the scores with the best vertex 
     JetVertexAssociation association(vertexScores, *jj, bestVertex, algoName);     
     jetVertexAssociationCollection.push_back(association);         
   }     

   iEvent.put(jetVertexAssociationCollection, outputLabel_)
}

// ------------ method called once each job just before starting event loop  ------------
void 
DisplacedJetSVAssociator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DisplacedJetSVAssociator::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
DisplacedJetSVAssociator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DisplacedJetSVAssociator::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DisplacedJetSVAssociator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DisplacedJetSVAssociator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisplacedJetSVAssociator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedJetSVAssociator);
