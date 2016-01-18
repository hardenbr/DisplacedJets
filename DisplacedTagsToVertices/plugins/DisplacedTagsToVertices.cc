// -*- C++ -*-
//
// Package:    DisplacedJets/DisplacedTagsToVertices
// Class:      DisplacedTagsToVertices
// 
/**\class DisplacedTagsToVertices DisplacedTagsToVertices.cc DisplacedJets/DisplacedTagsToVertices/plugins/DisplacedTagsToVertices.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Robert Hardenbrook
//         Created:  Tue, 16 Dec 2014 13:56:03 GMT
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

#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPData.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "DataFormats/BTauReco/interface/VertexTypes.h"


#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


//
// class declaration
//

class DisplacedTagsToVertices : public edm::EDProducer {
public:
  explicit DisplacedTagsToVertices(const edm::ParameterSet&);
  ~DisplacedTagsToVertices();
  
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  

  edm::InputTag tag_secondaryVertexTagInfo_;   
  edm::InputTag tag_genParticles_;   
  std::string outputLabel_;
  float jetPtCut_;
  int debug_;
  bool isSignalMC_;
  bool doGenMatch_;

      
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
DisplacedTagsToVertices::DisplacedTagsToVertices(const edm::ParameterSet& iConfig)
{

  tag_secondaryVertexTagInfo_ = iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertexTagInfo");
  tag_genParticles_	      = iConfig.getUntrackedParameter<edm::InputTag>("genParticleTag");
  jetPtCut_		      = iConfig.getUntrackedParameter<double>("jetPtCut");
  outputLabel_		      = iConfig.getUntrackedParameter<std::string>("outputLabel");
  isSignalMC_		      = iConfig.getUntrackedParameter<bool>("isSignalMC");
  doGenMatch_		      = iConfig.getUntrackedParameter<bool>("doGenMatch");
  debug_		      = iConfig.getUntrackedParameter<int>("debug");

  consumes<reco::SecondaryVertexTagInfoCollection>(tag_secondaryVertexTagInfo_);
  consumes<reco::GenParticleCollection>(tag_genParticles_);

  produces<reco::VertexCollection>(outputLabel_);  
}


DisplacedTagsToVertices::~DisplacedTagsToVertices()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DisplacedTagsToVertices::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::SecondaryVertexTagInfoCollection> secondaryVertexTagInfo;
   iEvent.getByLabel(tag_secondaryVertexTagInfo_, secondaryVertexTagInfo);

   std::auto_ptr<reco::VertexCollection> displacedSvVertexCollectionResult(new reco::VertexCollection);                                                                                          
   reco::VertexCollection displacedSvVertexCollection;       

   const reco::SecondaryVertexTagInfoCollection & sv = *(secondaryVertexTagInfo.product());
   reco::SecondaryVertexTagInfoCollection::const_iterator svinfo = sv.begin();
   // loop over each jet 
   for(; svinfo != sv.end(); ++svinfo){

     //apply a jet threshold
     if ( svinfo->jet()->pt()  < jetPtCut_){
      continue;
     }

     edm::Handle<reco::GenParticleCollection> genParticles;

     bool keepSV = false;
     // generator matching to truely displaced jets
     if (isSignalMC_ && doGenMatch_) {
       iEvent.getByLabel("genParticles", genParticles);    

       // get the jet attributes 
       float calopt = svinfo->jet()->pt();
       float caloeta = svinfo->jet()->eta();
       float calophi = svinfo->jet()->phi();

       if(debug_ > 0 && calopt > 40.0) { 
	   std::cout << "[SV DEBUG JET]"  << " pt " << calopt << " eta " << caloeta  <<  " phi " << calophi << std::endl;      	 
       }

       // check all gen particles for a match to the jet
       for(size_t pp = 0; pp < genParticles->size(); ++pp) {
	 const reco::GenParticle & part = (*genParticles)[pp];
	 int id = part.pdgId();
	 int st = part.status();  
	 
	 if (st != 3 || fabs(id) > 6 ) continue;     
	 
	 // const reco::Candidate * mom = part.mother();
	 double genpt = part.pt(), geneta = part.eta(), genphi = part.phi();
	 float dr = reco::deltaR( geneta, genphi, caloeta, calophi);
	 float dpt = fabs(calopt - genpt) / genpt;
	 
	 // found a match
	 if (dr < .5 && dpt < .2) {
	   std::cout << "[SV GEN MATCHED JET] id " << id << " status " << st << " pt " << genpt << " eta " << geneta  <<  " phi " << genphi << std::endl;      
	   if (debug_ > 0 ) std::cout << "[SV GEN MATCHED DEBUG] NUMBER OF SV: " << svinfo->nVertices() << std::endl;
	   keepSV = true;
	   break;
	 }
       } // end loop over gen particles
     } // end if statement for gen match

    int nSV = svinfo->nVertices();    

    // loop over each SVvertex that jet has
    for(int vv = 0; vv < nSV; vv++) {
     reco::Vertex svVertex = svinfo->secondaryVertex(vv);       
     if (keepSV || !doGenMatch_){
       if (debug_ > 0) std::cout << "[SV DEBUG] STORING SV" << std::endl;
       displacedSvVertexCollection.push_back(svVertex);                                                                    }
    }
   }

   *displacedSvVertexCollectionResult = displacedSvVertexCollection;
   iEvent.put(displacedSvVertexCollectionResult, outputLabel_);

/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::unique_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(std::move(pOut));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
DisplacedTagsToVertices::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DisplacedTagsToVertices::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
DisplacedTagsToVertices::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DisplacedTagsToVertices::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DisplacedTagsToVertices::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DisplacedTagsToVertices::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisplacedTagsToVertices::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedTagsToVertices);
