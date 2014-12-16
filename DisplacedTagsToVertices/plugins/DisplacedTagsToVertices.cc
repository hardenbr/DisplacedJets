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
  std::string outputLabel_;
  float jetPtCut_;
      
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
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed

  tag_secondaryVertexTagInfo_ = iConfig.getUntrackedParameter<edm::InputTag>("secondaryVertexTagInfo");
  jetPtCut_ = iConfig.getUntrackedParameter<double>("jetPtCut");
  outputLabel_ = iConfig.getUntrackedParameter<std::string>("outputLabel");

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

     // apply a jet threshold
     if (svinfo->jet()->pt()  > jetPtCut_){
       continue;
     }

    int nSV = svinfo->nVertices();     

    // loop over each SVvertex that jet has
    for(int vv = 0; vv < nSV; vv++) {
     reco::Vertex svVertex = svinfo->secondaryVertex(vv);       
     displacedSvVertexCollection.push_back(svVertex);                                                                                                                                            
    }
   }

   *displacedSvVertexCollectionResult = displacedSvVertexCollection;
   iEvent.put(displacedSvVertexCollectionResult, outputLabel_);//, "displacedSVertexCollection");            

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
