// -*- C++ -*-
//
// Package:    DisplacedJets/DisplacedSecondaryVertexNoPV
// Class:      DisplacedSecondaryVertexNoPV
// 
/**\class DisplacedSecondaryVertexNoPV DisplacedSecondaryVertexNoPV.cc DisplacedJets/DisplacedSecondaryVertexNoPV/plugins/DisplacedSecondaryVertexNoPV.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Robert Hardenbrook
//         Created:  Thu, 19 Mar 2015 08:43:57 GMT
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

// Taken from RecoBtag/SecondaryVertex/TemplatedSecondaryVertexProducer/ 
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "RecoBTag/SecondaryVertex/interface/TemplatedSecondaryVertex.h"

//interface
#include "RecoBTag/SecondaryVertex/interface/TrackSelector.h"
#include "RecoBTag/SecondaryVertex/interface/TrackSorting.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoBTag/SecondaryVertex/interface/VertexFilter.h"
#include "RecoBTag/SecondaryVertex/interface/VertexSorting.h"

// transient track related 
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/CandidatePtrTransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

using namespace reco;

GlobalVector flightDirection(const reco::Vertex & pv, const reco::Vertex & sv) {
return  GlobalVector(sv.x() - pv.x(), sv.y() - pv.y(),sv.z() - pv.z());
}
GlobalVector flightDirection(const reco::Vertex & pv, const reco::VertexCompositePtrCandidate & sv) {
return  GlobalVector(sv.vertex().x() - pv.x(), sv.vertex().y() - pv.y(),sv.vertex().z() - pv.z());
}

//
// class declaration
//

template <class IPTI,class VTX>
class DisplacedSecondaryVertexNoPV : public edm::EDProducer {
public:
  //type declarations
  typedef std::vector<TemplatedSecondaryVertexTagInfo<IPTI,VTX> > SVTagInfoCollection;
  typedef TemplatedSecondaryVertex<VTX> SecondaryVertex;
  typedef typename std::vector<reco::btag::IndexedTrackData> TrackDataVector;
  typedef typename IPTI::input_container input_container;

  explicit DisplacedSecondaryVertexNoPV(const edm::ParameterSet&);
  ~DisplacedSecondaryVertexNoPV();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //tokens to get
  edm::EDGetTokenT<reco::BeamSpot> token_BeamSpot;
  edm::EDGetTokenT<std::vector<IPTI> > token_trackIPTagInfo;

  //track sorting
  reco::btag::SortCriteria	sortCriterium;

  //vertex reconstruction parameters
  edm::ParameterSet		vtxRecoPSet;

  //sv vertex ordering criterium 
  VertexSorting<SecondaryVertex>	vertexSorting; 
 
  //one parameter function which takes a transient vertex and returns a secondary vertex
  struct SVBuilder :
    public std::unary_function<const VTX&, SecondaryVertex> {
    
    SVBuilder(const reco::Vertex &pv,
	      const GlobalVector &direction,
		          bool withPVError,
	      double minTrackWeight) :
      pv(pv), direction(direction),
      withPVError(withPVError),
      minTrackWeight(minTrackWeight) {}
    SecondaryVertex operator () (const TransientVertex &sv) const;
    
    SecondaryVertex operator () (const VTX &sv) const
    { return SecondaryVertex(pv, sv, direction, withPVError); }
        
    const Vertex	&pv;
    const GlobalVector	&direction;
    bool		withPVError; // not currently used
    double 		minTrackWeight; // not currently used
  };    

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
template <class IPTI,class VTX>
DisplacedSecondaryVertexNoPV::DisplacedSecondaryVertexNoPV(const edm::ParameterSet& iConfig):
  vtxRecoPSet(iConfig.getParameter<edm::ParameterSet>("vertexReco")),
  vertexSorting(iConfig.getParameter<edm::ParameterSet>("vertexSelection"))
{

  produces<SVTagInfoCollection>
}

template <class IPTI,class VTX>
DisplacedSecondaryVertexNoPV::~DisplacedSecondaryVertexNoPV()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
template <class IPTI,class VTX>w
void DisplacedSecondaryVertexNoPV::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::ESHandle<TransientTrackBuilder> trackBuilder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",
	                                   trackBuilder);
   // tag infos from jets
   edm::Handle<std::vector<IPTI>> trackIPTagInfos;
   iEvent.getByToken(token_trackIPTagInfo, trackIPTagInfos);

   std::auto_ptr<ConfigurableVertexReconstructor> vertexReco;

   edm::Handle<BeamSpot> beamSpot;
   event.getByToken(token_BeamSpot,beamSpot);

   // product to put into the event 
   std::auto_ptr<SVTagInfoCollection> tagInfos(new SVTagInfoCollection);

   //loop over each jet contained in the tag infos
   for(typename std::vector<IPTI>::const_iterator iterJets = trackIPTagInfos->begin(); 
       iterJets != trackIPTagInfos->end();
       ++iterJets) {
     
     //container for the track information
     TrackDataVector trackData;

     // PV associated with the jet
     const Vertex &pv = *iterJets->primaryVertex();

     // get the nominal jet direction 
     edm::RefToBase<Jet> jetRef = iterJets->jet();
     GlobalVector jetDir(jetRef->momentum().x(),
			 jetRef->momentum().y(),
			 jetRef->momentum().z());

     // extract ip information
     const std::vector<reco::btag::TrackIPData> &ipData = iterJets->impactParameterData();
     input_container trackRefs = iterJets->sortedTracks(indices);

     //transient track containers 
     std::vector<TransientTrack> fitTracks;
     // TransientTrackMap primariesMap;
     // TransientTrackMap::const_iterator pos = primariesMap.find(reco::btag::toTrack((trackRef)));

     reco::TrackRefVector jetTracks = iterJets->selectedTracks();

     for(unsigned int tt = 0; tt < indices.size(); tt++) {
       typedef typename TemplatedSecondaryVertexTagInfo<IPTI,VTX>::IndexedTrackData IndexedTrackData;
       
       //const input_item &trackRef = trackRefs[tt]; //normally used for the track selection

       trackData.push_back(IndexedTrackData());
       trackData.back().first = indices[tt];
       // every matched track is currently used for the vertex fit 
       trackData.back().second.svStatus = TemplatedSecondaryVertexTagInfo<IPTI,VTX>::TrackData::trackUsedForVertexFit;
     }

     // loop over the tracks to add them to be fit
     for(unsigned int ii = 0; ii < jetTracks.size(); ii++) {
       TransientTrack fitTrack = trackBuilder->(jetTracks[ii]);
       fitTracks.push_back(fitTrack);
     }

     // reconstruct the vertex using the fitter
     std::vector<TransientVertex> fittedSVs =  vertexReco->vertices(fitTracks);
     std::vector<SecondaryVertex> SVs;

     // converted the transient vertex into a secondary vertex object (keep everything)
     for(int ss = 0; ss < fitedSVs.size(); ss++ ) {
       SVs.push_back(DisplacedSecondaryVertexNoPV::SVBuilder(fittedSVs[ss])
     }
     
     //collect the sv data

     std::vector<unsigned int> vtxIndices = vertexSorting(SVs);
     std::vector<typename TemplatedSecondaryVertexTagInfo<IPTI, VTX>::VertexData> svData;
     svData.resize(vtxIndices.size());
     for(unsigned int vv = 0; vv < vtxIndices.size(); vv++) {
       const SecondaryVertex &sv = SVs[vtxIndices[vv]];
       
       svData[vv].vertex = sv;
       svData[vv].dist2d = sv.dist2d();
       svData[vv].dist3d = sv.dist3d();
       svData[vv].direction = flightDirection(pv, sv);
       // mark tracks successfully used in vertex fit
       markUsedTracks(trackData, trackRefs, sv, vv);
     }

     // final tagInfos push back
     tagInfos->push_back(
			 TemplatedSecondaryVertexTagInfo<IPTI, VTX>(
								   trackData, svData, SVs.size(),
								   edm::Ref<std::vector<IPTI> >(trackIPTagInfos,
												iterJets - trackIPTagInfos->begin())));

   }

   event.put(tagInfos);
} 


// ------------ method called once each job just before starting event loop  ------------
void 
DisplacedSecondaryVertexNoPV::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DisplacedSecondaryVertexNoPV::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
DisplacedSecondaryVertexNoPV::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DisplacedSecondaryVertexNoPV::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DisplacedSecondaryVertexNoPV::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DisplacedSecondaryVertexNoPV::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisplacedSecondaryVertexNoPV::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedSecondaryVertexNoPV);
