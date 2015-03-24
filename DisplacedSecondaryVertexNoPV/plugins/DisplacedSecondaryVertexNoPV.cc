// -*- C++ -*-
//
// Package:    DisplacedJets/DisplacedSecondaryVerxtexNoPV
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

//#define DEBUG 

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
  typedef typename IPTI::input_container::value_type input_item;


  explicit DisplacedSecondaryVertexNoPV(const edm::ParameterSet&);
  ~DisplacedSecondaryVertexNoPV();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void produce(edm::Event&, const edm::EventSetup&) override;

  //tokens to get
  edm::EDGetTokenT<reco::BeamSpot> token_BeamSpot;
  edm::EDGetTokenT<std::vector<IPTI> > token_trackIPTagInfo;

  //track sorting
  reco::btag::SortCriteria	sortCriterium;

  //vertex reconstruction parameters
  edm::ParameterSet		vtxRecoPSet;

  //sv vertex ordering criterium 
  VertexSorting<SecondaryVertex>	vertexSorting; 
 
  //method for marking tracks used in the vertex fit
  void markUsedTracks(TrackDataVector & trackData, const input_container & trackRefs, const SecondaryVertex & sv,size_t idx);

  //one parameter function which takes a transient vertex and returns a secondary vertex
  struct SVBuilder: public std::unary_function<const VTX&, SecondaryVertex> {
    
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
DisplacedSecondaryVertexNoPV<IPTI,VTX>::DisplacedSecondaryVertexNoPV(const edm::ParameterSet& iConfig):
  sortCriterium(TrackSorting::getCriterium(iConfig.getParameter<std::string>("trackSort"))),
  vtxRecoPSet(iConfig.getParameter<edm::ParameterSet>("vertexReco")),
  vertexSorting(iConfig.getParameter<edm::ParameterSet>("vertexSelection"))
{

  token_trackIPTagInfo =  consumes<std::vector<IPTI> >(iConfig.getParameter<edm::InputTag>("trackIPTagInfos"));
  token_BeamSpot = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"));

  produces<SVTagInfoCollection>();
}

template <class IPTI,class VTX>
DisplacedSecondaryVertexNoPV<IPTI,VTX>::~DisplacedSecondaryVertexNoPV()
{

}

//
// member functions
//
template <class IPTI,class VTX>
void DisplacedSecondaryVertexNoPV<IPTI,VTX>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

#ifdef DEBUG
   std::cout << "[DSVNoPV] Beginning Production....." << std::endl;
#endif


   edm::ESHandle<TransientTrackBuilder> trackBuilder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

   // tag infos from jets
   edm::Handle<std::vector<IPTI> > trackIPTagInfos;
   iEvent.getByToken(token_trackIPTagInfo, trackIPTagInfos);

   std::auto_ptr<ConfigurableVertexReconstructor> vertexReco;
   vertexReco.reset( new ConfigurableVertexReconstructor(vtxRecoPSet));


   edm::Handle<BeamSpot> beamSpot;
   iEvent.getByToken(token_BeamSpot, beamSpot);

   // product to put into the event 
   std::auto_ptr<SVTagInfoCollection> tagInfos(new SVTagInfoCollection);

#ifdef DEBUG
   std::cout << "[DSVNoPV] Looping Tag Info..." << std::endl;
#endif


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
     //const std::vector<reco::btag::TrackIPData> &ipData = iterJets->impactParameterData();
     std::vector<std::size_t> indices = iterJets->sortedIndexes(sortCriterium);

     input_container trackRefs = iterJets->sortedTracks(indices);

     //transient track containers 
     std::vector<TransientTrack> fitTracks;
     // TransientTrackMap primariesMap;
     // TransientTrackMap::const_iterator pos = primariesMap.find(reco::btag::toTrack((trackRef)));

     reco::TrackRefVector jetTracks = iterJets->selectedTracks();

#ifdef DEBUG
   std::cout << "[DSVNoPV] Looping Tracks in Tag Info..." << std::endl;
#endif

     for(unsigned int tt = 0; tt < indices.size(); tt++) {
       typedef typename TemplatedSecondaryVertexTagInfo<IPTI,VTX>::IndexedTrackData IndexedTrackData;
       
       //const input_item &trackRef = trackRefs[tt]; //normally used for the track selection
       trackData.push_back(IndexedTrackData());
       trackData.back().first = indices[tt];

       // every matched track is currently used for the vertex fit 
       trackData.back().second.svStatus = TemplatedSecondaryVertexTagInfo<IPTI,VTX>::TrackData::trackUsedForVertexFit;
     }

#ifdef DEBUG
   std::cout << "[DSVNoPV] Add tracks together for Fit..." << std::endl;
#endif

     // loop over the tracks to add them to be fit
     for(unsigned int ii = 0; ii < jetTracks.size(); ii++) {
       TransientTrack fitTrack = trackBuilder->build(jetTracks[ii]);
       fitTracks.push_back(fitTrack);
     }

#ifdef DEBUG
   std::cout << "[DSVNoPV] Performing Fit ..." << std::endl;
#endif
     
     // reconstruct the vertex using the fitter
     std::vector<TransientVertex> fittedSVs = vertexReco->vertices(fitTracks);
     std::vector<SecondaryVertex> SVs;

#ifdef DEBUG
   std::cout << "[DSVNoPV] Push Back SVs from fit.." << std::endl;
#endif


     SVBuilder svBuilder(pv, jetDir, false, 0); // pv, jetdir, withPVError, minTrackWeight
     // converted the transient vertex into a secondary vertex object (keep everything)
     for(int vv = 0; vv < fittedSVs.size(); vv++ ) {
       SVs.push_back(svBuilder(fittedSVs[vv]));
     }
     
#ifdef DEBUG
   std::cout << "[DSVNoPV] Extracting SV Info ..." << std::endl;
#endif

     //collect the sv data
     std::vector<unsigned int> vtxIndices = vertexSorting(SVs);
     std::vector<typename TemplatedSecondaryVertexTagInfo<IPTI, VTX>::VertexData> svData;
     svData.resize(vtxIndices.size());
     for(unsigned int vv = 0; vv < vtxIndices.size(); vv++) {
       const SecondaryVertex &sv = SVs[vtxIndices[vv]];
       
       //assign info from fit to the struct object
       svData[vv].vertex = sv;
       svData[vv].dist2d = sv.dist2d();
       svData[vv].dist3d = sv.dist3d();
       svData[vv].direction = flightDirection(pv, sv);

       // mark tracks successfully used in vertex fit
       //markUsedTracks(trackData, trackRefs, sv, vv);
     }

#ifdef DEBUG
   std::cout << "[DSVNoPV] Push Back Tag Infos ..." << std::endl;
#endif

     // final tagInfos push back
     tagInfos->push_back(TemplatedSecondaryVertexTagInfo<IPTI, VTX>( trackData, svData, SVs.size(), edm::Ref<std::vector<IPTI> >(trackIPTagInfos, iterJets - trackIPTagInfos->begin())));
   }
   
   iEvent.put(tagInfos);   
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
template <class IPTI, class VTX>
void DisplacedSecondaryVertexNoPV<IPTI,VTX>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  // desc.add<uint32_t>("vertexReco", 1);
  // desc.add<edm::InputTag>("trackIPTagInfos",edm::InputTag(" "," "));
  // desc.add<edm::InputTag>("beamSpotTag",edm::InputTag(" "," "));
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//Need specialized template because reco::Vertex iterators are TrackBase and it is a mess to make general
template <>
void  DisplacedSecondaryVertexNoPV<TrackIPTagInfo, reco::Vertex>::markUsedTracks(TrackDataVector & trackData, const input_container & trackRefs, const SecondaryVertex & sv, size_t idx)
{
	for(Vertex::trackRef_iterator iter = sv.tracks_begin();	iter != sv.tracks_end(); ++iter) {
		// if (sv.trackWeight(*iter) < minTrackWeight)
		// 	continue;

		typename input_container::const_iterator pos =
			std::find(trackRefs.begin(), trackRefs.end(),
					iter->castTo<input_item>());

		if (pos == trackRefs.end() ) {		  
		  throw cms::Exception("TrackNotFound")
		    << "Could not find track from secondary vertex in original tracks " << std::endl;
		} else {
			unsigned int index = pos - trackRefs.begin();
			trackData[index].second.svStatus =
				(btag::TrackData::Status)
				((unsigned int)btag::TrackData::trackAssociatedToVertex + idx);
		}
	}
}

// template <>
// void  DisplacedSecondaryVertexNoPV<CandIPTagInfo, reco::VertexCompositePtrCandidate>::markUsedTracks(TrackDataVector & trackData, const input_container & trackRefs, const SecondaryVertex & sv,size_t idx)
// {
// 	for(typename input_container::const_iterator iter = sv.daughterPtrVector().begin(); iter != sv.daughterPtrVector().end(); ++iter)
// 	{
// 		typename input_container::const_iterator pos =
// 			std::find(trackRefs.begin(), trackRefs.end(), *iter);

// 		if (pos != trackRefs.end() )
// 		{
// 			unsigned int index = pos - trackRefs.begin();
// 			trackData[index].second.svStatus =
// 				(btag::TrackData::Status)
// 				((unsigned int)btag::TrackData::trackAssociatedToVertex + idx);
// 		}
// 	}
// }

template <>
typename DisplacedSecondaryVertexNoPV<TrackIPTagInfo,reco::Vertex>::SecondaryVertex  
DisplacedSecondaryVertexNoPV<TrackIPTagInfo,reco::Vertex>::SVBuilder::operator () (const TransientVertex &sv) const {
	if(sv.originalTracks().size() > 0 && sv.originalTracks()[0].trackBaseRef().isNonnull())
		return SecondaryVertex(pv, sv, direction, withPVError);
	else {
		edm::LogError("UnexpectedInputs") << "Building from Candidates, should not happen!";
		return SecondaryVertex(pv, sv, direction, withPVError);
	}
}

template <>
typename DisplacedSecondaryVertexNoPV<CandIPTagInfo, reco::VertexCompositePtrCandidate>::SecondaryVertex
DisplacedSecondaryVertexNoPV<CandIPTagInfo, reco::VertexCompositePtrCandidate>::SVBuilder::operator () (const TransientVertex &sv) const
{
	if(sv.originalTracks().size()>0 && sv.originalTracks()[0].trackBaseRef().isNonnull())
	{
		edm::LogError("UnexpectedInputs") << "Building from Tracks, should not happen!";
		VertexCompositePtrCandidate vtxCompPtrCand;
		
		vtxCompPtrCand.setCovariance(sv.vertexState().error().matrix_new());
		vtxCompPtrCand.setChi2AndNdof(sv.totalChiSquared(), sv.degreesOfFreedom());
		vtxCompPtrCand.setVertex(Candidate::Point(sv.position().x(),sv.position().y(),sv.position().z()));
		
		return SecondaryVertex(pv, vtxCompPtrCand, direction, withPVError);
	}
	else
	{
		VertexCompositePtrCandidate vtxCompPtrCand;
		
		vtxCompPtrCand.setCovariance(sv.vertexState().error().matrix_new());
		vtxCompPtrCand.setChi2AndNdof(sv.totalChiSquared(), sv.degreesOfFreedom());
		vtxCompPtrCand.setVertex(Candidate::Point(sv.position().x(),sv.position().y(),sv.position().z()));
		
		Candidate::LorentzVector p4;
		for(std::vector<reco::TransientTrack>::const_iterator tt = sv.originalTracks().begin(); tt != sv.originalTracks().end(); ++tt)
		{
			// if (sv.trackWeight(*tt) < minTrackWeight)
			// 	continue;
			
			const CandidatePtrTransientTrack* cptt = dynamic_cast<const CandidatePtrTransientTrack*>(tt->basicTransientTrack());
			if ( cptt == 0 )
				edm::LogError("DynamicCastingFailed") << "Casting of TransientTrack to CandidatePtrTransientTrack failed!";
			else {
				p4 += cptt->candidate()->p4();
				vtxCompPtrCand.addDaughter(cptt->candidate());
			}
		}

		vtxCompPtrCand.setP4(p4);
		
		return SecondaryVertex(pv, vtxCompPtrCand, direction, withPVError);
	}
}


typedef DisplacedSecondaryVertexNoPV<TrackIPTagInfo, reco::Vertex> DisplacedSecondaryVertexProducer;
//typedef DisplacedSecondaryVertexNoPV<CandIPTagInfo, reco::VertexCompositePtrCandidate> DisplacedCandSecondaryVertexProducer;

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedSecondaryVertexProducer);
//DEFINE_FWK_MODULE(DisplacedCandSecondaryVertexProducer);
