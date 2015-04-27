#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BTauReco/interface/VertexTypes.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Used to match vertices to a jet based on some metric
// The metric is incoded in an association map for each jet which designates the 
// degree to which the vertex is matched to the jet
// The vertex with the best score is stored as bestVtx
class JetVertexAssociation {
 public:
  //  typedef edm::AssociationMap< edm::OneToValue<reco::VertexRef, float>> vertexFloatMap;
  typedef std::vector<std::pair<reco::Vertex, float>> vertexFloatMap;

  JetVertexAssociation(const vertexFloatMap &	vertexScores_,
		       const reco::CaloJet &	jet_,
		       const reco::Vertex &	bestVtx_,
		       const std::string &      algoName_ ){
    vertexScores = vertexScores_;
    jet	       = jet_;
    bestVtx      = bestVtx_;
    algoName     = algoName_;  
  }
      
  virtual JetVertexAssociation * clone(void) const {
    return new JetVertexAssociation(*this); 
  }

  const reco::Vertex & getBestVertex() const { return bestVtx; }

  const reco::CaloJet & getJet() const { return jet; }

  const std::string & getAssociationAlgo() const { return algoName; }

  const float getVertexScore(reco::Vertex vertex) {
    float score = -9999;

    vertexFloatMap::const_iterator finder = vertexScores.begin();
    for(; finder != vertexScores.end(); ++finder) {
      if (&finder->first == &vertex) {
	score = finder->second; 
	break;
      }
    }
    return score;
  }

 private:
  vertexFloatMap vertexScores; 
  reco::CaloJet jet;
  reco::Vertex	bestVtx;     
  std::string	algoName;
};




  
