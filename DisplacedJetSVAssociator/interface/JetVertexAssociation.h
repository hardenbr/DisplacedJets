#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Common/interface/AssociationMap.h"

// Used to match vertices to a jet based on some metric
// The metric is incoded in an association map for each jet which designates the 
// degree to which the vertex is matched to the jet
// The vertex with the best score is stored as bestVtx
class JetVertexAssocation {
 public:
  typedef edm::AssociationMap< edm::OneToValue<reco::VertexCollection, float>> vertexFloatMap;
  
  JetVertexAssocation(const edm::AssociationMap< edm::OneToValue<reco::VertexCollection, float>> &  vertexScores_,
		      const reco::CaloJet &							    jet_,
		      const reco::Vertex &							    bestVtx_,
		      const std::string &                                                           algo ){    
    
    jet		 = jet_;
    bestVtx	 = bestVtx_;
    vertexScores = vertexScores_;
    algoName	 = algoName_;      
  }
  
  JetVertexAssocation() {}
  
  virtual ~JetVertexAssocation() {}

  virtual JetVertexAssocation * clone(void) const {
    return new JetVertexAssociation(*this); 
  }

  const reco::Vertex & getBestVertex() const { return bestVtx; }

  const reco::CaloJet & getJet() const { return jet; }

  const std::string & getAssociationAlgo() const { return algoName; }

  const float getVertexScore(reco::Vertex & vertex) {
    vertexFloatMap::const_iterator finder = vertexScores.find(vertex);
    assert(finder != vertexScore.end());   

    return finder.val();
  }

 private:
  edm::AssociationMap< edm::OneToValue<reco::VertexCollection, float>>  vertexSecores_; 
  reco::CaloJet jet;
  reco::Vertex	bestVtx;     
};
  
