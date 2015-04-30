// jet 
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

// tracks
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/VertexTypes.h"

// kinematics
#include "DataFormats/Math/interface/deltaR.h"

class JetVertexAssociation {
 public:
  JetVertexAssociation(std::string nameString) { 
    name	    = nameString;
  }

  // fillers
  void addCaloJet( const reco::CaloJet & jet);
  void addVertex( const reco::Vertex & vertex);
  void setPrimaryVertex(const reco::Vertex & pv){ primaryVertex = pv; };
  
  // accessors
  int	getNVertices(){ return vertexCollection.size(); }
  int	getNJets() { return caloJetCollection.size(); } 
  float getBestVertexScore() { return bestVertexScore; }


  std::string getName() { return name; }
  // score
  const std::pair<const reco::Vertex, const float> getBestVertex(const reco::CaloJet&, const std::string&);
  float		getVertexJetScore(const reco::CaloJet&, const reco::Vertex & , const std::string &);    
  
 private:
  reco::Vertex		    primaryVertex;
  reco::CaloJetCollection   caloJetCollection;
  reco::VertexCollection    vertexCollection;    
  reco::Vertex   	    bestVertex;
  float			    bestVertexScore;  
  std::string		    name;
};

void JetVertexAssociation::addCaloJet(const reco::CaloJet & jet) { caloJetCollection.push_back(jet); }

void JetVertexAssociation::addVertex(const  reco::Vertex & vertex) { vertexCollection.push_back(vertex); }

float JetVertexAssociation::getVertexJetScore(const reco::CaloJet & jet, const reco::Vertex & vertex, const std::string & algo) {

  math::XYZTLorentzVectorD  p4	    = jet.detectorP4();
  float			    jet_eta = p4.eta();
  float			    jet_phi = p4.phi();
  float			    score   = 0;

  std::cout << "[DEBUG] [JVA] iterating over tracks" << std::endl;
  reco::Vertex::trackRef_iterator tt = vertex.tracks_begin();
  for(; tt != vertex.tracks_end(); ++tt) {
  std::cout << "[DEBUG] [JVA] Accessing kinematics" << std::endl;
    // float   track_outerPt  = (*tt)->outerPt();
  float   eta = (*tt)->eta();
  float   phi = (*tt)->phi();    
  std::cout << "[DEBUG] [JVA] DR Calc" << std::endl;
  float   dR		   = reco::deltaR(jet_eta, jet_phi, eta, phi);
  
  if (dR > 1.0) 
      continue;
  else 
    score += 1.0 / dR;	 	  	 
  }
  
  return score;
}

const std::pair<const reco::Vertex, const float> JetVertexAssociation::getBestVertex(const reco::CaloJet & jet, const std::string & algo) {
  // keep track of the best vertex
  float		bestScore = -1;
  reco::Vertex	bestVertex;

  std::cout << "[DEBUG] [JVA] iterating over Vertices" << std::endl;
  reco::VertexCollection::const_iterator ss = vertexCollection.begin();
  for(; ss != vertexCollection.end(); ++ss) {
    float score = getVertexJetScore(jet, *ss, algo);
    
    if (score > bestScore) {
      bestScore	 = score;
      bestVertex = *ss;
    }    
  }

  if (bestScore > 0) {
    const std::pair<const reco::Vertex, const float> vertexPair(bestVertex, bestScore);
    return vertexPair;  // we find at least one vertex with a track within dR = 1.0
  }

  const std::pair<const reco::Vertex, const float> vertexPair(primaryVertex, 0);
  return vertexPair;
}





  
