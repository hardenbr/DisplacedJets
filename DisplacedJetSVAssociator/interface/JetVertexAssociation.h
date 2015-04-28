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
  JetVertexAssociation(std::string nameString) { name = nameString;}

  // fillers
  void addCaloJet(const reco::CaloJet &);
  void addVertex(const reco::Vertex &);
  void setPrimaryVertex(const reco::Vertex & pv){ primaryVertex = pv; };
  
  // accessors
  int	getNVertices(){ return vertexCollection.size(); }
  int	getNJets() { return caloJetCollection.size(); } 
  std::string getName() { return name; }

  // score
  reco::Vertex	getBestVertex(const reco::CaloJet&, const std::string&);
  float		getVertexJetScore(const reco::CaloJet&, const reco::Vertex & , const std::string &);    
  
 private:
  reco::Vertex		    primaryVertex;
  reco::CaloJetCollection   caloJetCollection;
  reco::VertexCollection    vertexCollection;    
  std::string		    name;
};

void JetVertexAssociation::addCaloJet(const reco::CaloJet & jet) { caloJetCollection.push_back(jet); }

void JetVertexAssociation::addVertex(const reco::Vertex & vertex) { vertexCollection.push_back(vertex); }

float JetVertexAssociation::getVertexJetScore(const reco::CaloJet & jet, const reco::Vertex & vertex, const std::string & algo) {
  math::XYZTLorentzVectorD  p4	    = jet.detectorP4();
  //float			    jet_pt  = p4.pt();
  float			    jet_eta = p4.eta();
  float			    jet_phi = p4.phi();
  float			    score   = 0;

  std::vector<reco::TrackBaseRef>::const_iterator tt = vertex.tracks_begin();
  for(; tt != vertex.tracks_end(); ++tt) {
    float   track_outerEta = (*tt)->outerEta();
    float   track_outerPhi = (*tt)->outerPhi();    
    float   dR		   = reco::deltaR(jet_eta, jet_phi, track_outerEta, track_outerPhi);
      
    if (dR > 1.0) 
      continue;
    else 
      score += 1.0 / dR;	 	  	 
  }
  
  return score;
}

reco::Vertex JetVertexAssociation::getBestVertex(const reco::CaloJet & jet, const std::string & algo) {
  // keep track of the best vertex
  float	bestScore = -1;
  reco::Vertex	bestVertex;

  reco::VertexCollection::const_iterator ss = vertexCollection.begin();
  for(; ss != vertexCollection.end(); ++ss) {
    float score = getVertexJetScore(jet, *ss, algo);
    
    if (score > bestScore) {
      bestScore	 = score;
      bestVertex = *ss;
    }    
  }

  if (bestScore > 0) return bestVertex;  // we find at least one vertex with a track within dR = 1.0
  
  return primaryVertex;    
}





  
