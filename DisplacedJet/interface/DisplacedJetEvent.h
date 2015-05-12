
typedef std::vector<DisplacedJet>  DisplacedJetCollection;

class DisplacedJetEvent {
 public:

  DisplacedJetEvent(const bool&, const reco::CaloJetCollection&, const reco::VertexCollection&, const float&, const float&, const int&);

  // accessor
  int getNJets() { return djets.size(); } 
  
  // jet associated info mergers
  void mergeCaloIPTagInfo(const reco::TrackIPTagInfoCollection&) ;
  void mergeSVTagInfo(const reco::SecondaryVertexTagInfoCollection&);
  
  // associated collections
  void addIVFVertices(const reco::VertexCollection & vertices);
  
  // matching
  void doGenMatching(const reco::GenParticleCollection & genParticleCollection, 
		     const bool & doCaloJetMatch, const bool & doGenVtxMatch, const bool & doGenVtxID,
		     const float& ptMatch, const float& drMatch, const float& vtxMatchThreshold); 
  void doSimMatching(const edm::SimVertexContainer & simVertexCollection);
  
  // helper method
  DisplacedJet & findDisplacedJetByPtEtaPhi(const float& pt, const float& eta, const float& phi);
  DisplacedJetCollection & getDisplacedJets() { return djets; }
  
  float minPT;
  float minEta; 
  
 private:
  std::vector<DisplacedJet> djets;
  int jetIDCounter = 0;     
  int debug;  
};

DisplacedJetEvent::DisplacedJetEvent(const bool& isMC, const reco::CaloJetCollection & caloJets, const reco::VertexCollection & primaryVertices, const float& minPT_, const float& minEta_, const int & debug_) {

  minPT	 = minPT_;
  minEta = minEta_;
  debug	 = debug_;

  const reco::Vertex & firstPV = *primaryVertices.begin();

  if (debug > 1) std::cout << "[DEBUG] Constrcuted Event From Calo jets " << std::endl;
  reco::CaloJetCollection::const_iterator jetIter = caloJets.begin();
  for(; jetIter != caloJets.end(); ++jetIter) {    
    float pt = jetIter->pt(),  eta = jetIter->eta();
    if (pt < minPT || fabs(eta) > minEta) continue;
    
    DisplacedJet djet(*jetIter, firstPV, isMC, jetIDCounter, debug);
    
    djets.push_back(djet);
    jetIDCounter++;
  }  
}

void DisplacedJetEvent::addIVFVertices(const reco::VertexCollection & vertices) {
  if (debug > 1) std::cout << "[DEBUG] Adding IVF Vertices " << std::endl;
  std::vector<DisplacedJet>::iterator djetIter = djets.begin();
  for(; djetIter != djets.end(); ++djetIter) {
    djetIter->addIVFCollection(vertices);        
  }  
}

void DisplacedJetEvent::mergeSVTagInfo(const reco::SecondaryVertexTagInfoCollection& svTagInfoCollection) {
  if (debug > 1) std::cout << "[DEBUG] Merging SV Tag Info " << std::endl;
 
  reco::SecondaryVertexTagInfoCollection::const_iterator svinfo = svTagInfoCollection.begin();    
  int jj = 0;
  for(; svinfo != svTagInfoCollection.end(); ++svinfo, jj++){
    const reco::Jet * jet = svinfo->jet().get();    
    float pt = jet->pt(), eta = jet->eta(), phi = jet->phi();
    if (pt < minPT || fabs(eta) > minEta) continue;    

    // find the jet and corresponding track refrences
    DisplacedJet & djet = findDisplacedJetByPtEtaPhi(pt, eta, phi);
    
    djet.addSVTagInfo(*svinfo);
  }
}

// fills each jet with ip tag info variables 
// fills each jet track collection with tracks used for ip variables
void DisplacedJetEvent::mergeCaloIPTagInfo(const reco::TrackIPTagInfoCollection & ipTagInfo) {
  if (debug > 1) std::cout << "[DEBUG 1] Merging CALO IP Tag Info " << std::endl;

  reco::TrackIPTagInfoCollection::const_iterator ipInfoIter = ipTagInfo.begin(); 
  for(; ipInfoIter != ipTagInfo.end(); ++ipInfoIter) {
    
    // check the jet against the acceptance
    const reco::Jet jet = *ipInfoIter->jet();
    const float & pt = jet.pt(), eta = jet.eta(), phi = jet.phi();
    if (pt < minPT || fabs(eta) > minEta) continue;    

    // find the jet and corresponding track refrences
    DisplacedJet & djet = findDisplacedJetByPtEtaPhi(pt, eta, phi);
    const reco::TrackRefVector trackRefs = ipInfoIter->selectedTracks();        

    // add the track and ip info 
    djet.addIPTagInfo(*ipInfoIter);    
    djet.addCaloTrackInfo(trackRefs);
  }
}

DisplacedJet & DisplacedJetEvent::findDisplacedJetByPtEtaPhi(const float& pt, const float& eta, const float& phi) {
  if (debug > 2) std::cout << "[DEBUG] Finding Displaced Jet By PT ETA PHI " << std::endl;
  std::vector<DisplacedJet>::iterator djetIter = djets.begin();
  bool found = false;
  for(; djetIter != djets.end(); ++djetIter) {
    float djet_pt = djetIter->caloPt, djet_eta = djetIter->caloEta, djet_phi = djetIter->caloPhi;
    
    bool    pt_match  = djet_pt == pt;
    bool    eta_match = djet_eta == eta;
    bool    phi_match = djet_phi == phi;
    
    found = pt_match && eta_match && phi_match;

    if(found) {
      found = true;
      break;
    }          
  }

  assert(found);
  return *djetIter;
}

void DisplacedJetEvent::doGenMatching( const reco::GenParticleCollection& genParticles, 
				       const bool& doCaloJetMatch = true, const bool& doGenVtxMatch = true, const bool& doGenVtxID = true,
				       const float& ptMatch = 0.2, const float& dRMatch = 0.7,
				       const float& vtxMatchThreshold = 0.05) {

  // do the particle matching to calo jets
  // also do the vertex matching to gen particle (x,y,z)
  std::vector<DisplacedJet>::iterator djetIter = djets.begin();
  for(; djetIter != djets.end(); ++djetIter) {
    bool isCaloGenMatched = false, isIVFGenVertexMatched = false, isSVGenVertexMatched = false;

    if(doCaloJetMatch) {
       isCaloGenMatched = djetIter->doGenCaloJetMatching(ptMatch, dRMatch, genParticles);
    }
    
    if(doGenVtxMatch) {
      djetIter->doGenVertexJetMatching(vtxMatchThreshold, genParticles);
      isIVFGenVertexMatched = djetIter->ivfIsGenMatched;
      isSVGenVertexMatched  = djetIter->svIsGenMatched;      
    }    
    if(doGenVtxID) {
      djetIter->doGenVertexID(vtxMatchThreshold, genParticles);
    }
    if(debug > 1) std::cout <<  "[GEN MATCH] isCalomatch: " << isCaloGenMatched << " isIVFGenVertexMatched "
			  << isIVFGenVertexMatched << " isSVGEnVertexMatched " << isSVGenVertexMatched << std::endl;

  }  // end loop over displaced jets
}
