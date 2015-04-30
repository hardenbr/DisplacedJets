class DisplacedJetEvent {
 public:
  DisplacedJetEvent(const bool& isMC_, const reco::CaloJetCollection&, const float, const float);
  
  // jet associated info mergers
  void mergeCaloIPTagInfo(const reco::TrackIPTagInfoCollection&) ;
  void mergeSecondaryVertexTagInfo(const reco::SecondaryVertexTagInfoCollection&);
  
  // associated collections
  void addIVFVertices();
  
  // matching
  void doGenMatching(const reco::GenParticleCollection genParticles); 
  void doSimMatching(const edm::SimVertexContainer & simVertexCollection);
  
  DisplacedJet findDisplacedJetByPtEtaPhi(const float& pt, const float& eta, const float& phi);
  
  float minPT;
  float minEta; 
  
 private:
  std::vector<DisplacedJet> djets;
  int jetIDCounter = 0;     
};

DisplacedJetEvent::DisplacedJetEvent(const bool& isMC, const reco::CaloJetCollection & caloJets, const float minPT_, const float minEta_) {

  minPT = minPT_;
  minEta = minEta_;

  reco::CaloJetCollection::const_iterator jetIter = caloJets.begin();
  for(; jetIter != caloJets.end(); ++jetIter) {
    
    float pt = jetIter->pt(),  eta = jetIter->eta();
    if (pt < minPT || fabs(eta) > minEta) continue;
    
    DisplacedJet djet(*jetIter, isMC);
    
    djets.push_back(djet);
    jetIDCounter++;
  }  
}

// fills each jet with ip tag info variables 
// fills each jet track collection with tracks used for ip variables
void 
DisplacedJetEvent::mergeCaloIPTagInfo(const reco::TrackIPTagInfoCollection & ipTagInfo) {
  reco::TrackIPTagInfoCollection::const_iterator ipInfoIter = ipTagInfo.begin(); 
  for(; ipInfoIter != ipTagInfo.end(); ++ipInfoIter) {
    
    // check the jet against the acceptance
    const reco::Jet jet = *ipInfoIter->jet();
    const float & pt = jet.pt(), eta = jet.eta(), phi = jet.phi();
    if (pt < minPT || fabs(eta) > minEta) continue;    

    // find the jet and corresponding track refrences
    DisplacedJet djet = findDisplacedJetByPtEtaPhi(pt, eta, phi);
    const reco::TrackRefVector trackRefs = ipInfoIter->selectedTracks();        

    // add the track and ip info 
    djet.addIPTagInfo(*ipInfoIter);    
    djet.addCaloTrackInfo(trackRefs);
  }
}

DisplacedJet
DisplacedJetEvent::findDisplacedJetByPtEtaPhi(const float& pt, const float& eta, const float& phi) {

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

