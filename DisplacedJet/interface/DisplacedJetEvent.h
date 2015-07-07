
typedef std::vector<DisplacedJet>  DisplacedJetCollection;

class DisplacedJetEvent {
 public:
  DisplacedJetEvent(const bool&, const reco::CaloJetCollection&, const reco::VertexCollection&, const float&, const float&, const int&);

  // event analysis
  bool doesPassPreSelection();  

  void calcNIVFGenMatched(const float & metricThreshold, const reco::GenParticleCollection& genParticles);
  float genMatchMetric(const reco::GenParticle & particle, const reco::Vertex& vertex);


  // accessors
  int	getNJets() { return djets.size(); } 
  std::vector<int>  getNNoVertexTags() { return nNoVertexTagsVector; }
  std::vector<int>  getNShortTags()    { return nShortTagsVector; }
  std::vector<int>  getNMediumTags()   { return nMediumTagsVector; }
  std::vector<int>  getNLongTags()     { return nLongTagsVector; }
  std::vector<int>  getNTotalTags()    { return nTotalTagsVector; }
      
  // tag related
  void doJetTagging(std::vector<float> noVtxThres,
		    std::vector<float> shortThres,
		    std::vector<float> medThres,
		    std::vector<float> longThres, 
		    float shortDistance, float mediumDistance, float longDistance,
		    int dHTWP);
    
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
  
  // kinematic cuts for analysis
  float minPT;
  float minEta; 

  // kinematic variables
  float caloHT;
  std::vector<float> caloDHT;
  float caloMET;

  // ivf related 
  int nIVFReconstructed; 
  int nIVFGenMatched; 
  // n tag vectors indexed by working point threshold
  std::vector<int> nNoVertexTagsVector;
  std::vector<int> nShortTagsVector;
  std::vector<int> nMediumTagsVector;
  std::vector<int> nLongTagsVector;
  std::vector<int> nTotalTagsVector;
  
  // all jets that pass any tag
  std::vector<DisplacedJet> isTaggedVector;
  reco::VertexCollection ivfVertices;
  
  reco::VertexCollection pVertices;
  std::vector<reco::Track> caloMatchedTracks;
  std::vector<reco::Track> vtxMatchedTracks;

  reco::Vertex selPV;

private:
  static const int GEN_STATUS_CODE_MATCH     = 23;
  static const int GEN_STATUS_CODE_MATCH_MOM = 62;

  std::vector<DisplacedJet> djets;
  int jetIDCounter = 0;     
  int debug = 0;
};

// constructor designating the calojets, primary vertices, and kinematics cuts
DisplacedJetEvent::DisplacedJetEvent(const bool& isMC, const reco::CaloJetCollection & caloJets, const reco::VertexCollection & primaryVertices, const float& minPT_, const float& minEta_, const int & debug_) {

  minPT	 = minPT_;
  minEta = minEta_;
  debug	 = debug_;

  selPV = *primaryVertices.begin();
  //const reco::Vertex & firstPV = *primaryVertices.begin();
  caloHT = 0;
  caloMET = 0;    

  nIVFReconstructed = 0;
  nIVFGenMatched    = 0;

  // construct the empty displaced jet objects to merge info later
  if (debug > 1) std::cout << "[DEBUG 1] Constructing Event From Calo jets " << std::endl;
  reco::CaloJetCollection::const_iterator jetIter = caloJets.begin();
  for(; jetIter != caloJets.end(); ++jetIter) {    
    float pt = jetIter->pt(),  eta = jetIter->eta();
    if (pt > 40 && fabs(eta) < 2.4) caloHT += pt;
    if (pt < minPT || fabs(eta) > minEta) continue;
    
    DisplacedJet djet(*jetIter, selPV, isMC, jetIDCounter, debug);
    
    djets.push_back(djet);
    jetIDCounter++;
  }  

  // merge primary vertices into event
  if (debug > 1) std::cout << "[DEBUG 1] Storing Primary Vertices " << std::endl;
  reco::VertexCollection::const_iterator pvIter = primaryVertices.begin();
  for(; pvIter != primaryVertices.end(); ++pvIter) { 
    pVertices.push_back(*pvIter);
  }
}

// checks if the event passes the preselection after checking the preselection 
// on each jet
bool DisplacedJetEvent::doesPassPreSelection() {
  if (debug > 1) std::cout << "[DEBUG] Checking Event Preselection " << std::endl;

  int nJetsPass = 0;
  std::vector<DisplacedJet>::iterator djetIter = djets.begin();
  for(; djetIter != djets.end(); ++djetIter) {
    if(djetIter->doesPassPreSelection() && djetIter->caloPt > 80.0) nJetsPass++;
  }  

  bool didPass = nJetsPass >= 2;

  if (debug > 1) std::cout << "[DEBUG] Event Preselection pass? " << didPass << std::endl;
  return didPass;
}


// performs the tagging calculation and fills the NJet tag vectors
void DisplacedJetEvent::doJetTagging(std::vector<float> noVtxThres,
				     std::vector<float> shortThres,
				     std::vector<float> medThres,
				     std::vector<float> longThres,
				     float shortTagDist, float mediumTagDist, float longTagDist, 
				     int dHTWP) {
  float nWP = noVtxThres.size();
  // check that the same number of thresholds exists for each case
  assert(((noVtxThres.size() + shortThres.size() 
	   + medThres.size() + longThres.size()) / 4 ) == noVtxThres.size());

  // the displaced ht working point must be one of the tag working points

  assert(dHTWP <= nWP);

  // fill the arrays with default values
  for(int wp = 0; wp < nWP; ++wp) {
    nNoVertexTagsVector.push_back(0);
    nShortTagsVector.push_back(0);
    nMediumTagsVector.push_back(0);
    nLongTagsVector.push_back(0);    
    nTotalTagsVector.push_back(0); 
    caloDHT.push_back(0);
  }

  if (debug > 1) std::cout << "[DEBUG] Building nTags Arrays " << std::endl;

  // loop over each jet and check for passing each threshold 
  std::vector<DisplacedJet>::iterator djetIter = djets.begin();
  for(; djetIter != djets.end(); ++djetIter) {
    // vector of bools for passing each threshold in each distance category
    std::vector<bool>	noVtxPass  = djetIter->passNoVtxTag(noVtxThres); 
    std::vector<bool>	shortPass  = djetIter->passShortTag(shortThres, shortTagDist, mediumTagDist);
    std::vector<bool>	mediumPass = djetIter->passMediumTag(medThres, mediumTagDist, longTagDist);
    std::vector<bool>	longPass   = djetIter->passLongTag(longThres, longTagDist, 9999999999);
    
    // increment the threshold counter for each
    for(int wp = 0; wp < nWP; ++wp ) {
      if (debug > 3) std::cout << "[DEBUG] TAGS PASS WP? " <<  wp << " novtx " << 
		       noVtxPass[wp] << " short " << shortPass[wp] << " medium " << 
		       mediumPass[wp] << " long " << longPass[wp] <<  std::endl;

      bool passAny = noVtxPass[wp] || shortPass[wp] || mediumPass[wp] || longPass[wp];
      if(passAny) caloDHT[wp] += djetIter->caloPt;

      // arrays of integers count the number of passes
      nNoVertexTagsVector[wp] += noVtxPass[wp]  ? 1 : 0;
      nShortTagsVector[wp]    += shortPass[wp]  ? 1 : 0;
      nMediumTagsVector[wp]   += mediumPass[wp] ? 1 : 0;
      nLongTagsVector[wp]     += longPass[wp]   ? 1 : 0;      
      nTotalTagsVector[wp]    += nNoVertexTagsVector[wp] + nShortTagsVector[wp] + nMediumTagsVector[wp] + nLongTagsVector[wp];
    }  // loop: wp working point thresholds
  } // loop: displaced jets     
}

// metric for checking gen match
float DisplacedJetEvent::genMatchMetric(const reco::GenParticle & particle, const reco::Vertex& vertex) {
  float vx = vertex.x(), vy = vertex.y(), vz = vertex.z();
  float gx = particle.vx(), gy = particle.vy(), gz = particle.vz();
  float metric = std::sqrt(((gx-vx)*(gx-vx))/(gx*gx) + ((gy-vy)*(gy-vy))/(gy*gy) + ((gz-vz)*(gz-vz))/(gz*gz));
  return metric;
}

// check how many ivf vertices are gen matched
void DisplacedJetEvent::calcNIVFGenMatched(const float & metricThreshold, const reco::GenParticleCollection& genParticles) {
  if(debug > 1) std::cout <<  "[GEN MATCH] Calculating Number of IVF are gen matched" << std::endl;

  nIVFGenMatched = 0;
  reco::VertexCollection::const_iterator vtxIter = ivfVertices.begin();
  reco::GenParticleCollection::const_iterator genIter = genParticles.begin();
  for(; vtxIter != ivfVertices.end(); ++vtxIter ) {
    float foundMatch = false;
    for(; genIter != genParticles.end(); ++genIter) {
      if (genIter->status() != GEN_STATUS_CODE_MATCH) continue;  
      float ivfMatch_temp = genMatchMetric(*genIter, *vtxIter);      
      if (ivfMatch_temp < metricThreshold) { 
	foundMatch = true;
	break;
      }
    }
    if(foundMatch) {
      nIVFGenMatched++;
    }
  }
    
  if(debug > 1) std::cout <<  "[GEN MATCH] Total Matched:" << nIVFGenMatched <<  std::endl;
}

// add the ivf vertices for each jet
void DisplacedJetEvent::addIVFVertices(const reco::VertexCollection & vertices) {
  if (debug > 1) std::cout << "[DEBUG] Adding IVF Vertices " << std::endl;  
  nIVFReconstructed = vertices.size();
  ivfVertices = vertices;
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

  /* reco::VertexCollection::const_iterator pvIter = pVertices.begin(); */
  /* reco::GenParticleCollection::const_iterator genIter = genParticles.begin(); */
  /* for(; pvIter != pVertices.end(); ++pvIter) { */
  /*   for(; genIter != genParticles.end(); ++genIter) { */
  /*     if (genIter->status() != GEN_STATUS_CODE_MATCH_MOM) continue; */
      
  /*   } */

  /*   float x  = pvIter->x(), y = pvIter->y(), z = pvIter->z(); */
  /*   float dx = x - selPV.x() , dy = y - selPV.y(), dz = z - selPV.z(); */
  /*   float metric = std::sqrt(((gx-vx)*(gx-vx))/(gx*gx) + ((gy-vy)*(gy-vy))/(gy*gy) + ((gz-vz)*(gz-vz))/(gz*gz));     */
    
  /* } */

  // fill how many vertices are generator matched
  calcNIVFGenMatched(vtxMatchThreshold, genParticles);

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
      djetIter->doGenVertexID(vtxMatchThreshold, genParticles); // allow for twice the space for vertex ID matching (BKG MOTIVATED)
    }
    if(debug > 1) std::cout <<  "[GEN MATCH] isCalomatch: " << isCaloGenMatched << " isIVFGenVertexMatched "
			  << isIVFGenVertexMatched << " isSVGEnVertexMatched " << isSVGenVertexMatched << std::endl;

  }  // end loop over displaced jets
}
