
typedef std::vector<DisplacedJet>  DisplacedJetCollection;

class DisplacedJetEvent {
 public:
  // constructor designating the calojets, primary vertices, and kinematics cuts
 DisplacedJetEvent(const bool& isMC, const reco::CaloJetCollection & caloJets, 
		   const reco::VertexCollection & primaryVertices, 
		   const float& minPT_, const float& minEta_,
		   const edm::EventSetup& iSetup_, const int & debug_) : selPV(*(primaryVertices.begin())), minPT(minPT_), minEta(minEta_), debug(debug_), iSetup(iSetup_) {    

    caloHT	      = 0;
    caloMET	      = 0;    
    nIVFReconstructed = 0;
    nIVFGenMatched    = 0;

    // construct the empty displaced jet objects to merge info later
    if (debug > 1) std::cout << "[DEBUG 1] Constructing Event From Calo jets " << std::endl;
    reco::CaloJetCollection::const_iterator jetIter = caloJets.begin();
    caloLeadingJetPT = -1;
    caloSubLeadingJetPT = -1;
    int jj = 0;
    for(; jetIter != caloJets.end(); ++jj, ++jetIter) {    
      float pt = jetIter->pt(),  eta = jetIter->eta();

      if (pt > 40 && fabs(eta) < 3.0) caloHT += pt;
      if (pt < minPT || fabs(eta) > minEta) continue;
    
      // check leading sub-leading kinematic quantities
      // if this is the first jet
      if (pt > caloLeadingJetPT && caloLeadingJetPT == -1) {
	caloLeadingJetPT = pt;
      }
      // shift down the leading jet
      else if( pt > caloLeadingJetPT && caloLeadingJetPT > 0) {
	caloSubLeadingJetPT = caloLeadingJetPT;
	caloLeadingJetPT = pt;
      }
      else if (pt > caloSubLeadingJetPT) {
	caloSubLeadingJetPT = pt;
      }       
    
      DisplacedJet djet(*jetIter,  selPV, isMC, jetIDCounter, iSetup, debug);
    
      djets.push_back(djet);
      jetIDCounter++;
    }    

    // merge primary vertices into event
    if (debug > 1) std::cout << "[DEBUG 1] Storing Primary Vertices " << std::endl;
    reco::VertexCollection::const_iterator pvIter = primaryVertices.begin();
    for(; pvIter != primaryVertices.end(); ++pvIter) { 
      pVertices.push_back(*pvIter);
    }

    // Intialize the leading subleading track variables
    caloFewestPromptTracks       = 999;
    caloSubFewestPromptTracks    = 999;
    caloMostDispTracks	       = -1;
    caloSubMostDispTracks	       = -1;
    // HLT
    caloFewestPromptTracksHLT    = 999;
    caloSubFewestPromptTracksHLT = 999;
    caloMostDispTracksHLT	       = -1;
    caloSubMostDispTracksHLT     = -1;
    // by hadronic fraction for pt > 40 GeV                                     
    caloLeadingHadronicFraction  = -1;
  }   // end constructor

  // constructor intiatlized const values
  const reco::Vertex selPV;
  const float minPT;
  const float minEta; 
  const int debug;
  const edm::EventSetup& iSetup;

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

  // helper method
  DisplacedJet & findDisplacedJetByPtEtaPhi(const float& pt, const float& eta, const float& phi);
  DisplacedJetCollection & getDisplacedJets() { return djets; }
    
  // jet associated info mergers
  void mergeTrackAssociations(const reco::JetTracksAssociation::Container&, const reco::JetTracksAssociation::Container&);
  void mergeCaloIPTagInfo(const reco::TrackIPTagInfoCollection&, const reco::VertexCollection&);
  void mergeSVTagInfo(const reco::SecondaryVertexTagInfoCollection&);
  void fillLeadingSubleadingJets(const bool & isHLT);

  void doMultiJetClustering();
  
  // associated collections
  void addIVFVertices(const reco::VertexCollection & vertices);
  
  // matching
  void doGenMatching(const reco::GenParticleCollection & genParticleCollection, 
		     const bool & doCaloJetMatch, const bool & doGenVtxMatch, const bool & doGenVtxID,
		     const float& ptMatch, const float& drMatch, const float& vtxMatchThreshold); 

  void doSimMatching(const edm::SimVertexContainer & simVertexCollection);  
  
  // kinematic variables
  float caloHT;
  std::vector<float> caloDHT;
  float caloMET;

  // ordered quantities 
  float caloLeadingJetPT,  caloSubLeadingJetPT;
  // by pt for inclusive requirements
  float caloFewestPromptTracks, caloSubFewestPromptTracks;
  float caloMostDispTracks, caloSubMostDispTracks;
  // HLT
  float caloFewestPromptTracksHLT, caloSubFewestPromptTracksHLT;
  float caloMostDispTracksHLT, caloSubMostDispTracksHLT;
  // by hadronic fraction for pt > 40 GeV
  float caloLeadingHadronicFraction;
  // vbf numbers
  // Mqq for minimum dEta 3.0 max dEta 5 and min pt 20
  float caloLeadingMqq;

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

private:
  static const int GEN_STATUS_CODE_MATCH     = 23;
  static const int GEN_STATUS_CODE_MATCH_MOM = 62;

  std::vector<DisplacedJet> djets;
  int jetIDCounter = 0;     

};


// fill leading sub leading quantities for displaced jets
void DisplacedJetEvent::fillLeadingSubleadingJets(const bool & isHLT) {
  
  // temporary storage  
  float caloFewestPromptTracks_temp = 999, caloSubFewestPromptTracks_temp = 999;
  float caloMostDispTracks_temp	    = -1, caloSubMostDispTracks_temp = -1;

  std::vector<DisplacedJet>::iterator djetIter = djets.begin();
  for(; djetIter != djets.end(); ++djetIter) {    
    int	    nPrompt = djetIter->getNPromptTracks(isHLT);
    int	    nDisp   = djetIter->getNDispTracks(isHLT);
    float   pt	    = djetIter->caloPt;

    //Highest Hadronic Fraction 
    if (djetIter->isInclusive(true)  && pt > 40.0) {
      caloLeadingHadronicFraction = std::max(djetIter->caloHadEnergyFrac, caloLeadingHadronicFraction);
    }

    // Inclusive check
    if (pt > 40.0 && djetIter->caloHadEnergyFrac > 0.01){
      // first check
      if (nPrompt < caloFewestPromptTracks_temp && caloFewestPromptTracks_temp == 999) {
	caloFewestPromptTracks_temp = nPrompt;
      }
      // shift down the leading jet
      else if( nPrompt < caloFewestPromptTracks_temp && caloFewestPromptTracks_temp >= 0) {
	caloSubFewestPromptTracks_temp = caloFewestPromptTracks_temp;
	caloFewestPromptTracks_temp    = nPrompt;
      }
      // 2nd fewest
      else if (nPrompt < caloSubFewestPromptTracks_temp) {
	caloSubFewestPromptTracks_temp = nPrompt;
      }     
    }

    // Disp Track Check (for HLT level inclusive jets)
    // include an inclusive requirement using HLT iterations 0,1,2
    if (djetIter->isInclusive(true) && pt > 40 && djetIter->caloHadEnergyFrac > 0.01){
      // first iteration
      if (nDisp > caloMostDispTracks_temp && caloMostDispTracks_temp < 0) {
	caloMostDispTracks_temp = nDisp;
      }
      // shift down the leading jet
      else if( nDisp > caloMostDispTracks_temp && caloMostDispTracks_temp > 0) {
	caloSubMostDispTracks_temp = caloMostDispTracks_temp;
	caloMostDispTracks_temp = nDisp;
      }
      // 2nd highest
      else if (pt > caloSubMostDispTracks_temp) {
	caloSubMostDispTracks_temp = nDisp;
      }     
    }
  }

  if(isHLT) {
    caloFewestPromptTracksHLT	 = caloFewestPromptTracks_temp ;
    caloSubFewestPromptTracksHLT = caloSubFewestPromptTracks_temp;
    caloMostDispTracksHLT	 = caloMostDispTracks_temp;
    caloSubMostDispTracksHLT	 = caloSubMostDispTracks_temp;
  }
  else{
    caloFewestPromptTracks    = caloFewestPromptTracks_temp ;
    caloSubFewestPromptTracks = caloSubFewestPromptTracks_temp;
    caloMostDispTracks	      = caloMostDispTracks_temp;
    caloSubMostDispTracks     = caloSubMostDispTracks_temp;
  }  
}

// checks if the event passes the preselection after checking the preselection 
// on each jet
bool DisplacedJetEvent::doesPassPreSelection() {
  if (debug > 1) std::cout << "[DEBUG] Checking Event Preselection " << std::endl;

  int nJetsPass = 0;
  std::vector<DisplacedJet>::iterator djetIter = djets.begin();
  for(; djetIter != djets.end(); ++djetIter) {
    if(djetIter->doesPassPreSelection() && djetIter->caloPt > 40.0) nJetsPass++;
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

void DisplacedJetEvent::mergeTrackAssociations(const reco::JetTracksAssociation::Container& caloMatched, const reco::JetTracksAssociation::Container& vertexMatched) {
  std::vector<reco::JetBaseRef> caloJets   = reco::JetTracksAssociation::allJets(caloMatched);
  std::vector<reco::JetBaseRef> vertexJets = reco::JetTracksAssociation::allJets(vertexMatched);

  if(caloJets.size() != vertexJets.size()){
    throw cms::Exception("TrackAssociationJetMatchingFailure");
  }

  // loop over the jets contained in the caloJet association
  // add the calo matched track information stored in the assocation 
  std::vector<reco::JetBaseRef>::const_iterator jetIter = caloJets.begin();
  for(; jetIter != caloJets.end(); ++jetIter){
    float pt = (*jetIter)->pt(), eta = (*jetIter)->eta(), phi = (*jetIter)->phi();    
    if (pt < minPT || fabs(eta) > minEta) continue;    // dont look for jets outside of acceptance
    DisplacedJet & djet = findDisplacedJetByPtEtaPhi(pt, eta, phi);
    djet.addCaloTrackInfo(reco::JetTracksAssociation::getValue(caloMatched, *jetIter));
  }

  // loop over the jets contained in the vertexMatched association
  // add the  track information stored in the assocation 
  jetIter = vertexJets.begin();
  for(; jetIter != vertexJets.end(); ++jetIter){
    float pt = (*jetIter)->pt(), eta = (*jetIter)->eta(), phi = (*jetIter)->phi();    
    if (pt < minPT || fabs(eta) > minEta) continue;    // dont look for jets outside of acceptance
    DisplacedJet & djet = findDisplacedJetByPtEtaPhi(pt, eta, phi);
    djet.addVertexTrackInfo(reco::JetTracksAssociation::getValue(vertexMatched, *jetIter));
  }
}

// fills each jet with ip tag info variables 
// fills each jet track collection with tracks used for ip variables
void DisplacedJetEvent::mergeCaloIPTagInfo(const reco::TrackIPTagInfoCollection & ipTagInfo, const reco::VertexCollection & primaryVertices) {
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
    // the track colleciton above can be retrieved and passed to track angles below
    djet.addHitInfo(djet.getVertexMatchedTracks());    
    djet.addTrackAngles(djet.getVertexMatchedTracks());    
    // calculate the v0 candidates
    djet.addV0Info(djet.getVertexMatchedTrackRefs());    
    // calculate alpha for the vertces
    djet.calcJetAlpha(djet.getVertexMatchedTracks(), primaryVertices);
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

void DisplacedJetEvent::doMultiJetClustering() {
  // do clustering for each jet in the event
  std::vector<DisplacedJet>::iterator djetIter = djets.begin();
  for(; djetIter != djets.end(); ++djetIter) {
    // clone the list of jets and remove the current jet
    std::vector<DisplacedJet> temp;
    std::vector<DisplacedJet>::iterator neighbor = djets.begin();
    // add all the neighbor jets in the event for the clustering
    for(; neighbor != djets.end(); ++neighbor) {
      if (neighbor == djetIter) continue;
      temp.push_back(*neighbor);
    }
    // do the clustering comparing the vertices in the jet to the vertices in 
    // all other jets in the event
    djetIter->calcNJetClusterSize(djetIter->displacedV0VectorCleaned, temp, 2);    
  }  
}

