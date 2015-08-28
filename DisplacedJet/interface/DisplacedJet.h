
class DisplacedJet {
 public:
  DisplacedJet(const reco::CaloJet & jet_, const float alpha_, const reco::Vertex & primaryVertex, const bool & isMC_, const int& jetID_, const int & debug_) {
    debug = debug_;

    jet	  = jet_;
    isMC  = isMC_;
    jetID = jetID_;

    // primary vertex
    selPV = primaryVertex;

    // track association variables
    nTracks       = 0;

    // initialize calo related variables
    caloPt	  = jet.pt();
    caloPx        = jet.px();
    caloPy        = jet.py();
    caloPz        = jet.pz();
    caloEta	  = jet.eta();
    caloPhi	  = jet.phi();
    caloN60	  = 0; //jet.n60();
    caloN90	  = 0; //jet.n90();
    caloTowerArea = 0; //jet.towersArea();

    // alpha
    alpha = 0;
    
    // store quantites based on detector geometry (rather than vprimary vertex)
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> detP4 = jet.detectorP4();
    detPt  = detP4.pt();
    detEta = detP4.eta();
    detPhi = detP4.phi();
    
    // energy fractions for tagging very long lived jets
    caloEMEnergyFrac  = jet.emEnergyFraction();
    caloHadEnergyFrac = jet.energyFractionHadronic();
    
    // initialize ip related combination variables
    ipSigLogSum2D    = 0, ipSigLogSum3D = 0;
    ipLogSum2D	     = 0, ipLogSum3D = 0;
    jetDistLogSum    = 0, jetDistLogSum = 0;
    jetDistSigLogSum = 0, jetDistSigLogSum = 0;
    eipSigLogSum2D   = 0, eipSigLogSum3D = 0;
    eipSigSum2D	     = 0, eipSigSum3D = 0;

    // initalize ip distributional variable
    //mean
    meanIPSig2D	       = 0, meanIPSig3D = 0;
    meanIP2D	       = 0, meanIP3D = 0;
    meanIPLogSig2D     = 0, meanIPLogSig3D = 0;
    meanIPLog2D	       = 0, meanIPLog3D = 0;
    meanJetDist	       = 0, meanJetDistSig = 0;
    //median
    medianIPSig2D      = 0, medianIPSig3D = 0;  
    medianIP2D	       = 0, medianIP3D = 0;  
    medianIPLogSig2D   = 0, medianIPLogSig3D = 0;  
    medianIPLog2D      = 0, medianIPLog3D = 0;  
    medianJetDist      = 0, medianJetDistSig = 0;
    //variance
    varianceIPSig2D    = 0, varianceIPSig3D = 0;
    varianceIP2D       = 0, varianceIP3D = 0;
    varianceIPLogSig2D = 0, varianceIPLogSig3D = 0;
    varianceIPLog2D    = 0, varianceIPLog3D = 0;
    varianceJetDist    = 0, varianceJetDistSig = 0;    

    // ivf init
    ivfX = 0, ivfY = 0, ivfZ = 0;
    ivfXError = 0, ivfYError = 0, ivfZError = 0;
    ivfLxy = 0, ivfLxyz = 0;
    ivfLxySig = 0, ivfLxyzSig = 0;
    ivfMass = 0;
    ivfNTracks = 0;
    ivfIsGenMatched = false;
    ivfIsSimMatched = false;
    ivfGenVertexMatchMetric = FAKE_HIGH_NUMBER;
    selIVFIsPVScore = FAKE_HIGH_NUMBER;
    selIVFIsPV = true;
    
    // sv init
    svNVertex= 0;
    svX = 0, svY = 0, svZ = 0;
    svXError = 0, svYError = 0, svZError = 0;
    svChi2 = 0, svNDof = 0;
    svLxy = 0, svLxyz = 0;
    svLxySig = 0, svLxyzSig = 0;
    svMass = 0;
    svNTracks = 0;
    svIsGenMatched = false;
    svIsSimMatched = false;
    svGenVertexMatchMetric = FAKE_HIGH_NUMBER;    

    //angles
    sumTrackPt = 0;
    // pt weighted
    ptSumCosTheta2D	  = -99, ptSumCosTheta3D = -99;
    ptSumCosThetaDet2D	  = -99, ptSumCosThetaDet3D = -99;
    // aboslute sum
    sumCosTheta2D	  = -99, sumCosTheta3D = -99;
    sumCosThetaDet2D	  = -99, sumCosThetaDet3D = -99;
    // mean
    meanCosTheta2D	  = -99, meanCosTheta3D = -99;
    meanCosThetaDet2D	  = -99, meanCosThetaDet3D = -99;
    // median
    medianCosTheta2D	  = -99, medianCosTheta3D = -99;
    medianCosThetaDet2D	  = -99, medianCosThetaDet3D = -99;
    // variance
    varianceCosTheta2D	  = -99, varianceCosTheta3D = -99;
    varianceCosThetaDet2D = -99, varianceCosThetaDet3D = -99;

    trackSumMomCosTheta2D = -99;
    trackSumMomCosTheta3D = -99;
    trackSumMomMag2D	  = -99;
    trackSumMomMag3D	  = -99;

    ipPosSumMag3D	  = -99;
    ipPosSumMag2D	  = -99;
    // matching gen particles init
    caloGenPt		  = -999;
    caloGenEta		  = -999;
    caloGenPhi		  = -999;
    caloGenMass		  = -999;
    isCaloGenMatched	  = 0;
    // ivf mc id
    ivfHighestPtMomPt	  = 0;
    ivfHighestPtMomID	  = 0;

    // 
    // HIT RELATED
    //
    // overall distribution
    jetMedianInnerHitPos	   = -1;
    jetMedianOuterHitPos	   = -1;
    jetMeanInnerHitPos		   = -1;
    jetMeanOuterHitPos		   = -1;
    jetVarianceInnerHitPos	   = -1;
    jetVarianceOuterHitPos	   = -1;
    // distributions from inside the pixel layers
    jetMedianInnerHitPosInPixel	   = -1;
    jetMedianOuterHitPosInPixel	   = -1;
    jetMeanInnerHitPosInPixel	   = -1;
    jetMeanOuterHitPosInPixel	   = -1;
    jetVarianceInnerHitPosInPixel  = -1;
    jetVarianceOuterHitPosInPixel  = -1;
    // distributions outside the pixel layers
    jetMedianInnerHitPosOutPixel   = -1;
    jetMedianOuterHitPosOutPixel   = -1;
    jetMeanInnerHitPosOutPixel	   = -1;
    jetMeanOuterHitPosOutPixel	   = -1;
    jetVarianceInnerHitPosOutPixel = -1;
    jetVarianceOuterHitPosOutPixel = -1;
    // fraction valid hits
    jetMedianTrackValidHitFrac     = -1;
    jetMeanTrackValidHitFrac	   = -1;
    jetVarianceTrackValidHitFrac   = -1;
    // track counting
    jetNTracksNoPixel		   = -1;
    jetNTracksPixel		   = -1;
    jetPtSumTracksNoPixel	   = -1;
    jetPtSumTracksPixel		   = -1;
  }

  // tag related
  std::vector<bool> passNoVtxTag(const std::vector<float> thres);
  std::vector<bool> passShortTag(const std::vector<float> thres, float min, float max);
  std::vector<bool> passMediumTag(const std::vector<float> thres, float min, float max);
  std::vector<bool> passLongTag(const std::vector<float> thres, float min, float max);

  // selection related
  bool doesPassPreSelection() const;
  int  getNDispTracks(const bool & isHLT) const;
  int  getNPromptTracks(const bool & isHLT) const;
  bool isInclusive(const bool & isHLT) const;
  bool isDispTrack(const bool & isHLT) const;

  // jet info integration
  void addCaloTrackInfo(const reco::TrackRefVector&);
  void addVertexTrackInfo(const reco::TrackRefVector&);
  void addIVFCollection(const reco::VertexCollection&, const float& compatibilityScore);
  void addIPTagInfo(const reco::TrackIPTagInfo&);
  void addSVTagInfo(const reco::SecondaryVertexTagInfo&);
  
  // angular information relative to ejt
  void addTrackAngles(const reco::TrackCollection tracks, const edm::EventSetup& iSetup);
  void addHitInfo(const reco::TrackCollection tracks, const edm::EventSetup& iSetup);

  // generator matching
  bool doGenCaloJetMatching(const float& ptMatch, const float& dRMatch, const reco::GenParticleCollection& genParticles);
  bool doGenVertexJetMatching(const float& metricThreshold, const reco::GenParticleCollection& genParticles);
  void doGenVertexID(const float& metricThreshold, const reco::GenParticleCollection& genParticles);
  float genMatchMetric(const reco::GenParticle & particle, const reco::Vertex& vertex);
  
  // jet info extraction
  reco::TrackCollection getCaloMatchedTracks() { return caloMatchedTracks; }
  reco::TrackCollection getVertexMatchedTracks() { return vertexMatchedTracks; }
  reco::Vertex          getIVFVertexSelected() { return selIVF; }
  reco::Vertex          getSVVertex() { return selSV; }
  
  // jet distribution calculator
  float getJetMedian(const std::vector<float>&, bool);
  float getJetMean(const std::vector<float>&, bool);
  float getJetVariance(const std::vector<float>&, bool);

  //alpha calculation
  void calcJetAlpha(const reco::TrackCollection&, const reco::VertexCollection&);

  //////////////CALO INFORMATION////////////
  bool isMC;
  int jetID;

  // track association variables
  int nTracks;
  float sumTrackPt;  
  // calo related variables
  float caloPt, caloEta, caloPhi;
  float caloPx, caloPy, caloPz;
  float caloN60, caloN90;  
  float caloTowerArea;
  float detPt, detEta, detPhi;
  float caloEMEnergyFrac, caloHadEnergyFrac;

  float alpha, alphaMax;

  //////////////JET IP VARIABLES/////////////

  // ip related combined variables
  float ipSigLogSum2D, ipSigLogSum3D;
  float ipLogSum2D, ipLogSum3D;  
  float jetDistLogSum, jetDistSigLogSum;
  float eipSigLogSum2D, eipSigLogSum3D;
  float eipSigSum2D, eipSigSum3D;

  // ip distribution variables
  // mean 
  float meanIPSig2D, meanIPSig3D;
  float meanIP2D, meanIP3D;
  float meanIPLogSig2D, meanIPLogSig3D;
  float meanIPLog2D, meanIPLog3D;
  float meanJetDist, meanJetDistSig;
  // median
  float medianIPSig2D, medianIPSig3D;  
  float medianIP2D, medianIP3D;
  float medianIPLogSig2D, medianIPLogSig3D;  
  float medianIPLog2D, medianIPLog3D;
  float medianJetDist, medianJetDistSig;
  // variance
  float varianceIPSig2D, varianceIPSig3D;  
  float varianceIP2D, varianceIP3D;
  float varianceIPLogSig2D, varianceIPLogSig3D;  
  float varianceIPLog2D, varianceIPLog3D;
  float varianceJetDist, varianceJetDistSig;

  /////////////HIT RELATED////////////////

  float jetMedianInnerHitPos;
  float jetMedianOuterHitPos;
  float jetMeanInnerHitPos;
  float jetMeanOuterHitPos;
  float jetVarianceInnerHitPos;
  float jetVarianceOuterHitPos;
  // distributions from inside the pixel layers
  float jetMedianInnerHitPosInPixel;
  float jetMedianOuterHitPosInPixel;
  float jetMeanInnerHitPosInPixel;
  float jetMeanOuterHitPosInPixel;
  float jetVarianceInnerHitPosInPixel;
  float jetVarianceOuterHitPosInPixel;
  // distributions outside the pixel layers
  float jetMedianInnerHitPosOutPixel;
  float jetMedianOuterHitPosOutPixel;
  float jetMeanInnerHitPosOutPixel;
  float jetMeanOuterHitPosOutPixel;
  float jetVarianceInnerHitPosOutPixel;
  float jetVarianceOuterHitPosOutPixel;
  // fraction valid hits
  float jetMedianTrackValidHitFrac;
  float jetMeanTrackValidHitFrac;
  float jetVarianceTrackValidHitFrac ;
  // track counting
  float jetNTracksNoPixel;
  float jetNTracksPixel;
  float jetPtSumTracksNoPixel;
  float jetPtSumTracksPixel;

  //////////////VERTEX VARIABLES//////////////

  // ivf related variables
  // position
  float ivfX, ivfY, ivfZ;
  float ivfXError, ivfYError, ivfZError;  
  float ivfLxy, ivfLxyz;
  float ivfLxySig, ivfLxyzSig;
  // qualities
  float ivfMass;
  float ivfNTracks;
  // matching
  bool  ivfIsGenMatched;
  bool  ivfIsSimMatched;
  float ivfMatchingScore;
  float ivfGenVertexMatchMetric;

  // sv related variables
  // multiple vertices
  int svNVertex;
  // position
  float svX, svY, svZ;
  float svXError, svYError, svZError;  
  float svChi2, svNDof;
  float svLxy, svLxyz;
  float svLxySig, svLxyzSig;
  // qualities
  float svMass;
  float svNTracks;
  // matching
  bool  svIsGenMatched;
  bool  svIsSimMatched;   
  float svGenVertexMatchMetric;

  //////////////TRACK ANGLE VARIABLES
  

  // pt weighted
  float ptSumCosTheta2D, ptSumCosTheta3D;
  float ptSumCosThetaDet2D, ptSumCosThetaDet3D;
  // aboslute sum
  float sumCosTheta2D, sumCosTheta3D;
  float sumCosThetaDet2D, sumCosThetaDet3D;
  // mean
  float meanCosTheta2D, meanCosTheta3D;
  float meanCosThetaDet2D, meanCosThetaDet3D;
  // median
  float medianCosTheta2D, medianCosTheta3D;
  float medianCosThetaDet2D, medianCosThetaDet3D;
  // variance
  float varianceCosTheta2D, varianceCosTheta3D;
  float varianceCosThetaDet2D, varianceCosThetaDet3D;

  
  float trackSumMomCosTheta2D, trackSumMomCosTheta3D;
  float trackSumMomMag2D, trackSumMomMag3D;

  float ipPosSumMag3D, ipPosSumMag2D;

  //////////////GEN MATCH VARIABLES///////////

  bool isCaloGenMatched;
  bool isCaloPVGenMatched;
  // matched gen particle kinematics
  float caloGenPt, caloGenEta, caloGenPhi, caloGenMass;
  
  ////////////// IVF ID VARAIBLES /////////

  float	ivfHighestPtMomPt;
  int ivfHighestPtMomID;
  std::vector<std::pair<const int, const float>> genMomVector; 
  std::vector<std::pair<const int, const float>> genSonVector; 

  std::vector<bool> noVertexTagsVector;
  std::vector<bool> shortTagsVector;
  std::vector<bool> mediumTagsVector;
  std::vector<bool> longTagsVector;

 private: 
  int debug;
  static const int GEN_STATUS_CODE_MATCH = 23; 
  const float FAKE_HIGH_NUMBER = 999999999;
  // calo jet the displaced jet is built upon

  reco::CaloJet jet;

  // related vertices
  reco::Vertex selIVF;
  bool selIVFIsPV;
  float selIVFIsPVScore;
  reco::Vertex selSV;
  reco::Vertex selPV;

  // related track collections
  reco::TrackCollection caloMatchedTracks; 
  reco::TrackCollection vertexMatchedTracks; 
  std::vector<reco::btag::TrackIPData> lifetimeIPData; 
  std::vector<float> ip3dVector, ip3dsVector, ip2dVector, ip2dsVector;
  std::vector<int>   trackAlgo;
  std::vector<float> ipLog3dVector, ipLog3dsVector, ipLog2dVector, ipLog2dsVector;
  std::vector<float> jetAxisDistVector, jetAxisDistSigVector; 

  // cos related
  std::vector<float> cosTheta2DVector, cosThetaDet2DVector, cosTheta3DVector, cosThetaDet3DVector; 
};

float DisplacedJet::genMatchMetric(const reco::GenParticle & particle, const reco::Vertex& vertex) {
  float vx = vertex.x(), vy = vertex.y(), vz = vertex.z();
  float gx = particle.vx(), gy = particle.vy(), gz = particle.vz();
  float metric = std::sqrt(((gx-vx)*(gx-vx))/(gx*gx) + ((gy-vy)*(gy-vy))/(gy*gy) + ((gz-vz)*(gz-vz))/(gz*gz));  
  return metric;
}

bool DisplacedJet::doGenVertexJetMatching(const float & metricThreshold, const reco::GenParticleCollection& genParticles) {
  svIsGenMatched = false;
  ivfIsGenMatched = false;

  // is IVF is consistent with the primary vertex return false
  if (selIVFIsPV) {
    return false;
  }

  // loop over each gen particle
  reco::GenParticleCollection::const_iterator genIter = genParticles.begin();
  for(; genIter != genParticles.end(); ++genIter) {
    if (genIter->status() != GEN_STATUS_CODE_MATCH) continue;
    
    // calculate the score
    float   ivfMatch_temp   = genMatchMetric(*genIter, selIVF);
    float   svMatch_temp    = genMatchMetric(*genIter, selSV);
    // check for a match
    bool    ivfMatch	    = fabs(ivfMatch_temp) < metricThreshold;
    bool    svMatch	    = fabs(svMatch_temp) < metricThreshold;
    // only one needs to match
    ivfIsGenMatched	    = ivfIsGenMatched || ivfMatch;
    svIsGenMatched	    = svIsGenMatched || svMatch;
    // store the best matching Metric
    ivfGenVertexMatchMetric = std::min(ivfMatch_temp, ivfGenVertexMatchMetric);
    svGenVertexMatchMetric  = std::min(svMatch_temp, svGenVertexMatchMetric);      

    if (ivfMatch || svMatch) {            
      return true; 
    }
  } // loop over gen particles

  return false; 
}

void DisplacedJet::doGenVertexID(const float & metricThreshold, const reco::GenParticleCollection& genParticles) {
  // loop over each gen particle
  reco::GenParticleCollection::const_iterator genIter = genParticles.begin();
  if (debug > 2) std::cout << "\n[GEN VTX MATCHING 2] STARTING NEW JET MATCHING FOR IVF / SV --- IVF LXY: " << ivfLxy << std::endl;
  if (debug > 2) std::cout << "[GEN VTX MATCHING 2] IVF vs PV CONSITENCY SCORE: " << selIVFIsPVScore << std::endl;
  for(; genIter != genParticles.end(); ++genIter) {
    const int& status = genIter->status();
    const int& charge = genIter->charge();
    const int& id     = genIter->pdgId();
    const float& pt   = genIter->pt();
    if (charge == 0) continue; //hardest subprocess or neutral particle

    const reco::Candidate * mom  = genIter->mother();
    int momid = 0;
    float mompt = 0;

    if (mom != NULL) {
      momid = mom->pdgId(); //check for null pointer
      mompt = mom->pt();
    }

    // dont do matching if there is no IVF (not consistent with the PV)
    if(selIVFIsPV) {
      if (debug > 2 && selIVFIsPVScore == 0 ) std::cout << "[GEN VTX MATCHING 3] IVF not found" <<  std::endl;
      if (debug > 2 && selIVFIsPVScore > 0 ) std::cout << "[GEN VTX MATCHING 3] IVF consistent with PV" <<  std::endl;
      break; // if there is no IVF selected (the PV is selected) then dont do matching
    }

    // calculate the score
    float   ivfMatch_temp   = genMatchMetric(*genIter, selIVF);
    float   svMatch_temp    = genMatchMetric(*genIter, selSV);
    // check for a match
    bool    ivfMatch	    = fabs(ivfMatch_temp) < metricThreshold;

    if (ivfMatch ) {
      if (debug > 2) std::cout << "[GEN VTX MATCHING 2] STATUS:" << 
		       status << " ID " << id  << " MOM ID: " << momid << 
		       " MATCHING METRIC: " << std::min(fabs(ivfMatch_temp), fabs(svMatch_temp)) << std::endl;

      if (mompt > ivfHighestPtMomPt && mom != NULL) { // keep the highest pt mom
	ivfHighestPtMomPt = mompt;
	ivfHighestPtMomID = momid;
      }

      // add the son (always)
      if (debug > 5) std::cout << "[GEN VTX MATCHING 5] PUSHING BACK SON  PAIR:" << id << " pt: " << pt << std::endl;
      genSonVector.push_back(std::make_pair(id, pt));

      // check if the mother is already in the vecetor
      bool found_my_mother = false;
      std::vector<std::pair<const int, const float>>::const_iterator momIter =  genMomVector.begin();
      for(; momIter != genMomVector.end(); ++momIter) {
	int	momid_temp = momIter->first;
	float	mompt_temp = momIter->second;

	if (momid_temp == momid && mompt_temp == mompt) { // found a match
	  found_my_mother = true;
	  break;
	}	  
      }
      if (found_my_mother) { // found a match, skip to next particle
	if (debug > 2) std::cout << "[GEN VTX MATCHING 2] SON ALREADY HAS A MOTHER " << std::endl;	
	continue; // 
      }
	
      if (debug > 2) std::cout << "[GEN VTX MATCHING 2] ADDING MOTHER: " << momid <<  " PT :" << pt << std::endl;	      
      genMomVector.push_back(std::make_pair(momid, mompt)); // add the id momentum pair for the mom

    } // close ivf match    
  } // loop over gen particles
}

bool DisplacedJet::doesPassPreSelection() const {
  if (debug > 2) std::cout << "[DEBUG] Checking Jet Preselection  " << std::endl;
  return isInclusive(true) && isDispTrack(true);
}

// returns number of tracks found in iter 0,1,2 (run at HLT)
int DisplacedJet::getNPromptTracks(const bool& isHLT) const {  
  int nPTracks = 0;
  std::vector<float>::const_iterator ip2DIter = ip2dVector.begin();
  int index = 0;
  for(; ip2DIter != ip2dVector.end(); ++ip2DIter, ++index) {
    int algo  = trackAlgo[index];
    // hlt only runs iter 0,1,2 for prompt track calculations
    // algo != iteration
    bool pass_algo = (algo <= 6 || !isHLT);
    if (fabs(*ip2DIter) < 0.1 && pass_algo) nPTracks++;
  }  
  return  nPTracks;
}

//returns the number of tracks with ip significance >5 
//tracks must be found in iter 0,1,2, or 4 
int DisplacedJet::getNDispTracks(const bool & isHLT) const {
  int	nDTracks	    = 0;
  std::vector<float>::const_iterator ip2DSigIter = ip2dsVector.begin();
  int index = 0;
  for(; ip2DSigIter != ip2dsVector.end(); ++ip2DSigIter, ++index) {
    int algo  = trackAlgo[index];
    // if this is for HLT then only use iter 0,1,2,4
    // algo != iteration
    bool pass_algo =  (algo != 7 && algo < 8) || !isHLT;
    if (fabs(*ip2DSigIter) > 5.0 && pass_algo) nDTracks++;
  }
  return nDTracks;
}

// returns true if the inclusive criteria is satisfied (at most 2 prompt tracks)
bool DisplacedJet::isInclusive(const bool & isHLT) const {
  return isHLT ? (getNPromptTracks(true) <= 2) : (getNPromptTracks(false) <= 2) ;
}

bool DisplacedJet::isDispTrack(const bool & isHLT) const {
  return isHLT ? (getNDispTracks(true) >= 2) : (getNDispTracks(false) >= 2) ;
}

bool DisplacedJet::doGenCaloJetMatching(const float& ptMatch, const float& dRMatch, const reco::GenParticleCollection& genParticles) {
  if (debug > 2) std::cout << "\n[GEN MATCHING] CALO JET INFO  " << " pt " << caloPt << " eta " << caloEta <<  " phi " << caloPhi << std::endl;
  reco::GenParticleCollection::const_iterator genIter = genParticles.begin();
  
  for(; genIter != genParticles.end(); ++genIter) {

    int id = genIter->pdgId();
    int st = genIter->status();
  
    if (st != GEN_STATUS_CODE_MATCH) continue;

    // const reco::Candidate * mom = genIter->mother();                                         
    double genPt = genIter->pt(), genEta = genIter->eta(), genPhi = genIter->phi(), genMass = genIter->mass();
    float   dr     = reco::deltaR(genEta, genPhi, caloEta, caloPhi);
    float   dpt    = fabs(caloPt - genPt) / genPt;
    if (debug > 2) std::cout << "[GEN MATCHING] Candidate ID " << id <<   " dR " << dr << " dpt" << dpt << " eta " << genEta << " phi " << genPhi << std::endl; 

    // found a match
    if (dr < dRMatch && dpt < ptMatch) {
      if(debug > 1) std::cout << "[GEN MATCHED] id " << id << " status " << st << " pt " << genPt << " eta " << genEta  <<  " phi " << genPhi << std::endl;

      isCaloGenMatched = true;
      caloGenPt	       = genPt;
      caloGenEta       = genEta;
      caloGenPhi       = genPhi;
      caloGenMass      = genMass;

      return true;
    }
  } // loop over gen particles 
  
  return false;
}

//keep the collection of tracks matched at the calo surface
void DisplacedJet::addCaloTrackInfo(const reco::TrackRefVector & trackRefs) {
  if (debug > 2) std::cout << "[DEBUG] Adding Track Info  " << std::endl;
    reco::TrackRefVector::const_iterator trackIter = trackRefs.begin();
    for(; trackIter != trackRefs.end(); ++trackIter) {      
      caloMatchedTracks.push_back(**trackIter);
    }    
}

// count the number of tracks based on the association at the vertex
void DisplacedJet::addVertexTrackInfo(const reco::TrackRefVector & trackRefs) {
  if (debug > 2) std::cout << "[DEBUG] Adding Track Info  " << std::endl;
  reco::TrackRefVector::const_iterator trackIter = trackRefs.begin();
  nTracks = 0;
  sumTrackPt = 0;
  for(; trackIter != trackRefs.end(); ++trackIter) {
    float pt = (*trackIter)->pt();
    if(pt > 1.0) {
      nTracks++;
      sumTrackPt += pt;
      vertexMatchedTracks.push_back(**trackIter);
    }
  }    
}

void DisplacedJet::addIVFCollection(const reco::VertexCollection & vertices, const float & pvCompatibilityScore = .05) {
  if (debug > 2) std::cout << "[DEBUG] Building Jet Vertex Association  " << std::endl;
  // build the jet vertex association
  JetVertexAssociation JVAIVF("IVF", selPV, debug);    
  reco::VertexCollection::const_iterator vtxIter = vertices.begin();
  for(; vtxIter != vertices.end(); ++vtxIter ) {
    JVAIVF.addVertex(*vtxIter);
  }

  if (debug > 2) std::cout << "[DEBUG 2] Extracting Best Association  " << std::endl;
  // best vertex selection from JVA 
  const std::pair<reco::Vertex, float>    bestVertexPair  = JVAIVF.getBestVertex(jet, "oneOverR");
  const reco::Vertex                      bestVertex      = bestVertexPair.first;
  const float                             bestVertexScore = bestVertexPair.second;

  // set the global IVF to the best vertex
  selIVF = bestVertex;

  if (debug > 2) std::cout << "[DEBUG 2] Best IVF Vertex Score:" << bestVertexScore << std::endl;
  //flight distance from the firstPV
  float x  = bestVertex.x(), y = bestVertex.y(), z = bestVertex.z();    
  float dx = x - selPV.x() , dy = y - selPV.y(), dz = z - selPV.z();

  // set the compatibility score
  selIVFIsPVScore = std::sqrt((dx/x)*(dx/x) + (dy/y)*(dy/y) + (dz/z)*(dz/z));  
  selIVFIsPV = selIVFIsPVScore < pvCompatibilityScore; // default 5% consitency check 
  
  //build the total error
  float svxE = bestVertex.xError(), svyE = bestVertex.yError(), svzE = bestVertex.zError();
  float pvxE = selPV.xError(), pvyE = selPV.yError(), pvzE = selPV.zError();
  float xE   = std::sqrt(svxE * svxE + pvxE * pvxE), yE = std::sqrt(svyE * svyE + pvyE * pvyE), zE = std::sqrt(svzE * svzE + pvzE * pvzE);

  ivfX	     = selIVFIsPV ? 0 : x;
  ivfY	     = selIVFIsPV ? 0 : y;
  ivfZ	     = selIVFIsPV ? 0 : z;
  ivfXError  = selIVFIsPV ? 0 : svxE;
  ivfYError  = selIVFIsPV ? 0 : svyE;
  ivfZError  = selIVFIsPV ? 0 : svzE;  
  ivfNTracks = selIVFIsPV ? 0 : bestVertex.nTracks();  

  ivfMass    = selIVFIsPV ? 0 : bestVertex.p4().mass();    
  ivfLxySig  = selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy ) / std::sqrt(xE * xE + yE * yE);
  ivfLxyzSig = selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy + dz * dz) / std::sqrt(xE * xE + yE * yE + zE * zE);
  ivfLxy     = selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy );
  ivfLxyz    = selIVFIsPV ? 0 : std::sqrt( dx * dx + dy * dy + dz * dz );

  // matching score
  ivfMatchingScore = bestVertexScore;
}
// compute variables related to the track angles and the calo jet 
void DisplacedJet::addHitInfo(const reco::TrackCollection tracks, const edm::EventSetup& iSetup) {
  reco::TrackCollection::const_iterator tIter = tracks.begin();
  // iterative of the tracks
  std::vector<float> fractionValidHits;
  std::vector<float> innerHitPos;
  std::vector<float> outerHitPos;
  std::vector<float> innerHitPosInPixel;
  std::vector<float> innerHitPosOutPixel;
  std::vector<float> outerHitPosInPixel;
  std::vector<float> outerHitPosOutPixel;

  jetNTracksNoPixel	= 0;
  jetNTracksPixel	= 0;
  jetPtSumTracksNoPixel = 0;
  jetPtSumTracksPixel   = 0;

  for(; tIter != tracks.end(); ++tIter) {
    static GetTrackTrajInfo getTrackTrajInfo;
    const reco::Track & const_track = *tIter;
    std::vector<GetTrackTrajInfo::Result> trajInfo = getTrackTrajInfo.analyze(iSetup, const_track);

    float trackNValidHits = 0;
    std::vector<GetTrackTrajInfo::Result>::const_iterator ti = trajInfo.begin();
    for( ; ti != trajInfo.end(); ++ti) {
      if ((*ti).valid) trackNValidHits += 1;
    }
    
    fractionValidHits.push_back(float(trackNValidHits) / float(trajInfo.size()));   

    // check the inner and outer hit
    if (trajInfo[0].valid && trajInfo.back().valid) {

      // get the detector layer
      const DetLayer &  detLayerInner	 = *(trajInfo[0].detLayer); 
      const DetLayer &  detLayerOuter	 = *(trajInfo.back().detLayer); 
      
      //const bool	isBarInner	 = detLayerInner.isBarrel();
      //const bool	isBarOuter	 = detLayerOuter.isBarrel();
      GeomDetEnumerators::SubDetector subDetLayerInner = detLayerInner.subDetector();
      GeomDetEnumerators::SubDetector subDetLayerOuter = detLayerOuter.subDetector();
      
      // get the state on the surface
      const TrajectoryStateOnSurface&   tsosInnerHit  = trajInfo[0].detTSOS;
      const TrajectoryStateOnSurface&   tsosOuterHit  = trajInfo.back().detTSOS;
      // position of the hit
      const GlobalPoint &               innerPos      = tsosInnerHit.globalPosition();
      const GlobalPoint &               outerPos      = tsosOuterHit.globalPosition();

      float pt = tIter->pt(); 

      TVector3 outerv3( outerPos.x(), outerPos.y(), outerPos.z());
      TVector3 outerv2( outerPos.x(), outerPos.y(), 0);
      TVector3 ipv3( tIter->vx(), tIter->vy(), tIter->vz());
      TVector3 ipv2( tIter->vx(), tIter->vy(), 0);
      TVector3 innerv3( innerPos.x(), innerPos.y(), innerPos.z());
      TVector3 innerv2( innerPos.x(), innerPos.y(), 0);            

      float ri = innerv2.Mag();
      float ro = outerv2.Mag();
      
      innerHitPos.push_back(ri);
      outerHitPos.push_back(ro);

      bool hasNoPixelInner = subDetLayerInner!=GeomDetEnumerators::PixelBarrel && subDetLayerInner!=GeomDetEnumerators::PixelEndcap;
      bool hasNoPixelOuter = subDetLayerOuter!=GeomDetEnumerators::PixelBarrel && subDetLayerOuter!=GeomDetEnumerators::PixelEndcap;

      // if the inner hit is in the pixel layers
      if (!hasNoPixelInner) {
	innerHitPosInPixel.push_back(ri);
	jetNTracksPixel	    += 1;	
	jetPtSumTracksPixel += pt;
      }
      else {
	innerHitPosOutPixel.push_back(ri);
	jetNTracksNoPixel     += 1;
	jetPtSumTracksNoPixel += pt;	  
      }

      if (!hasNoPixelOuter) outerHitPosInPixel.push_back(ro);
      else outerHitPosOutPixel.push_back(ro);     
    } // end hits valid
  } //end loop filling vectors of hits

  // fill the jet parameters
  // overall distribution
  jetMedianInnerHitPos		 = getJetMedian(innerHitPos, false); 
  jetMedianOuterHitPos		 = getJetMedian(outerHitPos, false);
  jetMeanInnerHitPos		 = getJetMean(innerHitPos, false);
  jetMeanOuterHitPos		 = getJetMean(outerHitPos, false);
  jetVarianceInnerHitPos	 = getJetVariance(innerHitPos, false);
  jetVarianceOuterHitPos	 = getJetVariance(outerHitPos, false);
  // distributions from inside the pixel layers
  jetMedianInnerHitPosInPixel	 = getJetMedian(innerHitPosInPixel, false);
  jetMedianOuterHitPosInPixel	 = getJetMedian(outerHitPosInPixel, false);
  jetMeanInnerHitPosInPixel	 = getJetMean(innerHitPosInPixel, false);
  jetMeanOuterHitPosInPixel	 = getJetMean(outerHitPosInPixel, false);
  jetVarianceInnerHitPosInPixel  = getJetVariance(innerHitPosInPixel, false);
  jetVarianceOuterHitPosInPixel  = getJetVariance(outerHitPosInPixel, false);
  // distributions outside the pixel layers
  jetMedianInnerHitPosOutPixel   = getJetMedian(innerHitPosOutPixel, false);
  jetMedianOuterHitPosOutPixel   = getJetMedian(outerHitPosOutPixel, false);
  jetMeanInnerHitPosOutPixel	 = getJetMean(innerHitPosOutPixel, false);
  jetMeanOuterHitPosOutPixel	 = getJetMean(outerHitPosOutPixel, false);
  jetVarianceInnerHitPosOutPixel = getJetVariance(innerHitPosOutPixel, false);
  jetVarianceOuterHitPosOutPixel = getJetVariance(outerHitPosOutPixel, false);
  // fraction valid hits
  jetMedianTrackValidHitFrac     = getJetMedian(fractionValidHits, false);
  jetMeanTrackValidHitFrac	 = getJetMean(fractionValidHits,false);
  jetVarianceTrackValidHitFrac   = getJetVariance(fractionValidHits,false);  
}


// compute variables related to the track angles and the calo jet 
void DisplacedJet::addTrackAngles(const reco::TrackCollection tracks, const edm::EventSetup& iSetup) {
  reco::TrackCollection::const_iterator tIter = tracks.begin();
  
  //  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> jetP4 = jet.physicsP4();  
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> jetDetP4 = jet.detectorP4();  
  float jdpx = jetDetP4.px(), jdpy = jetDetP4.py(), jdpz = jetDetP4.pz();
  
  // jet vectors 2d and 3d
  TVector3 jetv3(jet.px(), jet.py(), jet.pz());
  TVector3 jetv2(jet.px(), jet.py(), 0);
  TVector3 jetv3_xyz(jet.px(), jet.py(), jet.pz());
  TVector3 jetv2_xy(jet.px(), jet.py(), 0);
  // unit jet vectors
  TVector3 jetv3_unit = jetv3_xyz.Unit();
  TVector3 jetv2_unit = jetv2_xy.Unit();

  // detector jet vectors 2d and 3d
  TVector3 jetdv3(jdpx, jdpy, jdpz);
  TVector3 jetdv2(jdpx, jdpy, 0);
  // detector unit jet vectors
  TVector3 jetdv3_unit = jetdv3.Unit();
  TVector3 jetdv2_unit = jetdv2.Unit();

  //  std::cout << "Jet px = " << jet.px() << " py = " << jet.py() << std::endl;
  // pt weighted
  ptSumCosTheta2D	  = 0, ptSumCosTheta3D = 0;
  ptSumCosThetaDet2D	  = 0, ptSumCosThetaDet3D = 0;
  // aboslute sum
  sumCosTheta2D	          = 0, sumCosTheta3D = 0;
  sumCosThetaDet2D	  = 0, sumCosThetaDet3D = 0;   

  TVector3 trackMomSumV3(0,0,0);
  TVector3 trackMomSumV2(0,0,0);
  TVector3 ipPosSumV3(0,0,0);
  TVector3 ipPosSumV2(0,0,0);

  // find the angle between the jet momentum and track momentum
  float sumTrackPtValid = 0;
  // Create the transient track builder
  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  for(; tIter != tracks.end(); ++tIter) {
    
    // build the track momentum sum for the jet
    TVector3 trackMom3d(tIter->px(), tIter->py(), tIter->pz());
    TVector3 trackMom2d(tIter->px(), tIter->py(), 0);
    trackMomSumV3 += trackMom3d;
    trackMomSumV2 += trackMom2d;
    
    // uses recotrack debug to approximate hit positions
    static GetTrackTrajInfo getTrackTrajInfo;
    const reco::Track & const_track = *tIter;
    std::vector<GetTrackTrajInfo::Result> trajInfo = getTrackTrajInfo.analyze(iSetup, const_track);
    if (trajInfo[0].valid && trajInfo.back().valid) {

      // get the inner and outer hit of the track for valid hits
      const TrajectoryStateOnSurface&   tsosInnerHit = trajInfo[0].detTSOS;
      const TrajectoryStateOnSurface&   tsosOuterHit = trajInfo.back().detTSOS;
      // parse the position from the tsos
      const GlobalPoint &               innerPos     = tsosInnerHit.globalPosition();
      const GlobalPoint &               outerPos     = tsosOuterHit.globalPosition();
      // const GlobalVector &              innerMom     = tsosInnerHit.globalMomentum();
      // const GlobalVector &              outerMom     = tsosOuterHit.globalMomentum();

      // build the transient track for manipulations
      reco::TransientTrack transientTrack = builder->build(const_track);
      // create the extraoplator using the transient track
      //      TransverseImpactPointExtrapolator extrapolator(transientTrack.field());
      // extrapolate to the closest point to the primary vertex
      //      TrajectoryStateOnSurface tsosClosestToPV2D = extrapolator.extrapolate(transientTrack.impactPointState(), RecoVertex::convertPos(selPV.position()));
      TrajectoryStateClosestToPoint tscpPV = transientTrack.trajectoryStateClosestToPoint(RecoVertex::convertPos(selPV.position()));
      // get the position from the tsos
      //      const GlobalPoint & closestToPV = tsosClosestToPV2D.globalPosition();
      const GlobalPoint & closestToPV = tscpPV.position();
      const GlobalVector & closestMom = tscpPV.momentum();

      TVector3 innerMomv3( closestMom.x(), closestMom.y(), closestMom.z());
      TVector3 innerMomv2( closestMom.x(), closestMom.y(), 0);

      
      TVector3 outerv3( outerPos.x() - selPV.x(), outerPos.y() - selPV.y(), outerPos.z() - selPV.z());
      TVector3 outerv2( outerPos.x() - selPV.x(), outerPos.y() - selPV.y(), 0);

      TVector3 ipv3( closestToPV.x() - selPV.x(), closestToPV.y() - selPV.y(), closestToPV.z() - selPV.z());
      TVector3 ipv2( closestToPV.x() - selPV.x(), closestToPV.y() - selPV.y(), 0);

      //      std::cout << "Closest to PV x = " << closestToPV.x() - selPV.x() << " y = " << closestToPV.y() - selPV.y() << " z = " <<closestToPV.y() - selPV.y() << std::endl;

      ipPosSumV3 += ipv3;
      ipPosSumV2 += ipv2;

      TVector3 innerv3( innerPos.x(), innerPos.y(), innerPos.z());
      TVector3 innerv2( innerPos.x(), innerPos.y(), 0);

      float pt = tIter->pt(); 

      /* float dx = outerPos.x() - innerPos.x(); //tIter->outerX() - tIter->innerPosition().X(); */
      /* float dy = outerPos.y() - innerPos.y(); //tIter->outerY() - tIter->innerPosition().Y();  */
      /* float dz = outerPos.z() - innerPos.z();//1 - tIter->vz(); //tIter->outerZ() - tIter->innerPosition().Z(); */
      /* 

      /* std::cout << "TRACK INNER " << innerPos.x() << " " << innerPos.y() << " " << innerPos.z() << std::endl; */
      /* std::cout << "TRACK OUTER " << outerPos.x() << " " << outerPos.y() << " " << outerPos.z() << std::endl; */
      /* std::cout << "TRACK Di " << dx << " " << dy << " " << dz << std::endl; */

      /* // track vectors */
      /* TVector3	trkv3(dx, dy, dz); */
      /* TVector2	trkv2(dx, dy); */
      /* // track unit vectors */
      /* TVector3	trkv3_unit = trkv3.Unit(); */
      /* TVector2	trkv2_unit = trkv2.Unit(); */

      /* std::cout << "TRACK UNIT  3D " << trkv3_unit.X() << " " << trkv3_unit.Y() << trkv3_unit.Z() << std::endl; */
      /* std::cout << "TRACK UNIT  2D " << trkv2_unit.X() << " " << trkv2_unit.Y() << std::endl; */
    
      // track cosine with jet unit
      float   cosTheta3D	  = (jetv3.Unit().Cross((innerMomv3).Unit())).Mag();//outerv3.Unit().Cross((outerv3 - ipv3).Unit()).Mag(); //trkv3_unit * jetv3_unit;
      float   cosTheta2D	  = (jetv2.Unit().Cross((innerMomv2).Unit())).Mag();//outerv2.Unit().Cross((outerv2 - ipv2).Unit()).Mag(); //trkv2_unit * jetv2_unit;
      //      cosTheta3D = cosTheta3D ? log(cosTheta3D) : -20;
      //      cosTheta2D = cosTheta2D ? log(cosTheta2D) : -20;

      /* std::cout << "COSTHETA3D " << cosTheta3D << std::endl; */
      /* std::cout << "COSTHETA2D " << cosTheta2D << std::endl; */
      /* std::cout << " -----------------" << std::endl; */

      // track cosine with jet unit detecotr
      float   cosThetaDet3D = 0; //trkv3_unit * jetdv3_unit;
      float   cosThetaDet2D = 0; //trkv2_unit * jetdv2_unit;

      cosTheta2DVector.push_back(cosTheta2D);
      cosTheta3DVector.push_back(cosTheta3D);
      cosThetaDet2DVector.push_back(cosThetaDet2D);
      cosThetaDet3DVector.push_back(cosThetaDet3D);
    
      // increment sumpt
      sumTrackPtValid	 += pt;

      // pt weighted
      ptSumCosTheta2D    += cosTheta2D *    pt;
      ptSumCosTheta3D    += cosTheta3D *    pt;
      ptSumCosThetaDet2D += cosThetaDet2D * pt;
      ptSumCosThetaDet3D += cosThetaDet3D * pt;
      // not pt weighted
      sumCosTheta2D      += cosTheta2D;
      sumCosTheta3D      += cosTheta3D;
      sumCosThetaDet2D   += cosThetaDet2D;
      sumCosThetaDet3D   += cosThetaDet3D;
    }

    // mean
    meanCosTheta2D	  = getJetMean(cosTheta2DVector, true);
    meanCosTheta3D	  = getJetMean(cosTheta3DVector, true);
    meanCosThetaDet2D	  = getJetMean(cosThetaDet2DVector, true);
    meanCosThetaDet3D	  = getJetMean(cosThetaDet3DVector, true);
    // variance
    varianceCosTheta2D	  = getJetVariance(cosTheta2DVector, true);
    varianceCosTheta3D	  = getJetVariance(cosTheta3DVector, true);
    varianceCosThetaDet2D = getJetVariance(cosThetaDet2DVector, true);
    varianceCosThetaDet3D = getJetVariance(cosThetaDet3DVector, true);
    // median
    medianCosTheta2D	  = getJetMedian(cosTheta2DVector, true);
    medianCosTheta3D	  = getJetMedian(cosTheta3DVector, true);
    medianCosThetaDet2D	  = getJetMedian(cosThetaDet2DVector, true);
    medianCosThetaDet3D	  = getJetMedian(cosThetaDet3DVector, true);    
  }

  // track sum
  trackSumMomCosTheta2D	 = jetv2_unit * trackMomSumV2.Unit();
  trackSumMomCosTheta3D	 = jetv3_unit * trackMomSumV3.Unit();
  trackSumMomMag2D	 = trackMomSumV2.Mag();
  trackSumMomMag3D	 = trackMomSumV3.Mag();
  // vector sum of IP
  //std::cout << "Sin(theta,jet,ipVectorSum) " << (jetv2_unit.Cross(ipPosSumV2.Unit())).Mag() << " |ipVectorSum| " << ipPosSumV2.Mag() << std::endl;

  ipPosSumMag3D		 = (jetv3_unit.Cross(ipPosSumV3.Unit())).Mag();
  ipPosSumMag2D		 = (jetv2_unit.Cross(ipPosSumV2.Unit())).Mag();
  // normalize by the sum pt
  ptSumCosTheta2D	/= sumTrackPtValid ? sumTrackPtValid : 1;
  ptSumCosTheta3D	/= sumTrackPtValid ? sumTrackPtValid : 1;
  ptSumCosThetaDet2D	/= sumTrackPtValid ? sumTrackPtValid : 1;
  ptSumCosThetaDet3D	/= sumTrackPtValid ? sumTrackPtValid : 1;

}

void DisplacedJet::addSVTagInfo(const reco::SecondaryVertexTagInfo& svTagInfo) {
  if (debug > 1) std::cout << "[DEBUG 1] Adding SV Tag Info To Jet  " << std::endl;
  // number of vertices reconstructed
  svNVertex = svTagInfo.nVertices();
  if (svNVertex == 0) return;

  // track the best vertex reconstructed
  int svi = 0;
  int mostTracks = 0;
  int tieBreaker = 0;

  if (debug > 1) std::cout << "[DEBUG 1] Choosing best SV Vertex  " << std::endl;
  // loop over all the reconstructed vertices and pick 1
  for(int vv = 0; vv < svNVertex; vv++) {
    reco::Vertex vtx = svTagInfo.secondaryVertex(vv);   
    float pt = vtx.p4().pt();
    int nTracksSV = vtx.nTracks();

    // take the vertex with the most tracks, tie breaker is the sum pt of vertex
    if ( (nTracksSV > mostTracks) || (nTracksSV == mostTracks && pt > tieBreaker) ){
      mostTracks = nTracksSV;
      tieBreaker = pt;
      svi = vv;
    }     
  }
  
  // pick the selected vertex
  const reco::Vertex & selVertex = svTagInfo.secondaryVertex(svi);

  // set the global vertex to the selected vertex
  selSV = selVertex;

  if (debug > 2) std::cout << "[DEBUG 2] Filling SV Vertex Information into DJet  " << std::endl;
  // with the one vertex, store the quantities
  // position
  svX = selVertex.x();
  svY = selVertex.y();
  svZ = selVertex.z();

  svNTracks = selSV.nTracks();

  svChi2 = selVertex.chi2();
  svNDof = selVertex.ndof();
  
  svXError = selVertex.xError();
  svYError = selVertex.yError();
  svZError = selVertex.zError();

  float pvxE = selPV.xError(), pvyE = selPV.yError(), pvzE = selPV.zError();
  float xE   = std::sqrt(svXError * svXError + pvxE * pvxE), 
    yE = std::sqrt(svYError * svYError + pvyE * pvyE), 
    zE = std::sqrt(svZError * svZError + pvzE * pvzE);
  float dx = svX - selPV.x() , dy = svY - selPV.y(), dz = svZ - selPV.z();

  svMass    = selVertex.p4().mass();    
  svLxySig  = std::sqrt( dx * dx + dy * dy ) / std::sqrt(xE * xE + yE * yE);
  svLxyzSig = std::sqrt( dx * dx + dy * dy + dz * dz) / std::sqrt(xE * xE + yE * yE + zE * zE);
  svLxy     = std::sqrt( dx * dx + dy * dy );
  svLxyz    = std::sqrt( dx * dx + dy * dy + dz * dz );  
  if (debug > 2) std::cout << "[DEBUG 2] Filling SV Vertex Information into DJet -- Complete  " << std::endl;
}

void DisplacedJet::addIPTagInfo(const reco::TrackIPTagInfo & ipTagInfo) {
  if (debug > 2) std::cout << "[DEBUG 2] Adding SV IP Info into DJet  " << std::endl;
  

  // pull the impact parameter data
  lifetimeIPData = ipTagInfo.impactParameterData();
  
  // loop over each tracks ip information
  std::vector<reco::btag::TrackIPData>::const_iterator ipIter = lifetimeIPData.begin();
  for(; ipIter != lifetimeIPData.end(); ++ipIter)  {    

    if (debug > 3) std::cout << "[DEBUG] Filing IP INFO  " << std::endl;
    // ip and jet distance information
    float ip3d = ipIter->ip3d.value(), ip3ds = ipIter->ip3d.significance();
    float ip2d = ipIter->ip2d.value(), ip2ds = ipIter->ip2d.significance();
    float jetAxisDist = ipIter->distanceToJetAxis.value(), jetAxisDistSig = ipIter->distanceToJetAxis.significance();

    // fill the vectors
    ip3dVector.push_back(ip3d);
    ip2dVector.push_back(ip2d);
    ip3dsVector.push_back(ip3ds);
    ip2dsVector.push_back(ip2ds);

    // log quantities
    float   ipSigLog2D = (ip3d ? log(fabs(ip2ds)) : 0);
    float   ipSigLog3D = (ip3ds? log(fabs(ip3ds)) : 0);    
    float   ipLog2D    = (ip2d ? log(fabs(ip2d)) : 0);
    float   ipLog3D    = (ip3d ? log(fabs(ip3d)) : 0);

    ipLog2dsVector.push_back(ipSigLog2D);
    ipLog3dsVector.push_back(ipSigLog3D);
    ipLog3dVector.push_back(ipLog3D);
    ipLog2dVector.push_back(ipLog2D);

    jetAxisDistVector.push_back(jetAxisDist);
    jetAxisDistSigVector.push_back(jetAxisDistSig);

    // ip log sums 3d and 2d
    if (debug > 4) std::cout << "[DEBUG 4] ipsiglogsum2d: " << ipSigLogSum2D << std::endl;
    if (debug > 4) std::cout << "[DEBUG 4] jetdist: " << jetAxisDist << std::endl;
    if (debug > 4) std::cout << "[DEBUG 4] jetdistsig: " << jetAxisDistSig << std::endl;
    ipSigLogSum2D    += (ip3d ? log(fabs(ip2ds)) : 0);
    ipSigLogSum3D    += (ip3ds? log(fabs(ip3ds)) : 0);    
    ipLogSum2D	     += (ip2d ? log(fabs(ip2d)) : 0);
    ipLogSum3D	     += (ip3d ? log(fabs(ip3d)) : 0);
    jetDistSigLogSum += fabs(jetAxisDistSig);
    jetDistLogSum    += fabs(jetAxisDist);
  }

  // loop over associated tracks
  reco::TrackRefVector selectedTracks = ipTagInfo.selectedTracks();   
  reco::TrackRefVector::const_iterator trackIter = selectedTracks.begin();
  for(; trackIter != selectedTracks.end(); ++trackIter) {    
    trackAlgo.push_back((*trackIter)->algo());
  }

  if (debug > 2) std::cout << "[DEBUG 2] Deriving Distributional Quantities of IP Info  " << std::endl;  
  // distributional quantities 
  // mean
  meanIPSig2D	   = getJetMean(ip2dsVector, false);
  meanIPSig3D	   = getJetMean(ip3dsVector, false);
  meanIP2D	   = getJetMean(ip2dVector, false);
  meanIP3D	   = getJetMean(ip3dVector, false);
  meanJetDist	   = getJetMean(jetAxisDistVector, false);
  meanJetDistSig   = getJetMean(jetAxisDistSigVector, false);
  meanIPLogSig2D   = getJetMean(ipLog2dsVector, true);	//signed value matters
  meanIPLogSig3D   = getJetMean(ipLog3dsVector, true);	//signed value matters
  meanIPLog2D      = getJetMean(ipLog2dVector, true);	//signed value matters
  meanIPLog3D      = getJetMean(ipLog3dVector, true);	//signed value matters  

  // median
  medianIPSig2D	   = getJetMedian(ip2dsVector, true);
  medianIPSig3D	   = getJetMedian(ip3dsVector, true);
  medianIP2D	   = getJetMedian(ip2dVector, false);
  medianIP3D	   = getJetMedian(ip3dVector, false);
  medianJetDist	   = getJetMedian(jetAxisDistVector, false);
  medianJetDistSig = getJetMedian(jetAxisDistSigVector, false);
  medianIPLogSig2D = getJetMedian(ipLog2dsVector, true);	//signed value matters
  medianIPLogSig3D = getJetMedian(ipLog3dsVector, true);	//signed value matters
  medianIPLog2D    = getJetMedian(ipLog2dVector, true);	//signed value matters
  medianIPLog3D    = getJetMedian(ipLog3dVector, true);	//signed value matters

  // variance
  varianceIPSig2D    = getJetVariance(ip2dsVector, false);
  varianceIPSig3D    = getJetVariance(ip3dsVector, false);
  varianceIP2D	     = getJetVariance(ip2dVector, false);
  varianceIP3D	     = getJetVariance(ip3dVector, false);
  varianceJetDist    = getJetVariance(jetAxisDistVector, false);
  varianceJetDistSig = getJetVariance(jetAxisDistSigVector, false);
  varianceIPLogSig2D = getJetVariance(ipLog2dsVector, true);	//signed value matters
  varianceIPLogSig3D = getJetVariance(ipLog3dsVector, true);	//signed value matters
  varianceIPLog2D    = getJetVariance(ipLog2dVector, true);	//signed value matters
  varianceIPLog3D    = getJetVariance(ipLog3dVector, true);	//signed value matters
}

// no IVF found
std::vector<bool> DisplacedJet::passNoVtxTag(const std::vector<float> thres){
  std::vector<bool> passResults;
  std::vector<float>::const_iterator thresIter = thres.begin();
  for(; thresIter != thres.end(); ++thresIter) {
    bool didPass = false;
    if(selIVFIsPV && nTracks == 0 && caloHadEnergyFrac > 0.95) didPass = true;
    passResults.push_back(didPass);
  }
  noVertexTagsVector = passResults;
  return passResults;
}

std::vector<bool> DisplacedJet::passShortTag(const std::vector<float> thres, float min, float max){
  std::vector<bool> passResults;
  std::vector<float>::const_iterator thresIter = thres.begin();
  for(; thresIter != thres.end(); ++thresIter) {
    bool didPass = false;
    if(selIVFIsPV && nTracks > 0 && medianIPLogSig2D > *thresIter  && (alphaMax / sumTrackPt) < .05) didPass = true;
    passResults.push_back(didPass);
  }
  shortTagsVector = passResults;
  return passResults;  
}

std::vector<bool> DisplacedJet::passMediumTag(const std::vector<float> thres, float min, float max){
  std::vector<bool> passResults;
  std::vector<float>::const_iterator thresIter = thres.begin();
  for(; thresIter != thres.end(); ++thresIter) {
    bool didPass = false;
    if(!selIVFIsPV && medianIPLogSig2D > *thresIter && ivfLxy > min && ivfLxy < max && (alphaMax / sumTrackPt) < .05) didPass = true;
    passResults.push_back(didPass);
  }
  mediumTagsVector = passResults;
  return passResults;  
}

std::vector<bool> DisplacedJet::passLongTag(const std::vector<float> thres, float min, float max){
  std::vector<bool> passResults;
  std::vector<float>::const_iterator thresIter = thres.begin();
  for(; thresIter != thres.end(); ++thresIter) {
    bool didPass = false;
    if(!selIVFIsPV && medianIPLogSig2D > *thresIter && ivfLxy > min && ivfLxy < max && (alphaMax / sumTrackPt) < .05) didPass = true;
    passResults.push_back(didPass);
  }
  longTagsVector = passResults;
  return passResults;  
}


float DisplacedJet::getJetMedian(const std::vector<float>& values, bool is_signed) {   
  if (values.size() == 0) return -999;
  
  int size = values.size();

  if(size == 0) return 0; //no tracks

  float * sorted = new float[size];
  for (int ii = 0; ii < size; ++ii) {
    sorted[ii] = is_signed ? values[ii] : fabs(values[ii]);
  }

  // sort the array
  for (int i = size - 1; i > 0; --i) {
    for (int j = 0; j < i; ++j) {
      if (sorted[j] > sorted[j+1]) {
	float dTemp = sorted[j];
	sorted[j] = sorted[j+1];
	sorted[j+1] = dTemp;
      }
    }
  }

  // Middle or average of middle values in the sorted array.
  float  dMedian = 0.0;
  if ((size % 2) == 0) {
    dMedian = (sorted[size/2] + sorted[(size/2) - 1])/2.0;
  } else {
    dMedian = sorted[size/2];
  }

  delete [] sorted;

  //safety check
  if(!is_signed ) { 
    if (dMedian < 0 ) std::cout << "ASSERT FAIL: "  << dMedian << std::endl;
    assert( dMedian >= 0 ); 
  }

  if (debug > 4) std::cout << "[DEBUG] jet median:  " << dMedian << std::endl;
  return dMedian;
}

float DisplacedJet::getJetMean(const std::vector<float> & values, bool is_signed) {
  if (values.size() == 0) return 0.0;

  float sum = 0;
  std::vector<float>::const_iterator val = values.begin();
  for (; val != values.end(); ++val) sum += is_signed ? *val : fabs(*val);  

  float mean = sum / float(values.size());
  if (debug > 4) std::cout << "[DEBUG] jet mean:  " << mean << std::endl;
  return mean;
}

float DisplacedJet::getJetVariance(const std::vector<float>& values, bool is_signed) {
  if (values.size() == 0 || values.size() == 1) return 0.0;

  float sum = 0;
  float mean = getJetMean(values, is_signed);

  std::vector<float>::const_iterator val = values.begin();
  for (; val != values.end(); ++val) {
    float temp =  (*val - mean) * (*val - mean);
    /* if (debug > 5) std::cout << "[DEBUG 5] value:" << *val << std::endl; */
    /* if (debug > 5) std::cout << "[DEBUG 5] mean: " << mean << std::endl; */
    /* if (debug > 5) std::cout << "[DEBUG 5] variance element: " << temp << std::endl; */
    sum += temp;
  }

  //safety check
  if(!is_signed ) { assert( sum >= 0 ); }
  float variance = std::sqrt(sum / float((values.size() - 1.0)));

  if (debug > 4) std::cout << "[DEBUG] jet variance:  " << variance << std::endl;
  return variance;
}

// calculate the jet variable alpha: the ratio of (vertex tracks sum pt matching the jet)
// divided by ( general tracks sum pt matching the jet)
void DisplacedJet::calcJetAlpha(const reco::TrackCollection& tracks, const reco::VertexCollection& primaryVertices) { 
  // calculate the sum track pt based on the tracks matched  

  // Take the scalar sum pt of tracks from the primary vertices relative to the total pt matched to the jets
  float sumJetPt	       = 0;
  float sumJetPtMax	       = 0;
  float sumJetPt_temp	       = 0;
  // highest vertex sum pt^2 
  float leadingVertexPtSq      = 0;
  float leadingVertexPtSq_temp = 0;

  // loop over vertices in the event (use beam spot constraint vertices
  reco::VertexCollection::const_iterator vtxIter = primaryVertices.begin();
  for(; vtxIter != primaryVertices.end(); ++vtxIter ) {
    // these are the sums for this specific vertex
    sumJetPt_temp = 0; // sum for matching within the jet
    leadingVertexPtSq_temp  = 0; // sum for the whole vertex

    std::vector<reco::TrackBaseRef>::const_iterator tIter = vtxIter->tracks_begin();
    // loop over tracks in the vertex
    for(; tIter != vtxIter->tracks_end(); ++tIter) {
      float pt = (*tIter)->pt();
      
      if((*tIter)->pt() < 1.0) continue; //apply a cut on the track pt 
      leadingVertexPtSq_temp += pt * pt; // sort by highest pt^2 

      // check for matching to the jet
      float dr = reco::deltaR((*tIter)->eta(), (*tIter)->phi(), caloEta, caloPhi);      
      if (dr < 0.4 ) sumJetPt_temp += (*tIter)->pt() ;
    }

    // pick the vertex with the highest sum pt^2
    if(leadingVertexPtSq_temp > leadingVertexPtSq) {
      leadingVertexPtSq  = leadingVertexPtSq_temp;
      sumJetPt = sumJetPt_temp;
    } // loop over tracks from vertex 

    // also keep the highest possible sum
    sumJetPtMax = std::max(sumJetPtMax, sumJetPt_temp);

  } // loop over vertices

  alpha	   = sumJetPt;
  alphaMax = sumJetPtMax;
}
