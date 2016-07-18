typedef std::vector<DisplacedTrack> DisplacedTrackCollection;
typedef std::vector<Displaced2TrackVertex> DisplacedV0Collection;

class DisplacedJet {
 public:
 DisplacedJet(const reco::CaloJet & jet_, 
	      const reco::Vertex & primaryVertex, 
	      const bool & isMC_, const int& jetID_, 
	      const edm::EventSetup& iSetup_, const int & debug_)
   :  debug(debug_), iSetup(iSetup_), jet(jet_), isMC(isMC_), jetID(jetID_), selPV(primaryVertex){

    if (debug > 2) std::cout << "[DEBUG] Creating a Displaced Jet with  " << jet_.pt() << " " << jet_.eta() << " " << jet_.phi() << std::endl;

    // track association variables
    nTracks		 = 0;
    nTracksPrompt	 = 0;
    nTracksDisp		 = 0;
    // regional tracking counts
    nTracksRegPrompt	 = -1;
    nTracksRegDisp	 = -1;
    nTracksRegPromptUp	 = -1;
    nTracksRegPromptDn	 = -1;
    nTracksRegDispUp	 = -1;
    nTracksRegDispDn	 = -1;
    nTracksRegDisp	 = -1;
    nTracksReg0124	 = -1;
    nTracksReg012	 = -1;
    nTracksReg4		 = -1;
    // pass the HLT requirements?
    passHLTPrompt	 = false;
    passHLTDisp		 = false;
    // up and down varied systematics
    passHLTPromptUp	 = false;
    passHLTDispUp	 = false;
    passHLTPromptDn	 = false;
    passHLTDispDn	 = false;
    passHLTPromptAndDisp = false;

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

    // jet vertex fraction variable
    alpha = 0, alphaMax = 0;
    
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
    svAngle3D = 0;
    svAngle2D = 0;
    svPt = 0;
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
    // NUCLEAR INTERACTION RELATED
    jetOneTrackNuclearCount	   = 0;
    jetTwoTrackNuclearCount	   = 0;
    jetTwoTrackInnerHitFake        = 0;
    jetVertexNearBPIX1		   = 0;
    jetVertexNearBPIX2		   = 0;
    jetVertexNearBPIX3		   = 0;
    jetVertexNearBPIX		   = 0;
    jetTightNuclear		   = 0;
    jetLooseNuclear		   = 0;
    jetNV0HitBehindVertex	   = 0;
    jetNV0NoHitBehindVertex	   = 0;
    jetNV0KShort		   = 0;
    jetNV0Lambda		   = 0;
    // GRAPH RELATED
    jetV0HIndex			   = -1;
    // CLUSTER RELATED
    jetV0ClusterSize		   = 0;
    jetV0ClusterLxy		   = 0;
    jetV0ClusterLxySig		   = 0;
    jetV0ClusterLxyzSig		   = 0;
    jetV0ClusterLxyz		   = 0;
    jetV0ClusterX		   = 0;
    jetV0ClusterY		   = 0;
    jetV0ClusterZ		   = 0;
    jetV0ClusterChi2		   = 0;
    jetV0ClusterIntercept	   = 0;
    jetV0ClusterAngle		   = 0;
    jetV0ClusterAngleMom           = 0;
    jetV0ClusterNTracks		   = 0;
    // N JET CLUSTERS
    jetV0NJetClusterSize	   = 0;
    jetV0NJetClusterLxy		   = 0;
    jetV0NJetClusterLxySig	   = 0;
    jetV0NJetClusterLxyzSig	   = 0;
    jetV0NJetClusterLxyz	   = 0;
    jetV0NJetClusterX		   = 0;
    jetV0NJetClusterY		   = 0;
    jetV0NJetClusterZ		   = 0;
    jetV0NJetClusterChi2	   = -1;
    jetV0NJetClusterIntercept	   = 0;
    jetV0NJetClusterAngle	   = 0;
    jetV0NJetClusterAngleMom	   = 0;
    jetV0NJetClusterNTracks	   = 0;
  } // class declaration end

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
  void addRegionalTracks(const reco::TrackRefVector&, const int & collectionID, const float & smear_2dip, const float & smear_2dipsig);
  void addVertexTrackInfo(const reco::TrackRefVector&);
  void addIVFCollection(const reco::VertexCollection&, const float& compatibilityScore);
  void addIPTagInfo(const reco::TrackIPTagInfo&);
  void addSVTagInfo(const reco::SecondaryVertexTagInfo&);
  
  // angular information relative to ejt
  void addTrackAngles(const DisplacedTrackCollection & tracks);
  // information about track inner and outer hits
  void addHitInfo(const reco::TrackCollection tracks);
  // v0 information for pairs of tracks
  void addV0Info(const reco::TrackRefVector tracks);

  /////////////////CLIQUE RELATED///////////////

  typedef std::vector<std::vector<int> > vertexGraph; 
  // clique vertexing 
  vertexGraph buildInputGraph(const std::vector<std::pair<int,int> > & graph_edges, const std::vector<int> & uniq_tracks, const bool & print);
  void findVertexCliques();
  int findRefTrack(const reco::TrackRef ref, const reco::TrackRefVector vector);
  int calcHIndex(const vertexGraph & graph);
  void calcClusterSize(const DisplacedV0Collection& vertices, const float & errorWindow);
  void calcNJetClusterSize(const DisplacedV0Collection& vertices, std::vector<DisplacedJet>& djets, const float & errorWindow);  

  //////////////////////////////////////////////

  // generator matching
  bool doGenCaloJetMatching(const float& ptMatch, const float& dRMatch, const reco::GenParticleCollection& genParticles);
  bool doGenVertexJetMatching(const float& metricThreshold, const reco::GenParticleCollection& genParticles);
  void doGenVertexID(const float& metricThreshold, const reco::GenParticleCollection& genParticles);
  float genMatchMetric(const reco::GenParticle & particle, const reco::Vertex& vertex);
  
  // jet info extraction
  DisplacedTrackCollection getDisplacedTracks() { return displacedTracks; }
  reco::TrackCollection getCaloMatchedTracks() { return caloMatchedTracks; }
  reco::TrackCollection getVertexMatchedTracks() { return vertexMatchedTracks; }
  reco::TrackRefVector  getVertexMatchedTrackRefs() { return vertexMatchedTrackRefs; }
  reco::Vertex          getIVFVertexSelected() { return selIVF; }
  reco::Vertex          getSVVertex() { return selSV; }
  
  // jet distribution calculator
  float getJetMedian(const std::vector<float>&, bool);
  float getJetMean(const std::vector<float>&, bool);
  float getJetVariance(const std::vector<float>&, bool);

  //alpha calculation
  void calcJetAlpha(const reco::TrackCollection&, const reco::VertexCollection&);

  // initialized varaibles (order matters)
  const int		    debug;
  const edm::EventSetup&    iSetup;
  const reco::CaloJet	    jet;
  const bool		    isMC;
  const int		    jetID;
  const reco::Vertex	    selPV;  

  //////////////CALO INFORMATION////////////
  // track association variables
  int nTracks;
  int nTracksPrompt;
  int nTracksDisp;
  int nTracksRegPrompt;
  int nTracksRegDisp;
  int nTracksRegPromptUp, nTracksRegPromptDn;
  int nTracksRegDispUp, nTracksRegDispDn;
  float smear_2dipsig, smear_2dip;
  int nTracksReg0124,nTracksReg012, nTracksReg4;

  bool passHLTPrompt, passHLTDisp, passHLTPromptAndDisp;
  bool passHLTPromptUp, passHLTDispUp;
  bool passHLTPromptDn, passHLTDispDn;
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

  ////////////// NUCLEAR INTERACTION /////////

  int jetOneTrackNuclearCount, jetTwoTrackNuclearCount;
  int jetTwoTrackInnerHitFake;
  int jetVertexNearBPIX1, jetVertexNearBPIX2, jetVertexNearBPIX3, jetVertexNearBPIX;
  int jetTightNuclear, jetLooseNuclear;
  int jetNV0NoHitBehindVertex, jetNV0HitBehindVertex;
  int jetNV0KShort, jetNV0Lambda;
  int jetV0HIndex;

  // cluster related variables
  int   jetV0ClusterSize;
  float jetV0ClusterLxy;
  float jetV0ClusterLxySig;
  float jetV0ClusterLxyzSig;
  float jetV0ClusterLxyz;
  float jetV0ClusterX;
  float jetV0ClusterY;
  float jetV0ClusterZ;
  float jetV0ClusterChi2;
  float jetV0ClusterIntercept;
  float jetV0ClusterAngle;
  float jetV0ClusterAngleMom;
  int   jetV0ClusterNTracks;

  int   jetV0NJetClusterSize;
  float jetV0NJetClusterLxy;
  float jetV0NJetClusterLxySig;
  float jetV0NJetClusterLxyzSig;
  float jetV0NJetClusterLxyz;
  float jetV0NJetClusterX;
  float jetV0NJetClusterY;
  float jetV0NJetClusterZ;
  float jetV0NJetClusterChi2;
  float jetV0NJetClusterIntercept;
  float jetV0NJetClusterAngle;
  float jetV0NJetClusterAngleMom;
  int   jetV0NJetClusterNTracks;

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
  // kinematics
  float svPt;
  float svEta, svPhi;
  float svAngle2D, svAngle3D;
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
  
  // combinatorial vertices for the jet
  vertexGraph v0Graph;
  DisplacedV0Collection displacedV0Vector;
  DisplacedV0Collection displacedV0VectorCleaned;  
  std::vector<TransientVertex> transientV0Vector;
  std::vector<TransientVertex> transientV0VectorCleaned;  
  
  // related vertices
  bool selIVFIsPV;
  float selIVFIsPVScore;
  reco::Vertex selIVF;
  reco::Vertex selSV;

  // clusters
  DisplacedCluster *v0Cluster = NULL;
  DisplacedCluster *v0NJetCluster = NULL;

  std::vector<float> cosTheta2DVector, cosThetaDet2DVector, cosTheta3DVector, cosThetaDet3DVector; 
  std::vector<float> ip3dVector, ip3dsVector, ip2dVector, ip2dsVector;
  std::vector<float> trEtaVector, trPhiVector, trPtVector;

  // related regional track collections
  DisplacedTrackCollection displacedTracks;
  DisplacedTrackCollection regionalTracks0124;
  DisplacedTrackCollection regionalTracks012;
  DisplacedTrackCollection regionalTracks4;

 private: 

  static const int GEN_STATUS_CODE_MATCH = 23; 
  const float FAKE_HIGH_NUMBER = 999999999;
  // calo jet the displaced jet is built upon


  // related  track collections  
  reco::TrackCollection caloMatchedTracks; 
  reco::TrackCollection vertexMatchedTracks; 
  reco::TrackRefVector  vertexMatchedTrackRefs; 
  std::vector<reco::btag::TrackIPData> lifetimeIPData; 

  std::vector<int>   trackAlgo;
  std::vector<float> ipLog3dVector, ipLog3dsVector, ipLog2dVector, ipLog2dsVector;
  std::vector<float> jetAxisDistVector, jetAxisDistSigVector; 

  // cos related


  // simplify quadrature calculations
  float metric2D(float x, float y) {
    return std::sqrt(x*x + y*y);
  }
  float metric3D(float x, float y, float z) {
    return std::sqrt(x*x + y*y + z*z);
  }  
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
  return isInclusive(false) && isDispTrack(false);
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
    if (fabs(*ip2DIter) < 0.05 && pass_algo) nPTracks++;
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
  if (debug > 2) std::cout << "[DEBUG] Not Adding Calorimeter Matched Track Info  " << std::endl;
    /* reco::TrackRefVector::const_iterator trackIter = trackRefs.begin(); */
    /* for(; trackIter != trackRefs.end(); ++trackIter) {       */
    /*   caloMatchedTracks.push_back(**trackIter); */
    /* }     */
}

void DisplacedJet::addRegionalTracks(const reco::TrackRefVector & trackRefs, const int & collectionID, const float & smear_2dip, const float & smear_2dipsig) {	     
  if (debug > 2) std::cout << "[DEBUG] Adding Vertex Matched Track Info  " << std::endl;
  reco::TrackRefVector::const_iterator trackIter = trackRefs.begin();
  
  // loop ver all the track associations and buid the displaced track collecitons first
  for(; trackIter != trackRefs.end(); ++trackIter) {
    // build the displacedTrack
    DisplacedTrack dTrack(*trackIter, selPV, iSetup, debug);
    // based on the collection ID add the displaced track to the correct set of tracks
    if(collectionID == 0) regionalTracks0124.push_back(dTrack);
    else if(collectionID == 1) regionalTracks012.push_back(dTrack);
    else if(collectionID == 2) regionalTracks4.push_back(dTrack);
    else {
      std::cout << "INVALID REGIONAL TRACKING COLLECTION ID: " << collectionID << std::endl;
    }   
  } // end loop over references from jet track association

  // make sure this colleciton hasnt been merged before
  // count the number of prompt tracks from the regional tracking

  if(collectionID == 1) {
    int ntrack012 = regionalTracks012.size();
    if(debug > 3) std::cout << "n tracks 012: " << ntrack012 << std::endl;
    nTracksRegPrompt   = 0;
    nTracksReg012      = 0;
    nTracksRegPromptUp = 0;
    nTracksRegPromptDn = 0;
    
    for(int tt = 0; tt < ntrack012; ++tt) {
      const DisplacedTrack & track = regionalTracks012[tt];
      if(debug > 3) std::cout << "012: checking track pt: " << track.pt << std::endl;
      if (track.pt  < 1.0) continue;
      float ip2d = track.ip2d;
      nTracksReg012++;
      if(fabs(ip2d) < 0.1) nTracksRegPrompt++;
      if(fabs(ip2d)*(1+smear_2dip) < 0.1) nTracksRegPromptUp++;
      if(fabs(ip2d)*(1-smear_2dip) < 0.1) nTracksRegPromptDn++;
    }
  }// close collectoinID 1  
  // make sure this colleciton hasnt been merged before
  // count the number of displaced tracks from the regional tracking
  else if(collectionID == 0) {
    nTracksRegDisp   = 0;
    nTracksReg0124   = 0;
    nTracksRegDispUp = 0;
    nTracksRegDispDn = 0;

    int ntrack0124 = regionalTracks0124.size();
    if(debug > 3)  std::cout << "n tracks 0124: " << ntrack0124 << std::endl;
    for(int tt = 0; tt < ntrack0124; ++tt) {
      const DisplacedTrack & track = regionalTracks0124[tt];
      if(debug > 3) std::cout << "0124: checking track pt: " << track.pt << std::endl;
      if(track.pt < 1.0) continue;
      float ip2dSig = track.ip2dSig;
      float ip2d    = track.ip2d;
      nTracksReg0124++;
      if(fabs(ip2dSig) > 5.0 && fabs(ip2d) > 0.05) nTracksRegDisp++;
      if((fabs(ip2dSig)*(1+smear_2dipsig)) > 5.0 && (fabs(ip2d)*(1+smear_2dip)) > 0.05) nTracksRegDispUp++;
      if((fabs(ip2dSig)*(1-smear_2dipsig)) > 5.0 && (fabs(ip2d)*(1-smear_2dip)) > 0.05) nTracksRegDispDn++;
    }
  }   // close collectionID 0 if  
  // check this collection has not yet been counted
  // only cout tracks with pt > 1.0
  else if(collectionID == 2) {
    int ntrack4 = regionalTracks4.size();
    if(debug > 3) std::cout << "n tracks 4: " << ntrack4 << std::endl;
    nTracksReg4 = 0;
    for(int tt = 0; tt < ntrack4; ++tt) {
      const DisplacedTrack & track = regionalTracks4[tt];
      if(debug > 3) std::cout << "4: checking track pt: " << track.pt << std::endl;
      if (track.pt  < 1.0) continue;
      nTracksReg4++;
    }    
  }
  else {
    std::cout << "INVALID COLLECTION ID FOR REGIONAL TRACK MERGER...EXITING" << std::endl;
    exit(1);
  }

  // decided whether the jet passed the event or not
  if(nTracksRegPrompt <= 2 && collectionID == 1) passHLTPrompt = true; //can only test this after checking the 012 iters
  if(nTracksRegDisp >= 1 && collectionID == 0) passHLTDisp = true; //can only test this after checking the 0124 iters
  // up and down systematics
  // up varied ip and ip2dsig
  if(nTracksRegPromptUp <= 2 && collectionID == 1) passHLTPromptUp = true; //can only test this after checking the 012 iters
  if(nTracksRegPromptDn <= 2 && collectionID == 1) passHLTPromptDn = true; //can only test this after checking the 012 iters
  // down varied ip and ip2dsig
  if(nTracksRegDispUp >= 1 && collectionID == 0) passHLTDispUp = true; //can only test this after checking the 0124 iters
  if(nTracksRegDispDn >= 1 && collectionID == 0) passHLTDispDn = true; //can only test this after checking the 0124 iters

  passHLTPromptAndDisp			  = false; // this is meaningless as it stands because it can only be decided after all tracking collections have been checked

  if(debug > 3) std::cout << "colID: " << collectionID << "jet status of track counting 0124 " << nTracksReg0124 << " 012: " << nTracksReg012 << " 4: " << nTracksReg4 << std::endl;
  
}

// count the number of tracks based on the association at the vertex
void DisplacedJet::addVertexTrackInfo(const reco::TrackRefVector & trackRefs) {
  vertexMatchedTrackRefs = trackRefs;

  if (debug > 2) std::cout << "[DEBUG] Adding Vertex Matched Track Info  " << std::endl;
  reco::TrackRefVector::const_iterator trackIter = trackRefs.begin();
  nTracks	= 0;
  sumTrackPt	= 0;
  nTracksPrompt = 0;
  nTracksDisp	= 0;
  for(; trackIter != trackRefs.end(); ++trackIter) {
    // build the displacedTrack
    DisplacedTrack dTrack(*trackIter, selPV, iSetup, debug);
    displacedTracks.push_back(dTrack);

    // apply a pt cut for the vertexMatchTrack Collection
    float pt = (*trackIter)->pt();
    if(pt > 1.0) {
      nTracks++;
      sumTrackPt += pt;
      vertexMatchedTracks.push_back(**trackIter);

      // incre
      if(dTrack.ip2d < 0.1) nTracksPrompt++;
      if(dTrack.ip2d > 0.05 && fabs(dTrack.ip2dSig) > 5.0) nTracksDisp++; 
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

int DisplacedJet::findRefTrack(const reco::TrackRef ref, const reco::TrackRefVector vector) {
  reco::TrackRefVector::const_iterator iter = vector.begin();
  int ii = 0;
  bool found = false;
  for(; iter != vector.end(); ++iter, ++ii) {    
    bool    pt_match  = (**iter).pt() == (*ref).pt();
    bool    eta_match = (**iter).eta() == (*ref).eta();
    bool    phi_match = (**iter).phi() == (*ref).phi();
    
    if(pt_match && eta_match && phi_match)   return ii;
  }     

  assert(!found);
  return -999;
}

DisplacedJet::vertexGraph DisplacedJet::buildInputGraph(const std::vector<std::pair<int,int> > & graph_edges, const std::vector<int> & uniq_tracks_, const bool & print) {
  if (debug > 4) std::cout << "[DEBUG 4]  Building Input Graph for h index" << std::endl;  

  std::vector<int> uniq_tracks(uniq_tracks_.begin(), uniq_tracks_.end());
  int ntracks = uniq_tracks.size();

  // sort the vector
  std::sort(uniq_tracks.begin(), uniq_tracks.end());    

  assert(int(uniq_tracks.size()) == ntracks);

  // fill the rows with zeros
  std::vector<std::vector<int>> rows;
  for(int ii = 0; ii < ntracks; ii++ ) {
    std::vector<int> row;
    for(int jj = 0; jj < ntracks; jj++ ) row.push_back(0);     
    rows.push_back(row);
  }

  // fill the individual positions with 1 = filled or 0 un-matched
  std::vector<std::pair<int,int> >::const_iterator edgeIter = graph_edges.begin();
  for(; edgeIter != graph_edges.end(); ++edgeIter) {    
    std::vector<int>::iterator  pos_row = std::find(uniq_tracks.begin(), uniq_tracks.end(), edgeIter->first);
    std::vector<int>::iterator  pos_col = std::find(uniq_tracks.begin(), uniq_tracks.end(), edgeIter->second);

    //std::cout << "row " << edgeIter->first << " col " << edgeIter->second << std::endl;
    // get the row and col from the pointer
    int row = std::distance(uniq_tracks.begin(), pos_row);
    int col = std::distance(uniq_tracks.begin(), pos_col);

    assert(row < ntracks && col < ntracks);

    // encode the column and row into the graph
    rows[row][col] = 1;
    rows[col][row] = 1;        
  }

  // print the result
  if(print) {
    std::cout << " ntracks " << ntracks << std::endl;
    std::cout << "ordering of the uniq tracks: ";
    for(int ii = 0; ii < ntracks; ii++ ) {
      std::cout << uniq_tracks[ii] << " ";
    }  

    for(int ii = 0; ii < ntracks; ii++ ) {
      std::cout << "\n";
      for(int jj = 0; jj < ntracks; jj++ )  {      
	std::cout << rows[ii][jj] << " ";
      }
    }
    std::cout << std::endl;
  }

  return rows;
}

// build a graph of track fit associations based on cleaned transient vertices
void DisplacedJet::findVertexCliques() {
  if (debug > 4) std::cout << "[DEBUG 4]  Finding Vertex Cliques" << std::endl;  
  // store the two tracks of the vertex as a graph edge
  std::vector<std::pair<int,int> > graph_edges;
  std::vector<int> uniq_tracks;
  
  // loop over each vertex in the jet
  DisplacedV0Collection::const_iterator vtxIter = displacedV0VectorCleaned.begin();  
  for(; vtxIter != displacedV0VectorCleaned.end(); ++vtxIter) {

    reco::TrackRef	    track_ref1 = (vtxIter->track1).trackRef;
    reco::TrackRef	    track_ref2 = (vtxIter->track2).trackRef;

    // interpret the index in the array
    int refIndex1 = findRefTrack(track_ref1, vertexMatchedTrackRefs); 
    int refIndex2 = findRefTrack(track_ref2, vertexMatchedTrackRefs); 

    // add the edge to the graph
    graph_edges.push_back(std::pair<int, int>(refIndex1, refIndex2)) ;     

    // check for the tracks inside of the uniq tracks vector 
    std::vector<int>::iterator findTrack1 = std::find(uniq_tracks.begin(), uniq_tracks.end(), refIndex1);
    if (findTrack1 == uniq_tracks.end()) uniq_tracks.push_back(refIndex1);
    std::vector<int>::iterator findTrack2 = std::find(uniq_tracks.begin(), uniq_tracks.end(), refIndex2);
    if (findTrack2 == uniq_tracks.end()) uniq_tracks.push_back(refIndex2);              

  } // end for loop over cleaned transient vertices 

  std::vector<int>::const_iterator		    trackIter = uniq_tracks.begin();
  std::vector<std::pair<int, int> >::const_iterator edgeIter  = graph_edges.begin();

  // print the uniq tracks in the combinatorial vertices for the jet
  // and the edges in the vertex graph

  /* std::cout << "\n Jet Median 2DIP log significance: " << medianIPLogSig2D << std::endl; */
  /* std::cout << " UNIQ TRACKS IN JET "; */
  /* for(; trackIter != uniq_tracks.end(); ++trackIter) { */
  /*   std::cout << *trackIter << " "; */
  /* } */
  /* std::cout << "\n EDGES IN GRAPH "; */
  /* for(; edgeIter != graph_edges.end(); ++edgeIter) { */
  /*   std::cout << "("<< edgeIter->first << "," << edgeIter->second << ")"; */
  /* } */

  // build the graph of track pairs
  v0Graph = buildInputGraph(graph_edges, uniq_tracks, false); // dont print
  // 
  jetV0HIndex = calcHIndex(v0Graph);
}


void DisplacedJet::calcNJetClusterSize(const DisplacedV0Collection& vertices, std::vector<DisplacedJet>& djets, const float & errorWindow) {
  if (debug > 4) std::cout << "[DEBUG 4]  Finding N Jet Cluster Size" << std::endl;

  if ( vertices.size() == 0 ) return;

  DisplacedCluster * maxCluster = NULL;

  // use each vertex in "this" jet as a center for a cluster, then check all possible vertices
  // from the jet collection passed to the function
  const int nVtx = vertices.size();
  for(int center = 0; center < nVtx; ++center) {
    Displaced2TrackVertex centerVtx = vertices[center];
    if(!centerVtx.isValid) continue;

    // create a candidate cluster
    DisplacedCluster cluster_temp(centerVtx, selPV, iSetup, debug);
    
    // check all the vertices in "this jet"
    for(int neighbor = 0; neighbor < nVtx; ++neighbor) {
      if (center == neighbor) continue; // dont count the center itself

      Displaced2TrackVertex neighVtx = vertices[neighbor];
      if(!neighVtx.isValid) continue;            
      // check the distance from the center
      float dx	     = centerVtx.x - neighVtx.x, dy = centerVtx.y - neighVtx.y;	//, dz = cz - nz;
      float distance = metric2D(dx, dy);
      float sig	     = distance / metric2D(neighVtx.tot_xyE, centerVtx.tot_xyE);

      if(sig < errorWindow) cluster_temp.addVertex(neighVtx);
    }// end loop over neighbors in this jet

    // check all the vertices in the neighboring jets
    std::vector<DisplacedJet>::const_iterator djIter = djets.begin();
    for(; djIter != djets.end(); ++djIter) {
      // loop over all the vertices in the current neighbor jet
      int nVtx = djIter->displacedV0VectorCleaned.size();
      for(int neighbor = 0; neighbor < nVtx; ++neighbor) {	
	Displaced2TrackVertex neighVtx = djIter->displacedV0VectorCleaned[neighbor];
	if(!neighVtx.isValid) continue;            

	// check the distance from the center
	float dx	 = centerVtx.x - neighVtx.x, dy = centerVtx.y - neighVtx.y; //, dz = cz - nz;
	float distance   = metric2D(dx, dy);
	float sig	 = distance / metric2D(neighVtx.tot_xyE, centerVtx.tot_xyE);

	if(sig < errorWindow) cluster_temp.addVertex(neighVtx);
      } // loop over all the candidate neighbors in the jet
    }  // loop over other displaced jets in the event  

    // check if the resulting cluster has the highest multiplicity
    if(maxCluster == NULL) { 
      maxCluster = &cluster_temp;
      continue;
    }
    if(cluster_temp.nV0 > maxCluster->nV0) maxCluster = &cluster_temp;         
  } // end loop over center vertices in "this" jet
  
  // store quatities related to the jet cluster
  maxCluster->buildClusterQuantities();
  // set the pointer for the jet to the cluster
  v0NJetCluster = maxCluster;
  // fill the related quantities
  jetV0NJetClusterSize	    = maxCluster->nV0;
  jetV0NJetClusterLxy	    = maxCluster->center.lxy;
  jetV0NJetClusterLxySig    = maxCluster->center.lxySig;
  jetV0NJetClusterLxyzSig   = maxCluster->center.lxyzSig;
  jetV0NJetClusterLxyz	    = maxCluster->center.lxyz;
  jetV0NJetClusterX	    = maxCluster->center.x;
  jetV0NJetClusterY	    = maxCluster->center.y;
  jetV0NJetClusterZ	    = maxCluster->center.z;
  jetV0NJetClusterChi2	    = maxCluster->fitChi2;
  jetV0NJetClusterIntercept = maxCluster->interceptPv;
  jetV0NJetClusterAngle	    = maxCluster->cosAngleToFit;  
  jetV0NJetClusterAngleMom  = maxCluster->cosAngleToMomentum;  
  jetV0NJetClusterNTracks   = maxCluster->nTracks;  
  if (debug > 4) std::cout << "[DEBUG 4]  End NJet Cluster Calculation" << std::endl;  
}

void DisplacedJet::calcClusterSize(const DisplacedV0Collection& vertices, const float& errorWindow) {
  if (debug > 4) std::cout << "[DEBUG 4]  Finding Cluster Size" << std::endl;  
  
  const int nVtx = vertices.size();
  // nothing to do if there are no vertices
  if (nVtx == 0) return;

  // start with a NULL cluster
  DisplacedCluster * maxCluster = NULL;

  // check every vertex in a jet for the largest number of vertices within an error window
  // call the vertexing being checked the center
  if (debug > 4) std::cout << "[DEBUG 4]  Looping Center" << std::endl;  
  for(int center = 0; center < nVtx; ++center) {
    Displaced2TrackVertex centerVtx = vertices[center];
    if(!centerVtx.isValid) continue;

    // candidate cluster
    DisplacedCluster cluster_temp(centerVtx, selPV, iSetup, debug);

    if (debug > 4) std::cout << "[DEBUG 4]  Looping Neighbors" << std::endl;  
    // loop over all the possible neighbors
    for(int neighbor = 0; neighbor < nVtx; ++neighbor) {
      if (center == neighbor) continue; // dont count the center itself
      Displaced2TrackVertex neighVtx = vertices[neighbor];

      if(!neighVtx.isValid) continue;            
      // check the distance from the center
      float dx	     = centerVtx.x - neighVtx.x, dy = centerVtx.y - neighVtx.y;	//, dz = cz - nz;
      float distance = metric2D(dx, dy);
      float sig	     = distance / metric2D(neighVtx.tot_xyE, centerVtx.tot_xyE);

      if(sig < errorWindow) cluster_temp.addVertex(neighVtx);
      /* std::cout << " KEEP NEIGHBOR: " << std::endl; */
      /* std::cout << "neighb num" << neighbor << std::endl; */
      /* std::cout << "neighb position x: " << nx << " y: " << ny << " z: " << nz << std::endl; */
      /* std::cout << "neighb error xE: " << nxE << " yE: " << nyE << " zE: " << nzE << std::endl; */
      /* std::cout << "distance : " << distance << " significance: " << sig << " cluster temp size: " << tempClusterSize << std::endl; */      
    } // end neighbor loop

    // set the max cluster for the first iteration
    if(maxCluster == NULL) { 
      maxCluster = &cluster_temp;
      continue;
    }
    // otherwise, it needs to have more vertices
    if(cluster_temp.nV0 > maxCluster->nV0) maxCluster = &cluster_temp;

  } // end center loop
  
  /* std::cout << "@@@@@@@@@@@@@@@@@@@@" << std::endl; */
  /* std::cout << "max cluster size " << maxClusterSize << std::endl; */
  /* std::cout << "max cluster x " << maxX << " y: " << maxY << " z: " << maxZ << std::endl; */
  /* std::cout << "@@@@@@@@@@@@@@@@@@@@\n" << std::endl; */
  if (debug > 4) std::cout << "[DEBUG 4]  Assigning Max Vertex Information" << std::endl;  

  // build the related cluster quantities
  maxCluster->buildClusterQuantities();
  // set quantities related to the max cluster
  v0Cluster             = maxCluster;
  jetV0ClusterSize	= maxCluster->nV0;
  jetV0ClusterLxy	= maxCluster->center.lxy;
  jetV0ClusterLxyz	= maxCluster->center.lxyz;
  jetV0ClusterLxySig	= maxCluster->center.lxySig;
  jetV0ClusterLxyzSig	= maxCluster->center.lxyzSig;
  jetV0ClusterX		= maxCluster->center.x;
  jetV0ClusterY		= maxCluster->center.y;
  jetV0ClusterZ		= maxCluster->center.z;
  jetV0ClusterChi2	= maxCluster->fitChi2;
  jetV0ClusterIntercept = maxCluster->interceptPv;
  jetV0ClusterAngle	= maxCluster->cosAngleToFit;  
  jetV0ClusterAngleMom	= maxCluster->cosAngleToMomentum;  
  jetV0ClusterNTracks   = maxCluster->nTracks;

  if (debug > 4) std::cout << "[DEBUG 4]  End NJet Cluster Calculation" << std::endl;  
}

int DisplacedJet::calcHIndex(const DisplacedJet::vertexGraph & v0Graph) {
  if (debug > 4) std::cout << "[DEBUG 4]  calculating h index for graph" << std::endl;  

  int nTracks = v0Graph.size();
  std::vector<int> trackOccurances;

  // count the number of vertices per track
  for(int ii = 0; ii < nTracks; ++ii ) {
    trackOccurances.push_back(0);
    for(int jj = 0; jj < nTracks; ++jj ) {
      if(v0Graph[ii][jj] == 1) trackOccurances[ii] += 1;
    }
  }

  // check each possible hIndex up to the number of tracks
  // loop over the possible indices
  int graph_hindex = 0;
  for(int hindex = 0; hindex <= nTracks; ++hindex ) {
    // number of tracks that have at least hindex occurances
    int track_count_with_hindex = 0;
    for(int jj = 0; jj < nTracks; ++jj) {
      // track occurances satisified
      if (trackOccurances[jj] > hindex) track_count_with_hindex+=1; 
    }
    // hindex condition satisfied
    if (track_count_with_hindex >= hindex) graph_hindex = hindex;
  }

  return graph_hindex;
}

// check pairs of tracks for v0 candidates
void DisplacedJet::addV0Info(const reco::TrackRefVector trackRefs) {

  // get the transient track builder
  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  if (debug > 5) std::cout << "[DEBUG 5]  V0 candidate Tracks to Build: " << trackRefs.size() << std::endl;

  // loop over all pairs of tracks
  int nTracks = displacedTracks.size();
  for(int tt1 = 0; tt1 < nTracks; ++tt1) {
    // build the first track
    // loop over second track
    for(int tt2 = 0; tt2 < nTracks; ++tt2) {
      if(tt1 == tt2 || tt2 < tt1) continue;

      DisplacedTrack track1 = displacedTracks[tt1];
      DisplacedTrack track2 = displacedTracks[tt2];
      if (!track1.isValid || !track2.isValid) continue;

      // require 2d ip significance of 2 for each track
      if (fabs(track1.ip2dSig) < 2 || fabs(track2.ip2dSig) < 2) continue;
						      
      Displaced2TrackVertex vertex(track1, track2, selPV, iSetup, debug);   

      if(vertex.chi2 > 20 || !vertex.isValid) continue;                           
      
      displacedV0Vector.push_back(vertex);
    } // end loop over second track
  } // end loop over first track


  int v0VectorSize = displacedV0Vector.size();
  if (debug > 5) std::cout << "[DEBUG 5]  looping over reconstructed vertices" << std::endl;
  // loop over the reconstructed vertices
  for(int vv = 0; vv < v0VectorSize; ++vv){

    Displaced2TrackVertex vertex = displacedV0Vector[vv];

    // checks that the inner hits are not consistent with nuclear interactions
    // the vertex should not be close to the inner hit in the layer
    static const float MIN_DISTANCE = 0.05;    
    bool    oneTrackNuclear	 = (vertex.dInnerHitTrack1 < MIN_DISTANCE) || (vertex.dInnerHitTrack2 < MIN_DISTANCE);
    bool    twoTrackNuclear	 = (vertex.dInnerHitTrack1 < MIN_DISTANCE) && (vertex.dInnerHitTrack2 < MIN_DISTANCE);
    //bool    twoTrackInnerHitFake = vertex.dInnerHitTrack1Track2 < MIN_DISTANCE;

    // make sure the track inner hits are in the pixel barrel and the the vertex is 
    // transversely close to a given layer 
    bool    vertexNearBeamPipe = vertex.lxy < 2.3 && vertex.lxy > 2.1 && vertex.tot_xyE < 0.03;
    bool    vertexNearBPIX1    = vertex.lxy < 5 && vertex.lxy > 4 && vertex.tot_xyE < 0.03;
    bool    vertexNearBPIX2    = vertex.lxy > 6.8 && vertex.lxy < 7.8 && vertex.tot_xyE < 0.03;
    bool    vertexNearBPIX3    = false;
    bool    vertexNearBPIX     = vertexNearBeamPipe || vertexNearBPIX1 || vertexNearBPIX2 || vertexNearBPIX3;

    if(oneTrackNuclear) jetOneTrackNuclearCount++;
    if(twoTrackNuclear) jetTwoTrackNuclearCount++;
    if(vertexNearBPIX1) jetVertexNearBPIX1++;
    if(vertexNearBPIX2) jetVertexNearBPIX2++;
    if(vertexNearBPIX3) jetVertexNearBPIX3++;
    if(vertexNearBPIX)  jetVertexNearBPIX++;
    if(vertexNearBPIX && twoTrackNuclear) jetTightNuclear++;
    if(vertexNearBPIX && oneTrackNuclear) jetLooseNuclear++;
    // restrict the inner hit behind calls to tracks that can enter the medianip significance calcualtion
    if((vertex.innerHitBehindVertex1)  && (vertex.innerHitBehindVertex2)) jetNV0HitBehindVertex++;
    //if(twoTrackInnerHitFake) 
    // count all vertices (including the ones from tracks which do not enter medianip significance caluclation)
    // require the vertices arent fake with hits behind the position
    // must have 2D significance of at least 1 as to not be consistent with the beamspot
    // must be outside .5 mm 
    if(!vertex.innerHitBehindVertex1 && !vertex.innerHitBehindVertex2 && vertex.chi2 < 20 && vertex.lxy > 0.05 && fabs(vertex.lxySig) > 3) {
      jetNV0NoHitBehindVertex++;      
      displacedV0VectorCleaned.push_back(vertex);
    }

    // particle comparisons
    jetNV0KShort += vertex.isKShort ? 1 : 0;
    jetNV0Lambda += vertex.isLambda ? 1 : 0;
  } // end loop on v0 candiate vertices  

  // create the graph of vertex association and find the h index
  //  findVertexCliques();
  // fill information for the clsutering of the vertices
  calcClusterSize(displacedV0VectorCleaned, 2);  
}



// compute variables related to the track angles and the calo jet 
void DisplacedJet::addHitInfo(const reco::TrackCollection tracks) {
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
      //      const DetLayer &  detLayerInner	 = *(trajInfo[0].detLayer); 
      // const DetLayer &  detLayerOuter	 = *(trajInfo.back().detLayer); 
      
      //const bool	isBarInner	 = detLayerInner.isBarrel();
      //const bool	isBarOuter	 = detLayerOuter.isBarrel();
      /* GeomDetEnumerators::SubDetector subDetLayerInner = detLayerInner.subDetector(); */
      /* GeomDetEnumerators::SubDetector subDetLayerOuter = detLayerOuter.subDetector(); */
      
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

      bool hasPixelHits = tIter->hitPattern().numberOfValidPixelHits() > 0;
      // if the inner hit is in the pixel layers
      if (hasPixelHits) {
	innerHitPosInPixel.push_back(ri);
	jetNTracksPixel	    += 1;	
	jetPtSumTracksPixel += pt;
      }
      else {
	innerHitPosOutPixel.push_back(ri);
	jetNTracksNoPixel     += 1;
	jetPtSumTracksNoPixel += pt;	  
      }

      if (hasPixelHits) outerHitPosInPixel.push_back(ro);
      else outerHitPosOutPixel.push_back(ro);     
    } // end hits valid
  } //end loop filling vectors of hits

  bool check_medians = false;
  if(check_medians) {
    // test that the median calculation is working
    float even[]		= {1, 3.2 , 4, 5};
    float odd[]		= {1, 3.2 , 4, 5, 9};
    float inverted_odd[]	= {3.2, 5 , 4, 5, 9};
    float inverted_even[] = {3.2 ,5 ,4,1};
    float large_nums[]	= {3,4,7,6,5,4,7,5,2,3,4,5,6,7,8,9,10,10,10,10,7,6,5,4,7,5,7,6,5,4,7,5};
    float one[]		= {3};
    float pos_neg[]	= {-1.3, 3.2, 5, -4, 1};

    // vectors
    std::vector<float> even_vec (even, even + sizeof(even) / sizeof(int));
    std::vector<float> odd_vec (odd, odd + sizeof(odd) / sizeof(int));
    std::vector<float> inverted_odd_vec (inverted_odd, inverted_odd + sizeof(inverted_odd) / sizeof(int));
    std::vector<float> inverted_even_vec (inverted_even, inverted_even + sizeof(inverted_even) / sizeof(int));
    std::vector<float> large_nums_vec (large_nums, large_nums + sizeof(large_nums) / sizeof(int));
    std::vector<float> one_vec (one, one + sizeof(one) / sizeof(int));
    std::vector<float> pos_neg_vec(pos_neg, pos_neg + sizeof(pos_neg) / sizeof(int));  

    // check the jet medians
    float even_diff	   = fabs(getJetMedian(even_vec, false) - 3.6);
    float odd_diff	   = fabs(getJetMedian(odd_vec, false) - 4);
    float inverted_odd_diff  = fabs(getJetMedian(inverted_odd_vec, false) - 5);
    float inverted_even_diff = fabs(getJetMedian(inverted_even_vec, false) - 3.6);
    float large_nums_diff    = fabs(getJetMedian(large_nums_vec, false) - 6);
    float one_diff	   = fabs(getJetMedian(one_vec, false) - 3);
    float pos_neg_diff	   = fabs(getJetMedian(pos_neg_vec, true) - 1);

    // check the jet medians
    if(even_diff > 0.001) std::cout << "[jetMedian TEST] even array diff val: " << even_diff << " median= " << getJetMedian(even_vec, false)<< std::endl;
    if(odd_diff > 0.001) std::cout << "[jetMedian TEST]  odd  array diff val: " << odd_diff <<  " median= " << getJetMedian(odd_vec, false) << std::endl;
    if(inverted_odd_diff > 0.001) std::cout << "[jetMedian TEST] inverted odd array diff val: " << inverted_odd_diff <<  " median= " << getJetMedian(inverted_odd_vec, false) << std::endl;
    if(inverted_even_diff > 0.001) std::cout << "[jetMedian TEST] inverted even array diff val: " <<  inverted_even_diff <<  " median= " << getJetMedian(inverted_even_vec, false) << std::endl;
    if(large_nums_diff > 0.001) std::cout << "[jetMedian TEST] large nums array diff val: " << large_nums_diff <<  " median= " << getJetMedian(large_nums_vec, false) << std::endl;
    if(one_diff > 0.001) std::cout << "[jetMedian TEST] one array diff val: " << one_diff << " median= " << getJetMedian(one_vec, false) << std::endl;
    if(pos_neg_diff > 0.001) std::cout << "[jetMedian TEST] pos neg diff diff val: " << pos_neg_diff  << " median= " << getJetMedian(pos_neg_vec, false) << std::endl;
  }
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
void DisplacedJet::addTrackAngles(const DisplacedTrackCollection & tracks) {
  DisplacedTrackCollection::const_iterator tIter = tracks.begin();
  
  // pt weighted
  ptSumCosTheta2D	  = 0, ptSumCosTheta3D = 0;
  ptSumCosThetaDet2D	  = 0, ptSumCosThetaDet3D = 0;
  // aboslute sum
  sumCosTheta2D	          = 0, sumCosTheta3D = 0;
  sumCosThetaDet2D	  = 0, sumCosThetaDet3D = 0;   

  // find the angle between the jet momentum and track momentum
  float sumTrackPtValid = 0;  
  for(; tIter != tracks.end(); ++tIter) {    
    if (tIter->isValid) {
      float pt = tIter->pt; 

      if(pt < 1) continue;

      // track cosine with respect to the momentum at the rference point
      // momentum is already defined at the reference point
      float   cosTheta3D    = tIter->angleMomentumAndPVAtInnerHit3D;
      float   cosTheta2D    = tIter->angleMomentumAndPVAtInnerHit2D;
      // track cosine with respect to the inner hit 
      float   cosThetaDet3D = tIter->angleMomentumAndPVAtOuterHit3D;
      float   cosThetaDet2D = tIter->angleMomentumAndPVAtOuterHit2D;    //trkv2_unit * jetdv2_unit;
      /* std::cout << "COSTHETA3D " << cosTheta3D << std::endl;  */
      /* std::cout << "COSTHETA2D " << cosTheta2D << std::endl;  */
      /* std::cout << " -----------------" << std::endl;  */

      // store the track kinematics in the jet
      trEtaVector.push_back(tIter->eta);      
      trPhiVector.push_back(tIter->phi);      
      trPtVector.push_back(pt);      

      cosTheta2DVector.push_back(cosTheta2D);
      cosTheta3DVector.push_back(cosTheta3D);
      cosThetaDet2DVector.push_back(cosThetaDet2D);
      cosThetaDet3DVector.push_back(cosThetaDet3D);
    
      // increment sumpt
      sumTrackPtValid	 += pt;

      // pt weighted
      ptSumCosTheta2D    += sin(cosTheta2D) *    pt;
      ptSumCosTheta3D    += sin(cosTheta3D) *    pt;
      ptSumCosThetaDet2D += sin(cosThetaDet2D) * pt;
      ptSumCosThetaDet3D += sin(cosThetaDet3D) * pt;
      // not pt weighted
      sumCosTheta2D      += sin(cosTheta2D);
      sumCosTheta3D      += sin(cosTheta3D);
      sumCosThetaDet2D   += sin(cosThetaDet2D);
      sumCosThetaDet3D   += sin(cosThetaDet3D);
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
  trackSumMomCosTheta2D	 = 0;
  trackSumMomCosTheta3D	 = 0;
  trackSumMomMag2D	 = 0;
  trackSumMomMag3D	 = 0;
  // vector sum of IP
  ipPosSumMag3D		 = 0;
  ipPosSumMag2D		 = 0;
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
  float dx = selPV.x() - svX , dy = selPV.y() - svY, dz = selPV.z() - svZ;

  TVector3 pvVector3D(selPV.x(), selPV.y(), selPV.z());
  TVector3 pvVector2D(selPV.x(), selPV.y(), 0);
  TVector3 svVector3D(svX, svY, svZ);
  TVector3 svVector2D(svX, svY, 0);

  // line pointing form the primary vertex through the sceondary vertex
  TVector3 svMom3D( selVertex.p4().x(), selVertex.p4().y(), selVertex.p4().z());
  TVector3 svMom2D( selVertex.p4().x(), selVertex.p4().y(), 0);

  // you want the negative when the momentum and sv are in the same 
  // direction relative to the PV
  // this makes sure the angle is not pi when the vertex is fit behind 
  // the primary vertex 
  float sign2D =  (svMom2D * (svVector2D - pvVector2D)) > 0 ? -1: 1;
  float sign3D =  (svMom3D * (svVector3D - pvVector3D)) > 0 ? -1: 1;

  TVector3 pvToVertex3D( sign3D * dx, sign3D * dy, sign3D * dz);
  TVector3 pvToVertex2D( sign2D * dx, sign2D * dy, 0);

  svAngle3D = pvToVertex3D.Angle(svMom3D);
  svAngle2D = pvToVertex2D.Angle(svMom2D);
  
  svMass    = selVertex.p4().mass();    
  svPt      = selVertex.p4().pt();
  svEta     = selVertex.p4().eta();
  svPhi     = selVertex.p4().phi();
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
    if(jetV0ClusterSize > (3 + 3*(*thresIter)) && jetV0ClusterAngleMom > 0.2 && jetV0ClusterAngleMom < 2.8) didPass = true;
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
    if(nTracks > 0 && medianIPLogSig2D > 0  && (alphaMax / sumTrackPt) < (.8 - (*thresIter * .14)) ) didPass = true;
    //    if( medianCosThetaDet2D > 0.05 && nTracks > 2 && medianIPLogSig2D > (*thresIter)
    //	&& (alphaMax / sumTrackPt) < 0.1) didPass = true;
    passResults.push_back(didPass);
  }
  shortTagsVector = passResults;
  return passResults;  
}

std::vector<bool> DisplacedJet::passMediumTag(const std::vector<float> thres, float min, float max){
  std::vector<bool> passResults;
  std::vector<float>::const_iterator thresIter = thres.begin();
  for(; thresIter != thres.end(); ++thresIter) {
    bool result = nTracks > 0 && medianIPLogSig2D > 2 && (alphaMax / sumTrackPt) < 0.1 && medianCosThetaDet2D > 0.05;
    passResults.push_back(result);
    //if (nTracks > 0 && medianIPLogSig2D > (.5 * (*thresIter)) && (alphaMax / sumTrackPt) < .8) didPass = true;
    passResults.push_back(result);
  }
  mediumTagsVector = passResults;
  return passResults;  
}

std::vector<bool> DisplacedJet::passLongTag(const std::vector<float> thres, float min, float max){
  std::vector<bool> passResults;
  std::vector<float>::const_iterator thresIter = thres.begin();
  for(; thresIter != thres.end(); ++thresIter) {
    bool result = nTracks > 0 && medianIPLogSig2D > 1.0 && (alphaMax / sumTrackPt) < 0.5 && medianCosThetaDet2D > 0.02;
    //    if( medianCosThetaDet2D > 0.05 && nTracks > 0 && medianIPLogSig2D > (*thresIter)
    //	&& (alphaMax / sumTrackPt) < 0.05) didPass = true;
    passResults.push_back(result);
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

  //  if (debug > 4) std::cout << "[DEBUG] jet median:  " << dMedian << std::endl;
  return dMedian;
}

float DisplacedJet::getJetMean(const std::vector<float> & values, bool is_signed) {
  if (values.size() == 0) return 0.0;

  float sum = 0;
  std::vector<float>::const_iterator val = values.begin();
  for (; val != values.end(); ++val) sum += is_signed ? *val : fabs(*val);  

  float mean = sum / float(values.size());
  //  if (debug > 4) std::cout << "[DEBUG] jet mean:  " << mean << std::endl;
  return mean;
}

float DisplacedJet::getJetVariance(const std::vector<float>& values, bool is_signed) {
  return 0;
  /* if (values.size() == 0 || values.size() == 1) return 0.0; */

  /* float sum = 0; */
  /* float mean = getJetMean(values, is_signed); */

  /* std::vector<float>::const_iterator val = values.begin(); */
  /* for (; val != values.end(); ++val) { */
  /*   float temp =  (*val - mean) * (*val - mean); */
  /*   sum += temp; */
  /* } */

  /* //safety check */
  /* if(!is_signed ) { assert( sum >= 0 ); } */
  /* float variance = std::sqrt(sum / float((values.size() - 1.0))); */

  /* //if (debug > 4) std::cout << "[DEBUG] jet variance:  " << variance << std::endl; */
  /* return variance; */
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
      if (dr < 0.4 ) sumJetPt_temp += (*tIter)->pt();
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


