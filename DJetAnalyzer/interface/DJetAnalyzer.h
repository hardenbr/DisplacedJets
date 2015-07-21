class DJetAnalyzer : public edm::EDAnalyzer {

 public:
  explicit  DJetAnalyzer(const edm::ParameterSet&);
  ~DJetAnalyzer();
  
  static void	fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
 private:
  void fillHandles(const edm::Event &);

  // gnerator
  void dumpGenInfo(const reco::GenParticleCollection &); 
  void dumpSimInfo(const edm::SimVertexContainer &);
  void dumpPreSelection(DisplacedJetEvent&);

  // trigger information
  void fillTriggerInfo(const edm::Event &  iEvent, const edm::TriggerResults & trigResults);

  // tree dumping displaced jet quantities
  void dumpCaloInfo(DisplacedJetEvent&);
  void dumpSVTagInfo(DisplacedJetEvent&);
  void dumpIPInfo(DisplacedJetEvent&);
  void dumpIVFInfo(DisplacedJetEvent&);
  void dumpDJTags(DisplacedJetEvent&);
  void dumpPVInfo(DisplacedJetEvent &, const reco::VertexCollection &);
  //tree dumping track quantities
  void dumpTrackInfo(DisplacedJetEvent&, const reco::TrackCollection &, const int & collectionID, const edm::EventSetup & iSetup);
  
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //file configuration tags
  std::string	outputFileName_;
  std::string	jetTreeName_;
  std::string	trackTreeName_;
  std::string	vertexTreeName_;
  std::string	genTreeName_;

  TFile *outputFile_;
  // analysis to dos
  bool	 doGenMatch_;
  bool	 doSimMatch_;
  bool	 applyEventPreSelection_;
  bool	 applyJetPreSelection_;
  bool   dumpGeneralTracks_;
  // keep trees
  bool   writeJetTree_;
  bool   writeTrackTree_;
  bool   writeEventTree_;
  bool   writeGenTree_;
  bool   writeVertexTree_;
  // event flags
  bool	 isMC_;
  bool	 isSignalMC_;
  
  //trigger tags
  std::string triggerResultPath_;
  edm::InputTag tag_triggerResults_;

  //tracking tags
  edm::InputTag tag_generalTracks_;
  edm::InputTag tag_trackIPTagInfoCollection_;
  edm::InputTag tag_lifetimeIPTagInfo_;
  edm::InputTag tag_secondaryVertexTagInfo_; 

  //vertex tags
  edm::InputTag tag_secondaryVertices_;
  edm::InputTag tag_inclusiveVertexCandidates_;
  edm::InputTag tag_inclusiveSecondaryVertices_;
  edm::InputTag tag_offlinePrimaryVertices_;

  //jet tags
  edm::InputTag tag_ak4CaloJets_;
  
  //gen info
  edm::InputTag tag_genParticles_;
  edm::InputTag tag_ak4GenJets_;
  edm::InputTag tag_genMetCalo_;
  
  //sim Tags
  edm::InputTag tag_simVertex_;    

  // tag threshold classifications
  float shortTagThresDist, mediumTagThresDist, longTagThresDist;
  int dHTWorkingPoint;
  //cuts
  float cut_jetPt, cut_jetEta;
  
  //output related
  TTree*    trackTree_;   
  TTree*    jetTree_;   
  TTree*    vertexTree_;
  TTree*    genTree_;
  TTree*    eventTree_; 
  TTree*    runStatTree_;

  static const Int_t    SIM_STATUS_CODE_MATCH = 0; 
  static const Int_t    GEN_STATUS_CODE_MATCH = 23; 
  const float		VERTEX_MATCH_METRIC   = 0.05;
  static const Int_t	MAX_TRIGGERS	      = 1000;
  static const Int_t	MAX_TRACKS	      = 5000;
  static const Int_t	MAX_JETS	      = 40;
  static const Int_t	MAX_VTX		      = 100;
  static const Int_t	MAX_CAT		      = 100; // max number of tagging categories
  static const Int_t	MAX_GEN		      = 500;
  static const Int_t	FAKE_HIGH_VAL	      = 9999;

  Int_t	debug = 0; 
 
  //bookkeeping
  Int_t run   = -1;
  Int_t lumi  = -1;
  Int_t event = -1;

  int	evNum = 0;
  Int_t jetid = 0;

  //tree variables
  Int_t nLiTracks = 0;
  Int_t nCaloJets = 0;
  Int_t nTracks	  = 0; // general tracks

  ////////////////// EVENT TREE SPECIFIC MEMBERS//////////

  Float_t eventCaloHT;
  Float_t eventCaloDHT[MAX_CAT];
  Float_t eventCaloMET;

  Int_t eventPassEventPreSelection;
  Int_t eventNJetsPassPreSelection;

  Int_t eventNIVFReco;
  Int_t eventNIVFRecoGenMatch;

  // number of working points for each category
  // loose, medium, tight --> nWP = 3
  Int_t nWP;

  // the decay regimes (split by lxyz of vtx)
  Int_t eventNNoVertexTags[MAX_CAT];
  Int_t eventNShortTags[MAX_CAT];
  Int_t eventNMediumTags[MAX_CAT];
  Int_t eventNLongTags[MAX_CAT];
  Int_t eventNTotalTags[MAX_CAT];

  ////////////////// CALO JETS ////////////

  // jet kinematics
  Int_t	    caloJetID[MAX_JETS];
  // ordered flat numbers
  // by pt
  Float_t   caloLeadingJetPT;
  Float_t   caloSubLeadingJetPT;
  // All Iterations
  Float_t   caloFewestPromptTracks;
  Float_t   caloSubFewestPromptTracks;
  Float_t   caloMostDispTracks;
  Float_t   caloSubMostDispTracks;
  // HLT Iterations
  Float_t   caloFewestPromptTracksHLT;
  Float_t   caloSubFewestPromptTracksHLT;
  Float_t   caloMostDispTracksHLT;
  Float_t   caloSubMostDispTracksHLT;

  // by hadronic fraction for pt > 40 GeV
  Float_t   caloLeadingHadronicFraction;
  // vbf numbers
  // Mqq for minimum dEta 3.0 max dEta 5 and min pt 20
  Float_t   caloLeadingMqq;

  Float_t   caloJetPt[MAX_JETS];
  Float_t   caloJetEta[MAX_JETS];
  Float_t   caloJetPhi[MAX_JETS];
  // size
  Float_t   caloJetN90[MAX_JETS];
  Float_t   caloJetN60[MAX_JETS];
  Float_t   caloJetTowerArea[MAX_JETS];
  // energy contribution
  Float_t   caloJetHfrac[MAX_JETS];
  Float_t   caloJetEfrac[MAX_JETS];

  ///////////////////// GEN MATCHED ////////////////////

  Int_t genMatch[MAX_JETS];
  Int_t genMatchTrack[MAX_TRACKS];

  Int_t     caloGenMatch[MAX_JETS];
  Float_t   caloGenPt[MAX_JETS];
  Float_t   caloGenEta[MAX_JETS];
  Float_t   caloGenPhi[MAX_JETS];
  Float_t   caloGenM[MAX_JETS];


  ///////////////////// TRACK MATCHING ////////////////////

  Int_t   jetNTracks[MAX_JETS];  

  //////////////////// LIFETIME TAG /////////////////  

  // track kinematics
  Float_t   liTrackPt[MAX_TRACKS];
  Float_t   liTrackEta[MAX_TRACKS];
  Float_t   liTrackPhi[MAX_TRACKS];
  
  // lifetime jet kinematics
  Int_t	    nLiJets = 0;
  Int_t	    liJetID[MAX_JETS];
  Float_t   liJetPt[MAX_JETS];
  Float_t   liJetEta[MAX_JETS];
  Float_t   liJetPhi[MAX_JETS];
  Int_t	    liJetNSelTracks[MAX_JETS];

  // lifetime track infoTags
  Float_t   liTrackIP3D[MAX_TRACKS];
  Float_t   liTrackIP2D[MAX_TRACKS];
  Float_t   liTrackIPSig3D[MAX_TRACKS];
  Float_t   liTrackIPSig2D[MAX_TRACKS];
  Float_t   liTrackDistanceJetAxis[MAX_TRACKS];
  Float_t   liTrackDistanceJetAxisSig[MAX_TRACKS];
  Float_t   liTrackChi2[MAX_TRACKS];
  
  //tracks to jet comparisons
  Int_t	    liTrackJetID[MAX_TRACKS];
  Float_t   liJetTrackDR[MAX_TRACKS];

  //////////////////// SV TAG ////////////////

  Int_t nSV = 0;

  Int_t	    svVertexJetID[MAX_VTX]; 
  Int_t	    nSvJets = 0;
  Float_t   svJetPt[MAX_JETS];
  Float_t   svJetEta[MAX_JETS];
  Float_t   svJetPhi[MAX_JETS];
  Float_t   svJetID[MAX_JETS];

  //kinematics
  Float_t   svMass[MAX_VTX];
  Float_t   svPx[MAX_VTX];
  Float_t   svPy[MAX_VTX];
  Float_t   svPt[MAX_VTX];
  Float_t   svEta[MAX_VTX];
  Float_t   svPhi[MAX_VTX];

  //position
  Float_t   svX[MAX_VTX];
  Float_t   svY[MAX_VTX];
  Float_t   svZ[MAX_VTX];

  //position error
  Float_t   svXErr[MAX_VTX];
  Float_t   svYErr[MAX_VTX];
  Float_t   svZErr[MAX_VTX];

  //quality
  Float_t   svChi2[MAX_VTX];
  Float_t   svNChi2[MAX_VTX];
  Float_t   svNDof[MAX_VTX];
  Int_t	    svIsValid[MAX_VTX];
  
  //tracking
  Int_t	    svNTracks[MAX_VTX];
  //  Int_t svTrackVertexID[MAX_VTX];   
  Float_t   svTotalCharge[MAX_VTX];

  // Flight Information
  Float_t   svFlight[MAX_VTX];
  Float_t   svFlightErr[MAX_VTX];
  Float_t   svFlight2D[MAX_VTX];
  Float_t   svFlight2DErr[MAX_VTX];

  // various DR 
  Float_t   svDRFlightJet[MAX_VTX];
  Float_t   svDRTrackJet[MAX_VTX];
  Float_t   svDRTrackFlight[MAX_VTX];

  ///////////////// JET TREE SPECIFIC MEMBERS /////////////////

  // analysis tracking
  Int_t jetPassPreSelection[MAX_JETS];

  //book keeping
  Int_t nJetWithSv; 

  // significance and aboslute IP weighted track energy
  Float_t   jetEIPSig2D[MAX_JETS];
  Float_t   jetEIPSigLog2D[MAX_JETS];
  Float_t   jetEIPSig3D[MAX_JETS];
  Float_t   jetEIPSigLog3D[MAX_JETS];

  // significance log IP weighted track energy
  Float_t   jetELogIPSig2D[MAX_JETS];
  Float_t   jetELogIPSig3D[MAX_JETS];
  
  // absolute IP sums
  Float_t   jetIPSum2D[MAX_JETS];
  Float_t   jetIPSum3D[MAX_JETS];
  // IP significance sums
  Float_t   jetIPSigSum2D[MAX_JETS];
  Float_t   jetIPSigSum3D[MAX_JETS];
  // IP significance log sums
  Float_t   jetIPSigLogSum2D[MAX_JETS];
  Float_t   jetIPSigLogSum3D[MAX_JETS];
  Float_t   jetIPLogSum2D[MAX_JETS];
  Float_t   jetIPLogSum3D[MAX_JETS];

  Float_t   jetDistLogSum[MAX_JETS];
  Float_t   jetDistSigLogSum[MAX_JETS];

  // IP significance averages
  // means
  Float_t   jetMeanIPSig2D[MAX_JETS];
  Float_t   jetMeanIPSig3D[MAX_JETS];
  Float_t   jetMeanIPLogSig2D[MAX_JETS];
  Float_t   jetMeanIPLogSig3D[MAX_JETS];
  Float_t   jetMeanJetDistSig[MAX_JETS];
  // median
  Float_t   jetMedianIPSig2D[MAX_JETS];
  Float_t   jetMedianIPSig3D[MAX_JETS];
  Float_t   jetMedianIPLogSig2D[MAX_JETS];
  Float_t   jetMedianIPLogSig3D[MAX_JETS];
  Float_t   jetMedianJetDistSig[MAX_JETS];
  // variance
  Float_t   jetVarianceIPSig2D[MAX_JETS];
  Float_t   jetVarianceIPSig3D[MAX_JETS];
  Float_t   jetVarianceIPLogSig2D[MAX_JETS];
  Float_t   jetVarianceIPLogSig3D[MAX_JETS];
  Float_t   jetVarianceJetDistSig[MAX_JETS];

  // Absolute IP averages
  // mean
  Float_t   jetMeanIP2D[MAX_JETS];
  Float_t   jetMeanIP3D[MAX_JETS];
  Float_t   jetMeanIPLog2D[MAX_JETS];
  Float_t   jetMeanIPLog3D[MAX_JETS];
  Float_t   jetMeanJetDist[MAX_JETS];
  // median
  Float_t   jetMedianIP2D[MAX_JETS];
  Float_t   jetMedianIP3D[MAX_JETS];
  Float_t   jetMedianIPLog2D[MAX_JETS];
  Float_t   jetMedianIPLog3D[MAX_JETS];
  Float_t   jetMedianJetDist[MAX_JETS];
  // variance
  Float_t   jetVarianceIP2D[MAX_JETS];
  Float_t   jetVarianceIP3D[MAX_JETS];
  Float_t   jetVarianceIPLog2D[MAX_JETS];
  Float_t   jetVarianceIPLog3D[MAX_JETS];
  Float_t   jetVarianceJetDist[MAX_JETS];

  // SV information
  Int_t	    jetNSv[MAX_JETS];
  Float_t   jetSvMass[MAX_JETS];
  Float_t   jetSvLxy[MAX_JETS];
  Float_t   jetSvLxySig[MAX_JETS];
  Float_t   jetSvLxyz[MAX_JETS];
  Float_t   jetSvLxyzSig[MAX_JETS];
  Int_t	    jetSvNTrack[MAX_JETS];  //vertex track multiplicty

  // SV position
  Float_t   jetSvX[MAX_JETS];
  Float_t   jetSvY[MAX_JETS];
  Float_t   jetSvZ[MAX_JETS];

  Float_t   jetSvZErr[MAX_JETS];
  Float_t   jetSvYErr[MAX_JETS];
  Float_t   jetSvXErr[MAX_JETS];

  // SV quality
  Float_t   jetSvChi2[MAX_JETS];
  Float_t   jetSvNChi2[MAX_JETS];
  Float_t   jetSvNDof[MAX_JETS];
  Int_t	    jetSvIsValid[MAX_JETS];

  // gen vertex matching and score
  Int_t	    jetSvGenVertexMatched[MAX_JETS];
  Float_t   jetSvGenVertexMatchMetric[MAX_JETS];
  // sim vtx matching
  Int_t	    jetSvSimVertexMatched[MAX_JETS];
  Float_t   jetSvSimVertexMatchMetric[MAX_JETS];

  // total number matched
  //n gen matched
  Int_t	    jetSvNGenMatched;
  Int_t	    jetSvNGenFake;
  // n sim matched
  Int_t	    jetSvNSimMatched;
  Int_t	    jetSvNSimFake;

  ///////////////////////////// IVF ///////////////////////

  // IVF Information
  Float_t   jetIVFMass[MAX_JETS];
  Float_t   jetIVFLxy[MAX_JETS];
  Float_t   jetIVFLxySig[MAX_JETS];
  Float_t   jetIVFLxyz[MAX_JETS];
  Float_t   jetIVFLxyzSig[MAX_JETS];
  Int_t	    jetIVFNTrack[MAX_JETS];  

  // IVF position
  Float_t   jetIVFX[MAX_JETS];
  Float_t   jetIVFY[MAX_JETS];
  Float_t   jetIVFZ[MAX_JETS];

  Float_t   jetIVFZErr[MAX_JETS];
  Float_t   jetIVFYErr[MAX_JETS];
  Float_t   jetIVFXErr[MAX_JETS];

  // IVF matching score
  Float_t   jetIVFMatchingScore[MAX_JETS];
  // mother id
  Int_t	    jetIVFVertexIDNMom;
  Float_t   jetIVFVertexIDMomPt[MAX_VTX];
  Int_t	    jetIVFVertexIDMomJetID[MAX_VTX];
  Int_t	    jetIVFVertexIDMomHighestPtID[MAX_JETS]; // jet index 
  Float_t   jetIVFVertexIDMomHighestPt[MAX_JETS]; // jet indexed
  Int_t	    jetIVFVertexIDMom[MAX_VTX];
  // son id
  Int_t	    jetIVFVertexIDNSon;
  Int_t	    jetIVFVertexIDSon[MAX_VTX];
  Float_t   jetIVFVertexIDSonPt[MAX_VTX];
  Int_t	    jetIVFVertexIDSonJetID[MAX_VTX];
  // IVF gen matching
  Int_t	    jetIVFGenVertexMatched[MAX_JETS];
  Float_t   jetIVFGenVertexMatchMetric[MAX_JETS];
  Int_t	    jetIVFSimVertexMatched[MAX_JETS];
  Float_t   jetIVFSimVertexMatchMetric[MAX_JETS];

  ///////////////// VERTEX TREE SPECIFIC MEMBERS /////////////////
  

  // Secondary Vertex Candidates from standalone
  Int_t	    vtxN;
  Float_t   vtxIsFake[MAX_VTX];
  Int_t	    vtxNTracks[MAX_VTX];
  Float_t   vtxChi2[MAX_VTX];
  Int_t	    vtxNDof[MAX_VTX];

  Float_t   vtxX[MAX_VTX];
  Float_t   vtxY[MAX_VTX];
  Float_t   vtxZ[MAX_VTX];
  Float_t   vtxLxy[MAX_VTX];
  Float_t   vtxLxyz[MAX_VTX];

  Float_t   vtxXSig[MAX_VTX];
  Float_t   vtxYSig[MAX_VTX];
  Float_t   vtxZSig[MAX_VTX];
  Float_t   vtxLxySig[MAX_VTX];

  // Inclusive Vertex Candidates 
  Int_t	    vtxIncCandN;
  Float_t   vtxIncCandIsFake[MAX_VTX];
  Int_t	    vtxIncCandNTracks[MAX_VTX];
  Float_t   vtxIncCandChi2[MAX_VTX];
  Int_t	    vtxIncCandNDof[MAX_VTX];

  Float_t   vtxIncCandX[MAX_VTX];
  Float_t   vtxIncCandY[MAX_VTX];
  Float_t   vtxIncCandZ[MAX_VTX];
  Float_t   vtxIncCandLxy[MAX_VTX];
  Float_t   vtxIncCandLxyz[MAX_VTX];

  Float_t   vtxIncCandXSig[MAX_VTX];
  Float_t   vtxIncCandYSig[MAX_VTX];
  Float_t   vtxIncCandZSig[MAX_VTX];
  Float_t   vtxIncCandLxySig[MAX_VTX];
  Float_t   vtxIncCandLxyzSig[MAX_VTX];

  Int_t	    vtxIncCandNGenMatched;
  Int_t	    vtxIncCandNGenFake;
  Int_t	    vtxIncCandGenMatched[MAX_VTX];
  Float_t   vtxIncCandGenMatchMetric[MAX_VTX];
  Int_t	    vtxIncCandNSimMatched;
  Int_t	    vtxIncCandNSimFake;
  Int_t	    vtxIncCandSimMatched[MAX_VTX];
  Float_t   vtxIncCandSimMatchMetric[MAX_VTX];

  // Inclusive Secondary  (post merge and arbitration of IVF Candidates)
  Int_t     vtxIncSecN;
  Float_t   vtxIncSecIsFake[MAX_VTX];
  Int_t	    vtxIncSecNTracks[MAX_VTX];
  Float_t   vtxIncSecChi2[MAX_VTX];
  Int_t	    vtxIncSecNDof[MAX_VTX];

  Float_t   vtxIncSecX[MAX_VTX];
  Float_t   vtxIncSecY[MAX_VTX];
  Float_t   vtxIncSecZ[MAX_VTX];
  Float_t   vtxIncSecLxy[MAX_VTX];
  Float_t   vtxIncSecLxyz[MAX_VTX];

  Float_t   vtxIncSecXSig[MAX_VTX];
  Float_t   vtxIncSecYSig[MAX_VTX];
  Float_t   vtxIncSecZSig[MAX_VTX];
  Float_t   vtxIncSecLxySig[MAX_VTX];  
  Float_t   vtxIncSecLxyzSig[MAX_VTX];  

  Int_t	    vtxIncSecNGenMatched;
  Int_t	    vtxIncSecNGenFake;
  Int_t	    vtxIncSecGenMatched[MAX_VTX];
  Float_t   vtxIncSecGenMatchMetric[MAX_VTX];
  Int_t	    vtxIncSecNSimMatched;
  Int_t	    vtxIncSecNSimFake;
  Int_t	    vtxIncSecSimMatched[MAX_VTX];
  Float_t   vtxIncSecSimMatchMetric[MAX_VTX];

  ///////////////// GENERATOR TREE SPECIFIC MEMBERS /////////////////

  // Qualities
  Int_t genPartN;
  Int_t genPartPID[MAX_GEN];
  Int_t genPartStatus[MAX_GEN];

  // Position
  Float_t   genPartPt[MAX_GEN];
  Float_t   genPartEta[MAX_GEN];
  Float_t   genPartPhi[MAX_GEN];

  // Vertex
  Float_t   genPartVX[MAX_GEN];
  Float_t   genPartVY[MAX_GEN];
  Float_t   genPartVZ[MAX_GEN];
  Float_t   genPartVLxy[MAX_GEN];
  Float_t   genPartVLxyz[MAX_GEN];

  // mother quantities
  Int_t	    genMomStatus[MAX_GEN];
  // Position
  Float_t   genMomPt[MAX_GEN];
  Float_t   genMomEta[MAX_GEN];
  Float_t   genMomPhi[MAX_GEN];
  Int_t	    genMomPID[MAX_GEN];
  Float_t   genMomBeta[MAX_GEN];
  Float_t   genMomGamma[MAX_GEN];
  Float_t   genMomLxy[MAX_GEN];
  Float_t   genMomLz[MAX_GEN];
  Float_t   genMomLxyz[MAX_GEN];
  Float_t   genMomCTau0[MAX_GEN];
  Float_t   genMom1CTau0 ;
  Float_t   genMom2CTau0 ;
  Float_t   genMom1Lxy ;
  Float_t   genMom2Lxy ;
  Float_t   genMom1Lxyz ;
  Float_t   genMom2Lxyz ;
  Float_t   genMom1Lz ;
  Float_t   genMom2Lz ;
  Float_t   genMom1Pt ;
  Float_t   genMom2Pt ;


  // Sim Vertex Information
  Int_t simVtxN;

  Int_t simVtxProcType[MAX_GEN];
  Int_t simVtxID[MAX_GEN];

  Float_t   simVtxTOF[MAX_GEN];
  Float_t   simVtxX[MAX_GEN];
  Float_t   simVtxY[MAX_GEN];
  Float_t   simVtxZ[MAX_GEN];
  Float_t   simVtxLxy[MAX_GEN];
  Float_t   simVtxLxyz[MAX_GEN];

  /////////////// PRIMARY VERTEX SPECIFIC MEMBERS ////////////////
  Int_t     pvN;
  Float_t   pvMass[MAX_VTX];
  Float_t   pvLxy[MAX_VTX];
  Float_t   pvLxySig[MAX_VTX];
  Float_t   pvLxyz[MAX_VTX];
  Float_t   pvLxyzSig[MAX_VTX];
  Float_t   pvChi2[MAX_VTX];
  Int_t	    pvNTrack[MAX_VTX];  
  Float_t   pvSumPtSq[MAX_VTX];
  // PV  position
  Float_t   pvX[MAX_VTX];
  Float_t   pvY[MAX_VTX];
  Float_t   pvZ[MAX_VTX];
  Float_t   pvZErr[MAX_VTX];
  Float_t   pvYErr[MAX_VTX];
  Float_t   pvXErr[MAX_VTX];

  ///////////////////// TRACK INFORMATION ////////////////////

  Int_t   trCollectionID[MAX_TRACKS];
  // nominal kinematics
  Float_t   trCharge[MAX_TRACKS];
  Float_t   trQOverP[MAX_TRACKS];
  Float_t   trPt[MAX_TRACKS];
  Float_t   trPtError[MAX_TRACKS];
  Float_t   trEta[MAX_TRACKS];
  Float_t   trEtaError[MAX_TRACKS];
  Float_t   trPhi[MAX_TRACKS];
  Float_t   trPhiError[MAX_TRACKS];
  
  // tracking angles
  Float_t   trTheta[MAX_TRACKS];
  Float_t   trThetaError[MAX_TRACKS];
  Float_t   trThetaSig[MAX_TRACKS];
  Float_t   trLambda[MAX_TRACKS];
  Float_t   trLambdaError[MAX_TRACKS];
  Float_t   trLambdaSig[MAX_TRACKS];

  // impact parameter proxies
  Float_t   trDxy[MAX_TRACKS];
  Float_t   trDxyError[MAX_TRACKS];
  Float_t   trDxySig[MAX_TRACKS];
  Float_t   trDz[MAX_TRACKS];
  Float_t   trDzError[MAX_TRACKS];
  Float_t   trDzSig[MAX_TRACKS];
  Float_t   trDsz[MAX_TRACKS];
  Float_t   trDszError[MAX_TRACKS];
  Float_t   trDszSig[MAX_TRACKS];

  // reference point
  Float_t   trRefR2D[MAX_TRACKS];
  Float_t   trRefR3D[MAX_TRACKS];
  Float_t   trRefX[MAX_TRACKS];
  Float_t   trRefY[MAX_TRACKS];
  Float_t   trRefZ[MAX_TRACKS];

  // inner positions
  Float_t   trInnerR2D[MAX_TRACKS];
  Float_t   trInnerR3D[MAX_TRACKS];
  Float_t   trInnerX[MAX_TRACKS];
  Float_t   trInnerY[MAX_TRACKS];
  Float_t   trInnerZ[MAX_TRACKS];
  Float_t   trInnerEta[MAX_TRACKS];
  Float_t   trInnerPhi[MAX_TRACKS];
  // inner momentum
  Float_t   trInnerPt[MAX_TRACKS];
  Float_t   trInnerPx[MAX_TRACKS];
  Float_t   trInnerPy[MAX_TRACKS];
  Float_t   trInnerPz[MAX_TRACKS];
  Float_t   trInnerP[MAX_TRACKS];

  // outer positions
  Float_t   trOuterR2D[MAX_TRACKS];
  Float_t   trOuterR3D[MAX_TRACKS];
  Float_t   trOuterX[MAX_TRACKS];
  Float_t   trOuterY[MAX_TRACKS];
  Float_t   trOuterZ[MAX_TRACKS];
  Float_t   trOuterEta[MAX_TRACKS];
  Float_t   trOuterPhi[MAX_TRACKS];
  Float_t   trOuterRadius[MAX_TRACKS];

  // outer momentum
  Float_t   trOuterPt[MAX_TRACKS];
  Float_t   trOuterPx[MAX_TRACKS];
  Float_t   trOuterPy[MAX_TRACKS];
  Float_t   trOuterPz[MAX_TRACKS];
  Float_t   trOuterP[MAX_TRACKS];

  // quality
  Float_t	trChi2[MAX_TRACKS];
  Float_t	trNDoF[MAX_TRACKS];
  Float_t	trNChi2[MAX_TRACKS];
  Float_t	trValidFraction[MAX_TRACKS];
  Int_t		trNLost[MAX_TRACKS];
  Int_t		trNFound[MAX_TRACKS];
  std::string	trAlgo[MAX_TRACKS];
  Int_t		trAlgoInt[MAX_TRACKS];

  // trigger related
  Int_t		nTrig;
  std::vector<std::string>   triggerNames;
  Int_t         triggerPass[MAX_TRIGGERS];
  Int_t         passDisplacedOR5e33;
  Int_t         passDisplacedOR14e34;
  Int_t         passHTControl; 
  Int_t         passHT200; 
  Int_t         passHT250; 
  Int_t         passHT300; 
  Int_t         passHT350; 
  Int_t         passHT400; 
  Int_t         passDisplaced350_40; 
  Int_t         passDisplaced500_40;
  Int_t         passDisplaced550_40;
  Int_t         passVBFHadronic;
  Int_t         passVBFDispTrack;
  Int_t         passVBFTriple;
  Int_t         passBigOR; 
  Int_t         passHT800;
  Int_t         passPFMET170;
  Int_t         passPFMET170NC;

  ///////////////////HANDLES////////////////////
  // tracks
  edm::Handle<reco::TrackCollection>			gTracks;
  // jets
  edm::Handle<reco::CaloJetCollection>			ak4CaloJets;
  // vertices
  edm::Handle<reco::TrackIPTagInfoCollection>		lifetimeIPTagInfo;
  edm::Handle<reco::SecondaryVertexTagInfoCollection>	secondaryVertexTagInfo;
  edm::Handle<reco::VertexCollection>			secondaryVertices;
  edm::Handle<reco::VertexCollection>			inclusiveVertexCandidates;
  edm::Handle<reco::VertexCollection>			inclusiveSecondaryVertices;
  edm::Handle<reco::VertexCollection>			offlinePrimaryVertices;
  //SIM Compatible 
  edm::Handle<reco::GenParticleCollection>		genParticles;
  edm::Handle<edm::SimVertexContainer>			simVertices;
  edm::Handle<edm::TriggerResults>			triggerResults;

  // trigger related

  //edm::TriggerNames					tNames;
  //edm::TriggerResultsByName				triggerResultsByName;

};

