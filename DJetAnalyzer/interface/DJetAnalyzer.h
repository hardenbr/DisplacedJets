class DJetAnalyzer : public edm::EDAnalyzer {

 public:
  explicit  DJetAnalyzer(const edm::ParameterSet&);
  ~DJetAnalyzer();
  
  static void	fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
 private:
  void fillHandles(const edm::Event &);

  // gnerator
  void dumpGenInfo(DisplacedJetEvent&, const reco::GenParticleCollection &); 
  void dumpSimInfo(const edm::SimVertexContainer &);
  void dumpWeights(const edm::Event &);
  void dumpPreSelection(DisplacedJetEvent&);

  // trigger information
  void fillTriggerInfo(const edm::Event &  iEvent, const edm::TriggerResults & trigResults);

  // tree dumping displaced jet quantities
  void	dumpCaloInfo(DisplacedJetEvent&);
  void	dumpSVTagInfo(DisplacedJetEvent&);
  void	dumpIPInfo(DisplacedJetEvent&);
  void	dumpIVFInfo(DisplacedJetEvent&);
  void	dumpDJTags(DisplacedJetEvent&);
  void	dumpV0Info(DisplacedJetEvent&);
  void	dumpPVInfo(DisplacedJetEvent &, const reco::VertexCollection &);
  //tree dumping track quantities
  void	dumpTrackInfo(DisplacedJetEvent&, const reco::TrackCollection &, const int & collectionID, const edm::EventSetup& iSetup);
  void	dumpDisplacedTrackInfo(DisplacedJetEvent& , const edm::EventSetup& iSetup);
  void  dumpRegionalTrackInfo(DisplacedJetEvent& , const edm::EventSetup& iSetup);

  virtual void beginJob() override;
  //virtual void beginRun(edm::Run const& , edm::EventSetup const&) override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  //file configuration tags
  std::string	outputFileName_;
  std::string	jetTreeName_;
  std::string	trackTreeName_;
  std::string	dTrackTreeName_;
  std::string	vertexTreeName_;
  std::string	genTreeName_;

  TFile *outputFile_;
  // analysis to dos
  bool	 doGenMatch_;
  bool	 doSimMatch_;

  bool	 applyEventPreSelection_;
  bool	 applyJetPreSelection_;
  bool   addRegionalTracking_;  // ONLY POSSIBLE OF RE_PROCESSED RAW
  bool   dumpGeneralTracks_;
  bool   dumpDisplacedTracks_;
  bool   dumpRegionalTracks_;

  // keep trees
  bool   writeJetTree_;
  bool   writeV0Tree_;
  bool   writeTrackTree_;
  bool   writeDTrackTree_;
  bool   writeEventTree_;
  bool   writeGenTree_;
  bool   writeVertexTree_;
  // event flags
  bool	 isMC_;
  bool	 isSignalMC_;
  
  //trigger tags
  std::string triggerResultPath_;
  edm::InputTag tag_triggerResults_;

  edm::EDGetTokenT<edm::MergeableCounter > token_eventCounterTotal;
  edm::EDGetTokenT<edm::MergeableCounter > token_eventCounterFiltered;
  edm::InputTag tag_eventCounterTotal_;
  edm::InputTag tag_eventCounterFiltered_;

  //tracking tags
  edm::InputTag tag_generalTracks_;
  edm::InputTag tag_trackIPTagInfoCollection_;
  edm::InputTag tag_lifetimeIPTagInfo_;
  edm::InputTag tag_secondaryVertexTagInfo_; 
  edm::InputTag tag_caloMatchedTracks_;
  edm::InputTag tag_vertexMatchedTracks_;
  // extra tracking collections for online regional tracking
  edm::InputTag tag_regionalTracksIter012_;
  edm::InputTag tag_regionalTracksIter0124_;
  edm::InputTag tag_regionalTracksIter4_;

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

  // alpha Tag
  edm::InputTag tag_alpha_;

  // tag threshold classifications
  float shortTagThresDist, mediumTagThresDist, longTagThresDist;
  int dHTWorkingPoint;
  //cuts
  float cut_jetPt, cut_jetEta;
  float smear_2dip, smear_2dipsig; 
  
  //output related
  TTree*    trackTree_;   
  TTree*    dTrackTree_;   
  TTree*    jetTree_;   
  TTree*    v0Tree_;   
  TTree*    vertexTree_;
  TTree*    genTree_;
  TTree*    eventTree_; 
  TTree*    runStatTree_;
  TTree*    filtCountTree_;

  static const Int_t    SIM_STATUS_CODE_MATCH = 0; 
  static const Int_t    GEN_STATUS_CODE_MATCH = 23; 
  const float		VERTEX_MATCH_METRIC   = 0.05;
  static const Int_t	MAX_TRIGGERS	      = 200;
  static const Int_t	MAX_TRACKS	      = 800;
  static const Int_t	MAX_V0	              = 1000;
  static const Int_t	MAX_JETS	      = 60;
  static const Int_t	MAX_VTX		      = 80;
  static const Int_t	MAX_CAT		      = 4; // max number of tagging categories
  static const Int_t	MAX_GEN		      = 500;
  static const Int_t	FAKE_HIGH_VAL	      = 9999;

  Int_t	debug = 0; 
 
  //bookkeeping
  Int_t run   = -1;
  Int_t lumi  = -1;
  Int_t event = -1;

  int	evNum = 0;
  Int_t jetid = 0;

  Int_t nEventsTotal;
  Int_t nEventsFiltered;

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

  Float_t   caloJetAlpha[MAX_JETS];
  Float_t   caloJetAlphaMax[MAX_JETS];
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
  Int_t   jetNTracksPrompt[MAX_JETS];    
  Int_t   jetNTracksDisp[MAX_JETS];    
  // regional tracking
  Int_t   jetNTracksRegPrompt[MAX_JETS];  // 2DIP > 0.05
  Int_t   jetNTracksRegDisp[MAX_JETS];  // 2DIPSig > 5? 
  // up
  Int_t   jetNTracksRegPromptUp[MAX_JETS];  // 2DIP > 0.05
  Int_t   jetNTracksRegDispUp[MAX_JETS];  // 2DIPSig > 5? 
  // dn
  Int_t   jetNTracksRegPromptDn[MAX_JETS];  // 2DIP > 0.05
  Int_t   jetNTracksRegDispDn[MAX_JETS];  // 2DIPSig > 5? 
  Int_t   jetNTracksReg0124[MAX_JETS];  
  Int_t   jetNTracksReg012[MAX_JETS];  
  Int_t   jetNTracksReg4[MAX_JETS];  
  // jet indexed passing of each requirements
  Int_t   jetPassRegHLTPrompt[MAX_JETS];
  Int_t   jetPassRegHLTDisp[MAX_JETS];
  Int_t   jetPassRegHLTPromptAndDisp[MAX_JETS];
  // up
  Int_t   jetPassRegHLTPromptUp[MAX_JETS];
  Int_t   jetPassRegHLTDispUp[MAX_JETS];
  Int_t   jetPassRegHLTPromptAndDispUp[MAX_JETS];
  // dn
  Int_t   jetPassRegHLTPromptDn[MAX_JETS];
  Int_t   jetPassRegHLTDispDn[MAX_JETS];
  Int_t   jetPassRegHLTPromptAndDispDn[MAX_JETS];



  // aggregate counts of jets passing the specific requirements
  Int_t   nJetsPassRegHLTPrompt; 
  Int_t   nJetsPassRegHLTDisp;
  Int_t   nJetsPassRegHLTPromptAndDisp;
  // up
  Int_t   nJetsPassRegHLTPromptUp; 
  Int_t   nJetsPassRegHLTDispUp;
  Int_t   nJetsPassRegHLTPromptAndDispUp;
  // dn
  Int_t   nJetsPassRegHLTPromptDn; 
  Int_t   nJetsPassRegHLTDispDn;
  Int_t   nJetsPassRegHLTPromptAndDispDn;



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

  // hit related set distributions
  // HIT RELATED
  Float_t jetMedianInnerHitPos[MAX_JETS];
  Float_t jetMedianOuterHitPos[MAX_JETS];
  Float_t jetMeanInnerHitPos[MAX_JETS];
  Float_t jetMeanOuterHitPos[MAX_JETS];
  Float_t jetVarianceInnerHitPos[MAX_JETS];
  Float_t jetVarianceOuterHitPos[MAX_JETS];
  // distributions from inside the pixel layers
  Float_t jetMedianInnerHitPosInPixel[MAX_JETS];
  Float_t jetMedianOuterHitPosInPixel[MAX_JETS];
  Float_t jetMeanInnerHitPosInPixel[MAX_JETS];
  Float_t jetMeanOuterHitPosInPixel[MAX_JETS];
  Float_t jetVarianceInnerHitPosInPixel[MAX_JETS];
  Float_t jetVarianceOuterHitPosInPixel[MAX_JETS];
  // distributions outside the pixel layers                                                                                                                                
  Float_t jetMedianInnerHitPosOutPixel[MAX_JETS];
  Float_t jetMedianOuterHitPosOutPixel[MAX_JETS];
  Float_t jetMeanInnerHitPosOutPixel[MAX_JETS];
  Float_t jetMeanOuterHitPosOutPixel[MAX_JETS];
  Float_t jetVarianceInnerHitPosOutPixel[MAX_JETS];
  Float_t jetVarianceOuterHitPosOutPixel[MAX_JETS];
  // fraction valid hits
  Float_t jetMedianTrackValidHitFrac[MAX_JETS];
  Float_t jetMeanTrackValidHitFrac[MAX_JETS];
  Float_t jetVarianceTrackValidHitFrac [MAX_JETS];
  // track counting
  Float_t jetNTracksNoPixel[MAX_JETS];
  Float_t jetNTracksPixel[MAX_JETS];
  Float_t jetPtSumTracksNoPixel[MAX_JETS];
  Float_t jetPtSumTracksPixel[MAX_JETS];

  // SV information
  Int_t	    jetNSv[MAX_JETS];
  Float_t   jetSvMass[MAX_JETS];
  Float_t   jetSvLxy[MAX_JETS];
  Float_t   jetSvLxySig[MAX_JETS];
  Float_t   jetSvLxyz[MAX_JETS];
  Float_t   jetSvLxyzSig[MAX_JETS];
  Int_t	    jetSvNTrack[MAX_JETS];  //vertex track multiplicty

  // kinematics
  Float_t   jetSvPt[MAX_JETS];
  Float_t   jetSvEta[MAX_JETS];
  Float_t   jetSvPhi[MAX_JETS];
  Float_t   jetSvAngle2D[MAX_JETS];
  Float_t   jetSvAngle3D[MAX_JETS];

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

  // jet indidivudal tagging variables
  Int_t noVertexTag[MAX_JETS];
  Int_t shortTag[MAX_JETS];
  Int_t mediumTag[MAX_JETS];
  Int_t longTag[MAX_JETS];
  Int_t anyTag[MAX_JETS];
  // tight variables
  Int_t noVertexTightTag[MAX_JETS];
  Int_t shortTightTag[MAX_JETS];
  Int_t mediumTightTag[MAX_JETS];
  Int_t longTightTag[MAX_JETS];
  Int_t anyTightTag[MAX_JETS];

  //////////////////////// V0 CANDIDATES ///////////
  //jet info
  Float_t v0JetEta[MAX_V0];
  Float_t v0JetPhi[MAX_V0];
  Float_t v0JetPt[MAX_V0];
  Float_t v0JetNV0[MAX_V0];
  Float_t v0JetNV0AboveP1[MAX_V0];
  Float_t v0JetMedianIPLogSig2D[MAX_V0];
  Float_t v0JetAlphaMax[MAX_V0];
  // cluster size for each vertex
  Int_t   v0JetClusterSize[MAX_V0];
  Int_t   v0JetNJetClusterSize[MAX_V0];
  Int_t   v0InCluster[MAX_V0];
  Int_t   v0InNJetCluster[MAX_V0];
  // info
  Int_t  v0SumValidHits[MAX_V0];
  Int_t  v0SumLostHits[MAX_V0];
  Int_t  nV0;
  Int_t  v0NTracks[MAX_V0];
  Float_t v0isOS[MAX_V0];
  Float_t v0Chi2[MAX_V0];
  Float_t v0NChi2[MAX_V0];
  Int_t   v0IsFake[MAX_V0];
  // kinematics
  Float_t v0Mass[MAX_V0];
  Float_t v0LambdaMass[MAX_V0];
  Float_t v0LambdaMassNoRefit[MAX_V0];
  Float_t v0Pt[MAX_V0];
  Float_t v0Px[MAX_V0];
  Float_t v0Py[MAX_V0];
  Float_t v0Pz[MAX_V0];
  // opening angle
  Float_t v0DR[MAX_V0];
  Float_t v0DRNoRefit[MAX_V0];
  // positions
  Float_t v0Eta[MAX_V0];
  Float_t v0Phi[MAX_V0];
  Float_t v0X[MAX_V0];
  Float_t v0Y[MAX_V0];
  Float_t v0Z[MAX_V0];
  Float_t v0XError[MAX_V0];
  Float_t v0YError[MAX_V0];
  Float_t v0ZError[MAX_V0];
  Float_t v0Lxy[MAX_V0];
  Float_t v0Lxyz[MAX_V0];
  Float_t v0LxySig[MAX_V0];
  Float_t v0LxyzSig[MAX_V0];
  Float_t v0Track1Chi2[MAX_V0];
  Float_t v0Track2Chi2[MAX_V0];
  Float_t v0Track1Pt[MAX_V0];
  Float_t v0Track2Pt[MAX_V0];
  Float_t v0Track1NoRefitPt[MAX_V0];
  Float_t v0Track2NoRefitPt[MAX_V0];
  // dxy
  Float_t v0Track1Dxy[MAX_V0];
  Float_t v0Track2Dxy[MAX_V0];
  Float_t v0Track1DxySig[MAX_V0];
  Float_t v0Track2DxySig[MAX_V0];
  // impact parameter
  // 2d
  Float_t v0Track1IP2D[MAX_V0];
  Float_t v0Track1IP2DSig[MAX_V0];
  Float_t v0Track2IP2D[MAX_V0];
  Float_t v0Track2IP2DSig[MAX_V0];
  // 3d
  Float_t v0Track1IP3D[MAX_V0];
  Float_t v0Track1IP3DSig[MAX_V0];
  Float_t v0Track2IP3D[MAX_V0];
  Float_t v0Track2IP3DSig[MAX_V0];


  ///////////////// V0 Related

  // NUCLEAR INTERACTIONS
  Int_t	    jetOneTrackNuclearCount[MAX_JETS];
  Int_t	    jetTwoTrackNuclearCount[MAX_JETS];
  Int_t	    jetTwoTrackInnerHitFake[MAX_JETS];
  Int_t	    jetVertexNearBPIX1[MAX_JETS];
  Int_t	    jetVertexNearBPIX2[MAX_JETS];
  Int_t	    jetVertexNearBPIX3[MAX_JETS];
  Int_t	    jetVertexNearBPIX[MAX_JETS];
  Int_t	    jetTightNuclear[MAX_JETS];
  Int_t	    jetLooseNuclear[MAX_JETS];  
  Int_t	    jetNV0HitBehindVertex[MAX_JETS];
  Int_t	    jetNV0NoHitBehindVertex[MAX_JETS];
  Int_t	    jetNV0KShort[MAX_JETS];
  Int_t	    jetNV0Lambda[MAX_JETS];
  Int_t	    jetV0HIndex[MAX_JETS];
  // CLUSTERS
  Int_t	    jetV0ClusterSize[MAX_JETS];  
  Float_t   jetV0ClusterLxy[MAX_JETS];
  Float_t   jetV0ClusterLxySig[MAX_JETS];
  Float_t   jetV0ClusterLxyz[MAX_JETS];
  Float_t   jetV0ClusterLxyzSig[MAX_JETS];
  Float_t   jetV0ClusterX[MAX_JETS];
  Float_t   jetV0ClusterY[MAX_JETS];
  Float_t   jetV0ClusterZ[MAX_JETS];
  Float_t   jetV0ClusterChi2[MAX_JETS];
  Float_t   jetV0ClusterIntercept[MAX_JETS];
  Float_t   jetV0ClusterAngle[MAX_JETS];
  Float_t   jetV0ClusterAngleMom[MAX_JETS];
  Int_t	    jetV0ClusterNTracks[MAX_JETS];
  // N JET CLUSTERS 
  Int_t	    jetV0NJetClusterSize[MAX_JETS];  
  Float_t   jetV0NJetClusterLxy[MAX_JETS];
  Float_t   jetV0NJetClusterLxySig[MAX_JETS];
  Float_t   jetV0NJetClusterLxyz[MAX_JETS];
  Float_t   jetV0NJetClusterLxyzSig[MAX_JETS];
  Float_t   jetV0NJetClusterX[MAX_JETS];
  Float_t   jetV0NJetClusterY[MAX_JETS];
  Float_t   jetV0NJetClusterZ[MAX_JETS];
  Float_t   jetV0NJetClusterChi2[MAX_JETS];
  Float_t   jetV0NJetClusterIntercept[MAX_JETS];
  Float_t   jetV0NJetClusterAngle[MAX_JETS];
  Float_t   jetV0NJetClusterAngleMom[MAX_JETS];
  Int_t	    jetV0NJetClusterNTracks[MAX_JETS];

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


  // TRACK ANGLE VARIABLES
  // track angle information
  Float_t   sumTrackPt[MAX_JETS];	      
  // pt weighted
  Float_t   ptSumCosTheta2D[MAX_JETS];	      
  Float_t   ptSumCosTheta3D[MAX_JETS];	      
  Float_t   ptSumCosThetaDet2D[MAX_JETS];	      
  Float_t   ptSumCosThetaDet3D[MAX_JETS];     
  // aboslute sum
  Float_t   sumCosTheta2D[MAX_JETS];	      
  Float_t   sumCosTheta3D[MAX_JETS];	      
  Float_t   sumCosThetaDet2D[MAX_JETS];	      
  Float_t   sumCosThetaDet3D[MAX_JETS];	      
  // mean
  Float_t   meanCosTheta2D[MAX_JETS];	      
  Float_t   meanCosTheta3D[MAX_JETS];	      
  Float_t   meanCosThetaDet2D[MAX_JETS];     
  Float_t   meanCosThetaDet3D[MAX_JETS];      
  // median
  Float_t   medianCosTheta2D[MAX_JETS];	      
  Float_t   medianCosTheta3D[MAX_JETS];	      
  Float_t   medianCosThetaDet2D[MAX_JETS];    
  Float_t   medianCosThetaDet3D[MAX_JETS];    
  // variance
  Float_t   varianceCosTheta2D[MAX_JETS];     
  Float_t   varianceCosTheta3D[MAX_JETS];     
  Float_t   varianceCosThetaDet2D[MAX_JETS];  
  Float_t   varianceCosThetaDet3D[MAX_JETS]; 

  // track vector momentum sum angles
  Float_t   trackSumMomCosTheta2D[MAX_JETS];
  Float_t   trackSumMomCosTheta3D[MAX_JETS];
  Float_t   trackSumMomMag2D[MAX_JETS];
  Float_t   trackSumMomMag3D[MAX_JETS];

  // ip vector sum magnitude
  Float_t ipPosSumMag3D[MAX_JETS];
  Float_t ipPosSumMag2D[MAX_JETS];

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

  // PV matching related
  Int_t	    hasMatchedGenPV;
  Int_t	    selectedPVIsMatched;
  Float_t   pvToGenPVDistance3D;
  Float_t   pvToGenPVDistance2D;
  Float_t   pvToGenPVDistanceZ;
  Float_t   bestPVDistance3D;
  Float_t   bestPVDistance2D;
  Float_t   bestPVDistanceZ;
  Float_t   bestPVX;
  Float_t   bestPVY;
  Float_t   bestPVZ;

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


  ///////////////////// REGIONAL TRACK INFORMATION ////////////

  // kinematics
  Int_t nDtrReg;
  Int_t dtrRegCollection[MAX_TRACKS];
  Int_t dtrRegJetIndex[MAX_TRACKS];
  Float_t dtrRegPt[MAX_TRACKS];
  Float_t dtrRegEta[MAX_TRACKS];
  Float_t dtrRegPhi[MAX_TRACKS];
  // tagging variables
  Float_t dtrRegTheta2D[MAX_TRACKS];
  Float_t dtrReg2DIPSig[MAX_TRACKS];
  Float_t dtrReg2DIP[MAX_TRACKS];
  Float_t dtrRegDxySig[MAX_TRACKS];
  Float_t dtrRegDxy[MAX_TRACKS];
  Int_t   dtrRegTagged[MAX_TRACKS];
  // associated jet information
  Float_t dtrRegJetPt[MAX_TRACKS];
  Float_t dtrRegJetEta[MAX_TRACKS];
  Float_t dtrRegJetPhi[MAX_TRACKS];
  // tracks multiplicities
  Int_t dtrRegJetNTracks[MAX_TRACKS];
  Int_t dtrRegJetNTracksReg012[MAX_TRACKS];
  Int_t dtrRegJetNTracksReg0124[MAX_TRACKS];
  Int_t dtrRegJetNTracksReg4[MAX_TRACKS];
  Int_t dtrRegJetNTracksPrompt[MAX_TRACKS];
  Int_t dtrRegJetNTracksPromptAndDisp[MAX_TRACKS];
  Int_t dtrRegJetNTracksDisp[MAX_TRACKS];
  // variations of the tracking systematics
  // up
  Int_t dtrRegJetNTracksDispUp[MAX_TRACKS];  
  Int_t dtrRegJetNTracksPromptUp[MAX_TRACKS];
  Int_t dtrRegJetNTracksPromptAndDispUp[MAX_TRACKS];
  // dn
  Int_t dtrRegJetNTracksDispDn[MAX_TRACKS];  
  Int_t dtrRegJetNTracksPromptDn[MAX_TRACKS];
  Int_t dtrRegJetNTracksPromptAndDispDn[MAX_TRACKS];

  Int_t dtrRegJetNTaggedTracks[MAX_TRACKS];
  // associated jet tagginginformation
  Float_t dtrRegJetMedian2DIPSig[MAX_TRACKS];
  Float_t dtrRegJetMedianTheta2D[MAX_TRACKS];
  Float_t dtrRegJetAlphaMax[MAX_TRACKS];
  ///////////////////// DISPLACED TRACK INFORMATION ////////////

  // kinematics
  Int_t nDtr;
  Int_t dtrJetIndex[MAX_TRACKS];
  Float_t dtrPt[MAX_TRACKS];
  Float_t dtrEta[MAX_TRACKS];
  Float_t dtrPhi[MAX_TRACKS];
  // tagging variables
  Float_t dtrTheta2D[MAX_TRACKS];
  Float_t dtr2DIPSig[MAX_TRACKS];
  Int_t   dtrTagged[MAX_TRACKS];
  // associated jet information
  Int_t dtrJetNTracks[MAX_TRACKS];
  Int_t dtrJetNTaggedTracks[MAX_TRACKS];


  Float_t dtrJetPt[MAX_TRACKS];
  Float_t dtrJetEta[MAX_TRACKS];
  Float_t dtrJetPhi[MAX_TRACKS];
  Float_t dtrJetAlphaMax[MAX_TRACKS];
  Float_t dtrJetMedian2DIPSig[MAX_TRACKS];
  Float_t dtrJetMedianTheta2D[MAX_TRACKS];
  
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

  // gen sim weights
  Float_t       genWeight;
  Float_t       genWeights[MAX_GEN];
  Float_t       genWeightsRel[MAX_GEN];
  Float_t       genWeightsRMS;
  Int_t         nWeights;
  Float_t       genPU;

  // trigger related
  Int_t		nTrig;
  std::vector<std::string>   triggerNames;
  Int_t         triggerPass[MAX_TRIGGERS];
  Int_t         passDisplacedOR5e33;
  Int_t         passDisplacedOR14e34;
  Int_t         passHTControl; 
  Int_t         passHT200; 
  Int_t         passHT275; 
  Int_t         passHT325; 
  Int_t         passHT425; 
  Int_t         passHT575; 
  Int_t         passDisplaced200_40; 
  Int_t         passDisplaced250_40; 
  Int_t         passDisplaced350_40;
  Int_t         passDisplaced350_40DT;
  Int_t         passDisplaced400_40; 
  Int_t         passDisplaced500_40;
  Int_t         passDisplaced550_40;
  Int_t         passVBFHadronic;
  Int_t         passVBFDispTrack;
  Int_t         passVBFTriple;
  Int_t         passBigOR; 
  Int_t         passHT800;
  Int_t         passPFMET170;
  Int_t         passPFMET170NC;
  Int_t         passMu20;

  ///////////////////HANDLES////////////////////
  // tracks
  edm::Handle<reco::TrackCollection>			           gTracks;
  edm::Handle<reco::JetTracksAssociation::Container>               caloMatchedTracks; 
  edm::Handle<reco::JetTracksAssociation::Container>               vertexMatchedTracks;
  edm::Handle<reco::JetTracksAssociation::Container>               regionalTracksIter012; 
  edm::Handle<reco::JetTracksAssociation::Container>               regionalTracksIter0124; 
  edm::Handle<reco::JetTracksAssociation::Container>               regionalTracksIter4; 

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

  edm::Handle<std::vector<double>>                      alpha;
  // trigger related


  edm::Handle<edm::MergeableCounter> nEventsTotalCounter;  
  edm::Handle<edm::MergeableCounter> nEventsFilteredCounter;
  //edm::TriggerNames					tNames;
  //edm::TriggerResultsByName				triggerResultsByName;

  int ipdf = 1;

};

namespace LHAPDF {
  void   initPDFSet(int nset, const std::string& filename, int member=0);
  int    numberPDF(int nset);
  void   usePDFMember(int nset, int member);
  double xfx(int nset, double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void   extrapolate(bool extrapolate=true);
}
