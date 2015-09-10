// container for tracks with all the necessary resources for impact parmameter 
// and hit related calculations 


class DisplacedTrack {
 public:    
 DisplacedTrack(const reco::TrackRef& ref, const reco::Vertex& pv,  const edm::EventSetup& iSetup, const int & debug_) :
  trackRef(ref), selPV(pv), track(*ref), debug(debug_) {

    // set positions for pv
    pvPos      = selPV.position();
    // set kinematics
    pt	       = track.pt();
    px	       = track.px();
    py	       = track.py();
    pz	       = track.pz();
    eta	       = track.eta();
    phi	       = track.phi();
    // quality information
    q	       = track.charge();
    qoverp     = track.qoverp();
    qoverpSig  = track.qoverp() / track.qoverpError();
    chi2       = track.chi2();
    nLostHits  = track.numberOfLostHits();
    nValidHits = track.numberOfValidHits();
    algoInt    = track.algo();    
    // impact related from pv
    dxy	       = track.dxy(pvPos);
    dz	       = track.dz(pvPos);
    dxySig     = track.dxy(pvPos) / track.dxyError();
    dzSig      = track.dz(pvPos) / track.dzError();

    // transient track builder required for vertexing
    edm::ESHandle<TransientTrackBuilder> builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
    transientTrack = builder->build(*ref);        
    // trajectory information for acessing hits
    static GetTrackTrajInfo getTrackTrajInfo;
    std::vector<GetTrackTrajInfo::Result> trajInfo = getTrackTrajInfo.analyze(iSetup, track);   

    ///@@@@@@@@@@@@@@@@@@@@@ HIT INFORMATION @@@@@@@@@@@@@@@@@@@@@@@@@@@
    // if the hits are valid fill the information
    if(trajInfo[0].valid  && trajInfo[1].valid && trajInfo.back().valid ) {
      validHits	    = true;
      // trajectory state on surface
      tsosInnerHit  = trajInfo[0].detTSOS;
      tsosNextHit   = trajInfo[1].detTSOS;            
      tsosLastHit   = trajInfo.back().detTSOS;            
      // get the detector layer from the trajectory information
      const DetLayer &  detLayerInner	 = *(trajInfo[0].detLayer); 
      const DetLayer &  detLayerOuter	 = *(trajInfo.back().detLayer); 
      // detector layers
      GeomDetEnumerators::SubDetector subDetLayerInner = detLayerInner.subDetector();
      GeomDetEnumerators::SubDetector subDetLayerOuter = detLayerOuter.subDetector();
      // check if the track has pixel hits
      hasPixelHits = subDetLayerInner == GeomDetEnumerators::PixelBarrel 
	|| subDetLayerInner == GeomDetEnumerators::PixelEndcap;
     
      // positions
      innerPos	    = tsosInnerHit.globalPosition();
      nextPos	    = tsosNextHit.globalPosition();
      lastPos	    = tsosLastHit.globalPosition();
      // momentum 
      /* innerPos	    = tsosInnerHit.globalMomentum(); */
      /* nextPos	    = tsosNextHit.globalMomentum(); */
      /* lastPos	    = tsosLastHit.globalMomentum(); */
      // coorindates
      innerX   = innerPos.x(), innerY = innerPos.y(), innerZ = innerPos.z();
      nextX    = nextPos.x(), nextY = nextPos.y(), nextZ = nextPos.z();
      lastX    = lastPos.x(), lastY = lastPos.y(), lastZ = innerPos.z();
      // distances
      innerR2D = metric2D(innerX, innerY);
      innerR3D = metric3D(innerX, innerY, innerZ);
      nextR2D  = metric2D(nextX, nextY);
      nextR3D  = metric3D(nextX, nextY, nextZ);
      lastR2D  = metric2D(lastX, lastY);
      lastR3D  = metric3D(lastX, lastY, lastZ);
    }
    else {
      validHits = false;
    }

    ///@@@@@@@@@@@@@@@@@@@@@ IP CALCULATION @@@@@@@@@@@@@@@@@@@@@@@@@@@
    if(transientTrack.isValid()) {
      // separate measurements for the transverse and 3d impact parameter measurements
      std::pair< bool, Measurement1D > ip2dMeasurement = IPTools::absoluteTransverseImpactParameter(transientTrack, selPV);
      std::pair< bool, Measurement1D > ip3dMeasurement = IPTools::absoluteImpactParameter3D(transientTrack, selPV);

      ip2d    = ip2dMeasurement.second.value();
      ip2dSig = ip2dMeasurement.second.significance();
      ip3d    = ip3dMeasurement.second.value();
      ip3dSig = ip3dMeasurement.second.significance();
    }      
  }  
 
  // info
  bool validHits; // has valid hits from hit pattern
    
  // kinematics
  float pt, eta, phi;
  float px, py, pz;
  float q, chi2;
  float qoverp, qoverpSig;
  int algoInt;
  // quality information
  float nLostHits, nValidHits;

  // detector related information
  bool hasPixelHits;

  // impact related
  float dxy, dz;
  float dxySig, dzSig;
  float ip2d, ip3d;    
  float ip2dSig, ip3dSig;

  // hit related
  float innerX, innerY, innerZ;
  float nextX, nextY, nextZ;
  float lastX, lastY, lastZ;
  float innerR2D, innerR3D;
  float nextR2D, nextR3D;
  float lastR2D, lastR3D;  

  // constructor initialized 
  const reco::TrackRef &    trackRef;
  const reco::Vertex &	    selPV;
  const reco::Track &	    track;
  const int		    debug;

  // track associated objects
  reco::Vertex::Point			pvPos;
  reco::TransientTrack			transientTrack ;
  std::vector<GetTrackTrajInfo::Result> trajInfo;

  // position of the hits
  GlobalPoint	innerPos;
  GlobalPoint	nextPos;
  GlobalPoint	lastPos;

  // trajectory state of the inner and next hit
  TrajectoryStateOnSurface tsosInnerHit ; 
  TrajectoryStateOnSurface tsosNextHit ;
  TrajectoryStateOnSurface tsosLastHit ;

  // simplify quadrature calculations
  float metric2D(float x, float y) {
    return std::sqrt(x*x + y*y);
  }
  float metric3D(float x, float y, float z) {
    return std::sqrt(x*x + y*y + z*z);
  }  
};
