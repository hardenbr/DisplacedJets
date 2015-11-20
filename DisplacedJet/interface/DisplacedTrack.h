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
    p          = metric3D(px,py,pz);
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
      isValid	    = true;
      // trajectory state on surface
      tsosInnerHit  = trajInfo[0].detTSOS;
      tsosNextHit   = trajInfo[1].detTSOS;            
      tsosLastHit   = trajInfo.back().detTSOS;            
      // get the detector layer from the trajectory information
      const DetLayer &  detLayerInner	 = *(trajInfo[0].detLayer); 
      const DetLayer &  detLayerOuter	 = *(trajInfo.back().detLayer); 
      // detector layers
      GeomDetEnumerators::SubDetector subDetLayerInner = detLayerInner.subDetector();
      //GeomDetEnumerators::SubDetector subDetLayerOuter = detLayerOuter.subDetector();
      // check if the track has pixel hits
      hasPixelHits = subDetLayerInner == GeomDetEnumerators::PixelBarrel 
	|| subDetLayerInner == GeomDetEnumerators::PixelEndcap;

      // positions
      innerPos	    = tsosInnerHit.globalPosition();
      nextPos	    = tsosNextHit.globalPosition();
      lastPos	    = tsosLastHit.globalPosition();
      // momentum 
      innerPosMom   = tsosInnerHit.globalMomentum();      
      lastPosMom    = tsosLastHit.globalMomentum();


      if(debug > 7) std::cout << "[DEBUG 7] pv X " << selPV.x() << " Y " <<  selPV.y() << " Z " << selPV.z() << std::endl;
      // pv vectors of global position
      TVector3 pvVector3D(selPV.x(), selPV.y(), selPV.z());
      TVector3 pvVector2D(selPV.x(), selPV.y(), 0);
      // global position of inner hit of tracks
      if(debug > 7) std::cout << "[DEBUG 7] INNER HIT X " << innerPos.x() << " Y " <<  innerPos.y() << " Z " << innerPos.z() << std::endl;
      TVector3 innerPos3D(innerPos.x(), innerPos.y(), innerPos.z());
      TVector3 innerPos2D(innerPos.x(), innerPos.y(), 0);
      // global position of outer hit
      TVector3 lastPos3D(lastPos.x(), lastPos.y(), lastPos.z());
      TVector3 lastPos2D(lastPos.x(), lastPos.y(), 0);
      // global position of reference point
      /* TVector3 refPoint3D(track.vx(), track.vy(), track.vz()); */
      /* TVector3 refPoint2D(track.vx(), track.vy(), 0); */
      /* // momentum vector at reference point */
      TVector3 refMom3D(px, py, pz); 
      TVector3 refMom2D(px, py, 0); 
      if(debug > 7) std::cout << "[DEBUG 7] Ref MOM X " << refMom3D.x() << " Y " <<  refMom3D.y() << " Z " << refMom3D.z() << std::endl;
      // momentum vector at inner Hit
      TVector3 innerMom3D(innerPosMom.x(), innerPosMom.y(), innerPosMom.z());
      TVector3 innerMom2D(innerPosMom.x(), innerPosMom.y(), 0);
      // momentum vector at outer hit
      TVector3 outerMom3D(lastPosMom.x(), lastPosMom.y(), lastPosMom.z());
      TVector3 outerMom2D(lastPosMom.x(), lastPosMom.y(), 0);

      if(debug > 7) std::cout << "[DEBUG 7] inner MOM X " << innerMom3D.x() << " Y " <<  innerMom3D.y() << " Z " << innerMom3D.z() << std::endl;

      // angles related to the primary vertex 
      angleMomentumAndPVAtInnerHit2D = (-1 * (pvVector2D - innerPos2D)).Angle((refMom2D)); 
      angleMomentumAndPVAtInnerHit3D = (-1 * (pvVector3D - innerPos3D)).Angle((refMom3D));
      angleMomentumAndPVAtOuterHit2D = (-1 * (pvVector2D - innerPos2D)).Angle((innerMom2D)); 
      angleMomentumAndPVAtOuterHit3D = (-1 * (pvVector3D - innerPos3D)).Angle((innerMom3D));

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
      isValid = false;
    }

    ///@@@@@@@@@@@@@@@@@@@@@ IP CALCULATION @@@@@@@@@@@@@@@@@@@@@@@@@@@
    if(transientTrack.isValid()) {
      // separate measurements for the transverse and 3d impact parameter measurements
      std::pair< bool, Measurement1D > ip2dMeasurement = IPTools::absoluteTransverseImpactParameter(transientTrack, selPV);
      std::pair< bool, Measurement1D > ip3dMeasurement = IPTools::absoluteImpactParameter3D(transientTrack, selPV);

      ip2d					       = ip2dMeasurement.second.value();
      ip2dSig					       = ip2dMeasurement.second.significance();
      ip3d					       = ip3dMeasurement.second.value();
      ip3dSig					       = ip3dMeasurement.second.significance();

      if(debug > 6) std::cout << "[DEBUG 6] 2dip " << ip2d << " 2dipsig " << ip2dSig << std::endl;
      if(debug > 6) std::cout << "[DEBUG 6] 3dip " << ip3d << " 3dipsig " << ip3dSig << std::endl;

    }
    else {
      isValid = false;
    }
  }  
 
  // info
  bool isValid; // has valid hits from hit pattern
    
  // kinematics
  float p;
  float pt, eta, phi;
  float px, py, pz;
  float q, chi2;
  float qoverp, qoverpSig;
  int algoInt;
  // quality information
  float nLostHits, nValidHits;

  // detector related information
  bool hasPixelHits;

  // angle related
  float angleMomentumAndPVAtInnerHit2D;
  float angleMomentumAndPVAtOuterHit2D;
  float angleMomentumAndPVAtInnerHit3D;
  float angleMomentumAndPVAtOuterHit3D;

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
  reco::TransientTrack			transientTrack;
  std::vector<GetTrackTrajInfo::Result> trajInfo;

  // position of the hits
  GlobalPoint	innerPos;
  GlobalPoint	nextPos;
  GlobalPoint	lastPos;
  
  GlobalVector	innerPosMom;
  GlobalVector	nextPosMom;
  GlobalVector	lastPosMom;

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
