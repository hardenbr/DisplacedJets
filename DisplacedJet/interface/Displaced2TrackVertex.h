
class Displaced2TrackVertex {
 public:    
 Displaced2TrackVertex(DisplacedTrack & track1_  , DisplacedTrack & track2_, const reco::Vertex& pv,  const edm::EventSetup& iSetup, const int & debug_) :
  track1(track1_), track2(track2_), selPV(pv), debug(debug_) {
    
    // build the pair
    std::vector<reco::TransientTrack> trackPair;    
    trackPair.push_back(track1.transientTrack);
    trackPair.push_back(track2.transientTrack);

    // fit the pair and check validity
    KalmanVertexFitter fitter(true); // refits the vertex
    if(track1.isValid && track2.isValid) {
      tVertex = fitter.vertex(trackPair);
      isValid = tVertex.isValid();
    }
    else {
      isValid = false;
    }

    // make sure the fit was valid
    if(isValid) {
      vertex = tVertex;		//conversion from tansientVertex to reco::VEretex      
      // quality
      chi2   = vertex.chi2();

      // kinematics
      mass   = vertex.p4().mass();
      charge = track1.q + track2.q;      
      eta    = vertex.p4().eta();
      phi    = vertex.p4().phi();
      pt     = vertex.p4().pt();
      px     = vertex.p4().px(), py = vertex.p4().py(), pz = vertex.p4().pz();
      dr     = reco::deltaR(track1.eta, track1.phi, track2.eta, track2.phi);
      
      // calculate the lambda hypothesis mass      
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> vec1;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> vec2;
      if (track1.pt > track2.pt && charge == 0 ){
	vec1.SetM(0.938272); // proton mass
	vec2.SetM(0.139570); // pion
      }
      else if (track1.pt < track2.pt && charge == 0) {
	vec1.SetM(0.139570); // pi 
	vec2.SetM(0.938272); // proton 	
      }
      else  { 
	vec1.SetM(0);
	vec2.SetM(0);
      }          

      // now set the momentum from the tracks
      //first track
      vec1.SetPx(track1.px);
      vec1.SetPy(track1.py);
      vec1.SetPz(track1.pz);
      //second
      vec2.SetPx(track2.px);
      vec2.SetPy(track2.py);
      vec2.SetPz(track2.pz);
     
      // build the mass candidate
      math::XYZTLorentzVectorD sumLambda;
      sumLambda += vec1;
      sumLambda += vec2;      
      massLambda = sumLambda.mass();      
      
      if(debug > 5) std::cout << "[DEBUG 5] Lambda Mass Hypothesis: " << massLambda << std::endl;

      // position calculations
      vtxPos   = vertex.position();
      x	       = vtxPos.x(), y = vtxPos.y(), z = vtxPos.z();
      xE       = vertex.xError(), yE = vertex.yError(), zE = vertex.zError();
      // distance from the primary vertex selected
      dx       = vtxPos.x() - selPV.x();
      dy       = vtxPos.y() - selPV.y();
      dz       = vtxPos.z() - selPV.z();
      // distances 2D and 3D from pv
      lxy      = metric2D(dx, dy);
      lxyz     = metric3D(dx, dy, dz);
      // total error 
      tot_xE   = metric2D(xE, selPV.xError());
      tot_yE   = metric2D(yE, selPV.yError());
      tot_zE   = metric2D(zE, selPV.zError());
      // 2 and 3d error
      tot_xyE  = metric2D(tot_xE, tot_yE);
      tot_xyzE = metric3D(tot_xE, tot_yE, tot_zE);
      // significances incuding the error
      lxySig   = lxy / tot_xyE;
      lxyzSig  = lxyz / tot_xyzE;
      
      // make the vectors pointing form the vertex to the inner and next hit
      TVector3 vertex_to_innerHit1(track1.innerX - x, track1.innerY - y, track1.innerZ - z);
      TVector3 vertex_to_innerHit2(track2.innerX - x, track2.innerY - y, track2.innerZ - z);
      // next hit after inner hit
      TVector3 vertex_to_nextHit1(track1.nextX - x, track1.nextY - y, track1.nextZ - z);
      TVector3 vertex_to_nextHit2(track2.nextX - x, track2.nextY - y, track2.nextZ - z);
      
      // inner product should be positive if vertex is real
      innerHitBehindVertex1 = (vertex_to_nextHit1 * vertex_to_innerHit1) < 0;
      innerHitBehindVertex2 = (vertex_to_nextHit2 * vertex_to_innerHit2) < 0;

      // tracks inner hit should not be near the vertex (likely a nuclear interaction
      dInnerHitTrack1 = metric3D( x - track1.innerX, y - track1.innerY, z - track1.innerZ);
      dInnerHitTrack2 = metric3D( x - track2.innerX, y - track2.innerY, z - track2.innerZ);
      // tracks inner hits should be separated by the time they reach their inner hit
      dInnerHitTrack1Track2 = metric3D( track1.innerX - track2.innerX, 
					track1.innerY - track2.innerY,
					track1.innerZ - track2.innerZ);          
      // total missing hits
      sumLostHits  = track1.nLostHits  + track2.nLostHits;
      sumValidHits = track1.nValidHits + track2.nValidHits;                
      
      // check kshort window
      isKShort = (fabs(mass - .493) < 0.02) && charge == 0 && track1.pt > 1 && track2.pt > 1 && lxySig > 50;
      isLambda = (fabs(massLambda - 1.1156) < 0.02) && charge == 0 && track1.pt > 1 && track2.pt > 1 && lxySig > 50;
      
    } // valid vertex initialization  
  } // end constructor


  // constructor initialized variables
  DisplacedTrack track1;
  DisplacedTrack track2;
  const reco::Vertex selPV;
  int debug;

  // associated objects
  reco::Vertex		vertex; 
  TransientVertex tVertex;
  math::XYZPoint vtxPos;

  // quality
  bool	isValid;
  float chi2;
  int	sumLostHits;
  int	sumValidHits;

  //kinematics
  int	charge;
  float eta, phi, pt;
  float px, py, pz;
  float mass;
  float massLambda;
  float dr;

  // position
  float x, y, z;
  float xE, yE, zE;  
  float dx, dy, dz;
  // errors combined with pv error
  float tot_xE, tot_yE, tot_zE;
  float tot_xyE, tot_xyzE;

  // positions
  float lxy, lxyz;
  float lxySig, lxyzSig;

  // relative position to tracks
  float dInnerHitTrack1;
  float dInnerHitTrack2;
  float dInnerHitTrack1Track2;

  // particle consistency
  bool isKShort, isLambda;

  // nuclear interaction checks
  bool    innerHitBehindVertex1;
  bool    innerHitBehindVertex2;

  // simplify quadrature calculations
  float metric2D(float x, float y) {
    return std::sqrt(x*x + y*y);
  }
  float metric3D(float x, float y, float z) {
    return std::sqrt(x*x + y*y + z*z);
  }
  
};
