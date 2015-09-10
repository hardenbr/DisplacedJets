
class Displaced2TrackVertex {
 public:
  Displaced2TrackVertex(const reco::Vertex vtx_, const edm::EventSetup& iSetup_) {
    iSetup = iSetup_;
    vtx	   = vtx_;
    vtxPt  = 0;
    vtxEta = 0;
  }
  
  edm::EventSetup& iSetup;
  reco::Vertex vtx;
  reco::Track tr1; 
  reco::Track tr2; 

  // vertex
  float vtxPt, vtxEta, vtxPhi; 
  float vtxRefitPt, vtxRefitEta, vtxRefitPhi; 
  float vtxChi2, vtxNChi2; 
  float vtxLxy, vtxLxyz;
  float vtxX, vtxY, vtxZ;
  float vtxErrorX, vtxErrorY, vtxErrorZ;
  // track 1
  float tr1Pt, tr1Eta, tr1Phi;
  float tr1RefitPt, tr1RefitEta, tr1RefitPhi;
  // track2
  float tr2Pt, tr2Eta, tr2Phi;  
  float tr2RefitPt, tr2RefitEta, tr2RefitPhi;  
}

  
