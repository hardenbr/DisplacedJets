class DisplacedJet {
 public:
  DisplacedJet(const reco::CaloJet & jet_, const reco::Vertex & primaryVertex, const bool & isMC_, const int& jetID_, const int & debug_) {
    debug = debug_;

    jet	  = jet_;
    isMC  = isMC_;
    jetID = jetID_;

    // primary vertex
    selPV = primaryVertex;

    // initialize calo related variables
    caloPt  = jet.pt();
    caloEta = jet.eta();
    caloPhi = jet.phi();

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

    // ivf initialization
    ivfX = 0, ivfY = 0, ivfZ = 0;
    ivfXError = 0, ivfYError = 0, ivfZError = 0;
    ivfLxy = 0, ivfLxyz = 0;
    ivfLxySig = 0, ivfLxyzSig = 0;
    ivfMass = 0;
    ivfNTracks = 0;
    ivfIsGenMatched = false;
    ivfIsSimMatched = false;
    ivfGenVertexMatchMetric = FAKE_HIGH_NUMBER;
    
    // sv initialization
    svX = 0, svY = 0, svZ = 0;
    svXError = 0, svYError = 0, svZError = 0;
    svLxy = 0, svLxyz = 0;
    svLxySig = 0, svLxyzSig = 0;
    svMass = 0;
    svNTracks = 0;
    svIsGenMatched = false;
    svIsSimMatched = false;
    svGenVertexMatchMetric = FAKE_HIGH_NUMBER;    
  }

  // jet info integration
  void addCaloTrackInfo(const reco::TrackRefVector&);
  void addIVFCollection(const reco::VertexCollection&);
  void addIPTagInfo(const reco::TrackIPTagInfo&);
  void addSVTagInfo(const reco::SecondaryVertexTagInfo&);

  // generator matching
  bool doGenCaloJetMatching(const float& ptMatch, const float& dRMatch, const reco::GenParticleCollection& genParticles);
  bool doGenVertexJetMatching(const float& metricThreshold, const reco::GenParticleCollection& genParticles);
  float genMatchMetric(const reco::GenParticle& gen, const reco::Vertex & vertex);

  // jet info extraction
  reco::TrackCollection getCaloMatchedTracks() { return caloMatchedTracks; }
  reco::TrackCollection getVertexMatchedTracks() { return vertexMatchedTracks; }
  reco::Vertex          getIVFVertexSelected() { return selIVF; }
  reco::Vertex          getSVVertex() { return selSV; }
  
  // jet distribution calculator
  float getJetMedian(const std::vector<float>&, bool);
  float getJetMean(const std::vector<float>&, bool);
  float getJetVariance(const std::vector<float>&, bool);

  //////////////CALO INFORMATION////////////
  bool isMC;
  int jetID;
  
  // calo related variables
  float caloPt, caloEta, caloPhi;
  float caloN60, caloN90;  
  float detPt, detEta, detPhi;
  float caloEMEnergyFrac, caloHadEnergyFrac;

  //////////////JET IP VARIABLES/////////////

  // ip related combined variables
  float ipSigLogSum2D, ipSigLogSum3D;
  float ipLogSum2D, ipLogSum3D;  
  float jetDistLogSum, jetDistSigLogSum;
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
  float svLxy, svLxyz;
  float svLxySig, svLxyzSig;
  // qualities
  float svMass;
  float svNTracks;
  // matching
  bool  svIsGenMatched;
  bool  svIsSimMatched;   
  float svGenVertexMatchMetric;
  
  //////////////GEN MATCH VARIABLES///////////

  bool isCaloGenMatched;
  bool isCaloPVGenMatched;
  // matched gen particle kinematics
  float caloGenPt, caloGenEta, caloGenPhi, caloGenMass;


 private: 
  int debug;
  static const int GEN_STATUS_CODE_MATCH = 23; 
  const float FAKE_HIGH_NUMBER = 999999999;
  // calo jet the displaced jet is built upon
  reco::CaloJet jet;

  // related vertices
  reco::Vertex selIVF;
  reco::Vertex selSV;
  reco::Vertex selPV;

  // related track collections
  reco::TrackCollection caloMatchedTracks; 
  reco::TrackCollection vertexMatchedTracks; 
  std::vector<reco::btag::TrackIPData> lifetimeIPData; 
  std::vector<float> ip3dVector, ip3dsVector, ip2dVector, ip2dsVector;
  std::vector<float> ipLog3dVector, ipLog3dsVector, ipLog2dVector, ipLog2dsVector;
  std::vector<float> jetAxisDistVector, jetAxisDistSigVector; 
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

    if (ivfMatch || svMatch ) return true; // if there was any match return true
  } // loop over gen particles

  return false; 
}

bool DisplacedJet::doGenCaloJetMatching(const float& ptMatch, const float& dRMatch, const reco::GenParticleCollection& genParticles) {
  
  reco::GenParticleCollection::const_iterator genIter = genParticles.begin();
  
  for(; genIter != genParticles.end(); ++genIter) {

    int id = genIter->pdgId();
    int st = genIter->status();
  
    if (st != GEN_STATUS_CODE_MATCH) continue;

    // const reco::Candidate * mom = genIter->mother();                                         
    double genPt = genIter->pt(), genEta = genIter->eta(), genPhi = genIter->phi(), genMass = genIter->mass();
    float   dr     = reco::deltaR(genEta, genPhi, caloEta, caloPhi);
    float   dpt    = fabs(caloPt - genPt) / genPt;

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

void DisplacedJet::addCaloTrackInfo(const reco::TrackRefVector & trackRefs) {
  if (debug > 2) std::cout << "[DEBUG] Adding Track Info  " << std::endl;
    reco::TrackRefVector::const_iterator trackIter = trackRefs.begin();
    for(; trackIter != trackRefs.end(); ++trackIter) {
      caloMatchedTracks.push_back(**trackIter);
    }    
}

void DisplacedJet::addIVFCollection(const reco::VertexCollection & vertices) {
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
    
  //build the total error
  float svxE = bestVertex.xError(), svyE = bestVertex.yError(), svzE = bestVertex.zError();
  float pvxE = selPV.xError(), pvyE = selPV.yError(), pvzE = selPV.zError();
  float xE   = std::sqrt(svxE * svxE + pvxE * pvxE), yE = std::sqrt(svyE * svyE + pvyE * pvyE), zE = std::sqrt(svzE * svzE + pvzE * pvzE);

  ivfX = x;
  ivfY = y;
  ivfZ = z;
  ivfXError = svxE;
  ivfYError = svyE;
  ivfZError = svzE;  

  ivfNTracks = bestVertex.nTracks();  
  ivfMass    = bestVertex.p4().mass();    
  ivfLxySig  = std::sqrt( dx * dx + dy * dy ) / std::sqrt(xE * xE + yE * yE);
  ivfLxyzSig = std::sqrt( dx * dx + dy * dy + dz * dz) / std::sqrt(xE * xE + yE * yE + zE * zE);
  ivfLxy     = std::sqrt( dx * dx + dy * dy );
  ivfLxyz    = std::sqrt( dx * dx + dy * dy + dz * dz );

  // matching score
  ivfMatchingScore = bestVertexScore;
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
    int nTracks = vtx.nTracks();

    // take the vertex with the most tracks, tie breaker is the sum pt of vertex
    if ( (nTracks > mostTracks) || (nTracks == mostTracks && pt > tieBreaker) ){
      mostTracks = nTracks;
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
    float ip3d = ipIter->ip3d.value(), ip3ds = ipIter->ip3d.significance();
    float ip2d = ipIter->ip2d.value(), ip2ds = ipIter->ip2d.significance();
    float jetAxisDist = ipIter->distanceToJetAxis.value(), jetAxisDistSig = ipIter->distanceToJetAxis.significance();

    // fill the vectors
    ip3dVector.push_back(ip3d);
    ip2dVector.push_back(ip2d);
    ip3dsVector.push_back(ip3ds);
    ip2dsVector.push_back(ip2ds);

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
    if (debug > 4) std::cout << "[DEBUG 4] jetdistsig: " << jetAxisDistSig << std::endl;
    ipSigLogSum2D    += (ip3d ? log(fabs(ip2ds)) : 0);
    ipSigLogSum3D    += (ip3ds? log(fabs(ip3ds)) : 0);    
    ipLogSum2D	     += (ip2d ? log(fabs(ip2d)) : 0);
    ipLogSum3D	     += (ip3d ? log(fabs(ip3d)) : 0);
    jetDistSigLogSum += fabs(jetAxisDistSig);
    jetDistLogSum    += fabs(jetAxisDist);
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
  medianIPSig2D	   = getJetMedian(ip2dsVector, false);
  medianIPSig3D	   = getJetMedian(ip3dsVector, false);
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
    if (debug > 5) std::cout << "[DEBUG 5] value:" << *val << std::endl;
    if (debug > 5) std::cout << "[DEBUG 5] mean: " << mean << std::endl;
    if (debug > 5) std::cout << "[DEBUG 5] variance element: " << temp << std::endl;
    sum += temp;
  }

  //safety check
  if(!is_signed ) { assert( sum >= 0 ); }
  float variance = std::sqrt(sum / float((values.size() - 1.0)));

  if (debug > 4) std::cout << "[DEBUG] jet variance:  " << variance << std::endl;
  return variance;
}
