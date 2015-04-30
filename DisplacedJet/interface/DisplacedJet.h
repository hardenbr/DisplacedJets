class DisplacedJet {
 public:
  DisplacedJet(const reco::CaloJet & jet, const bool & isMC_) {
    // initialize calo related variables
    caloPt  = jet.pt();
    caloEta = jet.eta();
    caloPhi = jet.phi();

    // store quantites based on detector geometry (rather than vprimary vertex)
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> detP4 = jet.detectorP4();
    detPt  = detP4.pt();
    detEta = detP4.eta();
    detPhi = detP4.phi();
    
    // energy fractions fro very long lived jets
    caloEMEnergyFrac  = jet.emEnergyFraction();
    caloHadEnergyFrac = jet.energyFractionHadronic();
    
    // initialize ip related combination variables
    ipSigLogSum2D    = 0, ipSigLogSum3D = 0;
    ipLogSum2D	     = 0, ipLogSum3D = 0;
    jetDistLogSum    = 0, jetDistLogSum = 0;
    jetDistSigLogSum = 0, jetDistSigLogSum = 0;

    // initalize ip distributional variable
    //mean
    meanIPSig2D	    = 0, meanIPSig3D = 0;
    meanIP2D	    = 0, meanIP3D = 0;
    meanJetDist	    = 0, meanJetDistSig = 0;
    //median
    medianIPSig2D   = 0, medianIPSig3D = 0;  
    medianIP2D      = 0, medianIP3D = 0;  
    medianJetDist   = 0, medianJetDistSig = 0;
    //variance
    varianceIPSig2D = 0, varianceIPSig3D = 0;
    varianceIP2D    = 0, varianceIP3D = 0;
    varianceJetDist = 0, varianceJetDistSig = 0;
    
    isMC    = isMC_;
  }

  // jet info integration
  void addCaloTrackInfo(const reco::TrackRefVector&);
  void addSVTagInfo();
  void addIVFCollection();
  void addIPTagInfo(const reco::TrackIPTagInfo&);

  // jet distribution calculator
  float getJetMedian(const std::vector<float>, bool);
  float getJetMean(const std::vector<float>, bool);
  float getJetVariance(const std::vector<float>, bool);

  //////////////CALO INFORMATION////////////
  bool isMC;
  
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
  float meanJetDist, meanJetDistSig;
  // median
  float medianIPSig2D, medianIPSig3D;  
  float medianIP2D, medianIP3D;
  float medianJetDist, medianJetDistSig;
  // variance
  float varianceIPSig2D, varianceIPSig3D;  
  float varianceIP2D, varianceIP3D;
  float varianceJetDist, varianceJetDistSig;

  //////////////VERTEX VARIABLES//////////////

  // ivf related variables

  // sv related variables

  ////////////////////////////////////////////
  
  // matching quantities 
  bool isGenMatched;
  bool isPVMatched;
  bool hasSimMatchedVertex;

 private: 
  // related vertices
  reco::Vertex selIVF;
  reco::Vertex selSV;
  reco::Vertex selPV;

  // related track collections
  reco::TrackCollection caloMatchedTracks; 
  reco::TrackCollection vertexMatchedTracks; 
  std::vector<reco::btag::TrackIPData> lifetimeIPData; 
  std::vector<float> ip3dVector, ip3dsVector, ip2dVector, ip2dsVector;
  std::vector<float> jetAxisDistVector, jetAxisDistSigVector; 
};

void
DisplacedJet::addCaloTrackInfo(const reco::TrackRefVector & trackRefs) {
    reco::TrackRefVector::const_iterator trackIter = trackRefs.begin();
    for(; trackIter != trackRefs.end(); ++trackIter) {
      caloMatchedTracks.push_back(**trackIter);
    }    
}

void
DisplacedJet::addIPTagInfo(const reco::TrackIPTagInfo & ipTagInfo) {
  // pull the impact parameter data
  lifetimeIPData = ipTagInfo.impactParameterData();

  // loop over each tracks ip information
  std::vector<reco::btag::TrackIPData>::const_iterator ipIter = lifetimeIPData.begin();
  for(; ipIter != lifetimeIPData.end(); ++ipIter)  {
    float ip3d = ipIter->ip3d.value(), ip3ds = ipIter->ip3d.significance();
    float ip2d = ipIter->ip2d.value(), ip2ds = ipIter->ip2d.significance();
    float jetAxisDist = ipIter->distanceToJetAxis.value(), jetAxisDistSig = ipIter->distanceToJetAxis.value();

    // fill the vectors
    ip3dVector.push_back(ip3d);
    ip2dVector.push_back(ip2d);
    ip3dsVector.push_back(ip3ds);
    ip2dsVector.push_back(ip2ds);
    jetAxisDistVector.push_back(jetAxisDist);
    jetAxisDistSigVector.push_back(jetAxisDistSig);

    // ip log sums 3d and 2d
    ipSigLogSum2D    += (ip3d ? log(fabs(ip2ds)) : 0);
    ipSigLogSum3D    += (ip3ds? log(fabs(ip3ds)) : 0);    
    ipLogSum2D	     += (ip2d ? log(fabs(ip2d)) : 0);
    ipLogSum3D	     += (ip3d ? log(fabs(ip3d)) : 0);
    jetDistSigLogSum += (jetAxisDist ? log(jetAxisDistSig) : 0);
    jetDistLogSum    += (jetAxisDist ? log(jetAxisDist) : 0);
  }
  
  // distributional quantities 
  // mean
  meanIPSig2D	     = getJetMean(ip2dsVector, false);
  meanIPSig3D	     = getJetMean(ip3dsVector, false);
  meanIP2D	     = getJetMean(ip2dVector, false);
  meanIP3D	     = getJetMean(ip3dVector, false);
  meanJetDist	     = getJetMean(jetAxisDistVector, false);
  meanJetDistSig     = getJetMean(jetAxisDistSigVector, false);

  // median
  medianIPSig2D	     = getJetMedian(ip2dsVector, false);
  medianIPSig3D	     = getJetMedian(ip3dsVector, false);
  medianIP2D	     = getJetMedian(ip2dVector, false);
  medianIP3D	     = getJetMedian(ip3dVector, false);
  medianJetDist	     = getJetMedian(jetAxisDistVector, false);
  medianJetDistSig   = getJetMedian(jetAxisDistSigVector, false);

  // variance
  varianceIPSig2D    = getJetVariance(ip2dsVector, false);
  varianceIPSig3D    = getJetVariance(ip3dsVector, false);
  varianceIP2D	     = getJetVariance(ip2dVector, false);
  varianceIP3D	     = getJetVariance(ip3dVector, false);
  varianceJetDist    = getJetVariance(jetAxisDistVector, false);
  varianceJetDistSig = getJetVariance(jetAxisDistSigVector, false);
}


float DisplacedJet::getJetMedian(std::vector<float> values, bool is_signed) {   
  
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
  
  return dMedian;
}

float DisplacedJet::getJetMean(std::vector<float> values, bool is_signed) {
  float sum = 0;
  std::vector<float>::const_iterator val = values.begin();
  for (; val != values.end(); ++val) sum += is_signed ? *val : fabs(*val);  

  return sum / float(values.size());
}

float DisplacedJet::getJetVariance(std::vector<float> values, bool is_signed) {
  float sum = 0;
  float mean = getJetMean(values, is_signed);

  std::vector<float>::const_iterator val = values.begin();
  for (; val != values.end(); ++val) {
    sum += (*val - mean) * (*val - mean);
  }

  //safety check
  if(!is_signed ) { assert( sum >= 0 ); }

  return sum / (values.size() - 1.0);
}
