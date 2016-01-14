// container for Displaced2TrackVertices contained in a combinatorial Cluster
typedef std::vector<Displaced2TrackVertex> DisplacedV0Collection;

class DisplacedCluster {
 public:    
  DisplacedCluster(const Displaced2TrackVertex center_, const reco::Vertex& pv,  const edm::EventSetup& iSetup, const int & debug_) :
  center(center_), selPV(pv), debug(debug_) {
    // recalculated lxy and lxySig
    // for information related to the center, access the center
    // directly
    nV0		       = 0; 
    nTracks	       = 0;
    // linear fit related
    fitChi2	       = 0;
    slope	       = 0;
    interceptPv	       = -999999;   // using t
    interceptZero      = -999999;
    cosAngleToFit      = -10;    
    cosAngleToMomentum = -10;    
    meanX	       = -9999;
    meanY	       = -9999;
    meanZ	       = -9999;
    sumPX	       = -9999;
    sumPY	       = -9999;
    sumPZ	       = -9999;
    
  }

  // decare methods
  void addVertex(Displaced2TrackVertex& vtx);
  void buildClusterQuantities() ;
  void findUniqueTracks();
  void fitVertices();
  bool containsVertex(const Displaced2TrackVertex& vtx);

  // quantities related to the vertex at the center
  // of the cluster
  float fitChi2;
  float slope;
  // intercept relative to the selected PV
  float interceptPv;
  // intercept with respect to (0,0,0)
  float interceptZero; 
  // angle between the line  from pv to center of cluster
  // and line from the fit
  float cosAngleToFit;
  float cosAngleToMomentum; 

  // number of combinatorial vertices contained in the cluster
  int nV0;
  // number of unique tracks contained in the cluster
  int nTracks;      

  // average over the cluster consitutients 
  float meanX, meanY, meanZ;
  float sumPX, sumPY, sumPZ;

  // constructor initialized quantities
  Displaced2TrackVertex 	center;
  const reco::Vertex &		selPV;
  const int			debug;

  // container of vertices
  std::vector<Displaced2TrackVertex> v0s; 
  std::vector<DisplacedTrack> uniq_tracks; 
};

// add a single v0
void DisplacedCluster::addVertex(Displaced2TrackVertex& vtx) { 
  v0s.push_back(vtx);
  nV0++;
}

// add v0 from a collection
/* void DisplacedCluster::addVertices(const DisplacedV0Collection& vertices) {      */
/*   DisplacedV0Collection::const_iterator vtxIter = vertices.begin(); */
/*   for(; vtxIter != vertices.end(); ++vtxIter) v0s.push_back(*vtxIter);        */
/*   nV0++; */
/* } */

// fit the vertices currently contained in the cluster
void DisplacedCluster::fitVertices() {
  if(debug>4) std::cout << "[DEBUG 4] FITTING VERTICES FOR CLUSTER"  << std::endl;
  int nVertices = v0s.size();

  // vector pointing to the center of the cluster
  TVector3 pvToCenter( -1 * (selPV.x() - meanX), -1 * (selPV.y() - meanY), 0);
  TVector3 vertexMomentum(center.px, center.py, 0);

  
  if(nVertices == 0) return;

  // use the direction of the vertex momtnum sum (has redundant tracks)
  cosAngleToMomentum = pvToCenter.Angle(vertexMomentum);

  if(nVertices < 2) return;

  // arrays to fil the tgraph
  float x[nVertices], y[nVertices];
  float xE[nVertices], yE[nVertices];

  // loop over all the vertices and fill the arrays for the TGraph
  // relative to the primary vertex system
  for(int vv = 0; vv < nVertices; ++vv) {    
    //if(fabs(v0s[vv].x) > 1000 || fabs(v0s[vv].y) > 1000  || fabs(selPV.x) > 1000 || fabs(selPV.y) > 1000) {          
    x[vv]  = v0s[vv].x - meanX + .000001 * vv;
    y[vv]  = v0s[vv].y - meanY + .000001 * vv;
    xE[vv] = v0s[vv].xE;
    yE[vv] = v0s[vv].yE;
    if(debug>4) std::cout << "[DEBUG 4] x:"  << x[vv] << " y " << y[vv] << " xE " << xE[vv] << " yE " << yE[vv] <<std::endl;    
  }  

  // build the graph  
  TGraphErrors graph(nVertices, x, y, 0, 0);  
  if(debug > 5) std::cout << "[DEBUG 5] Starting Fit" << std::endl;
  graph.Fit("pol1","Q");
  if(debug > 5) std::cout << "[DEBUG 5] Fit Complete" << std::endl;

  // extract the parameters
  TF1  *fit = graph.GetFunction("pol1"); 
  if (fit == NULL) return;
  float p0  = fit->GetParameter(0);
  if(debug > 5) std::cout << "[DEBUG 5] p0 intercept:" << p0 << std::endl;
  float p1  = fit->GetParameter(1);
  if(debug > 5) std::cout << "[DEBUG 5] p1 slope:" << p1 << std::endl;
  
  fitChi2 = fit->GetChisquare(); 

  // vector pointing along the fit line to the center of the cluster
  TVector3 fitExtrapPlus(1, p0 + p1, 0);
  TVector3 fitExtrapMinus(-1, p0 - p1, 0);
  
  TVector3 vertexSumMomentum(sumPX, sumPY, 0);

  // pick the direction in the direction of the vertex momentum 
  if(vertexMomentum * fitExtrapPlus > 0) {
  // calculate the cosine of the angle between of the angle 
    if(debug>4) std::cout << "[DEBUG 4] Picking Fit Extrapolation Plus"  << std::endl;
    cosAngleToFit = pvToCenter.Angle(fitExtrapPlus);
    interceptPv	  = pvToCenter.Mag() *	sin(pvToCenter.Angle(fitExtrapPlus)); 
  }
  else if (vertexMomentum * fitExtrapMinus > 0) {
    if(debug>4) std::cout << "[DEBUG 4] Picking Fit Extrapolation Minus"  << std::endl;
    cosAngleToFit = pvToCenter.Angle(fitExtrapMinus);
    interceptPv	  = pvToCenter.Mag() *	sin(pvToCenter.Angle(fitExtrapMinus)); 
  }
  else { // this shouldnt happen, but pick the smaller angle 
  if(debug>4) std::cout << "[DEBUG 4] Picking Smaller Cluster Angle with Fit"  << std::endl;
    cosAngleToFit = pvToCenter.Angle(fitExtrapMinus) < pvToCenter.Angle(fitExtrapPlus) ? pvToCenter.Angle(fitExtrapMinus) : pvToCenter.Angle(fitExtrapPlus);
    interceptPv	  = pvToCenter.Angle(fitExtrapMinus) < pvToCenter.Angle(fitExtrapPlus) ? 
      pvToCenter.Mag() * sin(pvToCenter.Angle(fitExtrapMinus)) : pvToCenter.Mag() * sin(pvToCenter.Angle(fitExtrapPlus)); 
  }
}

// since a cluster can contain combinatorial vertices, there are redundant
// tracks. This finds the unique tracks by reference and stores them in the
// uniq_tracks vector
void DisplacedCluster::findUniqueTracks() {
  if(debug>4) std::cout << "[DEBUG 4] Finding Unique Tracks"  << std::endl;
  // loop over all the vertices in the cluster and chceck both tracks
  if(nV0 == 0 || v0s.size() == 0) return;
  if(debug>4) std::cout << "[DEBUG 4] Looping Vertices"  << std::endl;
  DisplacedV0Collection::const_iterator vtxIter = v0s.begin();
  for(; vtxIter != v0s.end(); ++vtxIter) {    
    bool found_match1 = false, found_match2 = false;
    // check if each track in the vertex is already in the uniq track list
    if(!vtxIter->isValid) continue;

    if(debug>4) std::cout << "[DEBUG 4] Looping Uniq Tracks"  << std::endl;
    std::vector<DisplacedTrack>::const_iterator trackIter = uniq_tracks.begin();
    for(; trackIter != uniq_tracks.end(); ++trackIter) {
      if(debug > 4) std::cout << "[DEBUG 4] Checking Size"  << std::endl;
      if(uniq_tracks.size() == 0) break;
      if(debug > 4) std::cout << "[DEBUG 4] Checking Reference"  << std::endl;

      // compare the track references stored in the displacedtrack object
      DisplacedTrack tr1 = vtxIter->track1;
      DisplacedTrack tr2 = vtxIter->track2;
      // compare by reference
      bool  refMatch1 = trackIter->trackRef == tr1.trackRef;
      bool  refMatch2 = trackIter->trackRef == tr2.trackRef;
      // update the matching
      found_match1    = found_match1 || refMatch1;
      found_match2    = found_match2 || refMatch2;

      // if both are found we are done
      if(found_match1 && found_match2) continue;
    }    
    if(debug>4) std::cout << "[DEBUG 4] Finished Looping Uniq Tracks"  << std::endl;
    // if theres no match after check all trakcs in the list add it to the list
    DisplacedTrack tr1 = vtxIter->getTrack1();
    DisplacedTrack tr2 = vtxIter->getTrack2();
    if(debug>4) std::cout << "[DEBUG 4] Pushing Back Tracks "  << tr1.pt << " " << tr2.pt << std::endl;
    if(!found_match1) uniq_tracks.push_back(tr1);
    if(!found_match2) uniq_tracks.push_back(tr2);
  }  // end vertex loop 
  if(debug>4) std::cout << "[DEBUG 4] Finished Looping Vertices"  << std::endl;
  // fil the n tracks variable
  nTracks = uniq_tracks.size();
  if(debug>4) std::cout << "[DEBUG 4] End Unique Track "  << std::endl;
}

// check that a vertex is contained in a cluster
bool DisplacedCluster::containsVertex(const Displaced2TrackVertex& vtx) {
  if(debug>4) std::cout << "[DEBUG 4] Checking if vertex in cluster "  << std::endl;
  if (nV0 == 0) return false;
  bool found_match = false;
  if(debug>4) std::cout << "[DEBUG 4] interating on cluster "  << std::endl;
  DisplacedV0Collection::const_iterator vtxIter = v0s.begin();
  for(; vtxIter != v0s.end(); ++vtxIter) {    
    bool  refMatch1 = (vtx.track1).trackRef == vtxIter->track1.trackRef;
    bool  refMatch2 = (vtx.track2).trackRef == vtxIter->track2.trackRef;
    found_match = found_match || (refMatch1 && refMatch2);
  }

  return found_match;
}

// call to calculate all relevant quantities
// calls the necessary methods in the correct order  
void DisplacedCluster::buildClusterQuantities() {
  if(nV0 == 0 || v0s.size() == 0) return;  

  if(debug>4) std::cout << "[DEBUG 4] Build Cluster Quantities"  << std::endl;
  // find the average position
  meanX = 0, meanY = 0, meanZ = 0;
  if(debug>4) std::cout << "[DEBUG 4] Looping Vertices For mean Calculation"  << std::endl;
  DisplacedV0Collection::const_iterator vtxIter = v0s.begin();
  for(; vtxIter != v0s.end(); ++vtxIter) {    
    meanX += vtxIter->x / float(nV0);
    meanY += vtxIter->y / float(nV0);
    meanZ += vtxIter->z / float(nV0);
    // kinematic sum 
    sumPX += vtxIter->px;
    sumPY += vtxIter->py;
    sumPZ += vtxIter->pz;
  }  

  fitVertices();  
  //findUniqueTracks();
}
