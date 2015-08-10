// Associate jets with tracks by simple "dR" criteria
// Fedor Ratnikov (UMd), Aug. 28, 2007

#include "JetTracksAssociationInnerHit.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/Vector3D.h"

// inputs related to track inner hit position
#include "FWCore/Framework/interface/EventSetup.h"
#include "RecoTracker/DebugTools/interface/GetTrackTrajInfo.h"

JetTracksAssociationInnerHit::JetTracksAssociationInnerHit (double fDr) 
: mDeltaR2Threshold (fDr*fDr){}

void JetTracksAssociationInnerHit::produce (reco::JetTracksAssociation::Container* fAssociation, 
					    const std::vector <edm::RefToBase<reco::Jet> >& fJets,
					    const std::vector <reco::TrackRef>& fTracks,
					    const edm::EventSetup& fSetup) const 
{
  //std::cout << " call to produce" << std::endl;
  // cache tracks kinematics
  std::vector <math::RhoEtaPhiVector> trackP3s;
  // vector of the jet index where the track is assigned
  std::vector<unsigned> trackJetAssigned;
  trackP3s.reserve (fTracks.size());
  trackJetAssigned.reserve (fTracks.size());
  for (unsigned i = 0; i < fTracks.size(); ++i) {
    const reco::Track* track = &*(fTracks[i]);
    trackP3s.push_back (math::RhoEtaPhiVector (track->p(),track->eta(), track->phi())); 
    // default the jet assigned to the track to -1
  }

  // loop over all tracks and associate it to the jet with the lowest dR 
  // but within some minimum value

  //std::cout << " loopping tracks" << std::endl;
  for (unsigned t = 0; t < fTracks.size(); ++t) {
    double  lowestDR2Jet	     = 99999;
    unsigned  lowestDR2JetIndex	     = 99999;

    //std::cout << " gabbing inner hit" << std::endl;
    // grab the track inner hit
    // if(fTracks[t] == NULL)  {
    //   //std::cout << "pointer is null? size: " << fTracks.size() << std::endl;
    //   continue;
    // }
    const reco::Track & const_track = *fTracks[t];
    static GetTrackTrajInfo		    getTrackTrajInfo;
    //std::cout << " getting tracjectory state" << std::endl;
    std::vector<GetTrackTrajInfo::Result>   trajInfo	 = getTrackTrajInfo.analyze(fSetup, const_track);
    //std::cout << "trajInfo size: " << trajInfo.size() << std::endl;
    if(trajInfo.size() == 0) {
      //std::cout << "WARNING: no inner hit? size " << trajInfo.size() << std::endl;
      continue;
    }

    //std::cout << "getting tsos" << std::endl;
    const TrajectoryStateOnSurface& tsosInnerHit	   = trajInfo[0].detTSOS;
    //std::cout << "getting global position" << std::endl;

    if(!tsosInnerHit.isValid()) {
      //std::cout << "hit is not valid....continuing" << std::endl;
      continue;
    }
    //std::cout << "getting inner position" << std::endl;
    const GlobalPoint &		    innerPos		   = tsosInnerHit.globalPosition();
    //const GlobalVector &                    innerMom = tsosInnerHit.globalMomentum();
    
    //std::cout << " looping jets" << std::endl;
    for (unsigned j = 0; j < fJets.size(); ++j) {
      const reco::Jet & jet = *(fJets[j]); 

      //std::cout << " accessing inner hit point" << std::endl;
      // IMPORTANT: correct the jet eta, phi to the inner hit of the track     
      // math::XYZPoint innerHitPoint = math::XYZPoint(innerPos.x(), innerPos.y(), innerPos.z());
      reco::Candidate::Point innerHitPoint = math::XYZPoint(innerPos.x(), innerPos.y(), innerPos.z());
      reco::Candidate::Point jetPoint = math::XYZPoint(jet.vx(), jet.vy(), jet.vz());
      //std::cout << "inner point " << innerPos.x() << " " << innerPos.y() << " "  << innerPos.z() << std::endl;
      //std::cout << "jet point " << jet.vx() << " " << jet.vy() << " "  << jet.vz() << std::endl;
      //std::cout << "modifying jet p4" << std::endl;
      //std::cout << "eta phi  index" << j  << " eta=" << jet.eta() << " phi=" << jet.phi() << std::endl;
      const math::XYZTLorentzVectorD &	innerHitP4 = jet.physicsP4(innerHitPoint, jet, jetPoint);
      double				jetNewEta  = innerHitP4.eta();
      double				jetNewPhi  = innerHitP4.phi();

      // CaloPoint3D<Point> caloPoint(innerHitPoint, const_track.momentum()); 
      // static const Point np(0,0,0);
      // Vector detectorDir = caloPoint.caloPoint() - np;
      // double trackP = const_track.momentum().r();
      // Vector trackP3 = trackP * detectorDir.unit();
      // LorentzVector track_vector(trackP3.x(), trackP3.y(), trackP3.z(), const_track.energy());


      //std::cout << "eta phi index" << j  << " eta=" << innerHitP4.eta() << " phi=" << innerHitP4.phi() << " after " << std::endl;
      //double				dR2	   = deltaR2(jetNewEta, jetNewPhi, trackP3s[t].eta(), trackP3s[t].phi());
      double	dR2    = deltaR2(jetNewEta, jetNewPhi, trackP3s[t].eta(), trackP3s[t].phi());
      //      double	dR2Old = deltaR2(jet.eta(), jet.phi(), trackP3s[t].eta(), trackP3s[t].phi());
      //      double	dRDiff = fabs(dR2 - dR2Old);
      //std::cout << " jet track DR2 old " << dR2Old  << std::endl;
      //std::cout << " jet track DR2 " << dR2  << std::endl;
            
      //      if(dR2Old > mDeltaR2Threshold && dR2 < mDeltaR2Threshold) std::cout << "TRACK  IN: DR diff: " << std::sqrt(dRDiff) << " old DR " << std::sqrt(dR2Old) << " dxy " << const_track.dxy() << std::endl;

      //      if(dR2Old < mDeltaR2Threshold && dR2 > mDeltaR2Threshold) std::cout << "TRACK OUT: DR diff: " << std::sqrt(dRDiff) << " old DR " << std::sqrt(dR2Old) << " dxy " << const_track.dxy() << std::endl;

      // find the closest jet within the threshold
      if (dR2 < mDeltaR2Threshold && dR2 < lowestDR2Jet) {
	lowestDR2JetIndex    = j;
	lowestDR2Jet	     = dR2;
      }
    }
    // assign the track to the jet with the lowest DR2 within the threshold
    //std::cout << " lowest index " << lowestDR2JetIndex  << std::endl;
    trackJetAssigned.push_back(lowestDR2JetIndex);
  }
  
  //std::cout << " building jet association to tracks" << std::endl;
  // loop over the jets again and make an associator for each jet of all the tracks matched
  for (unsigned j = 0; j < fJets.size(); ++j) {
    reco::TrackRefVector assoTracks;      
    for (unsigned t = 0; t < trackJetAssigned.size(); ++t) {
      // if the index matches add it to the association
      if (trackJetAssigned[t] >  fJets.size()) continue; 
      if (trackJetAssigned[t] == j) {
	assoTracks.push_back(fTracks[t]);
      }
      //std::cout << " setting value of association" << std::endl;
    }
    reco::JetTracksAssociation::setValue (fAssociation, fJets[j], assoTracks);
  }

  // //loop on jets and associate
  // for (unsigned j = 0; j < fJets.size(); ++j) {
  //   reco::TrackRefVector assoTracks;
  //   const reco::Jet* jet = &*(fJets[j]); 
  //   double jetEta = jet->eta();
  //   double jetPhi = jet->phi();
  //   for (unsigned t = 0; t < fTracks.size(); ++t) {
  //     double dR2 = deltaR2 (jetEta, jetPhi, trackP3s[t].eta(), trackP3s[t].phi());
  //     if (dR2 < mDeltaR2Threshold)  assoTracks.push_back (fTracks[t]);
  //   }
  //   reco::JetTracksAssociation::setValue (fAssociation, fJets[j], assoTracks);
  // }
}
