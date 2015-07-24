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
  std::cout << " call to produce" << std::endl;
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
    trackJetAssigned.push_back(99999);
  }

  // loop over all tracks and associate it to the jet with the lowest dR 
  // but within some minimum value

  std::cout << " loopping tracks" << std::endl;
  for (unsigned t = 0; t < fTracks.size(); ++t) {
    double  lowestDR2Jet	     = 99999;
    unsigned  lowestDR2JetIndex	     = -1;

    // grab the track inner hit
    const reco::Track & const_track = *fTracks[t];
    static GetTrackTrajInfo		    getTrackTrajInfo;
    std::vector<GetTrackTrajInfo::Result>   trajInfo	 = getTrackTrajInfo.analyze(fSetup, const_track);
    const TrajectoryStateOnSurface&	    tsosInnerHit = trajInfo[0].detTSOS;
    const GlobalPoint &			    innerPos     = tsosInnerHit.globalPosition();
   //    const GlobalVector &                    innerMom     = tsosInnerHit.globalMomentum();
    
    std::cout << " loopping jets" << std::endl;

    for (unsigned j = 0; j < fJets.size(); ++j) {
      const reco::Jet & jet = *(fJets[j]); 

      std::cout << " accessing inner hit point" << std::endl;
      // IMPORTANT: correct the jet eta, phi to the inner hit of the track     
      // math::XYZPoint innerHitPoint = math::XYZPoint(innerPos.x(), innerPos.y(), innerPos.z());
      reco::Candidate::Point innerHitPoint = math::XYZPoint(innerPos.x(), innerPos.y(), innerPos.z());
      reco::Candidate::Point jetPoint = math::XYZPoint(jet.vx(), jet.vy(), jet.vz());
      std::cout << "inner point " << innerPos.x() << " " << innerPos.y() << " "  << innerPos.z() << std::endl;
      std::cout << "jet point " << jet.vx() << " " << jet.vy() << " "  << jet.vz() << std::endl;
      std::cout << " modifying jet p4" << std::endl;
      const math::XYZTLorentzVectorD &	innerHitP4 = jet.physicsP4(innerHitPoint, jet, jetPoint);
      double				jetNewEta  = innerHitP4.eta();
      double				jetNewPhi  = innerHitP4.phi();
      double				dR2	   = deltaR2(jetNewEta, jetNewPhi, trackP3s[t].eta(), trackP3s[t].phi());
      std::cout << " jet track DR2 " << dR2  << std::endl;

      // find the closest jet within the threshold
      if (dR2 < mDeltaR2Threshold && dR2 < lowestDR2Jet) {
	lowestDR2JetIndex    = j;
	lowestDR2Jet	     = dR2;
      }
    }
    // assign the track to the jet with the lowest DR2 within the threshold
    std::cout << " lowest index " << lowestDR2JetIndex  << std::endl;
    trackJetAssigned[t] = lowestDR2JetIndex;
  }
  
  std::cout << " building jet association to tracks" << std::endl;
  // loop over the jets again and make an associator for each jet of all the tracks matched
  for (unsigned j = 0; j < fJets.size(); ++j) {
    reco::TrackRefVector assoTracks;      
    for (unsigned t = 0; t < fTracks.size(); ++t) {
      // if the index matches add it to the association
      if (trackJetAssigned[t] == j) {
	assoTracks.push_back(fTracks[t]);
      }
      std::cout << " setting value of association" << std::endl;
      reco::JetTracksAssociation::setValue (fAssociation, fJets[j], assoTracks);
    }
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
