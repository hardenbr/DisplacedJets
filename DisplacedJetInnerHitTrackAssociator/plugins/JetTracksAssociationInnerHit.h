// \class JetTracksAssociationDRVertex
// Associate jets with tracks by simple "delta R" criteria
// Fedor Ratnikov (UMd), Aug. 28, 2007

#ifndef JetTracksAssociationInnerHit_h
#define JetTracksAssociationInnerHit_h

#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "FWCore/Framework/interface/EventSetup.h"


class JetTracksAssociationInnerHit {
 public:
  JetTracksAssociationInnerHit (double fDr);
  ~JetTracksAssociationInnerHit () {}

  void produce (reco::JetTracksAssociation::Container* fAssociation, 
		const std::vector <edm::RefToBase<reco::Jet> >& fJets,
		const std::vector <reco::TrackRef>& fTracks,
		const edm::EventSetup& fSetup) const;
 private:
  /// fidutial dR between track in the vertex and jet's reference direction
  double mDeltaR2Threshold;
};

#endif
