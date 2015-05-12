import FWCore.ParameterSet.Config as cms

displacedInclusiveVertexFinder  = cms.EDProducer("InclusiveVertexFinder",
       beamSpot = cms.InputTag("offlineBeamSpot"),
       primaryVertices = cms.InputTag("offlinePrimaryVertices"),
       tracks = cms.InputTag("displacedAssocToTracksCaloFace","displacedAssocToTracksCaloFace","ANA"),
#       tracks = cms.InputTag("displacedAssocToTracksCaloFace"), #djets: generalTracks -> displacedAk4JetTracksAssociatorAtCaloFace
       minHits = cms.uint32(0), #djet 8 -> 0 AOD produciton has problems with nhits
       maximumLongitudinalImpactParameter = cms.double(9999), #djet  .3 -> infty
       minPt = cms.double(0.8), #djet .8 -> 1 
       maxNTracks = cms.uint32(50), #djet 30 -> 50

       clusterizer = cms.PSet(
           seedMax3DIPSignificance = cms.double(9999.),
           seedMax3DIPValue = cms.double(9999.),
           seedMin3DIPSignificance = cms.double(5.0), 
           seedMin3DIPValue = cms.double(0.05),
           clusterMaxDistance = cms.double(.5), #500um #dijet .05 -> .5
           clusterMaxSignificance = cms.double(99999), #4.5 sigma  #djet  4.5 ---> infty
           distanceRatio = cms.double(20), # was cluster scale = 1 / density factor =0.05 
           clusterMinAngleCosine = cms.double(-99999), # only forward decays   #djet accept backward decays (unboosted topologies) .5 -> -9999
       ),

       vertexMinAngleCosine = cms.double(.6), # scalar prod direction of tracks and flight dir  #djet accept backward decays .95 -> .6
       vertexMinDLen2DSig = cms.double(2.5), #2.5 sigma 
       vertexMinDLenSig = cms.double(0.5), #0.5 sigma
       fitterSigmacut =  cms.double(3),
       fitterTini = cms.double(256),
       fitterRatio = cms.double(0.25),
       useDirectVertexFitter = cms.bool(True),
       useVertexReco  = cms.bool(True),
       vertexReco = cms.PSet(
               finder = cms.string('avr'),
               primcut = cms.double(1.0),
               seccut = cms.double(3),
               smoothing = cms.bool(True)
       )

)


