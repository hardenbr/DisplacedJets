import FWCore.ParameterSet.Config as cms

displacedOfflinePrimaryVertices = cms.EDProducer("PrimaryVertexProducer", #name change

    verbose = cms.untracked.bool(False),
    TrackLabel = cms.InputTag("generalTracks"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                        
    TkFilterParameters = cms.PSet(
        algorithm=cms.string('filter'),
        maxNormalizedChi2 = cms.double(20.0),
        minPixelLayersWithHits=cms.int32(2), #displaced 2->0
        minSiliconLayersWithHits = cms.int32(5), #displaced 5->0
        maxD0Significance = cms.double(5.0),  #displaced 5->infinity
        minPt = cms.double(0.0), #displaced 0 -> 1
        trackQuality = cms.string("any")
    ),

    TkClusParameters = cms.PSet(
        algorithm   = cms.string("DA"),
        TkDAClusParameters = cms.PSet(
            coolingFactor = cms.double(0.6),  #  moderate annealing speed
            Tmin = cms.double(4.),            #  end of annealing
            vertexSize = cms.double(0.01),    #  ~ resolution / sqrt(Tmin)
            d0CutOff = cms.double(3.), # displaced 3 -> infinity      # downweight high IP tracks
            dzCutOff = cms.double(4.)  # displaced 4 -> infinity      # outlier rejection after freeze-out (T<Tmin)
        )
    ),

    vertexCollections = cms.VPSet(
     [cms.PSet(label=cms.string(""),
               algorithm=cms.string("AdaptiveVertexFitter"),
               minNdof=cms.double(0.0),
               useBeamConstraint = cms.bool(False),
               maxDistanceToBeam = cms.double(1.0)   #displaced 1.0 -> infity
               ),
      cms.PSet(label=cms.string("WithBS"),
               algorithm = cms.string('AdaptiveVertexFitter'), 
               minNdof=cms.double(2.0),
               useBeamConstraint = cms.bool(True),
               maxDistanceToBeam = cms.double(1.0)  #displaced 1.0 -> infty
               )
      ]
    )
)

# These lines should be uncommented only for the gcc4X builds, with X>=6
displacedOfflinePrimaryVertices.TkClusParameters.algorithm = cms.string("DA_vect" )
displacedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.use_vdt = cms.untracked.bool( True )

