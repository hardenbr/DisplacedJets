import FWCore.ParameterSet.Config as cms

# generic track building tool
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

# iptaginfos with relaxed track cuts
from DisplacedJets.DisplacedImpactParameter.displacedImpactParameterTagInfos_cfi import *

# track counting
from DisplacedJets.DisplacedImpactParameter.trackCountingDJTag_cfi import * #applies track counting using IP signficance of first track

# es modules containing jet tag computers 
from DisplacedJets.DisplacedImpactParameter.displacedTrackCounting2D1st_cfi import * # IP significance of 1st track
from DisplacedJets.DisplacedImpactParameter.promptTrackCountingDJTag_cfi import * # number of prompt tracks


