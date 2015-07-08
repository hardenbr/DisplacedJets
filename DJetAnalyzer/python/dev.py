import FWCore.ParameterSet.Config as cms

base = '/afs/cern.ch/user/h/hardenbr/2014/LL_DIJET/TRACKING_STUDIES/CMSSW_7_4_6_patch2/src/DisplacedJets/'

# output options (to be appended to the file name outputted)
appendLifetime    = "gun300mm5k"
#appendLifetime   = "600_30mm"
#appendLifetime   = "dsusy500_10mm"
#appendLifetime   = "emerge"
appendBkg         = "470_600"
outputfile_string = "qcd.root"

############ FLAGS #############
#globaltags
#gtag =  "FALL1374_25V4"
gtag = "MCRUN2_74_V9"
#gtag = "PHYS14_25_V1"
# run related
nevents             = 5000
debugLevel          = 0
doedm               = False
# sample related
isSignalMC          = True
isMC                = True
# analysis related
doEventPreSelection = False
doJetPreSelection   = False
doApplyTrigger      = False
# matching
doGenMatch          = True
doSimVtxMatch       = False
# analysis cuts
cut_jetPt           = 40
cut_jetEta          = 2.0
# tag cateogires
shortTagThreshold    = 0.0
mediumTagThreshold  = 10
longTagThreshold    = 30
dHTWorkingPoint     = 2
######### signal input lists ####
input_file_list = None
#input_file_list = '/SignalMCLists/740patch1/mx300_100mm_aod.list'
#input_file_list = '/SignalMCLists/740patch1/mx300_300mm_aod.list'
#input_file_list = '/SignalMCLists/740patch1/mx300_30mm_aod.list'
#input_file_list = '/SignalMCLists/740patch1/mx600_30mm_aod.list'
#input_file_list = '/SignalMCLists/740patch1/mx100_30mm_aod.list'

# other samples
#input_file_list = '/SignalMCLists/DISPLACED_SUSY/displaced_susy.list'
#input_file_list = '/SignalMCLists/DISPLACED_SUSY/displaced_susy500_10.list'
#input_file_list = '/SignalMCLists/EMERGING_JETS/emerging_jets.txt'

# gun samples
input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau30mm.list'
input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau0mm.list'
input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau300mm.list'

# parse the input files to the file list
myfilelist = cms.untracked.vstring()
if input_file_list != None:
   list_from_input_list = open(base + input_file_list, "r") 
   lines = list_from_input_list.readlines()
   stripped_lines = map(lambda x: x.rstrip("\n"), lines)
   for line in stripped_lines:
      myfilelist.extend([line])

#fillter for the file list
if isSignalMC and input_file_list == None :
   myfilelist = cms.untracked.vstring()
   myfilelist.extend(['file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/MC_PRODUCTION/XXTo4J_M-300_CTau_30mm/AODSIM/XXTo4J_M-300_CTau-30mm_reco_102_1_ne5.root',
                      'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/MC_PRODUCTION/XXTo4J_M-300_CTau_30mm/AODSIM/XXTo4J_M-300_CTau-30mm_reco_105_1_1MO.root' ])
if not isSignalMC:
   myfilelist = cms.untracked.vstring()
#   myfilelist.extend(['root://xrootd-cms.infn.it//store/mc/Spring14dr/QCD_Pt-470to600_Tune4C_13TeV_pythia8/AODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/000F701A-95C8-E311-89D0-002618943983.root'])   
#   myfilelist.extend(['root://xrootd-cms.infn.it//store/mc/Fall13dr/QCD_Pt-470to600_Tune4C_13TeV_pythia8/AODSIM/castor_tsg_PU20bx25_POSTLS162_V2-v1/00000/005646AA-EB79-E311-A788-00304867924E.root'])
#   myfilelist.extend(['root://xrootd-cms.infn.it//store/mc/Phys14DR/QCD_Pt-470to600_Tune4C_13TeV_pythia8/AODSIM/AVE20BX25_tsg_castor_PHYS14_25_V3-v1/00000/1E538A06-988E-E411-8F36-0025905B8562.root'])

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")

#standard sequences (from 740x driver command)
process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag = cms.ESSource( "PoolDBESSource",
    globaltag = cms.string( gtag ), 
    RefreshEachRun = cms.untracked.bool( True ),
    RefreshOpenIOVs = cms.untracked.bool( False ),
    toGet = cms.VPSet( 
      cms.PSet(  record = cms.string( "JetCorrectionsRecord" ),
        tag = cms.string( "JetCorrectorParametersCollection_HLT_V1_AK4Calo" ),
        connect = cms.untracked.string( "frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS" ),
        label = cms.untracked.string( "AK4CaloHLT" )
      ),
      cms.PSet(  record = cms.string( "JetCorrectionsRecord" ),
        tag = cms.string( "JetCorrectorParametersCollection_HLT_trk1B_V1_AK4PF" ),
        connect = cms.untracked.string( "frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS" ),
        label = cms.untracked.string( "AK4PFHLT" )
      )
    ),
    DBParameters = cms.PSet( 
      authenticationPath = cms.untracked.string( "." ),
      connectionRetrialTimeOut = cms.untracked.int32( 60 ),
      idleConnectionCleanupPeriod = cms.untracked.int32( 10 ),
      messageLevel = cms.untracked.int32( 0 ),
      enablePoolAutomaticCleanUp = cms.untracked.bool( False ),
      enableConnectionSharing = cms.untracked.bool( True ),
      enableReadOnlySessionOnUpdateConnection = cms.untracked.bool( False ),
      connectionTimeOut = cms.untracked.int32( 0 ),
      connectionRetrialPeriod = cms.untracked.int32( 10 )
    ),
    RefreshAlways = cms.untracked.bool( False ),
    connect = cms.string( "frontier://(proxyurl=http://localhost:3128)(serverurl=http://localhost:8000/FrontierOnProd)(serverurl=http://localhost:8000/FrontierOnProd)(retrieve-ziplevel=0)/CMS_CONDITIONS" ),
    ReconnectEachRun = cms.untracked.bool( True ),
    BlobStreamerName = cms.untracked.string( "TBufferBlobStreamingService" )
)

# Deliver the missing payloads fro the globaltag (TO BE REMOVED)
if 'GlobalTag' in process.__dict__:
    from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
    process.GlobalTag           = customiseGlobalTag(process.GlobalTag, globaltag = gtag)
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_CONDITIONS'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    for pset in process.GlobalTag.toGet.value():
        pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
    # fix for multi-run processing
    process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
    process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nevents))

process.source = cms.Source("PoolSource", fileNames = myfilelist )

################################################################################################

#configure the analyzers
process.analyzerVTX = cms.EDAnalyzer('DJetAnalyzer')
process.analyzerCALO = cms.EDAnalyzer('DJetAnalyzer')

process.analyzerVTX.debugLevel  = cms.untracked.int32(debugLevel)
process.analyzerCALO.debugLevel = cms.untracked.int32(debugLevel)

#output configuration
if isSignalMC:
    process.analyzerVTX.outputFileName = cms.untracked.string('signalVTX%s.root' % appendLifetime)
    process.analyzerCALO.outputFileName = cms.untracked.string('signalCALO%s.root' % appendLifetime)
else:
#    process.analyzerVTX.outputFileName = cms.untracked.string('qcdVTX%s.root' % appendBkg)
#    process.analyzerCALO.outputFileName = cms.untracked.string('qcdCALO%s.root' % appendBkg)
   process.analyzerCALO.outputFileName = cms.untracked.string(outputfile_string)

#tree names
process.analyzerCALO.jetTreeName    = cms.untracked.string('jets')
process.analyzerCALO.trackTreeName  = cms.untracked.string('tracks')
process.analyzerCALO.vertexTreeName = cms.untracked.string('vtx')
process.analyzerCALO.genTreeName    = cms.untracked.string('genp')

process.analyzerVTX.jetTreeName     = cms.untracked.string('jets')
process.analyzerVTX.trackTreeName   = cms.untracked.string('tracks')
process.analyzerVTX.vertexTreeName  = cms.untracked.string('vtx')
process.analyzerVTX.genTreeName     = cms.untracked.string('genp')

# analysis dependent flags
#  applyEventPreSelection_ = iConfig.getUntrackedParameter<bool>("applyEventPreSelection");
#  applyJetPreSelection_   = iConfig.getUntrackedParameter<bool>("applyJetPreSelection");
process.analyzerCALO.applyEventPreSelection = cms.untracked.bool(doEventPreSelection)
process.analyzerCALO.applyJetPreSelection   = cms.untracked.bool(doJetPreSelection)

# MC dependent flags
process.analyzerVTX.isMC        = cms.untracked.bool(isMC)
process.analyzerCALO.isMC       = cms.untracked.bool(isMC)
process.analyzerVTX.isSignalMC  = cms.untracked.bool(isSignalMC)
process.analyzerCALO.isSignalMC = cms.untracked.bool(isSignalMC)
process.analyzerVTX.doGenMatch  = cms.untracked.bool(doGenMatch)
process.analyzerCALO.doGenMatch = cms.untracked.bool(doGenMatch)
process.analyzerVTX.doSimMatch  = cms.untracked.bool(doSimVtxMatch)
process.analyzerCALO.doSimMatch = cms.untracked.bool(doSimVtxMatch)

#tags
process.analyzerVTX.generalTracks  = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerVTX.ak4CaloJets    = cms.untracked.InputTag('ak4CaloJets', '', '')
process.analyzerVTX.genParticles   = cms.untracked.InputTag('genParticles', '', '')
process.analyzerCALO.generalTracks = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerCALO.ak4CaloJets   = cms.untracked.InputTag('ak4CaloJets', '', '')
process.analyzerCALO.genParticles  = cms.untracked.InputTag('genParticles', '', '')

# jet tagging categories
process.analyzerCALO.shortTagThreshold  = cms.untracked.double(shortTagThreshold)
process.analyzerCALO.mediumTagThreshold = cms.untracked.double(mediumTagThreshold)
process.analyzerCALO.longTagThreshold   = cms.untracked.double(longTagThreshold)
process.analyzerCALO.dHTWorkingPoint   = cms.untracked.int32(dHTWorkingPoint)

# vertex matched ip info
process.analyzerVTX.secondaryVertexTagInfo   = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPV', '', 'ANA')
process.analyzerVTX.lifetimeIPTagInfo        = cms.untracked.InputTag('displacedLifetimeTagInfos', '', 'ANA')

# calo matched ip and sv info
process.analyzerCALO.secondaryVertexTagInfo   = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPVCaloFace', '', 'ANA')
process.analyzerCALO.lifetimeIPTagInfo        = cms.untracked.InputTag('displacedLifetimeTagInfosCaloFace', '', 'ANA')

##reconstructed vertex information
process.analyzerCALO.secondaryVertex          = cms.untracked.InputTag('displacedTagsToVerticesCaloFace', '', 'ANA')
process.analyzerCALO.inclusiveVertexCand      = cms.untracked.InputTag('displacedInclusiveVertexFinder', '', 'ANA')
process.analyzerCALO.inclusiveVertexSecondary = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')

process.analyzerVTX.secondaryVertex          = cms.untracked.InputTag('displacedTagsToVertices', '', 'ANA')
process.analyzerVTX.inclusiveVertexCand      = cms.untracked.InputTag('displacedInclusiveVertexFinder', '', 'ANA')
process.analyzerVTX.inclusiveVertexSecondary = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')

# primary vertex
process.analyzerCALO.offlinePrimaryVertices = cms.untracked.InputTag('offlinePrimaryVerticesWithBS')
process.analyzerVTX.offlinePrimaryVertices  = cms.untracked.InputTag('offlinePrimaryVerticesWithBS')

#cuts
process.analyzerVTX.jetPt   = cms.untracked.double(cut_jetPt)
process.analyzerVTX.jetEta  = cms.untracked.double(cut_jetEta)
process.analyzerCALO.jetPt  = cms.untracked.double(cut_jetPt)
process.analyzerCALO.jetEta = cms.untracked.double(cut_jetEta)

# Tags related to the monte carlo
if isMC:
    # Gen Information
    process.analyzerVTX.ak5GenJets    = cms.untracked.InputTag('ak5GenJets', '', '')
    process.analyzerVTX.genMetCalo    = cms.untracked.InputTag('genMetCalo', '', '')
    process.analyzerVTX.genParticles  = cms.untracked.InputTag('genParticles', '', '')
    
    process.analyzerCALO.ak5GenJets   = cms.untracked.InputTag('ak5GenJets', '', '')
    process.analyzerCALO.genMetCalo   = cms.untracked.InputTag('genMetCalo', '', '')
    process.analyzerCALO.genParticles = cms.untracked.InputTag('genParticles', '', '')

    process.analyzerCALO.simVertices = cms.untracked.InputTag('g4SimHits', '', '')
    process.analyzerVTX.simVertices  = cms.untracked.InputTag('g4SimHits', '', '')        

#output anything produced in the ANA process
process.test_output = cms.OutputModule( "PoolOutputModule",
                                        fileName = cms.untracked.string( "/afs/cern.ch/work/h/hardenbr/edmoutput.root" ),
                                        fastCloning = cms.untracked.bool( False ),
                                        dataset = cms.untracked.PSet(filterName = cms.untracked.string( "" ),
                                                                      dataTier = cms.untracked.string( "RAW" )
                                                                      ),
                                        outputCommands = cms.untracked.vstring('keep *_*_*_*'))

# run the displaced jet tags
process.load('DisplacedJets.Configuration.RecoDJTag_cff')
#process.load('DisplacedJets.Configuration.AdditionalPATSequences_cff')
#process.load('DisplacedJets/DisplacedTriggerFilters/displacedTriggers_cff')

#create the main path to run
process.p = cms.Path()

# from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection

# # add the pat jets                                                                                                                                                          
# addJetCollection(
#    process,
#    labelName = 'AK4Calo',
#    jetSource = cms.InputTag('ak4CaloJets'),
#    jetCorrections = ('AK4Calo', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-1'), 
#    btagDiscriminators = ['None'] 
#    )

# #add the pat related sequences
# process.p *= process.makePatMETs

if doApplyTrigger: #apply the triggers and run dj tagging
   process.p *= process.InclusiveTrigger * process.DisplacedTracktrigger * process.djtagging
else: #just run the tagging sequence
   process.p *= process.djtagging

if doedm: #just dump the edm output of the djtagging sequence
    process.btag_output = cms.EndPath(process.test_output)
else: # run the analysis and tree dumper 
    process.p *= process.analyzerCALO


# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs                                                                  
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs                                                           
process = customisePostLS1(process)


process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool( True )
)
