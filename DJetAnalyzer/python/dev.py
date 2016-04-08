
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

# output directories
base                = '/afs/cern.ch/user/h/hardenbr/2014/LL_DIJET/TRACKING_STUDIES/CMSSW_7_6_3/src/DisplacedJets/'
outputDir           = "/afs/cern.ch/work/h/hardenbr/2015/DIJET/DJANALYSIS/JetHT"
#outputDir           = ""

# output options (to be appended to the file name outputted)
appendSignal       = ""
appendData         = ""
appendBkg          = ""
############ FLAGS #############
debugLevel         = 0
reportEveryNEvents = 100
isSignalMC         = False
isMC               = False or isSignalMC
isData             = not isMC
doedm              = False
nevents            = 1000

#-------------- globaltags
#gtag              = "74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0" #spring 15 25ns
#gtag              = "74X_mcRun2_asymptotic_realisticBS_v1" #realistic beam spot for 74x
#gtag              = "FALL1374_25V4"
#gtag              = "PHYS14_25_V1"
#gtag              = "MCRUN2_74_V9" #guns
#gtag              = "74X_dataRun2_Prompt_v2"  #data
gtag               = "76X_dataRun2_v15"
#gtag               = "76X_mcRun2_asymptotic_v12"

## -------------json
#JSON               = 'json_DCSONLY_Run2015B.txt'
#JSON               = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'
#JSON               = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
#JSON               = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258714_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#silver json
JSON                = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt'

#--------------trigger
trigger_process     = "HLT" 
#--------------analysis todos
#$$$$$$$$$trigger fixes
is76XTriggers       = True
doApplyTrigger      = False if isData else False #not isSignalMC
##mu related
doApplySingleMu     = False
isOnlyMu            = False
###
dispTriggersOnly    = False
removePFHT800       = False #NOT USED
###################
doEventPreSelection = False
doJetPreSelection   = False
dumpGeneralTracks   = False
dumpDisplacedTracks = True
#regional tracking
addRegionalTracking = True #regional tracking will be added into the jetTree
dumpRegionalTracks  = True #ONLY FOR RAWAOD++
# trees to write
writeTrackTree      = False
writeDTrackTree     = True  #displaced tracks
writeV0Tree         = False
writeEventTree      = True
writeJetTree        = True
writeVertexTree     = False
writeGenTree        = isSignalMC 
#----------- matching
doGenMatch          = isMC #and isSignalMC
doSimVtxMatch       = False
#----------- analysis cuts
cut_jetPt           = 60
cut_jetEta          = 2.0
#----------- tag categories
shortTagThreshold   = 0.0
mediumTagThreshold  = 0.0
longTagThreshold    = 30
dHTWorkingPoint     = 2

######### input lists #########
input_file_list     = None
#----xx4j samples
#input_file_list    = '/SignalMCLists/740patch1/mx300_1000mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx300_300mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx300_100mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx300_30mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx300_1mm_aod.list'
#input_file_list     = '/SignalMCLists/740patch1/mx300_0p1mm_aod.list'

#---- xx4j mass variants 
#input_file_list    = '/SignalMCLists/740patch1/mx600_30mm_aod.list'
#input_file_list    = '/SignalMCLists/740patch1/mx100_30mm_aod.list'

#----other samples
#input_file_list = '/SignalMCLists/DISPLACED_SUSY/displaced_susy.list'
#input_file_list = '/SignalMCLists/DISPLACED_SUSY/displaced_susy500_10.list'
#input_file_list = '/SignalMCLists/EMERGING_JETS/emerging_jets.txt'

#-----gun samples
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau30mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau0mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau300mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau10mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_flat_1mm_1000mm.list'
#input_file_list = 'SignalMCLists/DIJET_GUN/dijet_gun_m300_ctau0mm_bbar.list'

#---RAWAOD++ samples
input_file_list = 'RAW_CONTENT_QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8.list'

#----data samples
input_file_list = 'DataSampleLists/JetHT_AVR5_RAW_CONTENT.list'
#input_file_list  = 'DataSampleLists/PD_DisplacedJet_Run2015D_Sept21.txt'
#input_file_list = 'DataSampleLists/PD_DisplacedJet_Jul17AOD.txt'
#input_file_list = 'DataSampleLists/PD_JetHT_Jul17AOD.txt'
#input_file_list = 'DataSampleLists/PD_HTMHT_Jul17AOD.txt'
#input_file_list = 'DataSampleLists/PD_SingleMuon_Jul23AOD.txt'

# parse the input files to the file list
myfilelist = cms.untracked.vstring()
if input_file_list != None:
   list_from_input_list = open(base + input_file_list, "r") 
   lines = list_from_input_list.readlines()
   stripped_lines = map(lambda x: x.rstrip("\n"), lines)
   for line in stripped_lines:
      myfilelist.extend([line])

#filter for the file list
if isSignalMC and input_file_list == None :
   myfilelist = cms.untracked.vstring()
#   myfilelist.extend(["file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/PRODUCTION_74XGSDR/XX4J_M-600_CTAU30_NOPU/CRAB_PrivateMC/crab_XX4J_M-600_30mm_NOPU/600_30_nopu.root"])
   myfilelist.extend(["file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/PRODUCTION_74XGSDR/XX4J_M-600_CTAU30_NOPU/CRAB_PrivateMC/crab_XX4J_M-600_30mm_NOPU/160223_103026/0000/xx4j_600GeV_30mm_nopu_1kEv.root"])
#   myfilelist.extend(["/store/mc/RunIIFall15DR76/XXTo4J_M-300_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/020A2CB2-AFC4-E511-8C17-002590A2CCF2.root"])
#   myfilelist.extend(["/store/mc/RunIIFall15DR76/XXTo4J_M-500_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/00160C77-CCC6-E511-AB01-0025904C7DF0.root"])
   #print "NO SIGNAL INPUT?" 
#   exit(1)
#   myfilelist.extend(['file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/MC_PRODUCTION/XXTo4J_M-300_CTau_30mm/AODSIM/XXTo4J_M-300_CTau-30mm_reco_102_1_ne5.root',
 #                     'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/MC_PRODUCTION/XXTo4J_M-300_CTau_30mm/AODSIM/XXTo4J_M-300_CTau-30mm_reco_105_1_1MO.root' ])
if not isSignalMC and input_file_list == None and not isData:
   myfilelist = cms.untracked.vstring()
   #   myfilelist.extend(['file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/data/Run2015B/DisplacedJet/AOD/PromptReco-v1/000/251/562/00000/F6834634-9A2A-E511-9F6F-02163E012402.root'])   
#   qcd_files = ['/store/mc/RunIIFall15DR76/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/30000/00348C6E-599F-E511-B51D-02163E00F4BF.root'] 
#   qcd_files = ['file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/RAW_CONTENT_FILES/QCD_MAR31_REGIONALTRACKING_NOTRIGSUM_NOLOCALITY_FIXOUT/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/crab_qcd470_600/160331_115244/0000/edmoutput_18.root']
   qcd_files = ['file:/afs/cern.ch/user/h/hardenbr//eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/RAW_CONTENT_FILES/QCD_MAR31_REGIONALTRACKING_NOTRIGSUM_NOLOCALITY_FIXOUT/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/crab_qcd170_300/160331_115231/0000/edmoutput_55.root']
   # qcd_files = ['/store/mc/RunIIFall15DR76/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/40000/0AD1B83E-DCA0-E511-AA14-0025905A60F2.root',
   #              '/store/mc/RunIIFall15DR76/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/40000/18EEA885-D9A0-E511-9529-0025905A6138.root',
   #              '/store/mc/RunIIFall15DR76/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/40000/1AD1A49F-D9A0-E511-9484-0CC47A4D7618.root',
   #              '/store/mc/RunIIFall15DR76/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/40000/1CBB563A-F1A0-E511-A729-002590D9D89C.root']
   for ff in qcd_files: myfilelist.extend([ff])   
   # qcd_files = ['/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/5E794759-B9FB-E411-99F3-001E67397E90.root',
   #               '/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/5EADB739-17FB-E411-9D85-0025905B858E.root',
   #               '/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/5EC9A83E-FBFA-E411-87B4-002618943925.root',
   #               '/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/62F42BF1-EBFA-E411-8A7E-002590D0B0B6.root',
   #               '/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/6468E4FD-19FB-E411-83B5-0025905A6090.root',
   #               '/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/66AC1E2F-2CFB-E411-8596-002618FDA259.root',
   #               '/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/6A602B0D-0CFB-E411-9A94-00248C55CC3C.root',
   #               '/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/6A65D77E-07FB-E411-95C1-002590D0B002.root',
   #               '/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/70AD84BC-DCFB-E411-BF96-0025904CDDF8.root',
   #               '/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/7444B6BD-DCFB-E411-BD11-0025905C96E8.root']
                
#   myfilelist.extend(['/store/mc/RunIISpring15DR74/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/00D76158-CCFC-E411-89EA-AC853DA06B56.root'])
if isData and input_file_list == None:
   these_files = ['file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/RAW_CONTENT_FILES/JetHT_AVR5_RAW_CONTENT/JetHT/crab_JetHT/160405_130026/JetHT_0000_head200.root']
   for mfile in these_files: 
      print mfile
      myfilelist.extend([mfile])
#   myfilelist.extend(['file:pickevents.root'])

process = cms.Process("ANA", eras.Run2_25ns)

process.load("FWCore.MessageService.MessageLogger_cfi")
#standard sequences (from 740x driver command)
process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = reportEveryNEvents
process.MessageLogger.warnings.suppressInfo = cms.untracked.vstring()
supressWarnings = ["TwoTrackMinimumDistance", "displacedInclusiveVertexFinder","inclusiveVertexFinder","InclusiveVertexFinder","BasicTrajectoryState"]
for warning in supressWarnings: process.MessageLogger.suppressWarning.extend([warning])

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff') #new for navigation
#process.load('Configuration.StandardSequences.GeometryExtended_cff') #new for navigation
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff') #new for navigation
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
#global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.load('JetMETCorrections.Configuration.JetCorrectorsForReco_cff')
#process.load('JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('JetMETCorrections.Configuration.CorrectedJetProducers_cff')

process.GlobalTag.globaltag = gtag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nevents))
process.source = cms.Source("PoolSource", fileNames = myfilelist )

# add a JSON for data
import FWCore.PythonUtilities.LumiList as LumiList
if isData:
   process.source.lumisToProcess = LumiList.LumiList(filename = JSON).getVLuminosityBlockRange()

################################################################################################

#jet corrections taken from the global tag
process.correctJets = cms.Sequence( process.ak4CaloL2L3CorrectorChain * process.ak4CaloJetsL2L3)

#configure the analyzers
process.analyzerVTX  = cms.EDAnalyzer('DJetAnalyzer')
process.analyzerCALO = cms.EDAnalyzer('DJetAnalyzer')

process.analyzerVTX.debugLevel  = cms.untracked.int32(debugLevel)
process.analyzerCALO.debugLevel = cms.untracked.int32(debugLevel)

#output configuration
if isSignalMC:
    process.analyzerVTX.outputFileName = cms.untracked.string('signalVTX%s.root' % appendSignal)
    process.analyzerCALO.outputFileName = cms.untracked.string('%ssignal%s.root' % (outputDir, appendSignal))
elif not isSignalMC and isMC:
#    process.analyzerVTX.outputFileName = cms.untracked.string('qcdVTX%s.root' % appendBkg)
#    process.analyzerCALO.outputFileName = cms.untracked.string('qcdCALO%s.root' % appendBkg)
   process.analyzerCALO.outputFileName = cms.untracked.string("%sdjana.root" % outputDir)
else:
   process.analyzerVTX.outputFileName  = cms.untracked.string('dataVTX%s.root' % appendData)
   process.analyzerCALO.outputFileName = cms.untracked.string('%sdata%s.root' % (outputDir, appendData))   


#tree names
process.analyzerCALO.jetTreeName    = cms.untracked.string('jets')
process.analyzerCALO.trackTreeName  = cms.untracked.string('tracks')
process.analyzerCALO.dTrackTreeName  = cms.untracked.string('dtracks')
process.analyzerCALO.vertexTreeName = cms.untracked.string('vtx')
process.analyzerCALO.genTreeName    = cms.untracked.string('genp')

process.analyzerVTX.jetTreeName     = cms.untracked.string('jets')
process.analyzerVTX.trackTreeName   = cms.untracked.string('tracks')
process.analyzerVTX.vertexTreeName  = cms.untracked.string('vtx')
process.analyzerVTX.genTreeName     = cms.untracked.string('genp')

# analysis dependent flags
#  applyEventPreSelection_                  = iConfig.getUntrackedParameter<bool>("applyEventPreSelection");
#  applyJetPreSelection_                    = iConfig.getUntrackedParameter<bool>("applyJetPreSelection");
process.analyzerCALO.applyEventPreSelection = cms.untracked.bool(doEventPreSelection)
process.analyzerCALO.applyJetPreSelection   = cms.untracked.bool(doJetPreSelection)
process.analyzerCALO.dumpGeneralTracks      = cms.untracked.bool(dumpGeneralTracks)
process.analyzerCALO.dumpDisplacedTracks    = cms.untracked.bool(dumpDisplacedTracks)
process.analyzerCALO.dumpRegionalTracks     = cms.untracked.bool(dumpRegionalTracks)
process.analyzerCALO.addRegionalTracking    = cms.untracked.bool(addRegionalTracking)


# what to write out 
process.analyzerCALO.writeTrackTree  = cms.untracked.bool(writeTrackTree)
process.analyzerCALO.writeDTrackTree  = cms.untracked.bool(writeDTrackTree)
process.analyzerCALO.writeEventTree  = cms.untracked.bool(writeEventTree)
process.analyzerCALO.writeJetTree    = cms.untracked.bool(writeJetTree)
process.analyzerCALO.writeV0Tree    = cms.untracked.bool(writeV0Tree)
process.analyzerCALO.writeVertexTree = cms.untracked.bool(writeVertexTree)
process.analyzerCALO.writeGenTree    = cms.untracked.bool(writeGenTree)

# MC dependent flags
process.analyzerVTX.isMC        = cms.untracked.bool(isMC)
process.analyzerCALO.isMC       = cms.untracked.bool(isMC)
process.analyzerVTX.isSignalMC  = cms.untracked.bool(isSignalMC)
process.analyzerCALO.isSignalMC = cms.untracked.bool(isSignalMC)
process.analyzerVTX.doGenMatch  = cms.untracked.bool(doGenMatch)
process.analyzerCALO.doGenMatch = cms.untracked.bool(doGenMatch)
process.analyzerVTX.doSimMatch  = cms.untracked.bool(doSimVtxMatch)
process.analyzerCALO.doSimMatch = cms.untracked.bool(doSimVtxMatch)

#trigger
process.analyzerVTX.triggerResultPath  = cms.untracked.string(trigger_process)
process.analyzerCALO.triggerResultPath = cms.untracked.string(trigger_process)
process.analyzerVTX.triggerResults     = cms.untracked.InputTag('TriggerResults', '', '')
process.analyzerCALO.triggerResults    = cms.untracked.InputTag('TriggerResults', '', trigger_process)

# collection tags
process.analyzerVTX.generalTracks          = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerVTX.ak4CaloJets            = cms.untracked.InputTag('ak4CaloJets', '', '')
process.analyzerVTX.genParticles           = cms.untracked.InputTag('genParticles', '', '')

# regional tracking rom the HLT
process.analyzerVTX.regionalTracksIter012         = cms.untracked.InputTag('displacedAk4JetTracksAssociatorAtVertexRegionalIter012', '', '')
process.analyzerVTX.regionalTracksIter0124        = cms.untracked.InputTag('displacedAk4JetTracksAssociatorAtVertexRegionalIter0124', '', '')
process.analyzerVTX.regionalTracksIter4           = cms.untracked.InputTag('displacedAk4JetTracksAssociatorAtVertexRegionalIter4', '', '')
# process.analyzerVTX.regionalTracksIter012  = cms.untracked.InputTag('hltIter2MergedForBTag', '', '')
# process.analyzerVTX.regionalTracksIter0124 = cms.untracked.InputTag('hltIter4MergedWithIter012DisplacedJets', '', '')
# process.analyzerVTX.regionalTracksIter4    = cms.untracked.InputTag('hltDisplacedhltIter4PFlowTrackSelectionHighPurity', '', '')

process.analyzerCALO.generalTracks                 = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerCALO.ak4CaloJets                   = cms.untracked.InputTag('ak4CaloJetsL2L3', '', '')
process.analyzerCALO.genParticles                  = cms.untracked.InputTag('genParticles', '', '')
process.analyzerCALO.caloMatchedTrackAssociation   = cms.untracked.InputTag('displacedAk4JetTracksAssociatorAtCaloFace','','ANA')
process.analyzerCALO.vertexMatchedTrackAssociation = cms.untracked.InputTag('displacedAk4JetTracksAssociatorAtVertex','','ANA')
#regional tracking
process.analyzerCALO.regionalTracksIter012         = cms.untracked.InputTag('displacedAk4JetTracksAssociatorAtVertexRegionalIter012', '', '')
process.analyzerCALO.regionalTracksIter0124        = cms.untracked.InputTag('displacedAk4JetTracksAssociatorAtVertexRegionalIter0124', '', '')
process.analyzerCALO.regionalTracksIter4           = cms.untracked.InputTag('displacedAk4JetTracksAssociatorAtVertexRegionalIter4', '', '')

# jet tagging categories
process.analyzerCALO.shortTagThreshold  = cms.untracked.double(shortTagThreshold)
process.analyzerCALO.mediumTagThreshold = cms.untracked.double(mediumTagThreshold)
process.analyzerCALO.longTagThreshold   = cms.untracked.double(longTagThreshold)
process.analyzerCALO.dHTWorkingPoint    = cms.untracked.int32(dHTWorkingPoint)

# jet alpha 
#process.analyzerCALO.jetVertexAssociation = cms.untracked.InputTag("displacedJetVertexAssociation","Var","ANA")

# vertex matched ip info
process.analyzerVTX.secondaryVertexTagInfo    = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPV', '', 'ANA')
process.analyzerVTX.lifetimeIPTagInfo         = cms.untracked.InputTag('displacedLifetimeTagInfos', '', 'ANA')
# calo matched ip and sv info
process.analyzerCALO.secondaryVertexTagInfo   = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPV', '', 'ANA')

# use the tracks matched at the inner hit for the ip tags
process.analyzerCALO.lifetimeIPTagInfo        = cms.untracked.InputTag('displacedLifetimeTagInfos', '', 'ANA') #change to final associaton later
##reconstructed vertex information
process.analyzerCALO.secondaryVertex          = cms.untracked.InputTag('displacedTagsToVertices', '', 'ANA') 
process.analyzerCALO.inclusiveVertexCand      = cms.untracked.InputTag('displacedInclusiveVertexFinder', '', 'ANA')
process.analyzerCALO.inclusiveVertexSecondary = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')
process.analyzerVTX.secondaryVertex           = cms.untracked.InputTag('displacedTagsToVertices', '', 'ANA')
process.analyzerVTX.inclusiveVertexCand       = cms.untracked.InputTag('displacedInclusiveVertexFinder', '', 'ANA')
process.analyzerVTX.inclusiveVertexSecondary  = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')

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
#process.displacedLifetimeTagInfos.jetTracks   = cms.InputTag( "displacedAk4JetTracksAssociatorAtInnerHit" )
process.displacedLifetimeTagInfos.jetTracks   = cms.InputTag( "displacedAk4JetTracksAssociatorAtVertex" )
#process.load('DisplacedJets.Configuration.AdditionalPATSequences_cff')
#process.load('DisplacedJets/DisplacedTriggerFilters/displacedTriggers_cff')

#create the main path to run
process.p = cms.Path()


process.nEventsTotal                      = cms.EDProducer("EventCountProducer")
process.nEventsFiltered                   = cms.EDProducer("EventCountProducer")
process.analyzerCALO.eventCounter         = cms.untracked.InputTag('nEventsTotal')
process.analyzerCALO.eventCounterFiltered = cms.untracked.InputTag('nEventsFiltered')

process.p *= process.nEventsTotal

data_triggers =       [ 'HLT_HT350_DisplacedDijet40_DisplacedTrack_v*',
      'HLT_HT350_DisplacedDijet80_DisplacedTrack_v*',
      'HLT_HT500_DisplacedDijet40_Inclusive_v*',
      'HLT_HT550_DisplacedDijet40_Inclusive_v*',
      'HLT_HT650_DisplacedDijet80_Inclusive_v*',
      'HLT_HT750_DisplacedDijet80_Inclusive_v*']
#      'HLT_VBF_DisplacedJet40_Hadronic_v*',
#      'HLT_VBF_DisplacedJet40_DisplacedTrack_v*',
#      'HLT_VBF_DisplacedJet40_DisplacedTrack_v*',                                                  
#      'HLT_VBF_DisplacedJet40_TightID_DisplacedTrack_v*',                                          
#      'HLT_VBF_DisplacedJet40_TightID_Hadronic_v*',                                                
#      'HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_v*',                                         
#      'HLT_VBF_DisplacedJet40_VTightID_Hadronic_v*',                                               
#      'HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_v*',                                        
#      'HLT_VBF_DisplacedJet40_VVTightID_Hadronic_v*',   
#      'HLT_PFHT800_v*',
#      'HLT_L1_TripleJet_VBF_v*',                                                                  
#      'HLT_PFMET170_v*',                                                                  
 #     'HLT_PFMET170_NoiseCleaned_v*']

mc_triggers = [ 'HLT_HT250_DisplacedDijet40_DisplacedTrack_v*', #only 76x
      'HLT_HT350_DisplacedDijet40_DisplacedTrack_v*',
      'HLT_HT350_DisplacedDijet80_DisplacedTrack_v*', 
      'HLT_HT400_DisplacedDijet40_Inclusive_v*',   #only 76x
      'HLT_HT500_DisplacedDijet40_Inclusive_v*',
      'HLT_HT550_DisplacedDijet40_Inclusive_v*',
      'HLT_HT650_DisplacedDijet80_Inclusive_v*',
      'HLT_HT750_DisplacedDijet80_Inclusive_v*']
      # 'HLT_VBF_DisplacedJet40_Hadronic_v*',
      # 'HLT_VBF_DisplacedJet40_DisplacedTrack_v*',
      # 'HLT_VBF_DisplacedJet40_DisplacedTrack_v*',                                                  
      # 'HLT_VBF_DisplacedJet40_Hadronic_v*',                                                        
      # 'HLT_VBF_DisplacedJet40_TightID_DisplacedTrack_v*',                                          
      # 'HLT_VBF_DisplacedJet40_TightID_Hadronic_v*']
#      'HLT_L1_TripleJet_VBF_v*',                                                                  
#      'HLT_Mu20_v*',                                                                  
#      'HLT_PFMET170_v*',                                                                  
#      'HLT_PFMET170_NoiseCleaned_v*'                                                                  
#      'HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack_v*',                                         
#      'HLT_VBF_DisplacedJet40_VTightID_Hadronic_v*',                                               
#      'HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack_v*',                                        
#      'HLT_VBF_DisplacedJet40_VVTightID_Hadronic_v*',   
#      'HLT_PFHT200_v*',
#      'HLT_PFHT250_v*',
#      'HLT_PFHT300_v*',
#      'HLT_PFHT350_v*',
#      'HLT_PFHT400_v*',
#      'HLT_PFHT800_v*',
#      'HLT_L1_TripleJet_VBF_v*'                                                    

# triggers only in runD and the 76X re-reco              
mc_triggers_76x = ['HLT_HT400_DisplacedDijet40_Inclusive_v*', 'HLT_HT250_DisplacedDijet40_DisplacedTrack_v*', 'HLT_HT325_v*', 'HLT_HT425_v*', 'HLT_HT575_v*', 'HLT_HT275_v*']
data_triggers_76x = ['HLT_HT400_DisplacedDijet40_Inclusive_v*', 'HLT_HT250_DisplacedDijet40_DisplacedTrack_v*','HLT_HT275_v*','HLT_HT325_v*','HLT_HT425_v*','HLT_HT575_v*' ]

# add the extra triggers if they are available into the filter
if is76XTriggers: 
   for trigger in data_triggers_76x: data_triggers.append(trigger)
   for trigger in mc_triggers_76x:   mc_triggers.append(trigger)

if removePFHT800: 
   mc_triggers.remove("HLT_PFHT800_v*")
   data_triggers.remove("HLT_PFHT800_v*")

# only filter on the displaced jet triggers to reduce file size
if dispTriggersOnly:
   mc_triggers = ['HLT_HT500_DisplacedDijet40_Inclusive_v*', 'HLT_HT350_DisplacedDijet40_DisplacedTrack_v*']
   data_triggers = ['HLT_HT500_DisplacedDijet40_Inclusive_v*', 'HLT_HT350_DisplacedDijet40_DisplacedTrack_v*']
   if is76XTriggers: 
      mc_triggers.append('HLT_HT400_DisplacedDijet40_Inclusive_v*') 
      mc_triggers.append('HLT_HT250_DisplacedDijet40_DisplacedTrack_v*')
      data_triggers.append('HLT_HT400_DisplacedDijet40_Inclusive_v*') 
      data_triggers.append('HLT_HT250_DisplacedDijet40_DisplacedTrack_v*')

# trigger bits to keep
process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring(*data_triggers ),
    hltResults = cms.InputTag( "TriggerResults", "", trigger_process ),
    l1tResults = cms.InputTag( "" ),
    l1tIgnoreMask = cms.bool( False ),
    l1techIgnorePrescales = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True )
)


process.mcTriggerSelection = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring(*mc_triggers),
    hltResults = cms.InputTag( "TriggerResults", "", trigger_process ),
    l1tResults = cms.InputTag( "" ),
    l1tIgnoreMask = cms.bool( False ),
    l1techIgnorePrescales = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True )
)
# trigger bits to keep
process.singleMuTrigger = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring(
      'HLT_IsoMu20_v*'                                                                 
      ),
    hltResults = cms.InputTag( "TriggerResults", "", trigger_process ),
    l1tResults = cms.InputTag( "" ),
    l1tIgnoreMask = cms.bool( False ),
    l1techIgnorePrescales = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True )
)

if doApplyTrigger: #apply the triggers and run dj tagging
   if doApplySingleMu: 
      process.p *= process.singleMuTrigger

   if not (isOnlyMu and doApplySingleMu): # if we want more than just the single mu20 trigger
      if isMC:
         process.p *= process.mcTriggerSelection 
      else:
         process.p *= process.triggerSelection 

#count the events filtered 
process.p *= process.nEventsFiltered

#always add the jet corrections and tagging
process.p *=  process.correctJets * process.djtagging

#add in the regional track jet associations if we are running on RAWAOD+
if addRegionalTracking:  process.p *= process.regionalTrackAssocations 

if doedm: #just dump the edm output of the djtagging sequence no analyzer
    process.btag_output = cms.EndPath(process.test_output)
else: # run the analysis and tree dumper using the djtagging sequence output
    process.p *= process.analyzerCALO

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
#from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
#process = customisePostLS1(process)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool( True )
)
