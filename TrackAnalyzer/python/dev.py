import FWCore.ParameterSet.Config as cms

# output options (to be appended to the file name outputted)
appendLifetime = "single30"
appendBkg      = "120_170"

# flags for running
nevents       = 2000
isSignalMC    = True
doGenMatch    = True
doSimVtxMatch = True
isMC          = True
doedm         = False

# analysis cuts
cut_jetPt = 80
cut_jetEta = 2.0


input_file = None
if isSignalMC:
#    input_file = 'file:/afs/cern.ch/work/h/hardenbr/2015/DIJET/GEN_SIGNAL_TEST/dijet_700_300_ctau30.root'
#    input_file = 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU1000_40BX25_AOD/dijet_700_300_ctau1000_94_1_3mG.root'
    input_file = 'file:/afs/cern.ch/user/h/hardenbr/2014/LL_DIJET/SIGNAL_GENERATION/PAIR_XX/CMSSW_7_2_0/src/Configuration/GenProduction/python/ThirteenTeV/XXTo4J_M-1500_CTau-30mm_step3.root'
#   input_file = 'file:/afs/cern.ch/work/h/hardenbr/2015/DIJET/GEN_SIGNAL_TEST/dijet_700_300_ctau300.root'
#   input_file = 'file:/afs/cern.ch/work/h/hardenbr/2015/DIJET/GEN_SIGNAL_TEST/dijet_700_300_ctau3000.root'
#   input_file = 'file:/afs/cern.ch/work/h/hardenbr/2015/DIJET/GEN_SIGNAL_TEST/dijet_700_300_ctau3.root'
#   input_file = 'file:/afs/cern.ch/user/t/tkolberg/public/hepmcreco_RAW2DIGI_RECO.root'
else:
    input_file = 'file:/afs/cern.ch/work/h/hardenbr/TEST_FILES/QCD_Pt-120to170_Tune4C_13TeV_pythia8_castor_tsg_PU40bx25_POSTLS162_V2-v1.root'
#    "file:/afs/cern.ch/work/h/hardenbr/QCD_Pt-50to80_Tune4C_13TeV_pythia8_AOD.root'
#    input_file =  'file:/afs/cern.ch/work/h/hardenbr/QCD_470_600_AOD_40bx25.root'

process = cms.Process("ANA")
proc_label = "DIGI2RAW"

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

#geometry and global tag sequences
#process.load('Configuration.StandardSequences.Services_cff')
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'MCRUN2_72_V1A::All'

process.GlobalTag = cms.ESSource( "PoolDBESSource",
    globaltag = cms.string( "auto:run2_mc" ),
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
#    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'MCRUN2_72_V4A')
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'auto:run2_mc',conditions='TrackerAlignmentExtendedError_2011Realistic_v1_mc,TrackerAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+MuonDTAPEObjectsExtended_v0_mc,DTAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+MuonCSCAPEObjectsExtended_v0_mc,CSCAlignmentErrorExtendedRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalSamplesCorrelation_mc,EcalSamplesCorrelationRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalPulseShapes_mc,EcalPulseShapesRcd,frontier://FrontierProd/CMS_CONDITIONS+EcalPulseCovariances_mc,EcalPulseCovariancesRcd,frontier://FrontierProd/CMS_CONDITIONS')
#    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'MCRUN2_72_V1A')
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_CONDITIONS'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    for pset in process.GlobalTag.toGet.value():
        pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
    # fix for multi-run processing
    process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
    process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nevents))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_100_1_Eai.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_101_1_wzw.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_102_1_iv7.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_103_1_Oe5.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_104_1_gEW.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_105_1_Yzo.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_106_1_HWy.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_107_1_4KN.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_108_1_8Uv.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_109_1_iar.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_10_1_YLA.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_110_1_ZE2.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_111_1_YZK.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_112_1_Nsn.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_113_1_hpe.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_114_1_fgo.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_115_1_AyY.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_116_1_9H8.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_117_1_UCw.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_118_1_uTP.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_119_1_7rV.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_11_1_BDX.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_120_1_d9M.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_12_1_LnF.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_13_1_Qxo.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_14_1_IV4.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_15_1_BSH.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_16_1_5u4.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_17_1_44N.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_18_1_t15.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_19_1_LaR.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_1_1_gKa.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_20_1_zf6.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_21_1_PAJ.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_22_1_YGk.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_23_1_9mz.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_24_1_gFF.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_25_1_Gqa.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_26_1_YGd.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_27_1_VvU.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_28_1_0Uf.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_29_1_f54.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_2_1_ywf.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_30_1_eqv.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_31_1_05F.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_32_1_Z1x.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_33_1_NfH.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_34_1_GfL.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_35_1_IP7.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_36_1_BhU.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_37_1_KF2.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_38_1_Dd7.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_39_1_LkR.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_3_1_9T0.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_40_1_mHy.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_41_1_Wms.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_42_1_S6V.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_43_1_bkd.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_44_1_94t.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_45_1_MRq.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_46_1_CJq.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_47_1_qw1.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_48_1_tdT.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_49_1_GNO.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_4_1_kpz.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_50_1_9cf.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_51_1_dvA.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_52_1_Kb1.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_53_1_xrR.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_54_1_g9O.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_55_1_NaM.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_56_1_UD0.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_57_1_EMY.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_58_1_LHW.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_59_1_Fqf.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_5_1_f6t.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_60_1_yip.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_61_1_YB7.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_62_1_RtL.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_63_1_sbD.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_64_1_eLT.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_65_1_rvK.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_66_1_u4E.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_67_1_ovL.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_68_1_TmD.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_69_1_Lo2.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_6_1_Jsr.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_70_1_uNj.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_71_1_cYF.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_72_1_QJH.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_73_1_3Ck.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_74_1_XKW.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_75_1_2WB.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_76_1_NnF.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_77_1_rw8.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_78_1_Usn.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_79_1_Jqf.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_7_1_xkG.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_80_1_xIU.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_81_1_ddK.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_82_1_o6l.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_83_1_WqX.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_84_1_EM7.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_85_1_wPm.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_86_1_7cd.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_87_1_JvT.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_88_1_Hjl.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_89_1_4hO.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_8_1_Nhi.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_90_1_Cvu.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_91_1_PuN.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_92_1_Az3.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_93_1_rYM.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_94_1_qcb.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_95_1_FF8.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_96_1_iOj.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_97_1_emZ.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_98_1_aIn.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_99_1_OCm.root',
# 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/group/phys_susy/razor/josh/RAZOR_DIJET/DIJET_MH700_MX300_CTAU300_40BX25_AOD/dijet_700_300_ctau300_9_1_u6w.root'
        # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_0.root',
        # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_1.root',
        # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_2.root',
        # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_3.root',
        # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_4.root',
        # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_5.root',
        # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_6.root',
        # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_7.root',
        # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_8.root',
        # 'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/user/mwalker/SIGNAL/Displaced/LL_sbottom_500_100.0mm/step2_LL_sbottom_500_100.0mm_9.root'
#        'file:/afs/cern.ch/work/h/hardenbr/2015/DIJET/SIGNAL_SAMPLES/TEST_SAMPLES/DISPLACED_SUSSY/EXO-Phys14DR-00001_step2.root'
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_100_1_3ZU_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_101_1_6TZ_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_102_1_gm6_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_103_1_Wo4_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_104_1_FqX_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_105_1_d3J_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_106_1_X70_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_107_1_ufI_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_108_1_w0s_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_109_1_VAv_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_10_1_5Kw_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_110_1_Qlg_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_111_1_Nry_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_112_1_umx_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_113_1_5BV_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_114_1_a8z_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_115_1_pN5_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_116_1_Rjy_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_117_1_Xyc_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_118_1_eEC_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_119_1_7aP_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_11_1_jua_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_120_1_ENk_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_12_1_iN5_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_13_1_RY4_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_14_1_xXA_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_15_1_gDg_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_16_1_aH7_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_17_1_W14_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_18_1_ysT_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_1_1_It9_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_20_1_ber_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_21_1_mSO_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_22_1_vVX_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_23_1_H6o_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_24_1_nSJ_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_26_1_76D_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_27_1_pEF_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_28_1_EWK_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_29_1_e5i_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_2_1_3TK_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_30_1_P2A_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_32_1_5f9_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_33_1_5TZ_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_34_1_xN1_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_35_1_xMn_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_36_1_CQV_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_37_1_70q_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_38_1_iwv_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_39_1_wzj_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_3_1_JOv_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_40_1_Wnx_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_41_1_ioC_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_42_1_kXz_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_43_1_K2c_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_44_1_r75_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_45_1_Ob9_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_46_1_mrG_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_47_1_DTU_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_48_1_yoc_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_4_1_KgE_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_50_1_HF7_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_51_1_OVk_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_53_1_1KW_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_55_1_m4L_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_57_1_e5y_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_5_1_Bf4_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_61_1_UjV_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_62_1_WlE_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_63_1_oc9_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_67_1_vEr_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_68_1_FQc_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_6_1_fCq_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_70_1_hvv_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_71_1_n7f_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_72_1_Skw_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_73_1_CVZ_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_74_1_67S_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_75_1_ZEa_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_76_1_TyT_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_78_1_NMd_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_79_1_YtE_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_7_1_ylH_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_80_1_t0H_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_81_1_REX_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_82_1_C8V_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_83_1_VNy_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_85_1_wQl_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_86_1_sxt_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_88_1_HWa_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_8_1_GYa_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_90_1_E5x_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_92_1_HzY_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_93_1_a1c_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_96_1_aRb_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_97_1_RZJ_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_98_1_zmW_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_99_1_UcV_step3.root',
'file:/afs/cern.ch/user/h/hardenbr/eos/cms/store/caf/user/hardenbr/DIJET/FOR_EXO/XXTo4J_MX1500_1000mm/XXTo4J_M-1500_CTau-1000mm_step2_9_1_B59_step3.root'
      #        'file:/afs/cern.ch/work/h/hardenbr/QCD_Pt-50to80_Tune4C_13TeV_pythia8_AOD.root'
      #        'file:/afs/cern.ch/work/h/hardenbr/HTo2LongLivedTo4L_MH_700_MFF_300_CTau30_TSG_PU40BX25_AODSIM.root'
      #       'file:/afs/cern.ch/work/h/hardenbr/HTo2LongLivedTo4L_MH_700_MFF_300_CTau30_TSG_PU40BX25_AODSIM_v6.root'
      #        'file:/afs/cern.ch/work/h/hardenbr/QCD_470_600_AOD_40bx25.root'
      #        'file:/afs/cern.ch/work/h/hardenbr/TEST_FILES/HTo2LongLivedTo4L_MH_700_MFF_300_CTau30_TSG_PU40BX25_AODSIM_10ev.root'
      )
                            )



################################################################################################

#configure the analyzers
process.analyzerVTX = cms.EDAnalyzer('TrackAnalyzer')
process.analyzerCALO = cms.EDAnalyzer('TrackAnalyzer')

process.analyzerVTX.debugLevel  = cms.untracked.int32(0)
process.analyzerCALO.debugLevel = cms.untracked.int32(0)

#output configuration
if isSignalMC:
    process.analyzerVTX.outputFileName = cms.untracked.string('signalVTX%s.root' % appendLifetime)
    process.analyzerCALO.outputFileName = cms.untracked.string('signalCALO%s.root' % appendLifetime)

else:
    process.analyzerVTX.outputFileName = cms.untracked.string('qcdVTX%s.root' % appendBkg)
    process.analyzerCALO.outputFileName = cms.untracked.string('qcdCALO%s.root' % appendBkg)


process.analyzerCALO.jetTreeName = cms.untracked.string('jets')
process.analyzerCALO.trackTreeName = cms.untracked.string('tracks')
process.analyzerCALO.vertexTreeName = cms.untracked.string('vtx')
process.analyzerCALO.genTreeName = cms.untracked.string('gen')

process.analyzerVTX.jetTreeName = cms.untracked.string('jets')
process.analyzerVTX.trackTreeName = cms.untracked.string('tracks')
process.analyzerVTX.vertexTreeName = cms.untracked.string('vtx')
process.analyzerVTX.genTreeName = cms.untracked.string('gen')

process.analyzerVTX.isMC  = cms.untracked.bool(isMC)
process.analyzerCALO.isMC = cms.untracked.bool(isMC)

process.analyzerVTX.doGenMatch  = cms.untracked.bool(doGenMatch)
process.analyzerCALO.doGenMatch = cms.untracked.bool(doGenMatch)

process.analyzerVTX.doSimMatch  = cms.untracked.bool(doSimVtxMatch)
process.analyzerCALO.doSimMatch  = cms.untracked.bool(doSimVtxMatch)

process.analyzerVTX.isSignalMC  = cms.untracked.bool(isSignalMC)
process.analyzerCALO.isSignalMC = cms.untracked.bool(isSignalMC)

#tags
process.analyzerVTX.generalTracks  = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerVTX.ak5CaloJets    = cms.untracked.InputTag('ak4CaloJets', '', '')
process.analyzerVTX.genParticles   = cms.untracked.InputTag('genParticles', '', '')
process.analyzerCALO.generalTracks = cms.untracked.InputTag('generalTracks', '', '')
process.analyzerCALO.ak5CaloJets   = cms.untracked.InputTag('ak4CaloJets', '', '')
process.analyzerCALO.genParticles  = cms.untracked.InputTag('genParticles', '', '')

# vertex matched ip info
process.analyzerVTX.secondaryVertexTagInfo   = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPV', '', 'ANA')
process.analyzerVTX.lifetimeIPTagInfo        = cms.untracked.InputTag('displacedLifetimeTagInfos', '', 'ANA')

# calo matched ip and sv info
process.analyzerCALO.secondaryVertexTagInfo   = cms.untracked.InputTag('displacedSecondaryVertexTagInfosNoPVCaloFace', '', 'ANA')
process.analyzerCALO.lifetimeIPTagInfo        = cms.untracked.InputTag('displacedLifetimeTagInfosCaloFace', '', 'ANA')

##reconstructed vertex information
process.analyzerCALO.secondaryVertex          = cms.untracked.InputTag('displacedTagsToVerticesCaloFace', '', 'ANA')
process.analyzerCALO.inclusiveVertexCand      = cms.untracked.InputTag('displacedInclusiveVertexFinderJetMatchedTracksCaloFace', '', 'ANA')
process.analyzerCALO.inclusiveVertexSecondary = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')

process.analyzerVTX.secondaryVertex          = cms.untracked.InputTag('displacedTagsToVertices', '', 'ANA')
process.analyzerVTX.inclusiveVertexCand      = cms.untracked.InputTag('displacedInclusiveVertexFinderJetMatchedTracks', '', 'ANA')
process.analyzerVTX.inclusiveVertexSecondary = cms.untracked.InputTag('displacedInclusiveSecondaryVertices', '', 'ANA')

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
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
                                      outputCommands = cms.untracked.vstring(
        'keep *_*_*_*'))#        'keep *_*_*_ANA'))

#run the displaced jet tags
process.load('DisplacedJets.Configuration.RecoDJTag_cff')

#config gen matching for the output displaced vertices 
#if isSignalMC:
#    process.displacedTagsToVertices.isSignalMC = cms.untracked.bool(True)
#    process.displacedTagsToVertices.doGenMatch = cms.untracked.bool(False)

if doedm:
    process.p = cms.Path(process.djtagging)    
    process.btag_output = cms.EndPath( process.test_output)
else:
    process.p = cms.Path(process.djtagging + process.analyzerVTX + process.analyzerCALO)


# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs                                                                  
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs                                                           
process = customisePostLS1(process)
