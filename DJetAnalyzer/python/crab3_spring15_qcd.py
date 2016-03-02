from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName     = 'QCD_SPRING15_JAN7'
config.General.workArea        = 'CRAB3_QCD_SPRING15_JAN7'
config.General.transferOutputs = True
config.General.transferLogs    = False

config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'dev.py'
config.JobType.outputFiles = ['djana.root']

#config.Data.inputDataset = '/MinimumBias/Commissioning2015-v1/RAW'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 4
config.Data.totalUnits = -1 #config.Data.unitsPerJob * NJOBS
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.runRange = '193093-193999' # '193093-194075'

config.Data.outLFNDirBase  = '/store/group/phys_susy/razor/josh/RAZOR_DIJET/DJANALYSIS/QCD_SPRING15_JAN7_TRIGGERAPPLIED'
config.Data.publication = False
config.Data.ignoreLocality   = True #seth zenz fix 
#config.Data.publishDataName = 'CRAB3_tutorial_May2015_Data_analysis'
config.Site.storageSite = 'T2_CH_CERN'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    #config.General.workArea = 'crab_projects'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################


    # config.General.requestName = 'qcd5_10'
    # config.Data.inputDataset = '/QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM'
    # # config.Data.unitsPerJob = 2
    # # =config.Data.totalUnits = 4
    # submit(config)

    # config.General.requestName = 'qcd10_15'
    # config.Data.inputDataset = '/QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM'
    # # config.Data.unitsPerJob = 2
    # # =config.Data.totalUnits = 4
    # submit(config)

    # config.General.requestName = 'qcd15_30'
    # config.Data.inputDataset = '/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM'
    # # config.Data.unitsPerJob = 2
    # # =config.Data.totalUnits = 4
    # submit(config)

    # config.General.requestName = 'qcd30_50'
    # config.Data.inputDataset = '/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM'
    # # config.Data.unitsPerJob = 2
    # #config.Data.totalUnits = 4
    # submit(config)

    config.General.requestName = 'qcd50_80'
    config.Data.inputDataset = '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd80_120'
    config.Data.inputDataset = '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd120_170'
    config.Data.inputDataset = '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd170_300'
    config.Data.inputDataset = '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd300_470'
    config.Data.inputDataset = '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd470_600'
    config.Data.inputDataset = '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd600_800'
    config.Data.inputDataset = '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd800_1000'
    config.Data.inputDataset = '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd1000_1400'
    config.Data.inputDataset = '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd1400_1800'
    config.Data.inputDataset = '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd1800_2400'
    config.Data.inputDataset = '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd2400_3200'
    config.Data.inputDataset = '/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

    config.General.requestName = 'qcd3200'
    config.Data.inputDataset = '/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # config.Data.unitsPerJob = 2
    #config.Data.totalUnits = 4
    submit(config)

# /QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW

# /QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW

# /QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7_ext1-v1/GEN-SIM-RAW

# /QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW
# /QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW


# /QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW


# /QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v2/GEN-SIM-RAW
# /QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW
