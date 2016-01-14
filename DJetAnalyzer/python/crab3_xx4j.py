from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'XX4J_CENTRAL_DEC28_v2'
config.General.workArea = 'CRAB3_XX4J_CENTRAL_DEC28_v2'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'dev.py'
config.JobType.outputFiles = ['signal.root']

#config.Data.inputDataset = '/MinimumBias/Commissioning2015-v1/RAW'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1 #config.Data.unitsPerJob * NJOBS

config.Data.outLFNDirBase = '/store/group/phys_susy/razor/josh/RAZOR_DIJET/DJANALYSIS/XX4J_CENTRAL_DEC28'
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_May2015_Data_analysis'
config.Site.storageSite = 'T2_CH_CERN'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
     ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
     #############################################################################################

    # config.General.requestName  = 'xx4j_m1000_ctau1000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1000_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1000_ctau10mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1000_CTau-10mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1000_ctau2000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1000_ctau300mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1000_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1000_ctau30mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1000_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1000_ctau3mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1000_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1000_ctau1mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1000_CTau-1mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 


    config.General.requestName = 'xx4j_m100_ctau1000mm'
    config.Data.inputDataset   = '/XXTo4J_M-100_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    submit(config)	 
    config.General.requestName = 'xx4j_m100_ctau100mm'
    config.Data.inputDataset   = '/XXTo4J_M-100_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    submit(config)	 
    config.General.requestName = 'xx4j_m100_ctau10mm'
    config.Data.inputDataset   = '/XXTo4J_M-100_CTau-10mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    submit(config)	 
    config.General.requestName = 'xx4j_m100_ctau1mm'
    config.Data.inputDataset   = '/XXTo4J_M-100_CTau-1mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    submit(config)	 
    config.General.requestName = 'xx4j_m100_ctau2000mm'
    config.Data.inputDataset   = '/XXTo4J_M-100_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    submit(config)	 
    config.General.requestName = 'xx4j_m100_ctau300mm'
    config.Data.inputDataset   = '/XXTo4J_M-100_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    submit(config)	 
    config.General.requestName = 'xx4j_m100_ctau30mm'
    config.Data.inputDataset   = '/XXTo4J_M-100_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    submit(config)	 
    config.General.requestName = 'xx4j_m100_ctau3mm'
    config.Data.inputDataset   = '/XXTo4J_M-100_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    submit(config)	 
    
    # config.General.requestName = 'xx4j_m1500_ctau1000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1500_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1500_ctau100mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1500_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1500_ctau10mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1500_CTau-10mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1500_ctau3mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1500_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1500_ctau1mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1500_CTau-1mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1500_ctau2000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1500_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1500_ctau300mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1500_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1500_ctau30mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1500_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m1500_ctau3mm'
    # config.Data.inputDataset   = '/XXTo4J_M-1500_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 

    # config.General.requestName = 'xx4j_m3000_ctau1000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-3000_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m3000_ctau10mm'
    # config.Data.inputDataset   = '/XXTo4J_M-3000_CTau-10mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m3000_ctau1mm'
    # config.Data.inputDataset   = '/XXTo4J_M-3000_CTau-1mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m3000_ctau2000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-3000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m3000_ctau300mm'
    # config.Data.inputDataset   = '/XXTo4J_M-3000_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m3000_ctau30mm'
    # config.Data.inputDataset   = '/XXTo4J_M-3000_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m3000_ctau3mm'
    # config.Data.inputDataset   = '/XXTo4J_M-3000_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 

    # config.General.requestName = 'xx4j_m300_ctau1000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-300_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m300_ctau100mm'
    # config.Data.inputDataset   = '/XXTo4J_M-300_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m300_ctau10mm'
    # config.Data.inputDataset   = '/XXTo4J_M-300_CTau-10mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m300_ctau1mm'
    # config.Data.inputDataset   = '/XXTo4J_M-300_CTau-1mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m300_ctau2000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-300_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m300_ctau300mm'
    # config.Data.inputDataset   = '/XXTo4J_M-300_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m300_ctau30mm'
    # config.Data.inputDataset   = '/XXTo4J_M-300_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m300_ctau3mm'
    # config.Data.inputDataset   = '/XXTo4J_M-300_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 

    # config.General.requestName = 'xx4j_m500_ctau1000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-500_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m500_ctau100mm'
    # config.Data.inputDataset   = '/XXTo4J_M-500_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m500_ctau10mm'
    # config.Data.inputDataset   = '/XXTo4J_M-500_CTau-10mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m500_ctau1mm'
    # config.Data.inputDataset   = '/XXTo4J_M-500_CTau-1mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m500_ctau2000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-500_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m500_ctau300mm'
    # config.Data.inputDataset   = '/XXTo4J_M-500_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m500_ctau30mm'
    # config.Data.inputDataset   = '/XXTo4J_M-500_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m500_ctau3mm'
    # config.Data.inputDataset   = '/XXTo4J_M-500_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m500_ctau1mm'
    # config.Data.inputDataset   = '/XXTo4J_M-500_CTau-1mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 


    # config.General.requestName = 'xx4j_m50_ctau2000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-50_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m50_ctau1000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-50_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m50_ctau300mm'
    # config.Data.inputDataset   = '/XXTo4J_M-50_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m50_ctau100mm'
    # config.Data.inputDataset   = '/XXTo4J_M-50_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m50_ctau30mm'
    # config.Data.inputDataset   = '/XXTo4J_M-50_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m50_ctau10mm'
    # config.Data.inputDataset   = '/XXTo4J_M-50_CTau-10mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m50_ctau1mm'
    # config.Data.inputDataset   = '/XXTo4J_M-50_CTau-1mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m50_ctau2000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-50_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m50_ctau3mm'
    # config.Data.inputDataset   = '/XXTo4J_M-50_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 

    # config.General.requestName = 'xx4j_m700_ctau1000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-700_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m700_ctau100mm'
    # config.Data.inputDataset   = '/XXTo4J_M-700_CTau-100mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m700_ctau10mm'
    # config.Data.inputDataset   = '/XXTo4J_M-700_CTau-10mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m700_ctau1mm'
    # config.Data.inputDataset   = '/XXTo4J_M-700_CTau-1mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m700_ctau2000mm'
    # config.Data.inputDataset   = '/XXTo4J_M-700_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m700_ctau300mm'
    # config.Data.inputDataset   = '/XXTo4J_M-700_CTau-300mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m700_ctau30mm'
    # config.Data.inputDataset   = '/XXTo4J_M-700_CTau-30mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
    # config.General.requestName = 'xx4j_m700_ctau3mm'
    # config.Data.inputDataset   = '/XXTo4J_M-700_CTau-3mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    # submit(config)	 
