from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'LNuAJJ_EWK_test94X'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles = ['Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFchs.txt','Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFPuppi.txt','Summer16_07Aug2017_V11_MC_L2Relative_AK4PFPuppi.txt','Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFPuppi.txt']
config.JobType.psetName    = 'analysis_mc.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
###config.Data.outputPrimaryDataset = 'test_VBS_WGAMMA_94X'
config.Data.inputDataset = '/LNuAJJ_EWK_MJJ-120_13TeV-madgraph-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.outputDatasetTag = 'LNuAJJ_EWK_test94X'

config.section_("Site")
config.Site.storageSite = 'T2_CN_Beijing'
