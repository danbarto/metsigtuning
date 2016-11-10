import os
from CRABClient.UserUtilities import getUsernameFromSiteDB
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs = True
config.General.requestName = 'MetSig_Data_2016G_PromptReco'
config.General.workArea = config.General.requestName

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'makentuple_cfg.py'
config.JobType.outputFiles = ['ntuple.root']
config.JobType.inputFiles = ['Spring16_25nsV6_DATA.db','Spring16_25nsV6_MC.db']

config.section_("Data")
config.Data.lumiMask = 'json/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON.txt'

#config.Data.inputDataset = '/DoubleMuon/Run2016B-PromptReco-v2/MINIAOD'
#config.Data.runRange = '272007-275376'

#config.Data.inputDataset = '/DoubleMuon/Run2016C-PromptReco-v2/MINIAOD'
#config.Data.runRange = '275657-276283'

#config.Data.inputDataset = '/DoubleMuon/Run2016D-PromptReco-v2/MINIAOD'
#config.Data.runRange = '276315-276811'

#config.Data.inputDataset = '/DoubleMuon/Run2016E-PromptReco-v1/MINIAOD'

config.Data.inputDataset = '/DoubleMuon/Run2016G-PromptReco-v1/MINIAOD'
config.Data.runRange = '278820-280385'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.ignoreLocality = False

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.unitsPerJob = 50
#config.Data.totalUnits = 100

config.section_("Site")
#config.Site.blacklist = ['T2_US_Purdue', 'T2_US_Nebraska', 'T2_US_MIT', 'T2_US_Caltech']
config.Site.storageSite = 'T2_AT_Vienna'

config.section_("User")

config.section_("Debug")

