import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys
from os import environ
import os

options = VarParsing ('analysis')

options.setDefault( 'outputFile',
      'ntuple.root'
      )

options.register( 'globalTag',
      #'80X_mcRun2_asymptotic_2016_TrancheIV_v8',
      #'80X_dataRun2_2016SeptRepro_v7',

      #'80X_mcRun2_asymptotic_2016_miniAODv2_v1',
      #'80X_mcRun2_asymptotic_v20',
      #'80X_dataRun2_Prompt_ICHEP16JEC_v0',
      #'80X_mcRun2_asymptotic_2016_miniAODv2_v1',
      #'80X_mcRun2_asymptotic_v17',
      #'80X_dataRun2_v18',
      #'80X_dataRun2_Prompt_v11',
      '80X_dataRun2_Candidate_2016_09_02_10_26_48',
      #'74X_dataRun2_Prompt_v1',
      #'MCRUN2_74_V9',
      #'GR_P_V56',
      VarParsing.multiplicity.singleton,
      VarParsing.varType.string,
      "CMS Global Tag"
      )

#MC switch
options.register( 'runOnMC',
      True,
      VarParsing.multiplicity.singleton,
      VarParsing.varType.bool,
      "mc or data"
      )

options.register( 'mfTag',
      'RECO',
      VarParsing.multiplicity.singleton,
      VarParsing.varType.string,
      "for MET filters"
      )

options.parseArguments()

process = cms.Process("test")

#process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("metsigtuning.MakeNtuple.makentuple_cfi")
#process.load("RecoMET/METProducers.METSignificanceObjects_cfi")

process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')#
process.load('Configuration.StandardSequences.MagneticField_cff')#
#process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.GlobalTag.globaltag = ( options.globalTag )

JECdata = 'Summer16_23Sep2016AllV4_DATA'
JECMC   = 'Summer16_23Sep2016V4_MC'
if options.runOnMC:
    JECdb = JECMC
else: JECdb = JECdata


process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_'+JECdb+'_AK4PFchs'),
            # tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs')
            ),
      ),
      connect = cms.string('sqlite:'+JECdb+'.db')
      #connect = cms.string('sqlite:Fall15_V2_DATA.db')
)

## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #'/store/data/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/00039A2E-D7A7-E511-98EE-3417EBE64696.root'
       #'/store/data/Run2016G/DoubleMuon/MINIAOD/23Sep2016-v1/100000/00993A51-DF90-E611-A4EE-7845C4FC3650.root'
       #'/store/data/Run2016G/DoubleMuon/MINIAOD/23Sep2016-v1/100000/084F88CC-548F-E611-BEED-549F35AD8B7B.root'
       #'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/00000/0251DBB7-201B-E611-8653-0CC47A4F1C2E.root'
       #'/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/02A85EE9-70BA-E511-A0A2-0CC47A4D7678.root'
       #'file:00039A2E-D7A7-E511-98EE-3417EBE64696.root'
       #'/store/relval/CMSSW_8_0_20/RelValTTbar_13/MINIAODSIM/PU25ns_80X_mcRun2_asymptotic_2016_TrancheIV_v4_Tr4GT_v4-v1/00000/A8C282AE-D37A-E611-8603-0CC47A4C8ECE.root'
       #'/store/data/Run2016G/DoubleMuon/MINIAOD/PromptReco-v1/000/278/819/00000/E40327C9-9F63-E611-80EE-FA163E3490D9.root'
       #'/store/data/Run2016G/DoubleMuon/MINIAOD/PromptReco-v1/000/278/820/00000/227B551D-AD64-E611-A12B-FA163E951746.root'
       #'/store/data/Run2016C/DoubleMuon/MINIAOD/PromptReco-v2/000/275/601/00000/4423F253-7B3A-E611-8707-02163E013706.root'
       #'/store/data/Run2016C/DoubleMuon/MINIAOD/PromptReco-v2/000/275/657/00000/3460EDF8-7F3B-E611-9318-02163E01461C.root'
       '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/0EA60289-18C4-E611-8A8F-008CFA110AB4.root'
       #'/store/mc/RunIISummer16MiniAODv2/WW_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/0449B17C-BAD7-E611-8430-0025905B85EC.root'
       #'/store/data/Run2016G/DoubleMuon/MINIAOD/23Sep2016-v1/100000/00DD00F8-008C-E611-8CD0-00266CFFC9C4.root'
       #'/store/mc/RunIISpring16MiniAODv1/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/00000/06B334D9-4DFE-E511-B1A1-001E67A40604.root'
       #'/store/data/Run2016G/DoubleMuon/MINIAOD/03Feb2017-v1/100000/04996D94-19EB-E611-8EC6-001E67E0061C.root' 
    )
)

# print statistics
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)


# update JECs
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

corrections = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
if not options.runOnMC: # Do not forget 'L2L3Residual' on data!
   corrections += ['L2L3Residual']

updateJetCollection(
      process,
      jetSource = cms.InputTag('slimmedJets'),
      labelName = 'UpdatedJEC',
      jetCorrections = ('AK4PFchs', corrections, 'None') 
      )

#update MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
      isData =  not options.runOnMC,
      )

process.test = cms.EDAnalyzer('MakeNtuple',
      src                       = cms.InputTag("packedPFCandidates"),
      #jets                      = cms.InputTag("slimmedJets"), #same as selectedUpdatedPatJetsUpdatedJEC
      #jets                      = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
      jets                      = cms.InputTag("cleanedPatJets"),
      leptons                   = cms.VInputTag("slimmedElectrons", "slimmedMuons", "slimmedPhotons"),
      met                       = cms.InputTag("slimmedMETs"),
      metPF                     = cms.InputTag("patPFMet"),
      metPFT1                   = cms.InputTag("patPFMetT1"),
      metPFT1Smear              = cms.InputTag("patPFMetT1Smear"),
      metPFT1SmearJetResUp      = cms.InputTag("patPFMetT1SmearJetResUp"),
      metPFT1SmearJetResDown    = cms.InputTag("patPFMetT1SmearJetResDown"),
      metPFT1JetResUp           = cms.InputTag("patPFMetT1JetResUp"),
      metPFT1JetResDown         = cms.InputTag("patPFMetT1JetResDown"),
      metPFT1JetEnUp            = cms.InputTag("patPFMetT1JetEnUp"),
      metPFT1JetEnDown          = cms.InputTag("patPFMetT1JetEnDown"),
      metPFT1UnclusteredEnUp    = cms.InputTag("patPFMetT1UnclusteredEnUp"),
      metPFT1UnclusteredEnDown  = cms.InputTag("patPFMetT1UnclusteredEnDown"),
      muons                     = cms.InputTag("slimmedMuons"),
      electrons                 = cms.InputTag("slimmedElectrons"),
      vertices                  = cms.InputTag("offlineSlimmedPrimaryVertices"),
      runOnMC                   = cms.untracked.bool(options.runOnMC),
      addPileupInfo             = cms.InputTag("slimmedAddPileupInfo"),
      rho                       = cms.InputTag("fixedGridRhoAll"),
      generator                 = cms.InputTag("generator"),
      lheprod                   = cms.InputTag("externalLHEProducer"),
      pileup                    = cms.untracked.InputTag("slimmedAddPileupInfo"),
      srcJetSF                  = cms.string('AK4PFchs'),
      srcJetResPt               = cms.string('AK4PFchs_pt'),
      srcJetResPhi              = cms.string('AK4PFchs_phi'),
)

process.load('Configuration.StandardSequences.Services_cff')
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *

JERdataMC = 'MC'
if not options.runOnMC:
   JERdataMC = 'DATA'

process.jer = cms.ESSource("PoolDBESSource",
      CondDBSetup,
      toGet = cms.VPSet(
         # Pt Resolution
         cms.PSet(
            record = cms.string('JetResolutionRcd'),
            tag    = cms.string('JR_Spring16_25nsV6_'+JERdataMC+'_PtResolution_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs_pt')
            ),
         # Phi Resolution
         cms.PSet(
            record = cms.string('JetResolutionRcd'),
            tag    = cms.string('JR_Spring16_25nsV6_'+JERdataMC+'_PhiResolution_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs_phi')
            ),
         # Scale factors
         cms.PSet(
            record = cms.string('JetResolutionScaleFactorRcd'),
            tag    = cms.string('JR_Spring16_25nsV6_'+JERdataMC+'_SF_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs')
            ),
         ),
      connect = cms.string('sqlite:Spring16_25nsV6_'+JERdataMC+'.db')
      )

process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

# trigger filter                
#trigger_paths = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v']
trigger_paths = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v', 'HLT_IsoMu24_v', 'HLT_IsoTkMu24_v']
trigger_pattern = [path+"*" for path in trigger_paths]
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.triggerSelection = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::HLT",
      HLTPaths = trigger_pattern,
      throw=False
      )

#
## MET filters
#
process.HBHENoiseFilter = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::"+options.mfTag,
      HLTPaths = cms.vstring('Flag_HBHENoiseFilter'),
      throw=False
      )
process.HBHENoiseIsoFilter = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::"+options.mfTag,
      HLTPaths = cms.vstring('Flag_HBHENoiseIsoFilter'),
      throw=False
      )
process.globalTightHalo2016Filter = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::"+options.mfTag,
      HLTPaths = cms.vstring('Flag_globalTightHalo2016Filter'),
      throw=False
      )
process.EcalDeadCellTriggerPrimitiveFilter = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::"+options.mfTag,
      HLTPaths = cms.vstring('Flag_EcalDeadCellTriggerPrimitiveFilter'),
      throw=False
      )
process.goodVertices = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::"+options.mfTag,
      HLTPaths = cms.vstring('Flag_goodVertices'),
      throw=False
      )
process.eeBadScFilter = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::"+options.mfTag,
      HLTPaths = cms.vstring('Flag_eeBadScFilter'),
      throw=False
      )


process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")


if options.runOnMC:
    process.p = cms.Path(
      #process.triggerSelection *
      #process.HBHENoiseFilter *
      #process.HBHENoiseIsoFilter *
      #process.globalTightHalo2016Filter *
      #process.EcalDeadCellTriggerPrimitiveFilter *
      #process.goodVertices *
      #process.eeBadScFilter *
      #process.BadPFMuonFilter *
      #process.BadChargedCandidateFilter *
      process.test
      )

else:
    process.p = cms.Path(
      process.triggerSelection *
      process.HBHENoiseFilter *
      process.HBHENoiseIsoFilter *
      process.globalTightHalo2016Filter *
      process.EcalDeadCellTriggerPrimitiveFilter *
      process.goodVertices *
      process.eeBadScFilter *
      process.BadPFMuonFilter *
      process.BadChargedCandidateFilter *
      process.test
      )

#process.out = cms.OutputModule( "PoolOutputModule"
#      , fileName = cms.untracked.string( "test.root" )
#      , SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('p')
#        )
#      #, outputCommands = cms.untracked.vstring('keep *_selectedUpdatedPatJets*_*_*')
#      )
#from Configuration.EventContent.EventContent_cff import MINIAODSIMEventContent
#process.out.outputCommands = MINIAODSIMEventContent.outputCommands
#process.out.outputCommands.append('keep *_selectedUpdatedPatJets*_*_*')

#process.outpath = cms.EndPath( process.out )

