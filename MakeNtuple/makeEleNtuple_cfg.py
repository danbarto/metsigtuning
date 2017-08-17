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
      '80X_mcRun2_asymptotic_2016_TrancheIV_v8',
      #'80X_dataRun2_2016SeptRepro_v7',

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
process.load("metsigtuning.MakeNtuple.makeEleNtuple_cfi")
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


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #'/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00EFB35A-A8D5-E611-A740-001E67A3E872.root'
       #'/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/0EA60289-18C4-E611-8A8F-008CFA110AB4.root'
       #'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/0694F787-BCD7-E611-87BE-1CB72C0A3A5D.root'
       #'/store/mc/RunIISummer16MiniAODv2/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/10A6C6DC-9ECE-E611-B6CD-0CC47AC08816.root'
       '/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/0094B6D9-F8EA-E611-A809-0CC47AD9914A.root'
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

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

el_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
#add them to the VID producer
for idmod in el_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.test = cms.EDAnalyzer('MakeEleNtuple',
      src                       = cms.InputTag("packedPFCandidates"),
      #jets                      = cms.InputTag("slimmedJets"),
      #jets                      = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
      jets                      = cms.VInputTag("cleanedPatJets"),
      leptons                   = cms.VInputTag("slimmedElectrons", "slimmedMuons", "slimmedPhotons"),
      met                       = cms.InputTag("slimmedMETs"),
      electronVetoIdMap         = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
      electronLooseIdMap        = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
      electronMediumIdMap       = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
      electronTightIdMap        = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
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
#trigger_paths = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v', 'HLT_IsoMu24_v', 'HLT_IsoTkMu24_v']
trigger_paths = ['HLT_Ele27_eta2p1_WPTight_Gsf_v']
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
      process.egmGsfElectronIDSequence * 
      #process.photonIDValueMapProducer *
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
      process.egmGsfElectronIDSequence *
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

