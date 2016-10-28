import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import sys

options = VarParsing ('analysis')

options.setDefault( 'outputFile',
      'ntuple.root'
      )

options.register( 'globalTag',
      '80X_dataRun2_Candidate_2016_09_02_10_26_48',
      #'74X_dataRun2_Prompt_v1',
      #'MCRUN2_74_V9',
      #'GR_P_V56',
      VarParsing.multiplicity.singleton,
      VarParsing.varType.string,
      "CMS Global Tag"
      )

options.register( 'runOnMC',
      False,
      VarParsing.multiplicity.singleton,
      VarParsing.varType.bool,
      "mc or data"
      )

options.register( 'mfTag',
      'PAT',
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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# print statistics
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #'/store/data/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/00039A2E-D7A7-E511-98EE-3417EBE64696.root'
       '/store/data/Run2016G/DoubleMuon/MINIAOD/23Sep2016-v1/100000/00993A51-DF90-E611-A4EE-7845C4FC3650.root'
       #'/store/data/Run2016G/DoubleMuon/MINIAOD/23Sep2016-v1/100000/084F88CC-548F-E611-BEED-549F35AD8B7B.root'
       #'/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/02A85EE9-70BA-E511-A0A2-0CC47A4D7678.root'
       #'file:00039A2E-D7A7-E511-98EE-3417EBE64696.root'
    )
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
      src = cms.InputTag("packedPFCandidates"),
      #jets = cms.InputTag("slimmedJets"),
      jets = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
      leptons = cms.VInputTag("slimmedElectrons", "slimmedMuons", "slimmedPhotons"),
      met = cms.InputTag("slimmedMETs"),
      muons = cms.InputTag("slimmedMuons"),
      vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
      runOnMC = cms.untracked.bool(options.runOnMC),
      addPileupInfo = cms.InputTag("slimmedAddPileupInfo"),
      rho = cms.InputTag("fixedGridRhoAll"),
      generator = cms.InputTag("generator"),
      lheprod = cms.InputTag("externalLHEProducer"),
      pileup = cms.untracked.InputTag("slimmedAddPileupInfo"),
      srcJetSF             = cms.string('AK4PFchs'),
      srcJetResPt          = cms.string('AK4PFchs_pt'),
      srcJetResPhi         = cms.string('AK4PFchs_phi'),
)

process.load('Configuration.StandardSequences.Services_cff')
#process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *

JERdataMC = 'MC'
if not options.runOnMC:
   JERdataMC = 'DATA'

#process.jer = cms.ESSource("PoolDBESSource",
#      CondDBSetup,
#      toGet = cms.VPSet(
#         # Pt Resolution
#         cms.PSet(
#            record = cms.string('JetResolutionRcd'),
#            tag    = cms.string('JR_Fall15_25nsV2_'+JERdataMC+'_PtResolution_AK4PFchs'),
#            label  = cms.untracked.string('AK4PFchs_pt')
#            ),
#         # Phi Resolution
#         cms.PSet(
#            record = cms.string('JetResolutionRcd'),
#            tag    = cms.string('JR_Fall15_25nsV2_'+JERdataMC+'_PhiResolution_AK4PFchs'),
#            label  = cms.untracked.string('AK4PFchs_phi')
#            ),
#         # Scale factors
#         cms.PSet(
#            record = cms.string('JetResolutionScaleFactorRcd'),
#            tag    = cms.string('JR_Fall15_25nsV2_'+JERdataMC+'_SF_AK4PFchs'),
#            label  = cms.untracked.string('AK4PFchs')
#            ),
#         ),
#      connect = cms.string('sqlite:Fall15_25nsV2_'+JERdataMC+'.db')
#      )
#
#process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

# trigger filter                
trigger_paths = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v']
trigger_pattern = [path+"*" for path in trigger_paths]
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.triggerSelection = hltHighLevel.clone(
      TriggerResultsTag = "TriggerResults::HLT",
      HLTPaths = trigger_pattern,
      throw=False
      )


##___________________________CSC_Halo_Filter__________________________________||
process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
process.load('RecoMET.METFilters.CSCTightHalo2015Filter_cfi')
process.load('RecoMET.METFilters.CSCTightHaloTrkMuUnvetoFilter_cfi')

process.p = cms.Path(process.CSCTightHaloFilter)
process.p = cms.Path(process.CSCTightHalo2015Filter)
process.p = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)

##___________________________Global_Halo_Filter__________________________________||
process.load('RecoMET.METFilters.globalTightHalo2016Filter_cfi')
process.load('RecoMET.METFilters.globalSuperTightHalo2016Filter_cfi')

process.p = cms.Path(process.globalTightHalo2016Filter)
process.p = cms.Path(process.globalSuperTightHalo2016Filter)


##___________________________HCAL_Noise_Filter________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)

process.p = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseFilter)
process.p = cms.Path(process.HBHENoiseFilterResultProducer * process.HBHENoiseIsoFilter)

#process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('RecoMET.METFilters.HcalStripHaloFilter_cfi')
#process.p = cms.Path(process.HcalStripHaloFilter)



process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)
process.p = cms.Path(process.goodVertices)

process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.p = cms.Path(process.trackingFailureFilter)

process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.p = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)

process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.p = cms.Path(process.eeBadScFilter)

#process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
#process.p = cms.Path(process.BadChargedCandidateFilter)
#process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

#process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#process.BadPFMuonFilter.taggingMode = cms.bool(True)
##process.BadPFMuonFilter.debug = cms.bool(True)

#process.load('RecoMET.METFilters.BadChargedCandidateSummer16Filter_cfi')
#process.BadChargedCandidateSummer16Filter.taggingMode = cms.bool(True)
#
#process.load('RecoMET.METFilters.BadPFMuonSummer16Filter_cfi')
#process.BadPFMuonSummer16Filter.taggingMode = cms.bool(True)
##process.BadPFMuonFilter.debug = cms.bool(True)

#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
#                                           minimumNDOF = cms.uint32(4) ,
#                                           maxAbsZ = cms.double(24),
#                                           maxd0 = cms.double(2)
#                                           )
#process.p = cms.Path(process.primaryVertexFilter)

process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
##process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(True)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB=cms.vint32(12, 13, 14)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE=cms.vint32(12, 13, 14)
process.p = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)



process.p = cms.Path(
      process.triggerSelection *
      process.HBHENoiseFilter *
      process.HBHENoiseIsoFilter *
#      process.CSCTightHaloFilter *
#      process.CSCTightHalo2015Filter *
#      process.CSCTightHaloTrkMuUnvetoFilter *
      process.globalTightHalo2016Filter *
#      process.globalSuperTightHalo2016Filter *
#      process.goodVertices *
#      process.trackingFailureFilter *
#      process.primaryVertexFilter *
#      process.EcalDeadCellBoundaryEnergyFilter *
#      process.HBHENoiseFilterResultProducer *
#      process.HcalStripHaloFilter *
      process.EcalDeadCellTriggerPrimitiveFilter *
#      process.eeBadScFilter *
#      process.chargedHadronTrackResolutionFilter *
#      process.muonBadTrackFilter *
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

