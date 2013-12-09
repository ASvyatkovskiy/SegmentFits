import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process('Rec2')

#from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.ReMixingSeeds_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_5_2_0_pre6/RelValSingleMuPt1000/GEN-SIM-DIGI-RAW-HLTDEBUG/START52_V2-v3/0110/8C36B814-7360-E111-81FE-0026189438B1.root',
 )
)

import SimGeneral.MixingModule.mixNoPU_cfi
process.muonMix = SimGeneral.MixingModule.mixNoPU_cfi.mix.clone()

process.muonMix.playback = True
process.GlobalTag.globaltag = 'START52_V2::All'

# This is to load new CondDB
from CondCore.DBCommon.CondDBSetup_cfi import *

# Reading APEs 
process.muonAPEs = cms.ESSource("PoolDBESSource", CondDBSetup, 
 connect = cms.string('sqlite_file:/uscms/home/asvyatko/nobackup/initialMuonAlignment_DT-Aug11_CSC-Aug12_SetAPE.db'), 
 toGet = cms.VPSet(cms.PSet(record = cms.string("DTAlignmentRcd"), 
 tag = cms.string("DTAlignmentRcd")), 
 cms.PSet(record = cms.string("DTAlignmentErrorRcd"), 
 tag = cms.string("DTAlignmentErrorRcd")), 
 cms.PSet(record = cms.string("CSCAlignmentRcd"), 
 tag = cms.string("CSCAlignmentRcd")), 
 cms.PSet(record = cms.string("CSCAlignmentErrorRcd"), 
 tag = cms.string("CSCAlignmentErrorRcd")))) 
process.es_prefer_muonAPEs = cms.ESPrefer("PoolDBESSource","muonAPEs")

if 'debugdumps' in sys.argv:
    process.MessageLogger.cerr.threshold = 'DEBUG'
    process.MessageLogger.debugModules = cms.untracked.vstring('tevMuons')
    process.MessageLogger.categories.append('TrackFitters')
    process.MessageLogger.categories.append('GlobalMuonRefitter')

process.GLBsegsInFits = cms.EDAnalyzer('SegsInFits',
                                       srcs = cms.VInputTag(cms.InputTag("globalMuons", "", "Rec2"),
                                                            cms.InputTag("globalMuons", "", "Rec2")
                                                            ),
                                       is_StA_run =  cms.bool(False),
                                       doPhiCorr = cms.bool(False),
                                       signSelection = cms.double(1.),
                                       )
process.TFileService = cms.Service('TFileService', fileName=cms.string('histos.root'))

process.p = cms.Path(process.RawToDigi*process.reconstruction*process.GLBsegsInFits)
