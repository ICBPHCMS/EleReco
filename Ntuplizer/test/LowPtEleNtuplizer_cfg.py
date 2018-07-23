import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from EleReco.Ntuplizer.files_RERECO import *

process = cms.Process("Ntuplizer")

options = VarParsing.VarParsing('analysis')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
options.parseArguments()

sample = 4 #@@ choose preexisting samples here

files = [
    options.inputFiles,
    BToKee_Seed2p0,
    BToKee_Seed1p0,
    BToKee_Seed0p5,
    BToKee_Gsf0p5,
    BToKmm_Seed2p0,
    BToKmm_Seed0p5,
    ][sample]
output = [options.outputFile,
          'output_BToKee_Seed2p0.root',
          'output_BToKee_Seed1p0.root',
          'output_BToKee_Seed0p5.root',
          'output_BToKee_Gsf0p5.root',
          'output_BToKmm_Seed2p0.root'
          'output_BToKmm_Seed0p5.root'
          ][sample]

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(*files)
                            )

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(10)
    input = cms.untracked.int32(options.maxEvents)
)

process.load('EleReco.Ntuplizer.LowPtEleNtuplizer_cfi')
process.p = cms.Path(process.LowPtEleNtuplizer)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100 # silence output

process.TFileService=cms.Service('TFileService',
                                 fileName=cms.string(output)
                                 )
