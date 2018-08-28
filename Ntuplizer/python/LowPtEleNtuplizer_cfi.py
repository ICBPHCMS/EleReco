import FWCore.ParameterSet.Config as cms

LowPtEleNtuplizer = cms.EDAnalyzer("LowPtEleNtuplizer",
                                   hepMCProduct = cms.InputTag("generatorSmeared"),
                                   pileup = cms.InputTag("addPileupInfo"),
                                   vertices = cms.InputTag("offlinePrimaryVertices"),
                                   genParticles = cms.InputTag("genParticles"),
                                   generalTracks = cms.InputTag("generalTracks"),
                                   # recoPreIds_trackerDrivenElectronSeeds_preid_RECO
                                   preIds = cms.InputTag("trackerDrivenElectronSeeds","preid"),
                                   # recoElectronSeeds_electronMergedSeeds__RECO
                                   # recoElectronSeeds_electronMergedSeeds__RECO
                                   electronSeeds = cms.InputTag("electronMergedSeeds"), 
                                   gsfTracks = cms.InputTag("electronGsfTracks"),
                                   gsfElectronCores = cms.InputTag("gedGsfElectronCores"),
                                   gsfElectronsTmp = cms.InputTag("gedGsfElectronsTmp"),
                                   gsfElectrons = cms.InputTag("gedGsfElectrons"),
                                   recoMuons = cms.InputTag("muons"),
                                   )


