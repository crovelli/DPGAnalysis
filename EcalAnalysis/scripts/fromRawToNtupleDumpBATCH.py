import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO2")

# import of standard configurations        
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
 
# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'MCRUN2_74_V7::All'

# Logs
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 200 )
process.MessageLogger.cerr.threshold = 'ERROR' 
 
# start from RAW format for more flexibility
process.raw2digi_step = cms.Sequence(process.RawToDigi)

# events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('pippo.root'
                                                              )

process.ecalRtoD = cms.Sequence( process.raw2digi_step )

process.TFileService = cms.Service("TFileService",fileName = cms.string("OUTPUT"))
process.load('DPGAnalysis.EcalAnalysis.pulseDump_cfi')
#Setting to null value avoids reading mc infos  
process.pulseDump.mcProducer = cms.untracked.string('')  

process.mydump = cms.Sequence(process.pulseDump)

process.myPath = cms.Path( process.ecalRtoD * process.mydump )

