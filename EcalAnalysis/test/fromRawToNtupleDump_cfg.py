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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )

# MC
#path = '/store/relval/CMSSW_7_4_0_pre8/RelValRSGravitonToGaGa_13TeV/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_74_V7-v1/00000/'
#process.source = cms.Source("PoolSource",
#                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
#                            fileNames = cms.untracked.vstring(path+'143BB331-32BD-E411-A793-0025905A6118.root',
#                                                              path+'58708635-32BD-E411-957B-003048FFD7BE.root',
#                                                              path+'7099AF75-64BD-E411-AEC1-0025905A6088.root',
#                                                              path+'7EE5C235-32BD-E411-A576-0025905A607A.root',
#                                                              path+'8254F832-32BD-E411-AA58-0025905A6118.root',
#                                                              path+'A026D7BF-58BD-E411-A947-0025905A60B0.root',
#                                                              path+'B2B52A74-64BD-E411-B151-0025905A6092.root',
#                                                              path+'B83EE07E-3DBD-E411-AEDA-003048FFD75C.root',
#                                                              path+'D22236BF-58BD-E411-B25D-0025905A60B6.root'
#                                                              ))

# Data
from DPGAnalysis.EcalAnalysis.DoubleEle2012CZElectron_cff import *                                                      
process.source = cms.Source("PoolSource",                                                                                 
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = readFiles                                                                        
                            )                                    
#process.source = cms.Source("PoolSource",                                                                                                           # 
#                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),                                                          # 
#                            fileNames = cms.untracked.vstring('file:F6E5DF7C-6168-E211-B612-002618943829.root')
#                            )

process.ecalRtoD = cms.Sequence( process.raw2digi_step )

process.TFileService = cms.Service("TFileService",fileName = cms.string("pulseTree.root"))
process.load('DPGAnalysis.EcalAnalysis.pulseDump_cfi')
#Setting to null value avoids reading mc infos  
process.pulseDump.mcProducer = cms.untracked.string('')  

process.mydump = cms.Sequence(process.pulseDump)

process.myPath = cms.Path( process.ecalRtoD * process.mydump )

