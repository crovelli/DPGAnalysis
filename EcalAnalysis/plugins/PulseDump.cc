#include "PulseDump.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CommonTools/Utils/interface/StringToEnumValue.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <iostream>
#include <TLorentzVector.h>

using namespace std;
using namespace edm;
using namespace reco;

// NB:
// We assume MC is GEN-SIM-DIGI-RAW-HLTDEBUG
// and data is RAW-RECO

PulseDump::PulseDump(const edm::ParameterSet& conf)
{
  // Gen particles collection for MC, empty for data
  mcProducer_ = conf.getUntrackedParameter<std::string>("mcProducer");

  // Taking collections / parameters from cfg
  ebDigiCollectionToken_ = consumes<EBDigiCollection>(InputTag("ecalDigis","ebDigis"));
  eeDigiCollectionToken_ = consumes<EEDigiCollection>(InputTag("ecalDigis","eeDigis"));
  ecalHitEBToken_ = consumes<EcalRecHitCollection>(conf.getParameter<edm::InputTag>("reducedBarrelRecHitCollection"));
  ecalHitEEToken_ = consumes<EcalRecHitCollection>(conf.getParameter<edm::InputTag>("reducedEndcapRecHitCollection"));

  // Traslate string representation of flagsMapDBReco into enum values
  const ParameterSet & p=conf.getParameter< ParameterSet >("flagsMapDBReco");
  std::vector<std::string> recoflagbitsStrings = p.getParameterNames();
  v_DB_reco_flags_.resize(32);

  for (unsigned int i=0;i!=recoflagbitsStrings.size();++i){
    EcalRecHit::Flags recoflagbit = (EcalRecHit::Flags)
      StringToEnumValue<EcalRecHit::Flags>(recoflagbitsStrings[i]);
    std::vector<std::string> dbstatus_s =
      p.getParameter<std::vector<std::string> >(recoflagbitsStrings[i]);
    std::vector<uint32_t> dbstatuses;
    for (unsigned int j=0; j!= dbstatus_s.size(); ++j){
      EcalChannelStatusCode::Code dbstatus = (EcalChannelStatusCode::Code)
	StringToEnumValue<EcalChannelStatusCode::Code>(dbstatus_s[j]);
      dbstatuses.push_back(dbstatus);
    }
    v_DB_reco_flags_[recoflagbit]=dbstatuses;
  }

  flagmask_=0;
  flagmask_|= 0x1<<EcalRecHit::kNeighboursRecovered;
  flagmask_|= 0x1<<EcalRecHit::kTowerRecovered;
  flagmask_|= 0x1<<EcalRecHit::kDead;
  flagmask_|= 0x1<<EcalRecHit::kKilled;
  flagmask_|= 0x1<<EcalRecHit::kTPSaturated;
  flagmask_|= 0x1<<EcalRecHit::kL1SpikeFlag;
  flagmask_|= 0x1<<EcalRecHit::kNoisy;
  flagmask_|= 0x1<<EcalRecHit::kFaultyHardware;


  Service<TFileService> fs;
  _tree = fs->make<TTree>("pulse_tree","");
  _tree->Branch("run",&run_);
  _tree->Branch("event",&event_);
  _tree->Branch("barrel",&barrel_);
  _tree->Branch("gain",&gain_,"gain[10]/I");     
  _tree->Branch("pulse",&pulse_,"pulse[10]/D");
  _tree->Branch("firstGZ",&firstGZ_);
  _tree->Branch("firstG1",&firstG1_);
  _tree->Branch("firstG6",&firstG6_);
  _tree->Branch("lastUnsatSample",&lastUnsatSample_);
  _tree->Branch("ietaix",&ietaix_);
  _tree->Branch("iphiiy",&iphiiy_);
  _tree->Branch("iz",&iz_);
  _tree->Branch("mcEta",&mcEta_);
  _tree->Branch("mcEt", &mcEt_);
  _tree->Branch("kSaturated",&kSaturated_);      
  _tree->Branch("kLErecovered",&kLErecovered_);  
  _tree->Branch("kOutOfTime",&kOutOfTime_);  
  _tree->Branch("hitTime",&hitTime_);  
  _tree->Branch("swissCross",&swissCross_);

  // to count
  entry=0;
}

void PulseDump::analyze(const Event& e, const EventSetup& es) {

  entry++;
  if (entry%500==0) cout << entry << endl;

  // run and event
  run_   = e.eventAuxiliary().run();
  event_ = e.eventAuxiliary().event();

  // channel status to get flags
  es.get<EcalChannelStatusRcd>().get(chStatus);

  // geometry
  ESHandle<CaloGeometry> geoHandle;
  es.get<CaloGeometryRecord>().get(geoHandle);
  const CaloGeometry& geometry = *geoHandle;
  const CaloSubdetectorGeometry *geometry_pEB;
  const CaloSubdetectorGeometry *geometry_pEE;
  geometry_pEB = geometry.getSubdetectorGeometry(DetId::Ecal, EcalBarrel); 
  geometry_pEE = geometry.getSubdetectorGeometry(DetId::Ecal, EcalEndcap); 

  // topology
  edm::ESHandle<CaloTopology> theCaloTopology;
  es.get<CaloTopologyRecord>().get(theCaloTopology);
  theSubdetTopologyEB_ = theCaloTopology->getSubdetectorTopology(DetId::Ecal,EcalBarrel);
  
  // For data only: rechits
  const EBRecHitCollection* hitsEB=0;
  const EBRecHitCollection* hitsEE=0;
  if (mcProducer_.empty()) {
    Handle< EcalRecHitCollection > EcalBarrelRecHits;
    Handle< EcalRecHitCollection > EcalEndcapRecHits;
    e.getByToken(ecalHitEBToken_, EcalBarrelRecHits);
    e.getByToken(ecalHitEEToken_, EcalEndcapRecHits);
    hitsEB = EcalBarrelRecHits.product();
    hitsEE = EcalEndcapRecHits.product();
  }

  // MC truth study
  float mc1px = -999; float mc1py = -999.; float mc1pz = -999.; float mc1ene = -999.;
  float mc2px = -999; float mc2py = -999.; float mc2pz = -999.; float mc2ene = -999.;
  if (!mcProducer_.empty()) {
    
    bool foundFirst  = false;
    bool foundSecond = false;

    Handle<GenParticleCollection> genPcHandle;
    e.getByLabel(mcProducer_, genPcHandle);
    const GenParticleCollection *genPcColl = genPcHandle.product();
    
    for (GenParticleCollection::const_iterator p = genPcColl->begin(); p != genPcColl->end(); ++p) {

      if (foundFirst && foundSecond) break;

      int mothId = -100;
      const reco::Candidate & cand = *p;
      GenParticleCollection::const_iterator p2;
      for (p2 = genPcColl->begin(); p2 != genPcColl->end(); ++p2) {
	const reco::Candidate *mom = cand.mother();
	if(&(*p2)==&(*mom)) {
	  mothId = p2->pdgId();
	  break;
	}
      }

      if (p->pdgId()==22 && mothId==5100039) {
	if (!foundFirst) {
	  mc1px = (*p).px();
	  mc1py = (*p).py();
	  mc1pz = (*p).pz();
	  mc1ene = (*p).energy();
	  foundFirst = true;
	} else {
	  mc2px = (*p).px();
	  mc2py = (*p).py();
	  mc2pz = (*p).pz();
	  mc2ene = (*p).energy();
	  foundSecond = true;
	}
      }
    }
  }

  // taking digis and do the analysis
  Handle< EBDigiCollection > pEBDigis;
  Handle< EEDigiCollection > pEEDigis;
  e.getByToken( ebDigiCollectionToken_, pEBDigis);
  e.getByToken( eeDigiCollectionToken_, pEEDigis);
  for(EBDigiCollection::const_iterator itdg = pEBDigis->begin(); itdg != pEBDigis->end(); ++itdg) {
    if (!mcProducer_.empty()) FillRecHitMC(*itdg,mc1px,mc1py,mc1pz,mc1ene,mc2px,mc2py,mc2pz,mc2ene,geometry_pEB);
    else FillRecHitData(*itdg,hitsEB,hitsEE);
  }
  for(EEDigiCollection::const_iterator itdg = pEEDigis->begin(); itdg != pEEDigis->end(); ++itdg) {
    if (!mcProducer_.empty()) FillRecHitMC(*itdg,mc1px,mc1py,mc1pz,mc1ene,mc2px,mc2py,mc2pz,mc2ene,geometry_pEE);
    else FillRecHitData(*itdg,hitsEB,hitsEE);
  }
}

void PulseDump::FillRecHitData(const EcalDataFrame& dataFrame,const EBRecHitCollection* hitsEB,const EERecHitCollection* hitsEE) {

  const unsigned int nsample = EcalDataFrame::MAXSAMPLES;

  // event selection: at least one sample with G0 [ or G1 or G6 ]
  firstGZ_ = 100;
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    if (sample.gainId()==0) { 
      firstGZ_ = iSample;
      break;
    }
  }
  // if (firstGZ_>10) return;

  firstG1_ = 100;
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    if (sample.gainId()==3) { 
      firstG1_ = iSample;
      break;
    }
  }
  firstG6_ = 100;
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    if (sample.gainId()==2) { 
      firstG6_ = iSample;
      break;
    }
  }
  if (firstGZ_>10 && firstG1_>10 && firstG6_>10) return;

  lastUnsatSample_ = dataFrame.lastUnsaturatedSample();

  // detId for this digi
  DetId detid = dataFrame.id();
  barrel_ = detid.subdetId()==EcalBarrel;

  // check for channels to be excluded from reconstruction  
  EcalChannelStatusMap::const_iterator chit = chStatus->find(detid);
  EcalChannelStatusCode::Code dbstatus = chit->getStatusCode();
  if ( v_chstatus_.size() > 0) {
    std::vector<int>::const_iterator res =
      std::find( v_chstatus_.begin(), v_chstatus_.end(), dbstatus );
    if ( res != v_chstatus_.end() ) return;
  }
  uint32_t flagBits = setFlagBits(v_DB_reco_flags_, dbstatus);
  if (flagmask_ & flagBits) return;

  // dumping amplitude and gain for each sample 
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    double amplitude = (double)(sample.adc());
    gain_[iSample]  = sample.gainId();  
    pulse_[iSample] = amplitude;
  }

  if (barrel_) {
    EBDetId ebid(detid);
    ietaix_ = ebid.ieta();
    iphiiy_ = ebid.iphi();
    iz_ = 0;
  }
  else {
    EEDetId eeid(detid);
    ietaix_ = eeid.ix(); 
    iphiiy_ = eeid.iy();
    iz_ = eeid.zside();
  }

  // rechit flags and time
  EcalRecHitCollection::const_iterator itRH;
  if (barrel_) 
    itRH = hitsEB->find( detid );
  else
    itRH = hitsEE->find( detid );
  kSaturated_   = itRH->checkFlag(EcalRecHit::kSaturated);
  kLErecovered_ = itRH->checkFlag(EcalRecHit::kLeadingEdgeRecovered);
  kOutOfTime_   = itRH->checkFlag(EcalRecHit::kOutOfTime);
  hitTime_      = itRH->time();

  // swiss cross, barrel only
  swissCross_ = -999.;
  if (barrel_) {

    int iNeigh=0;
    CaloNavigator<DetId> cursorE = CaloNavigator<DetId>( detid, theSubdetTopologyEB_ );
    
    float sXsum  = 0.;
    float theMax = 0.;
    for(int ix=-1; ix<2; ++ix) {
      for(int iy=-1; iy<2; ++iy) {
	cursorE.home();
	cursorE.offsetBy( ix, iy );
	DetId cryId = cursorE.pos();

	if(cryId.subdetId()!=EcalBarrel) continue; 
	if (abs(ix)==1 && abs(iy)==1)    continue;
	EcalRecHitCollection::const_iterator itneigh = hitsEB->find( cryId );
	if( itneigh != hitsEB->end() ) { 
	  if (ix!=0 || iy!=0) sXsum  = sXsum + itneigh->energy();
	  if (ix==0 && iy==0) theMax = itneigh->energy(); 
	  // cout << "ix = " << ix << ", iy = " << iy << ", sXsum = " << sXsum << ", theMax = " << theMax << endl;
	  iNeigh++;
	}
      }
    }
    if (iNeigh!=5) cout << "problem: not 5 crystals!  ==> " << iNeigh << endl;
    if (theMax>0 && iNeigh==5) swissCross_ = 1. - (sXsum/theMax);
  }
      
  // MC variables dummy
  mcEta_ = -999.;
  mcEt_  = -999.;

  // Filling the tree
  _tree->Fill();
}

void PulseDump::FillRecHitMC(const EcalDataFrame& dataFrame,float mc1px,float mc1py,float mc1pz, float mc1ene,float mc2px,float mc2py,float mc2pz, float mc2ene,const CaloSubdetectorGeometry*geometry) {

  const unsigned int nsample = EcalDataFrame::MAXSAMPLES;

  // event selection: at least one sample with G0 [ or G1 or G6 ] 
  firstGZ_ = 100;
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    if (sample.gainId()==0) { 
      firstGZ_ = iSample;
      break;
    }
  }
  // if (firstGZ_>10) return;

  firstG1_ = 100;
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    if (sample.gainId()==3) { 
      firstG1_ = iSample;
      break;
    }
  }
  firstG6_ = 100;
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    if (sample.gainId()==2) { 
      firstG6_ = iSample;
      break;
    }
  }
  if (firstGZ_>10 && firstG1_>10 && firstG6_>10) return;

  // DetId for this digi
  DetId detid = dataFrame.id();
  barrel_ = detid.subdetId()==EcalBarrel;

  // check for channels to be excluded from reconstruction  
  EcalChannelStatusMap::const_iterator chit = chStatus->find(detid);
  EcalChannelStatusCode::Code dbstatus = chit->getStatusCode();
  if ( v_chstatus_.size() > 0) {
    std::vector<int>::const_iterator res =
      std::find( v_chstatus_.begin(), v_chstatus_.end(), dbstatus );
    if ( res != v_chstatus_.end() ) return;
  }
  uint32_t flagBits = setFlagBits(v_DB_reco_flags_, dbstatus);
  if (flagmask_ & flagBits) return;

  // dumping amplitude and gain for each sample
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    double amplitude = (double)(sample.adc());
    gain_[iSample]  = sample.gainId();  
    pulse_[iSample] = amplitude;
  }

  if (barrel_) {
    EBDetId ebid(detid);
    ietaix_ = ebid.ieta();
    iphiiy_ = ebid.iphi();
    iz_ = 0;
  }
  else {
    EEDetId eeid(detid);
    ietaix_ = eeid.ix(); 
    iphiiy_ = eeid.iy();
    iz_ = eeid.zside();
  }

  // closest mc truth gamma
  mcEta_ = -999.;
  mcEt_  = -999.;

  // position in CMS coordinates to compare with mc truth
  const CaloCellGeometry *this_cell;
  this_cell = (*geometry).getGeometry(dataFrame.id());
  GlobalPoint position = this_cell->getPosition();
  float cellEta = position.eta();
  float cellPhi = position.phi();
  
  // mc truth
  TLorentzVector theTLV1( mc1px, mc1py, mc1pz, mc1ene );
  TLorentzVector theTLV2( mc2px, mc2py, mc2pz, mc2ene );
  float mc1Eta = theTLV1.Eta();
  float mc1Phi = theTLV1.Phi();
  float mc2Eta = theTLV2.Eta();
  float mc2Phi = theTLV2.Phi();
  
  // choosing the closest
  float deltaEta1 = fabs(cellEta - mc1Eta);
  float deltaPhi1 = fabs(cellPhi - mc1Phi);
  if (deltaPhi1 > Geom::pi())  deltaPhi1 -= 2.*Geom::pi();
  if (deltaPhi1 < -Geom::pi()) deltaPhi1 += 2.*Geom::pi();
  
  float deltaEta2 = fabs(cellEta - mc2Eta);
  float deltaPhi2 = fabs(cellPhi - mc2Phi);
  if (deltaPhi2 > Geom::pi())  deltaPhi2 -= 2.*Geom::pi();
  if (deltaPhi2 < -Geom::pi()) deltaPhi2 += 2.*Geom::pi();
  
  float deltaR1 = sqrt(deltaEta1*deltaEta1+deltaPhi1*deltaPhi1);
  float deltaR2 = sqrt(deltaEta2*deltaEta2+deltaPhi2*deltaPhi2);
  
  bool match1 = false;
  bool match2 = false;
  if (deltaR1<deltaR2 && deltaR1<0.3) match1 = true;
  if (deltaR2<deltaR1 && deltaR2<0.3) match2 = true;
  
  if (match1) {
    mcEta_ = theTLV1.Eta();
    mcEt_  = theTLV1.Pt();
  } else if (match2) {
    mcEta_ = theTLV2.Eta();
    mcEt_  = theTLV2.Pt();
  }

  // data variables, dummy
  kSaturated_   = 0;
  kLErecovered_ = 0;
  kOutOfTime_   = 0;
  hitTime_      = -999999.;

  _tree->Fill();
}

// Take our association map of dbstatuses-> recHit flagbits and return the apporpriate flagbit word
uint32_t PulseDump::setFlagBits(const std::vector<std::vector<uint32_t> >& map,
				const uint32_t& status ){
  for (unsigned int i = 0; i!=map.size(); ++i){
    if (std::find(map[i].begin(), map[i].end(),status)!= map[i].end())
      return 0x1 << i;
  }
  return 0;
}
