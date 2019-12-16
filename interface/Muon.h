#ifndef Analysis_Core_Muon_h
#define Analysis_Core_Muon_h

#include "TLorentzVector.h"
#ifndef STANDALONE
#include "HMuMu/Core/interface/Track.h"
#else
#include "Track.h"
#endif

#include <math.h>       /* fabs */

namespace analysis
{
namespace core
{
class Muon : public Track
{
  public:
	Muon() : Track() { this->reset(); }

	virtual void reset()
	{
		_isTracker = 0;
		_isStandAlone = 0;
		_isGlobal = 0;
		_isTight = 0;
		_isMedium = 0;
		_isLoose = 0;

		_isPF = 0;

		// _pt_kinfit = 0; 
		// _ptErr_kinfit = 0;
		// _d0PV_kinfit = 0;
		// _dzPV_kinfit = 0;
		// _chi2_kinfit = 0;
		// _ndf_kinfit = 0;

		fsrP4.SetPxPyPzE(0,0,0,0);

		_pt_PF = 0;
		_pterr_PF = 0;
		_eta_PF = 0;
		_phi_PF = 0;

		_normChi2 = 0;
		_d0BS = 0;
		_dzBS = 0;
		_d0PV = 0;
		_dzPV = 0;
		_ip3d = 0;
		_sip3d = 0;

		_nPLs = 0;
		_nTLs = 0;
		_nSLs = 0;

		_vfrTrk = 0;
		_nvMHits = 0;
		_nvPHits = 0;
		_nvTHits = 0;
		_nvSHits = 0;
		_nSegMts = 0;
		_nMtsStations = 0;

		_trackIsoSumPt = 0;
		_trackIsoSumPtCorr = 0;

		_hIso = 0;
		_eIso = 0;
		_relCombIso = 0;
		_pfIso = 0;
		_miniIso = 0;
		
		Track::reset();
		_track.reset();
		_isHLTMatched.clear();

		_segmentCompatibility = 0;
		_combinedQChi2LocalPosition = 0;
		_combinedQTrkKink = 0;

		_sumChargedHadronPtR03 = 0;
		_sumChargedParticlePtR03 = 0;
		_sumNeutralHadronEtR03 = 0;
		_sumPhotonEtR03 = 0;
		_sumPUPtR03 = 0;
		_sumChargedHadronPtR04 = 0;
		_sumChargedParticlePtR04 = 0;
		_sumNeutralHadronEtR04 = 0;
		_sumPhotonEtR04 = 0;
		_sumPUPtR04 = 0;
		_roccCor = 0;
		_corrPT = 0;
		_geoPT = 0;
		
	}
	virtual ~Muon() {}

	bool _isTracker;
	Track _track;
	bool _isStandAlone;
	bool _isGlobal;
	bool _isTight;
	bool _isMedium;
	bool _isLoose;

	bool _isPF;
	

	float _normChi2;
	float _d0BS;
	float _dzBS;

	float _d0PV;
	float _dzPV;
	float _ip3d;
	float _sip3d;


	int _nPLs;
	int _nTLs;
	int _nSLs;

	float _vfrTrk;
	int _nvMHits;
	int _nvPHits;
	int _nvTHits;
	int _nvSHits;
	int _nSegMts;
	int _nMtsStations;

	float _trackIsoSumPt;
	float _trackIsoSumPtCorr;

	float _hIso;
	float _eIso;
	float _relCombIso;
	float _pfIso;
	float _miniIso;

	float _segmentCompatibility;
	float _combinedQChi2LocalPosition;
	float _combinedQTrkKink;
	std::vector<bool> _isHLTMatched;

	float _sumChargedHadronPtR03;
	float _sumChargedParticlePtR03;
	float _sumNeutralHadronEtR03;
	float _sumPhotonEtR03;
	float _sumPUPtR03;
	float _sumChargedHadronPtR04;
	float _sumChargedParticlePtR04;
	float _sumNeutralHadronEtR04;
	float _sumPhotonEtR04;
	float _sumPUPtR04;
	float _roccCor;
	float _corrPT;
	float _geoPT;

	// float _d0PV_kinfit;
	// float _dzPV_kinfit;
	// double _pt_kinfit; 
	// double _ptErr_kinfit; 
	// float _chi2_kinfit;
	// int _ndf_kinfit;

	double _pt_PF;
	double _pterr_PF;
	double _eta_PF;
	double _phi_PF;

	TLorentzVector fsrP4;

	bool operator==(const Muon& m2)
	{
		return (fabs(this->_pt - m2._pt) < this->_pt * .0001 && 
				fabs(this->_eta - m2._eta) < this->_eta * .0001 && 
				fabs(this->_phi - m2._phi) < this->_phi * .0001 && 
				fabs(this->_charge - m2._charge) < this->_charge * .0001);
	}

#ifdef STANDALONE
	ClassDef(Muon, 1)
#endif
};

typedef std::vector<analysis::core::Muon> Muons;
} // namespace core
} // namespace analysis

#ifdef STANDALONE
ClassImpUnique(analysis::core::Muon, Muon)
#endif

#endif
