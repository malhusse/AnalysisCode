#ifndef Analysis_Core_Electron_h
#define Analysis_Core_Electron_h

#ifndef STANDALONE
#include "HMuMu/Core/interface/Track.h"
#else
#include "Track.h"
#endif

namespace analysis
{
namespace core
{
class Electron : public Track
{
  public:
    Electron() : Track() { this->reset(); }
    virtual ~Electron() {}
    virtual void reset()
    {
        Track::reset();
        _sumChargedHadronPt = 0;
        _sumNeutralHadronEt = 0;
        _sumPhotonEt = 0;
        _sumPUPt = 0;
        _sumChargedParticlePt = 0;

        _dzPV = 0;
        _d0PV = 0;
        _ip3d = 0;
        _sip3d = 0;

        // _relCombIso = 0;
		_miniIso = 0;

        _isPF = false;
        _convVeto = false;

        _isTight = 0;
	    _isMedium = 0;
	    _isLoose = 0;
        _isVeto = 0;
    }

    float _sumChargedHadronPt;
    float _sumNeutralHadronEt;
    float _sumPhotonEt;
    float _sumPUPt;
    float _sumChargedParticlePt;

    float _dzPV;
    float _d0PV;
    float _ip3d;
    float _sip3d;

	// float _relCombIso;
	float _miniIso;
    
    bool _isPF;
    bool _convVeto;

    bool _isTight;
	bool _isMedium;
	bool _isLoose;
    bool _isVeto;

#ifdef STANDALONE
    ClassDef(Electron, 1)
#endif
};

typedef std::vector<analysis::core::Electron> Electrons;
} // namespace core
} // namespace analysis

#ifdef STANDALONE
ClassImpUnique(analysis::core::Electron, Electron)
#endif

#endif
