#ifndef Analysis_Core_Event_h
#define Analysis_Core_Event_h

#ifndef STANDALONE
#include "HMuMu/Core/interface/Object.h"
#else
#include "Object.h"
#endif

namespace analysis
{
namespace core
{
class Event : public Object
{
  public:
	Event() : Object() { this->reset(); }

	virtual void reset()
	{
		_run = 0;
		_lumi = 0;
		_event = 0;
		_bx = 0;
		_orbit = 0;
	}
	virtual ~Event() {}

	int _run;
	int _lumi;
	long long int _event;
	int _bx;
	int _orbit;

#ifdef STANDALONE
	ClassDef(Event, 1)
#endif
};

class EventAuxiliary : public Object
{
  public:
	EventAuxiliary() : Object() { this->reset(); }
	virtual ~EventAuxiliary() {}

	virtual void reset()
	{
		_nPU = 0;
		_bareMCWeight = 1;
		_genWeight = 1;
		_hasHLTFired.clear();
		_metFilterBits.clear();
		_passedMetFilters = true;
		_prefiringweight = 1;
		_prefiringweightup = 1;
		_prefiringweightdown = 1;
		_trigEffSF = 1;
		_trigEffSF_up = 1;
		_trigEffSF_down = 1;
		_idSF = 1;
		_idSF_up = 1;
		_idSF_down = 1;
		_isoSF = 1;
		_isoSF_up = 1;
		_isoSF_down = 1;
	}	

	int _nPU;
	int _nvtx;
	float _bareMCWeight;
	int _genWeight;
	bool _passedMetFilters;
	double _prefiringweight;
	double _prefiringweightup;
	double _prefiringweightdown;
	float _trigEffSF;
	float _trigEffSF_up;
	float _trigEffSF_down;
	float _idSF;
	float _idSF_up;
	float _idSF_down;
	float _isoSF;
	float _isoSF_up;
	float _isoSF_down;

	std::vector<bool> _hasHLTFired;
	std::map<string, int> _metFilterBits;



#ifdef STANDALONE
	ClassDef(EventAuxiliary, 1)
#endif
};

typedef std::vector<analysis::core::EventAuxiliary> EventAuxiliaries;
typedef std::vector<analysis::core::Event> Events;
} // namespace core
} // namespace analysis

#ifdef STANDALONE
ClassImpUnique(analysis::core::Event, Event)
ClassImpUnique(analysis::core::EventAuxiliary, EventAuxiliary)
#endif

#endif
