#ifndef Analysis_Core_Vertex_h
#define Analysis_Core_Vertex_h

#ifndef STANDALONE
#include "Object.h"
#else
#include "Object.h"
#endif

namespace analysis
{
namespace core
{
class Vertex : public Object
{
  public:
	Vertex() : Object() { this->reset(); }

	virtual void reset()
	{
		_isValid = 0;
		_x = 0;
		_y = 0;
		_z = 0;
		_xerr = 0;
		_yerr = 0;
		_zerr = 0;
		_chi2 = 0;
		_ndf = 0;
		_normChi2 = 0;
	}

	virtual ~Vertex() {}

	int _isValid;
	float _x;
	float _y;
	float _z;
	float _xerr;
	float _yerr;
	float _zerr;
	float _chi2;
	float _ndf;
	float _normChi2;

#ifdef STANDALONE
	ClassDef(Vertex, 1)
#endif
};

typedef std::vector<analysis::core::Vertex> Vertices;
} // namespace core
} // namespace analysis

#ifdef STANDALONE
ClassImpUnique(analysis::core::Vertex, Vertex)
#endif

#endif
