#ifndef ECHMET_ELMIG_PARAMS_FITTER_H
#define ECHMET_ELMIG_PARAMS_FITTER_H

#include <echmetelems.h>
#include <echmetsyscomp.h>
#include <new>

namespace ECHMET {
namespace ElmigParamsFitter {

ECHMET_ST_ENUM(RetCode) {
	OK = 0,
	E_INVALID_ARGUMENT = 0x1,
	E_NOT_ENOUGH_MEASUREMENTS = 0x2,
	E_BUFFERS_EXPDATA_MISMATCH = 0x3,
	E_INVALID_BUFFER = 0x4,
	E_NO_MEMORY = 0x5,
	E_REGRESSOR_INITIALIZATION = 0x6,
	E_REGRESSOR_NO_SOLUTION = 0x7,
	E_REGRESSOR_INTERNAL_ERROR = 0x8,
	E_REGRESSOR_PARAMETERS_NOT_SANE = 0x9,
	ENUM_FORCE_INT32_SIZE(ElmigParamsFitterRetCode)
};

ECHMET_ST_ENUM(FixedParameterType) {
	FPT_MOBILITY = 0,
	FPT_PKA = 0x1,
	ENUM_FORCE_INT32_SIZE(ElmigParamsFitterFixedParameterType)
};

class InBuffer {
public:
	SysComp::InConstituentVec *composition;
	RealVec *concentrations;
	double uEffExperimental;
};
IS_POD(InBuffer)

typedef Vec<InBuffer> InBufferVec;

class InSystem {
public:
	InBufferVec *buffers;
	SysComp::InConstituent analyte;
	NonidealityCorrections corrections;
};
IS_POD(InSystem)

typedef Vec<bool> InIsFixedVec;

class ExpectedCurvePoint {
public:
	double pH;
	double experimental;
	double expected;
};
IS_POD(ExpectedCurvePoint)

typedef Vec<ExpectedCurvePoint> ExpectedCurvePointVec;

class FittedParameter {
public:
	int charge;
	double value;
	double stDev;
};
IS_POD(FittedParameter)

typedef Vec<FittedParameter> FittedParameterVec;

class FitResults {
public:
	FittedParameterVec *mobilities;
	FittedParameterVec *pKas;
};
IS_POD(FitResults)

class ParametersFixer {
public:
	virtual RetCode ECHMET_CC add(const FixedParameterType type, const int id, const double value) ECHMET_NOEXCEPT = 0;
	virtual void ECHMET_CC destroy() ECHMET_NOEXCEPT = 0;
	virtual bool ECHMET_CC isFixed(const FixedParameterType type, const int id) const ECHMET_NOEXCEPT = 0;
	virtual void ECHMET_CC remove(const FixedParameterType type, const int id) ECHMET_NOEXCEPT = 0;
	virtual double ECHMET_CC value(const FixedParameterType type, const int id) const ECHMET_NOEXCEPT = 0;

protected:
	virtual ECHMET_CC ~ParametersFixer() ECHMET_NOEXCEPT = 0;
};

extern "C" {

ECHMET_API bool checkSanity(const SysComp::InConstituent &analyte) ECHMET_NOEXCEPT;
ECHMET_API InBufferVec * ECHMET_CC createInBufferVec() ECHMET_NOEXCEPT;
ECHMET_API ParametersFixer * ECHMET_CC createParametersFixer() ECHMET_NOEXCEPT;
ECHMET_API const char * ECHMET_CC EMPFerrorToString(const RetCode tRet) ECHMET_NOEXCEPT;
ECHMET_API RetCode ECHMET_CC expectedCurve(const InSystem &system, const FitResults &results, ExpectedCurvePointVec *&curve) ECHMET_NOEXCEPT;
ECHMET_API double ECHMET_CC mobilityLowerBound() ECHMET_NOEXCEPT;
ECHMET_API double ECHMET_CC mobilityUpperBound() ECHMET_NOEXCEPT;
ECHMET_API RetCode ECHMET_CC process(const InSystem &system, const ParametersFixer *fixer, FitResults &results) ECHMET_NOEXCEPT;
ECHMET_API void ECHMET_CC releaseInBuffer(const InBuffer &buffer) ECHMET_NOEXCEPT;
ECHMET_API void ECHMET_CC releaseInBufferVec(const InBufferVec *vec) ECHMET_NOEXCEPT;
ECHMET_API void ECHMET_CC releaseInSystem(InSystem &system) ECHMET_NOEXCEPT;
ECHMET_API void ECHMET_CC releaseResults(FitResults &results) ECHMET_NOEXCEPT;

}

} // namespace ElmigParamsFitter
} // namespace ECHMET

#endif // ECHMET_ELMIG_PARAMS_FITTER_H
