#include <echmetelmigparamsfitter.h>

#include "mobdissocregressor.h"
#include "parametersfixerimpl.h"

#include <cassert>
#include <fstream>

#define USE_ECHMET_CONTAINERS

#define ECHMET_IMPORT_INTERNAL
#include <containers/echmetvec_p.h>
#undef ECHMET_IMPORT_INTERNAL

#define M_STRINGIFY(input) #input
#define ERROR_CODE_CASE(erCase) case RetCode::erCase: return M_STRINGIFY(erCase)

namespace ECHMET {
namespace ElmigParamsFitter {

class TraceDumper {
public:
	~TraceDumper()
	{
		const auto log = TRACER_INSTANCE<RegressTracing>().logged();

		std::ofstream ofs("dump.txt");
		ofs << log;
	}
};

static
std::tuple<double> statDataFromVariance(const double variance) noexcept
{
	const double stDev = std::sqrt(variance);

	return { stDev };
}

static
std::vector<double> realVecToSTL(const RealVec *vec)
{
	std::vector<double> stlVec(vec->size());

	for (size_t idx = 0; idx < vec->size(); idx++)
		stlVec[idx] = vec->at(idx);

	return stlVec;
}

static
RetCode fillResults(const MobDissocRegressor::ParamsVector &params, const MobDissocRegressor::YTVector &variances, const SysComp::InConstituent &analyte, FitResults &results)
{
	const int chargeSpan = analyte.chargeHigh - analyte.chargeLow;
	FittedParameterVec *mobilities = nullptr;
	FittedParameterVec *pKas = nullptr;
	size_t idx = 0;
	int pKaNum = analyte.chargeLow;

	mobilities = createECHMETVec<FittedParameter, false>(chargeSpan + 1);
	pKas = createECHMETVec<FittedParameter, false>(chargeSpan);

	if (!(mobilities && pKas))
		goto err_out;

	for (int charge = analyte.chargeLow; charge <= analyte.chargeHigh; charge++) {
		if (charge == 0)
			continue;

		const auto statData = statDataFromVariance(variances(idx));

		FittedParameter p{charge, params.at(idx),
				  std::get<0>(statData)};

		if (mobilities->push_back(p) != ECHMET::RetCode::OK)
			goto err_out;
		idx++;
	}

	for (size_t jdx = 0; jdx < analyte.pKas->size(); jdx++, idx++, pKaNum++) {
		if (pKaNum == 0)
			pKaNum++;

		const auto statData = statDataFromVariance(variances(idx));

		FittedParameter p{pKaNum, params.at(idx),
				  std::get<0>(statData)};

		if (pKas->push_back(p) != ECHMET::RetCode::OK)
			goto err_out;
	}

	results.mobilities = mobilities;
	results.pKas = pKas;

	return RetCode::OK;

err_out:
	if (mobilities)
		mobilities->destroy();
	if (pKas)
		pKas->destroy();

	return RetCode::E_NO_MEMORY;
}

static
void fixRegressParameters(MobDissocRegressor &regressor, const ParametersFixer *fixer,
			  const SysComp::InConstituent &analyte) noexcept
{
	const auto pcnt = regressor.GetPTotal();
	MobDissocRegressor::index_type pidx = 0;

	for (int charge = analyte.chargeLow; charge <= analyte.chargeHigh; charge++) {
		if (charge == 0)
			continue;

		assert(pidx < pcnt);

		if (fixer->isFixed(FixedParameterType::FPT_MOBILITY, charge))
			regressor.FixParameter(pidx,
					       fixer->value(FixedParameterType::FPT_MOBILITY, charge));
		pidx++;
	}

	int charge = analyte.chargeLow;
	while (pidx < pcnt) {
		if (charge == 0) {
			charge++;
			continue;
		}

		if (fixer->isFixed(FixedParameterType::FPT_PKA, charge))
			regressor.FixParameter(pidx,
					       fixer->value(FixedParameterType::FPT_PKA, charge));

		pidx++;
		charge++;
	}
}

static
RetCode fit(BufferSystemVec bufSysVec, std::vector<double> expDataVec,
	    const SysComp::InConstituent &analyte,
	    const NonidealityCorrections corrs,
	    const ParametersFixer *fixer,
	    FitResults &results)
{
	TraceDumper dumper{};

	TRACER_INSTANCE<RegressTracing>().enableAllTracepoints();

	MobDissocRegressor regressor(std::move(bufSysVec), InConstituentWrapper(analyte), corrs);

	if (fixer != nullptr)
		fixRegressParameters(regressor, fixer, analyte);

	try {
		if (!regressor.CheckSanity())
			return RetCode::E_REGRESSOR_PARAMETERS_NOT_SANE;
		if (!regressor.Initialize(expDataVec))
			return RetCode::E_REGRESSOR_INITIALIZATION;
	} catch (const RegressCore::RegressException &ex) {
		return RetCode::E_REGRESSOR_INTERNAL_ERROR;
	}

	try {
		if (!regressor.Regress())
			return RetCode::E_REGRESSOR_NO_SOLUTION;
	} catch (const RegressCore::RegressException &) {
		return RetCode::E_REGRESSOR_INTERNAL_ERROR;
	}

	return fillResults(regressor.GetParameters(), regressor.GetVariances(), analyte, results);
}

static
std::tuple<BufferSystemVec,
	   std::vector<double>> makeBuffersVector(const InBufferVec *buffers, const NonidealityCorrections corrs)
{
	BufferSystemVec bVec{};
	std::vector<double> eVec{};

	bVec.reserve(buffers->size());
	eVec.reserve(buffers->size());

	for (size_t idx = 0; idx < buffers->size(); idx++) {
		const auto &buf = buffers->at(idx);

		std::vector<double> concs = realVecToSTL(buf.concentrations);

		bVec.emplace_back(buf.composition, std::move(concs), nonidealityCorrectionIsSet(corrs, NonidealityCorrectionsItems::CORR_DEBYE_HUCKEL));
		eVec.emplace_back(buf.uEffExperimental);
	}

	return { bVec, eVec };
}

bool ECHMET_CC checkSanity(const SysComp::InConstituent &analyte) noexcept
{
	return MobDissocRegressor::CheckAnalyteSanity(analyte);
}

InBufferVec * ECHMET_CC createInBufferVec() noexcept
{
	return ECHMET::createECHMETVec<InBuffer, false>(0);
}

ParametersFixer * ECHMET_CC createParametersFixer() noexcept
{
	return new (std::nothrow) ParametersFixerImpl();
}

const char * ECHMET_CC EMPFerrorToString(const RetCode tRet) noexcept
{
	switch (tRet) {
	ERROR_CODE_CASE(OK);
	ERROR_CODE_CASE(E_INVALID_ARGUMENT);
	ERROR_CODE_CASE(E_NOT_ENOUGH_MEASUREMENTS);
	ERROR_CODE_CASE(E_BUFFERS_EXPDATA_MISMATCH);
	ERROR_CODE_CASE(E_INVALID_BUFFER);
	ERROR_CODE_CASE(E_NO_MEMORY);
	ERROR_CODE_CASE(E_REGRESSOR_INITIALIZATION);
	ERROR_CODE_CASE(E_REGRESSOR_NO_SOLUTION);
	ERROR_CODE_CASE(E_REGRESSOR_INTERNAL_ERROR);
	ERROR_CODE_CASE(E_REGRESSOR_PARAMETERS_NOT_SANE);
	default:
		return "Unknown error";
	}
}

RetCode ECHMET_CC expectedCurve(const InSystem &system, const FitResults &results, ExpectedCurvePointVec *&curve) noexcept
{
	auto curveWrap = std::unique_ptr<ExpectedCurvePointVec,
					 std::function<void (ExpectedCurvePointVec *)>>(ECHMET::createECHMETVec<ExpectedCurvePoint, false>(0),
											[](ExpectedCurvePointVec *vec) { vec->destroy(); });
	if (curveWrap == nullptr)
		return RetCode::E_NO_MEMORY;

	InConstituentWrapper fixedAnalyte{system.analyte};

	for (size_t idx = 0; idx < results.pKas->size(); idx++)
		fixedAnalyte().pKas->operator[](idx) = results.pKas->at(idx).value;
	{
		size_t idx = 0, jdx = 0;
		auto charge = system.analyte.chargeLow;
		while (idx < results.mobilities->size()) {
			if (charge == 0)
				jdx++;

			fixedAnalyte().mobilities->operator[](jdx) = results.mobilities->at(idx).value;

			charge++;
			idx++;
			jdx++;
		}
	}

	for (size_t idx = 0; idx < system.buffers->size(); idx++) {
		const auto &buffer = system.buffers->at(idx);

		try {
			const auto sol = MobDissocRegressor::SolveBuffer(buffer.composition, fixedAnalyte(),
									 buffer.concentrations, { false, 0.0 },
									 system.corrections);
			auto tRet = curveWrap->push_back({ sol.second, buffer.uEffExperimental, sol.first });
			if (tRet != ECHMET::RetCode::OK)
				return RetCode::E_NO_MEMORY;
		} catch (const RegressCore::RegressException &ex) {
			return RetCode::E_REGRESSOR_INTERNAL_ERROR;
		}

	}

	curve = curveWrap.release();

	return RetCode::OK;
}

double ECHMET_CC mobilityLowerBound() noexcept
{
	return MobDissocRegressor::MOBILITY_LOWER_BOUND;
}

double ECHMET_CC mobilityUpperBound() noexcept
{
	return MobDissocRegressor::MOBILITY_UPPER_BOUND;
}

RetCode ECHMET_CC process(const InSystem &system, const ParametersFixer *fixer, FitResults &results) noexcept
{
	if (system.analyte.mobilities->size() < 1 ||
	    (system.analyte.pKas->size() != system.analyte.mobilities->size() - 1))
		return RetCode::E_INVALID_ARGUMENT;

	if (system.buffers->size() <
	    (system.analyte.mobilities->size() + system.analyte.pKas->size()))
		return RetCode::E_NOT_ENOUGH_MEASUREMENTS;

	try {
		auto dataPack = makeBuffersVector(system.buffers, system.corrections);

		return fit(std::move(std::get<0>(dataPack)),
			   std::move(std::get<1>(dataPack)),
			   system.analyte, system.corrections, fixer, results);
	} catch (const std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	} catch (const BufferSystem::BufferSystemException &) {
		return RetCode::E_INVALID_BUFFER;
	}
}

void ECHMET_CC releaseInBuffer(const InBuffer &buffer) noexcept
{
	SysComp::releaseInputData(buffer.composition);

	if (buffer.concentrations != nullptr)
		buffer.concentrations->destroy();
}

void ECHMET_CC releaseInBufferVec(const InBufferVec *vec) noexcept
{
	if (vec == nullptr)
		return;

	for (size_t idx = 0; idx < vec->size(); idx++) {
		const auto &buf = vec->at(idx);

		releaseInBuffer(buf);
	}

	vec->destroy();
}

void ECHMET_CC releaseInSystem(InSystem &system) noexcept
{
	releaseInBufferVec(system.buffers);
	SysComp::releaseInConstituent(system.analyte);
}

void ECHMET_CC releaseResults(FitResults &results) noexcept
{
	if (results.mobilities != nullptr) {
		results.mobilities->destroy();
		results.mobilities = nullptr;
	}

	if (results.pKas != nullptr) {
		results.pKas->destroy();
		results.pKas = nullptr;
	}
}

} // namespace ElmigParamsFitter
} // namespace ECHMET
