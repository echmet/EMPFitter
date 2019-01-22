#include "mobdissocregressor.h"

#include "totalequilibrium.h"
#include "types.h"

#include <echmetphchconsts.h>
#include <echmetcaes.h>
#include <echmetionprops.h>
#include <echmetsyscomp.h>

namespace ECHMET {

#ifdef ECHMET_MOBDISSOCREGRESSOR_DEBUG_DUMP

#include <iostream>

inline
void dumpAppliedParameters(const SysComp::InConstituent &ctuent)
{
	std::cout << "Current params\n";
	std::cout << "Mobilities:\t";
	for (size_t idx = 0; idx < ctuent.mobilities->size(); idx++)
		std::cout << ctuent.mobilities->at(idx) << " ";

	std::cout << "\npKas:\t";
	for (size_t idx = 0; idx < ctuent.pKas->size(); idx++)
		std::cout << ctuent.pKas->at(idx) << " ";
	std::cout << "\n";
}

inline
void debugDumpEquilibrium(const SysComp::ChemicalSystem &chemSystem, const SysComp::CalculatedProperties &calcProps)
{
	std::cout << "--- Composition ---\n";
	std::cout << "pH: " << IonProps::calculatepH_direct(calcProps.ionicConcentrations->at(0), calcProps.ionicStrength) << "\n";
	std::cout << "IS: " << calcProps.ionicStrength * 1000.0 << " (mM)\n";
	for (size_t idx = 2; idx < chemSystem.ionicForms->size(); idx++) {
		const auto iF = chemSystem.ionicForms->at(idx);
		std::cout << iF->name->c_str() << " " << calcProps.ionicConcentrations->at(iF->ionicConcentrationIndex) << "(mM)\n";
	}
	std::cout << "\n";
}

inline
void debugDumpMobilities(const SysComp::ChemicalSystem &chemSystem, const SysComp::CalculatedProperties &calcProps)
{
	std::cout << "--- Mobilities ---\n";
	for (size_t idx = 2; idx < chemSystem.ionicForms->size(); idx++) {
		const auto iF = chemSystem.ionicForms->at(idx);
		std::cout << iF->name->c_str() << " " << calcProps.ionicMobilities->at(iF->ionicMobilityIndex) << "\n";
	}
	std::cout << "\n";
}

inline
void debugDumpUEff(const double uEff, const double cH, const double ionicStrength)
{
	std::cout << "uEff at pH " << IonProps::calculatepH_direct(cH, ionicStrength) << ": " << uEff << "\n";
}

inline
void debugDumpInput(const size_t NMobs, const size_t NpKas, const size_t NParams)
{
	std::cout << "NMobilities: " << NMobs << ", NpKas: " << NpKas << ", NParams: " << NParams << "\n";
}
#else
inline
void dumpAppliedParameters(const SysComp::InConstituent &)
{
}

inline
void debugDumpEquilibrium(const SysComp::ChemicalSystem &, const SysComp::CalculatedProperties &)
{
}

inline
void debugDumpMobilities(const SysComp::ChemicalSystem &, const SysComp::CalculatedProperties &)
{
}

inline
void debugDumpUEff(const double, const double, const double)
{
}

inline
void debugDumpInput(const size_t, const size_t, const size_t)
{
}
#endif // ECHMET_MOBDISSOCREGRESSOR_DEBUG_DUMP


inline
std::vector<double> MakeActivityCoefficients(const SysComp::InConstituent &analyte, const double ionicStrength)
{
	const int max = std::max(std::abs(analyte.chargeLow), std::abs(analyte.chargeHigh));

	std::vector<double> actCoeffs{1.0};
	for (int ch = 1; ch <= max; ch++)
		actCoeffs.emplace_back(PhChConsts::activityCoefficient(ionicStrength, ch));

	return actCoeffs;
}

inline
CAES::TotalEquilibrium<double, true> MakeTotalEquilibrium(const SysComp::InConstituent &analyte)
{
	const auto pBs = [analyte]() {
		std::vector<double> vec{};

		for (size_t idx = 0; idx < analyte.pKas->size(); idx++)
			vec.emplace_back(analyte.pKas->at(idx));

		return vec;
	}();
	return { analyte.chargeLow, analyte.chargeHigh, pBs };
}

MobDissocRegressor::MobDissocRegressor(const BufferSystemVec &bufSysVec, const InConstituentWrapper &analyte,
				       const NonidealityCorrections corrections) :
	RFType(CountNumberOfParams(analyte())),
	m_bufSysVec(bufSysVec),
	m_analyte(analyte),
	m_solverCorrections(corrections),
	m_NMobilities(CountNumberOfMobilities(analyte())),
	m_NParams(CountNumberOfParams(analyte()))
{
	debugDumpInput(m_NMobilities, analyte().pKas->size(), m_NParams);
}

void MobDissocRegressor::AAssign(const RFType &otherBase)
{
	const auto &other = dynamic_cast<const MobDissocRegressor &>(otherBase);

	const_cast<size_t&>(m_NMobilities) = other.m_NMobilities;
	const_cast<size_t&>(m_NParams) = other.m_NParams;
	m_bufSysVec = other.m_bufSysVec;
	m_analyte = other.m_analyte;
	m_solverCorrections = other.m_solverCorrections;
}

double MobDissocRegressor::ACalculateDerivative(const BufferSystem &x, const ParamsVector &params, const typename RFType::index_type param_idx, msize_t idx) const
{
	(void)idx;

	static const double H{1.0e-11};

	auto shiftedParams = params;
	/* Low */
	shiftedParams[param_idx] -= H;
	const double low = CalculateEffectiveMobility(x, shiftedParams);

	/* High */
	shiftedParams[param_idx] += 2.0 * H;
	const double high = CalculateEffectiveMobility(x, shiftedParams);

	const double der = (high - low) / (2.0 * H);

	return der;
}

double MobDissocRegressor::ACalculateFx(const BufferSystem &x, const ParamsVector &params, msize_t idx) const
{
	(void)idx;

	return CalculateEffectiveMobility(x, params);
}

bool MobDissocRegressor::ACheckSanity(const ParamsVector &params) const
{
	return CheckSanityInternal(params, m_analyte(), m_NMobilities);
}

MobDissocRegressor * MobDissocRegressor::ACreate() const
{
	/* Calling this with no parameters makes no sense. Can this break stuff? */
	return new MobDissocRegressor(m_bufSysVec, m_analyte, m_solverCorrections);
}

bool MobDissocRegressor::ASetInitialParameters(ParamsVector &params, const XTVector &x, const YTVector &y)
{
	(void)x; (void)y;

	SetParameters(params, m_analyte());

	return this->ACheckSanity(params);
}

bool MobDissocRegressor::CheckAnalyteSanity(const SysComp::InConstituent &analyte)
{
	const auto NParams = CountNumberOfParams(analyte);

	if (analyte.chargeLow == 0 && analyte.chargeHigh == 0)
		return false;

	ParamsVector params{};
	params.resize(NParams);

	SetParametersRaw(params, analyte, NParams);

	return CheckSanityInternal(params, analyte, CountNumberOfMobilities(analyte));
}

bool MobDissocRegressor::CheckSanity()
{
	if (!CheckAnalyteSanity(m_analyte()))
		return false;

	return GetPCount() > 0;
}

bool MobDissocRegressor::Initialize(const std::vector<double> &inYVec)
{
	static_assert(std::numeric_limits<YTVector::Index>::max() <=
		      std::numeric_limits<std::vector<double>::size_type>::max(),
		      "Incompatible vector size types");

	YTVector yVec = YTVector(inYVec.size());

	/* Copy-ya-hey */
	for (size_t idx = 0; idx < inYVec.size(); idx++)
		yVec(idx) = inYVec.at(idx);

	/* WARNING: Derived type must make a manual call to the parent's
	 *          Initialize() function. This is pretty inconvenient
	 *          but it would required a change of the base regressor
	 *          initialization logic to get around this.
	 */
	return RFType::Initialize(m_bufSysVec, yVec, 1.0e-8, 200, true);
}

std::pair<double, double> MobDissocRegressor::SolveBuffer(const SysComp::InConstituentVec *composition,
							  const SysComp::InConstituent &analyte,
							  const RealVec *concentrations,
							  const std::pair<bool, double> cH,
							  const NonidealityCorrections corrs)
{
	SysComp::ChemicalSystem chemSystem{};
	SysComp::CalculatedProperties calcProps{};

	auto compositionFull = composition->duplicate();
	compositionFull->push_back(analyte);
	auto concVec = concentrations->duplicate();
	concVec->push_back(0);

	assert(compositionFull->size() == concVec->size());

	auto tRet = SysComp::makeComposition(chemSystem, calcProps, compositionFull);
	if (tRet != RetCode::OK) {
		concVec->destroy();
		compositionFull->destroy();

		throw RegressCore::RegressException("Unable to create system composition: " + std::string(errorToString(tRet)));
	}

	CAES::SolverContext *ctx{};
	CAES::Solver *solver{};

	tRet = CAES::createSolverContext(ctx, chemSystem);
	if (tRet != RetCode::OK) {
		concVec->destroy();
		compositionFull->destroy();
		SysComp::releaseChemicalSystem(chemSystem);
		SysComp::releaseCalculatedProperties(calcProps);

		throw RegressCore::RegressException("Unable to create solver context");
	}

	solver = CAES::createSolver(ctx, CAES::Solver::defaultOptions(), corrs);
	if (solver == nullptr) {
		concVec->destroy();
		ctx->destroy();
		compositionFull->destroy();
		SysComp::releaseChemicalSystem(chemSystem);
		SysComp::releaseCalculatedProperties(calcProps);

		throw RegressCore::RegressException("Unable to create solver");
	}

	if (cH.first) {
		tRet = solver->estimateDistributionFast(cH.second, concVec, calcProps);
		if (tRet != RetCode::OK)
			tRet = solver->estimateDistributionSafe(concVec, calcProps);
	} else
		tRet = solver->estimateDistributionSafe(concVec, calcProps);

	if (tRet != RetCode::OK) {
		concVec->destroy();
		solver->destroy();
		ctx->destroy();
		compositionFull->destroy();
		SysComp::releaseChemicalSystem(chemSystem);
		SysComp::releaseCalculatedProperties(calcProps);

		throw RegressCore::RegressException("Failed to estimate distribution");
	}

	debugDumpEquilibrium(chemSystem, calcProps);

	/* TODO: Figure out how to handle zero concentration of the analyte */
	/*tRet = solver->solve(concVec, calcProps, 1000, nullptr);
	if (tRet != RetCode::OK) {
		concVec->destroy();
		ctx->destroy();
		compositionFull->destroy();

		throw RegressException("Failed to solve equilibrium: " + std::string(errorToString(tRet)));
	}

	debugDumpEquilibrium(chemSystem, calcProps);*/

	auto ionCtx = IonProps::makeComputationContext(chemSystem, IonProps::ComputationContext::defaultOptions());
	if (ionCtx == nullptr) {
		concVec->destroy();
		solver->destroy();
		ctx->destroy();
		compositionFull->destroy();
		SysComp::releaseChemicalSystem(chemSystem);
		SysComp::releaseCalculatedProperties(calcProps);

		throw RegressCore::RegressException("Failed to create IonProps context");
	}

	tRet = IonProps::correctMobilities(ionCtx, corrs, concVec, calcProps);
	if (tRet != RetCode::OK) {
		concVec->destroy();
		ionCtx->destroy();
		solver->destroy();
		ctx->destroy();
		compositionFull->destroy();
		SysComp::releaseChemicalSystem(chemSystem);
		SysComp::releaseCalculatedProperties(calcProps);

		throw RegressCore::RegressException("Failed to correct ionic mobilities");
	}

	debugDumpMobilities(chemSystem, calcProps);

	auto te = MakeTotalEquilibrium(analyte);
	auto actCoeffs = MakeActivityCoefficients(analyte, calcProps.ionicStrength);
	auto dist = te.distribution(calcProps.ionicConcentrations->at(0), actCoeffs);
	assert(dist.size() == analyte.mobilities->size());

	/* Now calculate the actual mobility */
	double uEff{0.0};
	{
		const auto sgn = [](const double v) {
			return (v > 0) - (v < 0);
		};

		const auto &csAnalyte = chemSystem.constituents->back(); /* Our analyte is always in the back */

		for (size_t idx = 0; idx < csAnalyte->ionicForms->size(); idx++) {
			const auto iF = csAnalyte->ionicForms->at(idx);
			const double mob = calcProps.ionicMobilities->at(iF->ionicMobilityIndex);

			uEff += sgn(iF->totalCharge) * mob * dist.at(idx);
		}
	}
	debugDumpUEff(uEff, calcProps.ionicConcentrations->at(0), calcProps.ionicStrength);

	const double pH = [&]() {
		if (cH.first)
			return cH.second;
		return IonProps::calculatepH_direct(calcProps.ionicConcentrations->at(0),
						    calcProps.ionicStrength);
	}();

	concVec->destroy();
	ionCtx->destroy();
	solver->destroy();
	ctx->destroy();
	compositionFull->destroy();
	SysComp::releaseChemicalSystem(chemSystem);
	SysComp::releaseCalculatedProperties(calcProps);

	return { uEff, pH };
}

void MobDissocRegressor::ApplyParameters(SysComp::InConstituent &analyte, const ParamsVector &params) const
{
	size_t idx{0};
	PVSize param_idx{0};
	for (int charge = m_analyte().chargeLow; charge <= m_analyte().chargeHigh; charge++, idx++) {
		if (charge == 0)
			continue;

		analyte.mobilities->operator[](idx) = params.at(param_idx);
		param_idx++;
	}

	for (idx = 0; param_idx < m_NParams; idx++, param_idx++)
		analyte.pKas->operator[](idx) = params.at(param_idx);

	dumpAppliedParameters(analyte);
}

double MobDissocRegressor::CalculateEffectiveMobility(const BufferSystem &x, const ParamsVector &params) const
{
	InConstituentWrapper analyte = m_analyte;

	ApplyParameters(analyte(), params);

	return SolveBuffer(x.composition(), analyte(), x.concentrationsRVec(),
			   { true, x.cH() }, m_solverCorrections).first;
}

bool MobDissocRegressor::CheckSanityInternal(const ParamsVector &params, const SysComp::InConstituent &analyte,
					     const size_t NMobilities) noexcept
{
	static const auto isMobilitiesReasonable = [](const double uPrev, const double uCurr,
						      const double uBase) {
		const double lowerBound = uPrev + uBase * MOBILITY_LOWER_BOUND;
		const double upperBound = uPrev + uBase * MOBILITY_UPPER_BOUND;
		return (uCurr > lowerBound) && (uCurr < upperBound);
	};

	const auto mobilityIndex = [&analyte](int charge) {
		const PVSize uZeroIdx = static_cast<PVSize>(-analyte.chargeLow);	/* "Imaginary" index of mobility of the zero charge parameter */

		if (charge > 0)
			charge--;
		return uZeroIdx + charge;
	};

	PVSize idx{0};

	/* Elementary sanity checks.
	 * 1) All mobilities except for the mobility of zero charge
	 * must be positive. */
	for (int charge = analyte.chargeLow; charge <= analyte.chargeHigh; charge++) {
		if (charge == 0)
			continue;
		else if (charge != 0 && params.at(idx) <= 0.0)
			return false;

		idx++;
	}

	/* 2) pKa values must be descending */
	idx = NMobilities + 1;
	for (size_t jdx = 1; jdx < analyte.pKas->size(); idx++, jdx++) {
		double pKaPrev = params.at(idx - 1);
		double pKaCurr = params.at(idx);

		if (pKaCurr >= pKaPrev)
			return false;
	}

	/* More advanced sanity checks.
	 * 1) pKas values must be reasonable.
	 */
	idx = NMobilities;
	for (size_t jdx = 0; jdx < analyte.pKas->size(); idx++, jdx++) {
		const double pKa = params.at(idx);

		if (pKa < 0.0 || pKa > 14.0)
			return false;

	}

	/* 2)
	 * Each additional charge can increase the overall mobility
	 * only by the the amount similar to mobility of the analyte
	 * with +/-1 charge. Margin of flexibility is defined by
	 * <tt>MOBILITY_LOWER_BOUND</tt> and <tt>MOBILITY_UPPER_BOUND</tt>
	 * constants.
	 */
	int fromCharge = analyte.chargeHigh > -2 ? -2 : analyte.chargeHigh - 1;
	idx = mobilityIndex(fromCharge);
	auto baseIdx = mobilityIndex(fromCharge + 1);
	for (int charge = fromCharge; charge >= analyte.chargeLow; charge--, idx--) {
		const double uBase = params.at(baseIdx) / std::abs(fromCharge + 1);
		double uPrev = params.at(idx + 1);
		double uCurr = params.at(idx);


		if (!isMobilitiesReasonable(uPrev, uCurr, uBase))
			return false;
	}

	fromCharge = analyte.chargeLow < 2 ? 2 : analyte.chargeLow + 1;
	idx = mobilityIndex(fromCharge);
	baseIdx = mobilityIndex(fromCharge - 1);
	for (int charge = fromCharge; charge <= analyte.chargeHigh; charge++, idx++) {
		const double uBase = params.at(baseIdx) / std::abs(fromCharge - 1);
		double uPrev = params.at(idx - 1);
		double uCurr = params.at(idx);

		if (!isMobilitiesReasonable(uPrev, uCurr, uBase))
			return false;
	}

	return true;

}

MobDissocRegressor::PVSize MobDissocRegressor::CountNumberOfMobilities(const SysComp::InConstituent &analyte) noexcept
{
	auto sgn = [](const auto v) {
		using T = decltype(v);
		return (v > T{0}) - (v < T{0});
	};

	PVSize count = analyte.chargeHigh - analyte.chargeLow;
	if (sgn(analyte.chargeHigh) == sgn(analyte.chargeLow))
		count++;

	return count;
}

MobDissocRegressor::PVSize MobDissocRegressor::CountNumberOfParams(const SysComp::InConstituent &analyte) noexcept
{
	return CountNumberOfMobilities(analyte) + analyte.pKas->size();
}

void MobDissocRegressor::SetParameters(ParamsVector &params, const SysComp::InConstituent &analyte)
{
	assert(params.size() == m_NParams);

	/*
	 * Set vector of parameters.
	 * Mobilities come first in the vector of parameters
	 * followed by pKa constants.
	 * Data order:
	 * u(lowCharge) -> u(highCharge) -> pKa(lowCharge) -> pKa(highCharge), mobility of zero charge is excluded
		 */
	size_t idx{0};
	PVSize param_idx{0};
	for (int charge = analyte.chargeLow; charge <= analyte.chargeHigh; charge++, idx++) {
		if (charge == 0)
			continue;

		this->SetParam(params, param_idx, analyte.mobilities->at(idx));
		param_idx++;
	}

	for (idx = 0; param_idx < m_NParams; idx++, param_idx++)
		this->SetParam(params, param_idx, analyte.pKas->at(idx));
}

void MobDissocRegressor::SetParametersRaw(ParamsVector &params, const SysComp::InConstituent &analyte,
					  const PVSize NParams)
{
	assert(params.size() == NParams);

	/*
	 * Set vector of parameters.
	 * Mobilities come first in the vector of parameters
	 * followed by pKa constants.
	 * Data order:
	 * u(lowCharge) -> u(highCharge) -> pKa(lowCharge) -> pKa(highCharge), mobility of zero charge is excluded
	 */
	size_t idx{0};
	PVSize param_idx{0};
	for (int charge = analyte.chargeLow; charge <= analyte.chargeHigh; charge++, idx++) {
		if (charge == 0)
			continue;

		params[param_idx] = analyte.mobilities->at(idx);
		param_idx++;
	}

	for (idx = 0; param_idx < NParams; idx++, param_idx++)
		params[param_idx] = analyte.pKas->at(idx);
}

} // namespace ECHMET
