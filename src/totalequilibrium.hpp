#ifndef ECHMET_CAES_TOTALEQUILIBRIUM_HPP
#define ECHMET_CAES_TOTALEQUILIBRIUM_HPP

namespace ECHMET {
namespace CAES {

template <typename CAESReal>
static ECHMET_FORCE_INLINE
void calculateTs(std::vector<CAESReal> &ts, const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X, const std::vector<CAESReal> &Ls, const int numLow)
{
	const size_t len = ts.size();
	X = 1.0;
	ts[0] = 1.0;

	const CAESReal activityLow = activityCoefficients[std::abs(numLow)];
	const CAESReal vMul = v * activityCoefficients[1];
	int num = numLow + 1;
	CAESReal vPow = vMul;
	for (size_t idx = 1; idx < len; idx++) {
		const CAESReal T = Ls[idx] * vPow * (activityLow / activityCoefficients[std::abs(num++)]);

		vPow *= vMul;

		ts[idx] = T;
		X += T;
	}
}

template <typename CAESReal>
static ECHMET_FORCE_INLINE
void calculatedTsdV(std::vector<CAESReal> &dts, const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X, const std::vector<CAESReal> &Ls, const int numLow)
{
	const size_t len = dts.size();

	X = 0.0;
	dts[0] = 0.0;

	const CAESReal activityLow = activityCoefficients[std::abs(numLow)];
	const CAESReal vMul = v * activityCoefficients[1];
	int num = numLow + 1;
	CAESReal vPow = 1.0;
	for (size_t idx = 1; idx < len; idx++) {
		const CAESReal T = idx * Ls[idx] * vPow * (activityLow / activityCoefficients[std::abs(num++)]);

		vPow *= vMul;

		dts[idx] = T;
		X += T;
	}
}

/*!
 * TotalEquilibrium c-tor.
 *
 * @param[in] numLow Lowest equilibrium index.
 * @param[in] numHigh Highest equilibirium index.
 * @param[in] pBs Vector of consecutive equilibrium constants.
 * @param[in] concentration Analytical concentration of the constituent.
 */
template <typename CAESReal>
TotalEquilibrium<CAESReal, true>::TotalEquilibrium(const int numLow, const int numHigh, const std::vector<CAESReal> &pBs, const size_t concentrationIndex) :
	concentrationIndex(concentrationIndex),
	Ls(calculateLs(pBs)),
	numLow(numLow),
	numHigh(numHigh)
{
}

/*!
 * TotalEquilibrium c-tor.
 *
 * @param[in] numLow Lowest equilibrium index.
 * @param[in] numHigh Highest equilibirium index.
 * @param[in] pBs Vector of consecutive equilibrium constants.
 * @param[in] concentration Analytical concentration of the constituent.
 */
template <typename CAESReal>
TotalEquilibrium<CAESReal, false>::TotalEquilibrium(const int numLow, const int numHigh, const std::vector<CAESReal> &pBs, const size_t concentrationIndex) :
	concentrationIndex(concentrationIndex),
	Ls(calculateLs(pBs)),
	numLow(numLow),
	numHigh(numHigh)
{
	m_concentrations.resize(Ls.size());
	m_dDistdV.resize(Ls.size());
	m_distribution.resize(Ls.size());
	m_dTsdV.resize(Ls.size());
	m_Ts.resize(Ls.size());
}

template <typename CAESReal>
TotalEquilibrium<CAESReal, true>::TotalEquilibrium(const TotalEquilibrium &other) :
	concentrationIndex(other.concentrationIndex),
	Ls(std::move(other.Ls)),
	numLow(other.numLow),
	numHigh(other.numHigh)
{
}

template <typename CAESReal>
TotalEquilibrium<CAESReal, false>::TotalEquilibrium(const TotalEquilibrium &other) :
	concentrationIndex(other.concentrationIndex),
	Ls(std::move(other.Ls)),
	numLow(other.numLow),
	numHigh(other.numHigh),
	m_concentrations(other.m_concentrations),
	m_dDistdV(other.m_dDistdV),
	m_distribution(other.m_distribution),
	m_dTsdV(other.m_dTsdV),
	m_Ts(other.m_Ts)
{
}


/*!
 * TotalEquilibrium move c-tor.
 *
 * @param[in] other TotalEquilibrium object being moved.
 */
template <typename CAESReal>
TotalEquilibrium<CAESReal, true>::TotalEquilibrium(TotalEquilibrium &&other) :
	concentrationIndex(other.concentrationIndex),
	Ls(std::move(other.Ls)),
	numLow(other.numLow),
	numHigh(other.numHigh)
{
}

/*!
 * TotalEquilibrium move c-tor.
 *
 * @param[in] other TotalEquilibrium object being moved.
 */
template <typename CAESReal>
TotalEquilibrium<CAESReal, false>::TotalEquilibrium(TotalEquilibrium &&other) :
	concentrationIndex(other.concentrationIndex),
	Ls(std::move(other.Ls)),
	numLow(other.numLow),
	numHigh(other.numHigh),
	m_concentrations(std::move(other.m_concentrations)),
	m_dDistdV(std::move(other.m_dDistdV)),
	m_distribution(std::move(other.m_distribution)),
	m_dTsdV(std::move(other.m_dTsdV)),
	m_Ts(std::move(other.m_Ts))
{
}

template <typename CAESReal>
std::vector<CAESReal> TotalEquilibrium<CAESReal, true>::concentrations(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, const RealVec *analyticalConcentrations) const
{
	CAESReal X = 0.0;
	std::vector<CAESReal> concentrations{};
	const auto ts = Ts(v, activityCoefficients, X); /* Recalculate Ts */

	concentrations.resize(ts.size());
	size_t ctr = 0;
	const CAESReal &c = analyticalConcentrations->elem(concentrationIndex);
	for (const CAESReal &T : ts) {
		const CAESReal fC = c * T / X;

		concentrations[ctr] = fC;

		ctr++;
	}

	return concentrations;
}

template <typename CAESReal>
const std::vector<CAESReal> & TotalEquilibrium<CAESReal, false>::concentrations(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, const RealVec *analyticalConcentrations)
{
	CAESReal X = 0.0;
	Ts(v, activityCoefficients, X); /* Recalculate Ts */

	size_t ctr = 0;
	const CAESReal &c = analyticalConcentrations->elem(concentrationIndex);
	for (const CAESReal &T : m_Ts) {
		const CAESReal fC = c * T / X;

		m_concentrations[ctr] = fC;

		ctr++;
	}

	return m_concentrations;
}

template <typename CAESReal>
typename TotalEquilibrium<CAESReal, true>::DDPack TotalEquilibrium<CAESReal, true>::distributionAnddDistdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients) const
{
	const size_t len = Ls.size();

	CAESReal X = 0.0;
	CAESReal dX = 0.0;
	const auto ts = Ts(v, activityCoefficients, X);		/* Recalculate Ts */
	const auto dts = dTsdV(v, activityCoefficients, dX);	/* Recalculate dTsdV */

	std::vector<CAESReal> distribution;
	std::vector<CAESReal> dDistdV;

	distribution.resize(len);
	dDistdV.resize(len);

	size_t ctr = 0;

	for (size_t idx = 0; idx < len; idx++) {
		const CAESReal &T = ts[idx];
		const CAESReal &dT = dts[idx];

		/* Distribution */
		const CAESReal fC = T / X;
		distribution[ctr] = fC;

		/* dDistdV */
		const CAESReal fD = (dT * X - T * dX) / (X * X);
		dDistdV[idx] = fD;
	}

	return { distribution, dDistdV };
}

template <typename CAESReal>
typename TotalEquilibrium<CAESReal, false>::DDPack TotalEquilibrium<CAESReal, false>::distributionAnddDistdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients)
{
	const size_t len = Ls.size();

	CAESReal X = 0.0;
	CAESReal dX = 0.0;
	Ts(v, activityCoefficients, X);		/* Recalculate Ts */
	dTsdV(v, activityCoefficients, dX);	/* Recalculate dTsdV */

	size_t ctr = 0;

	for (size_t idx = 0; idx < len; idx++) {
		const CAESReal &T = m_Ts[idx];
		const CAESReal &dT = m_dTsdV[idx];

		/* Distribution */
		const CAESReal fC = T / X;
		m_distribution[ctr] = fC;

		/* dDistdV */
		const CAESReal fD = (dT * X - T * dX) / (X * X);
		m_dDistdV[idx] = fD;
	}

	return { m_distribution, m_dDistdV };
}

template <typename CAESReal>
std::vector<CAESReal> TotalEquilibrium<CAESReal, true>::distribution(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients) const
{
	CAESReal X = 0.0;
	const auto ts = Ts(v, activityCoefficients, X); /* Recalculate Ts */
	std::vector<CAESReal> distribution;

	distribution.resize(ts.size());
	size_t ctr = 0;
	for (const CAESReal &T : ts) {
		const CAESReal fC = T / X;

		distribution[ctr] = fC;

		ctr++;
	}

	return distribution;
}

template <typename CAESReal>
const std::vector<CAESReal> & TotalEquilibrium<CAESReal, false>::distribution(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients)
{
	CAESReal X = 0.0;
	Ts(v, activityCoefficients, X); /* Recalculate Ts */

	size_t ctr = 0;
	for (const CAESReal &T : m_Ts) {
		const CAESReal fC = T / X;

		m_distribution[ctr] = fC;

		ctr++;
	}

	return m_distribution;
}

template <typename CAESReal>
std::vector<CAESReal> TotalEquilibrium<CAESReal, true>::dTsdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X) const
{
	std::vector<CAESReal> dts{};

	const size_t len = Ls.size();

	assert(len == static_cast<size_t>(numHigh - numLow) + 1);

	dts.resize(len);
	calculatedTsdV(dts, v, activityCoefficients, X, Ls, numLow);

	return dts;
}

template <typename CAESReal>
std::vector<CAESReal> & TotalEquilibrium<CAESReal, false>::dTsdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X)
{
	assert(m_dTsdV.size() == static_cast<size_t>(numHigh - numLow) + 1);
	calculatedTsdV(m_dTsdV, v, activityCoefficients, X, Ls, numLow);

	return m_dTsdV;
}

template <typename CAESReal>
std::vector<CAESReal> TotalEquilibrium<CAESReal, true>::Ts(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X) const
{
	std::vector<CAESReal> ts{};

	const size_t len = Ls.size();

	assert(len == static_cast<size_t>(numHigh - numLow + 1));

	ts.resize(len);
	calculateTs(ts, v, activityCoefficients, X, Ls, numLow);

	return ts;
}

template <typename CAESReal>
const std::vector<CAESReal> & TotalEquilibrium<CAESReal, false>::Ts(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X)
{
	assert(m_Ts.size() == static_cast<size_t>(numHigh - numLow + 1));

	calculateTs(m_Ts, v, activityCoefficients, X, Ls, numLow);

	return m_Ts;
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_TOTALEQUILIBRIUM_HPP


