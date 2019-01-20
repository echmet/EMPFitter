#ifndef ECHMET_CAES_TOTALEQUILIBRIUM_H
#define ECHMET_CAES_TOTALEQUILIBRIUM_H

#include <cassert>
#include <echmetcaes.h>
#include <vector>
#include "funcs.h"

namespace ECHMET {
namespace CAES {

/*!
 * Calculates total equilibrium constants from constecutive constants
 *
 * @param[in] pBs Consecutive equilibrium constants
 *
 * @return Vector of total equilibrium constants
 */
template <typename CAESReal>
std::vector<CAESReal> calculateLsBase(const std::vector<CAESReal> &pBs)
{
	std::vector<CAESReal> _Ls;
	const size_t len = pBs.size() + 1;

	_Ls.resize(len);
	_Ls[0] = 1.0;

	for (size_t idx = 1; idx < len; idx++)
		_Ls[idx] = _Ls[idx - 1] / X10(pBs[idx - 1] - 3.0); /* Use 1/L values so we can do vPow * L instead of vPow / L later on */

	return _Ls;
}

class TotalEquilibriumBase {
public:
	virtual ~TotalEquilibriumBase() {}
};

/*!
 * Description of total distribution equilibria for a given constituent
 */
template <typename CAESReal, bool ThreadSafe>
class TotalEquilibrium : public TotalEquilibriumBase
{};

template <typename CAESReal>
class TotalEquilibrium<CAESReal, true> : public TotalEquilibriumBase {
public:
	typedef std::tuple<std::vector<CAESReal>, std::vector<CAESReal>> DDPack;

	TotalEquilibrium();

	/*!
	 * TotalEquilibrium c-tor
	 *
	 * @param[in] numLow Lowest equilibrium index
	 * @param[in] numHigh Highest equilibirium index
	 * @param[in] pBs Vector of consecutive equilibrium constants
	 * @param[in] concentrationIndex Index of analytical concentration of the constituent
	 */
	template <typename T = CAESReal, typename std::enable_if<!std::is_same<ECHMETReal, T>::value>::type>
	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<ECHMETReal> &pBs, const size_t concentrationIndex) :
		concentrationIndex(concentrationIndex),
		Ls(calculateLs(pBs)),
		numLow(numLow),
		numHigh(numHigh)
	{
	}

	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<CAESReal> &pBs, const size_t concentrationIndex = SIZE_MAX);
	TotalEquilibrium(const TotalEquilibrium &other);
	TotalEquilibrium(TotalEquilibrium &&other);
	std::vector<CAESReal> concentrations(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, const RealVec *analyticalConcentrations) const;
	DDPack distributionAnddDistdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients) const;
	std::vector<CAESReal> distribution(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients) const;
	std::vector<CAESReal> dTsdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X) const;
	std::vector<CAESReal> Ts(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X) const;

	const size_t concentrationIndex;	/*!< Analytical concentration index of the constituent */
	const std::vector<CAESReal> Ls;		/*!< Vector of total equilibirum constants */
	const int numLow;			/*!< Lowest equilibrium index */
	const int numHigh;			/*!< Highest equilibrium index */

private:
	/*!
	 * Calculates total equilibrium constants from constecutive constants
	 *
	 * @param[in] pBs Consecutive equilibrium constants
	 *
	 * @return Vector of total equilibrium constants
	 */
	template <typename T = CAESReal, typename std::enable_if<!std::is_same<ECHMETReal, T>::value>::type>
	std::vector<CAESReal> calculateLs(const std::vector<ECHMETReal> &pBs)
	{
		std::vector<CAESReal> _pBs;
		_pBs.reserve(pBs.size());

		for (const ECHMETReal &d : pBs)
			_pBs.emplace_back(d);

		return calculateLsBase<CAESReal>(_pBs);
	}

	std::vector<CAESReal> calculateLs(const std::vector<CAESReal> &pBs)
	{
		return calculateLsBase<CAESReal>(pBs);
	}
};

template <typename CAESReal>
class TotalEquilibrium<CAESReal, false> : public TotalEquilibriumBase {
public:
	typedef std::tuple<const std::vector<CAESReal> &, const std::vector<CAESReal> &> DDPack;

	TotalEquilibrium();

	/*!
	 * TotalEquilibrium c-tor
	 *
	 * @param[in] numLow Lowest equilibrium index
	 * @param[in] numHigh Highest equilibirium index
	 * @param[in] pBs Vector of consecutive equilibrium constants
	 * @param[in] concentrationIndex Index of analytical concentration of the constituent
	 */
	template <typename T = CAESReal, typename std::enable_if<!std::is_same<ECHMETReal, T>::value>::type>
	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<ECHMETReal> &pBs, const size_t concentrationIndex) :
		concentrationIndex(concentrationIndex),
		Ls(calculateLs(pBs)),
		numLow(numLow),
		numHigh(numHigh)
	{
		m_concentrations.resize(Ls.size());
		m_distribution.resize(Ls.size());
		m_dDistdV.resize(Ls.size());
		m_dTsdV.resize(Ls.size());
		m_Ts.resize(Ls.size());
	}

	TotalEquilibrium(const int numLow, const int numHigh, const std::vector<CAESReal> &pBs, const size_t concentrationIndex = SIZE_MAX);
	TotalEquilibrium(const TotalEquilibrium &other);
	TotalEquilibrium(TotalEquilibrium &&other);
	const std::vector<CAESReal> & concentrations(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, const RealVec *analyticalConcentrations);
	DDPack distributionAnddDistdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients);
	const std::vector<CAESReal> & distribution(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients);
	std::vector<CAESReal> & dTsdV(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X);
	const std::vector<CAESReal> & Ts(const CAESReal &v, const std::vector<CAESReal> &activityCoefficients, CAESReal &X);

	const size_t concentrationIndex;	/*!< Analytical concentration index of the constituent */
	const std::vector<CAESReal> Ls;		/*!< Vector of total equilibirum constants */
	const int numLow;			/*!< Lowest equilibrium index */
	const int numHigh;			/*!< Highest equilibrium index */

private:
	/*!
	 * Calculates total equilibrium constants from constecutive constants
	 *
	 * @param[in] pBs Consecutive equilibrium constants
	 *
	 * @return Vector of total equilibrium constants
	 */
	template <typename T = CAESReal, typename std::enable_if<!std::is_same<ECHMETReal, T>::value>::type>
	std::vector<CAESReal> calculateLs(const std::vector<ECHMETReal> &pBs)
	{
		std::vector<CAESReal> _pBs;
		_pBs.reserve(pBs.size());

		for (const ECHMETReal &d : pBs)
			_pBs.emplace_back(d);

		return calculateLs<CAESReal>(_pBs);
	}

	std::vector<CAESReal> calculateLs(const std::vector<CAESReal> &pBs)
	{
		return calculateLsBase<CAESReal>(pBs);
	}


	std::vector<CAESReal> m_concentrations;
	std::vector<CAESReal> m_distribution;
	std::vector<CAESReal> m_dDistdV;
	std::vector<CAESReal> m_dTsdV;
	std::vector<CAESReal> m_Ts;
};

}
}

#include "totalequilibrium.hpp"

#endif // ECHMET_CAES_TOTALEQUILIBRIUM_H
