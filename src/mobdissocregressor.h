#ifndef ECHMET_MOBDISSOCREGRESSOR_H
#define ECHMET_MOBDISSOCREGRESSOR_H

#include "types.h"

#include <regress.h>

namespace ECHMET {

class MobDissocRegressor : public RegressCore::RegressFunction<BufferSystem, double> {
public:
	using RFType = RegressFunction<BufferSystem, double>;

	using YTMatrix = RFType::YTMatrix;
	using XTVector = RFType::XTVector;
	using YTVector = RFType::YTVector;
	using ParamsVector = RFType::ParamsVector;
	using PVSize = ParamsVector::size_type;
	using msize_t = RegressCore::msize_t;

	MobDissocRegressor(const BufferSystemVec &bufSysVec, const InConstituentWrapper &analyte,
			   const NonidealityCorrections corrections);

	/*!
	 * Assignment operator analog.
	 *
	 * @param otherBase Object to copy from.
	 */
	void AAssign(const RFType &otherBase) override;

	/*!
	 * Calculates derivative of the function by the given parameter.
	 *
	 * @param x X-value of the function
	 * @param params Set of function's parameters
	 * @param param_idx Index of the parameter to calculate derivative for.
	 * @param idx UNKNOWN
	 */
	double ACalculateDerivative(const BufferSystem &x, const ParamsVector &params, const typename RFType::index_type param_idx, msize_t idx) const override;

	/*!
	 * Calculates Y value of the function for given X and parameters.
	 *
	 * @param x X value
	 * @param params Function's parameters
	 * @param idx UNKNOWN
	 *
	 * @return Calculated Y value
	 */
	double ACalculateFx(const BufferSystem &x, const ParamsVector &params, msize_t idx) const override;

	/*!
	 * Checks if the function's parameters are sane.
	 *
	 * @return true of the parameters are sane, false otherwise
	 */
	bool ACheckSanity(const ParamsVector &params) const override;

	/*!
	 * Creates a new instance of the regressor.
	 *
	 * @return Pointer to a new instance of the RegressFunction object or \p nullptr
	 *         if the new instance creation failed.
	 */
	MobDissocRegressor * ACreate() const override;

	/*!
	 * Sets initial parameters for the fit to start with.
	 * Called by parent's internal Initialize().
	 *
	 * @param params Reference to array of parameters used by RegressFunction
	 * @param x Vector of input X values
	 * @param y Vector of input Y values
	 *
	 * @return true of the initial estimate is sensible, false otherwise
	 */
	bool ASetInitialParameters(ParamsVector &params, const XTVector &x, const YTVector &y) override;

	bool CheckSanity();

	bool Initialize(const std::vector<double> &inYVec);

	static bool CheckAnalyteSanity(const SysComp::InConstituent &analyte);

	/*!
	 * Calculates effective mobility and optionally pH for a given combination
	 * of buffer and analyte.
	 *
	 * @param[in] composition Composition of the background buffer.
	 * @param[in] analyte Physical-chemical properties of the analyte whose mobility is to be calculated.
	 * @param[in] concentrations Vector of analytical concentrations of background buffer components.
	 * @param[in] cH <tt>first</tt> field controls whether the cH value is set to a valid value,
	 *               <tt>second</tt> field contains the actual value of <tt>cH</tt>.
	 * @param[in] corrs Nonideality corrections that shall be applied.
	 *
	 * @return Effective mobility and pH as <tt>std::pair&lt;double, double&gt;</tt>.
	 *         Effective mobility is first, pH second.
	 */
	static
	std::pair<double, double> SolveBuffer(const SysComp::InConstituentVec *composition,
					      const SysComp::InConstituent &analyte,
					      const RealVec *concentrations,
					      const std::pair<bool, double> cH,
					      const NonidealityCorrections corrs);

	constexpr static const double MOBILITY_LOWER_BOUND{0.55};
	constexpr static const double MOBILITY_UPPER_BOUND{1.2};

private:
	/*!
	 * Applies set of function parameter to the analyte object.
	 *
	 * @param analyte Analyte object to operate upon.
	 * @param params Set of parameters to be applied.
	 */
	void ApplyParameters(SysComp::InConstituent &analyte, const ParamsVector &params) const;

	/*!
	 * Calculates effective mobility of analyte for a given set of paramers.
	 *
	 * @param x Buffer system to use in the calculation.
	 * @param params Physical-chemical properties of the analyte.
	 *
	 * @return Effective mobility
	 */
	double CalculateEffectiveMobility(const BufferSystem &x, const ParamsVector &params) const;

	void SetParameters(ParamsVector &params, const SysComp::InConstituent &analyte);

	static bool CheckSanityInternal(const ParamsVector &params, const SysComp::InConstituent &analyte,
					const size_t NMobilities) noexcept;

	static PVSize CountNumberOfMobilities(const SysComp::InConstituent &analyte) noexcept;

	static PVSize CountNumberOfParams(const SysComp::InConstituent &analyte) noexcept;

	static void SetParametersRaw(ParamsVector &params, const SysComp::InConstituent &analyte,
				     const PVSize NParams);

	BufferSystemVec m_bufSysVec;
	InConstituentWrapper m_analyte;

	NonidealityCorrections m_solverCorrections;

	const PVSize m_NMobilities;
	const PVSize m_NParams;
};

} // namespace ECHMET

#endif // ECHMET_MOBDISSOCREGRESSOR_H
