#ifndef ECHMET_ELMIG_PARAMS_FITTER_H
#define ECHMET_ELMIG_PARAMS_FITTER_H

#include <echmetelems.h>
#include <echmetsyscomp.h>
#include <new>

namespace ECHMET {

/*!
 * Package that implements electromigration parameters fitting calculator
 */
namespace ElmigParamsFitter {

/*!
 * Enumeration of allowed return codes
 */
ECHMET_ST_ENUM(RetCode) {
	OK = 0,					/*!< Success */
	E_INVALID_ARGUMENT = 0x1,		/*!< Invalid argument passed to a function */
	E_NOT_ENOUGH_MEASUREMENTS = 0x2,	/*!< Number of conducted measurements is insufficient
						     to fit all parameters */
	E_INVALID_BUFFER = 0x4,			/*!< Processed system contains invalid buffer */
	E_NO_MEMORY = 0x5,			/*!< Insufficient memory */
	E_REGRESSOR_INITIALIZATION = 0x6,	/*!< Failed to initialize regressor with the given
						     set of input date */
	E_REGRESSOR_NO_SOLUTION = 0x7,		/*!< Regressor failed to converge */
	E_REGRESSOR_INTERNAL_ERROR = 0x8,	/*!< Unspecified regressor error */
	E_REGRESSOR_PARAMETERS_NOT_SANE = 0x9,	/*!< Initial estimates are not sane */
	ENUM_FORCE_INT32_SIZE(ElmigParamsFitterRetCode)
};

/*!
 * Fixable parameters of fit
 */
ECHMET_ST_ENUM(FixedParameterType) {
	FPT_MOBILITY = 0,			/*!< Fix mobility */
	FPT_PKA = 0x1,				/*!< Fix pKa */
	ENUM_FORCE_INT32_SIZE(ElmigParamsFitterFixedParameterType)
};

/*!
 * Fit options
 */
ECHMET_WK_ENUM(FitOptions) {
	FO_DISABLE_MOB_CONSTRAINTS = (1 << 0),	/*!< Disables constraning of limit mobilities to physically
						     sensible values */
	FO_UNSCALED_STDERRS = (1 << 1),		/*!< Use unscaled standard errors */
	ENUM_FORCE_INT32_SIZE(ElmigParamsFitterFitOptions)
};

template <typename T>
bool isOptionSet(const T &item, const T &options)
{
#ifdef ECHMET_HAVE_CPP11
	using UT = std::underlying_type<T>::type;

	return static_cast<UT>(item) & static_cast<UT>(options);
#else
	return static_cast<int32_t>(item) & static_cast<int32_t>(options);
#endif // ECHMET_HAVE_CPP11
}

/*!
 * Description of a tracepoint.
 */
class TracepointInfo {
public:
	int32_t id;				/*!< Internal ID of the tracepoint. Used to
						     set the tracepoint state. */
	FixedString *description;		/*!< Human-readable description of the tracepoint. */
};
IS_POD(TracepointInfo)

/*!
 * <tt>Vec</tt> of <tt>TracepointInfo</tt>s
 */
typedef Vec<TracepointInfo> TracepointInfoVec;

/*!
 * Input experimental conditions for a given background buffer
 */
class InBuffer {
public:
	SysComp::InConstituentVec *composition;	/*!< Buffer composition */
	RealVec *concentrations;		/*!< Concentrations of buffer constituents.
						     Order of concentrations shall be the same
						     as the order of constituents in <tt>composition</tt>
						     vector. */
	double uEffExperimental;		/*!< Measured mobility of the analyte */
};
IS_POD(InBuffer)

/*!
 * <tt>Vec</tt> of <tt>InBuffer</tt>s
 */
typedef Vec<InBuffer> InBufferVec;

/*!
 * Input description of the system to be processed
 */
class InSystem {
public:
	InBufferVec *buffers;			/*!< Vector of buffers where the measurements were conducted. */
	SysComp::InConstituent analyte;		/*!< Analyte with estimated mobility and pKa values */
	NonidealityCorrections corrections;	/*!< Enabled ionic effects corrections */
};
IS_POD(InSystem)

/*!
 * Point on the resulting fitted curve
 */
class ExpectedCurvePoint {
public:
	double pH;				/*!< Independent variable */
	double experimental;			/*!< Measured (input) experimental mobility */
	double expected;			/*!< Expected mobility computed from the fitted parameters */
};
IS_POD(ExpectedCurvePoint)

/*!
 * <tt>Vec</tt> of <tt>ExpectedCurvePoint</tt>s
 */
typedef Vec<ExpectedCurvePoint> ExpectedCurvePointVec;

/*!
 * Computed fitted parameter
 */
class FittedParameter {
public:
	int charge;				/*!< Respective analyte charge */
	double value;				/*!< Value of the parameter */
	double stdErr;				/*!< Standard error of the value. */
};
IS_POD(FittedParameter)

/*!
 * <tt>Vec</tt> of <tt>FittedParameter</tt>s
 */
typedef Vec<FittedParameter> FittedParameterVec;

/*!
 * Results of a successful fit.
 * Values for zero charge are excluded as the mobility of
 * uncharged state is always expected to be zero.
 */
class FitResults {
public:
	FittedParameterVec *mobilities;		/*!< Fitted mobilities */
	FittedParameterVec *pKas;		/*!< Fitted pKas */
	double rSquared;			/*!< Coefficient of determination */
};
IS_POD(FitResults)

/*!
 * Helper gadget for fixing parameters of fit
 */
class ParametersFixer {
public:
	/*!
	 * Fixes a parameter.
	 *
	 * @param[in] type Type of parameter to fix. Can be either <tt>FPT_MOBILITy</tt> or <tt>FPT_PKA</tt>.
	 * @param[in] id Analyte charge for which to fix the parameter.
	 * @param[in] value Fixed value
	 *
	 * @retval RetCode::OK Success
	 * @retval RetCode::E_NO_MEMORY Insufficient memory
	 */
	virtual RetCode ECHMET_CC add(const FixedParameterType type, const int id, const double value) ECHMET_NOEXCEPT = 0;

	/*!
	 * Destroy the object
	 */
	virtual void ECHMET_CC destroy() ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns whether a parameter is fixed.
	 *
	 * @param[in] type Type of parameter to fix. Can be either <tt>FPT_MOBILITy</tt> or <tt>FPT_PKA</tt>.
	 * @param[in] id Analyte charge for which to fix the parameter.
	 *
	 * @return true if the parameter is fixed, false otherwise
	 */
	virtual bool ECHMET_CC isFixed(const FixedParameterType type, const int id) const ECHMET_NOEXCEPT = 0;

	/*!
	 * Unfixes a fixed parameter.
	 *
	 * @param[in] type Type of parameter to fix. Can be either <tt>FPT_MOBILITy</tt> or <tt>FPT_PKA</tt>.
	 * @param[in] id Analyte charge for which to fix the parameter.
	 */
	virtual void ECHMET_CC remove(const FixedParameterType type, const int id) ECHMET_NOEXCEPT = 0;

	/*!
	 * Returns value of a fixed parameter.
	 *
	 * @param[in] type Type of parameter to fix. Can be either <tt>FPT_MOBILITy</tt> or <tt>FPT_PKA</tt>.
	 * @param[in] id Analyte charge for which to fix the parameter.
	 *
	 * @return Value of a fixed parameter. If the given parameter is not fixed, zero is returned.
	 */
	virtual double ECHMET_CC value(const FixedParameterType type, const int id) const ECHMET_NOEXCEPT = 0;

protected:
	virtual ECHMET_CC ~ParametersFixer() ECHMET_NOEXCEPT = 0;
};

extern "C" {

/*!
 * Checks whether analyte estimates are sane.
 *
 * Sane estimates must meet the following requirements:
 *
 * - All pKa values must be descending
 * - All charged forms must have non-zero mobility.
 * - Uncharged form must have zero mobility
 * - Mobility of consecutive forms must fall into range given by
 *   <tt>mobilityLowerBound()</tt> and <tt>mobilityUpperBound()</tt>.
 *   This requirement may be overriden by setting <tt>FO_DISABLE_MOB_CONSTRAINTS</tt>
 *   fit option.
 *
 * @param[in] analyte Analyte to be checked.
 * @param[in] options Additional options that modify behavior of the fit algorithm.
 *
 * @return true if estimates are sane, false otherwise
 */
ECHMET_API bool ECHMET_CC checkSanity(const SysComp::InConstituent &analyte, const FitOptions options) ECHMET_NOEXCEPT;

/*!
 * Creates empty vector if input buffers.
 *
 * @return Pointer to the buffers vector or <tt>NULL</tt> if the operation fails.
 */
ECHMET_API InBufferVec * ECHMET_CC createInBufferVec() ECHMET_NOEXCEPT;

/*!
 * Creates parameter fixing gadget.
 *
 * @return Pointer to the gadget or <tt>NULL</tt> if the operation fails.
 */
ECHMET_API ParametersFixer * ECHMET_CC createParametersFixer() ECHMET_NOEXCEPT;

/*!
 * Default fit options.
 *
 * @return Default fit options
 */
ECHMET_API FitOptions ECHMET_CC defaultFitOptions() ECHMET_NOEXCEPT;

/*!
 * Converts return code to respective human-readable error string.
 *
 * @param[in] tRet Return code to convert
 *
 * @return Human-readable error string
 */
ECHMET_API const char * ECHMET_CC EMPFerrorToString(const RetCode tRet) ECHMET_NOEXCEPT;

/*!
 * Computes expected pH vs. mobility curve.
 *
 * @param[in] system Description of the processed system.
 * @param[in] results Computed results.
 * @param[in,out] curve Pointer to vector of <tt>ExpectedCurvePoint</tt>s.
 *
 * @retval RetCode::E_MO_MEMORY Insufficient memory
 * @retval RetCode::E_REGRESSOR_INTERNAL_ERROR Unable to compute expected mobility
 * @retval RetCode::E_INVALID_BUFFER Invalid input buffer
 */
ECHMET_API RetCode ECHMET_CC expectedCurve(const InSystem &system, const FitResults &results, ExpectedCurvePointVec *&curve) ECHMET_NOEXCEPT;

/*!
 * Value of relative lower mobility bound
 *
 * @return Value of relative lower mobility bound
 */
ECHMET_API double ECHMET_CC mobilityLowerBound() ECHMET_NOEXCEPT;

/*!
 * Value of relative upper mobility bound
 *
 * @return Value of relative upper mobility bound
 */
ECHMET_API double ECHMET_CC mobilityUpperBound() ECHMET_NOEXCEPT;

/*!
 * Performs the fit.
 *
 * @param[in] system Input description of the system to be processed.
 * @param[in] fixer Parameter fixing gadget. If the value is <tt>NULL</tt> all parameters are free.
 * @param[in] options Additional options that modify behavior of the fit algorithm
 * @param[in,out] results Fitted parameters.
 */
ECHMET_API RetCode ECHMET_CC process(const InSystem &system, const ParametersFixer *fixer, const FitOptions options,
				     FitResults &results) ECHMET_NOEXCEPT;

/*!
 * Convenience function to properly release resources claimed by <tt>InBuffer</tt> object.
 *
 * @param[in] buffer <tt>InBuffer</tt> to be released.
 */
ECHMET_API void ECHMET_CC releaseInBuffer(const InBuffer &buffer) ECHMET_NOEXCEPT;

/*!
 * Convenience function to properly release all resources claimed by a vector
 * of <tt>InBuffer</tt>s. Contents of the vector are released as well.
 *
 * @param[in] vec <tt>InBufferVec</tt> to be released.
 */
ECHMET_API void ECHMET_CC releaseInBufferVec(const InBufferVec *vec) ECHMET_NOEXCEPT;

/*!
 * Convenience function to release all resources claimed by <tt>InSystem</tt> object.
 *
 * @param[in] system <tt>InSystem</tt> to be released.
 */
ECHMET_API void ECHMET_CC releaseInSystem(InSystem &system) ECHMET_NOEXCEPT;

/*!
 * Convenience function to release all resources claimed by <tt>Resutls</tt> object.
 *
 * @param[in] results <tt>Results</tt> object to be released,
 */
ECHMET_API void ECHMET_CC releaseResults(FitResults &results) ECHMET_NOEXCEPT;

/*!
 * Sets all tracepoints to the given state.
 *
 * @param[in] state If \p true all tracepoints will be enabled and vice versa.
 */
ECHMET_API void ECHMET_CC toggleAllTracepoints(const bool state) ECHMET_NOEXCEPT;

/*!
 * Set state of one tracepoint.
 *
 * @param[in] TPID Internal ID of the tracepoint to set.
 * @param[in] state If \p true the tracepoint will be enabled and vice versa.
 */
ECHMET_API void ECHMET_CC toggleTracepoint(const int32_t TPID, const bool state) ECHMET_NOEXCEPT;

/*!
 * Returns the complete trace.
 *
 * @param[in] dontClear If \p true the trace log will not be cleared.
 *
 * @return String containing the whole trace.
 */
ECHMET_API FixedString * ECHMET_CC trace(const bool dontClear = false) ECHMET_NOEXCEPT;

/*!
 * Returns information about available tracepoints.
 *
 * @retval Pointer to a vector of all available tracepoints. May be \p NULL if
 *         the operation fails or if no tracepoins are available.
 */
ECHMET_API TracepointInfoVec * ECHMET_CC tracepointInfo() ECHMET_NOEXCEPT;

/*!
 * Returns the state of a given tracepoint.
 *
 * @param[in] TPID Internal ID of the tracepoint whose state is requested.
 *
 * @retval true if the tracepoint is enabled and vice versa.
 */
ECHMET_API bool ECHMET_CC tracepointState(const int32_t TPID) ECHMET_NOEXCEPT;

/*!
 * Returns library version ID string
 * The string is formatted as "MAJOR.MINOR.PATCH"
 * where MAJOR, MINOR and PATCH are positive integers.
 */
ECHMET_API const char * ECHMET_CC versionString() ECHMET_NOEXCEPT;

}

} // namespace ElmigParamsFitter
} // namespace ECHMET

#endif // ECHMET_ELMIG_PARAMS_FITTER_H
