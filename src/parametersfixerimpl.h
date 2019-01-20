#ifndef ECHMET_PARAMETERS_FIXER_IMPL_H
#define ECHMET_PARAMETERS_FIXER_IMPL_H

#include <echmetelmigparamsfitter.h>
#include <map>

namespace ECHMET {
namespace ElmigParamsFitter {

class ParametersFixerImpl : public ParametersFixer {
public:
	~ParametersFixerImpl() noexcept override;
	RetCode ECHMET_CC add(const FixedParameterType type, const int id, const double value) noexcept override;
	void ECHMET_CC destroy() noexcept override;
	bool ECHMET_CC isFixed(const FixedParameterType type, const int id) const noexcept override;
	void ECHMET_CC remove(const FixedParameterType type, const int id) noexcept override;
	double ECHMET_CC value(const FixedParameterType type, const int id) const noexcept override;

private:
	std::map<int, double> m_fixedMobilities;
	std::map<int, double> m_fixedpKas;
};

} // namespace ElmigParamsFitter
} // namespace ECHMET

#endif // ECHMET_PARAMETERS_FIXER_IMPL_H
