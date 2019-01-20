#include "parametersfixerimpl.h"

#include <cassert>

namespace ECHMET {
namespace ElmigParamsFitter {

ParametersFixerImpl::~ParametersFixerImpl() noexcept
{
}

RetCode ECHMET_CC ParametersFixerImpl::add(const FixedParameterType type, const int id, const double value) noexcept
{
	try {
		switch (type) {
		case FixedParameterType::FPT_PKA:
			m_fixedpKas[id] = value;
			return RetCode::OK;
		case FixedParameterType::FPT_MOBILITY:
			m_fixedMobilities[id] = value;
			return RetCode::OK;
		}
	} catch (const std::bad_alloc &) {
		return RetCode::E_NO_MEMORY;
	}

	assert(false);
}

void ECHMET_CC ParametersFixerImpl::destroy() noexcept
{
	delete this;
}

bool ECHMET_CC ParametersFixerImpl::isFixed(const FixedParameterType type, const int id) const noexcept
{
	switch (type) {
	case FixedParameterType::FPT_PKA:
		return m_fixedpKas.find(id) != m_fixedpKas.end();
	case FixedParameterType::FPT_MOBILITY:
		return m_fixedMobilities.find(id) != m_fixedMobilities.end();
	}

	assert(false);
}

void ECHMET_CC ParametersFixerImpl::remove(const FixedParameterType type, const int id) noexcept
{
	switch (type) {
	case FixedParameterType::FPT_PKA:
		m_fixedpKas.erase(id);
		return;
	case FixedParameterType::FPT_MOBILITY:
		m_fixedMobilities.erase(id);
		return;
	}

	assert(false);
}

double ECHMET_CC ParametersFixerImpl::value(const FixedParameterType type, const int id) const noexcept
{
	try {
		switch (type) {
		case FixedParameterType::FPT_PKA:
			return m_fixedpKas.at(id);
		case FixedParameterType::FPT_MOBILITY:
			return m_fixedMobilities.at(id);
		}
	} catch (const std::out_of_range &) {
		return 0.0;
	}

	assert(false);
}

ParametersFixer::~ParametersFixer() noexcept
{
}

} // namespace ElmigParamsFitter
} // namespace ECHMET
