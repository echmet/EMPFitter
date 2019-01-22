#ifndef ECHMET_TT_TYPES_H
#define ECHMET_TT_TYPES_H

#include <cassert>
#include <echmetsyscomp.h>
#include <echmetcaes.h>
#include <echmetionprops.h>
#include <ostream>
#include <vector>

namespace ECHMET {

typedef std::pair<std::vector<double>, std::vector<double>> ExperimentalData;

class BufferSystem {
public:
	class BufferSystemException : public std::runtime_error {
	public:
		using std::runtime_error::runtime_error;
	};

	BufferSystem(const SysComp::InConstituentVec *composition, std::vector<double> _concentrations,
		     const bool debyeHuckel) :
		concentrations(std::move(_concentrations)),
		m_debyeHuckel(debyeHuckel),
		m_movedAway(false)
	{
		m_composition = SysComp::duplicateInConstituentVec(composition);
		if (composition == nullptr)
			throw std::bad_alloc();

		try {
			m_concentrationsRVec = makeConcentrationsVec(this->concentrations);
		} catch (const std::bad_alloc &) {
			SysComp::releaseInputData(m_composition);

			throw;
		}

		calculateEquilibrium();
	}

	BufferSystem(SysComp::InConstituentVec * &&composition, std::vector<double> _concentrations,
		     const bool debyeHuckel) :
		concentrations(std::move(_concentrations)),
		m_composition(composition),
		m_debyeHuckel(debyeHuckel),
		m_movedAway(false)
	{
		try {
			m_concentrationsRVec = makeConcentrationsVec(this->concentrations);
		} catch (const std::bad_alloc &) {
			SysComp::releaseInputData(m_composition);

			throw;
		}

		calculateEquilibrium();
	}

	BufferSystem(const BufferSystem &other) :
		concentrations(other.concentrations),
		m_cH(other.m_cH),
		m_pH(other.m_pH),
		m_debyeHuckel(other.m_debyeHuckel),
		m_movedAway(other.m_movedAway)
	{
		m_composition = SysComp::duplicateInConstituentVec(other.m_composition);
		if (m_composition == nullptr)
			throw std::bad_alloc();

		m_concentrationsRVec = other.m_concentrationsRVec->duplicate();
		if (m_concentrationsRVec == nullptr) {
			m_composition->destroy();

			throw std::bad_alloc();
		}
	}

	BufferSystem(BufferSystem &&other) noexcept :
		concentrations(std::move(other.concentrations)),
		m_composition(other.m_composition),
		m_concentrationsRVec(other.m_concentrationsRVec),
		m_cH(other.m_cH),
		m_pH(other.m_pH),
		m_debyeHuckel(other.m_debyeHuckel),
		m_movedAway(other.m_movedAway)
	{
		other.m_movedAway = true;
	}

	~BufferSystem()
	{
		if (!m_movedAway) {
			SysComp::releaseInputData(m_composition);
			m_concentrationsRVec->destroy();
		}
	}

	double cH() const noexcept
	{
		return m_cH;
	}

	const SysComp::InConstituentVec * composition() const noexcept
	{
		return m_composition;
	}

	const RealVec * concentrationsRVec() const noexcept
	{
		return m_concentrationsRVec;
	}

	double pH() const noexcept
	{
		return m_pH;
	}

	std::ostream & operator<<(std::ostream &os) const
	{
		dump(os, *this);

		return os;
	}

	BufferSystem & operator=(const BufferSystem &other)
	{
		auto concs = other.concentrations;
		auto comp = SysComp::duplicateInConstituentVec(other.m_composition);
		if (comp == nullptr)
			throw std::bad_alloc();

		auto cRVec = other.m_concentrationsRVec->duplicate();
		if (cRVec == nullptr) {
			SysComp::releaseInputData(comp);

			throw std::bad_alloc();
		}

		const_cast<std::vector<double>&>(concentrations) = std::move(concs);
		m_composition = comp;
		m_concentrationsRVec = cRVec;
		m_cH = other.m_cH;
		m_pH = other.m_pH;
		m_debyeHuckel = other.m_debyeHuckel;
		m_movedAway = other.m_movedAway;

		return *this;
	}

	BufferSystem & operator=(BufferSystem &&other) noexcept
	{
		const_cast<std::vector<double>&>(concentrations) = std::move(other.concentrations);
		m_composition = other.m_composition;
		m_concentrationsRVec = other.m_concentrationsRVec;
		m_cH = other.m_cH;
		m_pH = other.m_pH;
		m_debyeHuckel = other.m_debyeHuckel;
		m_movedAway = other.m_movedAway;

		other.m_movedAway = true;

		return *this;
	}

	friend std::ostream & operator<<(std::ostream &os, const BufferSystem &bufSys);

	const std::vector<double> concentrations;

private:
	void calculateEquilibrium()
	{
		SysComp::ChemicalSystem chemSystem;
		SysComp::CalculatedProperties calcProps;

		NonidealityCorrections corrections = defaultNonidealityCorrections();
		if (m_debyeHuckel)
			nonidealityCorrectionSet(corrections, NonidealityCorrectionsItems::CORR_DEBYE_HUCKEL);

		auto tRet = SysComp::makeComposition(chemSystem, calcProps, m_composition);
		if (tRet != RetCode::OK)
			throw BufferSystemException("Unable to create system composition: " + std::string(errorToString(tRet)));

		CAES::SolverContext *ctx;
		CAES::Solver *solver;

		tRet = CAES::createSolverContext(ctx, chemSystem);
		if (tRet != RetCode::OK) {
			SysComp::releaseChemicalSystem(chemSystem);
			SysComp::releaseCalculatedProperties(calcProps);

			throw BufferSystemException("Unable to create solver context");
		}

		solver = CAES::createSolver(ctx, CAES::Solver::defaultOptions(), corrections);
		if (solver == nullptr) {
			ctx->destroy();
			SysComp::releaseChemicalSystem(chemSystem);
			SysComp::releaseCalculatedProperties(calcProps);

			throw BufferSystemException("Unable to create solver");
		}

		tRet = solver->estimateDistributionSafe(m_concentrationsRVec, calcProps);
		if (tRet != RetCode::OK) {
			solver->destroy();
			ctx->destroy();
			SysComp::releaseChemicalSystem(chemSystem);
			SysComp::releaseCalculatedProperties(calcProps);

			throw BufferSystemException("Failed to estimate distribution");
		}

		m_cH = calcProps.ionicConcentrations->at(0);
		m_pH = IonProps::calculatepH_direct(m_cH, calcProps.ionicStrength);

		solver->destroy();
		ctx->destroy();
		SysComp::releaseChemicalSystem(chemSystem);
		SysComp::releaseCalculatedProperties(calcProps);
	}

	static void dump(std::ostream &os, const BufferSystem &bufSys)
	{
		assert(bufSys.m_composition->size() == bufSys.concentrations.size());

		os << std::string("---\n");
		os << "pH: " << bufSys.m_pH << "\n";
		for (size_t idx = 0; idx < bufSys.m_composition->size(); idx++) {
			const auto &ctuent = bufSys.m_composition->at(idx);

			os << ctuent.name->c_str() << ", " << bufSys.concentrations.at(idx) << "(mM)\n";
		}
		os << std::string("---\n");
	}

	RealVec * makeConcentrationsVec(const std::vector<double> &src)
	{
		RealVec * rVec = createRealVec(src.size());
		if (rVec == nullptr)
			throw std::bad_alloc();

		for (double d : src) {
			if (rVec->push_back(d) != RetCode::OK) {
				rVec->destroy();

				throw std::bad_alloc();
			}
		}

		return rVec;
	}

	SysComp::InConstituentVec *m_composition;
	RealVec *m_concentrationsRVec;

	double m_cH;
	double m_pH;
	bool m_debyeHuckel;

	bool m_movedAway;
};
typedef std::vector<BufferSystem> BufferSystemVec;

inline
std::ostream & operator<<(std::ostream &os, const BufferSystem &bufSys)
{
	BufferSystem::dump(os, bufSys);

	return os;
}

class InConstituentWrapper {
public:
	class NotInitializedException : public std::runtime_error {
	public:
		using std::runtime_error::runtime_error;
	};

	InConstituentWrapper() :
		m_initialized(false),
		m_movedAway(false)
	{}

	InConstituentWrapper(const SysComp::InConstituent &inC) :
		m_initialized(true),
		m_movedAway(false)
	{
		copyContent(inC);
	}

	InConstituentWrapper(SysComp::InConstituent &&inC) noexcept :
		m_constituent(inC),
		m_initialized(true),
		m_movedAway(false)
	{
	}

	InConstituentWrapper(const InConstituentWrapper &other) :
		m_initialized(other.m_initialized),
		m_movedAway(false)
	{
		if (m_initialized)
			copyContent(other.m_constituent);
	}

	InConstituentWrapper(InConstituentWrapper &&other) noexcept :
		m_initialized(other.m_initialized),
		m_movedAway(false)
	{
		moveContent(std::move(other));
	}

	~InConstituentWrapper()
	{
		if (!m_movedAway && m_initialized)
			SysComp::releaseInConstituent(m_constituent);
	}

	InConstituentWrapper & operator=(const InConstituentWrapper &other)
	{
		m_initialized = other.m_initialized;

		if (m_initialized)
			copyContent(other.m_constituent);

		m_movedAway = false;

		return *this;
	}

	InConstituentWrapper & operator=(InConstituentWrapper &&other) noexcept
	{
		m_initialized = other.m_initialized;

		moveContent(std::move(other));
		m_movedAway = false;

		return *this;
	}

	SysComp::InConstituent & operator()()
	{
		if (!m_initialized)
			throw NotInitializedException("InConstituent not initialized");

		return m_constituent;
	}

	const SysComp::InConstituent & operator()() const
	{
		if (!m_initialized)
			throw NotInitializedException("InConstituent not initialized");

		return m_constituent;
	}

private:
	void copyContent(const SysComp::InConstituent &inC)
	{
		if (SysComp::duplicateInConstituent(m_constituent, inC) != RetCode::OK)
			throw std::bad_alloc();
	}

	void moveContent(InConstituentWrapper &&other) noexcept
	{
		m_constituent = other.m_constituent;
		other.m_movedAway = true;
	}

	SysComp::InConstituent m_constituent;

	bool m_initialized;
	bool m_movedAway;
};

} // namespace ECHMET

#endif // ECHMET_TT_TYPES_H
