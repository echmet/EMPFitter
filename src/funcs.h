#include <cmath>

template<typename T>
static inline T pX(const T &t)
{
	return -std::log10(t);
}

template<typename T>
static inline T X10(const T &t)
{
	return std::pow<T>(10.0, -t);
}

template<>
inline double X10(const double &t)
{
	return std::exp(-t * M_LN10);
}
