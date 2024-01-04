/**
 * !file omp_utils.h
 * \brief OpenMP utilities
 */

#ifndef OMP_UTILS_H
#define OMP_UTILS_H

#ifdef USEOPENMP
#include <omp.h>
#endif // USEOPENMP

namespace NBody
{

#ifdef USEOPENMP
int get_available_threads()
{
	int nthreads;
	#pragma omp parallel
	#pragma omp single
	{
		nthreads = omp_get_num_threads();
	}
	return nthreads;
}
#else
int get_available_threads()
{
	return 1;
}
#endif // USEOPENMP

#ifdef USEOPENMP
/// Returns whether OpenMP nesting is enabled or not
static bool _omp_get_nested()
{
#if _OPENMP >= 200805
	return omp_get_max_active_levels() > 1;
#else
	return omp_get_nested();
#endif
}

/// Set OpenMP to enabled (or not) nested parallelism
static void _omp_set_nested(bool enable)
{
	static constexpr int MAX_OPENMP_ACTIVE_LEVELS = 20;
#if _OPENMP >= 200805
	omp_set_max_active_levels(enable ? MAX_OPENMP_ACTIVE_LEVELS : 1);
#else
	omp_set_nested(int(enable));
#endif
}
#endif // USEOPENMP

#ifdef USEOPENMP
/**
 * A class that enables OpenMP nested calls upon construction, if requested,
 * and restores the previous behavior upon destruction.
 */
class OmpNestedEnabler
{
public:
	OmpNestedEnabler(bool enable)
	  : nested_previously_enabled(_omp_get_nested()),
	    _available_threads(get_available_threads())
	{
		if (nested_previously_enabled || _available_threads == 1) {
			return;
		}
		_omp_set_nested(enable);
	}

	~OmpNestedEnabler()
	{
		_omp_set_nested(nested_previously_enabled);
	}

	int available_threads()
	{
		return _available_threads;
	}

private:
	bool nested_previously_enabled;
	int _available_threads;
};
#else
class OmpNestedEnabler
{
public:
	OmpNestedEnabler(bool enable)
	{
	}

	int available_threads()
	{
		return get_available_threads();
	}
};
#endif

}

#endif // OMP_UTILS_H