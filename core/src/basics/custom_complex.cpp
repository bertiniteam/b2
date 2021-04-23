#include "bertini2/custom_complex.hpp"


namespace bertini{
	#ifdef USE_THREAD_LOCAL
		mpfr_float thread_local custom_complex::temp_[8]{};
	#else
		mpfr_float custom_complex::temp_[8]{};
	#endif
}
