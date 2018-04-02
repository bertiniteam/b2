#include "bertini2/mpfr_complex.hpp"

#if 0

namespace bertini{
	#ifdef USE_THREAD_LOCAL
		mpfr_float thread_local complex::temp_[8]{};
	#else
		mpfr_float complex::temp_[8]{};
	#endif
}

#endif // remove me -- testing
