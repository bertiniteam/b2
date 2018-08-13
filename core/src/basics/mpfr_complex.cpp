#include "bertini2/mpfr_complex.hpp"

#if USE_BMP_COMPLEX
	//nothing to see here, folks.  move along.
#else

namespace bertini{
	#ifdef USE_THREAD_LOCAL
		mpfr_float thread_local custom_complex::temp_[8]{};
	#else
		mpfr_float custom_complex::temp_[8]{};
	#endif
}

#endif
