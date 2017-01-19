#include "bertini2/system/slice.hpp"

namespace bertini {
	
	std::ostream& operator<<(std::ostream& out, LinearSlice const& s)
	{
		out << "linear slice on " << s.NumVariables() << " variables:\n";
		for (auto& v : s.sliced_vars_)
			out << *v << " ";

		out << "\n\ncoefficient matrix:\n\n";
		out << std::get<Mat<dbl> >(s.coefficients_working_) << "\n\n";

		
		if (s.is_homogeneous_)
			out << "slice is homogeneous";
		else
		{
			out << "slice is not homogeneous, with constants\n";
			out << std::get<Vec<dbl> >(s.constants_working_) << "\n";
		}
		
		return out;
	}
}