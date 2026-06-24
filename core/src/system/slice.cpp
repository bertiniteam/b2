#include "bertini2/system/slice.hpp"
#include "bertini2/system/system.hpp"

namespace bertini {

	void Slice::AddTo(System & s) const
	{
		s.AddBlock(block_);
	}

	System Slice::AsSystem() const
	{
		System s;
		s.AddVariableGroup(sliced_vars_);
		AddTo(s);
		return s;
	}

	std::ostream& operator<<(std::ostream& out, Slice const& s)
	{
		out << "linear slice on " << s.NumVariables() << " variables:\n";
		for (auto& v : s.sliced_vars_)
			out << *v << " ";

		out << "\n\naugmented coefficient matrix (last column is the constant term):\n\n";
		out << s.Coefficients() << "\n\n";

		out << (s.is_homogeneous_ ? "slice is homogeneous" : "slice is not homogeneous") << "\n";

		return out;
	}
}
