#include <sstream>

#include "bertini2/system/slice.hpp"
#include "bertini2/system/system.hpp"

namespace bertini {

	void Slice::AddTo(System & s) const
	{
		// The block's columns are indexed by the system's variable ordering; if the counts disagree
		// the matrix-vector product is meaningless.  Fail here with a clear message rather than
		// producing a silently-wrong (or out-of-bounds) evaluation downstream.
		if (NumVariables() != s.NumVariables())
		{
			std::stringstream ss;
			ss << "cannot add a slice on " << NumVariables() << " variables to a system on "
			   << s.NumVariables() << " variables; build the slice over the system's variables";
			throw std::runtime_error(ss.str());
		}
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
