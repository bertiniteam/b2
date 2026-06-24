#include <sstream>

#include "bertini2/system/slice.hpp"
#include "bertini2/system/system.hpp"

namespace bertini {

	void Slice::AddTo(System & s) const
	{
		// If the system has been homogenized (it has minted homogenizing variables), a freshly-added
		// block is NOT homogenized by the system, so we must homogenize the slice ourselves to match
		// -- folding its constant term onto the homogenizing variable, exactly as System::Homogenize
		// did to the system's other blocks.  Otherwise (an affine system, or a born-projective system
		// whose homogeneous variable group carries no separate homogenizing variable) the slice is
		// added as authored.
		if (s.NumHomVariables() == 0)
		{
			// The block's columns are indexed by the system's variable ordering; if the counts
			// disagree the matrix-vector product is meaningless.  Fail with a clear message rather
			// than a silently-wrong (or out-of-bounds) evaluation downstream.
			if (NumVariables() != s.NumVariables())
			{
				std::stringstream ss;
				ss << "cannot add a slice on " << NumVariables() << " variables to a system on "
				   << s.NumVariables() << " variables; build the slice over the system's variables";
				throw std::runtime_error(ss.str());
			}
			s.AddBlock(block_);
			return;
		}

		// Homogenized system: fold the slice's constant onto the homogenizing variable, then add.
		// LinearFormsBlock::Homogenize supports a single affine variable group, so mirror that limit.
		if (s.NumVariableGroups() != 1)
			throw std::runtime_error("adding a slice to a homogenized system with multiple affine "
			                         "variable groups is not yet supported");
		if (NumVariables() != s.NumNaturalVariables())
		{
			std::stringstream ss;
			ss << "cannot add a slice on " << NumVariables() << " variables to a homogenized system on "
			   << s.NumNaturalVariables() << " natural variables; build the slice over the system's "
			   << "(affine) variables before it was homogenized";
			throw std::runtime_error(ss.str());
		}

		blocks::LinearFormsBlock homogenized_block = block_;
		homogenized_block.Homogenize(s.VariableGroups()[0], s.HomogenizingVariables()[0]);
		s.AddBlock(homogenized_block);
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
