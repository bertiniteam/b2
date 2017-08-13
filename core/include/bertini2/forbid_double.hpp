

#pragma once

namespace boost{
	namespace multiprecision
	{
template <class Backend>
struct is_compatible_arithmetic_type<double, number<Backend> > : public mpl::false_ {};


	}
}

