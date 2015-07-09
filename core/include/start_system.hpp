
#include "system.hpp"



namespace bertini 
{
	namespace start_system{
		class TotalDegree : public System
		{
		public:
			TotalDegree() = default;
			~TotalDegree() = default;
			
			/**
			 Constructor for making a total degree start system from a polynomial system
			*/
			TotalDegree(System const& s);


			mpfr RandomValue(size_t index)
			{
				return random_values_(index);
			}


			Vec<mpfr> const& RandomValues()
			{
				return random_values_;
			}


		private:

			Vec<mpfr> random_values_; ///< stores the random values for the start functions.  x^d-r, where r is stored in this vector.

		};
	}
}
