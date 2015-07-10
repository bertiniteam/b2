
#include "system.hpp"



namespace bertini 
{
	namespace start_system{

		
		template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
		template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
			

		class StartSystem : public System
		{

		public:
			

			virtual size_t NumStartPoints() const = 0;

			template<typename T>
			Vec<T> StartPoint(size_t index) const
			{
				return GenerateStartPoint(T(),index);
			}

		private:
			virtual Vec<dbl> GenerateStartPoint(dbl,size_t index) const = 0;
			virtual Vec<mpfr> GenerateStartPoint(mpfr,size_t index) const = 0;

		};


		class TotalDegree : public StartSystem
		{
		public:
			TotalDegree() = default;
			virtual ~TotalDegree() = default;
			
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

			size_t NumStartPoints() const override;


		private:


			Vec<dbl> GenerateStartPoint(dbl,size_t index) const override;
			Vec<mpfr> GenerateStartPoint(mpfr,size_t index) const override;

			Vec<mpfr> random_values_; ///< stores the random values for the start functions.  x^d-r, where r is stored in this vector.
			std::vector<int> degrees_;
		};
	}
}
