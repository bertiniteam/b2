


#pragma once
#include <functional>

namespace bertini{

	namespace nag_datatype{

		namespace policy{

			template<typename T>
			struct Copy
			{
				using HeldT = T;


				static
				const T & AtGet(HeldT const& t)
				{
					return t;
				}

			};


			template<typename T>
			struct Reference
			{
				using HeldT = std::reference_wrapper<T>;


				static
				const T & AtGet(HeldT const& t)
				{
					return t;
				}

			};
			

			template<typename T>
			struct SharedPtr
			{
				using HeldT = std::shared_ptr<T>;


				static
				const T & AtGet(HeldT const& t)
				{
					return *t;
				}

			};

		}// policy
	}// nag_datatype
}//bertini