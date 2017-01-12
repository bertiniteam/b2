


#pragma once
#include <functional>

namespace bertini{

	namespace nag_datatype{

		namespace policy{

			template<typename T>
			struct RefToGiven
			{
				using HeldT = std::reference_wrapper<T>;

				static 
				HeldT AtSet(T const& t)
				{
					return std::ref(t);
				}

				static
				const T & AtGet(HeldT const& t)
				{
					return t;
				}

			};

			template<typename T>
			struct CopyGiven
			{
				using HeldT = T;

				static 
				HeldT AtSet(T const& t)
				{
					return t;
				}

				static
				const T & AtGet(HeldT const& t)
				{
					return t;
				}

			};

		}// policy
	}// nag_datatype
}//bertini