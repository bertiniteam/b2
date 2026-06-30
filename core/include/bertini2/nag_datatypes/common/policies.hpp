


#pragma once
#include <functional>

namespace bertini{

	namespace nag_datatype{

		namespace policy{

			/// \brief Storage policy that holds objects by value (a copy).
			template<typename T>
			struct Copy
			{
				using HeldT = T;  ///< The held type (the object itself).


				/// \brief Access the held object.
				static
				const T & AtGet(HeldT const& t)
				{
					return t;
				}

			};


			/// \brief Storage policy that holds a reference to an externally-owned object.
			template<typename T>
			struct Reference
			{
				using HeldT = std::reference_wrapper<T>;  ///< The held type (a reference wrapper).


				/// \brief Access the referenced object.
				static
				const T & AtGet(HeldT const& t)
				{
					return t;
				}

			};


			/// \brief Storage policy that holds objects via std::shared_ptr.
			template<typename T>
			struct SharedPtr
			{
				using HeldT = std::shared_ptr<T>;  ///< The held type (a shared pointer).


				/// \brief Access the pointed-to object.
				static
				const T & AtGet(HeldT const& t)
				{
					return *t;
				}

			};

		}// policy
	}// nag_datatype
}//bertini