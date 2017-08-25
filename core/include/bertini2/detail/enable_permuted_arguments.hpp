//This file is part of Bertini 2.
//
//enable_permuted_arguments.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//enable_permuted_arguments.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with enable_permuted_arguments.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
// additionally, this file contains code
// based on and refined from SO post 
// https://stackoverflow.com/questions/19329297/constructor-permutations-for-passing-parameters-in-arbitrary-order
// by user Daniel Frey.

/**
\file enable_permuted_arguments.hpp

Contains helper templates, based on and refined from SO post 
https://stackoverflow.com/questions/19329297/constructor-permutations-for-passing-parameters-in-arbitrary-order
by user Daniel Frey.
*/

#ifndef BERTINI_ENABLE_PERMUTED_ARGUMENTS
#define BERTINI_ENABLE_PERMUTED_ARGUMENTS

#include <type_traits>
#include <tuple>


namespace bertini {

	
	namespace detail {
	// find the index of a type in a list of types,
	// return sizeof...(Remaining) if T is not found
	template< typename T, typename... Remaining >
	struct IndexByType : std::integral_constant< std::size_t, 0 > 
	{};

	// template specialization, enforcing the equality of the number of variadic arguments, and the corresponding index.
	// if two types are the same, this template will be preferred since it's a closer match, and the static_assert will fail, preventing compilation.
	template< typename T, typename... RemainingT >
	struct IndexByType< T, T, RemainingT... > : std::integral_constant< std::size_t, 0 >
	{
	   static_assert( IndexByType< T, RemainingT... >::value == sizeof...( RemainingT ), "No duplicates are allowed when enabling permuted arguments" );
	};

	// template specialization on two types.  looks for the index of T in the latter arguments.  If T1==T2, this template will not be used, but instead the above will.
	template< typename T1, typename T2, typename... RemainingT >
	struct IndexByType< T1, T2, RemainingT... > : std::integral_constant< std::size_t, IndexByType< T1, RemainingT... >::value + 1 > 
	{};

	// a free function which gets a type from a tuple with no duplicates.  Remember, duplicates are not allowed.
	template< std::size_t I, std::size_t J, 
			typename T, typename... PresentT, typename... GivenT >
	auto GetByIndex( const std::tuple< PresentT... >&, 
	                 const std::tuple< GivenT... >& absent_args )
	   -> typename std::enable_if< I == sizeof...( PresentT ), const T& >::type
	{
	   return std::get< J >( absent_args );
	}

	// the other half of the get by index function
	template< std::size_t I, std::size_t J, 
	        typename T, typename... PresentT, typename... GivenT >
	auto GetByIndex(const std::tuple< PresentT... >& present_args, 
	                const std::tuple< GivenT... >& ) // these arguments will have to be default constructed
	   -> typename std::enable_if< I != sizeof...( PresentT ), const T& >::type /* remember, `I != sizeof...( PresentT )` means it found T in PresentT */
	{
	   return std::get< I >( present_args );
	}




	// helper to validate that all Us are in Ts...
	template< bool > struct AreValidArguments;

	// the false one is never instantiated, so the 'false' template parameter will fail, deliberately
	template<> struct AreValidArguments< true > : std::true_type 
	{};

	// a proxy function which to ensure that the types are distinct, and that we can indeed get the index by type.
	template< std::size_t... > 
	void ValidateTypes() 
	{}



	// default construct the objects not present in the arguments, but required by the unpermute function
	template< typename T >
	struct DefaultConstruct
	{ static const T value; };

	template< typename T >
	const T DefaultConstruct< T >::value
	{};



	} // namespace detail

	// helper template which reorders parameters
	template< typename... UnpermutedT, 	// these must be declared when using the function
	          typename... PermutedPresentT > // these are inferred, cannot be declared
	std::tuple< const UnpermutedT&... > Unpermute( const PermutedPresentT&... present_permuted_args )
	{
		using namespace detail;
	   auto present_args = std::tie( present_permuted_args... );
	   auto absent_args  = std::tie( DefaultConstruct< UnpermutedT >::value... );

	   //next we validate that the input arguments are valid for unpermuting.  must have all distinct types.
	   ValidateTypes< AreValidArguments< IndexByType< const PermutedPresentT&, 
	                                            const UnpermutedT&... >::value != sizeof...( UnpermutedT ) 
	                                            >::value... 
	                 >(); // this will fail if AreValidArguments are not valid ('true')

	   // finally, we return a tie of the unpermuted types, default constructed if necessary
	   return std::tie( GetByIndex< IndexByType< const UnpermutedT&, const   PermutedPresentT&... >::value,
	                                IndexByType< const UnpermutedT&, const UnpermutedT&... >::value, UnpermutedT >
	                                ( present_args, absent_args )... );
	}


} // re: namespace bertini

#endif
