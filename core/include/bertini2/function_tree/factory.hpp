//This file is part of Bertini 2.
//
//include/bertini2/function_tree/factory.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//include/bertini2/function_tree/factory.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with include/bertini2/function_tree/factory.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire


/**
\file include/bertini2/function_tree/factory.hpp

\brief Methods for producing nodes of various types. 

These wrapper functions are intended to be re-written to use something like a workspace, sometime in the future.  This will allow one to refer to things like the `1` node, or the `pi` node, or the `f2` node.  Other desired functionality is things like inspecting all declared symbols or variables, or having multiple workspaces.

*/

#pragma once

#include <memory>
#include <vector>
#include "bertini2/function_tree/forward_declares.hpp"

namespace bertini {
	
		// template <typename T>
		// struct Factory{

		// 	template <typename ... Args>
		// 	static
		// 	std::shared_ptr<T> Make (Args&& ... args)
		// 	{
		// 		return T::Make(args...);
		// 	}

		// };




	class Factory
	{


	public:
		using HeldType = std::shared_ptr<node::Node>;

		// static
		// HeldType NonPtrAdd(ObjT d)
		// {
		// 	// held_data_.push_back(PointerPolicy::FromObj(d));
		// 	held_data_.push_back(PointerPolicy::FromObj(d));
		// 	return held_data_.back();
		// }

		// static
		// HeldType PtrAdd(HeldType d)
		// {
		// 	held_data_.push_back(d);
		// 	return held_data_.back();
		// }

		template<typename NodeT, typename ... Args>
		static
		std::shared_ptr<NodeT> Make(Args&& ...args )
		{
			std::shared_ptr<NodeT> ptr = NodeT::Make(args...);
			held_data_.push_back(ptr);
			return ptr;
		}

		template<typename NodeT, typename ... Args>
		static
		std::shared_ptr<NodeT> MakeInPlace(NodeT * raw_ptr, Args&& ...args )
		{
			std::shared_ptr<NodeT> ptr = NodeT::MakeInPlace(raw_ptr, args...);
			held_data_.push_back(ptr);
			return ptr;
		}


		static
		void PurgeCache()
		{
			std::remove_if(held_data_.begin(), held_data_.end(), [](HeldType const& h){return h.use_count()==1;});
		}

	private:
		static
		std::vector<HeldType> held_data_;
	};

	



	// /**
	// \brief Make a variable.

	// \return A node for you!

	// \param name The name of the variable to make.  In principle can be anything, but stick to Bertini naming rules.  Don't be silly, you'll cause serialization-deserialization round tripping problems.

	// \return A shared pointer to a variable, the name of which you stated.
	// */
	// template<typename T>
	// std::shared_ptr<node::Variable> MakeVariable(T const& t)
	// {
	template <typename ... T>
	std::shared_ptr<node::Variable> MakeVariable(T&& ...t)
	{
		return Factory::Make<node::Variable>(t...);
	}
	// }


	// /**
	// \brief Make an integer.

	// \return A node for you!

	// \param name The integer to make.
	// */
	// template<typename T>
	// std::shared_ptr<node::Integer> MakeInteger(T const& t)
	// {
	template <typename ... T>
	std::shared_ptr<node::Integer> MakeInteger(T&& ... t)
	{
		return Factory::Make<node::Integer>(t...);
	}
	// }


	// template<typename ... T>
	// std::shared_ptr<node::Rational> MakeRational(T const& ...t)
	// {
	template <typename ... T>
	std::shared_ptr<node::Rational> MakeRational(T&& ... t)
	{
		return Factory::Make<node::Rational>(t...);
	}
	// }


	// template<typename ... T>
	// std::shared_ptr<node::Float> MakeFloat(T const& ...t)
	// {
	template <typename ... T>
	std::shared_ptr<node::Float> MakeFloat(T&& ... t)
	{
		return Factory::Make<node::Float>(t...);
	}
	// }

	// template<typename ... T>
	// std::shared_ptr<node::LinearProduct> MakeLinearProduct(T const& ...t)
	// {
	template <typename ... T>
	std::shared_ptr<node::LinearProduct> MakeLinearProduct(T&& ... t)
	{
		return Factory::Make<node::LinearProduct>(t...);
	}
	// }

	// template<typename ... T>
	// std::shared_ptr<node::Function> MakeFunction(T const& ...t)
	// {
	template <typename ... T>
	std::shared_ptr<node::Function> MakeFunction(T&& ... t)
	{
		return Factory::Make<node::Function>(t...);
	}
	// }

	// template<typename ... T>
	// std::shared_ptr<node::Differential> MakeDifferential(T const& ...t)
	// {
	template <typename ... T>
	std::shared_ptr<node::Differential> MakeDifferential(T&& ... t)
	{
		return Factory::Make<node::Differential>(t...);
	}
	// }

	// template<typename ... T>
	// std::shared_ptr<node::Jacobian> MakeJacobian(T const& ...t)
	// {
	template <typename ... T>
	std::shared_ptr<node::Jacobian> MakeJacobian(T&& ... t)
	{
		return Factory::Make<node::Jacobian>(t...);
	}

	template <typename ... T>
	std::shared_ptr<node::DiffLinear> MakeDiffLinear(T&& ... t)
	{
		return Factory::Make<node::DiffLinear>(t...);
	}


	template <typename ... T>
	std::shared_ptr<node::SumOperator> MakeSumOperator(T&& ... t)
	{
		return Factory::Make<node::SumOperator>(t...);
	}

	template <typename ... T>
	std::shared_ptr<node::MultOperator> MakeMultOperator(T&& ... t)
	{
		return Factory::Make<node::MultOperator>(t...);
	}
	// }

}// namespace bertini
