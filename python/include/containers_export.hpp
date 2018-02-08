//This file is part of Bertini 2.
//
//python/containers_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/containers_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/containers_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016-2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//  Danielle Brake
//  UWEC
//  Spring 2018
//
//  python/containers_export.hpp:  Exports all needed containers from Bertini 2.0 to python.

#pragma once
#ifndef BERTINI_PYTHON_CONTAINERS_EXPORT_HPP
#define BERTINI_PYTHON_CONTAINERS_EXPORT_HPP


#include <bertini2/function_tree.hpp>

#include "python_common.hpp"


namespace bertini{ namespace python{

template< typename T>
inline std::ostream& operator<<(std::ostream & out, const std::vector<T> & t)
{
	out << "[";
	for (int ii = 0; ii < t.size(); ++ii)
	{
		out << t[ii];
		if (ii!=t.size()-1)
		{
			out << ", ";
		}
	}
	out << "]";

	return out;
}


template< typename T>
inline std::ostream& operator<<(std::ostream & out, const std::deque<T> & t)
{
	out << "[";
	for (int ii = 0; ii < t.size(); ++ii)
	{
		out << t[ii];
		if (ii!=t.size()-1)
		{
			out << ", ";
		}
	}
	out << "]";

	return out;
}


/**
Adds functionality to iterable types
 */
template<typename ContT>
class ListVisitor: public def_visitor<ListVisitor<ContT> >
{
	friend class def_visitor_access;
	
public:
	template<class PyClass>
	void visit(PyClass& cl) const;
	
private:


	static std::string __str__(const object& obj)
	{
		std::ostringstream oss;
		const ContT& self=extract<ContT>(obj)();
		std::stringstream ss;
		ss << "[";
		for (int ii = 0; ii < self.size(); ++ii)
		{
			ss << self[ii];
			if (ii!=self.size()-1)
			{
				ss << ", ";
			}
		}
		ss << "]";
		return ss.str();
	};

	static std::string __repr__(const object& obj)
	{
		return __str__(obj);
		// std::ostringstream oss;
		// const ContT& self=extract<ContT>(obj)();
		// std::stringstream ss;
		// ss << self.str(0,std::ios::scientific);
		// return ss.str();
	};	

};// ListVisitor class



void ExportContainers();

}} // namespaces
#endif
