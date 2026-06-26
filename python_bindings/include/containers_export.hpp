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
// Copyright(C) Bertini2 Development Team
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
//  silviana amethyst
//  UWEC
//  Spring 2018
//
//  python/containers_export.hpp:  Exports all needed containers from Bertini 2.0 to python.

#pragma once
#ifndef BERTINI_PYTHON_CONTAINERS_EXPORT_HPP
#define BERTINI_PYTHON_CONTAINERS_EXPORT_HPP

#include <deque>
#include "python_common.hpp"

#include <bertini2/nag_algorithms/zero_dim_solve.hpp>
#include <bertini2/function_tree.hpp>
#include <boost/python/stl_iterator.hpp>




namespace bertini{ namespace python{

template< typename T>
inline std::ostream& operator<<(std::ostream & out, const std::vector<T> & t)
{
	out << "[";
	for (size_t ii = 0; ii < t.size(); ++ii)
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
	for (size_t ii = 0; ii < t.size(); ++ii)
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
	friend class ::boost::python::def_visitor_access;
	
public:
	template<class PyClass>
	void visit(PyClass& cl) const;
	
private:


	static std::string __str__(const object& obj)
	{
		std::ostringstream oss;
		ContT self=extract<ContT>(obj)();
		std::stringstream ss;
		ss << "[";
		for (size_t ii = 0; ii < self.size(); ++ii)
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






// This block of code lets us construct a container in C++ from a list of things in Python.
// i found this problem difficult.  
//
// fortunately, there were a number of questions and answers of varying quality about it, and the below
// worked readily.
//
// derived from https://stackoverflow.com/questions/56290774/boost-python-exposing-c-class-with-constructor-taking-a-stdlist

template<typename ContT>
std::shared_ptr<ContT> create_MyClass(boost::python::object const& iterable)
{
	using ContainedT = typename ContT::value_type;

    // Accept ANY Python iterable, not just a list -- so e.g. VariableGroup(linalg.variable_vector('x', 3))
    // works directly on the numpy object array, with no list() wrapping.
    ContT temp{ boost::python::stl_input_iterator<ContainedT>(iterable)
        , boost::python::stl_input_iterator<ContainedT>() };
    return std::make_shared<ContT>(temp);
}


template<typename ContT>
struct std_list_to_python
{
    static PyObject* convert(ContT const& l)
    {
        boost::python::list result;
        for (auto const& value : l) {
            result.append(value);
        }
        return boost::python::incref(result.ptr());
    }
};


template<typename ContT>
struct pylist_converter
{	
	using ContainedT = typename ContT::value_type;

    static void* convertible(PyObject* object)
    {
        // Accept any Python iterable (list, tuple, numpy object array, ...) -- so e.g. a
        // bertini.linalg.variable_vector converts to a VariableGroup with no list() wrapping.
        // Strings/bytes are iterable but are never a container of our element type, so reject them.
        if (object == Py_None || PyUnicode_Check(object) || PyBytes_Check(object))
            return nullptr;

        PyObject* iter = PyObject_GetIter(object);
        if (!iter) {
            PyErr_Clear();
            return nullptr;
        }
        Py_DECREF(iter);
        return object;
    }

    static void construct(PyObject* object, boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        typedef boost::python::converter::rvalue_from_python_storage<ContT> storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

        data->convertible = new (storage) ContT();

        ContT* l = (ContT*)(storage);

        boost::python::object iterable(boost::python::handle<>(boost::python::borrowed(object)));
        boost::python::stl_input_iterator<ContainedT> it(iterable), end;
        for (; it != end; ++it)
            l->push_back(*it);
    }
};






void ExportContainers();

}} // namespaces
#endif
