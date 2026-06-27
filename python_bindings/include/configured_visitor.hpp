//This file is part of Bertini 2.
//
//configured_visitor.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//configured_visitor.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with configured_visitor.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Wisconsin-Eau Claire
//
//
//  python/configured_visitor.hpp:  A reusable boost.python def_visitor that, for
//  any object deriving from bertini::detail::Configured<...> (every tracker and
//  nag_algorithm), exposes a uniform Python interface for getting and setting its
//  configuration structs.  It is driven entirely by the object's own
//  `Config::UsedConfigs` type-list, so a newly added algorithm gets this whole
//  interface from a single `.def(ConfiguredVisitor<Algo>())` -- no per-config
//  boilerplate, and nothing about the config list is duplicated in Python.

#pragma once

#include "python_common.hpp"

#include <bertini2/detail/typelist.hpp>

namespace bertini{
	namespace python{

		namespace config_detail{

			// The boost.python class object registered for the C++ type C, or None
			// if C has not (yet) been exposed.  Config structs are always exported
			// before the owners that use them, so this is populated by the time a
			// ConfiguredVisitor runs.
			template<typename C>
			boost::python::object RegisteredClass()
			{
				auto const* reg = boost::python::converter::registry::query(boost::python::type_id<C>());
				if (reg && reg->m_class_object)
					return boost::python::object(boost::python::handle<>(boost::python::borrowed(
						reinterpret_cast<PyObject*>(reg->m_class_object))));
				return boost::python::object(); // None
			}

			// Recursion over the config type-list: return a *copy* of the owner's
			// configuration struct whose registered Python class is `cls`.
			template<typename OwnerT>
			boost::python::object GetConfigDispatch(OwnerT&, boost::python::object const& /*cls*/)
			{
				PyErr_SetString(PyExc_KeyError,
					"this object has no configuration of the requested type; see config_types()");
				boost::python::throw_error_already_set();
				return boost::python::object();
			}

			template<typename OwnerT, typename C, typename... Rest>
			boost::python::object GetConfigDispatch(OwnerT& self, boost::python::object const& cls)
			{
				boost::python::object rc = RegisteredClass<C>();
				if (!rc.is_none() && rc.ptr() == cls.ptr())
					return boost::python::object(self.template Get<C>()); // copy out
				return GetConfigDispatch<OwnerT, Rest...>(self, cls);
			}

			// Bundles the type-list-driven free functions bound onto the owner.
			template<typename OwnerT, typename... Cs>
			struct Dispatch
			{
				static boost::python::object Get(OwnerT& self, boost::python::object cls)
				{
					return GetConfigDispatch<OwnerT, Cs...>(self, cls);
				}

				static boost::python::list Types(OwnerT&)
				{
					boost::python::list out;
					(out.append(RegisteredClass<Cs>()), ...);
					return out;
				}
			};

		} // namespace config_detail


		/**
		 A def_visitor exposing the configuration interface of any
		 bertini::detail::Configured<...> owner (trackers, nag_algorithms):

		   - set_config(cfg)     : store a configuration struct (overloaded per type)
		   - get_config(cls)     : get a copy of the stored config of class `cls`
		   - config_types()      : list the config classes this owner accepts

		 Pythonic conveniences (update/repr/configure/...) are layered on top of
		 these in pure Python (bertini.config).
		 */
		template<typename OwnerT>
		class ConfiguredVisitor : public def_visitor<ConfiguredVisitor<OwnerT> >
		{
			friend class ::boost::python::def_visitor_access;

			template<class PyClass, typename... Cs>
			static void DoVisit(PyClass& cl, bertini::detail::TypeList<Cs...>*)
			{
				// one set_config, overloaded by argument type, per config the owner holds
				(cl.def("set_config",
				        static_cast<void (OwnerT::*)(Cs const&)>(&OwnerT::template Set<Cs>),
				        (arg("self"), arg("config")),
				        "Store one of this object's configuration structs (dispatched by the config's type)."), ...);

				cl.def("get_config", &config_detail::Dispatch<OwnerT, Cs...>::Get,
				       (arg("self"), arg("config_type")),
				       "Return a copy of this object's stored configuration struct of the given class.");

				cl.def("config_types", &config_detail::Dispatch<OwnerT, Cs...>::Types,
				       (arg("self")),
				       "List the configuration struct classes this object accepts.");
			}

		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				DoVisit(cl, static_cast<typename OwnerT::Config::UsedConfigs*>(nullptr));
			}
		}; // ConfiguredVisitor


	}} // namespaces
