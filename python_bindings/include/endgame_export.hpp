//This file is part of Bertini 2.
//
//python/endgame_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/endgame_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/endgame_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Notre Dame
//  Summer 2016, Summer 2023, Fall 2024
//
//
//  python/endgame_export.hpp:  Header file for exposing endgames to python.

#pragma once

#include "python_common.hpp"

#include <bertini2/endgames.hpp>

namespace bertini{
	namespace python{

		using namespace bertini::tracking;

		/**
		 Abstract Endgame class
		 */
		template<typename EndgameT>
		class EndgameBaseVisitor: public def_visitor<EndgameBaseVisitor<EndgameT> >
		{
			friend class ::boost::python::def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const;

		private:

			using BCT = typename TrackerTraits<typename EndgameT::TrackerType>::BaseComplexT;
			using BRT = typename TrackerTraits<typename EndgameT::TrackerType>::BaseRealT;
			using BaseEGT = typename EndgameT::BaseEGT;

			static
			SuccessCode WrapRunDefaultTime(EndgameT & self, BCT t, Vec<BCT> const& s){
				return self.Run(t, s);
			}

			static
			 SuccessCode WrapRunCustomTime(EndgameT & self, BCT t, Vec<BCT> const& s, BCT const& u){
				return self.Run(t, s, u);
			}


			using unsigned_of_void = unsigned (BaseEGT::*)() const;
			static unsigned_of_void GetCycleNumberFn()
			{
				return &BaseEGT::CycleNumber;
			};

			template <typename T>
			static
			Vec<T> return_final_approximation(EndgameT const& self)
			{
				return self.template FinalApproximation<T>();
			}

		};// EndgameVisitor class



		/**
		 Particulars for the PowerSeries endgame.
		 */
		template<typename PowerSeriesT>
		class PowerSeriesVisitor: public def_visitor<PowerSeriesVisitor<PowerSeriesT> >
		{
			friend class ::boost::python::def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const;

		private:

			using BCT = typename TrackerTraits<typename PowerSeriesT::TrackerType>::BaseComplexT;
			using BRT = typename TrackerTraits<typename PowerSeriesT::TrackerType>::BaseRealT;


		};// CauchyVisitor class




		/**
		 Particulars for the Cauchy endgame.
		 */
		template<typename CauchyT>
		class CauchyVisitor: public def_visitor<CauchyVisitor<CauchyT> >
		{
			friend class ::boost::python::def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const;

		private:

			using BCT = typename TrackerTraits<typename CauchyT::TrackerType>::BaseComplexT;
			using BRT = typename TrackerTraits<typename CauchyT::TrackerType>::BaseRealT;


		};// CauchyVisitor class




		// now prototypes for expose functions defined in the .cpp files for the python bindings.

		/**
		The main function for exporting the bound endgames to Python.

		This should be the only function called from the main function defining the module, and should call all those functions exposing particular endgames.
		*/
		void ExportEndgames();


		/**
		export the power series endgame incarnations
		*/
		void ExportAMPPSEG();
		void ExportFDPSEG();
		void ExportFMPSEG();


		/**
		export the cauchy endgame incarnations
		*/
		void ExportAMPCauchyEG();
		void ExportFDCauchyEG();
		void ExportFMCauchyEG();



}}// re: namespaces

