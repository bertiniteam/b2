//This file is part of Bertini 2.
//
//bertini2/io/generators.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/generators.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/generators.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/io/generators.hpp 

\brief Provides the generators for bertini2.
*/

#pragma once


#include "bertini2/mpfr_complex.hpp"
#include "bertini2/eigen_extensions.hpp"

#include <boost/math/special_functions/modf.hpp> 

BOOST_MATH_STD_USING

#include <boost/fusion/adapted.hpp>
#include <boost/fusion/include/adapted.hpp>

#define BOOST_SPIRIT_USE_PHOENIX_V3 1

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/fusion/include/std_pair.hpp>


#include <boost/spirit/include/support_istream_iterator.hpp>



// BOOST_FUSION_ADAPT_ADT(
//     std::complex<double>,
//     (bool, bool, obj.imag() != 0, /**/)
//     (double, double, obj.real(), /**/)
//     (double, double, obj.imag(), /**/)
// )


// BOOST_FUSION_ADAPT_ADT(
//     bertini::complex,
//     (bool, bool, obj.imag() != 0, /**/)
//     (bertini::mpfr_float, bertini::mpfr_float, obj.real(), /**/)
//     (bertini::mpfr_float, bertini::mpfr_float, obj.imag(), /**/)
// )


namespace bertini{
	namespace generators{

		namespace karma = boost::spirit::karma;

		template <typename Num>
		struct BertiniNumPolicy : public karma::real_policies<Num>
		{
		    // we want the numbers always to be in scientific format
		    static int floatfield(Num n) { return std::ios_base::scientific; }

		    static unsigned int precision(Num) {
		        return std::numeric_limits<Num>::max_digits10;
		      }
		};

		template<>
		struct BertiniNumPolicy<mpfr_float>  : public karma::real_policies<mpfr_float>
		{
		    // we want the numbers always to be in scientific format
		    static int floatfield(mpfr_float n) { return std::ios_base::scientific; }

		    static unsigned int precision(mpfr_float const& x) {
		        return x.precision();
		      }
		};



		template<typename Num>
		using FullPrec = boost::spirit::karma::real_generator<Num, BertiniNumPolicy<Num> >;

		FullPrec<double> const full_prec_d = FullPrec<double>();
		FullPrec<mpfr_float> const full_prec_mp = FullPrec<mpfr_float>();


		struct Classic{

			template <typename OutputIterator>
			static bool generate(OutputIterator sink, double const& c)
			{
	            using boost::spirit::karma::omit;
	            using boost::spirit::karma::generate;

	            return generate(sink,

	                //  Begin grammar
	                (
	                   full_prec_d << " " << full_prec_d
	                ),
	                //  End grammar

	                c, 0     //  Data to output
	                );
			}


			template <typename OutputIterator>
			static bool generate(OutputIterator sink, std::complex<double> const& c)
			{
	            using boost::spirit::karma::omit;
	            using boost::spirit::karma::generate;

	            return generate(sink,

	                //  Begin grammar
	                (
	                   full_prec_d << " " << full_prec_d
	                ),
	                //  End grammar

	                c.real(), c.imag()     //  Data to output
	                );
			}



			template <typename OutputIterator>
			static bool generate(OutputIterator sink, mpfr_float const& c)
			{
	            using boost::spirit::karma::omit;
	            using boost::spirit::karma::generate;

	            return generate(sink,

	                //  Begin grammar
	                (
	                   boost::spirit::karma::string
	                ),
	                //  End grammar

	                c.str()     //  Data to output
	                );
			}


			template <typename OutputIterator>
			static bool generate(OutputIterator sink, mpfr const& c)
			{
	            using boost::spirit::karma::omit;
	            using boost::spirit::karma::generate;

	            return generate(sink,

	                //  Begin grammar
	                (
	                   boost::spirit::karma::string << " " << boost::spirit::karma::string
	                ),
	                //  End grammar

	                c.real().str(), c.imag().str()     //  Data to output
	                );
			}


			template <typename OutputIterator, typename T>
			static bool generate(OutputIterator sink, Vec<T> const& c)
			{
				auto s = c.size();
				for (decltype(s) ii=0; ii<s; ++ii)
				{
					bool b = generate(sink, c(ii));
					if (!b)
						return b;
					boost::spirit::karma::generate(sink, (boost::spirit::karma::char_), '\n');
				}
				return true;
			}

		};



		struct CPlusPlus{

			template <typename OutputIterator>
			static bool generate(OutputIterator sink, std::complex<double> const& c)
			{
	            using boost::spirit::karma::omit;
	            using boost::spirit::karma::generate;

	            return generate(sink,

	                //  Begin grammar
	                (
	                   !full_prec_d(0.0) << '(' << full_prec_d << ", " << full_prec_d << ')'
	                |   omit[full_prec_d] << full_prec_d
	                ),
	                //  End grammar

	                c.imag(), c.real(), c.imag()     //  Data to output
	                );
			}
		};



	}
}

