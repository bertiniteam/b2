//This file is part of Bertini 2.
//
//block.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//block.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with block.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file bertini2/system/blocks/block.hpp

\brief The evaluation-block contract for the block-composed System.

A System is composed of evaluation blocks held in a std::variant and dispatched by
std::visit -- closed set, no type erasure.  Each block contributes a contiguous range
of rows to the system's function vector / Jacobian.  Every block must provide, for each
numeric type T in {dbl, mpfr_complex}:

    template<typename T> void EvalInPlace(Eigen::Ref<Vec<T>> seg, Vec<T> const& vars, T const& t) const;
    template<typename T> void JacobianInPlace(Eigen::Ref<Mat<T>> blk, Vec<T> const& vars, T const& t) const;
    template<typename T> void TimeDerivInPlace(Eigen::Ref<Vec<T>> seg, Vec<T> const& vars, T const& t) const;

plus the non-templated metadata / precision surface:

    size_t  NumFunctions() const;          // rows this block contributes
    bool    DependsOnPathVariable() const; // false => TimeDerivInPlace writes zeros, and dH/dt can skip it
    unsigned Precision() const;
    void     Precision(unsigned) const;

`seg`/`blk` are caller-provided sub-ranges of the system's storage (the block writes its
own rows; the caller owns segmentation, the patch, etc.).  `vars` is the full system
variable vector; `t` is the path-variable value (ignored by blocks that don't depend on
it).  The contract is value-in / value-out so blocks are self-contained and testable
without a System.

We are on C++17, so the contract is enforced structurally (std::visit + the is_block
detection trait below + static_assert), not via a C++20 concept.
*/

#pragma once

#include <type_traits>
#include <variant>

#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"

#include "bertini2/system/blocks/products_of_linears_block.hpp"
#include "bertini2/system/blocks/blend_block.hpp"

namespace bertini {
namespace blocks {

/**
\brief Detection trait: does B satisfy the evaluation-block contract?

Checks the dbl instantiation of the templated eval entry points plus the metadata
methods.  Adding a non-conforming type to the block variant then fails a one-line
static_assert instead of producing a deep template error at the visit site.
*/
template <typename, typename = void>
struct is_block : std::false_type {};

template <typename B>
struct is_block<B, std::void_t<
	decltype(std::declval<const B&>().NumFunctions()),
	decltype(std::declval<const B&>().DependsOnPathVariable()),
	decltype(std::declval<const B&>().Precision()),
	decltype(std::declval<const B&>().template EvalInPlace<dbl>(
		std::declval<Eigen::Ref<Vec<dbl>>>(), std::declval<Vec<dbl> const&>(), std::declval<dbl const&>())),
	decltype(std::declval<const B&>().template JacobianInPlace<dbl>(
		std::declval<Eigen::Ref<Mat<dbl>>>(), std::declval<Vec<dbl> const&>(), std::declval<dbl const&>())),
	decltype(std::declval<const B&>().template TimeDerivInPlace<dbl>(
		std::declval<Eigen::Ref<Vec<dbl>>>(), std::declval<Vec<dbl> const&>(), std::declval<dbl const&>()))
>> : std::true_type {};

template <typename B>
inline constexpr bool is_block_v = is_block<B>::value;

// As blocks are added, extend this static_assert list (and, later, the Block variant).
static_assert(is_block_v<ProductsOfLinearsBlock>,
              "ProductsOfLinearsBlock must satisfy the evaluation-block contract");

} // namespace blocks

// Forward declaration: BlendBlock holds its operands by shared_ptr<const System>, and a
// System contains Blocks -- the template parameter keeps System a dependent name so the
// recursive type closes without System being complete here.
class System;

/// The closed set of evaluation blocks a System can be composed of.  Grows as block
/// types are added (eventually a PolynomialBlock so the polynomial part is uniform too).
using Block = std::variant<blocks::ProductsOfLinearsBlock, blocks::BlendBlock<System>>;

} // namespace bertini
