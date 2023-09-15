//This file is part of Bertini 2.
//
//trig.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//trig.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with trig.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2015
//
// silviana amethyst, university of wisconsin-eau claire
//
//
//
// trig.hpp:  Declares the trigonometric operator classes.

/**
\file trig.hpp

\brief Provides the abstract TrigOperator, and the concrete types such as SinOperator.

*/

#ifndef BERTINI_TRIGONOMETRIC_OPERATOR_HPP
#define BERTINI_TRIGONOMETRIC_OPERATOR_HPP

#include "bertini2/function_tree/operators/operator.hpp"
#include "bertini2/function_tree/operators/arithmetic.hpp"



namespace bertini {
namespace node{	
	/**
	\brief Abstract class for trigonometric Operator types.

	Abstract class for trigonometric Operator types.
	*/
	class TrigOperator: public virtual UnaryOperator
	{
	public:
		BERTINI_DEFAULT_VISITABLE()
		
		TrigOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};

		virtual ~TrigOperator() = default;

		/**
		 Compute the degree with respect to a single variable.
		 
		 For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
	protected:
		TrigOperator(){}
	private:
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<UnaryOperator>(*this);
		}


	};





	/**
	\brief Provides the sine Operator.

	This class represents the sine function.  FreshEval method
	is defined for sine and takes the sine of the child node.
	*/
	class SinOperator : public virtual TrigOperator, public virtual EnableSharedFromThisVirtual<SinOperator>
	{
	public:
		BERTINI_DEFAULT_VISITABLE()

		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;

		template<typename... Ts> 
		static 
		std::shared_ptr<SinOperator> Make(Ts&& ...ts){ 
			return std::shared_ptr<SinOperator>( new SinOperator(ts...) );
		}

	private:
		SinOperator(const std::shared_ptr<Node> & N) : TrigOperator(N), UnaryOperator(N)
		{};
		
	public:
		
		
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the sine function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;

		virtual ~SinOperator() = default;
		
	protected:
		
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		mpfr_complex FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		
		
		void FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		
	private:
		SinOperator() = default;
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};
	

	/**
	\brief Provides the inverse sine Operator.

	This class represents the inverse sine function.  FreshEval method
	is defined for arcsine and takes the sine of the child node.
	*/
	class ArcSinOperator : public  virtual TrigOperator, public virtual EnableSharedFromThisVirtual<ArcSinOperator>
	{
	public:
		BERTINI_DEFAULT_VISITABLE()

		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;
		
		template<typename... Ts> 
		static 
		std::shared_ptr<ArcSinOperator> Make(Ts&& ...ts){ 
			return std::shared_ptr<ArcSinOperator>( new ArcSinOperator(ts...) );
		}

	private:
		ArcSinOperator(const std::shared_ptr<Node> & N) : TrigOperator(N), UnaryOperator(N)
		{};

	public:
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the sine function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;

		virtual ~ArcSinOperator() = default;
		
	protected:
		
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		mpfr_complex FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;

		
		void FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

	private:
		ArcSinOperator() = default;
		friend class boost::serialization::access;
		

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};




	
	
	
	
	/**
	\brief  Provides the cosine Operator.

	This class represents the cosine function.  FreshEval method
	is defined for cosine and takes the cosine of the child node.
	*/
	class CosOperator : public  virtual TrigOperator, public virtual EnableSharedFromThisVirtual<CosOperator>
	{
	public:
		BERTINI_DEFAULT_VISITABLE()
		
		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;

		template<typename... Ts> 
		static 
		std::shared_ptr<CosOperator> Make(Ts&& ...ts){ 
			return std::shared_ptr<CosOperator>( new CosOperator(ts...) );
		}

	private:
		CosOperator(const std::shared_ptr<Node> & N) : TrigOperator(N), UnaryOperator(N)
		{};
		
	public:
		
		
		
		void print(std::ostream & target) const override;
		
		
		
		
		/**
		 Differentiates the cosine function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;

		virtual ~CosOperator() = default;
		
	protected:
		
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		mpfr_complex FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		
		
		void FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		
		
	private:
		CosOperator() = default;
		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};
	
	
	/**
	\brief Provides the arc cosine Operator.

	This class represents the inverse cosine function.  FreshEval method
	is defined for arccosine and takes the arccosine of the child node.
	*/
	class ArcCosOperator : public  virtual TrigOperator, public virtual EnableSharedFromThisVirtual<ArcCosOperator>
	{
	public:
		BERTINI_DEFAULT_VISITABLE()

		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;
		
		template<typename... Ts> 
		static 
		std::shared_ptr<ArcCosOperator> Make(Ts&& ...ts){ 
			return std::shared_ptr<ArcCosOperator>( new ArcCosOperator(ts...) );
		}

	private:
		ArcCosOperator(const std::shared_ptr<Node> & N) : TrigOperator(N), UnaryOperator(N)
		{};
	public:
		
		
		
		void print(std::ostream & target) const override;
		
		
		
		/**
		 Differentiates the cosine function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		virtual ~ArcCosOperator() = default;
		
	protected:
		
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		mpfr_complex FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		
		
		void FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		
	private:
		ArcCosOperator() = default;
		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};
	
	
	
	
	


	
	
	
	/**
	\brief Provides the tangent Operator.

	This class represents the tangent function.  FreshEval method
	is defined for tangent and takes the tangent of the child node.
	*/
	class TanOperator : public  virtual TrigOperator, public virtual EnableSharedFromThisVirtual<TanOperator>
	{
	public:
		BERTINI_DEFAULT_VISITABLE()
		
		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;
		
		template<typename... Ts> 
		static 
		std::shared_ptr<TanOperator> Make(Ts&& ...ts){ 
			return std::shared_ptr<TanOperator>( new TanOperator(ts...) );
		}

	private:
		TanOperator(const std::shared_ptr<Node> & N) : TrigOperator(N), UnaryOperator(N)
		{};
	public:
		
		
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the tangent function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		virtual ~TanOperator() = default;
		
	protected:
		
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		mpfr_complex FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		
		
		void FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		
	private:
		TanOperator() = default;
		friend class boost::serialization::access;
		

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};
	
	
	/**
	\brief Provides the inverse tangent Operator.

	This class represents the inverse tangent function.  FreshEval method
	is defined for arctangent and takes the arc tangent of the child node.
	*/
	class ArcTanOperator : public  virtual TrigOperator, public virtual EnableSharedFromThisVirtual<ArcTanOperator>
	{
	public:
		BERTINI_DEFAULT_VISITABLE()
		
		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;
		
		template<typename... Ts> 
		static 
		std::shared_ptr<ArcTanOperator> Make(Ts&& ...ts){ 
			return std::shared_ptr<ArcTanOperator>( new ArcTanOperator(ts...) );
		}

	private:
		ArcTanOperator(const std::shared_ptr<Node> & N) : TrigOperator(N), UnaryOperator(N)
		{};
	public:
		
		
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the tangent function.
		 */		
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		virtual ~ArcTanOperator() = default;
		
	protected:
		
		
		// Specific implementation of FreshEval for arctangent.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		mpfr_complex FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		
		
		void FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;
		
	private:
		ArcTanOperator() = default;
		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};
	
	
	// begin the overload of operators

	inline std::shared_ptr<Node> sin(const std::shared_ptr<Node> & N)
	{
		return SinOperator::Make(N);
	}
	
	inline std::shared_ptr<Node> asin(const std::shared_ptr<Node> & N)
	{
		return ArcSinOperator::Make(N);
	}
	


	inline std::shared_ptr<Node> cos(const std::shared_ptr<Node> & N)
	{
		return CosOperator::Make(N);
	}



	inline std::shared_ptr<Node> acos(const std::shared_ptr<Node> & N)
	{
		return ArcCosOperator::Make(N);
	}



	inline std::shared_ptr<Node> tan(const std::shared_ptr<Node> & N)
	{
		return TanOperator::Make(N);
	}


	
	inline std::shared_ptr<Node> atan(const std::shared_ptr<Node> & N)
	{
		return ArcTanOperator::Make(N);
	}
	
	
} // re: namespace node
} // re: namespace bertini


#endif
