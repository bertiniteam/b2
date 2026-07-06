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
// Copyright(C) Bertini2 Development Team
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
	class TrigOperator: public UnaryOperator
	{
	public:
		BERTINI_DEFAULT_VISITABLE()
		
		/// \brief Construct a trigonometric operator over a single operand node.
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
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<UnaryOperator>(*this);
		}


	};





	/**
	\brief Provides the sine Operator.

	This class represents the sine function.
	*/
	class SinOperator : public TrigOperator
	{
	public:
		BERTINI_DEFAULT_VISITABLE()

		std::shared_ptr<Node> Simplified() const override;
		std::shared_ptr<Node> Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const override;

		/// \brief Construct (and intern) a SinOperator node.
		template<typename... Ts> 
		static 
		std::shared_ptr<SinOperator> Make(Ts&& ...ts){ 
			return std::static_pointer_cast<SinOperator>(Intern(std::shared_ptr<Node>( new SinOperator(ts...) )));
		}

	private:
		SinOperator(const std::shared_ptr<Node> & N) : TrigOperator(N)
		{};
		
	public:
		
		
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the sine function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;

		virtual ~SinOperator() = default;
		
	protected:
		
		
		
		
		
		

		
	private:
		SinOperator() = default;
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};
	

	/**
	\brief Provides the inverse sine Operator.

	This class represents the inverse sine function.
	*/
	class ArcSinOperator : public TrigOperator
	{
	public:
		BERTINI_DEFAULT_VISITABLE()

		std::shared_ptr<Node> Simplified() const override;
		std::shared_ptr<Node> Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const override;
		
		/// \brief Construct (and intern) a ArcSinOperator node.
		template<typename... Ts> 
		static 
		std::shared_ptr<ArcSinOperator> Make(Ts&& ...ts){ 
			return std::static_pointer_cast<ArcSinOperator>(Intern(std::shared_ptr<Node>( new ArcSinOperator(ts...) )));
		}

	private:
		ArcSinOperator(const std::shared_ptr<Node> & N) : TrigOperator(N)
		{};

	public:
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the sine function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;

		virtual ~ArcSinOperator() = default;
		
	protected:
		
		
		

		
		

	private:
		ArcSinOperator() = default;
		friend class boost::serialization::access;
		

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};




	
	
	
	
	/**
	\brief  Provides the cosine Operator.

	This class represents the cosine function.
	*/
	class CosOperator : public TrigOperator
	{
	public:
		BERTINI_DEFAULT_VISITABLE()
		
		std::shared_ptr<Node> Simplified() const override;
		std::shared_ptr<Node> Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const override;

		/// \brief Construct (and intern) a CosOperator node.
		template<typename... Ts> 
		static 
		std::shared_ptr<CosOperator> Make(Ts&& ...ts){ 
			return std::static_pointer_cast<CosOperator>(Intern(std::shared_ptr<Node>( new CosOperator(ts...) )));
		}

	private:
		CosOperator(const std::shared_ptr<Node> & N) : TrigOperator(N)
		{};
		
	public:
		
		
		
		void print(std::ostream & target) const override;
		
		
		
		
		/**
		 Differentiates the cosine function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;

		virtual ~CosOperator() = default;
		
	protected:
		
		
		
		
		
		

		
		
	private:
		CosOperator() = default;
		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};
	
	
	/**
	\brief Provides the arc cosine Operator.

	This class represents the inverse cosine function.
	*/
	class ArcCosOperator : public TrigOperator
	{
	public:
		BERTINI_DEFAULT_VISITABLE()

		std::shared_ptr<Node> Simplified() const override;
		std::shared_ptr<Node> Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const override;
		
		/// \brief Construct (and intern) a ArcCosOperator node.
		template<typename... Ts> 
		static 
		std::shared_ptr<ArcCosOperator> Make(Ts&& ...ts){ 
			return std::static_pointer_cast<ArcCosOperator>(Intern(std::shared_ptr<Node>( new ArcCosOperator(ts...) )));
		}

	private:
		ArcCosOperator(const std::shared_ptr<Node> & N) : TrigOperator(N)
		{};
	public:
		
		
		
		void print(std::ostream & target) const override;
		
		
		
		/**
		 Differentiates the cosine function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		virtual ~ArcCosOperator() = default;
		
	protected:
		
		
		
		
		
		

		
	private:
		ArcCosOperator() = default;
		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};
	
	
	
	
	


	
	
	
	/**
	\brief Provides the tangent Operator.

	This class represents the tangent function.
	*/
	class TanOperator : public TrigOperator
	{
	public:
		BERTINI_DEFAULT_VISITABLE()
		
		std::shared_ptr<Node> Simplified() const override;
		std::shared_ptr<Node> Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const override;
		
		/// \brief Construct (and intern) a TanOperator node.
		template<typename... Ts> 
		static 
		std::shared_ptr<TanOperator> Make(Ts&& ...ts){ 
			return std::static_pointer_cast<TanOperator>(Intern(std::shared_ptr<Node>( new TanOperator(ts...) )));
		}

	private:
		TanOperator(const std::shared_ptr<Node> & N) : TrigOperator(N)
		{};
	public:
		
		
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the tangent function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		virtual ~TanOperator() = default;
		
	protected:
		
		
		
		
		
		

		
	private:
		TanOperator() = default;
		friend class boost::serialization::access;
		

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};
	
	
	/**
	\brief Provides the inverse tangent Operator.

	This class represents the inverse tangent function.
	*/
	class ArcTanOperator : public TrigOperator
	{
	public:
		BERTINI_DEFAULT_VISITABLE()
		
		std::shared_ptr<Node> Simplified() const override;
		std::shared_ptr<Node> Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const override;
		
		/// \brief Construct (and intern) a ArcTanOperator node.
		template<typename... Ts> 
		static 
		std::shared_ptr<ArcTanOperator> Make(Ts&& ...ts){ 
			return std::static_pointer_cast<ArcTanOperator>(Intern(std::shared_ptr<Node>( new ArcTanOperator(ts...) )));
		}

	private:
		ArcTanOperator(const std::shared_ptr<Node> & N) : TrigOperator(N)
		{};
	public:
		
		
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the tangent function.
		 */		
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		virtual ~ArcTanOperator() = default;
		
	protected:
		
		
		
		
		
		
		
	private:
		ArcTanOperator() = default;
		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<TrigOperator>(*this);
		}
	};
	
	
	// begin the overload of operators

	/// \brief Build a sine expression-tree node.
	inline std::shared_ptr<Node> sin(const std::shared_ptr<Node> & N)
	{
		return SinOperator::Make(N);
	}
	
	/// \brief Build a arcsine expression-tree node.
	inline std::shared_ptr<Node> asin(const std::shared_ptr<Node> & N)
	{
		return ArcSinOperator::Make(N);
	}
	


	/// \brief Build a cosine expression-tree node.
	inline std::shared_ptr<Node> cos(const std::shared_ptr<Node> & N)
	{
		return CosOperator::Make(N);
	}



	/// \brief Build a arccosine expression-tree node.
	inline std::shared_ptr<Node> acos(const std::shared_ptr<Node> & N)
	{
		return ArcCosOperator::Make(N);
	}



	/// \brief Build a tangent expression-tree node.
	inline std::shared_ptr<Node> tan(const std::shared_ptr<Node> & N)
	{
		return TanOperator::Make(N);
	}


	
	/// \brief Build a arctangent expression-tree node.
	inline std::shared_ptr<Node> atan(const std::shared_ptr<Node> & N)
	{
		return ArcTanOperator::Make(N);
	}
	
	
} // re: namespace node
} // re: namespace bertini


#endif
