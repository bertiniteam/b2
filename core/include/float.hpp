#ifndef BERTINI_FLOAT
#define BERTINI_FLOAT

#include <iostream>
#include <random>
#include <memory>
#include <sstream>
#include <cmath>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>


#include <random>
#include <cmath>
#include <exception>
#include <boost/random.hpp>




#include <eigen3/Eigen/Core>




namespace boost { namespace serialization {
	
	
	/**
	 Save a mpfr_float type to a boost archive.
	 */
	template <typename Archive>
	void save(Archive& ar, ::boost::multiprecision::backends::mpfr_float_backend<0> const& r, unsigned /*version*/)
	{
		std::string tmp = r.str(0, std::ios::fixed);
		ar & tmp;
	}
	
	/**
	 Load a mpfr_float type from a boost archive.
	 */
	template <typename Archive>
	void load(Archive& ar, ::boost::multiprecision::backends::mpfr_float_backend<0>& r, unsigned /*version*/)
	{
		std::string tmp;
		ar & tmp;
		r = tmp.c_str();
	}
	
} } // re: namespaces

BOOST_SERIALIZATION_SPLIT_FREE(::boost::multiprecision::backends::mpfr_float_backend<0>)







namespace bertini {
	

	
	using boost::multiprecision::mpfr_float;
	
	/**
	 \class the fundamental adaptive multiple precision floating point class for bertini 2.
	 
	 
	 Members of this class have pointers to a double and a multiple-precision number, which is populated as needed.  Operations are required to stay within the precision of two numbers, and (e.g.) addition of two bertini::Floats will result in a throw of a runtime_error.
	 
	 
	 
	 \note If trying to round-trip a bertini::Float to an ostream, you probably will have to use precision on the ostream higher than the Precision() of the Float.
	 
	 \note This class is Boost serialized.
	 
	 
	 @startuml{Float.ChangePrecision,scale=1/2}
	 
	 
	 namespace bertini{
	 
	 class Float {
	 
	 The Bertini class for floating point Adaptive Multiple Precision numbers.
	 ..
	 Will hold multiple precision and double type numbers, and have methods
	 for switching between the two, as well as changing the precision of the MP Float.
	 ==
	 data members
	 --
	 -std::unique_ptr<double> number_as_double_;
	 -std::unique_ptr<boost::multiprecision::mpfr_float> number_as_multiple_precision_;
	 -unsigned int current_precision_;
	 
	 
	 ==
	 setters
	 --
	 +void ChangePrecision(unsigned int num_bits)
	 sets the precision in terms of the number of bits desired.
	 switches between number types internally based on the desired precision.
	 if going from high to low, then allocate the double type, copy the value, and clear the mpfr type.
	 if going from low to high, then allocate the mpfr type to desired precision, copy the value, and clear the double.
	 ..
	 
	 
	 ==
	 getters
	 --
	 +unsigned int Precision() const; ///< get the current working precision in digits.  16 is double precision
	 
	 ==
	 operators
	 --
	 +Float operator = (const Float& right); // the assignment operator
	 
	 ==
	 friend operators
	 ..
	 +friend ostream & operator << (ostream & out, Float & n);
	 ..
	 +friend Float operator +(const Float & left, const Float & right);
	 ..
	 +friend Float operator -(const Float & left, const Float & right);
	 ..
	 +friend Float operator *(const Float & left, const Float & right);
	 ..
	 +friend Float operator /(const Float & left, const Float & right);
	 ..
	 +friend Float operator %(const Float & left, const Float & right);
	 
	 
	 ==
	 mpi methods
	 --
	 optionally compiled via a #ifdef construct
	 
	 int Send(unsigned int target)
	 ..
	 int Receive(unsigned int source)
	 ..
	 int Send(unsigned int target, communicator comm)
	 ..
	 int Receive(unsigned int source, communicator comm)
	 ..
	 int BcastSend()
	 ..
	 int BcastSend(communicator)
	 ..
	 int BcastReceive(source)
	 ..
	 int BcastReceive(source, communicator)
	 
	 
	 ==
	 serialization
	 --
	 see << >> operators
	 use Boost.Serialization for archiving
	 }
	 }
	 @enduml
	 
	 
	 
	 
	 
	 */
	class Float{
		
		
		
		
	private:
		
		/**
		 \brief The underlying multiple precision number.  
		 
		 is null if precision indicates double, otherwise populated.
		 */
		std::unique_ptr<mpfr_float> number_as_multiple_precision_;
		
		/**
		 \brief The underlying double precision number.
		 
		 is null if precision indicates high precision, otherwise populated.
		 */
		std::unique_ptr<double> number_as_double_;
		
		/**
		\brief The current precision, IN DIGITS, of the Float.
		 */
		unsigned int current_precision_;
		
		
		
		// let the boost serialization library have access to the private members of this class.
		friend class boost::serialization::access;
		

		
		/**
		 \brief save method for archiving a Float
		 */
		template<class Archive>
		void save(Archive & ar, const unsigned int version) const
		{
			// note, version is always the latest when saving
			
			ar & current_precision_;
			
			if (current_precision_<17) {
				ar & *number_as_double_;
			}
			else{
				ar & *number_as_multiple_precision_;
			}
		}
		
		
		/**
		 \brief save method for archiving a Float
		 */
		template<class Archive>
		void load(Archive & ar, const unsigned int version)
		{
			ar & current_precision_;
			if (current_precision_<17) {
				double temp;
				ar & temp;
				number_as_double_.reset(new double(temp));
				number_as_multiple_precision_.reset(nullptr);
			}
			else{
				number_as_multiple_precision_.reset(new boost::multiprecision::mpfr_float);
				number_as_multiple_precision_->precision(current_precision_);
				ar & *number_as_multiple_precision_;
				number_as_double_.reset(nullptr);
			}
		}
		
		
		
		
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		// this has to be here so that the Boost.Serialization library
		// knows that there are different methods for the serialize() method.
		
		
		
		
		
	public:
		
		
		
		
		/**
		 \brief Default constructor.
		 
		 The default constructor creates an AMP number with value 0, of precision equal to the current default precision, boost::multiprecision::mpfr_float::default_precision().
		 */
		Float() : current_precision_(mpfr_float::default_precision())
		{
			if (current_precision_<17) {
				number_as_double_.reset(new double(0));
				number_as_multiple_precision_.reset(nullptr);
			}
			else{
				number_as_double_.reset(nullptr);
				number_as_multiple_precision_.reset(new mpfr_float("0.0"));
			}
		}
		
		
		/**
		 \brief Construct Float from string, at default_precision.
		 
		 Construct a Float from a string.  make it a multiple precision, with precision equal to the current default precision.
		 */
		Float(std::string new_value) : current_precision_(mpfr_float::default_precision())
		{
			
			if (current_precision_<17) {
				number_as_double_.reset(new double(boost::lexical_cast<double>(new_value)));
				number_as_multiple_precision_.reset(nullptr);
			}
			else{
				number_as_double_.reset(nullptr);
				number_as_multiple_precision_.reset(new mpfr_float(new_value)); // this may fail if the string is not a valid number.  we let boost handle this fault.
			}
		}
		
		
		/**
		 \brief Construct a float from a string, at given precision (digits).
		 
		 Construct a Float from a string.  if it's short, make it a double.  Otherwise, make it a multiple precision, with precision equal to the current default precision.
		 
		 \param new_value The string from which to construct the number.
		 \param desired_precision The number of digits desired for the number.
		 
		 \note This function changes the default precision, and the resets it, threatening safety in a multithreaded environment.
		 */
		Float(std::string new_value, unsigned int desired_precision) : current_precision_(desired_precision)
		{
			unsigned int cached_precision = mpfr_float::default_precision();
			mpfr_float::default_precision(desired_precision);
			
			if (desired_precision<17) {
				char* end_of_num;
				number_as_double_.reset(new double( boost::lexical_cast<double>(new_value) ));
				number_as_multiple_precision_.reset(nullptr);// this may fail if the string is not a
			}
			else{
				
				number_as_double_.reset(nullptr);
				number_as_multiple_precision_.reset(new mpfr_float(new_value));// this may fail if the string is not a valid number.  we let boost handle this fault.
			}
			
			mpfr_float::default_precision(cached_precision);

		}
		
		
		
		/**
		 \brief Create a Float from a double.
		 */
		Float(double d) : current_precision_(mpfr_float::default_precision())
		{
			if (current_precision_<17) {
				number_as_double_.reset(new double(d));
				number_as_multiple_precision_.reset(nullptr);
			}
			else {
				number_as_double_.reset(nullptr);
				number_as_multiple_precision_.reset(new mpfr_float(d));
			}
		}
		
		
		
		/**
		 \brief Create a multiprecision Float from a boost::multiprecision::mpfr_float.
		 
		 The resulting number will have the precision equal to that of the input number.
		 */
		Float(const boost::multiprecision::mpfr_float & m) : current_precision_(m.precision())
		{
			number_as_multiple_precision_.reset(new mpfr_float(m));
			// the multiple precision field is left empty
		}
		
		
		
		
		/**
		 \brief Copy constructor.
		 */
		Float(const Float & new_value) : current_precision_(0)
		{
			this->ChangePrecision(new_value.Precision());
			if (current_precision_==0) {
				return;
			}
			else if (current_precision_<17) {
				*number_as_double_ = *(new_value.number_as_double_);
			}
			else{
				*number_as_multiple_precision_ = *(new_value.number_as_multiple_precision_);
			}
		}
		
		
		
		
		
		
		/**
		 \brief The swap function.
		 */
		friend void swap(Float & left, Float & right)
		{
			std::swap(left.current_precision_, right.current_precision_);
			std::swap(left.number_as_double_, right.number_as_double_);
			std::swap(left.number_as_multiple_precision_, right.number_as_multiple_precision_);
		}
		
		
		
		
		/**
		 \brief The assignment operator.
		 
		 This uses the copy-and-swap idiom.
		 
		 \return the result of the assignment
		 \param other Another Float from which to assign
		 */
				Float& operator = (Float other)
		{
			swap(*this,other);
			return *this;
		}
		
		
		
		/**
		 \brief The move constructor.
		 
		 Using the copy-and-swap idiom.
		 */
		Float(Float && other) : Float()
		{
			swap(*this, other);
		}
		
		
		
		
		
		/**
		 \brief The assignment operator, from a double.
		 
		 \return the result of the assignment
		 \param other a double from which to assign
		 
		 This function will throw a runtime_error if you are assigning from is not un-assigned, or is of a different precision.
		 */
		Float& operator = (double other)
		{
			if (this->Precision()!=16 && Precision()!=0) {
				throw std::runtime_error("trying to assign bertini::Float to be double, and previous precision is not 16.");
			}
			else{
				*number_as_double_ = other;
			}
			return *this;
		}
		
		
		
		/**
		 the assignment operator
		 
		 \return the result of the assignment
		 \param other a boost::multiprecision::mpfr_float from which to assign
		 
		 this operator uses the copy-and-swap idiom.
		 */
		Float& operator = (mpfr_float other)
		{
			if (this->Precision()!=other.precision() && this->Precision()!=0) {
				throw std::runtime_error("trying to assign bertini::Float to be boost::multiprecision::mpfr_float, and precisions do not match.");
			}
			else{
				std::swap(*number_as_multiple_precision_,other);
			}
			return *this;
		}
		
		
		
		
		/**
		 \brief Change the precision of a bertini::Float to that which you desire.
		 
		 \param new_precision The new precision, in digits, of the number.  if 0, will throw, because this is unacceptable.  Otherwise, sets the fields accordingly.  
		 
		 If the input precision is less than 17, the Float becomes a double.  Otherwise, it will be a multiple precision number.
		 
		 \note The new precision cannot be 0.
		 
		 
		 @startuml
				 start
				
				 if (new_precision==0) then (yes)
					:throw]
					 stop
				 endif
		 
	
				 if (this->Precision() <= 16) then (yes, currently double)
					if (new_precision <= 16) then (yes, new is double)
						stop
					else (new_precision>16)
						 :number_as_multiple_precision_.reset(new mp)]
						 :number_as_multiple_precision_.precision(new_precision)]
						 :number_as_multiple_precision_ = *number_as_double]
						 :number_as_double_.reset(nullptr)]
					endif
				 else (no, mp already)
					if (new_precision <= 16) then (yes, turn into a double)
						 :number_as_double_.reset(new double(*number_as_multiple_precision_))]
						 :number_as_multiple_precision_.reset(nullptr)]
					else
						:number_as_multiple_precision_.precision(new_precision)]
					endif
				 endif

				 :current_precision_ = new_precision]
				 stop
		 
		 @enduml
		 
		 
		 
		 */
		void ChangePrecision(unsigned int new_precision)
		{
			
			if (new_precision==0) {
				throw std::runtime_error("trying to set a bertini::Float with 0 precision.");
			}
			
			
			if (this->Precision()==0) {
				if (new_precision<17) {
					number_as_double_.reset(new double(0));
					number_as_multiple_precision_.release();
					new_precision = 16;
				}
				else{
					number_as_multiple_precision_.reset(new mpfr_float(0));
					number_as_multiple_precision_->precision(new_precision);
					number_as_double_.release();
				}
			}
			else
			{// if pre-existing precision is not 0; already set to something.
				
					if (this->Precision()<17) {
						if (new_precision<17) {
							return;
						}
						else{
							number_as_multiple_precision_.reset(new mpfr_float);
							number_as_multiple_precision_->precision(new_precision);
							*number_as_multiple_precision_ = *number_as_double_;
							number_as_double_.release();
						}
					}
					else{
						if (new_precision<17) {
							number_as_double_.reset(new double(*number_as_multiple_precision_));
							number_as_multiple_precision_.release();
							new_precision = 16;
						}
						else{
							//simply change the precision, as already mp, just not at same precision.
							number_as_multiple_precision_->precision(new_precision);
						}
					}
			}//re: if Precision==0
			 //done manually changing the fields to correct precision, now change the tracked number representing the precision.
			current_precision_ = new_precision;
		}
		
		
		/**
		 \brief Get the current precision of the number, in digits.
		 
		 \return the current number of DIGITS of the number
		 */
		unsigned int Precision() const
		{
			return current_precision_;
		}
		
		
		
		/**
		 \brief Turn the number into a randomly distributed number in [0,1).
		 
		 \note This calls a function which changes the default precision and back again, so is potentially unsafe in a multithreaded environment.
		 */
		Float MakeRandom()
		{
			if (Precision()==0) {
				throw std::runtime_error("trying to make a precision-0 number random.");
			}
			else if (Precision()<17) {
				*number_as_double_ = double(rand());
			}
			else{
				RandomLongNumberUniformUnitInterval(Precision(),*number_as_multiple_precision_);
			}
			return *this;
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		/***********************
		 **
		 **  comparitors
		 **
		 **********************/
		
		
		/**********
		 equality
		 ********/
		
		
		/**
		 \brief Test for equality.
		 */
		bool operator==(const Float & other) const
		{
			if (this->Precision()!= other.Precision()) {
				throw std::runtime_error("trying to compare == two Floats of differing precision");
			}
			else
			{ // the two precisions are the same
				if (Precision()<17) {
					if ( *(this->number_as_double_) == *(other.number_as_double_) ) { // dereference the pointers to compare the values
						return true;
					}
					else{
						return false;
					}
				}
				else {
					if ( *(this->number_as_multiple_precision_) == *(other.number_as_multiple_precision_) ) { // dereference the pointers to compare the values
						return true;
					}
					else{
						return false;
					}
				}
			}
		}
		
		
		
		/**
		 \brief Test for equality with a double, right hand side.
		 
		 Throws if the Float is not at precision 16.
		 */
		bool operator==(double rhs) const
		{
			if (this->Precision()!= 16) {
				throw std::runtime_error("trying to compare == non-double bertini::Float with double");
			}
			else
			{ // the two precisions are the same
				return (*number_as_double_==rhs);
			}
		}
		
		/**
		 \brief Test for equality with a double, left hand side.
		 
		 Throws if the Float is not at precision 16.
		 */
		friend bool operator==(double lhs, const Float & rhs)
		{
			if (rhs.Precision()!=16) {
				throw std::runtime_error("trying to compare == non-double bertini::Float with double");
			}
			else{
				return lhs== *rhs.number_as_double_;
			}
		}
		
		
		/**
		 \brief Test for equality with an boost::multiprecision::mpfr_float.
		 
		 Throws if the Float is not at same precision as other.
		 */
		bool operator==(const mpfr_float & other) const
		{
			if (this->Precision()!= other.precision()) {
				throw std::runtime_error("trying to compare == numbers of different precision or equality");
			}
			else
			{ // the two precisions are the same
				if (this->Precision()<17) {
					return (*number_as_double_==other);
				}
				else{
					return (*number_as_multiple_precision_==other);
				}
			}
		}
		
		/**
		 \brief Test for equality with an boost::multiprecision::mpfr_float.
		 
		 Throws if the Float is not at same precision as other
		 */
		friend bool operator==(const mpfr_float & lhs, const Float & rhs)
		{
			return (rhs==lhs);
		}
		
		
		
		
		
		/**
		 \brief Test for inequality.
		 */
		bool operator!=(const Float & other) const
		{
			return !(*this==other);
		}
		
		
		
		/**
		 \brief Test for inequality with a double.
		 */
		bool operator!=(double other) const
		{
			return !(*this==other);
		}
		
		/**
		 \brief Test for inequality with a double.
		 */
		friend bool operator!=(double lhs, const Float & rhs)
		{
			return !(rhs==lhs);
		}
		
		
		/**
		 \brief Test for inequality with an boost::multiprecision::mpfr_float.
		 */
		bool operator!=(const mpfr_float & other) const
		{
			return !(*this==other);
		}
		
		/**
		 \brief Test for inequality with an boost::multiprecision::mpfr_float.
		 */
		friend bool operator!=(const mpfr_float & lhs, const Float & rhs)
		{
			return !(rhs==lhs);
		}
		
		
		
		
		
		
		/**********
		 less than
		 ********/
		
		
		
		
		/**
		 \brief Test for less than.
		 
		 \note Permits comparison of numbers of differing precision.
		 
		 */
		bool operator<(const Float & other) const
		{
	
			if (this->Precision()<17 && other.Precision()<17) {
				return (*this->number_as_double_ < *other.number_as_double_);
			}
			else if (this->Precision()<17 && other.Precision()>=17) {
				return (*this->number_as_double_ < *other.number_as_multiple_precision_);
			}
			else if (this->Precision()>=17 && other.Precision()<17) {
				return (*this->number_as_multiple_precision_ < *other.number_as_double_);
			}
			else {
				return (*this->number_as_multiple_precision_ < *other.number_as_multiple_precision_);
			}

		}
		
		/**
		 \brief Test for less than with a double.
		 
		 \note Permits comparison of numbers of differing precision.
		 
		 */
		bool operator<(double other) const
		{
			if (this->Precision()<17) {
				return (*this->number_as_double_ < other);
			}
			else{
				return (*this->number_as_multiple_precision_ < other);
			}
		}
		
		/**
		 \brief Test for less than with an mpfr_float.
		 \note Permits comparison of numbers of differing precision.
		 */
		bool operator<(const mpfr_float & other) const
		{
			if (this->Precision()<17){
				return (*this->number_as_double_ < other);
			}
			else{
				return (*this->number_as_multiple_precision_ < other);
			}
		}
		
		
		/**
		 \brief Test for less than with a double on the left.
		 \note Permits comparison of numbers of differing precision.
		 */
		friend bool operator<(double lhs, const Float & rhs)
		{
			if (rhs.Precision()<17) {
				return lhs < *rhs.number_as_double_;
			}
			else{
				return lhs < *rhs.number_as_multiple_precision_;
			}
		}
		
		/**
		 \brief Test for less than with an mpfr_float on the left
		 \note Permits comparison of numbers of differing precision.
		 */
		friend bool operator<(const mpfr_float & lhs, const Float & rhs)
		{
			if (rhs.Precision()<17) {
				return lhs < *rhs.number_as_double_;
			}
			else{
				return lhs < *rhs.number_as_multiple_precision_;
			}
		}
		
		
		
		/*******
		 greater than
		 *******/
		
		/**
		 \brief Test for greater than.
		 \note Permits comparison of numbers of differing precision.
		 */
		bool operator>(const Float & other) const
		{
			if (this->Precision()<17 && other.Precision()<17) {
				return (*this->number_as_double_ > *other.number_as_double_);
			}
			else if (this->Precision()<17 && other.Precision()>=17) {
				return (*this->number_as_double_ > *other.number_as_multiple_precision_);
			}
			else if (this->Precision()>=17 && other.Precision()<17) {
				return (*this->number_as_multiple_precision_ > *other.number_as_double_);
			}
			else {
				return (*this->number_as_multiple_precision_ > *other.number_as_multiple_precision_);
			}
			
		}
		
		/**
		 \brief Test for greater than with a double.
		 \note Permits comparison of numbers of differing precision.
		 */
		bool operator>(double other) const
		{
			if (this->Precision()<17) {
				return (*this->number_as_double_ > other);
			}
			else{
				return (*this->number_as_multiple_precision_ > other);
			}
		}
		
		/**
		 \brief Test for greater than with an mpfr_float
		 \note Permits comparison of numbers of differing precision.
		 */
		bool operator>(const mpfr_float & other) const
		{
			if (this->Precision()<17){
				return (*this->number_as_double_ > other);
			}
			else{
				return (*this->number_as_multiple_precision_ > other);
			}
		}
		
		/**
		 \brief Test for greater than with a double on the left.
		 \note Permits comparison of numbers of differing precision.
		 */
		friend bool operator>(double lhs, const Float & rhs)
		{
			if (rhs.Precision()<17) {
				return lhs > *rhs.number_as_double_;
			}
			else{
				return lhs > *rhs.number_as_multiple_precision_;
			}
		}
		
		/**
		 \brief Test for greater than with an mpfr_float on the left.
		 \note Permits comparison of numbers of differing precision.
		 */
		friend bool operator>(const mpfr_float & lhs, const Float & rhs)
		{
			if (rhs.Precision()<17) {
				return lhs > *rhs.number_as_double_;
			}
			else{
				return lhs > *rhs.number_as_multiple_precision_;
			}
		}
		
		
		
		
		
		
		
		
		/**********
		 less than or equal to
		 ********/
		
		
		
		
		/**
		 \brief Test for less than or equal to.
		 \note Permits comparison of numbers of differing precision.
		 */
		bool operator<=(const Float & other) const
		{
			
			if (this->Precision()<17 && other.Precision()<17) {
				return (*this->number_as_double_ <= *other.number_as_double_);
			}
			else if (this->Precision()<17 && other.Precision()>=17) {
				return (*this->number_as_double_ <= *other.number_as_multiple_precision_);
			}
			else if (this->Precision()>=17 && other.Precision()<17) {
				return (*this->number_as_multiple_precision_ <= *other.number_as_double_);
			}
			else {
				return (*this->number_as_multiple_precision_ <= *other.number_as_multiple_precision_);
			}
			
		}
		
		/**
		 \brief Test for less than or equal to with a double.
		 \note Permits comparison of numbers of differing precision.
		 */
		bool operator<=(double other) const
		{
			if (this->Precision()<17) {
				return (*this->number_as_double_ <= other);
			}
			else{
				return (*this->number_as_multiple_precision_ <= other);
			}
		}
		
		/**
		 \brief Test for less than or equal to with an mpfr_float.
		 */
		bool operator<=(const mpfr_float & other) const
		{
			if (this->Precision()<17){
				return (*this->number_as_double_ <= other);
			}
			else{
				return (*this->number_as_multiple_precision_ <= other);
			}
		}
		
		friend bool operator<=(double lhs, const Float & rhs)
		{
			if (rhs.Precision()<17) {
				return lhs <= *rhs.number_as_double_;
			}
			else{
				return lhs <= *rhs.number_as_multiple_precision_;
			}
		}
		
		friend bool operator<=(const mpfr_float & lhs, const Float & rhs)
		{
			if (rhs.Precision()<17) {
				return lhs <= *rhs.number_as_double_;
			}
			else{
				return lhs <= *rhs.number_as_multiple_precision_;
			}
		}
		
		
		
		/*******
		 greater than or equal to
		 *******/
		
		/**
		 \brief Test for greater than than or equal to.
		 */
		bool operator>=(const Float & other) const
		{
			
			if (this->Precision()<17 && other.Precision()<17) {
				return (*this->number_as_double_ >= *other.number_as_double_);
			}
			else if (this->Precision()<17 && other.Precision()>=17) {
				return (*this->number_as_double_ >= *other.number_as_multiple_precision_);
			}
			else if (this->Precision()>=17 && other.Precision()<17) {
				return (*this->number_as_multiple_precision_ >= *other.number_as_double_);
			}
			else {
				return (*this->number_as_multiple_precision_ >= *other.number_as_multiple_precision_);
			}
			
		}
		
		
		bool operator>=(double other) const
		{
			if (this->Precision()<17) {
				return (*this->number_as_double_ >= other);
			}
			else{
				return (*this->number_as_multiple_precision_ >= other);
			}
		}
		
		bool operator>=(const mpfr_float & other) const
		{
			if (this->Precision()<17){
				return (*this->number_as_double_ >= other);
			}
			else{
				return (*this->number_as_multiple_precision_ >= other);
			}
		}
		
		friend bool operator>=(double lhs, const Float & rhs)
		{
			if (rhs.Precision()<17) {
				return lhs >= *rhs.number_as_double_;
			}
			else{
				return lhs >= *rhs.number_as_multiple_precision_;
			}
		}
		
		friend bool operator>=(const mpfr_float & lhs, const Float & rhs)
		{
			if (rhs.Precision()<17) {
				return lhs >= *rhs.number_as_double_;
			}
			else{
				return lhs >= *rhs.number_as_multiple_precision_;
			}
		}

		
		

		
		
		
		
		
		
		
		
		
		
		
		
		
		/******************
		 *
		 *  I/O streams
		 *
		 ****************/
		
		/**
		 \brief Send to an ostream.
		 */
		friend std::ostream & operator<<(std::ostream & out, const Float & printme)
		{
			if (printme.Precision()<17) {
				out << *printme.number_as_double_;
			}
			else{
				out << *printme.number_as_multiple_precision_;
			}
			
			return out;
		}
		
		/**
		 \brief Read from an istream.
		 */
		friend std::istream & operator>>(std::istream & in_stream, Float & writeontome)
		{
			if (writeontome.Precision()<=16) {
				in_stream >> *(writeontome.number_as_double_);
			}
			else{
				in_stream >> *writeontome.number_as_multiple_precision_;
			}
			return in_stream;
		}
		
		
		
		
		
		
		
		/*******************
		 **
		 **   basic arithmetic operations and operators
		 **
		 *****************/
		
		
		
		
		Float & operator+=(const Float & rhs)
		{
			if (this->Precision()==0) {
				throw std::runtime_error("trying to add two bertini::Floats and the left has 0 precision");
			}
			if (rhs.Precision()==0) {
				throw std::runtime_error("trying to add two bertini::Floats and the right has 0 precision");
			}
			
			if (this->Precision()!=rhs.Precision()) {
				throw std::runtime_error("trying to add two numbers of different precision");
			}
			else{
				if (this->Precision()<17){
					// both double, so add the double.
					*(this->number_as_double_) += *(rhs.number_as_double_);
				}
				else{
					*(this->number_as_multiple_precision_) += *(rhs.number_as_multiple_precision_);
				}
			}
			
			return *this;
		}
		
		
		
		
	

		
		
		
		
		
		
		
		
		
		
		
		/**
		 plain ol' subtraction
		 
		 Throws if the two numbers are not at the same precision.
		 */
		Float & operator-=(const Float & rhs)
		{
			if (this->Precision()==0) {
				throw std::runtime_error("trying to subtract bertini::Floats and the left has 0 precision");
			}
			if (rhs.Precision()==0) {
				throw std::runtime_error("trying to subtract bertini::Floats and the right has 0 precision");
			}
			
			if (this->Precision()!=rhs.Precision()) {
				throw std::runtime_error("trying to subtract numbers of different precision");
			}
			else{
				if (this->Precision()<17){
					// both double, so add the double.
					*(this->number_as_double_) -= *(rhs.number_as_double_);
				}
				else{
					*(this->number_as_multiple_precision_) -= *(rhs.number_as_multiple_precision_);
				}
			}
			
			return *this;
		}
		
		
		
		
		
		
		
		
		
		/**
		 \brief plain ol' negation
		 
		 Throws if the number is unset.
		 */
		Float operator-() const
		{
			if (this->Precision()==0) {
				throw std::runtime_error("trying to negate a bertini::Floats and it has 0 precision");
			}
			
			if (this->Precision()<17)
			{
				return bertini::Float(-*(this->number_as_double_));
			}
			else{
				return bertini::Float(-*(this->number_as_multiple_precision_));
			}
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		/**
		 \brief plain ol' multiplication
		 
		 \note Throws if the numbers are of differing precision.
		 */
		Float & operator*=(const Float & rhs)
		{
			if (this->Precision()==0) {
				throw std::runtime_error("trying to multiply  bertini::Floats and the left has 0 precision");
			}
			if (rhs.Precision()==0) {
				throw std::runtime_error("trying to multiply bertini::Floats and the right has 0 precision");
			}
			
			if (this->Precision()!=rhs.Precision()) {
				throw std::runtime_error("trying to multiply numbers of different precision");
			}
			else{
				if (this->Precision()<17){
					// both double, so add the double.
					*(this->number_as_double_) *= *(rhs.number_as_double_);
				}
				else{
					*(this->number_as_multiple_precision_) *= *(rhs.number_as_multiple_precision_);
				}
			}
			
			return *this;
		}
		
		
		
		

		
		
		
		
		/**
		 \brief plain ol' division
		 
		 \note Throws if the numbers are of differing precision.
		 */
		Float & operator/=(const Float & rhs)
		{
			if (this->Precision()==0) {
				throw std::runtime_error("trying to divide  bertini::Floats and the left has 0 precision");
			}
			if (rhs.Precision()==0) {
				throw std::runtime_error("trying to divide bertini::Floats and the right has 0 precision");
			}
			
			if (this->Precision()!=rhs.Precision()) {
				throw std::runtime_error("trying to divide numbers of different precision");
			}
			else{
				if (this->Precision()<17){
					// both double, so add the double.
					*(this->number_as_double_) /= *(rhs.number_as_double_);
				}
				else{
					*(this->number_as_multiple_precision_) /= *(rhs.number_as_multiple_precision_);
				}
			}
			
			return *this;
		}
		
		
		
		
		
		
		
		
		
		/**
		 \brief plain ol' addition
		 */
		inline friend Float operator+(Float lhs, const Float & rhs)
		{
			return lhs+=rhs;
		}
		
		
		
		
		
		
		
		
		/**
		 \brief plain ol' subtraction
		 */
		inline friend Float operator-(Float lhs, const Float & rhs)
		{
			return lhs-=rhs;
		}
		
		
		
		
		
		
		/**
		 \brief plain ol' multiplication
		 */
		inline friend Float operator*(Float lhs, const Float & rhs)
		{
			return lhs*=rhs;
		}
		
		
		
		
		
		
		/**
		 \brief plain ol' division
		 */
		inline friend Float operator/(Float lhs, const Float & rhs)
		{
			return lhs/=rhs;
		}
		

		
		
		
		
		
		/******************
		 **
		 **  non-elementary math operations
		 **
		 *******************/
		
		
		/**
		 \brief Support for the absolute value function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float abs(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take the absolute value of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( fabs(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take absolute value of bertini::Float at different precision from default.  Returned value would be of differing precision from input");
				}
				else
					return Float( fabs(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		

		
		/**
		 \brief Support for the square root function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float sqrt(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take the square root of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( sqrt(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take square root of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				else
					return Float( sqrt(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		
		/**
		 \brief Support for the floor function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float floor(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take the floor of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( floor(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take floor of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				else
					return Float( floor(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		
		/**
		 \brief Support for the ceiling function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float ceil(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take the ceiling of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( ceil(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take ceiling of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				else
					return Float( ceil(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		
		/**
		 \brief Support for the exponential function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float exp(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take exp of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( exp(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take exp of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				else
					return Float( exp(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		
		/**
		 \brief Support for the log function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float log(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take log of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( log(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take log of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				return Float( log(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		/**
		 \brief Support for the log10 function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float log10(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take log10 of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( log10(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take log10 of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				return Float( log10(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		/**
		 \brief Support for the cosine function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float cos(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take cosine of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( cos(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take cosine of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				return Float( cos(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		/**
		 \brief Support for the sine function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float sin(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take sine of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( sin(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take sine of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				
				return Float( sin(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		
		/**
		 \brief Support for the tangent function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float tan(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take tangent of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( tan(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take tangent of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				return Float( tan(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		/**
		 \brief Support for the arccosine function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float acos(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take arccosine of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( acos(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take arccosine of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				return Float( acos(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		/**
		 \brief Support for the arcsine function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float asin(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take arcsine of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( asin(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take arcsine of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				return Float( asin(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		
		/**
		 \brief Support for the arctangent function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float atan(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take arctangent of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( atan(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take arctangent of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				return Float( atan(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		
		/**
		 \brief Support for the hyperbolic cosine function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float cosh(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take hyperbolic cosine of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( cosh(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take hyperbolic cosine of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				return Float( cosh(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		
		/**
		 \brief Support for the hyperbolic sine function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float sinh(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take hyperbolic sine of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( sinh(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take hyperbolic sine of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				return Float( sinh(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		
		
		/**
		 \brief Support for the hyperbolic tangent function.
		 
		 \note Throws if the number is of multiple precision, and default precision is different from current precision.
		 */
		friend Float tanh(const Float & argument)
		{
			if (argument.Precision()==0) {
				throw std::runtime_error("trying to take hyperbolic tangent of a bertini::Float, and it has 0 precision");
			}
			else if (argument.Precision()<17){
				return Float( tanh(*(argument.number_as_double_)) );
			}
			else{
				if (argument.Precision()!=mpfr_float::default_precision()){
					throw std::runtime_error("trying to take hyperbolic tangent of bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
				}
				
				return Float( tanh(*(argument.number_as_multiple_precision_)) );
			}
		}
		
		
		
		/**
		 \brief Support for the power function.
		 
		 \note Throws if the base is of multiple precision, and default precision is different from current precision.
		 */
		friend Float pow(const Float & base, const Float & exponent)
		{
			if (base.Precision()==0) {
				throw std::runtime_error("trying to take pow of a bertini::Float, and base has 0 precision");
			}
			
			if (exponent.Precision()==0) {
				throw std::runtime_error("trying to take pow of a bertini::Float, and exponent has 0 precision");
			}
			
			if ( base.Precision()>16  &&  (base.Precision()!=mpfr_float::default_precision()) ){
				throw std::runtime_error("trying to take pow(,) of base bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
			}
			
			if (base.Precision()<17 && exponent.Precision()<17){
				return Float( std::pow(*base.number_as_double_,*exponent.number_as_double_) );
			}
			else if (base.Precision()<17 && exponent.Precision()>=17){
				return Float( pow(*base.number_as_double_,*exponent.number_as_multiple_precision_) );
			}
			else if (base.Precision()>=17 && exponent.Precision()<17){
				return Float( pow(*base.number_as_multiple_precision_,*exponent.number_as_double_) );
			}
			else{
				return Float( pow(*base.number_as_multiple_precision_,*exponent.number_as_multiple_precision_) );
			}
		}
		
		
		
		/**
		 \brief Support for the std pow function, with a double exponent
		 
		 \note Throws if the base is of different precision from current default.
		 */
		friend Float pow(const Float & base, double exponent)
		{
			if (base.Precision()==0) {
				throw std::runtime_error("trying to take pow of a bertini::Float, and base has 0 precision");
			}
			
			
			
			if ( base.Precision()>=17  &&  (base.Precision()!=mpfr_float::default_precision()) ){
				throw std::runtime_error("trying to take pow(,) of base bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
			}
			
			
			if (base.Precision()<17){
				return Float( std::pow(*base.number_as_double_,exponent) );
			}
			else{
				return Float( pow(*base.number_as_multiple_precision_,exponent) );
			}
		}
		
		
		
		
		/**
		 \brief Support for the std pow function, with an integer exponent
		 
		 \note Throws if the base is of different precision from current default.
		 */
		friend Float pow(const Float & base, int exponent)
		{
			if (base.Precision()==0) {
				throw std::runtime_error("trying to take pow of a bertini::Float, and base has 0 precision");
			}
			
			
			
			if ( base.Precision()>=17  &&  (base.Precision()!=mpfr_float::default_precision()) ){
				throw std::runtime_error("trying to take pow(,) of base bertini::Float at different precision from default.  Returned value would be of differing precision from input.  This is disallowed");
			}
			
			
			if (base.Precision()<17){
				return Float( std::pow(*base.number_as_double_,exponent) );
			}
			else{
				return Float( pow(*base.number_as_multiple_precision_,exponent) );
			}
		}
		
		
		
		
		/**
		 brief Support for atan2 operation
		 */
		friend Float atan2(const Float & left, const Float & right)
		{
			if (left.Precision()==0) {
				throw std::runtime_error("trying to take atan2 of a bertini::Float, and left has 0 precision");
			}
			
			if (right.Precision()==0) {
				throw std::runtime_error("trying to take atan2 of a bertini::Float, and right has 0 precision");
			}
			
			
			if (left.Precision()!=right.Precision()) {
				throw std::runtime_error("trying to take atan2 of two bertini::Floats at different precision");
			}
			
			
			if (left.Precision()<17){
				return Float( atan2(*left.number_as_double_,*right.number_as_double_) );
			}
			else{
				return Float( atan2(*left.number_as_multiple_precision_,*right.number_as_multiple_precision_) );
			}
		}
		
		
		
		
		
		
		
		
		
		/*******************
		 **
		 **   conversion operators
		 **
		 *****************/
		 
		 
		 
		
		
		
		
		friend Float copysign(const Float & lhs, const Float & rhs)
		{
			if (lhs.Precision()!=rhs.Precision()) {
				throw std::runtime_error("trying to mix precisions of bertini::Floats during copysign");
			}
			
			if (lhs.Precision()<17) {
				return Float(copysign(*lhs.number_as_double_,*rhs.number_as_double_));
			}
			else
			{
				return Float( abs(*lhs.number_as_multiple_precision_)* rhs.number_as_multiple_precision_->sign());
			}
			
		}
		
		
		/**
		 \brief query whether the number is nan
		 */
		friend bool isnan(const Float & x)
		{
			if (x.Precision()<=16){
				return std::isnan(*x.number_as_double_);
			}
			else{
				return boost::math::isnan(*x.number_as_multiple_precision_);
			}
		}
			
		/**
		 \brief query whether the number is infinity
		 */
		friend bool isinf(const Float & x)
		{
			if (x.Precision()<=16){
				return std::isinf(*x.number_as_double_);
			}
			else{
				return boost::math::isinf(*x.number_as_multiple_precision_);
			}
		}
			
			
		/**
		 \brief Convert Float to mpfr_float.
		 
		 returns nan if precision is 0.
		 */
		operator boost::multiprecision::mpfr_float() const
		{
			if (Precision()==0) {
				return mpfr_float(nan(""));
			}
			else if (Precision() <17){
				return mpfr_float(*number_as_double_);
			}
			else{
				return *number_as_multiple_precision_;
			}
			
		}
		
		
		
		/*****************
		 **
		 **   constants
		 **
		 ******************/
		
		
		
		/**
		 \brief get the number \f$\pi\f$, to default precision.
		 */
		static Float pi()
		{
			if (boost::multiprecision::mpfr_float::default_precision()<17) {
				return Float(acos(-1.0));
			}
			else{
				return Float(acos(boost::multiprecision::mpfr_float("-1.0")));
			}		}
		
		/**
		 \brief get the number \f$e\f$, to default precision.
		 */
		static Float e()
		{
			if (boost::multiprecision::mpfr_float::default_precision()<17) {
				return Float(exp(1.0));
			}
			else{
				return Float(exp(boost::multiprecision::mpfr_float("1.0")));
			}
			
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		/**
		 a templated function for producing random numbers in the unit interval, of a given number of digits.
		 
		 \tparam length_in_digits The length of the desired random number
		 \param number_to_make_random The number which you desire to populate with a random number.
		 
		 \note this function changes the default precision to the desired precision and back, threatening its use in a threaded environment which desires to have different threads create mpfr_floats of differing precision.
		 */
		template <unsigned int length_in_digits>
		static void RandomLongNumberSpecificDigitsUniformUnitInterval(boost::multiprecision::mpfr_float & number_to_make_random)
		{
			unsigned int prior_default_precision = boost::multiprecision::mpfr_float::default_precision();
			boost::multiprecision::mpfr_float::default_precision(number_to_make_random.precision());
			
			static boost::uniform_01<boost::multiprecision::mpfr_float> uf;
			static boost::random::independent_bits_engine<boost::random::mt19937, length_in_digits*1000L/301L, boost::multiprecision::mpz_int> gen;
			number_to_make_random = uf(gen);
			
			boost::multiprecision::mpfr_float::default_precision(prior_default_precision);
			
			
			return;
		}
		
		
		
		/**
		 \brief create a random number of a length at least a given number of digits
		 
		 \note this function calls the templated function RandomLongNumberSpecificDigitsUniformUnitInterval, which changes the default precision, then resets it.  This could possibly cause problems in a multithreaded environment.
		 
		 \param num_digits The desired minimum number of digits.
		 \param number_to_make_random The number whose contents you are overwriting with a random number.
		 */
		static void RandomLongNumberUniformUnitInterval(unsigned int num_digits, boost::multiprecision::mpfr_float & number_to_make_random)
		{
			if (num_digits<=50)
				return RandomLongNumberSpecificDigitsUniformUnitInterval<50>(number_to_make_random);
			else if (num_digits<=100)
				return RandomLongNumberSpecificDigitsUniformUnitInterval<100>(number_to_make_random);
			else if (num_digits<=200)
				return RandomLongNumberSpecificDigitsUniformUnitInterval<200>(number_to_make_random);
			else if (num_digits<=400)
				return RandomLongNumberSpecificDigitsUniformUnitInterval<400>(number_to_make_random);
			else if (num_digits<=800)
				return RandomLongNumberSpecificDigitsUniformUnitInterval<800>(number_to_make_random);
			else if (num_digits<=1600)
				return RandomLongNumberSpecificDigitsUniformUnitInterval<1600>(number_to_make_random);
			else if (num_digits<=3200)
				return RandomLongNumberSpecificDigitsUniformUnitInterval<3200>(number_to_make_random);
			else if (num_digits<=6400)
				return RandomLongNumberSpecificDigitsUniformUnitInterval<6400>(number_to_make_random);
			else
				throw std::out_of_range("requesting random long number, longer than provided by RandomLongNumber(d)");
		}
		
		
		
	}; // re: class Float
	
	
	
	
	
} //re: bertini namespaces

BOOST_CLASS_VERSION(bertini::Float, 0)


 namespace Eigen {
	 
	 /**
	  \brief this templated struct permits us to use the Float type in Eigen matrices.
	  */
	 template<> struct NumTraits<bertini::Float>  // permits to get the epsilon, dummy_precision, lowest, highest functions
	 {
		 typedef bertini::Float Real;
		 typedef bertini::Float NonInteger;
		 typedef bertini::Float Nested;
		 enum {
			 IsComplex = 0,
			 IsInteger = 0,
			 IsSigned = 1,
			 RequireInitialization = 1, // our default constructor is nonempty.  therefore, requires initialization.
			 ReadCost = 1,
			 AddCost = 3,
			 MulCost = 3
		 };
		 
		 
		 
		 //TODO verify this returned number
		 static Real epsilon()
		 {
			 return Real( pow( bertini::Float(10.0), 1-int(boost::multiprecision::mpfr_float::default_precision()) ));
		 }
		 
	 };
 }


#endif


