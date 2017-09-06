// a little bit of code to generate random real numbers, either integral or floating point.

#include <random>
#include <iostream>

#pragma once


namespace demo{
	
template<typename T>
struct Random
{
	static 
	T Generate()
	{
		static_assert(std::is_arithmetic<T>::value, "must use an arithmetic type");
		return Generate(std::is_integral<T>());
	}


private:

	static
	T Generate(std::true_type)
	{
		static std::random_device rd;
		static std::default_random_engine gen(rd());
		static std::uniform_int_distribution<T> dist(std::numeric_limits<T>::min(),std::numeric_limits<T>::max());
		return dist(gen); 
	} 

	static
	T Generate(std::false_type)
	{
		static std::random_device rd;
		static std::default_random_engine gen(rd());
		static std::uniform_real_distribution<T> dist(T{-1},T{1});
		return dist(gen); 
	} 

   
};

}

