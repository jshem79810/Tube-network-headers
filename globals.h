#ifndef GLOBALS_H
#define GLOBALS_H

#include <iostream>
#include <vector>

//general functions

inline void abort_on_failure()   //generic abort function
{
	exit(EXIT_FAILURE);
}

template <typename T>
inline size_t sizet_round(const T & num)   //round number to size_t
{
	size_t nfloor = size_t(num);
	if(num - ((T) nfloor) >= 0.5)
	{
		return nfloor + 1;
	}
	else
	{
		return nfloor;
	}

}

template <typename T>
std::vector<T> vectorise(const T * value_array, const int & N)    //turn c array into vector
{
	std::vector<T> v;
	for(int n = 0; n < N; n++)
	{
		v.push_back(value_array[n]);
	}
	return v;
}



#endif